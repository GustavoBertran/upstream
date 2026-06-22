"""htsprep — HTS preprocessing pipeline

Usage:
  htsprep qc         --samples samples.csv --outdir results/
  htsprep rnaseq     --samples samples.csv --star-index /ref/star --salmon-index /ref/salmon --outdir results/
  htsprep atacseq    --samples samples.csv --bowtie2-index /ref/bt2/hg38 --outdir results/
  htsprep chipseq    --samples samples.csv --bowtie2-index /ref/bt2/hg38 --peak-type narrow --outdir results/
  htsprep methylation --samples samples.csv --bismark-genome /ref/bismark --outdir results/

Explanations are shown by default. Use --no-explain to suppress them.

Samplesheet CSV format (samples.csv):
  name,group,r1,r2
  GSM123,disease,/path/sample1_R1.fastq.gz,/path/sample1_R2.fastq.gz
  GSM456,control,/path/sample2_R1.fastq.gz,/path/sample2_R2.fastq.gz
"""

from __future__ import annotations

import csv
import importlib.resources
import sys
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel

from . import checkpoints, downloader, obama, runner

app = typer.Typer(
    name="htsprep",
    help="HTS preprocessing pipeline — trim, align, quantify, and produce OBAMA-format matrices.",
    no_args_is_help=True,
    add_completion=False,
)
console = Console()

# ── Types (shorthand for Annotated options) ──────────────────────────────

SamplesOpt = Annotated[Path, typer.Option("--samples", help="Samplesheet CSV (name,group,r1,r2)")]
OutdirOpt  = Annotated[Path, typer.Option("--outdir",  help="Output directory (created if missing)")]
ThreadsOpt = Annotated[int,  typer.Option("--threads", help="CPU threads for tools that support it")]
ExplainOpt = Annotated[bool, typer.Option("--explain/--no-explain",
                                           help="Show educational explanations (on by default)")]


# ── Helpers ───────────────────────────────────────────────────────────────


def _read_samplesheet(path: Path, allow_input: bool = False) -> list[dict[str, str]]:
    if not path.exists():
        _die(f"Samplesheet not found: {path}")
    with path.open() as f:
        reader = csv.DictReader(f)
        missing_cols = {"name", "group", "r1", "r2"} - set(reader.fieldnames or [])
        if missing_cols:
            _die(f"Samplesheet is missing columns: {', '.join(sorted(missing_cols))}")
        rows = list(reader)

    valid_groups = {"disease", "control"} | ({"input"} if allow_input else set())
    for row in rows:
        if row["group"] not in valid_groups:
            extra = " or 'input'" if allow_input else ""
            _die(
                f"Invalid group '{row['group']}' for sample '{row['name']}'. "
                f"Must be 'disease' or 'control'{extra} (OBAMA requires the exact strings "
                "'disease'/'control')."
            )
        # r1 is always required; r2 is optional (empty = single-end).
        if not Path(row["r1"]).exists():
            _die(f"File not found for sample '{row['name']}': {row['r1']}")
        r2 = (row.get("r2") or "").strip()
        if r2 and not Path(r2).exists():
            _die(f"File not found for sample '{row['name']}': {r2}")
    return rows


def _is_paired(s: dict) -> bool:
    """A sample is paired-end when it has a non-empty r2 column."""
    return bool((s.get("r2") or "").strip())


def _fastp_cmd(
    s: dict, name: str, trim_dir: Path, threads: int, extra: list[str],
) -> tuple[list[str], Path, Optional[Path]]:
    """Build the fastp command for a sample, handling paired vs single-end.

    Returns (cmd, trimmed_r1, trimmed_r2_or_None).  `--detect_adapter_for_pe`
    in *extra* is dropped for single-end input (fastp auto-detects SE adapters).
    """
    paired = _is_paired(s)
    out1 = trim_dir / f"{name}_R1.fastq.gz"
    out2 = trim_dir / f"{name}_R2.fastq.gz" if paired else None

    cmd = ["fastp", "--in1", str(Path(s["r1"])), "--out1", str(out1)]
    if paired:
        cmd += ["--in2", str(Path(s["r2"])), "--out2", str(out2)]
    cmd += [
        "--json", str(trim_dir / f"{name}_fastp.json"),
        "--html", str(trim_dir / f"{name}_fastp.html"),
        "--thread", str(threads),
    ]
    cmd += [f for f in extra if not (f == "--detect_adapter_for_pe" and not paired)]
    return cmd, out1, out2


def _build_peak_outputs(
    results: list[tuple[str, str, Path, Path]],  # (name, group, peak_file, filtered_bam)
    outdir: Path,
    output_format: str,
    threads: int,
    explain: bool,
    export_doc: str,
) -> list[Path]:
    """Write peak-based outputs (ATAC-seq / ChIP-seq), returning files written.

    'obama' → peak-score matrix (samples × peaks). 'matrix' → consensus peak set +
    per-peak read counts (deeptools multiBamSummary) → counts_matrix.csv + coldata.csv
    for DESeq2/edgeR. 'both' → both.
    """
    written: list[Path] = []

    if output_format in ("obama", "both"):
        obama_path = outdir / "obama_matrix.csv"
        obama.build_atacseq([(n, g, pf) for n, g, pf, _b in results], obama_path)
        _check(*checkpoints.check_obama_format(obama_path))
        written.append(obama_path)

    if output_format in ("matrix", "both"):
        if explain:
            _explain(export_doc)
        consensus_bed = outdir / "consensus_peaks.bed"
        n_peaks = obama.build_atac_consensus([r[2] for r in results], consensus_bed)
        console.print(f"  Consensus peak set: {n_peaks:,} merged regions → {consensus_bed.name}")

        raw = outdir / "_multibamsummary_raw.tab"
        npz = outdir / "_multibamsummary.npz"
        rc = runner.run([
            "multiBamSummary", "BED-file",
            "--BED",          str(consensus_bed),
            "--bamfiles",     *[str(r[3]) for r in results],
            "--labels",       *[r[0] for r in results],
            "-p",             str(threads),
            "--outRawCounts", str(raw),
            "-o",             str(npz),
        ])
        if rc != 0:
            _die("multiBamSummary (deeptools) failed — are the filtered BAMs indexed?")
        counts_path = outdir / "counts_matrix.csv"
        coldata_path = outdir / "coldata.csv"
        obama.build_atac_counts_matrix(
            raw, [(n, g) for n, g, _pf, _b in results], counts_path, coldata_path,
        )
        npz.unlink(missing_ok=True)
        raw.unlink(missing_ok=True)
        _check(*checkpoints.check_counts_matrix(counts_path, coldata_path))
        written += [counts_path, coldata_path]

    return written


def _die(msg: str) -> None:
    console.print(f"[bold red]Error:[/bold red] {msg}")
    raise typer.Exit(1)


def _require_dir(path: Path, label: str) -> None:
    if not path.exists():
        _die(f"{label} not found: {path}")


def _check(ok: bool, msg: str) -> None:
    if ok:
        runner.ok(msg)
    else:
        runner.fail(msg)
        raise typer.Exit(1)


def _explain(content_file: str) -> None:
    try:
        pkg = importlib.resources.files("upstream") / "content" / content_file
        text = pkg.read_text(encoding="utf-8")
        console.print(Panel(Markdown(text), border_style="dim blue", padding=(1, 2)))
    except Exception:
        pass  # explanations are optional; missing file is not an error


# ── QC ────────────────────────────────────────────────────────────────────


@app.command()
def qc(
    samples: SamplesOpt,
    outdir: OutdirOpt,
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """FastQC + MultiQC quality control on all samples in the samplesheet."""
    sample_list = _read_samplesheet(samples)
    outdir.mkdir(parents=True, exist_ok=True)
    console.print(f"\n[bold]QC pipeline[/bold] — {len(sample_list)} sample(s) → {outdir}\n")

    fastqc_out = outdir / "fastqc"
    fastqc_out.mkdir(exist_ok=True)
    multiqc_out = outdir / "multiqc"
    multiqc_out.mkdir(exist_ok=True)

    runner.step_header("FastQC", 1, 2)
    if explain:
        _explain("qc_fastqc.md")
    all_fastq = []
    for s in sample_list:
        all_fastq.append(s["r1"])
        if _is_paired(s):
            all_fastq.append(s["r2"])
    rc = runner.run(["fastqc", "--outdir", str(fastqc_out), "--threads", str(threads)] + all_fastq)
    if rc != 0:
        _die("FastQC failed.")
    _check(*checkpoints.check_fastqc_output(fastqc_out))

    runner.step_header("MultiQC", 2, 2)
    if explain:
        _explain("qc_multiqc.md")
    rc = runner.run(["multiqc", str(fastqc_out), "--outdir", str(multiqc_out), "--force"])
    if rc != 0:
        _die("MultiQC failed.")
    _check(*checkpoints.check_multiqc_output(multiqc_out))

    console.print(f"\n[bold green]QC complete.[/bold green] Open {multiqc_out}/multiqc_report.html")


# ── RNA-seq ───────────────────────────────────────────────────────────────


@app.command()
def rnaseq(
    samples: SamplesOpt,
    outdir: OutdirOpt,
    aligner: Annotated[str, typer.Option(
        "--aligner",
        help="'salmon' (alignment-free, default) or 'star' (genome alignment + gene counts).",
    )] = "salmon",
    salmon_index: Annotated[Optional[Path], typer.Option(
        "--salmon-index", help="Pre-built Salmon index directory (required for --aligner salmon).",
    )] = None,
    star_index: Annotated[Optional[Path], typer.Option(
        "--star-index", help="Pre-built STAR genome index directory (required for --aligner star).",
    )] = None,
    output_format: Annotated[str, typer.Option(
        "--format",
        help="Output: 'obama' (default), 'matrix' (counts_matrix.csv + coldata.csv for "
             "DESeq2/edgeR/limma-voom), or 'both'.",
    )] = "obama",
    gtf: Annotated[Optional[Path], typer.Option(
        "--gtf", help="GTF (GENCODE/Ensembl) for Salmon transcript→gene aggregation → gene-level counts.",
    )] = None,
    tx2gene: Annotated[Optional[Path], typer.Option(
        "--tx2gene", help="2-column CSV (transcript_id,gene) for Salmon gene-level aggregation. Overrides --gtf.",
    )] = None,
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """RNA-seq: trim (fastp) → align/quantify (STAR or Salmon) → counts matrix.

    Use --aligner salmon (default) for alignment-free quantification via Salmon, or
    --aligner star for genome alignment with STAR (gene counts via --quantMode
    GeneCounts). Both produce RAW COUNTS — Salmon NumReads, STAR ReadsPerGene.

    Output is the OBAMA matrix by default; --format matrix (or both) also writes a
    feature x sample counts_matrix.csv + coldata.csv for DESeq2/edgeR/limma-voom.
    For gene-level Salmon counts, pass --gtf (or --tx2gene); otherwise Salmon output
    stays transcript-level.
    """
    if aligner not in ("salmon", "star"):
        _die("--aligner must be 'salmon' or 'star'.")
    if aligner == "salmon" and salmon_index is None:
        _die("--salmon-index is required when --aligner salmon.")
    if aligner == "star" and star_index is None:
        _die("--star-index is required when --aligner star.")
    if output_format not in ("obama", "matrix", "both"):
        _die("--format must be 'obama', 'matrix', or 'both'.")
    if gtf is not None and not gtf.exists():
        _die(f"GTF not found: {gtf}")
    if tx2gene is not None and not tx2gene.exists():
        _die(f"tx2gene file not found: {tx2gene}")

    sample_list = _read_samplesheet(samples)
    if aligner == "salmon":
        _require_dir(salmon_index, "Salmon index")
    else:
        _require_dir(star_index, "STAR index")
    outdir.mkdir(parents=True, exist_ok=True)

    TOTAL = 3
    console.print(
        f"\n[bold]RNA-seq pipeline[/bold] ({aligner}) — "
        f"{len(sample_list)} sample(s) → {outdir}\n"
    )

    result_dirs: list[tuple[str, str, Path]] = []

    for s in sample_list:
        name, group = s["name"], s["group"]
        paired = _is_paired(s)
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)

        layout = "paired-end" if paired else "single-end"
        console.rule(f"[bold white]Sample: {name}  ({group})  [{layout}][/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("rnaseq_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        cmd, t_r1, t_r2 = _fastp_cmd(s, name, trim_dir, threads, ["--detect_adapter_for_pe"])
        rc = runner.run(cmd)
        if rc != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        if aligner == "salmon":
            # 2 — Salmon quantification (alignment-free)
            runner.step_header("Quantify — Salmon", 2, TOTAL)
            if explain:
                _explain("rnaseq_quantify.md")
            quant_dir = sdir / "salmon"
            reads = (["--mates1", str(t_r1), "--mates2", str(t_r2)] if paired
                     else ["--unmatedReads", str(t_r1)])
            rc = runner.run([
                "salmon", "quant",
                "--index",    str(salmon_index),
                "--libType",  "A",
                *reads,
                "--threads",  str(threads),
                "--output",   str(quant_dir),
                "--validateMappings",
            ])
            if rc != 0:
                _die(f"Salmon quantification failed for sample '{name}'.")
            _check(*checkpoints.check_salmon_sf(quant_dir / "quant.sf"))
            result_dirs.append((name, group, quant_dir))

        else:
            # 2 — STAR alignment with built-in gene counts
            runner.step_header("Align — STAR", 2, TOTAL)
            if explain:
                _explain("rnaseq_align.md")
            aln_dir = sdir / "aligned"
            aln_dir.mkdir(exist_ok=True)
            read_files = [str(t_r1)] + ([str(t_r2)] if paired else [])
            rc = runner.run([
                "STAR",
                "--runThreadN",       str(threads),
                "--genomeDir",        str(star_index),
                "--readFilesIn",      *read_files,
                "--readFilesCommand", "zcat",
                "--outSAMtype",       "BAM", "SortedByCoordinate",
                "--outSAMattributes", "NH", "HI", "AS", "NM",
                "--outFileNamePrefix", str(aln_dir / f"{name}_"),
                "--quantMode",        "GeneCounts",
                "--runRNGseed",       "42",
            ])
            if rc != 0:
                _die(f"STAR alignment failed for sample '{name}'.")
            _check(*checkpoints.check_star_bam(
                aln_dir / f"{name}_Aligned.sortedByCoord.out.bam",
                aln_dir / f"{name}_Log.final.out",
            ))
            result_dirs.append((name, group, aln_dir / f"{name}_ReadsPerGene.out.tab"))

    # 3 — Build output matrices (OBAMA and/or DESeq2/edgeR/limma matrix)
    runner.step_header("Build matrix output", 3, TOTAL)
    if explain and output_format in ("matrix", "both"):
        _explain("rnaseq_export_formats.md")

    source = "salmon" if aligner == "salmon" else "star"
    tx2gene_map = None
    if source == "salmon":
        if gtf or tx2gene:
            tx2gene_map = obama.load_tx2gene(gtf, tx2gene)
            if tx2gene_map is None:
                _die("Could not build a transcript→gene map from the provided --gtf/--tx2gene.")
            matched, total = obama.salmon_tx2gene_match(result_dirs[0][2], tx2gene_map)
            pct = (matched / total * 100) if total else 0
            console.print(
                f"  Transcript→gene map: {len(tx2gene_map):,} mappings; matched "
                f"{matched:,}/{total:,} quantified transcripts ({pct:.0f}%)."
            )
            if pct < 50:
                console.print(
                    "[yellow]Warning:[/yellow] under half of quantified transcripts matched the "
                    "tx2gene map. Was the Salmon index built from the same GENCODE/Ensembl release "
                    "(versioned IDs)? Output may stay largely transcript-level."
                )
        elif output_format in ("matrix", "both"):
            console.print(
                "[yellow]Note:[/yellow] no --gtf/--tx2gene given — Salmon counts stay "
                "transcript-level. DESeq2/edgeR are usually run at gene level."
            )

    written = obama.write_rnaseq_outputs(result_dirs, outdir, output_format, source, tx2gene_map)

    if output_format in ("obama", "both"):
        _check(*checkpoints.check_obama_format(outdir / "obama_matrix.csv"))
    if output_format in ("matrix", "both"):
        _check(*checkpoints.check_counts_matrix(outdir / "counts_matrix.csv", outdir / "coldata.csv"))

    console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── ATAC-seq ──────────────────────────────────────────────────────────────


@app.command()
def atacseq(
    samples: SamplesOpt,
    bowtie2_index: Annotated[Path, typer.Option("--bowtie2-index", help="Bowtie2 index prefix (path/to/hg38)")],
    outdir: OutdirOpt,
    output_format: Annotated[str, typer.Option(
        "--format",
        help="Output: 'obama' (default), 'matrix' (counts_matrix.csv + coldata.csv — reads "
             "counted in a consensus peak set, for DESeq2/edgeR), or 'both'.",
    )] = "obama",
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """ATAC-seq: QC → trim (fastp) → align (Bowtie2) → filter → peaks (MACS2) → matrix.

    Default output is the OBAMA matrix (MACS2 peak scores). --format matrix (or both)
    instead builds a consensus peak set across samples and counts reads per peak with
    deeptools multiBamSummary, producing a peak × sample raw-count matrix + coldata for
    DESeq2/edgeR differential accessibility.
    """
    if output_format not in ("obama", "matrix", "both"):
        _die("--format must be 'obama', 'matrix', or 'both'.")
    sample_list = _read_samplesheet(samples)
    bt2_prefix = bowtie2_index  # e.g. /ref/bowtie2/hg38 (no .bt2 extension)
    outdir.mkdir(parents=True, exist_ok=True)

    TOTAL = 5
    console.print(f"\n[bold]ATAC-seq pipeline[/bold] — {len(sample_list)} sample(s) → {outdir}\n")

    # (name, group, peak_file, filtered_bam) per sample
    atac_results: list[tuple[str, str, Path, Path]] = []

    for s in sample_list:
        name, group = s["name"], s["group"]
        paired = _is_paired(s)
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)

        layout = "paired-end" if paired else "single-end"
        console.rule(f"[bold white]Sample: {name}  ({group})  [{layout}][/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("atacseq_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        cmd, t_r1, t_r2 = _fastp_cmd(s, name, trim_dir, threads, ["--detect_adapter_for_pe"])
        rc = runner.run(cmd)
        if rc != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        # 2 — Align
        runner.step_header("Align — Bowtie2", 2, TOTAL)
        if explain:
            _explain("atacseq_align.md")
        aln_dir = sdir / "aligned"
        aln_dir.mkdir(exist_ok=True)
        bam_raw = aln_dir / f"{name}.bam"

        # Paired-end uses concordant-pair flags (-X/--no-mixed/--no-discordant);
        # single-end aligns unpaired reads with -U.
        reads = (["-1", str(t_r1), "-2", str(t_r2),
                  "--no-mixed", "--no-discordant", "-X", "2000"] if paired
                 else ["-U", str(t_r1)])
        rc = runner.pipe(
            ["bowtie2",
             "--threads", str(threads),
             "-x", str(bt2_prefix),
             *reads,
             "--very-sensitive"],
            ["samtools", "sort", "-o", str(bam_raw), "-"],
        )
        if rc != 0:
            _die(f"Bowtie2/samtools alignment failed for sample '{name}'.")
        runner.run(["samtools", "index", str(bam_raw)])
        _check(*checkpoints.check_bowtie2_bam(bam_raw))

        # 3 — Filter
        runner.step_header("Filter — remove mito, dedup, MAPQ ≥ 30", 3, TOTAL)
        if explain:
            _explain("atacseq_filter.md")
        filt_dir = sdir / "filtered"
        filt_dir.mkdir(exist_ok=True)

        autosomes = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
        bam_nomito = filt_dir / f"{name}.nomito.bam"
        bam_dedup  = filt_dir / f"{name}.dedup.bam"
        bam_final  = filt_dir / f"{name}.filtered.bam"

        for cmd in [
            ["samtools", "view", "-b", str(bam_raw)] + autosomes + ["-o", str(bam_nomito)],
            ["samtools", "index", str(bam_nomito)],
            ["samtools", "markdup", "-r", str(bam_nomito), str(bam_dedup)],
            ["samtools", "index", str(bam_dedup)],
            ["samtools", "view", "-b", "-q", "30", str(bam_dedup), "-o", str(bam_final)],
            ["samtools", "index", str(bam_final)],
        ]:
            if runner.run(cmd) != 0:
                _die(f"Filtering failed at: {' '.join(str(c) for c in cmd[:3])}")
        _check(*checkpoints.check_atac_filtered_bam(bam_final))

        # 4 — Peaks
        runner.step_header("Peak calling — MACS2", 4, TOTAL)
        if explain:
            _explain("atacseq_peaks.md")
        peaks_dir = sdir / "peaks"
        peaks_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "macs2", "callpeak",
            "-t", str(bam_final),
            "-f", "BAM",
            "--nomodel", "--shift", "-100", "--extsize", "200",
            "-g", "hs",
            "--outdir", str(peaks_dir),
            "-n", name,
        ])
        if rc != 0:
            _die(f"MACS2 peak calling failed for sample '{name}'.")
        peak_file = peaks_dir / f"{name}_peaks.narrowPeak"
        _check(*checkpoints.check_macs2_peaks(peak_file))

        atac_results.append((name, group, peak_file, bam_final))

    # 5 — Build output matrices
    runner.step_header("Build matrix output", 5, TOTAL)
    written = _build_peak_outputs(atac_results, outdir, output_format, threads,
                                  explain, "atacseq_export_formats.md")
    console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── ChIP-seq ─────────────────────────────────────────────────────────────────


@app.command()
def chipseq(
    samples: SamplesOpt,
    bowtie2_index: Annotated[Path, typer.Option("--bowtie2-index", help="Bowtie2 index prefix (path/to/hg38)")],
    outdir: OutdirOpt,
    peak_type: Annotated[str, typer.Option(
        "--peak-type",
        help="'narrow' (default; TFs, H3K4me3, H3K27ac) or 'broad' (H3K27me3, H3K9me3, H3K36me3).",
    )] = "narrow",
    output_format: Annotated[str, typer.Option(
        "--format",
        help="Output: 'matrix' (default; consensus-peak counts + coldata for DESeq2/edgeR), "
             "'obama' (peak-score matrix — experimental for peaks: OBAMA's gene-based "
             "interpretation needs peak→gene annotation), or 'both'.",
    )] = "matrix",
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """ChIP-seq: fastp → Bowtie2 → filter → MACS2 (optional input control) → matrix.

    The samplesheet may add an optional 'control' column naming each ChIP sample's
    input by 'name'; rows with group=input are aligned to provide that control BAM
    but are not peak-called or placed in the matrix. --peak-type broad calls broad
    domains (histone marks). Default output is a consensus-peak count matrix for
    DESeq2/edgeR differential binding (the validated downstream for peak data;
    OBAMA's gene-centric modules need peak→gene annotation, so --format obama is
    experimental for ChIP).
    """
    if peak_type not in ("narrow", "broad"):
        _die("--peak-type must be 'narrow' or 'broad'.")
    if output_format not in ("obama", "matrix", "both"):
        _die("--format must be 'obama', 'matrix', or 'both'.")
    broad = peak_type == "broad"

    sample_list = _read_samplesheet(samples, allow_input=True)
    inputs = {s["name"]: s for s in sample_list if s["group"] == "input"}
    signal = [s for s in sample_list if s["group"] != "input"]
    if not signal:
        _die("No ChIP signal samples found (every row is group=input).")
    for s in signal:
        ctl = (s.get("control") or "").strip()
        if ctl and ctl not in inputs:
            _die(f"control '{ctl}' for sample '{s['name']}' must name a group=input row.")

    bt2_prefix = bowtie2_index
    outdir.mkdir(parents=True, exist_ok=True)

    TOTAL = 5
    n_in = len(inputs)
    console.print(
        f"\n[bold]ChIP-seq pipeline[/bold] ({peak_type}) — {len(signal)} ChIP sample(s)"
        + (f" + {n_in} input(s)" if n_in else " (no input control)")
        + f" → {outdir}\n"
    )

    # Align + filter every sample, inputs first so control BAMs exist before peak-calling.
    ordered = list(inputs.values()) + signal
    bam_for: dict[str, Path] = {}
    chip_results: list[tuple[str, str, Path, Path]] = []

    for s in ordered:
        name, group = s["name"], s["group"]
        paired = _is_paired(s)
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)
        layout = "paired-end" if paired else "single-end"
        role = "input/control" if group == "input" else group
        console.rule(f"[bold white]Sample: {name}  ({role})  [{layout}][/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("chipseq_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        cmd, t_r1, t_r2 = _fastp_cmd(s, name, trim_dir, threads, ["--detect_adapter_for_pe"])
        if runner.run(cmd) != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        # 2 — Align
        runner.step_header("Align — Bowtie2", 2, TOTAL)
        if explain:
            _explain("chipseq_align.md")
        aln_dir = sdir / "aligned"
        aln_dir.mkdir(exist_ok=True)
        bam_raw = aln_dir / f"{name}.bam"
        reads = (["-1", str(t_r1), "-2", str(t_r2), "--no-mixed", "--no-discordant", "-X", "2000"]
                 if paired else ["-U", str(t_r1)])
        rc = runner.pipe(
            ["bowtie2", "--threads", str(threads), "-x", str(bt2_prefix), *reads, "--very-sensitive"],
            ["samtools", "sort", "-o", str(bam_raw), "-"],
        )
        if rc != 0:
            _die(f"Bowtie2/samtools alignment failed for sample '{name}'.")
        runner.run(["samtools", "index", str(bam_raw)])
        _check(*checkpoints.check_bowtie2_bam(bam_raw))

        # 3 — Filter
        runner.step_header("Filter — dedup, MAPQ ≥ 30, remove mito", 3, TOTAL)
        if explain:
            _explain("chipseq_filter.md")
        filt_dir = sdir / "filtered"
        filt_dir.mkdir(exist_ok=True)
        autosomes = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
        bam_nomito = filt_dir / f"{name}.nomito.bam"
        bam_dedup  = filt_dir / f"{name}.dedup.bam"
        bam_final  = filt_dir / f"{name}.filtered.bam"
        for cmd in [
            ["samtools", "view", "-b", str(bam_raw)] + autosomes + ["-o", str(bam_nomito)],
            ["samtools", "index", str(bam_nomito)],
            ["samtools", "markdup", "-r", str(bam_nomito), str(bam_dedup)],
            ["samtools", "index", str(bam_dedup)],
            ["samtools", "view", "-b", "-q", "30", str(bam_dedup), "-o", str(bam_final)],
            ["samtools", "index", str(bam_final)],
        ]:
            if runner.run(cmd) != 0:
                _die(f"Filtering failed at: {' '.join(str(c) for c in cmd[:3])}")
        _check(*checkpoints.check_atac_filtered_bam(bam_final))
        bam_for[name] = bam_final

        if group == "input":
            console.print(f"  [dim]{name}: input/control — aligned + filtered, not peak-called.[/dim]")
            continue

        # 4 — Peaks: ChIP-correct MACS2 (model-based narrow, or --broad; optional -c input).
        # NOT ATAC's --nomodel/--shift/--extsize, which model the Tn5 cut site.
        runner.step_header("Peak calling — MACS2", 4, TOTAL)
        if explain:
            _explain("chipseq_peaks.md")
        peaks_dir = sdir / "peaks"
        peaks_dir.mkdir(exist_ok=True)
        ctl = (s.get("control") or "").strip()
        control_bam = bam_for.get(ctl) if ctl else None
        macs2_cmd = ["macs2", "callpeak", "-t", str(bam_final)]
        if control_bam:
            macs2_cmd += ["-c", str(control_bam)]
        macs2_cmd += ["-f", "BAMPE" if paired else "BAM", "-g", "hs",
                      "--outdir", str(peaks_dir), "-n", name]
        if broad:
            macs2_cmd += ["--broad", "--broad-cutoff", "0.1"]
        if runner.run(macs2_cmd) != 0:
            _die(f"MACS2 peak calling failed for sample '{name}'.")
        peak_file = peaks_dir / f"{name}_peaks.{'broadPeak' if broad else 'narrowPeak'}"
        _check(*checkpoints.check_chip_peaks(peak_file, broad))
        chip_results.append((name, group, peak_file, bam_final))

    # 5 — Build output matrices
    runner.step_header("Build matrix output", 5, TOTAL)
    written = _build_peak_outputs(chip_results, outdir, output_format, threads,
                                  explain, "chipseq_export_formats.md")
    console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── Methylation ────────────────────────────────────────────────────────────


@app.command()
def methylation(
    outdir: OutdirOpt,
    method: Annotated[str, typer.Option("--method",
                      help="wgbs (Bismark pipeline, default) or array (Illumina 450K/EPIC from GEO).")] = "wgbs",
    samples: Annotated[Optional[Path], typer.Option("--samples",
                       help="Samplesheet CSV (for --method wgbs).")] = None,
    bismark_genome: Annotated[Optional[Path], typer.Option("--bismark-genome",
                              help="Bismark genome directory (for --method wgbs).")] = None,
    betas: Annotated[Optional[Path], typer.Option("--betas",
                    help="Beta matrix CSV from GEO (for --method array).")] = None,
    metadata: Annotated[Optional[Path], typer.Option("--metadata",
                        help="Metadata CSV with geo_accession and disease.state (for --method array).")] = None,
    output_format: Annotated[str, typer.Option(
        "--format",
        help="Output: 'obama' (default), 'matrix' (mvalues_matrix.csv + coldata.csv for "
             "limma), or 'both'. M-values = log2(beta/(1-beta)).",
    )] = "obama",
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """Methylation: WGBS (Bismark pipeline) or Illumina array (450K/EPIC beta matrix from GEO).

    Output is the OBAMA matrix by default; --format matrix (or both) also writes an
    M-value matrix (features × samples) + coldata.csv for limma. limma models
    M-values, not beta/percent values.
    """
    if method not in ("wgbs", "array"):
        _die(f"Unknown --method: {method!r}. Use 'wgbs' or 'array'.")
    if output_format not in ("obama", "matrix", "both"):
        _die("--format must be 'obama', 'matrix', or 'both'.")

    outdir.mkdir(parents=True, exist_ok=True)

    if method == "array":
        if not betas:
            _die("--betas is required for --method array.")
        if not metadata:
            _die("--metadata is required for --method array.")
        if not betas.exists():
            _die(f"Beta matrix not found: {betas}")
        if not metadata.exists():
            _die(f"Metadata file not found: {metadata}")

        TOTAL = 1
        console.print(f"\n[bold]Methylation (array) pipeline[/bold] — building matrices\n")
        runner.step_header("Merge beta values with metadata", 1, TOTAL)
        if explain:
            _explain("methylation_array.md")
            if output_format in ("matrix", "both"):
                _explain("methylation_export_formats.md")

        try:
            written = obama.write_methylation_array_outputs(betas, metadata, outdir, output_format)
        except ValueError as e:
            _die(str(e))
        if output_format in ("obama", "both"):
            _check(*checkpoints.check_obama_format(outdir / "obama_matrix.csv"))
        if output_format in ("matrix", "both"):
            _check(*checkpoints.check_counts_matrix(outdir / "mvalues_matrix.csv", outdir / "coldata.csv"))
        console.print(f"\n[bold green]Done.[/bold green] Wrote: "
                      f"{', '.join(p.name for p in written)} → {outdir}")
        return

    # ── WGBS path ──
    if not samples:
        _die("--samples is required for --method wgbs.")
    if not bismark_genome:
        _die("--bismark-genome is required for --method wgbs.")
    sample_list = _read_samplesheet(samples)
    _require_dir(bismark_genome, "Bismark genome directory")

    TOTAL = 4
    console.print(f"\n[bold]Methylation pipeline[/bold] — {len(sample_list)} sample(s) → {outdir}\n")

    cx_reports: list[tuple[str, str, Path]] = []

    for s in sample_list:
        name, group = s["name"], s["group"]
        paired = _is_paired(s)
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)

        layout = "paired-end" if paired else "single-end"
        console.rule(f"[bold white]Sample: {name}  ({group})  [{layout}][/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("methylation_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        cmd, t_r1, t_r2 = _fastp_cmd(
            s, name, trim_dir, threads,
            ["--detect_adapter_for_pe", "--trim_poly_g", "--length_required", "36"],
        )
        rc = runner.run(cmd)
        if rc != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        # 2 — Bismark alignment
        runner.step_header("Align — Bismark", 2, TOTAL)
        if explain:
            _explain("methylation_align.md")
        bismark_dir = sdir / "bismark"
        bismark_dir.mkdir(exist_ok=True)
        # Bismark: paired-end uses -1/-2 and emits *_pe.bam / *_PE_report.txt;
        # single-end takes the read file directly and emits *.bam / *_SE_report.txt.
        reads = (["-1", str(t_r1), "-2", str(t_r2)] if paired else [str(t_r1)])
        rc = runner.run([
            "bismark",
            "--genome",     str(bismark_genome),
            *reads,
            "--output_dir", str(bismark_dir),
            "-p",           "2",
            "--basename",   name,
        ])
        if rc != 0:
            _die(f"Bismark alignment failed for sample '{name}'.")
        # Resolve Bismark's output by globbing rather than hardcoding: the exact
        # name depends on the Bismark version and whether --basename strips the
        # "_bismark_bt2" tag. PE emits *_pe.bam/*_PE_report.txt; SE *.bam/*_SE_report.txt.
        bam_glob = "*_pe.bam" if paired else "*.bam"
        rep_glob = "*_PE_report.txt" if paired else "*_SE_report.txt"
        bams = sorted(bismark_dir.glob(f"{name}{bam_glob}")) or sorted(bismark_dir.glob(bam_glob))
        reps = sorted(bismark_dir.glob(f"{name}*{rep_glob[1:]}")) or sorted(bismark_dir.glob(rep_glob))
        if not bams:
            _die(f"Bismark BAM not found in {bismark_dir} for sample '{name}'.")
        bam = bams[0]
        report = reps[0] if reps else None
        _check(*checkpoints.check_bismark_bam(bam, report))

        # 3 — Methylation extraction
        runner.step_header("Methylation extraction — bismark_methylation_extractor", 3, TOTAL)
        if explain:
            _explain("methylation_extract.md")
        methyl_dir = sdir / "methylation"
        methyl_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "bismark_methylation_extractor",
            "--paired-end" if paired else "--single-end",
            "--comprehensive",
            "--CX_context",
            "--cytosine_report",
            "--genome_folder", str(bismark_genome),
            "--output",        str(methyl_dir),
            str(bam),
        ])
        if rc != 0:
            _die(f"bismark_methylation_extractor failed for sample '{name}'.")

        cx_files = list(methyl_dir.glob(f"{name}*.CX_report.txt"))
        if not cx_files:
            cx_files = list(methyl_dir.glob("*.CX_report.txt"))
        if not cx_files:
            _die(f"CX cytosine report not found in {methyl_dir} for sample '{name}'.")
        cx_report = cx_files[0]
        _check(*checkpoints.check_methylation_cx_report(cx_report))

        cx_reports.append((name, group, cx_report))

    # 4 — Build output matrices (OBAMA percent-methylation and/or limma M-values)
    runner.step_header("Build matrix output", 4, TOTAL)
    if explain and output_format in ("matrix", "both"):
        _explain("methylation_export_formats.md")
    written = obama.write_methylation_outputs(cx_reports, outdir, output_format)
    if output_format in ("obama", "both"):
        _check(*checkpoints.check_obama_format(outdir / "obama_matrix.csv"))
    if output_format in ("matrix", "both"):
        _check(*checkpoints.check_counts_matrix(outdir / "mvalues_matrix.csv", outdir / "coldata.csv"))

    console.print(f"\n[bold green]Done.[/bold green] Wrote: "
                  f"{', '.join(p.name for p in written)} → {outdir}")


# ── Download ───────────────────────────────────────────────────────────────


@app.command()
def download(
    track: Annotated[Optional[str], typer.Option(
        "--track",
        help="Track to set up: rnaseq, atacseq, chipseq, methylation, qc. Omit to list all.",
    )] = None,
    outdir: Annotated[Path, typer.Option(
        "--outdir",
        help="Directory to write samples.csv (and where to look for locally downloaded files).",
    )] = Path("data"),
    shared_data: Annotated[Path, typer.Option(
        "--shared-data",
        envvar="UPSTREAM_SHARED_DATA",
        help="Path to the server's shared FASTQ directory (default: /data/upstream/shared).",
    )] = Path("/data/upstream/shared"),
    instructions: Annotated[bool, typer.Option(
        "--instructions",
        help="Show manual download instructions instead of writing samples.csv.",
    )] = False,
) -> None:
    """Locate example data and write a samples.csv ready for the pipeline.

    Checks the shared server data directory first (set up by the course admin
    via scripts/prepare_sample_data.sh), then looks in --outdir for files the
    user downloaded themselves.  If neither is found, prints instructions for
    obtaining the data independently.

    Set UPSTREAM_SHARED_DATA to override the default shared-data path.
    """
    if track is None:
        downloader.print_catalog()
        return

    downloader.locate_or_download(
        track=track,
        outdir=outdir,
        shared_data=shared_data,
        instructions_only=instructions,
    )


# ── Serve ─────────────────────────────────────────────────────────────────


@app.command()
def serve(
    port: Annotated[int, typer.Option("--port", help="Port to listen on.")] = 8421,
    host: Annotated[str, typer.Option("--host", help="Host to bind (127.0.0.1 = localhost only).")] = "127.0.0.1",
) -> None:
    """Launch the web UI in a browser-accessible local server.

    Runs on localhost by default — access it at http://localhost:<port>.
    Over SSH, forward the port first:

      ssh -L <port>:localhost:<port> your_username@server

    Then open http://localhost:<port> in your local browser.
    """
    try:
        import uvicorn
    except ImportError:
        _die("uvicorn is not installed. Run: pip install uvicorn")

    from .server import app as web_app

    url = f"http://{host}:{port}"
    console.print(f"\n[bold]upstream web UI[/bold] → [bold cyan]{url}[/bold cyan]")
    if host == "127.0.0.1":
        console.print(
            f"[dim]SSH tunnel (if on a remote server):[/dim]\n"
            f"  [dim]ssh -L {port}:localhost:{port} your_username@server[/dim]\n"
        )
    uvicorn.run(web_app, host=host, port=port, log_level="warning")
