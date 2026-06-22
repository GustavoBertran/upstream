"""htsprep — HTS preprocessing pipeline

Usage:
  htsprep qc         --samples samples.csv --outdir results/
  htsprep rnaseq     --samples samples.csv --star-index /ref/star --salmon-index /ref/salmon --outdir results/
  htsprep atacseq    --samples samples.csv --bowtie2-index /ref/bt2/hg38 --outdir results/
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


def _read_samplesheet(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        _die(f"Samplesheet not found: {path}")
    with path.open() as f:
        reader = csv.DictReader(f)
        missing_cols = {"name", "group", "r1", "r2"} - set(reader.fieldnames or [])
        if missing_cols:
            _die(f"Samplesheet is missing columns: {', '.join(sorted(missing_cols))}")
        rows = list(reader)

    for row in rows:
        if row["group"] not in ("disease", "control"):
            _die(
                f"Invalid group '{row['group']}' for sample '{row['name']}'. "
                "Must be 'disease' or 'control' (OBAMA requires these exact strings)."
            )
        for key in ("r1", "r2"):
            if not Path(row[key]).exists():
                _die(f"File not found for sample '{row['name']}': {row[key]}")
    return rows


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
    all_fastq = [f for s in sample_list for f in (s["r1"], s["r2"])]
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
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """RNA-seq: trim (fastp) → align/quantify (STAR or Salmon) → OBAMA matrix.

    Use --aligner salmon (default) for alignment-free TPM quantification via Salmon.
    Use --aligner star for genome alignment with STAR; gene counts are produced via
    STAR's built-in --quantMode GeneCounts (no separate quantification tool needed).
    """
    if aligner not in ("salmon", "star"):
        _die("--aligner must be 'salmon' or 'star'.")
    if aligner == "salmon" and salmon_index is None:
        _die("--salmon-index is required when --aligner salmon.")
    if aligner == "star" and star_index is None:
        _die("--star-index is required when --aligner star.")

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
        r1, r2 = Path(s["r1"]), Path(s["r2"])
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)

        console.rule(f"[bold white]Sample: {name}  ({group})[/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("rnaseq_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "fastp",
            "--in1",  str(r1),
            "--in2",  str(r2),
            "--out1", str(trim_dir / f"{name}_R1.fastq.gz"),
            "--out2", str(trim_dir / f"{name}_R2.fastq.gz"),
            "--json", str(trim_dir / f"{name}_fastp.json"),
            "--html", str(trim_dir / f"{name}_fastp.html"),
            "--thread", str(threads),
            "--detect_adapter_for_pe",
        ])
        if rc != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        if aligner == "salmon":
            # 2 — Salmon quantification (alignment-free)
            runner.step_header("Quantify — Salmon", 2, TOTAL)
            if explain:
                _explain("rnaseq_quantify.md")
            quant_dir = sdir / "salmon"
            rc = runner.run([
                "salmon", "quant",
                "--index",    str(salmon_index),
                "--libType",  "A",
                "--mates1",   str(trim_dir / f"{name}_R1.fastq.gz"),
                "--mates2",   str(trim_dir / f"{name}_R2.fastq.gz"),
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
            rc = runner.run([
                "STAR",
                "--runThreadN",       str(threads),
                "--genomeDir",        str(star_index),
                "--readFilesIn",      str(trim_dir / f"{name}_R1.fastq.gz"),
                                      str(trim_dir / f"{name}_R2.fastq.gz"),
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

    # 3 — OBAMA matrix
    runner.step_header("OBAMA matrix", 3, TOTAL)
    matrix_path = outdir / "obama_matrix.csv"
    if aligner == "salmon":
        obama.build_rnaseq(result_dirs, matrix_path)
    else:
        obama.build_rnaseq_star(result_dirs, matrix_path)
    _check(*checkpoints.check_obama_format(matrix_path))

    console.print(f"\n[bold green]Done.[/bold green] OBAMA matrix → {matrix_path}")


# ── ATAC-seq ──────────────────────────────────────────────────────────────


@app.command()
def atacseq(
    samples: SamplesOpt,
    bowtie2_index: Annotated[Path, typer.Option("--bowtie2-index", help="Bowtie2 index prefix (path/to/hg38)")],
    outdir: OutdirOpt,
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """ATAC-seq: QC → trim (fastp) → align (Bowtie2) → filter → peaks (MACS2) → OBAMA matrix."""
    sample_list = _read_samplesheet(samples)
    bt2_prefix = bowtie2_index  # e.g. /ref/bowtie2/hg38 (no .bt2 extension)
    outdir.mkdir(parents=True, exist_ok=True)

    TOTAL = 5
    console.print(f"\n[bold]ATAC-seq pipeline[/bold] — {len(sample_list)} sample(s) → {outdir}\n")

    peak_files: list[tuple[str, str, Path]] = []

    for s in sample_list:
        name, group = s["name"], s["group"]
        r1, r2 = Path(s["r1"]), Path(s["r2"])
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)

        console.rule(f"[bold white]Sample: {name}  ({group})[/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("atacseq_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "fastp",
            "--in1",  str(r1), "--in2", str(r2),
            "--out1", str(trim_dir / f"{name}_R1.fastq.gz"),
            "--out2", str(trim_dir / f"{name}_R2.fastq.gz"),
            "--json", str(trim_dir / f"{name}_fastp.json"),
            "--html", str(trim_dir / f"{name}_fastp.html"),
            "--thread", str(threads),
            "--detect_adapter_for_pe",
        ])
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

        rc = runner.pipe(
            ["bowtie2",
             "--threads", str(threads),
             "-x", str(bt2_prefix),
             "-1", str(trim_dir / f"{name}_R1.fastq.gz"),
             "-2", str(trim_dir / f"{name}_R2.fastq.gz"),
             "--very-sensitive", "--no-mixed", "--no-discordant", "-X", "2000"],
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

        peak_files.append((name, group, peak_file))

    # 5 — OBAMA matrix
    runner.step_header("OBAMA matrix", 5, TOTAL)
    matrix_path = outdir / "obama_matrix.csv"
    obama.build_atacseq(peak_files, matrix_path)
    _check(*checkpoints.check_obama_format(matrix_path))

    console.print(f"\n[bold green]Done.[/bold green] OBAMA matrix → {matrix_path}")


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
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
) -> None:
    """Methylation: WGBS (Bismark pipeline) or Illumina array (450K/EPIC beta matrix from GEO)."""
    if method not in ("wgbs", "array"):
        _die(f"Unknown --method: {method!r}. Use 'wgbs' or 'array'.")

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
        console.print(f"\n[bold]Methylation (array) pipeline[/bold] — building OBAMA matrix\n")
        runner.step_header("OBAMA matrix — merge beta values with metadata", 1, TOTAL)
        if explain:
            _explain("methylation_array.md")

        matrix_path = outdir / "obama_matrix.csv"
        try:
            obama.build_methylation_array(betas, metadata, matrix_path)
        except ValueError as e:
            _die(str(e))
        _check(*checkpoints.check_obama_format(matrix_path))
        console.print(f"\n[bold green]Done.[/bold green] OBAMA matrix → {matrix_path}")
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
        r1, r2 = Path(s["r1"]), Path(s["r2"])
        sdir = outdir / name
        sdir.mkdir(exist_ok=True)

        console.rule(f"[bold white]Sample: {name}  ({group})[/bold white]", style="white")

        # 1 — Trim
        runner.step_header("Trim — fastp", 1, TOTAL)
        if explain:
            _explain("methylation_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "fastp",
            "--in1",  str(r1), "--in2", str(r2),
            "--out1", str(trim_dir / f"{name}_R1.fastq.gz"),
            "--out2", str(trim_dir / f"{name}_R2.fastq.gz"),
            "--json", str(trim_dir / f"{name}_fastp.json"),
            "--html", str(trim_dir / f"{name}_fastp.html"),
            "--thread", str(threads),
            "--detect_adapter_for_pe",
            "--trim_poly_g",
            "--length_required", "36",
        ])
        if rc != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        # 2 — Bismark alignment
        runner.step_header("Align — Bismark", 2, TOTAL)
        if explain:
            _explain("methylation_align.md")
        bismark_dir = sdir / "bismark"
        bismark_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "bismark",
            "--genome",     str(bismark_genome),
            "-1",           str(trim_dir / f"{name}_R1.fastq.gz"),
            "-2",           str(trim_dir / f"{name}_R2.fastq.gz"),
            "--output_dir", str(bismark_dir),
            "-p",           "2",
            "--basename",   name,
        ])
        if rc != 0:
            _die(f"Bismark alignment failed for sample '{name}'.")
        bam = bismark_dir / f"{name}_bismark_bt2_pe.bam"
        report = bismark_dir / f"{name}_bismark_bt2_PE_report.txt"
        _check(*checkpoints.check_bismark_bam(bam, report))

        # 3 — Methylation extraction
        runner.step_header("Methylation extraction — bismark_methylation_extractor", 3, TOTAL)
        if explain:
            _explain("methylation_extract.md")
        methyl_dir = sdir / "methylation"
        methyl_dir.mkdir(exist_ok=True)
        rc = runner.run([
            "bismark_methylation_extractor",
            "--paired-end",
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

    # 4 — OBAMA matrix
    runner.step_header("OBAMA matrix", 4, TOTAL)
    matrix_path = outdir / "obama_matrix.csv"
    obama.build_methylation(cx_reports, matrix_path)
    _check(*checkpoints.check_obama_format(matrix_path))

    console.print(f"\n[bold green]Done.[/bold green] OBAMA matrix → {matrix_path}")


# ── Download ───────────────────────────────────────────────────────────────


@app.command()
def download(
    track: Annotated[Optional[str], typer.Option(
        "--track",
        help="Track to set up: rnaseq, atacseq, methylation, qc. Omit to list all.",
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
