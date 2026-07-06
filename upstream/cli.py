"""htsprep — HTS preprocessing pipeline

Usage:
  htsprep qc          --samples samples.csv --outdir results/
  htsprep genomics    --samples samples.csv --reference /ref/hg38.fa --outdir results/
  htsprep rnaseq      --samples samples.csv --star-index /ref/star --salmon-index /ref/salmon --outdir results/
  htsprep atacseq     --samples samples.csv --bowtie2-index /ref/bt2/hg38 --outdir results/
  htsprep chipseq     --samples samples.csv --bowtie2-index /ref/bt2/hg38 --peak-type narrow --outdir results/
  htsprep methylation --samples samples.csv --bismark-genome /ref/bismark --outdir results/
  htsprep proteomics  --intensities proteinGroups.txt --metadata metadata.csv --outdir results/

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

from . import (
    __version__, checkpoints, diagnostics, downloader, indexhelp, obama,
    preflight, runner, variants,
)
from .samplesheet import scan_fastq_dir

app = typer.Typer(
    name="htsprep",
    help="HTS preprocessing pipeline — trim, align, quantify, and produce OBAMA-format matrices.",
    no_args_is_help=True,
    add_completion=False,
)
console = Console()


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"upstream {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: Annotated[bool, typer.Option(
        "--version", help="Show the upstream version and exit.",
        callback=_version_callback, is_eager=True,
    )] = False,
) -> None:
    """upstream — HTS preprocessing pipeline (trim, align, quantify → matrices)."""


# ── Types (shorthand for Annotated options) ──────────────────────────────

SamplesOpt = Annotated[Path, typer.Option("--samples", help="Samplesheet CSV (name,group,r1,r2)")]
OutdirOpt  = Annotated[Path, typer.Option("--outdir",  help="Output directory (created if missing)")]
ThreadsOpt = Annotated[int,  typer.Option("--threads", help="CPU threads for tools that support it")]
ExplainOpt = Annotated[bool, typer.Option("--explain/--no-explain",
                                           help="Show educational explanations (on by default)")]
DryRunOpt  = Annotated[bool, typer.Option("--dry-run",
                                          help="Print the commands each step would run, without executing.")]


# ── Helpers ───────────────────────────────────────────────────────────────


def _read_samplesheet(path: Path, allow_input: bool = False) -> list[dict[str, str]]:
    issues = preflight.validate_samplesheet(path, allow_input)
    if issues:
        _die("samplesheet problems found:\n  - " + "\n  - ".join(issues))
    with path.open() as f:
        return list(csv.DictReader(f))


def _require_tools(track: str, dry: bool = False, **opts) -> None:
    """Fail fast with a clear message if a track's tools aren't on PATH.

    Skipped under --dry-run, which previews commands without needing the tools.
    """
    if dry:
        return
    missing = preflight.missing_tools(preflight.required_tools(track, **opts))
    if missing:
        _die(
            f"required tool(s) not found on PATH: {', '.join(missing)}. "
            "Activate the environment first (conda activate upstream)."
        )


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
    if runner.DRY_RUN:
        console.print(f"[yellow][dry-run][/yellow] would build '{output_format}' output(s) "
                      f"in {outdir} (consensus peaks + multiBamSummary counts / OBAMA matrix).")
        return []

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


def _require_index(path: Optional[Path], tool: str, label: str) -> None:
    """Like _require_dir, but a missing index points the user at `index-help`."""
    if path is None or not Path(path).exists():
        console.print(f"[bold red]Error:[/bold red] {label} not found: {path}")
        console.print(f"[dim]Don't have this index yet? See how to build it:[/dim] "
                      f"[bold]upstream index-help --tool {tool} --genome <human|mouse>[/bold]")
        raise typer.Exit(1)


def _require_bwa_index(reference: Optional[Path], need_dict: bool = False) -> None:
    """Require a reference FASTA with its bwa (.bwt) and samtools faidx (.fai) sidecars.

    A bare FASTA without these would pass an existence check then fail mid-run, so we
    validate the companions up front and point at index-help. With need_dict (GATK), a
    sequence dictionary (.dict) is required too.
    """
    hint = ("[dim]Build it:[/dim] "
            "[bold]upstream index-help --tool bwa --genome <human|mouse>[/bold]")
    if reference is None or not Path(reference).is_file():
        console.print(f"[bold red]Error:[/bold red] reference FASTA not found: {reference}")
        console.print(hint)
        raise typer.Exit(1)
    ref = Path(reference)
    if not Path(str(ref) + ".bwt").exists():
        console.print(f"[bold red]Error:[/bold red] bwa index missing for {ref.name} "
                      f"(expected {ref.name}.bwt).")
        console.print(hint)
        raise typer.Exit(1)
    if not Path(str(ref) + ".fai").exists():
        console.print(f"[bold red]Error:[/bold red] FASTA index missing for {ref.name} "
                      f"(expected {ref.name}.fai). Run: [bold]samtools faidx {ref}[/bold]")
        raise typer.Exit(1)
    if need_dict and not ref.with_suffix(".dict").exists():
        console.print(f"[bold red]Error:[/bold red] GATK sequence dictionary missing "
                      f"(expected {ref.with_suffix('.dict').name}). Run: "
                      f"[bold]samtools dict {ref} -o {ref.with_suffix('.dict')}[/bold]")
        raise typer.Exit(1)


def _check(ok: bool, msg: str) -> None:
    if runner.DRY_RUN:
        return  # nothing was produced to validate
    if ok:
        runner.ok(msg)
    else:
        runner.fail(msg)
        raise typer.Exit(1)


def _next_steps(track: str, output_format: Optional[str]) -> list[str]:
    """Human-readable 'what to do next' lines for the run summary."""
    if track == "qc":
        return ["Open multiqc_report.html in a browser to review per-sample QC."]
    if track == "genomics":
        return [
            "Review per-sample variants/<name>.vcf.gz and variant_summary.csv; "
            "load cohort.vcf.gz (if merged) into IGV or a variant browser.",
            "Production gold standard: GATK HaplotypeCaller (GVCF) + GenotypeGVCFs.",
        ]
    if track == "proteomics":
        return [
            "limma: load proteins_matrix.csv + coldata.csv "
            "(see content/proteomics_export_formats.md). Values are normalized log2 "
            "intensities; missing values (blank) are tolerated.",
        ]
    steps: list[str] = []
    fmt = output_format or "obama"
    if fmt in ("obama", "both"):
        steps.append("OBAMA: load obama_matrix.csv into the OBAMA Shiny app "
                     "(geo_accession + disease.state columns).")
    if fmt in ("matrix", "both"):
        if track == "methylation":
            steps.append("limma: load mvalues_matrix.csv + coldata.csv "
                         "(see content/methylation_export_formats.md).")
        elif track in ("atacseq", "chipseq"):
            steps.append(f"DESeq2/edgeR: load counts_matrix.csv + coldata.csv "
                         f"(see content/{track}_export_formats.md).")
        else:
            steps.append("DESeq2/edgeR/limma-voom: load counts_matrix.csv + coldata.csv "
                         "(see content/rnaseq_export_formats.md).")
    return steps


def _write_run_summary(
    outdir: Path,
    track: str,
    written: list,
    output_format: Optional[str] = None,
    inputs: Optional[dict] = None,
) -> None:
    """Write run_summary.txt (skipped in dry-run): inputs, outputs, and next steps."""
    if runner.DRY_RUN:
        return
    from datetime import datetime as _dt
    lines = [
        f"upstream {track} — run summary",
        f"generated:        {_dt.now().isoformat(timespec='seconds')}",
        f"output directory: {Path(outdir).resolve()}",
        "",
    ]
    if inputs:
        lines.append("inputs:")
        lines += [f"  {k}: {v}" for k, v in inputs.items() if v]
        lines.append("")
    lines.append("outputs:")
    lines += [f"  - {Path(p).name}" for p in written] or ["  (none)"]
    steps = _next_steps(track, output_format)
    if steps:
        lines.append("")
        lines.append("next steps:")
        lines += [f"  - {s}" for s in steps]
    (Path(outdir) / "run_summary.txt").write_text("\n".join(lines) + "\n")
    console.print(f"[dim]Run summary → {Path(outdir) / 'run_summary.txt'}[/dim]")


def _auto_multiqc(outdir: Path):
    """Best-effort: aggregate the run's logs (fastp/STAR/Salmon/Bismark/…) into one
    MultiQC report. Skipped in dry-run or if multiqc isn't installed (it's a bonus,
    not a hard requirement). Returns the report path, or None."""
    import shutil as _sh
    if runner.DRY_RUN or _sh.which("multiqc") is None:
        return None
    mqc_dir = Path(outdir) / "multiqc"
    console.print("  Aggregating QC with MultiQC…")
    rc = runner.run(["multiqc", str(outdir), "--outdir", str(mqc_dir), "--force", "--quiet"])
    report = mqc_dir / "multiqc_report.html"
    return report if (rc == 0 and report.exists()) else None


def _explain(content_file: str) -> None:
    try:
        pkg = importlib.resources.files("upstream") / "content" / content_file
        text = pkg.read_text(encoding="utf-8")
        console.print(Panel(Markdown(text), border_style="dim blue", padding=(1, 2)))
    except Exception:
        pass  # explanations are optional; missing file is not an error


# ── QC ────────────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Quality control")
def qc(
    samples: SamplesOpt,
    outdir: OutdirOpt,
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
    dry_run: DryRunOpt = False,
) -> None:
    """FastQC + MultiQC quality control on all samples in the samplesheet."""
    sample_list = _read_samplesheet(samples)
    _require_tools("qc", dry=dry_run)
    runner.DRY_RUN = dry_run
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

    _write_run_summary(outdir, "qc",
                       [multiqc_out / "multiqc_report.html", fastqc_out], None,
                       {"samples": str(samples)})
    console.print(f"\n[bold green]QC complete.[/bold green] Open {multiqc_out}/multiqc_report.html")


# ── Genomics ──────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Genomics")
def genomics(
    samples: SamplesOpt,
    reference: Annotated[Path, typer.Option(
        "--reference",
        help="Reference genome FASTA. Its bwa index (.bwt/.amb/.ann/.pac/.sa) and "
             ".fai must sit alongside it — see index-help --tool bwa.",
    )],
    outdir: OutdirOpt,
    caller: Annotated[str, typer.Option(
        "--caller",
        help="Variant caller: 'bcftools' (default, lightweight) or 'gatk' "
             "(GATK HaplotypeCaller — the field standard; needs a .dict).",
    )] = "bcftools",
    merge: Annotated[bool, typer.Option(
        "--merge/--no-merge",
        help="After per-sample calling, merge per-sample VCFs into cohort.vcf.gz "
             "(bcftools merge). Auto-skipped with a single sample.",
    )] = True,
    min_mapq: Annotated[int, typer.Option(
        "--min-mapq", help="Minimum read mapping quality (bcftools mpileup -q).")] = 20,
    min_baseq: Annotated[int, typer.Option(
        "--min-baseq", help="Minimum base quality (bcftools mpileup -Q).")] = 20,
    threads: ThreadsOpt = 4,
    explain: ExplainOpt = True,
    dry_run: DryRunOpt = False,
) -> None:
    """Genomics: trim (fastp) → align (BWA-MEM) → mark duplicates → call variants.

    Germline short-variant calling. Each sample is aligned with BWA-MEM, duplicates
    are marked, and SNVs/indels are called per sample — with bcftools mpileup|call
    (default) or GATK HaplotypeCaller (--caller gatk, the field standard). With
    --merge (default) the per-sample VCFs are combined into cohort.vcf.gz. Output is
    VCF + variant_summary.csv — no matrix, no OBAMA. Joint genotyping (GATK GVCF →
    GenotypeGVCFs) is the further production step, described in the step explanations.
    """
    if caller not in ("bcftools", "gatk"):
        _die("--caller must be 'bcftools' or 'gatk'.")
    sample_list = _read_samplesheet(samples)
    _require_tools("genomics", dry=dry_run, caller=caller)
    if not dry_run:
        _require_bwa_index(reference, need_dict=(caller == "gatk"))
    runner.DRY_RUN = dry_run
    outdir.mkdir(parents=True, exist_ok=True)

    TOTAL = 5
    console.print(f"\n[bold]Genomics pipeline[/bold] ({caller}) — "
                  f"{len(sample_list)} sample(s) → {outdir}\n")

    results: list[tuple[str, str, Path, Path]] = []   # (name, group, vcf, stats)

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
            _explain("genomics_trim.md")
        trim_dir = sdir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        cmd, t_r1, t_r2 = _fastp_cmd(s, name, trim_dir, threads, ["--detect_adapter_for_pe"])
        if runner.run(cmd) != 0:
            _die(f"fastp failed for sample '{name}'.")
        _check(*checkpoints.check_fastp_trim(trim_dir / f"{name}_fastp.json"))

        # 2 — Align (BWA-MEM | samtools sort)
        runner.step_header("Align — BWA-MEM", 2, TOTAL)
        if explain:
            _explain("genomics_align.md")
        aln_dir = sdir / "aligned"
        aln_dir.mkdir(exist_ok=True)
        # bwa wants a literal '\t'-delimited @RG string; SM:<name> drives the VCF sample column.
        rg = f"@RG\\tID:{name}\\tSM:{name}\\tPL:ILLUMINA\\tLB:{name}"
        reads = [str(t_r1), str(t_r2)] if paired else [str(t_r1)]
        bam_sorted = aln_dir / f"{name}.sorted.bam"
        rc = runner.pipe(
            ["bwa", "mem", "-t", str(threads), "-R", rg, str(reference), *reads],
            ["samtools", "sort", "-@", str(threads), "-o", str(bam_sorted), "-"],
        )
        if rc != 0:
            _die(f"bwa mem/samtools sort failed for sample '{name}'.")
        runner.run(["samtools", "index", str(bam_sorted)])
        _check(*checkpoints.check_bwa_bam(bam_sorted))

        # 3 — Mark duplicates (germline FLAGS dups, does not remove them)
        runner.step_header("Mark duplicates — samtools markdup", 3, TOTAL)
        if explain:
            _explain("genomics_markdup.md")
        markdup_bam = aln_dir / f"{name}.markdup.bam"
        if paired:
            # markdup needs the ms tags from `fixmate -m`; chain collate→fixmate→sort→markdup.
            collate_bam = aln_dir / f"{name}.collate.bam"
            fixmate_bam = aln_dir / f"{name}.fixmate.bam"
            possort_bam = aln_dir / f"{name}.possort.bam"
            cmds = [
                ["samtools", "collate", "-@", str(threads), "-o", str(collate_bam), str(bam_sorted)],
                ["samtools", "fixmate", "-m", str(collate_bam), str(fixmate_bam)],
                ["samtools", "sort", "-@", str(threads), "-o", str(possort_bam), str(fixmate_bam)],
                ["samtools", "markdup", str(possort_bam), str(markdup_bam)],
                ["samtools", "index", str(markdup_bam)],
            ]
        else:
            # single-end: no mate, so markdup runs directly on the coord-sorted BAM.
            cmds = [
                ["samtools", "markdup", str(bam_sorted), str(markdup_bam)],
                ["samtools", "index", str(markdup_bam)],
            ]
        for cmd in cmds:
            if runner.run(cmd) != 0:
                _die(f"markdup chain failed at: {' '.join(str(c) for c in cmd[:2])}")
        _check(*checkpoints.check_markdup_bam(markdup_bam))

        # 4 — Call variants (bcftools mpileup|call, or GATK HaplotypeCaller), bgzipped per sample
        runner.step_header(
            f"Call variants — {'GATK HaplotypeCaller' if caller == 'gatk' else 'bcftools'}",
            4, TOTAL,
        )
        if explain:
            _explain("genomics_call.md")
        var_dir = sdir / "variants"
        var_dir.mkdir(exist_ok=True)
        vcf = var_dir / f"{name}.vcf.gz"
        if caller == "gatk":
            # HaplotypeCaller: local reassembly caller; needs the ref .fai + .dict, a
            # coord-sorted, duplicate-marked, indexed BAM with read groups (all present).
            rc = runner.run([
                "gatk", "HaplotypeCaller",
                "-R", str(reference),
                "-I", str(markdup_bam),
                "-O", str(vcf),
            ])
            if rc != 0:
                _die(f"GATK HaplotypeCaller failed for sample '{name}'.")
            runner.run(["bcftools", "index", "-f", "-t", str(vcf)])  # ensure a tabix index for merge
        else:
            rc = runner.pipe(
                ["bcftools", "mpileup", "-f", str(reference),
                 "-q", str(min_mapq), "-Q", str(min_baseq),
                 "-a", "FORMAT/AD,FORMAT/DP", "-Ou", str(markdup_bam)],
                ["bcftools", "call", "-mv", "-Oz", "-o", str(vcf)],
            )
            if rc != 0:
                _die(f"bcftools mpileup/call failed for sample '{name}'.")
            runner.run(["bcftools", "index", "-t", str(vcf)])  # tabix index (merge needs it)
        _check(*checkpoints.check_vcf(vcf))
        stats = var_dir / f"{name}.stats.txt"
        runner.run_capture(["bcftools", "stats", str(vcf)], stats)
        results.append((name, group, vcf, stats))

    # 5 — Merge (optional) + combined summary
    runner.step_header("Variant summary — bcftools stats", 5, TOTAL)
    if explain:
        _explain("genomics_summary.md")
    if dry_run:
        console.print(f"[yellow][dry-run][/yellow] would merge {len(results)} VCF(s) (if --merge) "
                      f"and write variant_summary.csv → {outdir}")
        console.print("\n[bold green]Dry run complete.[/bold green] No tools were executed.")
        return

    written: list[Path] = [v for _n, _g, v, _s in results]
    if merge and len(results) >= 2:
        cohort = outdir / "cohort.vcf.gz"
        rc = runner.run(["bcftools", "merge", "-Oz", "-o", str(cohort),
                         "--threads", str(threads), *[str(v) for _n, _g, v, _s in results]])
        if rc != 0:
            _die("bcftools merge failed.")
        runner.run(["bcftools", "index", "-t", str(cohort)])
        _check(*checkpoints.check_vcf(cohort))
        runner.run_capture(["bcftools", "stats", str(cohort)], outdir / "cohort.stats.txt")
        written.append(cohort)
    elif merge:
        console.print("  [dim]Single sample — skipping cohort merge (nothing to merge).[/dim]")

    summary = outdir / "variant_summary.csv"
    variants.write_variant_summary(results, summary)
    _check(*checkpoints.check_variant_summary(summary))
    written.append(summary)

    mqc = _auto_multiqc(outdir)
    if mqc:
        written = list(written) + [mqc]
    _write_run_summary(outdir, "genomics", written, None,
                       {"samples": str(samples), "reference": str(reference),
                        "caller": caller, "merge": str(merge)})
    console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── RNA-seq ───────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Transcriptomics")
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
    dry_run: DryRunOpt = False,
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
    _require_tools("rnaseq", dry=dry_run, aligner=aligner, output_format=output_format)
    if aligner == "salmon":
        _require_index(salmon_index, "salmon", "Salmon index")
    else:
        _require_index(star_index, "star", "STAR index")
    runner.DRY_RUN = dry_run
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
    if dry_run:
        console.print(f"[yellow][dry-run][/yellow] would build '{output_format}' output(s) "
                      f"from {len(result_dirs)} sample(s) → {outdir}")
        console.print("\n[bold green]Dry run complete.[/bold green] No tools were executed.")
        return
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

    mqc = _auto_multiqc(outdir)
    if mqc:
        written = list(written) + [mqc]
    _write_run_summary(outdir, "rnaseq", written, output_format,
                       {"samples": str(samples), "aligner": aligner})
    if dry_run:
        console.print("\n[bold green]Dry run complete.[/bold green] No tools were executed.")
    else:
        console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── ATAC-seq ──────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Epigenomics")
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
    dry_run: DryRunOpt = False,
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
    _require_tools("atacseq", dry=dry_run, output_format=output_format)
    bt2_prefix = bowtie2_index  # e.g. /ref/bowtie2/hg38 (no .bt2 extension)
    runner.DRY_RUN = dry_run
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
    mqc = _auto_multiqc(outdir)
    if mqc:
        written = list(written) + [mqc]
    _write_run_summary(outdir, "atacseq", written, output_format, {"samples": str(samples)})
    if dry_run:
        console.print("\n[bold green]Dry run complete.[/bold green] No tools were executed.")
    else:
        console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── ChIP-seq ─────────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Epigenomics")
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
    dry_run: DryRunOpt = False,
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

    _require_tools("chipseq", dry=dry_run, output_format=output_format)
    bt2_prefix = bowtie2_index
    runner.DRY_RUN = dry_run
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
    mqc = _auto_multiqc(outdir)
    if mqc:
        written = list(written) + [mqc]
    _write_run_summary(outdir, "chipseq", written, output_format,
                       {"samples": str(samples), "peak_type": peak_type})
    if dry_run:
        console.print("\n[bold green]Dry run complete.[/bold green] No tools were executed.")
    else:
        console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── Methylation ────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Epigenomics")
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
    dry_run: DryRunOpt = False,
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

    runner.DRY_RUN = dry_run
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

        if dry_run:
            console.print(f"[yellow][dry-run][/yellow] would merge betas + metadata into "
                          f"'{output_format}' output(s) → {outdir} (no external tools run).")
            console.print("\n[bold green]Dry run complete.[/bold green]")
            return
        try:
            written = obama.write_methylation_array_outputs(betas, metadata, outdir, output_format)
        except ValueError as e:
            _die(str(e))
        if output_format in ("obama", "both"):
            _check(*checkpoints.check_obama_format(outdir / "obama_matrix.csv"))
        if output_format in ("matrix", "both"):
            _check(*checkpoints.check_counts_matrix(outdir / "mvalues_matrix.csv", outdir / "coldata.csv"))
        _write_run_summary(outdir, "methylation", written, output_format,
                           {"method": "array", "betas": str(betas), "metadata": str(metadata)})
        console.print(f"\n[bold green]Done.[/bold green] Wrote: "
                      f"{', '.join(p.name for p in written)} → {outdir}")
        return

    # ── WGBS path ──
    if not samples:
        _die("--samples is required for --method wgbs.")
    if not bismark_genome:
        _die("--bismark-genome is required for --method wgbs.")
    sample_list = _read_samplesheet(samples)
    _require_tools("methylation", dry=dry_run, method="wgbs")
    _require_index(bismark_genome, "bismark", "Bismark genome directory")

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
        if dry_run:
            bam = bismark_dir / f"{name}_bismark_bt2{'_pe' if paired else ''}.bam"  # expected (not resolved)
        else:
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

        if dry_run:
            cx_report = methyl_dir / f"{name}.CX_report.txt"  # expected (not resolved)
        else:
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
    if dry_run:
        console.print(f"[yellow][dry-run][/yellow] would build '{output_format}' output(s) "
                      f"from {len(cx_reports)} sample(s) → {outdir}")
        console.print("\n[bold green]Dry run complete.[/bold green] No tools were executed.")
        return
    if explain and output_format in ("matrix", "both"):
        _explain("methylation_export_formats.md")
    written = obama.write_methylation_outputs(cx_reports, outdir, output_format)
    if output_format in ("obama", "both"):
        _check(*checkpoints.check_obama_format(outdir / "obama_matrix.csv"))
    if output_format in ("matrix", "both"):
        _check(*checkpoints.check_counts_matrix(outdir / "mvalues_matrix.csv", outdir / "coldata.csv"))

    mqc = _auto_multiqc(outdir)
    if mqc:
        written = list(written) + [mqc]
    _write_run_summary(outdir, "methylation", written, output_format,
                       {"method": "wgbs", "samples": str(samples)})
    console.print(f"\n[bold green]Done.[/bold green] Wrote: "
                  f"{', '.join(p.name for p in written)} → {outdir}")


# ── Proteomics ───────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Proteomics")
def proteomics(
    intensities: Annotated[Path, typer.Option(
        "--intensities",
        help="MaxQuant proteinGroups.txt OR a generic protein x sample matrix (CSV/TSV).")],
    metadata: Annotated[Path, typer.Option(
        "--metadata",
        help="Metadata CSV with columns 'sample' and 'group' (group = disease/control).")],
    outdir: OutdirOpt,
    input_type: Annotated[str, typer.Option(
        "--input-type",
        help="'maxquant' (proteinGroups.txt) or 'matrix' (generic protein x sample table).")] = "maxquant",
    intensity_col: Annotated[str, typer.Option(
        "--intensity-col",
        help="MaxQuant quant column: 'LFQ' (default), 'iBAQ', or 'Intensity'.")] = "LFQ",
    min_valid: Annotated[float, typer.Option(
        "--min-valid",
        help="Min fraction (0..1) of non-missing values required PER GROUP to keep a protein.")] = 0.5,
    normalize: Annotated[str, typer.Option(
        "--normalize",
        help="Normalization: 'median' (default), 'quantile', or 'none'.")] = "median",
    explain: ExplainOpt = True,
    dry_run: DryRunOpt = False,
) -> None:
    """Proteomics: intensity matrix → filter → log2 → normalize → limma matrix + coldata.

    Starts from a MaxQuant proteinGroups.txt (or a generic protein x sample matrix),
    drops contaminant/reverse hits, filters proteins by a per-group valid-value
    fraction, log2-transforms, normalizes, and writes proteins_matrix.csv + coldata.csv
    for limma. Pure Python — no external tools, no OBAMA. The raw-spectra → matrix step
    (MaxQuant/FragPipe) runs upstream and is out of scope.
    """
    if input_type not in ("maxquant", "matrix"):
        _die("--input-type must be 'maxquant' or 'matrix'.")
    if intensity_col not in ("LFQ", "iBAQ", "Intensity"):
        _die("--intensity-col must be 'LFQ', 'iBAQ', or 'Intensity'.")
    if normalize not in ("median", "quantile", "none"):
        _die("--normalize must be 'median', 'quantile', or 'none'.")
    if not intensities.exists():
        _die(f"Intensity table not found: {intensities}")
    if not metadata.exists():
        _die(f"Metadata file not found: {metadata}")

    runner.DRY_RUN = dry_run
    outdir.mkdir(parents=True, exist_ok=True)

    TOTAL = 1
    console.print("\n[bold]Proteomics pipeline[/bold] — building matrices\n")
    runner.step_header("Build protein expression matrix", 1, TOTAL)
    if explain:
        _explain("proteomics_quantify.md")
        _explain("proteomics_export_formats.md")

    if dry_run:
        console.print(f"[yellow][dry-run][/yellow] would parse {input_type} intensities, drop "
                      f"contaminant/reverse rows, valid-value filter (≥{min_valid}/group), log2, "
                      f"normalize ('{normalize}') → proteins_matrix.csv + coldata.csv → {outdir} "
                      "(no external tools).")
        console.print("\n[bold green]Dry run complete.[/bold green]")
        return

    try:
        written = obama.write_proteomics_outputs(
            intensities, metadata, outdir,
            input_type=input_type, intensity_col=intensity_col,
            min_valid=min_valid, normalize=normalize,
        )
    except ValueError as e:
        _die(str(e))
    _check(*checkpoints.check_counts_matrix(outdir / "proteins_matrix.csv", outdir / "coldata.csv"))
    _write_run_summary(outdir, "proteomics", written, None,
                       {"intensities": str(intensities), "metadata": str(metadata),
                        "input_type": input_type, "intensity_col": intensity_col,
                        "min_valid": str(min_valid), "normalize": normalize})
    console.print(f"\n[bold green]Done.[/bold green] Wrote: {', '.join(p.name for p in written)} → {outdir}")


# ── Samplesheet generation ───────────────────────────────────────────────────


@app.command(rich_help_panel="Utilities")
def samplesheet(
    directory: Annotated[Path, typer.Option("--dir", help="Folder of FASTQ files to scan.")],
    out: Annotated[Optional[Path], typer.Option(
        "--out", help="Where to write samples.csv (default: <dir>/samples.csv).")] = None,
) -> None:
    """Generate a starter samples.csv from a folder of FASTQ files.

    Pairs R1/R2 mates, infers paired vs single-end, and derives sample names. The
    'group' column is left blank — fill it with disease/control before running.
    """
    rows, warnings = scan_fastq_dir(directory)
    for w in warnings:
        console.print(f"[yellow]•[/yellow] {w}")
    if not rows:
        _die(f"No samples detected in {directory}.")

    out_path = out or (Path(directory) / "samples.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "group", "r1", "r2"])
        writer.writeheader()
        writer.writerows(rows)

    n_pe = sum(1 for r in rows if r["r2"])
    console.print(
        f"\n[bold green]✓ Wrote {out_path}[/bold green] — {len(rows)} sample(s) "
        f"({n_pe} paired-end, {len(rows) - n_pe} single-end)."
    )
    console.print("[dim]Next: fill in the 'group' column (disease/control), then run a track "
                  "(or 'upstream check').[/dim]")


# ── Index help ───────────────────────────────────────────────────────────────


@app.command(name="index-help", rich_help_panel="Utilities")
def index_help(
    tool: Annotated[str, typer.Option("--tool", help="salmon | star | bowtie2 | bismark | bwa")],
    genome: Annotated[str, typer.Option("--genome", help="human | mouse")] = "human",
) -> None:
    """Print recommended commands to obtain/build a reference index.

    This only prints guidance — it does not download or build anything (index builds
    are large and machine-specific). Copy the commands and run them where you have
    the RAM/disk (STAR needs ~30 GB RAM, Bismark ~100 GB disk for human).
    """
    if tool not in indexhelp.TOOLS:
        _die(f"--tool must be one of: {', '.join(indexhelp.TOOLS)}")
    if genome not in indexhelp.GENOMES:
        _die(f"--genome must be one of: {', '.join(indexhelp.GENOMES)}")
    snippet = indexhelp.recommend(tool, genome)
    console.print(Panel(Markdown("```bash\n" + snippet + "\n```"),
                        title=f"[bold]Build a {tool} index — {indexhelp.GENOMES[genome]['label']}[/bold]",
                        border_style="cyan", padding=(1, 2)))
    console.print("[dim]GENCODE release shown is a sensible default — bump it for a newer one.[/dim]")


# ── Bug report ───────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Utilities")
def bug(
    out: Annotated[Optional[Path], typer.Option(
        "--out", help="Also save the diagnostics to this file.")] = None,
) -> None:
    """Print environment diagnostics and where to report a bug.

    Paste the diagnostics into a new issue along with the command you ran and the
    error message you saw.
    """
    text = diagnostics.as_text()
    console.print(Panel(text, title="[bold]Environment diagnostics[/bold]", border_style="yellow"))
    console.print(f"\nReport issues at: [bold cyan]{diagnostics.ISSUES_NEW_URL}[/bold cyan]")
    console.print("Include: (1) the diagnostics above, (2) the command you ran, "
                  "(3) the full error message.")
    if out:
        out.write_text(text + "\n")
        console.print(f"[dim]Saved diagnostics → {out}[/dim]")


# ── Check (preflight) ────────────────────────────────────────────────────────


@app.command(rich_help_panel="Utilities")
def check(
    track: Annotated[str, typer.Option(
        "--track", help="Track to validate: qc, genomics, rnaseq, atacseq, chipseq, "
                        "methylation, proteomics.")],
    samples: Annotated[Optional[Path], typer.Option("--samples", help="Samplesheet to validate.")] = None,
    aligner: Annotated[str, typer.Option("--aligner", help="rnaseq aligner (salmon/star).")] = "salmon",
    method: Annotated[str, typer.Option("--method", help="methylation method (wgbs/array).")] = "wgbs",
    output_format: Annotated[str, typer.Option("--format", help="Output format you plan to use.")] = "obama",
    salmon_index: Annotated[Optional[Path], typer.Option("--salmon-index")] = None,
    star_index: Annotated[Optional[Path], typer.Option("--star-index")] = None,
    bowtie2_index: Annotated[Optional[Path], typer.Option("--bowtie2-index")] = None,
    bismark_genome: Annotated[Optional[Path], typer.Option("--bismark-genome")] = None,
    reference: Annotated[Optional[Path], typer.Option("--reference", help="genomics reference FASTA.")] = None,
    caller: Annotated[str, typer.Option("--caller", help="genomics variant caller (bcftools/gatk).")] = "bcftools",
    intensities: Annotated[Optional[Path], typer.Option("--intensities", help="proteomics intensity table.")] = None,
    betas: Annotated[Optional[Path], typer.Option("--betas")] = None,
    metadata: Annotated[Optional[Path], typer.Option("--metadata")] = None,
) -> None:
    """Validate a run's inputs WITHOUT running anything.

    Reports every problem at once — bad samplesheet rows, missing FASTQ files,
    invalid group labels, unmatched ChIP controls, tools not on PATH, and missing
    index/genome/reference files — so you can fix them before a long run starts.
    """
    issues = preflight.run_issues(
        track,
        samples=str(samples) if samples else None,
        aligner=aligner, method=method, output_format=output_format,
        salmon_index=str(salmon_index) if salmon_index else None,
        star_index=str(star_index) if star_index else None,
        bismark_genome=str(bismark_genome) if bismark_genome else None,
        reference=str(reference) if reference else None,
        caller=caller,
        intensities=str(intensities) if intensities else None,
        betas=str(betas) if betas else None,
        metadata=str(metadata) if metadata else None,
    )
    if issues:
        console.print(f"\n[bold red]✗ {len(issues)} problem(s) found:[/bold red]")
        for x in issues:
            console.print(f"  [red]•[/red] {x}")
        console.print("\n[dim]Fix these and re-run [bold]upstream check[/bold], or run the track directly.[/dim]")
        raise typer.Exit(1)
    console.print(f"\n[bold green]✓ All checks passed[/bold green] — '{track}' is ready to run.")


# ── Download ───────────────────────────────────────────────────────────────


@app.command(rich_help_panel="Utilities")
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


@app.command(rich_help_panel="Utilities")
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


if __name__ == "__main__":   # enables `python -m upstream.cli` (server _upstream_argv fallback)
    app()
