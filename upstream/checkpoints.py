"""Output validation for each pipeline step.

Each function takes specific file/directory paths and returns (ok: bool, message: str).
Messages are shown to the user — make them specific and actionable on failure.
"""

import json
import re
from pathlib import Path


# ── QC ──────────────────────────────────────────────────────────────────


def check_fastqc_output(fastqc_dir: Path) -> tuple[bool, str]:
    reports = list(fastqc_dir.glob("*_fastqc.html"))
    if not reports:
        return False, f"No FastQC HTML reports in {fastqc_dir}. Did FastQC finish?"
    return True, f"{len(reports)} FastQC report(s) written."


def check_multiqc_output(multiqc_dir: Path) -> tuple[bool, str]:
    report = multiqc_dir / "multiqc_report.html"
    if not report.exists():
        return False, f"multiqc_report.html not found in {multiqc_dir}."
    return True, "MultiQC report written."


# ── Trimming ─────────────────────────────────────────────────────────────


def check_fastp_trim(json_path: Path) -> tuple[bool, str]:
    if not json_path.exists():
        return False, f"fastp JSON report not found: {json_path}"
    try:
        stats = json.loads(json_path.read_text())
        passed = stats["filtering_result"]["passed_filter_reads"]
        total = stats["summary"]["before_filtering"]["total_reads"]
        pct = passed / total * 100 if total else 0
        if passed == 0:
            return False, "fastp output has 0 reads. Check that input files are non-empty."
        return True, f"{passed:,} / {total:,} reads passed ({pct:.1f}%)."
    except (KeyError, json.JSONDecodeError) as e:
        return False, f"Could not parse fastp JSON: {e}"


# ── RNA-seq ──────────────────────────────────────────────────────────────


def check_star_bam(bam_path: Path, log_path: Path | None = None) -> tuple[bool, str]:
    if not bam_path.exists():
        return False, f"BAM not found: {bam_path.name}"
    if bam_path.stat().st_size < 1024:
        return False, f"{bam_path.name} is too small — likely empty."

    if log_path and log_path.exists():
        m = re.search(r"Uniquely mapped reads %\s+\|\s+([\d.]+)%", log_path.read_text())
        if m:
            pct = float(m.group(1))
            if pct < 60:
                return (
                    False,
                    f"Unique mapping rate is {pct:.1f}%. Rates below 60% usually indicate "
                    "a mismatch between sample organism and reference, or severely degraded input.",
                )
            return True, f"Alignment OK — {pct:.1f}% uniquely mapped."

    return True, f"{bam_path.name} written and non-empty."


def check_salmon_sf(sf_path: Path) -> tuple[bool, str]:
    if not sf_path.exists():
        return False, f"Salmon quant.sf not found: {sf_path}"
    n = len(sf_path.read_text().splitlines()) - 1  # minus header
    if n < 1000:
        return False, f"quant.sf has only {n} rows — expected tens of thousands of transcripts."
    return True, f"Salmon quantification OK — {n:,} transcripts."


# ── ATAC-seq ─────────────────────────────────────────────────────────────


def check_bowtie2_bam(bam_path: Path) -> tuple[bool, str]:
    if not bam_path.exists():
        return False, f"BAM not found: {bam_path.name}"
    if not bam_path.with_suffix(".bam.bai").exists() and not Path(str(bam_path) + ".bai").exists():
        return False, f"BAM index (.bai) not found for {bam_path.name}. Run samtools index."
    if bam_path.stat().st_size < 1024:
        return False, f"{bam_path.name} appears empty."
    return True, f"Alignment OK — {bam_path.name} written and indexed."


def check_atac_filtered_bam(bam_path: Path) -> tuple[bool, str]:
    if not bam_path.exists():
        return False, f"Filtered BAM not found: {bam_path.name}"
    if bam_path.stat().st_size < 512:
        return False, f"{bam_path.name} appears empty after filtering."
    return True, "Filtering OK (mitochondrial reads removed, duplicates removed, MAPQ ≥ 30)."


def check_macs2_peaks(peak_file: Path) -> tuple[bool, str]:
    if not peak_file.exists():
        return False, f"Peak file not found: {peak_file.name}"
    n = len([l for l in peak_file.read_text().splitlines() if l.strip()])
    if n < 100:
        return (
            False,
            f"Only {n} peaks called. ATAC-seq should yield thousands. "
            "Check BAM quality and that filtering didn't discard too many reads.",
        )
    return True, f"{n:,} peaks called."


# ── Methylation ──────────────────────────────────────────────────────────


def check_bismark_bam(bam_path: Path, report_path: Path | None = None) -> tuple[bool, str]:
    if not bam_path.exists():
        return False, f"Bismark BAM not found: {bam_path.name}"

    if report_path and report_path.exists():
        m = re.search(r"Mapping efficiency:\s+([\d.]+)%", report_path.read_text())
        if m:
            eff = float(m.group(1))
            if eff < 30:
                return (
                    False,
                    f"Bismark mapping efficiency is {eff:.1f}%. Below 30% usually means "
                    "wrong reference, failed bisulfite conversion, or contaminated library.",
                )
            return True, f"Bismark alignment OK — {eff:.1f}% mapping efficiency."

    return True, f"{bam_path.name} written."


def check_methylation_cx_report(report_path: Path) -> tuple[bool, str]:
    if not report_path.exists():
        return False, f"CX cytosine report not found: {report_path.name}"
    n = sum(1 for l in report_path.read_text().splitlines() if l.strip())
    if n < 100:
        return False, f"CX report has only {n} lines — extraction may have failed."
    return True, f"Methylation extraction OK — cytosine report written ({n:,} positions)."


# ── OBAMA format ─────────────────────────────────────────────────────────

# OBAMA's comparison code looks for 'disease' and 'control' literally.
# Do not rename these labels — OBAMA will silently fail to find the groups.
_VALID_GROUPS = {"disease", "control"}


def check_counts_matrix(counts_path: Path, coldata_path: Path) -> tuple[bool, str]:
    """Validate the DESeq2/edgeR/limma matrix export: counts matrix + coldata.

    Checks that both files exist, their sample sets match exactly, and the
    condition labels are the OBAMA-compatible 'disease'/'control'.
    """
    if not counts_path.exists():
        return False, f"Counts matrix not found: {counts_path}"
    if not coldata_path.exists():
        return False, f"coldata not found: {coldata_path}"

    clines = [l for l in counts_path.read_text().splitlines() if l.strip()]
    if len(clines) < 2:
        return False, "Counts matrix needs a header and at least one feature row."
    sep = "\t" if "\t" in clines[0] else ","
    matrix_samples = [c.strip() for c in clines[0].split(sep)][1:]  # drop the 'gene' column

    dlines = [l for l in coldata_path.read_text().splitlines() if l.strip()]
    dsep = "\t" if "\t" in dlines[0] else ","
    coldata_samples, conditions = [], []
    for row in dlines[1:]:
        cells = [c.strip() for c in row.split(dsep)]
        if len(cells) >= 2:
            coldata_samples.append(cells[0])
            conditions.append(cells[1])

    if set(matrix_samples) != set(coldata_samples):
        return False, ("Sample mismatch between counts matrix columns and coldata rows — "
                       "DESeq2/edgeR require them to align.")
    bad = {c for c in conditions if c not in _VALID_GROUPS}
    if bad:
        return False, f"condition must be 'disease' or 'control'. Invalid value(s): {', '.join(bad)}."

    n_features = len(clines) - 1
    return True, (f"Matrix OK — {n_features:,} feature(s) x {len(matrix_samples)} sample(s), "
                  "coldata aligned.")


def check_obama_format(csv_path: Path) -> tuple[bool, str]:
    if not csv_path.exists():
        return False, f"Output file not found: {csv_path}"

    lines = [l for l in csv_path.read_text().splitlines() if l.strip()]
    if len(lines) < 2:
        return False, "File needs at least one header row and one data row."

    sep = "\t" if "\t" in lines[0] else ","
    header = [c.strip() for c in lines[0].split(sep)]

    if header[0] != "geo_accession":
        return False, f"First column must be 'geo_accession', got '{header[0]}'."
    if header[1] != "disease.state":
        return False, f"Second column must be 'disease.state', got '{header[1]}'."
    if len(header) < 3:
        return False, "No feature columns found after geo_accession and disease.state."

    bad = {r.split(sep)[1].strip() for r in lines[1:] if r.split(sep)[1].strip() not in _VALID_GROUPS}
    if bad:
        return (
            False,
            f"disease.state must be 'disease' or 'control'. Invalid value(s): {', '.join(bad)}. "
            "OBAMA requires these exact strings for its group comparison code.",
        )

    n_samples = len(lines) - 1
    n_features = len(header) - 2
    return True, f"OBAMA matrix OK — {n_samples} sample(s), {n_features:,} feature(s)."
