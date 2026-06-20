"""Build OBAMA-compatible output matrices from pipeline results.

OBAMA requires:
  - First column:  geo_accession  (sample identifier)
  - Second column: disease.state  (must be exactly 'disease' or 'control')
  - Remaining columns: one per feature (transcript, peak, or CpG position)

The labels 'disease' and 'control' are hardcoded expectations in OBAMA's comparison
code. Do not rename them even when the experiment has no actual disease/control
distinction — OBAMA will silently fail to find the groups otherwise.
"""

import csv
from pathlib import Path


def build_rnaseq(
    samples: list[tuple[str, str, Path]],  # (name, group, salmon_quant_dir)
    out_path: Path,
) -> None:
    """Merge Salmon quant.sf files (TPM column) into one OBAMA CSV."""
    transcripts: list[str] = []
    data: dict[str, dict] = {}

    for name, group, quant_dir in samples:
        sf = quant_dir / "quant.sf"
        data[name] = {"group": group, "tpm": {}}
        for row in csv.DictReader(sf.open(), delimiter="\t"):
            tx = row["Name"]
            data[name]["tpm"][tx] = float(row["TPM"])
            if name == samples[0][0]:
                transcripts.append(tx)

    _write(out_path, samples, transcripts, lambda name, feat: data[name]["tpm"].get(feat, 0.0))


def build_atacseq(
    samples: list[tuple[str, str, Path]],  # (name, group, peaks_narrowPeak_file)
    out_path: Path,
) -> None:
    """Use disease sample's peaks as the feature set; fill scores for all samples."""
    # Use the first sample's peaks as the reference feature set
    ref_sample_name, _, ref_peak_file = samples[0]
    peaks: list[str] = []
    data: dict[str, dict[str, float]] = {}

    for name, group, peak_file in samples:
        data[name] = {}
        for line in peak_file.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            peak_id = f"{parts[0]}:{parts[1]}-{parts[2]}"
            data[name][peak_id] = float(parts[4])  # MACS2 score column
            if name == ref_sample_name:
                peaks.append(peak_id)

    _write(out_path, samples, peaks, lambda name, feat: data[name].get(feat, 0.0))


def build_methylation(
    samples: list[tuple[str, str, Path]],  # (name, group, CX_report_path)
    out_path: Path,
) -> None:
    """Build CpG percent-methylation matrix from Bismark CX reports."""
    cpgs: list[str] = []
    data: dict[str, dict[str, float]] = {}

    for name, group, cx_report in samples:
        data[name] = {}
        for line in cx_report.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _, coverage, methylated, *ctx = parts
            if int(coverage) == 0:
                continue
            context = ctx[0] if ctx else ""
            if "CG" not in context:  # CpG only
                continue
            cpg_id = f"{chrom}:{pos}"
            pct = round(int(methylated) / int(coverage) * 100, 2)
            data[name][cpg_id] = pct
            if name == samples[0][0]:
                cpgs.append(cpg_id)

    _write(out_path, samples, cpgs, lambda name, feat: data[name].get(feat, 0.0))


def _write(
    out_path: Path,
    samples: list[tuple[str, str, Path]],
    features: list[str],
    value_fn,
) -> None:
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["geo_accession", "disease.state"] + features)
        for name, group, _ in samples:
            writer.writerow([name, group] + [value_fn(name, feat) for feat in features])
