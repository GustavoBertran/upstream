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


def build_rnaseq_star(
    samples: list[tuple[str, str, Path]],  # (name, group, ReadsPerGene.out.tab path)
    out_path: Path,
) -> None:
    """Merge STAR ReadsPerGene.out.tab files (unstranded counts) into one OBAMA CSV."""
    genes: list[str] = []
    data: dict[str, dict] = {}

    for name, group, tab_path in samples:
        data[name] = {"group": group, "counts": {}}
        for line in tab_path.read_text().splitlines():
            if line.startswith("N_"):   # skip STAR summary rows (N_unmapped etc.)
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene_id = parts[0]
            count   = float(parts[1])   # column 1 = unstranded read count
            data[name]["counts"][gene_id] = count
            if name == samples[0][0]:
                genes.append(gene_id)

    _write(out_path, samples, genes, lambda name, feat: data[name]["counts"].get(feat, 0.0))


def build_methylation(
    samples: list[tuple[str, str, Path]],  # (name, group, CX_report_path)
    out_path: Path,
    min_coverage: int = 1,
) -> None:
    """Build CpG percent-methylation matrix from Bismark CX reports.

    Only CpGs covered at >= min_coverage reads in EVERY sample are included.
    This intersection avoids imputing 0% methylation for uncovered sites, which
    would generate false differences between samples in OBAMA.
    """
    data: dict[str, dict[str, float]] = {}
    covered_per_sample: list[set[str]] = []

    for name, group, cx_report in samples:
        data[name] = {}
        covered: set[str] = set()
        for line in cx_report.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            chrom, pos, _, coverage, methylated, context, *_ = parts
            if int(coverage) < min_coverage:
                continue
            if "CG" not in context:  # CpG only
                continue
            cpg_id = f"{chrom}:{pos}"
            data[name][cpg_id] = round(int(methylated) / int(coverage) * 100, 2)
            covered.add(cpg_id)
        covered_per_sample.append(covered)

    # Restrict to sites covered in all samples; preserve order from sample[0].
    shared = covered_per_sample[0].intersection(*covered_per_sample[1:]) if len(covered_per_sample) > 1 else covered_per_sample[0]
    cpgs = [cpg for cpg in data[samples[0][0]] if cpg in shared]

    _write(out_path, samples, cpgs, lambda name, feat: data[name][feat])


def build_methylation_array(
    betas_path: Path,
    metadata_path: Path,
    out_path: Path,
) -> None:
    """Build OBAMA CSV from an Illumina 450K / EPIC beta-value matrix.

    Supports two layouts automatically:
    - Probes-as-rows: first column = probe IDs (header cell "ID_REF" or ""),
      remaining header cells = GSM accessions. Standard GEO supplementary format.
    - Samples-as-rows: first column = GSM accessions (unnamed header cell),
      remaining header cells = probe IDs (start with "cg"/"ch"). Typical R
      output after t(exprs(gse)).

    Probes with NA in any selected sample are dropped — imputing 0.0 creates false
    unmethylated signal in OBAMA (same rationale as the WGBS CpG intersection).
    """
    meta: dict[str, str] = {}
    with metadata_path.open() as f:
        for row in csv.DictReader(f):
            acc = row["geo_accession"].strip()
            meta[acc] = row["disease.state"].strip()

    with betas_path.open() as f:
        reader = csv.reader(f)
        header = next(reader)
        probe_prefix = ("cg", "ch")
        # If the second header cell starts with "cg"/"ch", samples are in rows
        samples_as_rows = len(header) > 1 and header[1].strip().startswith(probe_prefix)

        if samples_as_rows:
            # header = ["", "cg00000029", "cg00000108", ...]
            probe_names = [h.strip() for h in header[1:]]
            n = len(probe_names)
            sample_vals: dict[str, list] = {}
            samples_ordered: list[str] = []
            for line in reader:
                if not line:
                    continue
                acc = line[0].strip()
                if acc not in meta:
                    continue
                vals: list = []
                for raw in line[1: n + 1]:
                    raw = raw.strip()
                    if raw and raw not in ("NA", "NaN", "nan"):
                        try:
                            vals.append(float(raw))
                        except ValueError:
                            vals.append(None)
                    else:
                        vals.append(None)
                sample_vals[acc] = vals
                samples_ordered.append(acc)
            valid_idx = [
                i for i in range(n)
                if all(sample_vals[acc][i] is not None for acc in samples_ordered)
            ]
            valid_probes = [probe_names[i] for i in valid_idx]
            out_data = {acc: [sample_vals[acc][i] for i in valid_idx] for acc in samples_ordered}

        else:
            # header = ["ID_REF", "GSM123", "GSM456", ...]
            all_sid = [h.strip() for h in header[1:]]
            col_map = [(i + 1, sid) for i, sid in enumerate(all_sid) if sid in meta]
            samples_ordered = [sid for _, sid in col_map]
            probe_names_all: list[str] = []
            sample_vals = {sid: [] for sid in samples_ordered}
            valid_probes = []
            for line in reader:
                if not line:
                    continue
                probe = line[0].strip()
                row_vals: dict[str, float] = {}
                has_na = False
                for col, sid in col_map:
                    raw = line[col].strip() if col < len(line) else ""
                    if raw and raw not in ("NA", "NaN", "nan"):
                        try:
                            row_vals[sid] = float(raw)
                        except ValueError:
                            has_na = True
                            break
                    else:
                        has_na = True
                        break
                probe_names_all.append(probe)
                if not has_na:
                    valid_probes.append(probe)
                    for sid in samples_ordered:
                        sample_vals[sid].append(row_vals[sid])
            out_data = sample_vals

    if not samples_ordered:
        raise ValueError(
            "No samples in the beta matrix matched the metadata geo_accession values. "
            "Check for leading/trailing spaces in the metadata file."
        )
    if not valid_probes:
        raise ValueError(
            "No probes survived NA filtering. "
            "Verify that the metadata accessions match the beta matrix column headers."
        )

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["geo_accession", "disease.state"] + valid_probes)
        for acc in samples_ordered:
            writer.writerow([acc, meta[acc]] + out_data[acc])


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
