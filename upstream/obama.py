"""Build downstream-analysis matrices from pipeline results.

Two output layouts are supported (see write_rnaseq_outputs):

OBAMA format (default) — one CSV, samples as rows:
  - First column:  geo_accession  (sample identifier)
  - Second column: disease.state  (must be exactly 'disease' or 'control')
  - Remaining columns: one per feature (gene, peak, or CpG position)
  The labels 'disease'/'control' are hardcoded expectations in OBAMA's comparison
  code; do not rename them.

Matrix format (for DESeq2 / edgeR / limma-voom) — two files, features as rows:
  - counts_matrix.csv: first column = feature id, one column per sample
  - coldata.csv:       sample id + condition (disease/control)
  R tools read these directly; see content/rnaseq_export_formats.md.

RNA-seq values are RAW COUNTS (what OBAMA expects): Salmon NumReads (optionally
aggregated transcript -> gene via a tx2gene map) and STAR ReadsPerGene counts.
"""

import csv
import gzip
import re
from pathlib import Path
from typing import Optional


# ── tx2gene ────────────────────────────────────────────────────────────────


def _open_text(path: Path):
    """Open a text file, transparently handling gzip (.gz)."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_tx2gene(
    gtf_path: Optional[Path] = None,
    tx2gene_path: Optional[Path] = None,
) -> Optional[dict[str, str]]:
    """Return a {transcript_id: gene} map, or None if no source is given.

    A 2-column tx2gene CSV (transcript_id, gene) takes priority. Otherwise a GTF
    is parsed: each 'transcript' feature maps its transcript_id to gene_name
    (falling back to gene_id). GENCODE GTFs may be gzipped.
    """
    if tx2gene_path:
        m: dict[str, str] = {}
        with _open_text(tx2gene_path) as f:
            for row in csv.reader(f):
                if len(row) < 2:
                    continue
                tx, gene = row[0].strip(), row[1].strip()
                if tx.lower() in ("transcript_id", "tx", "name", "target_id"):
                    continue  # header
                if tx:
                    m[tx] = gene
        return m or None

    if gtf_path:
        m = {}
        tx_re = re.compile(r'transcript_id "([^"]+)"')
        gn_re = re.compile(r'gene_name "([^"]+)"')
        gid_re = re.compile(r'gene_id "([^"]+)"')
        with _open_text(gtf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "transcript":
                    continue
                tx = tx_re.search(parts[8])
                if not tx:
                    continue
                gene = gn_re.search(parts[8]) or gid_re.search(parts[8])
                if gene:
                    m[tx.group(1)] = gene.group(1)
        return m or None

    return None


# ── per-sample count readers ─────────────────────────────────────────────────


def _read_salmon_counts(quant_dir: Path, tx2gene: Optional[dict[str, str]]) -> dict[str, float]:
    """Read Salmon NumReads (estimated counts); aggregate to gene if tx2gene given."""
    counts: dict[str, float] = {}
    with (quant_dir / "quant.sf").open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n = float(row["NumReads"])
            if tx2gene:
                tx = row["Name"]
                gene = tx2gene.get(tx) or tx2gene.get(tx.split(".")[0])  # version fallback
                if gene is None:
                    gene = tx  # unmapped transcript: keep its own id rather than drop reads
                counts[gene] = counts.get(gene, 0.0) + n
            else:
                counts[row["Name"]] = n
    return counts


def _read_star_counts(tab_path: Path) -> dict[str, float]:
    """Read STAR ReadsPerGene.out.tab unstranded counts (column 1), gene-level."""
    counts: dict[str, float] = {}
    for line in tab_path.read_text().splitlines():
        if line.startswith("N_"):  # skip STAR summary rows
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        counts[parts[0]] = float(parts[1])
    return counts


def _collect_rnaseq_counts(
    samples: list[tuple[str, str, Path]],
    source: str,
    tx2gene: Optional[dict[str, str]],
) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Return (ordered feature list from the first sample, {sample: {feature: count}})."""
    data: dict[str, dict[str, float]] = {}
    features: list[str] = []
    for i, (name, _group, path) in enumerate(samples):
        counts = _read_salmon_counts(path, tx2gene) if source == "salmon" else _read_star_counts(path)
        data[name] = counts
        if i == 0:
            features = list(counts.keys())  # dict preserves insertion order
    return features, data


def write_rnaseq_outputs(
    samples: list[tuple[str, str, Path]],
    outdir: Path,
    fmt: str,
    source: str,
    tx2gene: Optional[dict[str, str]] = None,
) -> list[Path]:
    """Write RNA-seq results in the requested format(s); return the files written.

    fmt: 'obama' (default), 'matrix' (counts_matrix.csv + coldata.csv), or 'both'.
    source: 'salmon' or 'star'. Raw counts either way.
    """
    features, data = _collect_rnaseq_counts(samples, source, tx2gene)
    written: list[Path] = []

    if fmt in ("obama", "both"):
        obama_path = outdir / "obama_matrix.csv"
        _write(obama_path, samples, features, lambda name, feat: data[name].get(feat, 0.0))
        written.append(obama_path)

    if fmt in ("matrix", "both"):
        counts_path = outdir / "counts_matrix.csv"
        coldata_path = outdir / "coldata.csv"
        _write_counts_matrix(counts_path, coldata_path, samples, features, data)
        written += [counts_path, coldata_path]

    return written


def _write_counts_matrix(
    counts_path: Path,
    coldata_path: Path,
    samples: list[tuple[str, str, Path]],
    features: list[str],
    data: dict[str, dict[str, float]],
) -> None:
    """Write a feature x sample integer count matrix + a coldata table.

    Counts are rounded to integers so DESeqDataSetFromMatrix() accepts them
    directly (Salmon NumReads are fractional); STAR counts are already integers.
    """
    names = [name for name, _g, _p in samples]
    with counts_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene"] + names)
        for feat in features:
            writer.writerow([feat] + [int(round(data[name].get(feat, 0.0))) for name in names])

    with coldata_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "condition"])
        for name, group, _p in samples:
            writer.writerow([name, group])


def build_rnaseq(
    samples: list[tuple[str, str, Path]],  # (name, group, salmon_quant_dir)
    out_path: Path,
    tx2gene: Optional[dict[str, str]] = None,
) -> None:
    """OBAMA CSV from Salmon NumReads (raw counts), optionally gene-level."""
    features, data = _collect_rnaseq_counts(samples, "salmon", tx2gene)
    _write(out_path, samples, features, lambda name, feat: data[name].get(feat, 0.0))


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
    """OBAMA CSV from STAR ReadsPerGene.out.tab unstranded counts (gene-level)."""
    features, data = _collect_rnaseq_counts(samples, "star", None)
    _write(out_path, samples, features, lambda name, feat: data[name].get(feat, 0.0))


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
