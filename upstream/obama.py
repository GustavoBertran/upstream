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

This file also hosts the non-OBAMA matrix builders that share its CSV/coldata
utilities: the methylation-array import and the proteomics intensity-matrix
analysis (write_proteomics_outputs) — both produce a feature x sample matrix +
coldata for limma, never an OBAMA matrix.
"""

import csv
import gzip
import math
import re
import statistics
from pathlib import Path
from typing import Optional


# ── shared helpers for the matrix/coldata export ─────────────────────────────


def _write_coldata(path: Path, rows: list[tuple[str, str]]) -> None:
    """Write a coldata table: sample,condition (condition = disease/control)."""
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "condition"])
        writer.writerows(rows)


def _beta_to_mvalue(beta: float, eps: float = 1e-3) -> float:
    """M-value = log2(beta / (1 - beta)); beta clamped to [eps, 1-eps] to avoid ±inf."""
    b = min(max(beta, eps), 1.0 - eps)
    return round(math.log2(b / (1.0 - b)), 4)


def _counts_to_mvalue(methylated: int, unmethylated: int, alpha: float = 1.0) -> float:
    """WGBS M-value from read counts: log2((M + a) / (U + a)); pseudocount avoids ±inf."""
    return round(math.log2((methylated + alpha) / (unmethylated + alpha)), 4)


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


def salmon_tx2gene_match(quant_dir: Path, tx2gene: dict[str, str]) -> tuple[int, int]:
    """Return (matched, total) quant.sf transcripts that map via tx2gene.

    Lets callers detect a release/ID mismatch (e.g. a non-GENCODE index) that
    would silently leave output transcript-level despite a populated tx2gene map.
    """
    matched = total = 0
    with (quant_dir / "quant.sf").open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            tx = row["Name"]
            if tx in tx2gene or tx.split(".")[0] in tx2gene:
                matched += 1
    return matched, total


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

    _write_coldata(coldata_path, [(name, group) for name, group, _p in samples])


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


def _collect_wgbs_counts(
    samples: list[tuple[str, str, Path]],
    min_coverage: int,
) -> tuple[list[str], dict[str, dict[str, tuple[int, int]]]]:
    """Parse Bismark CX reports → (shared CpGs ordered, {sample: {cpg: (meth, unmeth)}}).

    CX report columns (tab): chrom, pos, strand, count_methylated, count_unmethylated,
    context, trinucleotide. Coverage = methylated + unmethylated. Only CpG-context
    sites covered at >= min_coverage in EVERY sample are kept (the intersection avoids
    imputing values for uncovered sites, which would fabricate between-sample
    differences).
    """
    data: dict[str, dict[str, tuple[int, int]]] = {}
    covered_per_sample: list[set[str]] = []

    for name, _group, cx_report in samples:
        data[name] = {}
        covered: set[str] = set()
        for line in cx_report.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            chrom, pos, _strand, meth, unmeth, context = parts[:6]
            if "CG" not in context:  # CpG only
                continue
            m, u = int(meth), int(unmeth)
            if m + u < min_coverage:
                continue
            cpg_id = f"{chrom}:{pos}"
            data[name][cpg_id] = (m, u)
            covered.add(cpg_id)
        covered_per_sample.append(covered)

    shared = (covered_per_sample[0].intersection(*covered_per_sample[1:])
              if len(covered_per_sample) > 1 else covered_per_sample[0])
    cpgs = [cpg for cpg in data[samples[0][0]] if cpg in shared]
    return cpgs, data


def build_methylation(
    samples: list[tuple[str, str, Path]],  # (name, group, CX_report_path)
    out_path: Path,
    min_coverage: int = 1,
) -> None:
    """OBAMA CpG percent-methylation matrix from Bismark CX reports."""
    cpgs, data = _collect_wgbs_counts(samples, min_coverage)

    def pct(name: str, cpg: str) -> float:
        m, u = data[name][cpg]
        cov = m + u
        return round(m / cov * 100, 2) if cov else 0.0

    _write(out_path, samples, cpgs, pct)


def write_methylation_outputs(
    samples: list[tuple[str, str, Path]],  # (name, group, CX_report_path)
    outdir: Path,
    fmt: str,
    min_coverage: int = 1,
) -> list[Path]:
    """Write WGBS results: OBAMA percent-methylation and/or M-value matrix + coldata.

    fmt: 'obama' (default), 'matrix' (mvalues_matrix.csv + coldata.csv for limma),
    or 'both'. M-values use log2((M+1)/(U+1)) from read counts.
    """
    cpgs, data = _collect_wgbs_counts(samples, min_coverage)
    names = [n for n, _g, _p in samples]
    written: list[Path] = []

    if fmt in ("obama", "both"):
        obama_path = outdir / "obama_matrix.csv"
        _write(obama_path, samples, cpgs,
               lambda name, cpg: (lambda m, u: round(m / (m + u) * 100, 2) if (m + u) else 0.0)(*data[name][cpg]))
        written.append(obama_path)

    if fmt in ("matrix", "both"):
        mpath = outdir / "mvalues_matrix.csv"
        cpath = outdir / "coldata.csv"
        with mpath.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["cpg"] + names)
            for cpg in cpgs:
                writer.writerow([cpg] + [_counts_to_mvalue(*data[name][cpg]) for name in names])
        _write_coldata(cpath, [(n, g) for n, g, _p in samples])
        written += [mpath, cpath]

    return written


def _parse_beta_matrix(
    betas_path: Path,
    metadata_path: Path,
) -> tuple[list[str], list[str], dict[str, list[float]], dict[str, str]]:
    """Parse a 450K/EPIC beta matrix + metadata → (samples, probes, betas, meta).

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

    return samples_ordered, valid_probes, out_data, meta


def _write_beta_obama(out_path, samples_ordered, valid_probes, out_data, meta) -> None:
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["geo_accession", "disease.state"] + valid_probes)
        for acc in samples_ordered:
            writer.writerow([acc, meta[acc]] + out_data[acc])


def build_methylation_array(
    betas_path: Path,
    metadata_path: Path,
    out_path: Path,
) -> None:
    """OBAMA CSV (samples x probes, beta values) from a 450K/EPIC beta matrix."""
    samples_ordered, valid_probes, out_data, meta = _parse_beta_matrix(betas_path, metadata_path)
    _write_beta_obama(out_path, samples_ordered, valid_probes, out_data, meta)


def write_methylation_array_outputs(
    betas_path: Path,
    metadata_path: Path,
    outdir: Path,
    fmt: str,
) -> list[Path]:
    """Write array results: OBAMA beta matrix and/or M-value matrix + coldata (limma).

    M-value = log2(beta/(1-beta)); limma is run on M-values, not betas.
    """
    samples_ordered, valid_probes, out_data, meta = _parse_beta_matrix(betas_path, metadata_path)
    written: list[Path] = []

    if fmt in ("obama", "both"):
        obama_path = outdir / "obama_matrix.csv"
        _write_beta_obama(obama_path, samples_ordered, valid_probes, out_data, meta)
        written.append(obama_path)

    if fmt in ("matrix", "both"):
        mpath = outdir / "mvalues_matrix.csv"
        cpath = outdir / "coldata.csv"
        with mpath.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["probe"] + samples_ordered)
            for i, probe in enumerate(valid_probes):
                writer.writerow([probe] + [_beta_to_mvalue(out_data[acc][i]) for acc in samples_ordered])
        _write_coldata(cpath, [(acc, meta[acc]) for acc in samples_ordered])
        written += [mpath, cpath]

    return written


# ── Proteomics (intensity matrix → limma matrix; non-OBAMA) ──────────────────

# MaxQuant proteinGroups quant-column prefixes, by --intensity-col choice.
_PROT_QUANT_PREFIX = {"LFQ": "LFQ intensity ", "iBAQ": "iBAQ ", "Intensity": "Intensity "}
# Rows flagged '+' in any present column below are dropped (decoys/contaminants).
_PROT_FLAG_COLS = ("Potential contaminant", "Contaminant", "Reverse", "Only identified by site")
_PROT_NA = {"", "NA", "NaN", "nan", "NAN", "#NUM!"}


def _read_proteomics_metadata(path: Path) -> dict[str, str]:
    """Read a proteomics metadata CSV → ordered {sample: group}.

    Requires columns 'sample' and 'group'; group must be 'disease' or 'control'.
    """
    meta: dict[str, str] = {}
    bad: list[str] = []
    with Path(path).open() as f:
        reader = csv.DictReader(f)
        cols = set(reader.fieldnames or [])
        if not {"sample", "group"} <= cols:
            raise ValueError("metadata CSV must have 'sample' and 'group' columns.")
        for row in reader:
            sample = (row.get("sample") or "").strip()
            group = (row.get("group") or "").strip()
            if not sample:
                continue
            if group not in ("disease", "control"):
                bad.append(f"{sample}={group or '(empty)'}")
            meta[sample] = group
    if bad:
        raise ValueError("group must be 'disease' or 'control'. Invalid: " + ", ".join(bad))
    if not meta:
        raise ValueError("metadata CSV has no sample rows.")
    return meta


def _to_intensity(raw: str) -> Optional[float]:
    """Parse one intensity cell: blank/NA/0 → None (0 = not detected in MaxQuant)."""
    raw = (raw or "").strip()
    if raw in _PROT_NA:
        return None
    try:
        v = float(raw)
    except ValueError:
        return None
    return v if v > 0 else None


def _read_maxquant_proteingroups(
    path: Path, intensity_col: str,
) -> tuple[list[str], list[str], dict[str, list[Optional[float]]]]:
    """Parse a MaxQuant proteinGroups.txt → (protein_ids, sample_names, {pid: [vals]}).

    Drops contaminant / reverse / 'only identified by site' rows (any present flag
    column == '+'). Protein id = first accession of 'Majority protein IDs'
    (fallback 'Protein IDs'). Quant columns are detected by the prefix for the
    chosen intensity_col; the sample name is the text after the prefix. Spurious
    matches (e.g. 'iBAQ peptides') are harmless — they are dropped later when the
    sample set is intersected with the metadata.
    """
    prefix = _PROT_QUANT_PREFIX[intensity_col]
    with Path(path).open() as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError("proteinGroups.txt appears empty.")
        idx = {name: i for i, name in enumerate(header)}
        id_col = idx.get("Majority protein IDs", idx.get("Protein IDs"))
        if id_col is None:
            raise ValueError("proteinGroups.txt has no 'Majority protein IDs' / 'Protein IDs' column.")
        flag_cols = [idx[c] for c in _PROT_FLAG_COLS if c in idx]
        quant = [(i, h[len(prefix):].strip()) for i, h in enumerate(header)
                 if h.startswith(prefix) and h[len(prefix):].strip()]
        if not quant:
            raise ValueError(
                f"No '{prefix.strip()}' quant columns found in proteinGroups.txt. "
                f"Check --intensity-col (looked for headers starting '{prefix}')."
            )
        sample_names = [s for _i, s in quant]

        protein_ids: list[str] = []
        data: dict[str, list[Optional[float]]] = {}
        seen: set[str] = set()
        for row in reader:
            if len(row) <= id_col:
                continue
            if any(row[c].strip() == "+" for c in flag_cols if c < len(row)):
                continue
            pid = row[id_col].split(";")[0].strip()
            if not pid or pid in seen:
                continue
            seen.add(pid)
            protein_ids.append(pid)
            data[pid] = [_to_intensity(row[i]) if i < len(row) else None for i, _s in quant]
    return protein_ids, sample_names, data


def _read_generic_matrix(
    path: Path,
) -> tuple[list[str], list[str], dict[str, list[Optional[float]]]]:
    """Parse a generic protein x sample table (CSV or TSV, auto-detected)."""
    text_lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    if len(text_lines) < 2:
        raise ValueError("intensity matrix needs a header and at least one protein row.")
    sep = "\t" if "\t" in text_lines[0] else ","
    header = [c.strip() for c in text_lines[0].split(sep)]
    sample_names = header[1:]
    if not sample_names:
        raise ValueError("intensity matrix has no sample columns after the protein column.")
    protein_ids: list[str] = []
    data: dict[str, list[Optional[float]]] = {}
    seen: set[str] = set()
    for line in text_lines[1:]:
        cells = [c.strip() for c in line.split(sep)]
        pid = cells[0]
        if not pid or pid in seen:
            continue
        seen.add(pid)
        protein_ids.append(pid)
        data[pid] = [_to_intensity(cells[i]) if i < len(cells) else None
                     for i in range(1, len(sample_names) + 1)]
    return protein_ids, sample_names, data


def _proteomics_normalize_median(
    logvals: dict[str, dict[str, Optional[float]]], proteins: list[str], samples: list[str],
) -> dict[str, dict[str, Optional[float]]]:
    """Median-center each sample to the global median of per-sample medians.

    Aligns sample medians (corrects loading differences) while preserving scale —
    after this every sample's median of observed values equals the common median.
    """
    sample_med: dict[str, float] = {}
    for s in samples:
        obs = [logvals[p][s] for p in proteins if logvals[p][s] is not None]
        sample_med[s] = statistics.median(obs) if obs else 0.0
    gmed = statistics.median(list(sample_med.values())) if sample_med else 0.0
    return {
        p: {s: (None if logvals[p][s] is None else logvals[p][s] - (sample_med[s] - gmed))
            for s in samples}
        for p in proteins
    }


def _proteomics_normalize_quantile(
    logvals: dict[str, dict[str, Optional[float]]], proteins: list[str], samples: list[str],
) -> dict[str, dict[str, Optional[float]]]:
    """NA-aware rank-based quantile normalization.

    Builds a reference distribution by interpolating each sample's sorted observed
    values onto a common n-point grid and averaging across samples, then maps every
    value back via its within-sample rank. With a complete (no-missing) matrix this
    reduces to classic quantile normalization (all columns share one sorted vector).
    """
    n = len(proteins)
    ref = [0.0] * n
    contributing = [0] * n

    def interp(sorted_obs: list[float], q: float) -> float:
        m = len(sorted_obs)
        if m == 1:
            return sorted_obs[0]
        idx = q * (m - 1)
        lo = math.floor(idx)
        hi = math.ceil(idx)
        return sorted_obs[lo] * (1 - (idx - lo)) + sorted_obs[hi] * (idx - lo)

    for s in samples:
        obs = sorted(v for v in (logvals[p][s] for p in proteins) if v is not None)
        if not obs:
            continue
        for k in range(n):
            q = k / (n - 1) if n > 1 else 0.0
            ref[k] += interp(obs, q)
            contributing[k] += 1
    ref = [ref[k] / contributing[k] if contributing[k] else 0.0 for k in range(n)]

    out: dict[str, dict[str, Optional[float]]] = {p: {} for p in proteins}
    for s in samples:
        pairs = sorted((logvals[p][s], p) for p in proteins if logvals[p][s] is not None)
        m = len(pairs)
        for rank, (_v, p) in enumerate(pairs):
            q = rank / (m - 1) if m > 1 else 0.0
            out[p][s] = interp(ref, q)
        for p in proteins:
            if logvals[p][s] is None:
                out[p][s] = None
    return out


def write_proteomics_outputs(
    intensities_path: Path,
    metadata_path: Path,
    outdir: Path,
    input_type: str = "maxquant",
    intensity_col: str = "LFQ",
    min_valid: float = 0.5,
    normalize: str = "median",
) -> list[Path]:
    """Proteomics intensity matrix → proteins_matrix.csv + coldata.csv (for limma).

    Real pure-Python analysis: drop contaminant/reverse hits (MaxQuant), keep only
    samples present in the metadata, filter proteins by a per-group valid-value
    fraction, log2-transform, normalize, and write a feature x sample matrix +
    coldata. No imputation (residual missing stays blank; limma tolerates NA).
    No OBAMA output.
    """
    meta = _read_proteomics_metadata(metadata_path)
    if input_type == "maxquant":
        prot_ids, file_samples, data = _read_maxquant_proteingroups(intensities_path, intensity_col)
    else:
        prot_ids, file_samples, data = _read_generic_matrix(intensities_path)

    file_idx = {s: i for i, s in enumerate(file_samples)}
    samples = [s for s in meta if s in file_idx]   # metadata order, intersected
    if not samples:
        raise ValueError(
            "No samples in the intensity table matched the metadata 'sample' values. "
            "Check that sample names match the quant-column suffixes."
        )

    vals = {pid: {s: data[pid][file_idx[s]] for s in samples} for pid in prot_ids}

    groups: dict[str, list[str]] = {}
    for s in samples:
        groups.setdefault(meta[s], []).append(s)

    kept: list[str] = []
    for pid in prot_ids:
        keep = True
        for gsamples in groups.values():
            present = sum(1 for s in gsamples if vals[pid][s] is not None)
            if present < min_valid * len(gsamples) - 1e-9:
                keep = False
                break
        if keep:
            kept.append(pid)
    if not kept:
        raise ValueError(
            "No proteins survived the valid-value filter. "
            f"Lower --min-valid (currently {min_valid}) or check the input."
        )

    logvals = {
        pid: {s: (math.log2(vals[pid][s]) if vals[pid][s] is not None else None) for s in samples}
        for pid in kept
    }
    if normalize == "median":
        logvals = _proteomics_normalize_median(logvals, kept, samples)
    elif normalize == "quantile":
        logvals = _proteomics_normalize_quantile(logvals, kept, samples)

    matrix_path = Path(outdir) / "proteins_matrix.csv"
    coldata_path = Path(outdir) / "coldata.csv"
    with matrix_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["protein"] + samples)
        for pid in kept:
            writer.writerow(
                [pid] + ["" if logvals[pid][s] is None else round(logvals[pid][s], 4) for s in samples]
            )
    _write_coldata(coldata_path, [(s, meta[s]) for s in samples])
    return [matrix_path, coldata_path]


def build_atac_consensus(peak_files: list[Path], out_bed: Path) -> int:
    """Merge all samples' narrowPeak intervals into a consensus BED; return n regions.

    Overlapping/abutting peaks across samples are merged per chromosome so every
    sample is counted over the same feature set (the input DESeq2/edgeR need).
    """
    intervals: dict[str, list[tuple[int, int]]] = {}
    for pf in peak_files:
        for line in Path(pf).read_text().splitlines():
            if not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 3:
                continue
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            intervals.setdefault(p[0], []).append((start, end))

    merged: list[tuple[str, int, int]] = []
    for chrom in sorted(intervals):
        ivs = sorted(intervals[chrom])
        cs, ce = ivs[0]
        for s, e in ivs[1:]:
            if s <= ce:                 # overlap or touch → extend
                ce = max(ce, e)
            else:
                merged.append((chrom, cs, ce))
                cs, ce = s, e
        merged.append((chrom, cs, ce))

    with out_bed.open("w") as f:
        for i, (chrom, s, e) in enumerate(merged, 1):
            f.write(f"{chrom}\t{s}\t{e}\tpeak_{i}\n")
    return len(merged)


def build_atac_counts_matrix(
    raw_counts_path: Path,
    samples: list[tuple[str, str]],   # (name, group) in --labels order
    counts_path: Path,
    coldata_path: Path,
) -> None:
    """Convert deeptools multiBamSummary --outRawCounts into counts_matrix + coldata.

    multiBamSummary preserves --bamfiles/--labels order, so count columns 4..N map
    positionally to `samples`. Counts are rounded to integers for DESeq2/edgeR.
    """
    lines = [l for l in raw_counts_path.read_text().splitlines() if l.strip()]
    if len(lines) < 2:
        raise ValueError("multiBamSummary produced no count rows.")
    n_data = len(lines[0].split("\t")) - 3
    if n_data != len(samples):
        raise ValueError(
            f"multiBamSummary returned {n_data} count column(s) for {len(samples)} sample(s)."
        )

    names = [n for n, _g in samples]
    with counts_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["peak"] + names)
        for line in lines[1:]:
            p = line.split("\t")
            if len(p) < 3 + len(samples):
                continue
            peak_id = f"{p[0]}:{p[1]}-{p[2]}"
            writer.writerow([peak_id] + [int(round(float(p[3 + i]))) for i in range(len(samples))])

    _write_coldata(coldata_path, list(samples))


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
