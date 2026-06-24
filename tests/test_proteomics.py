"""Tests for the proteomics intensity-matrix analysis in upstream.obama."""

import csv
import statistics
from pathlib import Path

import pytest

from upstream import checkpoints, obama


@pytest.fixture()
def tmp(tmp_path):
    return tmp_path


# ── fixture builders ──────────────────────────────────────────────────────


def _write_maxquant(path, samples, rows, contaminant_header="Potential contaminant",
                    prefix="LFQ intensity "):
    cols = ["Majority protein IDs", contaminant_header, "Reverse",
            "Only identified by site"] + [prefix + s for s in samples]
    lines = ["\t".join(cols)]
    for r in rows:
        cells = [r["id"], r.get("contaminant", ""), r.get("reverse", ""), r.get("site", "")]
        cells += [str(r["vals"][s]) for s in samples]
        lines.append("\t".join(cells))
    Path(path).write_text("\n".join(lines) + "\n")


def _write_matrix(path, samples, rows, sep=","):
    lines = [sep.join(["protein"] + samples)]
    for pid, vals in rows:
        lines.append(sep.join([pid] + [str(v) for v in vals]))
    Path(path).write_text("\n".join(lines) + "\n")


def _write_meta(path, mapping):
    Path(path).write_text("sample,group\n" + "\n".join(f"{s},{g}" for s, g in mapping) + "\n")


def _read_matrix(path):
    """Return (header_samples, {protein: [str cells]})."""
    rows = list(csv.reader(Path(path).open()))
    return rows[0][1:], {r[0]: r[1:] for r in rows[1:]}


# ── MaxQuant parsing ──────────────────────────────────────────────────────


def test_read_maxquant_proteingroups(tmp):
    p = tmp / "proteinGroups.txt"
    _write_maxquant(p, ["A", "B"], [
        {"id": "P1;P1b", "vals": {"A": 1000, "B": 2000}},
        {"id": "CON1", "contaminant": "+", "vals": {"A": 5, "B": 5}},
        {"id": "REV1", "reverse": "+", "vals": {"A": 5, "B": 5}},
        {"id": "P2", "vals": {"A": 0, "B": 800}},      # 0 -> None
    ])
    ids, samples, data = obama._read_maxquant_proteingroups(p, "LFQ")
    assert samples == ["A", "B"]
    assert ids == ["P1", "P2"]                          # contaminant + reverse dropped; first accession only
    assert data["P2"][0] is None and data["P2"][1] == 800.0


def test_read_maxquant_older_contaminant_header(tmp):
    p = tmp / "proteinGroups.txt"
    _write_maxquant(p, ["A", "B"], [
        {"id": "P1", "vals": {"A": 100, "B": 200}},
        {"id": "CON", "contaminant": "+", "vals": {"A": 1, "B": 1}},
    ], contaminant_header="Contaminant")
    ids, _samples, _data = obama._read_maxquant_proteingroups(p, "LFQ")
    assert ids == ["P1"]                                # older 'Contaminant' header still drops it


def test_read_generic_matrix_na_and_zero(tmp):
    p = tmp / "m.csv"
    _write_matrix(p, ["S1", "S2"], [("P1", ["NA", 5]), ("P2", [0, 9])])
    ids, samples, data = obama._read_generic_matrix(p)
    assert samples == ["S1", "S2"] and ids == ["P1", "P2"]
    assert data["P1"][0] is None and data["P1"][1] == 5.0
    assert data["P2"][0] is None                        # 0 -> None


def test_read_generic_matrix_tsv(tmp):
    p = tmp / "m.tsv"
    _write_matrix(p, ["S1", "S2"], [("P1", [10, 20])], sep="\t")
    _ids, samples, _data = obama._read_generic_matrix(p)
    assert samples == ["S1", "S2"]                       # auto-detects tab separator


# ── Full pipeline ─────────────────────────────────────────────────────────


def test_write_proteomics_outputs_maxquant(tmp):
    p = tmp / "proteinGroups.txt"
    _write_maxquant(p, ["A", "B", "C", "D"], [
        {"id": "P1", "vals": {"A": 1000, "B": 1100, "C": 900, "D": 950}},   # all -> kept
        {"id": "CON", "contaminant": "+", "vals": {"A": 1, "B": 1, "C": 1, "D": 1}},
        {"id": "Pdrop", "vals": {"A": 500, "B": 600, "C": 0, "D": 0}},      # control missing -> dropped
        {"id": "P2", "vals": {"A": 300, "B": 320, "C": 280, "D": 310}},     # all -> kept
    ])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "disease"), ("C", "control"), ("D", "control")])

    written = obama.write_proteomics_outputs(p, meta, tmp, input_type="maxquant",
                                             intensity_col="LFQ", min_valid=0.5, normalize="median")
    assert [w.name for w in written] == ["proteins_matrix.csv", "coldata.csv"]
    samples, rows = _read_matrix(tmp / "proteins_matrix.csv")
    assert samples == ["A", "B", "C", "D"]
    assert set(rows) == {"P1", "P2"}                    # contaminant + control-missing dropped
    # values are log2-scaled (P1 ~ log2(1000) ≈ 10), not the raw intensity
    assert all(0 < float(c) < 20 for c in rows["P1"])
    assert checkpoints.check_counts_matrix(tmp / "proteins_matrix.csv", tmp / "coldata.csv")[0] is True


def test_write_proteomics_outputs_generic(tmp):
    m = tmp / "m.csv"
    _write_matrix(m, ["A", "B", "C", "D"], [
        ("P1", [1000, 1100, 900, 950]),
        ("P2", [300, 320, 280, 310]),
    ])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "disease"), ("C", "control"), ("D", "control")])
    written = obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix", normalize="none")
    samples, rows = _read_matrix(tmp / "proteins_matrix.csv")
    assert samples == ["A", "B", "C", "D"] and set(rows) == {"P1", "P2"}


def test_proteomics_min_valid_per_group(tmp):
    m = tmp / "m.csv"
    _write_matrix(m, ["A", "B", "C", "D"], [
        ("P_keep", [100, 200, 150, 250]),       # present in both groups -> kept
        ("P_dropctrl", [100, 200, 0, 0]),        # missing across control -> dropped
        ("P_dropdis", [0, 0, 150, 250]),         # missing across disease -> dropped
    ])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "disease"), ("C", "control"), ("D", "control")])
    obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix", min_valid=0.5, normalize="none")
    _samples, rows = _read_matrix(tmp / "proteins_matrix.csv")
    assert set(rows) == {"P_keep"}


def _col_medians(matrix_path):
    samples, rows = _read_matrix(matrix_path)
    meds = []
    for j in range(len(samples)):
        vals = [float(cells[j]) for cells in rows.values() if cells[j] != ""]
        meds.append(statistics.median(vals))
    return meds


def test_proteomics_normalize_median_aligns(tmp):
    m = tmp / "m.csv"
    # samples deliberately on different scales (S2 ~10x, S3 ~0.5x of S1)
    _write_matrix(m, ["A", "B", "C", "D"], [
        ("P%d" % i, [100 * i, 1000 * i, 50 * i, 400 * i]) for i in range(1, 7)
    ])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "disease"), ("C", "control"), ("D", "control")])
    obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix", min_valid=0.5, normalize="median")
    meds = _col_medians(tmp / "proteins_matrix.csv")
    assert max(meds) - min(meds) < 1e-6           # medians aligned (NOT necessarily ~0)


def test_proteomics_normalize_quantile_invariant(tmp):
    m = tmp / "m.csv"
    _write_matrix(m, ["A", "B", "C", "D"], [
        ("P%d" % i, [100 * i, 1000 * i, 50 * i, 400 * i]) for i in range(1, 7)
    ])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "disease"), ("C", "control"), ("D", "control")])
    obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix", min_valid=0.5, normalize="quantile")
    samples, rows = _read_matrix(tmp / "proteins_matrix.csv")
    cols = [sorted(round(float(cells[j]), 4) for cells in rows.values()) for j in range(len(samples))]
    for c in cols[1:]:
        assert c == cols[0]                       # complete matrix → all columns share one sorted vector


# ── error paths ───────────────────────────────────────────────────────────


def test_proteomics_no_sample_match(tmp):
    m = tmp / "m.csv"
    _write_matrix(m, ["X", "Y"], [("P1", [10, 20])])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "control")])
    with pytest.raises(ValueError, match="matched"):
        obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix")


def test_proteomics_bad_group(tmp):
    m = tmp / "m.csv"
    _write_matrix(m, ["A", "B"], [("P1", [10, 20])])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "tumor"), ("B", "control")])
    with pytest.raises(ValueError, match="tumor"):
        obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix")


def test_proteomics_all_dropped(tmp):
    m = tmp / "m.csv"
    _write_matrix(m, ["A", "B", "C", "D"], [("P1", [100, 0, 0, 0])])
    meta = tmp / "meta.csv"
    _write_meta(meta, [("A", "disease"), ("B", "disease"), ("C", "control"), ("D", "control")])
    with pytest.raises(ValueError, match="survived"):
        obama.write_proteomics_outputs(m, meta, tmp, input_type="matrix", min_valid=0.5)
