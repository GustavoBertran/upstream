"""Tests for htsprep.checkpoints validation functions."""

import json
import textwrap
from pathlib import Path

import pytest

from upstream import checkpoints


@pytest.fixture()
def tmp(tmp_path):
    return tmp_path


# ── FastQC / MultiQC ─────────────────────────────────────────────────────


def test_fastqc_pass(tmp):
    d = tmp / "fastqc"
    d.mkdir()
    (d / "sample_R1_fastqc.html").write_text("<html/>")
    (d / "sample_R2_fastqc.html").write_text("<html/>")
    ok, msg = checkpoints.check_fastqc_output(d)
    assert ok and "2" in msg


def test_fastqc_empty_dir(tmp):
    d = tmp / "fastqc"
    d.mkdir()
    ok, _ = checkpoints.check_fastqc_output(d)
    assert not ok


def test_multiqc_pass(tmp):
    d = tmp / "multiqc"
    d.mkdir()
    (d / "multiqc_report.html").write_text("<html/>")
    ok, _ = checkpoints.check_multiqc_output(d)
    assert ok


def test_multiqc_missing(tmp):
    d = tmp / "multiqc"
    d.mkdir()
    ok, _ = checkpoints.check_multiqc_output(d)
    assert not ok


# ── fastp ─────────────────────────────────────────────────────────────────


def _fastp_json(path: Path, passed: int = 500_000, total: int = 600_000) -> None:
    data = {
        "filtering_result": {"passed_filter_reads": passed},
        "summary": {"before_filtering": {"total_reads": total}},
    }
    path.write_text(json.dumps(data))


def test_fastp_pass(tmp):
    j = tmp / "sample_fastp.json"
    _fastp_json(j)
    ok, msg = checkpoints.check_fastp_trim(j)
    assert ok and "500,000" in msg


def test_fastp_zero_reads(tmp):
    j = tmp / "sample_fastp.json"
    _fastp_json(j, passed=0)
    ok, msg = checkpoints.check_fastp_trim(j)
    assert not ok and "0 reads" in msg


def test_fastp_missing(tmp):
    ok, _ = checkpoints.check_fastp_trim(tmp / "missing.json")
    assert not ok


# ── STAR ──────────────────────────────────────────────────────────────────


def test_star_bam_pass(tmp):
    bam = tmp / "sample_Aligned.sortedByCoord.out.bam"
    bam.write_bytes(b"x" * 2048)
    log = tmp / "sample_Log.final.out"
    log.write_text("Uniquely mapped reads % |\t85.32%\n")
    ok, msg = checkpoints.check_star_bam(bam, log)
    assert ok and "85" in msg


def test_star_bam_low_mapping(tmp):
    bam = tmp / "sample_Aligned.sortedByCoord.out.bam"
    bam.write_bytes(b"x" * 2048)
    log = tmp / "sample_Log.final.out"
    log.write_text("Uniquely mapped reads % |\t42.00%\n")
    ok, msg = checkpoints.check_star_bam(bam, log)
    assert not ok and "42.0" in msg


def test_star_bam_missing(tmp):
    ok, _ = checkpoints.check_star_bam(tmp / "missing.bam")
    assert not ok


# ── Salmon ───────────────────────────────────────────────────────────────


def test_salmon_sf_pass(tmp):
    sf = tmp / "quant.sf"
    lines = ["Name\tLength\tTPM\tNumReads"] + [f"tx{i}\t1000\t1.0\t100" for i in range(5000)]
    sf.write_text("\n".join(lines))
    ok, msg = checkpoints.check_salmon_sf(sf)
    assert ok and "5,000" in msg


def test_salmon_sf_too_few(tmp):
    sf = tmp / "quant.sf"
    sf.write_text("Name\tTPM\ntx1\t1.0\n")
    ok, msg = checkpoints.check_salmon_sf(sf)
    assert not ok


def test_salmon_sf_missing(tmp):
    ok, _ = checkpoints.check_salmon_sf(tmp / "missing.sf")
    assert not ok


# ── OBAMA format ──────────────────────────────────────────────────────────


def _write_obama(path: Path, content: str) -> None:
    path.write_text(textwrap.dedent(content))


def test_obama_valid(tmp):
    p = tmp / "out.csv"
    _write_obama(p, """\
        geo_accession,disease.state,GENE1,GENE2
        s1,disease,1.0,2.0
        s2,control,3.0,4.0
        """)
    ok, msg = checkpoints.check_obama_format(p)
    assert ok and "2 sample(s)" in msg and "2 feature(s)" in msg


def test_obama_bad_first_col(tmp):
    p = tmp / "out.csv"
    _write_obama(p, "sample_id,disease.state,GENE1\ns1,disease,1.0\n")
    ok, msg = checkpoints.check_obama_format(p)
    assert not ok and "geo_accession" in msg


def test_obama_bad_second_col(tmp):
    p = tmp / "out.csv"
    _write_obama(p, "geo_accession,condition,GENE1\ns1,disease,1.0\n")
    ok, msg = checkpoints.check_obama_format(p)
    assert not ok and "disease.state" in msg


def test_obama_invalid_group_label(tmp):
    p = tmp / "out.csv"
    _write_obama(p, "geo_accession,disease.state,GENE1\ns1,tumor,1.0\n")
    ok, msg = checkpoints.check_obama_format(p)
    assert not ok and "tumor" in msg


def test_obama_no_features(tmp):
    p = tmp / "out.csv"
    _write_obama(p, "geo_accession,disease.state\ns1,disease\n")
    ok, _ = checkpoints.check_obama_format(p)
    assert not ok


def test_obama_missing_file(tmp):
    ok, _ = checkpoints.check_obama_format(tmp / "missing.csv")
    assert not ok
