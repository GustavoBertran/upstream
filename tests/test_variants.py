"""Tests for upstream.variants (bcftools stats parsing + variant summary)."""

from pathlib import Path

import pytest

from upstream import checkpoints, variants


@pytest.fixture()
def tmp(tmp_path):
    return tmp_path


_STATS = """\
# This file was produced by bcftools stats
SN\t0\tnumber of samples:\t1
SN\t0\tnumber of records:\t1234
SN\t0\tnumber of SNPs:\t1000
SN\t0\tnumber of indels:\t234
TSTV\t0\t800\t400\t2.05\t790\t395\t2.00
"""


def test_parse_bcftools_stats(tmp):
    p = tmp / "s.stats.txt"
    p.write_text(_STATS)
    s = variants.parse_bcftools_stats(p)
    assert s == {"records": 1234, "snps": 1000, "indels": 234, "ts_tv": 2.05}


def test_parse_bcftools_stats_missing_fields(tmp):
    p = tmp / "s.stats.txt"
    # no indels SN line, no TSTV line
    p.write_text("SN\t0\tnumber of records:\t10\nSN\t0\tnumber of SNPs:\t10\n")
    s = variants.parse_bcftools_stats(p)
    assert s["records"] == 10 and s["indels"] == 0 and s["ts_tv"] == 0.0


def test_parse_bcftools_stats_missing_file(tmp):
    s = variants.parse_bcftools_stats(tmp / "nope.txt")
    assert s == {"records": 0, "snps": 0, "indels": 0, "ts_tv": 0.0}


def test_write_variant_summary(tmp):
    (tmp / "a.stats.txt").write_text(_STATS)
    (tmp / "b.stats.txt").write_text(
        "SN\t0\tnumber of records:\t50\nSN\t0\tnumber of SNPs:\t40\n"
        "SN\t0\tnumber of indels:\t10\nTSTV\t0\t30\t10\t3.0\t.\t.\t.\n"
    )
    results = [
        ("tumor", "disease", tmp / "tumor.vcf.gz", tmp / "a.stats.txt"),
        ("normal", "control", tmp / "normal.vcf.gz", tmp / "b.stats.txt"),
    ]
    out = tmp / "variant_summary.csv"
    variants.write_variant_summary(results, out)

    rows = [l for l in out.read_text().splitlines() if l.strip()]
    assert rows[0] == "sample,group,total_records,snps,indels,ts_tv"
    assert rows[1] == "tumor,disease,1234,1000,234,2.05"
    assert rows[2] == "normal,control,50,40,10,3.0"
    # cross-check: the checkpoint accepts what we wrote
    assert checkpoints.check_variant_summary(out)[0] is True
