"""CLI smoke tests via --dry-run (no external tools needed).

Under --dry-run, runner.run/pipe/run_capture short-circuit and _require_tools is
skipped, so the genomics/proteomics commands preview their plan without bwa,
bcftools, fastp, etc. being installed.
"""

from pathlib import Path

import pytest
from typer.testing import CliRunner

from upstream.cli import app

runner = CliRunner()


@pytest.fixture()
def tmp(tmp_path):
    return tmp_path


def _samplesheet(tmp, paired=True):
    r1 = tmp / "s_R1.fastq.gz"; r1.write_bytes(b"")
    cols = "name,group,r1,r2\n"
    if paired:
        r2 = tmp / "s_R2.fastq.gz"; r2.write_bytes(b"")
        rows = f"s1,disease,{r1},{r2}\n"
    else:
        rows = f"s1,disease,{r1},\n"
    sheet = tmp / "samples.csv"
    sheet.write_text(cols + rows)
    return sheet


def _reference(tmp, with_index=True):
    fa = tmp / "genome.fa"
    fa.write_text(">chr1\nACGTACGT\n")
    if with_index:
        (tmp / "genome.fa.bwt").write_bytes(b"x")
        (tmp / "genome.fa.fai").write_text("chr1\t8\t6\t8\t9\n")
    return fa


# ── genomics ──────────────────────────────────────────────────────────────


def test_genomics_dry_run(tmp):
    sheet = _samplesheet(tmp, paired=True)
    ref = _reference(tmp)
    res = runner.invoke(app, ["genomics", "--samples", str(sheet), "--reference", str(ref),
                              "--outdir", str(tmp / "out"), "--dry-run", "--no-explain"])
    assert res.exit_code == 0, res.output
    out = res.output
    for token in ["bwa", "mem", "bcftools", "mpileup", "call", "1/5", "5/5"]:
        assert token in out, f"missing {token!r} in dry-run output"


def test_genomics_single_end_dry_run(tmp):
    sheet = _samplesheet(tmp, paired=False)
    ref = _reference(tmp)
    res = runner.invoke(app, ["genomics", "--samples", str(sheet), "--reference", str(ref),
                              "--outdir", str(tmp / "out"), "--dry-run", "--no-explain"])
    assert res.exit_code == 0, res.output
    # single-end markdup path skips collate/fixmate
    assert "fixmate" not in res.output and "collate" not in res.output


def test_genomics_missing_index(tmp, monkeypatch):
    # pretend the tools are installed so the bwa-index check (not the tool check) fires
    monkeypatch.setattr("upstream.preflight.missing_tools", lambda tools: [])
    sheet = _samplesheet(tmp, paired=True)
    ref = _reference(tmp, with_index=False)   # FASTA without .bwt
    res = runner.invoke(app, ["genomics", "--samples", str(sheet), "--reference", str(ref),
                              "--outdir", str(tmp / "out"), "--no-explain"])
    assert res.exit_code == 1
    assert "index-help --tool bwa" in res.output


def test_genomics_gatk_dry_run(tmp):
    sheet = _samplesheet(tmp, paired=True)
    ref = _reference(tmp)
    res = runner.invoke(app, ["genomics", "--samples", str(sheet), "--reference", str(ref),
                              "--outdir", str(tmp / "out"), "--caller", "gatk",
                              "--dry-run", "--no-explain"])
    assert res.exit_code == 0, res.output
    assert "gatk" in res.output and "HaplotypeCaller" in res.output
    assert "mpileup" not in res.output   # gatk path replaces bcftools calling


def test_genomics_bad_caller(tmp):
    sheet = _samplesheet(tmp, paired=True)
    ref = _reference(tmp)
    res = runner.invoke(app, ["genomics", "--samples", str(sheet), "--reference", str(ref),
                              "--outdir", str(tmp / "out"), "--caller", "bogus", "--no-explain"])
    assert res.exit_code == 1


# ── proteomics ──────────────────────────────────────────────────────────────


def test_proteomics_dry_run(tmp):
    # content need not be valid — the dry-run guard returns before parsing
    intens = tmp / "proteinGroups.txt"; intens.write_text("junk\n")
    meta = tmp / "meta.csv"; meta.write_text("junk\n")
    out = tmp / "out"
    res = runner.invoke(app, ["proteomics", "--intensities", str(intens), "--metadata", str(meta),
                              "--outdir", str(out), "--dry-run", "--no-explain"])
    assert res.exit_code == 0, res.output
    assert "would parse" in res.output
    assert not (out / "proteins_matrix.csv").exists()   # nothing written in dry-run


def test_proteomics_bad_normalize(tmp):
    intens = tmp / "proteinGroups.txt"; intens.write_text("junk\n")
    meta = tmp / "meta.csv"; meta.write_text("junk\n")
    res = runner.invoke(app, ["proteomics", "--intensities", str(intens), "--metadata", str(meta),
                              "--outdir", str(tmp / "out"), "--normalize", "bogus", "--no-explain"])
    assert res.exit_code == 1


# ── help grouping ─────────────────────────────────────────────────────────


def test_help_grouping():
    res = runner.invoke(app, ["--help"])
    assert res.exit_code == 0
    for panel in ["Genomics", "Transcriptomics", "Epigenomics", "Proteomics",
                  "Quality control", "Utilities"]:
        assert panel in res.output, f"missing help panel {panel!r}"
