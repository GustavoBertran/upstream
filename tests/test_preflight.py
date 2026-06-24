"""Tests for upstream.preflight for the new genomics + proteomics tracks."""

import pytest

from upstream import preflight


@pytest.fixture()
def tmp(tmp_path):
    return tmp_path


def test_required_tools_new_tracks():
    assert preflight.required_tools("genomics") == ["fastp", "bwa", "samtools", "bcftools"]
    assert preflight.required_tools("proteomics") == []


def test_genomics_missing_samples_and_reference():
    issues = preflight.run_issues("genomics", samples=None, reference=None)
    assert any("samplesheet" in i for i in issues)
    assert any("reference" in i and "not provided" in i for i in issues)


def test_genomics_reference_without_bwa_index(tmp):
    fa = tmp / "genome.fa"
    fa.write_text(">chr1\nACGT\n")          # FASTA exists, but no .bwt / .fai sidecars
    issues = preflight.run_issues("genomics", reference=str(fa))
    assert any(".bwt" in i for i in issues)
    assert any(".fai" in i for i in issues)


def test_proteomics_missing_files_no_samplesheet(tmp):
    # validate_files flags provided-but-missing paths (matching the methylation-array path)
    issues = preflight.run_issues("proteomics",
                                  intensities=str(tmp / "nope.txt"),
                                  metadata=str(tmp / "nometa.csv"))
    assert any("intensities" in i for i in issues)
    assert any("metadata" in i for i in issues)
    # proteomics must NOT demand a samplesheet
    assert not any("samplesheet" in i for i in issues)
