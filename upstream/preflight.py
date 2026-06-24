"""Pre-run validation — catch samplesheet, tool, and path problems up front.

Shared by the CLI (`upstream check` + each track's fail-fast preflight) and the web
UI (`/api/check`, called before a run starts). Validators COLLECT all problems and
return them as a list of human-readable strings, rather than dying on the first one.
"""
from __future__ import annotations

import csv
import shutil
from pathlib import Path
from typing import Optional


def required_tools(
    track: str,
    aligner: Optional[str] = None,
    method: Optional[str] = None,
    output_format: Optional[str] = None,
) -> list[str]:
    """External tools that must be on PATH to run *track* with the given options."""
    if track == "qc":
        return ["fastqc", "multiqc"]
    if track == "rnaseq":
        return ["fastp", "STAR" if aligner == "star" else "salmon"]
    if track == "genomics":
        return ["fastp", "bwa", "samtools", "bcftools"]
    if track == "proteomics":
        return []  # pure-Python intensity-matrix analysis
    if track in ("atacseq", "chipseq"):
        tools = ["fastp", "bowtie2", "samtools", "macs2"]
        if output_format in ("matrix", "both"):
            tools.append("multiBamSummary")  # deeptools
        return tools
    if track == "methylation":
        if method == "array":
            return []  # pure-Python beta-matrix reshaping
        return ["fastp", "bismark", "bismark_methylation_extractor", "samtools"]
    if track in ("download", "download_geo"):
        return ["prefetch", "fasterq-dump", "seqtk"]
    return []


def missing_tools(tools: list[str]) -> list[str]:
    return [t for t in tools if shutil.which(t) is None]


def validate_samplesheet(path, allow_input: bool = False) -> list[str]:
    """Return ALL problems with the samplesheet (empty list = valid)."""
    p = Path(path)
    if not p.exists():
        return [f"samplesheet not found: {path}"]
    with p.open() as f:
        reader = csv.DictReader(f)
        cols = set(reader.fieldnames or [])
        missing = {"name", "group", "r1", "r2"} - cols
        if missing:
            return [f"missing column(s): {', '.join(sorted(missing))}"]
        rows = list(reader)
    if not rows:
        return ["samplesheet has no data rows"]

    valid_groups = {"disease", "control"} | ({"input"} if allow_input else set())
    inputs = {(r.get("name") or "").strip()
              for r in rows if (r.get("group") or "").strip() == "input"}
    issues: list[str] = []
    seen: set[str] = set()

    for i, row in enumerate(rows, start=2):  # header is line 1
        name = (row.get("name") or "").strip()
        tag = name or f"row {i}"
        if not name:
            issues.append(f"row {i}: empty sample name")
        else:
            if name in seen:
                issues.append(f"{tag}: duplicate sample name")
            seen.add(name)
            if any(c.isspace() for c in name):
                issues.append(f"{tag}: name contains whitespace (use '_' instead)")

        grp = (row.get("group") or "").strip()
        if grp not in valid_groups:
            issues.append(f"{tag}: invalid group '{grp}' — must be "
                          f"{' / '.join(sorted(valid_groups))}")

        r1 = (row.get("r1") or "").strip()
        if not r1:
            issues.append(f"{tag}: r1 is empty")
        elif not Path(r1).exists():
            issues.append(f"{tag}: r1 file not found: {r1}")

        r2 = (row.get("r2") or "").strip()
        if r2 and not Path(r2).exists():
            issues.append(f"{tag}: r2 file not found: {r2}")

        if allow_input:
            ctl = (row.get("control") or "").strip()
            if ctl and ctl not in inputs:
                issues.append(f"{tag}: control '{ctl}' is not a group=input sample")

    if allow_input and not any(g != "input" for g in
                               [(r.get("group") or "").strip() for r in rows]):
        issues.append("every row is group=input — no ChIP signal samples to call peaks on")

    return issues


def validate_dirs(labeled_paths: dict[str, Optional[str]]) -> list[str]:
    """Each given path must exist as a directory. {label: path}; None/empty skipped."""
    issues = []
    for label, path in labeled_paths.items():
        if path and not Path(path).is_dir():
            issues.append(f"{label} directory not found: {path}")
    return issues


def validate_files(labeled_paths: dict[str, Optional[str]]) -> list[str]:
    issues = []
    for label, path in labeled_paths.items():
        if path and not Path(path).is_file():
            issues.append(f"{label} file not found: {path}")
    return issues


def run_issues(
    track: str,
    samples: Optional[str] = None,
    aligner: Optional[str] = None,
    method: Optional[str] = None,
    output_format: Optional[str] = None,
    salmon_index: Optional[str] = None,
    star_index: Optional[str] = None,
    bismark_genome: Optional[str] = None,
    betas: Optional[str] = None,
    metadata: Optional[str] = None,
    reference: Optional[str] = None,
    intensities: Optional[str] = None,
) -> list[str]:
    """All problems that would block a run of *track* — the shared preflight used by
    both `upstream check` and the web UI's pre-run gate. Empty list = ready to run."""
    issues: list[str] = []
    allow_input = track == "chipseq"
    needs_samplesheet = (
        track in ("qc", "rnaseq", "atacseq", "chipseq", "genomics")
        or (track == "methylation" and method != "array")
    )

    if needs_samplesheet:
        if samples:
            issues += [f"samplesheet: {x}" for x in validate_samplesheet(samples, allow_input)]
        else:
            issues.append("samplesheet: not provided")

    for tool in missing_tools(
        required_tools(track, aligner=aligner, method=method, output_format=output_format)
    ):
        issues.append(f"tool: '{tool}' not found on PATH (activate the environment?)")

    dirs: dict[str, Optional[str]] = {}
    if track == "rnaseq":
        if (aligner or "salmon") == "salmon":
            dirs["Salmon index"] = salmon_index
        else:
            dirs["STAR index"] = star_index
    if track == "methylation" and method != "array":
        dirs["Bismark genome"] = bismark_genome
    issues += validate_dirs(dirs)

    if track == "methylation" and method == "array":
        issues += validate_files({"betas": betas, "metadata": metadata})

    if track == "genomics":
        if reference:
            ref = Path(reference)
            if not ref.is_file():
                issues.append(f"reference FASTA not found: {reference}")
            else:
                if not Path(str(ref) + ".bwt").exists():
                    issues.append(f"reference: bwa index missing ({ref.name}.bwt) — "
                                  "see index-help --tool bwa")
                if not Path(str(ref) + ".fai").exists():
                    issues.append(f"reference: samtools faidx index missing ({ref.name}.fai) — "
                                  f"run samtools faidx {ref}")
        else:
            issues.append("reference: not provided")

    if track == "proteomics":
        issues += validate_files({"intensities": intensities, "metadata": metadata})

    return issues
