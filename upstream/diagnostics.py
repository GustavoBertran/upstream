"""Environment diagnostics for bug reports.

Collects the details that make an issue actionable — version, OS, Python, conda
env, and which tools are on PATH — so users can paste them into a GitHub issue.
"""
from __future__ import annotations

import os
import platform
import shutil
import sys

from . import __version__

REPO_URL = "https://github.com/GustavoBertran/upstream"
ISSUES_NEW_URL = REPO_URL + "/issues/new"

# Every external tool any track might use.
_ALL_TOOLS = [
    "fastqc", "multiqc", "fastp", "salmon", "STAR", "bowtie2", "samtools",
    "macs2", "bismark", "bismark_methylation_extractor", "seqtk",
    "prefetch", "fasterq-dump", "multiBamSummary",
]


def collect() -> dict:
    return {
        "upstream": __version__,
        "python": sys.version.split()[0],
        "platform": f"{platform.system()} {platform.release()} ({platform.machine()})",
        "conda_env": os.environ.get("CONDA_DEFAULT_ENV", "(none)"),
        "tools": {t: ("found" if shutil.which(t) else "MISSING") for t in _ALL_TOOLS},
    }


def as_text() -> str:
    d = collect()
    lines = [
        f"upstream:  {d['upstream']}",
        f"python:    {d['python']}",
        f"platform:  {d['platform']}",
        f"conda env: {d['conda_env']}",
        "tools on PATH:",
    ]
    lines += [f"  {t}: {status}" for t, status in d["tools"].items()]
    return "\n".join(lines)
