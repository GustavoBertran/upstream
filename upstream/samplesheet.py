"""Generate a starter samplesheet from a folder of FASTQ files.

Pairs R1/R2 mates, infers paired vs single-end, and derives a sample name from the
filename. The `group` column is left blank for the user to fill (disease/control).
"""
from __future__ import annotations

import re
from pathlib import Path

# <sample>_R1.fastq.gz, <sample>_R1_001.fastq.gz, <sample>_1.fq.gz, ...
_MATE_RE = re.compile(
    r"^(?P<sample>.+?)_(?P<mate>R?[12])(?:_\d+)?\.(?:fastq|fq)(?:\.gz)?$",
    re.IGNORECASE,
)
_FASTQ_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)


def scan_fastq_dir(directory) -> tuple[list[dict], list[str]]:
    """Return (rows, warnings). rows = [{name, group, r1, r2}] with group blank."""
    d = Path(directory)
    if not d.is_dir():
        return [], [f"not a directory: {directory}"]

    files = sorted(p for p in d.iterdir() if p.is_file() and _FASTQ_RE.search(p.name))
    paired: dict[str, dict[str, str]] = {}   # name -> {"1": path, "2": path}
    singles: dict[str, str] = {}             # name -> path (no R1/R2 token)
    warnings: list[str] = []

    for p in files:
        m = _MATE_RE.match(p.name)
        if m:
            name = m.group("sample")
            mate = m.group("mate").upper().lstrip("R")  # "1" or "2"
            slot = paired.setdefault(name, {})
            if mate in slot:
                warnings.append(
                    f"{name}: multiple R{mate} files — keeping {Path(slot[mate]).name}, "
                    f"ignoring {p.name} (merge lanes yourself if needed)"
                )
            else:
                slot[mate] = str(p.resolve())
        else:
            name = _FASTQ_RE.sub("", p.name)
            if name in singles:
                warnings.append(f"{name}: multiple files without an R1/R2 token — keeping the first")
            else:
                singles[name] = str(p.resolve())

    rows: list[dict] = []
    for name in sorted(paired):
        slot = paired[name]
        r1, r2 = slot.get("1"), slot.get("2", "")
        if not r1 and r2:
            warnings.append(f"{name}: has R2 but no R1 — skipped")
            continue
        rows.append({"name": name, "group": "", "r1": r1 or "", "r2": r2 or ""})

    for name in sorted(singles):
        if name in paired:
            continue  # already represented
        rows.append({"name": name, "group": "", "r1": singles[name], "r2": ""})

    if not rows:
        warnings.append("no FASTQ files found (looked for *.fastq / *.fq, optionally .gz)")
    return rows, warnings
