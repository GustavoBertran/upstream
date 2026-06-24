"""VCF/stats summarization for the genomics track — no matrices, no OBAMA.

The genomics track's natural product is a VCF, not a sample x feature matrix, so it
does not flow into OBAMA (OBAMA is transcriptomics-only). This module parses the
`bcftools stats` reports produced per sample into a small, human-readable
variant_summary.csv (one row per sample).
"""
from __future__ import annotations

import csv
from pathlib import Path


def parse_bcftools_stats(stats_path: Path) -> dict:
    """Parse a `bcftools stats` text report → {records, snps, indels, ts_tv}.

    Reads the SN ("summary numbers") and TSTV lines. Missing fields default to 0 /
    0.0 rather than raising — a tiny VCF may omit some lines.
    """
    out = {"records": 0, "snps": 0, "indels": 0, "ts_tv": 0.0}
    if not Path(stats_path).exists():
        return out

    for line in Path(stats_path).read_text().splitlines():
        if line.startswith("SN\t"):
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            key, val = parts[2].strip().lower(), parts[3].strip()
            try:
                n = int(val)
            except ValueError:
                continue
            if key == "number of records:":
                out["records"] = n
            elif key == "number of snps:":
                out["snps"] = n
            elif key == "number of indels:":
                out["indels"] = n
        elif line.startswith("TSTV\t"):
            # TSTV  id  ts  tv  ts/tv  ts(1st ALT)  tv(1st ALT)  ts/tv(1st ALT)
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    out["ts_tv"] = float(parts[4].strip())
                except ValueError:
                    pass
    return out


def write_variant_summary(
    results: list[tuple[str, str, Path, Path]],  # (name, group, vcf, stats_path)
    out_path: Path,
) -> None:
    """Write one summary row per sample (parsed from each sample's bcftools stats)."""
    with Path(out_path).open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "group", "total_records", "snps", "indels", "ts_tv"])
        for name, group, _vcf, stats_path in results:
            s = parse_bcftools_stats(stats_path)
            writer.writerow([name, group, s["records"], s["snps"], s["indels"], s["ts_tv"]])
