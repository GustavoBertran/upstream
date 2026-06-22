"""NCBI E-utilities helpers — GEO series metadata and SRA run resolution.

NCBI allows 3 unauthenticated requests/second.  All calls go through _get()
which enforces a 0.4 s delay (2.5 req/s) to stay safely under that limit.
"""
from __future__ import annotations

import json
import re
import time
import urllib.request
from typing import Optional

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_DELAY  = 0.4   # seconds between NCBI calls


def _lk(d: dict) -> dict:
    """Return d with all top-level keys lowercased — NCBI API key casing varies by dataset."""
    return {k.lower(): v for k, v in d.items()}


def _get(url: str) -> dict:
    time.sleep(_DELAY)
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.loads(r.read())


def fetch_series(gse: str) -> dict:
    """Return series metadata: title, organism, and [{gsm, title}] sample list.

    Raises ValueError if the accession is not found or the response is empty.
    """
    gse = gse.strip().upper()
    if not gse.startswith("GSE"):
        raise ValueError(f"Expected a GSE accession (e.g. GSE12345), got: {gse!r}")

    # Search the GEO DataSets database for this accession
    result = _get(f"{_EUTILS}/esearch.fcgi?db=gds&term={gse}[ACCN]&retmode=json")
    ids = result["esearchresult"]["idlist"]
    if not ids:
        raise ValueError(f"GEO series {gse} not found.")

    # Series UIDs start with "2" in GDS; sample UIDs start with "1"
    series_id = next((i for i in ids if i.startswith("2")), ids[0])

    result = _get(f"{_EUTILS}/esummary.fcgi?db=gds&id={series_id}&retmode=json")
    doc = _lk(result["result"].get(series_id, {}))
    if not doc:
        raise ValueError(f"Could not retrieve summary for {gse}.")

    samples = [
        {"gsm": _lk(s).get("accession", ""), "title": _lk(s).get("title", "")}
        for s in doc.get("samples", [])
    ]
    return {
        "gse":       gse,
        "title":     doc.get("title", ""),
        "organism":  doc.get("taxon", ""),
        "n_samples": doc.get("n_samples", len(samples)),
        "samples":   [s for s in samples if s["gsm"]],
    }


def fetch_srr(gsm: str) -> list[str]:
    """Return SRR run accessions linked to a GSM accession.

    Returns an empty list if no SRA links are found (e.g. microarray study).
    """
    gsm = gsm.strip().upper()

    # Get the GDS numeric UID for this GSM
    result = _get(f"{_EUTILS}/esearch.fcgi?db=gds&term={gsm}[ACCN]&retmode=json")
    ids = result["esearchresult"]["idlist"]
    if not ids:
        return []
    gsm_id = next((i for i in ids if i.startswith("1")), ids[0])

    # Find linked SRA experiments
    result = _get(
        f"{_EUTILS}/elink.fcgi?dbfrom=gds&db=sra&id={gsm_id}&retmode=json"
    )
    sra_ids: list[str] = []
    for linkset in result.get("linksets", []):
        for linksetdb in linkset.get("linksetdbs", []):
            if linksetdb.get("dbto") == "sra":
                sra_ids.extend(str(i) for i in linksetdb.get("links", []))
    if not sra_ids:
        return []

    # Fetch SRA experiment summaries and parse SRR accessions
    ids_str = ",".join(sra_ids[:10])
    result = _get(f"{_EUTILS}/esummary.fcgi?db=sra&id={ids_str}&retmode=json")

    srrs: list[str] = []
    for sra_id in sra_ids[:10]:
        doc = result.get("result", {}).get(sra_id, {})
        # The "runs" field contains XML snippets like: acc="SRR1234567"
        srrs.extend(re.findall(r'acc="(SRR\d+)"', doc.get("runs", "")))

    return list(dict.fromkeys(srrs))  # deduplicate, preserve order


def fetch_characteristics(gse: str) -> dict:
    """Return sample characteristics from the GEO series matrix file.

    Downloads the series matrix from NCBI FTP, decompresses it in streaming
    fashion, and stops parsing once the expression-matrix section begins —
    so only the compact header is read even for large series.

    Returns {gsm: {title, characteristics: {key: value}}}.
    Returns {} if the file is unavailable or cannot be parsed.
    """
    import io as _io
    import zlib as _zlib

    gse = gse.strip().upper()
    digits = gse[3:]
    prefix = digits[:-3] + "nnn" if len(digits) > 3 else digits
    base = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{prefix}/{gse}/matrix"
    )
    primary_url = f"{base}/{gse}_series_matrix.txt.gz"

    # Resolve the actual matrix filename (some series use platform suffixes)
    matrix_url: Optional[str] = None
    time.sleep(_DELAY)
    try:
        req = urllib.request.Request(primary_url, method="HEAD")
        urllib.request.urlopen(req, timeout=12)
        matrix_url = primary_url
    except Exception:
        try:
            with urllib.request.urlopen(f"{base}/", timeout=12) as r:
                listing = r.read().decode("utf-8", errors="replace")
            found = re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', listing)
            if found:
                # href values are bare filenames in NCBI FTP listings
                matrix_url = f"{base}/{found[0]}"
        except Exception:
            return {}

    if not matrix_url:
        return {}

    # Stream + decompress until we hit the expression matrix marker
    header_lines: list[str] = []
    time.sleep(_DELAY)
    try:
        decomp = _zlib.decompressobj(_zlib.MAX_WBITS | 32)  # 32 = auto-detect gzip
        leftover = b""
        stop = False
        with urllib.request.urlopen(matrix_url, timeout=120) as resp:
            while not stop:
                chunk = resp.read(131_072)
                if not chunk:
                    break
                try:
                    raw = decomp.decompress(chunk)
                except _zlib.error:
                    break
                combined = leftover + raw
                nl = combined.rfind(b"\n")
                if nl < 0:
                    leftover = combined
                    continue
                block, leftover = combined[: nl + 1], combined[nl + 1 :]
                for raw_line in block.split(b"\n"):
                    if not raw_line:
                        continue
                    line = raw_line.decode("utf-8", errors="replace").rstrip()
                    if line.startswith("!series_matrix_table_begin"):
                        stop = True
                        break
                    header_lines.append(line)
    except Exception:
        return {}

    # Parse SOFT-format header
    gsm_order: list[str] = []
    char_rows: dict[str, list[str]] = {}   # characteristic key → per-sample values
    sample_titles: dict[str, str] = {}

    for line in header_lines:
        if not line.startswith("!Sample_"):
            continue
        parts = line.split("\t")
        field = parts[0]
        vals = [v.strip('"') for v in parts[1:]]

        if field == "!Sample_geo_accession":
            gsm_order = vals
        elif field == "!Sample_title":
            for i, gsm in enumerate(gsm_order):
                if i < len(vals):
                    sample_titles[gsm] = vals[i]
        elif field == "!Sample_characteristics_ch1":
            # Determine the key from the first non-empty "key: value" entry
            key = next(
                (v.split(":", 1)[0].strip() for v in vals if ":" in v),
                None,
            )
            if key:
                char_rows[key] = [
                    v.split(":", 1)[1].strip() if ":" in v else ""
                    for v in vals
                ]

    return {
        gsm: {
            "title": sample_titles.get(gsm, ""),
            "characteristics": {
                k: (char_rows[k][i] if i < len(char_rows[k]) else "")
                for k in char_rows
            },
        }
        for i, gsm in enumerate(gsm_order)
    }
