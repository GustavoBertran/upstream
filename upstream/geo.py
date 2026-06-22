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
    doc = result["result"].get(series_id, {})
    if not doc:
        raise ValueError(f"Could not retrieve summary for {gse}.")

    samples = [
        {"gsm": s["accession"], "title": s["title"]}
        for s in doc.get("samples", [])
    ]
    return {
        "gse":       gse,
        "title":     doc.get("title", ""),
        "organism":  doc.get("taxon", ""),
        "n_samples": doc.get("n_samples", len(samples)),
        "samples":   samples,
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
