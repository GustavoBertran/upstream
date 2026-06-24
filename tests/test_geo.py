"""Tests for the GEO supplementary-file listing (matrix-track download)."""

import pytest

from upstream import geo


class _FakeResp:
    def __init__(self, data: bytes):
        self._d = data

    def read(self):
        return self._d

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


# A realistic NCBI Apache directory index for a GEO series suppl/ folder.
_LISTING = b"""<html><head><title>Index of /geo/series/GSE59nnn/GSE59685/suppl</title></head>
<body><h1>Index of /geo/series/GSE59nnn/GSE59685/suppl</h1>
<table>
<tr><th><a href="?C=N;O=D">Name</a></th><th>Size</th></tr>
<tr><td><a href="/geo/series/GSE59nnn/GSE59685/">Parent Directory</a></td></tr>
<tr><td><a href="GSE59685_RAW.tar">GSE59685_RAW.tar</a></td><td>1.2G</td></tr>
<tr><td><a href="GSE59685_betas.csv.gz">GSE59685_betas.csv.gz</a></td><td>40M</td></tr>
<tr><td><a href="GSE59685_signal_intensities.csv.gz">GSE59685_signal_intensities.csv.gz</a></td><td>80M</td></tr>
<tr><td><a href="filelist.txt">filelist.txt</a></td><td>2K</td></tr>
</table></body></html>"""


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    monkeypatch.setattr(geo.time, "sleep", lambda *a, **k: None)


def test_fetch_supplementary_files_parses_listing(monkeypatch):
    monkeypatch.setattr(geo.urllib.request, "urlopen", lambda url, timeout=30: _FakeResp(_LISTING))
    out = geo.fetch_supplementary_files("gse59685")
    names = [f["name"] for f in out["files"]]
    # data files kept; parent dir, sort link, and filelist.txt excluded
    assert names == ["GSE59685_RAW.tar", "GSE59685_betas.csv.gz", "GSE59685_signal_intensities.csv.gz"]
    # urls are absolute under the series suppl/ base
    betas = next(f for f in out["files"] if f["name"] == "GSE59685_betas.csv.gz")
    assert betas["url"] == out["base_url"] + "GSE59685_betas.csv.gz"
    assert out["base_url"].endswith("/GSE59685/suppl/")


def test_fetch_supplementary_files_bad_accession():
    with pytest.raises(ValueError, match="GSE"):
        geo.fetch_supplementary_files("PXD000001")


def test_fetch_supplementary_files_empty(monkeypatch):
    monkeypatch.setattr(geo.urllib.request, "urlopen",
                        lambda url, timeout=30: _FakeResp(b"<html><body>nothing here</body></html>"))
    with pytest.raises(ValueError, match="No supplementary files"):
        geo.fetch_supplementary_files("GSE000000")


def test_supplementary_base_url_prefix():
    assert geo.supplementary_base_url("GSE12345").endswith("/GSE12nnn/GSE12345/suppl/")
    assert geo.supplementary_base_url("GSE1").endswith("/GSE1/GSE1/suppl/")
