"""
Web UI server for the upstream teaching tool.

Students launch with `upstream serve`, then open http://localhost:<port> in
a browser.  The actual pipeline tools still run on the server; output is
streamed back to the browser in real-time via Server-Sent Events.

Architecture:
  POST /api/run   → spawns `upstream <track> ...` as a subprocess in a thread;
                    returns a run_id
  GET  /api/stream/<run_id> → SSE feed of captured lines
  POST /api/cancel/<run_id> → SIGTERM the subprocess
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import threading
import uuid
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import __version__, diagnostics, preflight
from .samplesheet import scan_fastq_dir

_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]")

app = FastAPI(title="upstream", docs_url=None, redoc_url=None)

# Serve profile logos (and any other static assets) from upstream/static/.
_STATIC_DIR = Path(__file__).parent / "static"
if _STATIC_DIR.is_dir():
    app.mount("/static", StaticFiles(directory=str(_STATIC_DIR)), name="static")

# run_id → {"status": "running"|"done"|"error"|"cancelled",
#            "lines": [{"text": str, "cls": str}],
#            "proc": Popen | None,
#            "exit_code": int | None}
_runs: dict[str, dict] = {}
_runs_lock = threading.Lock()


# ── Request model ─────────────────────────────────────────────────────────


class RunRequest(BaseModel):
    track: str
    samples: Optional[str] = None
    outdir: str = "results/"
    threads: int = 4
    explain: bool = True
    dry_run: bool = False
    # rnaseq aligner choice
    aligner: Optional[str] = None
    # rnaseq output format + gene-level aggregation
    output_format: Optional[str] = None   # obama | matrix | both
    gtf: Optional[str] = None
    tx2gene: Optional[str] = None
    # chipseq
    peak_type: Optional[str] = None       # narrow | broad
    # legacy catalog download
    data_track: Optional[str] = None
    # geo download
    geo_samples: Optional[list] = None
    reads: Optional[int] = None            # subsample target; 0/None for the default, <0 for all reads
    # geo supplementary-matrix download (proteomics / methylation array)
    suppl_url: Optional[str] = None
    suppl_name: Optional[str] = None
    # track-specific index paths
    star_index: Optional[str] = None
    salmon_index: Optional[str] = None
    bowtie2_index: Optional[str] = None
    bismark_genome: Optional[str] = None
    # methylation array
    method: Optional[str] = None
    betas: Optional[str] = None
    metadata_csv: Optional[str] = None
    # genomics
    reference: Optional[str] = None
    merge: Optional[bool] = None
    min_mapq: Optional[int] = None
    min_baseq: Optional[int] = None
    # proteomics
    intensities: Optional[str] = None
    input_type: Optional[str] = None
    intensity_col: Optional[str] = None
    min_valid: Optional[float] = None
    normalize: Optional[str] = None


# ── Routes ────────────────────────────────────────────────────────────────


_DEFAULT_PROFILES = {
    "profiles": [
        {"id": "aog", "name": "Applied Optimization Group (AOG)",
         "blurb": "Transcriptomics and the OBAMA pipelines — QC, RNA-seq, ATAC-seq, Methylation, and Download.",
         "image": "aog.png",
         "tracks": ["qc", "rnaseq", "atacseq", "methylation", "download"]},
        {"id": "all", "name": "Whole department",
         "blurb": "Every pipeline upstream provides, across all domains.",
         "tracks": ["*"]},
    ]
}


def load_profiles() -> dict:
    """Load lab-group profiles from upstream/profiles.json (department-editable).

    Falls back to a built-in default if the file is missing or malformed, so the UI
    always renders. Unknown track ids in a profile are filtered out client-side with
    a visible warning. Profiles are navigation only — they do not restrict the CLI.
    """
    path = Path(__file__).parent / "profiles.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        profs = data.get("profiles")
        if isinstance(profs, list):
            cleaned = [p for p in profs
                       if isinstance(p, dict) and p.get("id") and isinstance(p.get("tracks"), list)]
            if cleaned:
                return {"profiles": cleaned}
    except Exception:
        pass
    return _DEFAULT_PROFILES


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    # json.dumps is valid JS; escape "<" so a blurb can't break out of <script>.
    profiles_js = json.dumps(load_profiles()).replace("<", "\\u003c")
    return (_HTML
            .replace("{{VERSION}}", __version__)
            .replace("{{PROFILES}}", profiles_js))


@app.post("/api/run")
def start_run(req: RunRequest) -> dict:
    run_id = uuid.uuid4().hex[:8]
    with _runs_lock:
        _runs[run_id] = {"status": "running", "lines": [], "proc": None,
                         "procs": set(), "exit_code": None}
    threading.Thread(target=_run_pipeline, args=(run_id, req), daemon=True).start()
    return {"run_id": run_id}


@app.get("/api/stream/{run_id}")
def stream(run_id: str) -> StreamingResponse:
    if run_id not in _runs:
        raise HTTPException(404, "run not found")
    return StreamingResponse(_sse_generator(run_id), media_type="text/event-stream")


def _preflight_issues(req: "RunRequest") -> list[str]:
    if req.track == "download_suppl":
        return []  # pure HTTP fetch — nothing to preflight
    if req.track == "download_geo":
        return [f"tool: '{t}' not found on PATH (activate the environment?)"
                for t in preflight.missing_tools(preflight.required_tools("download_geo"))]
    return preflight.run_issues(
        req.track,
        samples=req.samples,
        aligner=req.aligner, method=req.method, output_format=req.output_format,
        salmon_index=req.salmon_index, star_index=req.star_index,
        bismark_genome=req.bismark_genome,
        reference=req.reference, intensities=req.intensities,
        betas=req.betas, metadata=req.metadata_csv,
    )


@app.post("/api/check")
def api_check(req: RunRequest) -> dict:
    """Pre-run validation for the web UI's Run gate. Returns {issues: [...]}."""
    return {"issues": _preflight_issues(req)}


@app.get("/api/diagnostics")
def api_diagnostics() -> dict:
    """Environment diagnostics for the 'Report a bug' link."""
    return {"text": diagnostics.as_text(), "issues_url": diagnostics.ISSUES_NEW_URL}


@app.get("/api/samplesheet")
def api_samplesheet(dir: str) -> dict:
    """Scan a folder of FASTQs and write a starter samples.csv into it."""
    import csv as _csv
    rows, warnings = scan_fastq_dir(dir)
    if not rows:
        return {"ok": False, "warnings": warnings, "path": None}
    out_path = Path(dir) / "samples.csv"
    with out_path.open("w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=["name", "group", "r1", "r2"])
        w.writeheader()
        w.writerows(rows)
    return {
        "ok": True,
        "path": str(out_path.resolve()),
        "n": len(rows),
        "paired": sum(1 for r in rows if r["r2"]),
        "warnings": warnings,
    }


@app.post("/api/cancel/{run_id}")
def cancel(run_id: str) -> dict:
    procs = []
    with _runs_lock:
        run = _runs.get(run_id)
        if run and run["status"] == "running":
            run["status"] = "cancelled"
            if run.get("proc"):
                procs.append(run["proc"])
            procs.extend(run.get("procs", ()))
    # terminate() outside the lock — it can block briefly
    for p in procs:
        try:
            p.terminate()
        except Exception:
            pass
    return {"ok": True}


# ── Pipeline runner ───────────────────────────────────────────────────────


def _upstream_argv() -> list[str]:
    exe = shutil.which("upstream")
    return [exe] if exe else [sys.executable, "-m", "upstream.cli"]


def _build_cmd(req: RunRequest) -> list[str]:
    if req.track == "download":
        cmd = _upstream_argv() + ["download", "--outdir", req.outdir]
        if req.data_track:
            cmd += ["--track", req.data_track]
        return cmd

    if req.track == "methylation" and req.method == "array":
        cmd = _upstream_argv() + [
            "methylation",
            "--method",   "array",
            "--betas",    req.betas or "",
            "--metadata", req.metadata_csv or "",
            "--outdir",   req.outdir,
        ]
        if req.output_format:
            cmd += ["--format", req.output_format]
        if not req.explain:
            cmd.append("--no-explain")
        if req.dry_run:
            cmd.append("--dry-run")
        return cmd

    # Proteomics takes an intensity matrix + metadata, NOT a samplesheet — handle it
    # before the generic samplesheet branch (which hard-requires --samples).
    if req.track == "proteomics":
        cmd = _upstream_argv() + [
            "proteomics",
            "--intensities", req.intensities or "",
            "--metadata",    req.metadata_csv or "",
            "--outdir",      req.outdir,
            "--input-type",  req.input_type or "maxquant",
            "--intensity-col", req.intensity_col or "LFQ",
            "--normalize",   req.normalize or "median",
        ]
        if req.min_valid is not None:
            cmd += ["--min-valid", str(req.min_valid)]
        if not req.explain:
            cmd.append("--no-explain")
        if req.dry_run:
            cmd.append("--dry-run")
        return cmd

    cmd = _upstream_argv() + [
        req.track,
        "--samples", req.samples,
        "--outdir",  req.outdir,
        "--threads", str(req.threads),
    ]
    if not req.explain:
        cmd.append("--no-explain")
    if req.aligner:
        cmd += ["--aligner", req.aligner]
    if req.method:
        cmd += ["--method", req.method]
    if req.star_index:
        cmd += ["--star-index", req.star_index]
    if req.salmon_index:
        cmd += ["--salmon-index", req.salmon_index]
    if req.bowtie2_index:
        cmd += ["--bowtie2-index", req.bowtie2_index]
    if req.bismark_genome:
        cmd += ["--bismark-genome", req.bismark_genome]
    if req.output_format and req.track in ("rnaseq", "methylation", "atacseq", "chipseq"):
        cmd += ["--format", req.output_format]
    if req.track == "rnaseq":
        if req.gtf:
            cmd += ["--gtf", req.gtf]
        if req.tx2gene:
            cmd += ["--tx2gene", req.tx2gene]
    if req.track == "chipseq" and req.peak_type:
        cmd += ["--peak-type", req.peak_type]
    if req.track == "genomics":
        if req.reference:
            cmd += ["--reference", req.reference]
        if req.merge is False:
            cmd.append("--no-merge")
        if req.min_mapq is not None:
            cmd += ["--min-mapq", str(req.min_mapq)]
        if req.min_baseq is not None:
            cmd += ["--min-baseq", str(req.min_baseq)]
    if req.dry_run:
        cmd.append("--dry-run")
    return cmd


def _run_pipeline(run_id: str, req: RunRequest) -> None:
    if req.track == "download_geo":
        _run_geo_download(run_id, req)
        return
    if req.track == "download_suppl":
        _run_suppl_download(run_id, req)
        return

    run = _runs[run_id]
    cmd = _build_cmd(req)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
    )
    run["proc"] = proc
    for raw in proc.stdout:
        text = _ANSI_RE.sub("", raw).rstrip()
        if text:
            with _runs_lock:
                run["lines"].append({"text": text, "cls": _classify(text)})
    proc.wait()
    with _runs_lock:
        run["exit_code"] = proc.returncode
        if run["status"] == "running":
            run["status"] = "done" if proc.returncode == 0 else "error"


def _run_geo_download(run_id: str, req: RunRequest) -> None:
    """Download selected GEO samples as subsampled FASTQ, then write samples.csv.

    Faster, more robust path than bare `fasterq-dump`:
      prefetch SRR (cloud mirror)  →  fasterq-dump on the local .sra  →  seqtk
    `prefetch` pulls from the AWS/GCP SRA mirrors and avoids the NCBI streaming
    throttling that bare `fasterq-dump` hits.  Samples download concurrently.
    Library layout (paired vs single) is detected from the extracted FASTQ.
    """
    import csv as _csv
    import concurrent.futures as _futures

    run     = _runs[run_id]
    outdir  = Path(req.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    samples = req.geo_samples or []
    # reads: <0 → keep all reads (needed for variant calling, where 1M reads is far too
    # little coverage); >0 → that many; None/0 → the teaching default.
    full_reads = req.reads is not None and req.reads < 0
    nreads = req.reads if (req.reads and req.reads > 0) else 1_000_000
    workers = min(4, len(samples)) or 1

    def append(text: str, cls: str = "out") -> None:
        with _runs_lock:
            run["lines"].append({"text": text, "cls": cls})

    def cancelled() -> bool:
        with _runs_lock:
            return run["status"] == "cancelled"

    needed = ["prefetch", "fasterq-dump"] + ([] if full_reads else ["seqtk"])
    missing = [t for t in needed if shutil.which(t) is None]
    if missing:
        append(f"✗ Required tool(s) not found on PATH: {', '.join(missing)}. "
               f"Activate the environment (conda activate upstream) and retry.", "error")
        with _runs_lock:
            run["status"] = "error"; run["exit_code"] = 1
        return

    def spawn(cmd: list[str], tag: str) -> int:
        """Run a subprocess, prefixing output with *tag* (samples run concurrently)."""
        proc = subprocess.Popen(
            [str(c) for c in cmd],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        with _runs_lock:
            run["procs"].add(proc)
        try:
            for raw in proc.stdout:
                t = _ANSI_RE.sub("", raw).rstrip()
                if t:
                    append(f"  [{tag}] {t}")
            proc.wait()
        finally:
            with _runs_lock:
                run["procs"].discard(proc)
        return proc.returncode

    def find_sra(base: Path, srr: str) -> Optional[Path]:
        for cand in (base / srr / f"{srr}.sra", base / f"{srr}.sra", base / srr / srr):
            if cand.exists():
                return cand
        hits = list(base.glob(f"{srr}*/*.sra")) + list(base.glob(f"{srr}*.sra"))
        return hits[0] if hits else None

    def subsample(src: Path, dest: Path) -> bool:
        if full_reads:                       # keep every read — just compress
            with dest.open("wb") as fh:
                p = subprocess.Popen(["gzip", "-c", str(src)], stdout=fh)
                p.wait()
            return p.returncode == 0
        p1 = subprocess.Popen(["seqtk", "sample", "-s", "42", str(src), str(nreads)],
                              stdout=subprocess.PIPE)
        with dest.open("wb") as fh:
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fh)
        p1.stdout.close(); p2.wait(); p1.wait()
        return p2.returncode == 0 and p1.returncode == 0

    def do_sample(s: dict) -> dict:
        srr, name, group = s["srr"], s["name"], s["group"]
        r1 = outdir / f"{name}_R1.fastq.gz"
        r2 = outdir / f"{name}_R2.fastq.gz"
        base = {"name": name, "group": group}

        if r1.exists():
            layout = "paired" if r2.exists() else "single"
            append(f"[{name}] files already present, skipping ({layout}).", "ok")
            return {**base, "r1": r1, "r2": (r2 if r2.exists() else None), "ok": True}
        if cancelled():
            return {**base, "ok": False}

        sra_dir = outdir / f".sra_{srr}"
        try:
            # 1 — prefetch the run from the cloud mirror
            append(f"[{name}] prefetch {srr} (cloud mirror)…")
            if spawn(["prefetch", srr, "-O", str(sra_dir), "--max-size", "100g"], name) != 0:
                if not cancelled():
                    append(f"[{name}] ✗ prefetch failed.", "error")
                return {**base, "ok": False}
            if cancelled():
                return {**base, "ok": False}
            sra = find_sra(sra_dir, srr)
            if sra is None:
                append(f"[{name}] ✗ prefetch produced no .sra file.", "error")
                return {**base, "ok": False}

            # 2 — extract locally (no network); layout detected from output files
            append(f"[{name}] extracting reads…")
            if spawn(["fasterq-dump", str(sra), "--split-files",
                      "--outdir", str(sra_dir)], name) != 0:
                if not cancelled():
                    append(f"[{name}] ✗ fasterq-dump failed.", "error")
                return {**base, "ok": False}
            if cancelled():
                return {**base, "ok": False}

            f1, f2 = sra_dir / f"{srr}_1.fastq", sra_dir / f"{srr}_2.fastq"
            f0 = sra_dir / f"{srr}.fastq"
            if f1.exists() and f2.exists():
                layout, jobs = "paired", [(f1, r1), (f2, r2)]
            elif f0.exists():
                layout, jobs = "single", [(f0, r1)]
            else:
                append(f"[{name}] ✗ no FASTQ produced by fasterq-dump.", "error")
                return {**base, "ok": False}

            # 3 — subsample + compress
            append(f"[{name}] " + (f"compressing full read set ({layout})…" if full_reads
                    else f"subsampling to {nreads:,} reads ({layout})…"))
            for src, dest in jobs:
                if cancelled():
                    return {**base, "ok": False}
                if not subsample(src, dest):
                    dest.unlink(missing_ok=True)
                    append(f"[{name}] ✗ seqtk/gzip failed.", "error")
                    return {**base, "ok": False}

            append(f"[{name}] ✓ ready ({layout}).", "ok")
            return {**base, "r1": r1, "r2": (r2 if layout == "paired" else None), "ok": True}
        except Exception as e:                       # noqa: BLE001 — surface, don't crash the run
            append(f"[{name}] ✗ {e}", "error")
            return {**base, "ok": False}
        finally:
            shutil.rmtree(sra_dir, ignore_errors=True)

    append(f"Downloading {len(samples)} sample(s), {workers} in parallel "
           f"(prefetch → fasterq-dump → "
           f"{'full reads' if full_reads else f'subsample {nreads:,}'})…", "step")

    results: list[dict] = []
    with _futures.ThreadPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(do_sample, samples))

    if cancelled():
        append("Run cancelled.", "error")
        with _runs_lock:
            run["exit_code"] = 1
        return

    failed = [r["name"] for r in results if not r.get("ok")]
    if failed:
        append(f"✗ Failed: {', '.join(failed)}. No samplesheet written.", "error")
        with _runs_lock:
            run["status"] = "error"; run["exit_code"] = 1
        return

    has_single = any(not r.get("r2") for r in results)
    csv_rows = [{"name": r["name"], "group": r["group"],
                 "r1": str(Path(r["r1"]).resolve()),
                 "r2": (str(Path(r["r2"]).resolve()) if r.get("r2") else "")}
                for r in results]

    csv_path = outdir / "samples.csv"
    with csv_path.open("w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=["name", "group", "r1", "r2"])
        w.writeheader(); w.writerows(csv_rows)

    append(f"✓ samples.csv written → {csv_path}", "ok")
    if has_single:
        append("Note: single-end sample(s) were written with an empty r2 column. "
               "The pipeline tracks currently assume paired-end input — running them "
               "on single-end data is a separate step.", "out")
    append("Done. Run the pipeline with:", "success")
    append(f"  upstream rnaseq --samples {csv_path} --outdir results/", "cmd")
    with _runs_lock:
        run["status"] = "done"; run["exit_code"] = 0


def _run_suppl_download(run_id: str, req: RunRequest) -> None:
    """Download a GEO supplementary file (or any direct URL) → outdir; gunzip if needed.

    The matrix-track download (proteomics intensity matrix, methylation-array beta
    matrix) — a plain HTTP fetch, no SRA tools. The matrix parsers aren't gzip-aware,
    so a .gz download is decompressed to the file the pipeline expects.
    """
    import gzip as _gzip
    import urllib.request as _u

    run = _runs[run_id]
    outdir = Path(req.outdir or "data")
    outdir.mkdir(parents=True, exist_ok=True)

    def append(text: str, cls: str = "out") -> None:
        with _runs_lock:
            run["lines"].append({"text": text, "cls": cls})

    def cancelled() -> bool:
        with _runs_lock:
            return run["status"] == "cancelled"

    url = (req.suppl_url or "").strip()
    if not url:
        append("✗ No file URL provided.", "error")
        with _runs_lock:
            run["status"] = "error"; run["exit_code"] = 1
        return

    name = (req.suppl_name or url.rsplit("/", 1)[-1] or "download").split("?")[0]
    dest = outdir / name
    append(f"Downloading {name} → {dest}…", "step")
    try:
        with _u.urlopen(url, timeout=120) as resp:
            total = int(resp.headers.get("Content-Length") or 0)
            got, last_pct = 0, -1
            with dest.open("wb") as fh:
                while True:
                    if cancelled():
                        append("Run cancelled.", "error")
                        dest.unlink(missing_ok=True)
                        with _runs_lock:
                            run["exit_code"] = 1
                        return
                    chunk = resp.read(262_144)
                    if not chunk:
                        break
                    fh.write(chunk)
                    got += len(chunk)
                    if total:
                        pct = int(got * 100 / total)
                        if pct >= last_pct + 10:
                            last_pct = pct
                            append(f"  {pct}%  ({got:,} / {total:,} bytes)")
        append(f"✓ Downloaded {got:,} bytes.", "ok")
    except Exception as e:                       # noqa: BLE001 — surface, don't crash the run
        append(f"✗ Download failed: {e}", "error")
        dest.unlink(missing_ok=True)
        with _runs_lock:
            run["status"] = "error"; run["exit_code"] = 1
        return

    final = dest
    if dest.suffix == ".gz":
        unzipped = dest.with_suffix("")          # drop the .gz extension
        append(f"Decompressing → {unzipped.name}…")
        try:
            with _gzip.open(dest, "rb") as gz, unzipped.open("wb") as out:
                shutil.copyfileobj(gz, out)
            dest.unlink(missing_ok=True)
            final = unzipped
            append("✓ Decompressed.", "ok")
        except Exception as e:                   # noqa: BLE001
            append(f"✗ Could not decompress (keeping the .gz): {e}", "error")
            final = dest

    append(f"✓ Saved → {final.resolve()}", "success")
    append("This is a data matrix — pair it with a metadata CSV and run, e.g.:", "out")
    append(f"  upstream proteomics  --intensities {final} --metadata meta.csv --outdir results/", "cmd")
    append(f"  upstream methylation --method array --betas {final} --metadata meta.csv --outdir results/", "cmd")
    with _runs_lock:
        run["status"] = "done"; run["exit_code"] = 0


def _classify(text: str) -> str:
    t = text.strip()
    if re.search(r"\d+/\d+ —", t):
        return "step"
    if t.startswith("✓"):
        return "ok"
    if t.startswith("✗") or re.search(r"[Ee]rror:", t):
        return "error"
    if t.startswith("$"):
        return "cmd"
    if re.match(r"\[bold green\]", t) or ("Done." in t and len(t) < 80):
        return "success"
    return "out"


# ── GEO metadata routes ───────────────────────────────────────────────────


@app.get("/api/geo/{gse}")
def geo_series(gse: str) -> dict:
    from . import geo
    try:
        return geo.fetch_series(gse)
    except ValueError as e:
        raise HTTPException(404, str(e))


@app.get("/api/geo/srr/{gsm}")
def geo_srr(gsm: str) -> dict:
    from . import geo
    return {"gsm": gsm, "srrs": geo.fetch_srr(gsm)}


@app.get("/api/geo/srr-all/{gse}")
def geo_srr_all(gse: str) -> dict:
    """Return {gsm: first_srr} for all samples in a series via batch API calls."""
    from . import geo
    try:
        return geo.fetch_srr_all(gse)
    except Exception as e:
        raise HTTPException(500, str(e))


@app.get("/api/geo/char/{gse}")
def geo_char(gse: str) -> dict:
    """Return sample characteristics from the GEO series matrix file."""
    from . import geo
    return geo.fetch_characteristics(gse)


@app.get("/api/geo/suppl/{gse}")
def geo_suppl(gse: str) -> dict:
    """List a GEO series' supplementary files (intensity/beta matrices)."""
    from . import geo
    try:
        return geo.fetch_supplementary_files(gse)
    except ValueError as e:
        raise HTTPException(404, str(e))


@app.get("/api/browse")
def browse(path: str = "~") -> dict:
    """List directory contents for the in-browser file picker."""
    import os
    p = Path(path).expanduser()
    if not p.exists():
        p = p.parent
    if not p.is_dir():
        p = p.parent
    if not p.exists():
        p = Path.home()

    entries = []
    try:
        for child in sorted(p.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
            if child.name.startswith("."):
                continue
            try:
                entries.append({"name": child.name, "path": str(child), "is_dir": child.is_dir()})
            except (PermissionError, OSError):
                pass
    except (PermissionError, OSError):
        pass

    parent = str(p.parent) if p != p.parent else None
    return {"path": str(p), "parent": parent, "entries": entries}


# ── SSE generator ─────────────────────────────────────────────────────────


async def _sse_generator(run_id: str):
    run = _runs[run_id]
    sent = 0
    import asyncio
    while True:
        await asyncio.sleep(0.05)
        with _runs_lock:
            batch = run["lines"][sent:]
            status = run["status"]
            exit_code = run["exit_code"]
        for line in batch:
            yield f"data: {json.dumps(line)}\n\n"
            sent += 1
        if status != "running" and sent >= len(run["lines"]):
            yield f"data: {json.dumps({'done': True, 'exit_code': exit_code, 'status': status})}\n\n"
            return


# ── HTML ──────────────────────────────────────────────────────────────────

_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>upstream</title>
<style>
:root{
  --bg:#1a1b26;--bg2:#24283b;--bg3:#1f2335;
  --border:#414868;--text:#a9b1d6;--bright:#c0caf5;
  --accent:#7aa2f7;--green:#9ece6a;--red:#f7768e;
  --yellow:#e0af68;--cyan:#7dcfff;--dim:#8a93bd;
}
*{box-sizing:border-box;margin:0;padding:0}
/* visible focus ring for keyboard users (mouse focus stays subtle) */
:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:3px}
/* respect users who prefer reduced motion */
@media (prefers-reduced-motion:reduce){
  *{transition:none!important}
  .spin{animation:none!important}
}
body{background:var(--bg);color:var(--text);font-family:system-ui,sans-serif;
     height:100vh;display:flex;overflow:hidden}
/* sidebar */
.sidebar{width:220px;min-width:220px;background:var(--bg2);
         border-right:1px solid var(--border);display:flex;
         flex-direction:column;padding:1.25rem .875rem;gap:.75rem;overflow-y:auto}
.sidebar h1{font-size:1rem;color:var(--accent);font-weight:700;
            letter-spacing:.05em;padding-bottom:.75rem;
            border-bottom:1px solid var(--border)}
.sidebar h2{font-size:.65rem;text-transform:uppercase;letter-spacing:.1em;
            color:var(--dim);margin-bottom:.3rem}
.track-domain{margin:.3rem 0 .25rem}
.track-domain:first-child{margin-top:0}
.track-btn{background:none;border:1px solid transparent;border-radius:6px;
           color:var(--text);cursor:pointer;padding:.45rem .65rem;
           text-align:left;width:100%;transition:all .15s}
.track-btn:hover{background:var(--bg3);border-color:var(--border)}
.track-btn.active{background:var(--bg3);border-color:var(--accent);color:var(--bright)}
.track-btn .tn{font-weight:600;font-size:.825rem}
.track-btn .td{font-size:.68rem;color:var(--dim);margin-top:.15rem;line-height:1.3}
.sidebar-footer{margin-top:auto;font-size:.68rem;color:var(--dim);line-height:1.6;
                padding-top:.75rem;border-top:1px solid var(--border)}
/* lab-group badge + switcher (sidebar) */
.lab-badge{display:flex;align-items:center;gap:.5rem;font-size:.7rem;color:var(--dim);
           margin:-.4rem 0 .2rem;flex-wrap:wrap}
.lab-badge #lab-name{color:var(--bright);font-weight:600}
.lab-logo{width:18px;height:18px;border-radius:4px;background:#fff;object-fit:contain;flex-shrink:0}
.lab-switch{background:none;border:none;color:var(--accent);cursor:pointer;font-size:.7rem;padding:0}
.lab-switch:hover{text-decoration:underline}
/* profile picker (landing) */
.profile-picker{position:fixed;inset:0;z-index:200;background:var(--bg);
                display:none;align-items:center;justify-content:center;padding:2rem;overflow-y:auto}
.profile-picker.open{display:flex}
body.picker-open .sidebar,body.picker-open .main{display:none}
.picker-box{width:100%;max-width:560px}
.picker-title{font-size:1.4rem;color:var(--accent);font-weight:700;letter-spacing:.05em;margin-bottom:.4rem}
.picker-sub{color:var(--dim);font-size:.9rem;margin-bottom:1.4rem}
.picker-cards{display:flex;flex-direction:column;gap:.7rem}
.profile-card{display:flex;align-items:center;gap:.95rem;width:100%;text-align:left;background:var(--bg2);
              border:1px solid var(--border);border-radius:8px;padding:1rem 1.1rem;
              cursor:pointer;color:var(--text);transition:border-color .15s,background .15s}
.profile-card:hover,.profile-card:focus-visible{border-color:var(--accent);background:var(--bg3)}
.profile-card .pc-logo{width:48px;height:48px;border-radius:8px;background:#fff;
                       object-fit:contain;flex-shrink:0}
.profile-card .pc-text{min-width:0}
.profile-card .pc-name{font-weight:700;font-size:1rem;color:var(--bright)}
.profile-card .pc-blurb{font-size:.8rem;color:var(--dim);margin-top:.3rem;line-height:1.5}
.profile-card .pc-meta{font-size:.7rem;color:var(--accent);margin-top:.5rem}
.profile-card[disabled]{opacity:.45;cursor:not-allowed}
.picker-warn{color:var(--yellow);font-size:.72rem;margin-top:1rem}
/* main */
.main{flex:1;display:flex;flex-direction:column;padding:1.75rem;gap:1.25rem;
      overflow-y:auto;min-width:0}
/* form */
.panel-title{font-size:.95rem;color:var(--bright);font-weight:600;margin-bottom:1rem}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:.875rem}
.field{display:flex;flex-direction:column;gap:.35rem}
.field.full{grid-column:1/-1}
.field label{font-size:.68rem;text-transform:uppercase;letter-spacing:.06em;color:var(--dim)}
.field input,.field select{
  background:var(--bg2);border:1px solid var(--border);border-radius:4px;
  color:var(--bright);font-family:monospace;font-size:.825rem;
  padding:.45rem .65rem;outline:none;transition:border-color .15s}
.field input:focus,.field select:focus{border-color:var(--accent)}
.hint{font-size:.68rem;color:var(--dim);margin-top:.25rem}
.hint code{background:var(--bg3);padding:.05rem .3rem;border-radius:3px}
.actions{display:flex;gap:.625rem;align-items:center;margin-top:1.125rem}
.btn{background:var(--accent);border:none;border-radius:6px;color:#1a1b26;
     cursor:pointer;font-size:.825rem;font-weight:700;padding:.5rem 1.25rem;
     transition:opacity .15s;white-space:nowrap}
.btn:hover{opacity:.85}
.btn.danger{background:var(--red)}
.btn.ghost{background:none;border:1px solid var(--border);color:var(--text)}
.btn:disabled{opacity:.4;cursor:not-allowed}
/* toggle */
.toggle-row{display:flex;align-items:center;gap:.625rem}
.toggle{position:relative;display:inline-block;width:32px;height:18px}
.toggle input{opacity:0;width:0;height:0}
.slider{position:absolute;inset:0;background:var(--border);border-radius:18px;
        cursor:pointer;transition:background .2s}
.slider::before{content:"";position:absolute;width:12px;height:12px;border-radius:50%;
                background:white;left:3px;top:3px;transition:transform .2s}
input:checked+.slider{background:var(--accent)}
input:checked+.slider::before{transform:translateX(14px)}
/* geo download panel */
#geo-panel{display:none}
.geo-fetch-row{display:flex;gap:.5rem;align-items:flex-end;margin-bottom:.875rem}
.geo-fetch-row .field{flex:1;max-width:300px}
#geo-msg{font-size:.74rem;min-height:1.1em;margin-bottom:.5rem}
#geo-series-info{font-size:.78rem;color:var(--dim);margin-bottom:.75rem}
.geo-table-wrap{overflow-x:auto;border:1px solid var(--border);border-radius:6px;margin-bottom:1rem}
.geo-table{width:100%;border-collapse:collapse;font-size:.775rem}
.geo-table th,.geo-table td{padding:.4rem .7rem;text-align:left;border-bottom:1px solid var(--border)}
.geo-table tr:last-child td{border-bottom:none}
.geo-table th{color:var(--dim);text-transform:uppercase;font-size:.65rem;
              letter-spacing:.07em;background:var(--bg3)}
.geo-table select{background:var(--bg2);border:1px solid var(--border);
                  border-radius:3px;color:var(--bright);font-size:.775rem;padding:.2rem .4rem}
.srr-cell{font-family:monospace;color:var(--cyan);font-size:.72rem}
.srr-err{color:var(--red);font-size:.72rem}
.geo-dl-row{display:flex;align-items:flex-end;gap:.875rem;flex-wrap:wrap}
/* bulk-assign + search toolbar */
.geo-bulk{display:flex;flex-wrap:wrap;gap:.4rem;align-items:center;margin-bottom:.6rem}
.geo-search{flex:1;min-width:180px;max-width:340px;background:var(--bg2);
            border:1px solid var(--border);border-radius:4px;color:var(--bright);
            font-size:.78rem;padding:.32rem .55rem;outline:none}
.geo-search:focus{border-color:var(--accent)}
.geo-bulk-lbl{font-size:.72rem;color:var(--dim);margin-left:.3rem}
.geo-bulk-count{font-size:.7rem;color:var(--dim);margin-left:auto}
.geo-bulk .btn{font-size:.72rem;padding:.28rem .6rem}
.geo-bulk .btn.grp-disease{border-color:var(--red);color:var(--red)}
.geo-bulk .btn.grp-control{border-color:var(--cyan);color:var(--cyan)}
/* log panel */
#log-panel{flex:1;display:flex;flex-direction:column;gap:.75rem;min-height:0}
.log-header{display:flex;justify-content:space-between;align-items:flex-start;flex-shrink:0}
.log-header h2{font-size:.95rem;color:var(--bright);font-weight:600}
.steps{display:flex;gap:.4rem;align-items:center;margin-top:.4rem}
.dot{width:8px;height:8px;border-radius:50%;background:var(--border);transition:background .25s}
.dot.done{background:var(--green)}
.dot.active{background:var(--accent)}
.dot.err{background:var(--red)}
.steps-legend{font-size:.66rem;color:var(--dim);margin-top:.3rem}
/* per-step timing */
.step-time{font-weight:normal;color:var(--bright);margin-left:.6rem;font-size:.72rem}
.step-time.running{color:var(--accent)}
.step-typical{color:var(--dim);font-size:.68rem;font-weight:normal;
              margin:.05rem 0 .15rem 1.4rem}
.log{flex:1;background:var(--bg2);border:1px solid var(--border);border-radius:6px;
     font-family:monospace;font-size:.775rem;overflow-y:auto;padding:.875rem;
     line-height:1.65;min-height:0}
.ll{white-space:pre-wrap;word-break:break-all}
.ll.step{color:var(--cyan);font-weight:bold;margin-top:.5rem}
.ll.ok{color:var(--green)}
.ll.error{color:var(--red)}
.ll.cmd{color:var(--dim)}
.ll.success{color:var(--green);font-weight:bold}
.status-bar{font-size:.72rem;color:var(--dim);display:flex;align-items:center;
            gap:.4rem;flex-shrink:0}
@keyframes spin{to{transform:rotate(360deg)}}
.spin{display:inline-block;animation:spin .9s linear infinite}
/* file browser */
.field-row{display:flex;gap:.4rem;align-items:stretch}
.field-row input{flex:1;min-width:0}
.browse-btn{background:var(--bg2);border:1px solid var(--border);border-radius:4px;
            color:var(--dim);cursor:pointer;padding:0 .6rem;font-size:.75rem;
            white-space:nowrap;transition:all .15s;flex-shrink:0}
.browse-btn:hover{border-color:var(--accent);color:var(--bright)}
.browser-modal{display:none;position:fixed;inset:0;background:rgba(0,0,0,.65);
               z-index:100;align-items:center;justify-content:center}
.browser-modal.open{display:flex}
.browser-box{background:var(--bg2);border:1px solid var(--border);border-radius:8px;
             width:560px;max-width:calc(100vw - 2rem);max-height:80vh;
             display:flex;flex-direction:column;overflow:hidden}
.browser-header{padding:.65rem 1rem;border-bottom:1px solid var(--border);
                display:flex;gap:.5rem;align-items:center}
.browser-path{font-size:.72rem;font-family:monospace;color:var(--bright);
              flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.browser-bar{padding:.4rem .75rem;border-bottom:1px solid var(--border)}
.browser-list{flex:1;overflow-y:auto;max-height:52vh;padding:.25rem 0}
.browser-row{display:flex;align-items:center;gap:.5rem;padding:.32rem .85rem;
             cursor:pointer;font-size:.78rem;color:var(--text)}
.browser-row:hover{background:var(--bg3);color:var(--bright)}
.browser-icon{flex-shrink:0;opacity:.7}
/* GEO table: sort + facet filters */
.geo-sort{cursor:pointer;user-select:none}
.geo-sort:hover{color:var(--bright)}
.sort-ind{font-size:.6rem;margin-left:.2rem;color:var(--accent);vertical-align:middle}
#geo-facets{display:none;flex-wrap:wrap;gap:.35rem;margin-bottom:.6rem}
.facet-group{position:relative;display:inline-block}
.facet-btn{background:var(--bg2);border:1px solid var(--border);border-radius:4px;
           color:var(--dim);cursor:pointer;font-size:.7rem;padding:.22rem .6rem;
           white-space:nowrap;transition:border-color .15s,color .15s}
.facet-btn:hover,.facet-btn.active{border-color:var(--accent);color:var(--bright)}
.facet-popup{display:none;position:absolute;top:calc(100% + 3px);left:0;z-index:60;
             min-width:155px;max-height:215px;overflow-y:auto;background:var(--bg2);
             border:1px solid var(--border);border-radius:5px;padding:.2rem 0;
             box-shadow:0 4px 14px rgba(0,0,0,.45)}
.facet-popup.open{display:block}
.facet-ctrl{display:flex;gap:.35rem;padding:.25rem .65rem .3rem;
            border-bottom:1px solid var(--border);margin-bottom:.1rem}
.facet-ctrl-btn{background:none;border:none;color:var(--accent);cursor:pointer;
                font-size:.68rem;padding:0}
.facet-ctrl-btn:hover{text-decoration:underline}
.facet-item{display:flex;align-items:center;gap:.4rem;padding:.2rem .65rem;font-size:.74rem}
.facet-item:hover{background:var(--bg3)}
.facet-count{margin-left:auto;font-size:.64rem;color:var(--dim)}
/* GEO selection summary */
.geo-summary-hdr{display:flex;align-items:center;justify-content:space-between;
                  margin:.8rem 0 .4rem;padding-top:.6rem;border-top:1px solid var(--border)}
.geo-summary-hdr span{font-size:.78rem;color:var(--dim)}
.grp-disease{color:var(--red)}
.grp-control{color:var(--cyan)}
/* CSV column selector */
.csv-export-hdr{display:flex;align-items:center;justify-content:space-between;
                padding-top:.55rem;margin-top:.6rem;border-top:1px solid var(--border)}
.csv-export-hdr>span{font-size:.76rem;color:var(--dim)}
#geo-csv-cols{display:flex;flex-wrap:wrap;gap:.3rem;margin-top:.4rem}
.csv-col-tag{display:inline-flex;align-items:center;gap:.28rem;background:var(--bg3);
             border:1px solid var(--border);border-radius:4px;
             padding:.18rem .5rem;font-size:.7rem;cursor:pointer;user-select:none;
             transition:border-color .12s,color .12s}
.csv-col-tag input[type=checkbox]{width:11px;height:11px;margin:0;cursor:pointer;
                                   accent-color:var(--accent)}
.csv-col-tag.checked{border-color:var(--accent);color:var(--bright)}
/* collapsible info note (paired vs single-end) */
.info-note{background:var(--bg3);border:1px solid var(--border);border-radius:6px;
           margin-bottom:.875rem;font-size:.76rem;overflow:hidden}
.info-note summary{cursor:pointer;padding:.5rem .75rem;color:var(--accent);
                   user-select:none;list-style:none;font-weight:600}
.info-note summary::-webkit-details-marker{display:none}
.info-note summary::before{content:"\\25B8";display:inline-block;margin-right:.45rem;
                           font-size:.7rem;transition:transform .15s}
.info-note[open] summary::before{transform:rotate(90deg)}
.info-note[open] summary{border-bottom:1px solid var(--border)}
.info-note-body{padding:.65rem .85rem;color:var(--text);line-height:1.65}
.info-note-body p{margin-bottom:.5rem}
.info-note-body p:last-child{margin-bottom:0}
.info-note-body strong{color:var(--bright)}
.info-note-body code{background:var(--bg2);padding:.05rem .3rem;border-radius:3px;
                     font-size:.72rem}
</style>
</head>
<body>

<!-- Lab-group picker (landing) — shown until a profile is chosen -->
<div class="profile-picker" id="profile-picker" role="dialog" aria-modal="true" aria-label="Select your lab group">
  <div class="picker-box">
    <h1 class="picker-title"><span aria-hidden="true">&#9889;</span> upstream</h1>
    <p class="picker-sub">Select your lab group to see its pipelines.</p>
    <div id="picker-cards" class="picker-cards"></div>
    <div id="picker-warn" class="picker-warn" role="alert"></div>
  </div>
</div>

<div class="sidebar">
  <h1><span aria-hidden="true">&#9889;</span> upstream</h1>
  <div class="lab-badge">
    <img id="lab-logo" class="lab-logo" alt="" style="display:none">
    <span id="lab-name"></span>
    <button class="lab-switch" onclick="switchLab()">&#8646; switch lab</button>
  </div>
  <nav aria-label="Pipeline track">
    <h2 id="track-heading">Track</h2>
    <div id="track-list" aria-labelledby="track-heading"
         style="display:flex;flex-direction:column;gap:.2rem;margin-top:.4rem"></div>
  </nav>
  <div class="sidebar-footer">
    Pipeline runs on the server.<br>
    Output streams here live.<br>
    <span style="opacity:.65">upstream v{{VERSION}}</span><br>
    <a href="#" onclick="reportBug();return false;"
       style="color:var(--accent);text-decoration:none">&#9873; Report a bug</a>
  </div>
</div>

<!-- File browser modal -->
<div class="browser-modal" id="browser-modal">
  <div class="browser-box" role="dialog" aria-modal="true" aria-label="File browser">
    <div class="browser-header">
      <span class="browser-path" id="browser-path" aria-label="Current folder"></span>
      <button class="btn ghost" onclick="closeBrowser()" style="padding:.3rem .65rem;font-size:.72rem">Cancel</button>
    </div>
    <div class="browser-bar" id="browser-bar">
      <button class="btn" onclick="selectBrowserDir()" style="padding:.3rem .75rem;font-size:.75rem">Select this folder</button>
    </div>
    <div class="browser-list" id="browser-list" role="listbox" aria-label="Folder contents"></div>
  </div>
</div>

<div class="main">

  <!-- Setup panel -->
  <div id="setup-panel">
    <div class="panel-title" id="form-title">Select a track to begin</div>

    <!-- Standard fields (hidden for download track) -->
    <div class="grid" id="standard-grid">
      <div class="field full">
        <label>Samples CSV</label>
        <div class="field-row">
          <input id="inp-samples" type="text" placeholder="/data/upstream/shared/fastq/rnaseq/samples.csv">
          <button class="browse-btn" aria-label="Browse for samples CSV file" onclick="openBrowser('inp-samples','file')">...</button>
          <button class="browse-btn" onclick="buildSamplesheet()"
                  title="Generate samples.csv from a folder of FASTQ files"
                  aria-label="Generate samples.csv from a folder of FASTQ files">&#9881; from folder</button>
        </div>
        <div class="hint">Generate with: <code>upstream download --track &lt;track&gt; --outdir .</code></div>
      </div>
      <div class="field">
        <label>Output directory</label>
        <div class="field-row">
          <input id="inp-outdir" type="text" placeholder="results/">
          <button class="browse-btn" aria-label="Browse for output directory" onclick="openBrowser('inp-outdir','dir')">...</button>
        </div>
      </div>
      <div class="field">
        <label>CPU threads</label>
        <input id="inp-threads" type="number" value="4" min="1" max="64">
      </div>
      <div id="extra-fields" style="display:contents"></div>
      <div class="field full">
        <div class="toggle-row">
          <label class="toggle">
            <input type="checkbox" id="inp-explain" aria-label="Show step explanations" checked>
            <span class="slider"></span>
          </label>
          <span style="font-size:.825rem">Show step explanations</span>
        </div>
      </div>
      <div class="field full">
        <div class="toggle-row">
          <label class="toggle">
            <input type="checkbox" id="inp-dryrun" aria-label="Dry run — print the commands, run nothing">
            <span class="slider"></span>
          </label>
          <span style="font-size:.825rem">Dry run — print the commands, run nothing</span>
        </div>
      </div>
    </div>

    <!-- GEO download panel (shown only for download track) -->
    <div id="geo-panel">
      <div class="field full" style="max-width:560px">
        <label>What to download</label>
        <select id="inp-dl-mode" onchange="updateDownloadMode()">
          <option value="fastq">Sequencing reads — FASTQ from SRA (genomics &middot; RNA-seq &middot; ATAC/ChIP &middot; WGBS)</option>
          <option value="matrix">Data matrix — GEO supplementary / URL (proteomics &middot; methylation array)</option>
        </select>
      </div>

      <div id="dl-fastq-section">
      <details class="info-note">
        <summary>Paired-end vs single-end reads</summary>
        <div class="info-note-body">
          <p><strong>Paired-end</strong> sequencing reads each DNA fragment from
          <em>both</em> ends, producing two FASTQ files per sample &mdash; <code>R1</code>
          (forward) and <code>R2</code> (reverse). Knowing both ends and the distance
          between them improves alignment accuracy and helps detect splice junctions,
          insertions/deletions, and structural rearrangements. SRA labels these
          <code>PAIRED</code>, and download tools split them into <code>_1.fastq</code>
          and <code>_2.fastq</code>.</p>
          <p><strong>Single-end</strong> sequencing reads each fragment from
          <em>one</em> end only, producing a single FASTQ file per sample. It is cheaper
          and faster and is perfectly adequate for straightforward read counting, such as
          standard RNA-seq quantification or basic peak calling. SRA labels these
          <code>SINGLE</code>.</p>
          <p>Neither is universally better: paired-end gives more positional information
          per fragment, while single-end costs less. The right choice depends on the
          experiment &mdash; the GEO/SRA record for a dataset tells you which was used.</p>
        </div>
      </details>
      <div class="geo-fetch-row">
        <div class="field">
          <label>GEO Series Accession</label>
          <input id="inp-gse" type="text" placeholder="GSE12345"
                 style="text-transform:uppercase"
                 onkeydown="if(event.key==='Enter')fetchGSE()">
        </div>
        <button class="btn ghost" id="btn-fetch" onclick="fetchGSE()">Fetch metadata</button>
      </div>
      <div id="geo-msg" aria-live="polite" style="color:var(--red)"></div>

      <div id="geo-results" style="display:none">
        <div id="geo-series-info"></div>
        <div id="geo-facets"></div>
        <div class="geo-bulk">
          <input id="geo-search" type="text" class="geo-search"
                 placeholder="Filter by GSM, title, or characteristic…"
                 aria-label="Filter samples" oninput="geoSearch(this.value)">
          <span class="geo-bulk-lbl">Assign shown:</span>
          <button class="btn ghost grp-disease" onclick="bulkAssign('disease')">&rarr; disease</button>
          <button class="btn ghost grp-control" onclick="bulkAssign('control')">&rarr; control</button>
          <button class="btn ghost" onclick="bulkAssign('skip')">&rarr; skip</button>
          <span id="geo-bulk-count" class="geo-bulk-count"></span>
        </div>
        <div class="geo-table-wrap">
          <table class="geo-table" id="geo-table">
            <thead id="geo-thead"></thead>
            <tbody id="geo-tbody"></tbody>
          </table>
        </div>
        <!-- Selection summary (appears once any group is assigned) -->
        <div id="geo-summary" style="display:none">
          <div class="geo-summary-hdr">
            <span id="geo-summary-title"></span>
            <button class="btn ghost" onclick="resetGeoGroups()"
                    style="font-size:.72rem;padding:.25rem .65rem">Reset groups</button>
          </div>
          <div class="geo-table-wrap">
            <table class="geo-table">
              <thead>
                <tr><th>Sample (GSM)</th><th>Title</th><th>Group</th><th>SRR</th></tr>
              </thead>
              <tbody id="geo-summary-tbody"></tbody>
            </table>
          </div>
          <!-- CSV export -->
          <div id="geo-csv-export" style="display:none">
            <div class="csv-export-hdr">
              <span>Export columns:</span>
              <button class="btn ghost" onclick="downloadSelectionCSV()"
                      style="font-size:.73rem;padding:.25rem .75rem">&#8595; Download CSV</button>
            </div>
            <div id="geo-csv-cols"></div>
          </div>
        </div>

        <div class="geo-dl-row">
          <div class="field" style="max-width:280px">
            <label>Output directory</label>
            <div class="field-row">
              <input id="inp-geo-outdir" type="text" value="data/">
              <button class="browse-btn" aria-label="Browse for output directory" onclick="openBrowser('inp-geo-outdir','dir')">...</button>
            </div>
          </div>
          <div class="field" style="max-width:175px">
            <label>Reads / sample</label>
            <input id="inp-geo-reads" type="text" value="1000000">
            <div class="hint"><code>all</code> = full coverage (use for genomics)</div>
          </div>
          <button class="btn" onclick="startGeoDownload()">Download selected</button>
          <span id="geo-dl-msg" role="alert" style="font-size:.74rem;color:var(--red)"></span>
        </div>
      </div>
      </div><!-- /dl-fastq-section -->

      <div id="dl-matrix-section" style="display:none">
        <details class="info-note">
          <summary>Downloading a data matrix (proteomics / methylation array)</summary>
          <div class="info-note-body">
            <p>Proteomics and microarray (450K/EPIC) data are not sequencing reads in
            SRA &mdash; they are deposited as <strong>supplementary files</strong> on GEO
            (a protein intensity matrix, or a beta-value matrix). Enter a GEO series to
            list and download one, or paste a direct URL (e.g. from PRIDE / ProteomeXchange
            or a journal). A <code>.gz</code> file is decompressed for you.</p>
          </div>
        </details>
        <div class="geo-fetch-row">
          <div class="field">
            <label>GEO Series Accession</label>
            <input id="inp-suppl-gse" type="text" placeholder="GSE59685"
                   style="text-transform:uppercase"
                   onkeydown="if(event.key==='Enter')listSuppl()">
          </div>
          <button class="btn ghost" id="btn-suppl" onclick="listSuppl()">List files</button>
        </div>
        <div id="suppl-msg" aria-live="polite" style="color:var(--red)"></div>
        <div id="suppl-list" style="display:none">
          <div class="geo-table-wrap">
            <table class="geo-table">
              <thead><tr><th>Supplementary file</th><th></th></tr></thead>
              <tbody id="suppl-tbody"></tbody>
            </table>
          </div>
        </div>
        <div class="field" style="max-width:560px;margin-top:.7rem">
          <label>Or paste a direct file URL</label>
          <input id="inp-suppl-url" type="text"
                 placeholder="https://… supplementary matrix (PRIDE / journal)">
        </div>
        <div class="geo-dl-row">
          <div class="field" style="max-width:280px">
            <label>Output directory</label>
            <div class="field-row">
              <input id="inp-suppl-outdir" type="text" value="data/">
              <button class="browse-btn" aria-label="Browse for output directory" onclick="openBrowser('inp-suppl-outdir','dir')">...</button>
            </div>
          </div>
          <button class="btn" onclick="downloadSupplUrl()">Download from URL</button>
          <span id="suppl-dl-msg" role="alert" style="font-size:.74rem;color:var(--red)"></span>
        </div>
      </div>
    </div>

    <!-- Run Pipeline button (hidden for download track) -->
    <div class="actions" id="run-actions" style="display:none">
      <button class="btn" onclick="startRun()">Run Pipeline</button>
      <span id="form-msg" role="alert" style="font-size:.74rem;color:var(--red)"></span>
    </div>
  </div>

  <!-- Log panel (hidden initially) -->
  <div id="log-panel" style="display:none">
    <div class="log-header">
      <div>
        <h2 id="log-title"></h2>
        <div class="steps" id="steps-row"></div>
        <div class="steps-legend" id="steps-legend" style="display:none">
          &#9201; live elapsed per step &middot; <em>typical</em> ranges are rough and assume full-size input
        </div>
      </div>
      <div class="actions">
        <button class="btn ghost" id="btn-newrun" onclick="newRun()" style="display:none">New run</button>
        <button class="btn danger" id="btn-stop" onclick="stopRun()">Stop</button>
      </div>
    </div>
    <div class="log" id="log" role="log" aria-live="polite" aria-label="Pipeline output"></div>
    <div class="status-bar" id="log-status" role="status" aria-live="polite"><span class="spin" aria-hidden="true">&#8635;</span>&nbsp;Running&hellip;</div>
  </div>

</div>

<script>
const TRACKS = {
  qc:          {label:"QC",           group:"Quality control", desc:"FastQC + MultiQC", nsteps:2, extra:[]},
  genomics:    {label:"Genomics",     group:"Genomics",
                desc:"fastp → BWA-MEM → markdup → bcftools → VCF + summary", nsteps:5,
                isGenomics:true, extra:["reference"]},
  rnaseq:      {label:"RNA-seq",      group:"Transcriptomics",
                desc:"fastp → STAR or Salmon → counts matrix (OBAMA / DESeq2)", nsteps:3, isRnaseq:true},
  atacseq:     {label:"ATAC-seq",     group:"Epigenomics",
                desc:"fastp → Bowtie2 → filter → MACS2 → matrix (OBAMA / DESeq2)", nsteps:5, extra:["bowtie2-index"]},
  chipseq:     {label:"ChIP-seq",     group:"Epigenomics",
                desc:"fastp → Bowtie2 → filter → MACS2 (±input) → counts (DESeq2/edgeR)", nsteps:5, extra:["bowtie2-index"]},
  methylation: {label:"Methylation",  group:"Epigenomics",
                desc:"WGBS (Bismark) or 450K/EPIC array (GEO beta matrix)", nsteps:4, isMethylation:true},
  proteomics:  {label:"Proteomics",   group:"Proteomics",
                desc:"intensity matrix → filter → log2 → normalize → limma matrix", nsteps:1, isProteomics:true},
  download:    {label:"Download Data",group:"Utilities", desc:"Browse GEO, pick conditions, fasterq-dump", isDownload:true},
};

// Sidebar grouping order (domains with no visible track are skipped).
const DOMAINS = ["Quality control", "Genomics", "Transcriptomics", "Epigenomics", "Proteomics", "Utilities"];

const EXTRA_INFO = {
  "star-index":    {label:"STAR index directory",     ph:"/data/upstream/shared/indices/star_hg38"},
  "salmon-index":  {label:"Salmon index directory",   ph:"/data/upstream/shared/indices/salmon_hg38"},
  "bowtie2-index": {label:"Bowtie2 index prefix",     ph:"/data/upstream/shared/indices/bowtie2_hg38/hg38"},
  "bismark-genome":{label:"Bismark genome directory", ph:"/data/upstream/shared/indices/bismark_hg38"},
  "reference":     {label:"Reference genome FASTA",   ph:"/data/upstream/shared/indices/bwa_hg38/hg38.fa", btype:"file"},
};

let track = null, runId = null, sse = null, stepIdx = 0, _runNsteps = 0;

// ── Lab-group profiles (department-editable, navigation only) ───────────────
const PROFILE_CFG = {{PROFILES}};
const PROFILES = (PROFILE_CFG && PROFILE_CFG.profiles) || [];
let _currentProfile = null;

function _profileById(id) {
  for (var i = 0; i < PROFILES.length; i++) if (PROFILES[i].id === id) return PROFILES[i];
  return null;
}

// A profile's logo: a full URL/absolute path as-is, else a file under /static/.
function _profileImgSrc(image) {
  if (!image) return null;
  if (image.indexOf("http://") === 0 || image.indexOf("https://") === 0 ||
      image.charAt(0) === "/") return image;
  return "/static/" + image;
}

// Resolve a profile's track ids to the ones that actually exist ("*" = all);
// unknown ids are skipped with a visible warning rather than silently vanishing.
function resolveTracks(profile) {
  var ids = (profile.tracks || []).slice();
  if (ids.indexOf("*") !== -1) ids = Object.keys(TRACKS);
  var valid = [], bad = [];
  ids.forEach(function(id){ (TRACKS[id] ? valid : bad).push(id); });
  if (bad.length) console.warn("profile '"+profile.id+"' references unknown track id(s): "+bad.join(", "));
  return {valid: valid, bad: bad};
}

// Build the sidebar scoped to a set of track ids, grouped by domain.
function buildSidebar(allowed) {
  const tl = document.getElementById("track-list");
  tl.innerHTML = "";
  for (const dom of DOMAINS) {
    const entries = Object.entries(TRACKS).filter(function(e){
      return e[1].group === dom && allowed.indexOf(e[0]) !== -1;
    });
    if (!entries.length) continue;
    const h = document.createElement("h2");
    h.className = "track-domain"; h.textContent = dom;
    tl.appendChild(h);
    for (const [id, t] of entries) {
      const b = document.createElement("button");
      b.className = "track-btn"; b.id = "tb-"+id;
      b.setAttribute("aria-pressed", "false");
      b.innerHTML = '<div class="tn">'+t.label+'</div><div class="td">'+t.desc+'</div>';
      b.onclick = () => selectTrack(id);
      tl.appendChild(b);
    }
  }
}

function showPicker() {
  document.body.classList.add("picker-open");
  document.getElementById("profile-picker").classList.add("open");
  var warn = document.getElementById("picker-warn"); if (warn) warn.textContent = "";
  var wrap = document.getElementById("picker-cards");
  wrap.innerHTML = "";
  if (!PROFILES.length && warn) warn.textContent = "No profiles configured (upstream/profiles.json).";
  PROFILES.forEach(function(p){
    var res = resolveTracks(p);
    var card = document.createElement("button");
    card.className = "profile-card"; card.type = "button";
    var n = res.valid.length;
    card.innerHTML = '<div class="pc-text"><div class="pc-name"></div>'
      + '<div class="pc-blurb"></div>'
      + '<div class="pc-meta">'+n+' track'+(n===1?'':'s')+'</div></div>';
    var src = _profileImgSrc(p.image);
    if (src) {  // build via DOM (not innerHTML) so the src can't break out of the attribute
      var img = document.createElement("img");
      img.className = "pc-logo"; img.src = src; img.alt = "";
      card.insertBefore(img, card.firstChild);
    }
    card.querySelector(".pc-name").textContent = p.name || p.id;
    card.querySelector(".pc-blurb").textContent = p.blurb || "";
    if (!n) {
      card.disabled = true;
      card.querySelector(".pc-meta").textContent = "no available tracks";
    } else {
      (function(pid){ card.onclick = function(){ enterProfile(pid); }; })(p.id);
      if (res.bad.length && warn) {
        warn.textContent = "Note: some configured track ids were not recognized and were skipped "
          + "(see the browser console).";
      }
    }
    wrap.appendChild(card);
  });
  var firstCard = wrap.querySelector(".profile-card:not([disabled])");
  if (firstCard) firstCard.focus();
}

function enterProfile(id) {
  var p = _profileById(id);
  if (!p) { showPicker(); return; }
  var res = resolveTracks(p);
  if (!res.valid.length) { showPicker(); return; }
  _currentProfile = id;
  try { localStorage.setItem("us:profile", id); } catch(e) {}
  buildSidebar(res.valid);
  document.getElementById("lab-name").textContent = p.name || p.id;
  var logo = document.getElementById("lab-logo");
  var lsrc = _profileImgSrc(p.image);
  if (logo) {
    if (lsrc) { logo.src = lsrc; logo.style.display = ""; }
    else { logo.removeAttribute("src"); logo.style.display = "none"; }
  }
  document.body.classList.remove("picker-open");
  document.getElementById("profile-picker").classList.remove("open");
  selectTrack(res.valid[0]);   // default to the profile's first track
}

function switchLab() { showPicker(); }

function selectTrack(id) {
  if (runId) return;
  track = id;
  document.querySelectorAll(".track-btn").forEach(b => {
    b.classList.remove("active"); b.setAttribute("aria-pressed", "false");
  });
  const activeBtn = document.getElementById("tb-"+id);
  activeBtn.classList.add("active"); activeBtn.setAttribute("aria-pressed", "true");
  document.getElementById("form-title").textContent = TRACKS[id].label;
  const fm = document.getElementById("form-msg"); if (fm) fm.textContent = "";

  const isDownload     = !!TRACKS[id].isDownload;
  const isRnaseq       = !!TRACKS[id].isRnaseq;
  const isMethylation  = !!TRACKS[id].isMethylation;
  const isGenomics     = !!TRACKS[id].isGenomics;
  const isProteomics   = !!TRACKS[id].isProteomics;

  document.getElementById("standard-grid").style.display = isDownload ? "none" : "";
  document.getElementById("run-actions").style.display   = isDownload ? "none" : "";
  document.getElementById("geo-panel").style.display     = isDownload ? "block" : "none";
  if (isDownload && typeof updateDownloadMode === "function") updateDownloadMode();

  // Reset the standard "Samples CSV" row: proteomics (and methylation-array, set later)
  // don't use a samplesheet. Resetting here also re-shows it when switching back from
  // a matrix-input track to a FASTQ track.
  const _samplesRow = (document.getElementById("inp-samples") || {}).closest
    ? document.getElementById("inp-samples").closest(".field") : null;
  if (_samplesRow) _samplesRow.style.display = isProteomics ? "none" : "";

  if (!isDownload) {
    const ef = document.getElementById("extra-fields");
    if (isRnaseq) {
      ef.innerHTML =
        '<div class="field full">' +
          '<label>Aligner / quantifier</label>' +
          '<select id="inp-aligner" onchange="updateAlignerField()">' +
            '<option value="salmon">Salmon — alignment-free (recommended)</option>' +
            '<option value="star">STAR — genome alignment + gene counts</option>' +
          '</select>' +
        '</div>' +
        '<div class="field full" id="aligner-idx-field"></div>' +
        '<div class="field full">' +
          '<label>Output format</label>' +
          '<select id="inp-format" onchange="updateFormatField()">' +
            '<option value="obama">OBAMA matrix (samples × features)</option>' +
            '<option value="matrix">DESeq2 / edgeR / limma — counts matrix + coldata</option>' +
            '<option value="both">Both</option>' +
          '</select>' +
        '</div>' +
        '<div class="field full" id="gtf-field"></div>';
      updateAlignerField();
      updateFormatField();
    } else if (isMethylation) {
      ef.innerHTML =
        '<div class="field full">' +
          '<label>Method</label>' +
          '<select id="inp-meth-method" onchange="updateMethylationMethod()">' +
            '<option value="wgbs">WGBS — whole-genome bisulfite (Bismark)</option>' +
            '<option value="array">450K / EPIC array — GEO beta matrix</option>' +
          '</select>' +
        '</div>' +
        '<div class="field full">' +
          '<label>Output format</label>' +
          '<select id="inp-format">' +
            '<option value="obama">OBAMA matrix (samples × features)</option>' +
            '<option value="matrix">limma — M-value matrix + coldata</option>' +
            '<option value="both">Both</option>' +
          '</select>' +
        '</div>' +
        '<div id="meth-extra-fields" style="display:contents"></div>';
      updateMethylationMethod();
    } else if (isProteomics) {
      ef.innerHTML =
        '<div class="field full">' +
          '<label>Input type</label>' +
          '<select id="inp-input-type">' +
            '<option value="maxquant">MaxQuant proteinGroups.txt</option>' +
            '<option value="matrix">Generic protein × sample matrix (CSV/TSV)</option>' +
          '</select>' +
        '</div>' +
        '<div class="field full">' +
          '<label>Intensity table</label>' +
          '<div class="field-row">' +
            '<input type="text" id="inp-intensities" placeholder="/data/proteinGroups.txt">' +
            '<button class="browse-btn" data-inp="inp-intensities" data-btype="file">...</button>' +
          '</div>' +
        '</div>' +
        '<div class="field full">' +
          '<label>Metadata CSV (sample, group columns)</label>' +
          '<div class="field-row">' +
            '<input type="text" id="inp-metadata-csv" placeholder="/data/metadata.csv">' +
            '<button class="browse-btn" data-inp="inp-metadata-csv" data-btype="file">...</button>' +
          '</div>' +
        '</div>' +
        '<div class="field">' +
          '<label>MaxQuant intensity column</label>' +
          '<select id="inp-intensity-col">' +
            '<option value="LFQ">LFQ — label-free (recommended)</option>' +
            '<option value="iBAQ">iBAQ</option>' +
            '<option value="Intensity">Intensity (raw)</option>' +
          '</select>' +
        '</div>' +
        '<div class="field">' +
          '<label>Normalization</label>' +
          '<select id="inp-normalize">' +
            '<option value="median">Median centering (recommended)</option>' +
            '<option value="quantile">Quantile</option>' +
            '<option value="none">None</option>' +
          '</select>' +
        '</div>' +
        '<div class="field">' +
          '<label>Min valid fraction / group</label>' +
          '<input type="number" id="inp-min-valid" value="0.5" min="0" max="1" step="0.1">' +
        '</div>' +
        '<div class="hint">Starts from a MaxQuant proteinGroups.txt (or generic matrix). ' +
          'The raw spectra → matrix step (MaxQuant/FragPipe) runs upstream.</div>';
      _wireBrowse(ef);
    } else {
      ef.innerHTML = (TRACKS[id].extra || []).map(function(k) {
        const i = EXTRA_INFO[k];
        return '<div class="field full"><label>'+i.label+'</label>' +
               '<div class="field-row">' +
               '<input type="text" id="inp-'+k+'" placeholder="'+i.ph+'">' +
               '<button class="browse-btn" data-inp="inp-'+k+'" data-btype="'+(i.btype||"dir")+'">...</button>' +
               '</div></div>';
      }).join("");
      if (id === "atacseq") {
        ef.innerHTML +=
          '<div class="field full">' +
            '<label>Output format</label>' +
            '<select id="inp-format">' +
              '<option value="obama">OBAMA matrix (samples × peaks)</option>' +
              '<option value="matrix">DESeq2 / edgeR — consensus-peak counts + coldata</option>' +
              '<option value="both">Both</option>' +
            '</select>' +
          '</div>';
      }
      if (id === "chipseq") {
        ef.innerHTML +=
          '<div class="field full">' +
            '<label>Peak type</label>' +
            '<select id="inp-peak-type">' +
              '<option value="narrow">Narrow — TFs, H3K4me3, H3K27ac</option>' +
              '<option value="broad">Broad — H3K27me3, H3K9me3, H3K36me3</option>' +
            '</select>' +
          '</div>' +
          '<div class="field full">' +
            '<label>Output format</label>' +
            '<select id="inp-format">' +
              '<option value="matrix">DESeq2 / edgeR — consensus-peak counts + coldata</option>' +
              '<option value="obama">OBAMA matrix (experimental — peak-coordinate features)</option>' +
              '<option value="both">Both</option>' +
            '</select>' +
          '</div>' +
          '<div class="hint">Samplesheet may add a <code>control</code> column naming each ' +
            'ChIP sample\\u2019s input; rows with <code>group=input</code> become the MACS2 control.</div>';
      }
      _wireBrowse(ef);
    }
  }
  if (typeof _restoreInputs === "function") _restoreInputs();
}

function _wireBrowse(el) {
  el.querySelectorAll(".browse-btn[data-inp]").forEach(function(btn) {
    btn.onclick = function() { openBrowser(btn.dataset.inp, btn.dataset.btype || "any"); };
    if (!btn.getAttribute("aria-label")) {
      btn.setAttribute("aria-label",
        btn.dataset.btype === "dir" ? "Browse for folder" : "Browse for file");
    }
  });
}

function updateAlignerField() {
  const sel = document.getElementById("inp-aligner");
  const div = document.getElementById("aligner-idx-field");
  if (!sel || !div) return;
  if (sel.value === "salmon") {
    div.innerHTML = '<label>Salmon index directory</label>' +
      '<div class="field-row">' +
        '<input type="text" id="inp-salmon-index" placeholder="'+EXTRA_INFO["salmon-index"].ph+'">' +
        '<button class="browse-btn" data-inp="inp-salmon-index" data-btype="dir">...</button>' +
      '</div>';
  } else {
    div.innerHTML = '<label>STAR index directory</label>' +
      '<div class="field-row">' +
        '<input type="text" id="inp-star-index" placeholder="'+EXTRA_INFO["star-index"].ph+'">' +
        '<button class="browse-btn" data-inp="inp-star-index" data-btype="dir">...</button>' +
      '</div>';
  }
  _wireBrowse(div);
  updateFormatField();   // GTF field is Salmon-only; refresh when aligner changes
  if (typeof _restoreInputs === "function") _restoreInputs();
}

function updateFormatField() {
  const div = document.getElementById("gtf-field");
  if (!div) return;
  const aligner = (document.getElementById("inp-aligner") || {}).value || "salmon";
  if (aligner === "salmon") {
    div.style.display = "";
    div.innerHTML =
      '<label>GTF for gene-level counts (optional)</label>' +
      '<div class="field-row">' +
        '<input type="text" id="inp-gtf" placeholder="/path/gencode.annotation.gtf.gz">' +
        '<button class="browse-btn" data-inp="inp-gtf" data-btype="file">...</button>' +
      '</div>' +
      '<div class="hint">Aggregates Salmon transcripts → genes (tximport-style). ' +
        'Without it, Salmon counts stay transcript-level. STAR is already gene-level.</div>';
    _wireBrowse(div);
  } else {
    div.style.display = "none";
    div.innerHTML = "";   // STAR output is already gene-level; no GTF needed
  }
}

function updateMethylationMethod() {
  const sel = document.getElementById("inp-meth-method");
  const div = document.getElementById("meth-extra-fields");
  if (!sel || !div) return;
  const samplesField = document.getElementById("inp-samples");
  const samplesRow   = samplesField && samplesField.closest(".field");
  if (sel.value === "array") {
    if (samplesRow) samplesRow.style.display = "none";
    div.innerHTML =
      '<div class="field full">' +
        '<label>Beta matrix CSV (probes as rows, GSM accessions as columns)</label>' +
        '<div class="field-row">' +
          '<input type="text" id="inp-betas" placeholder="/data/GSE59685_betas.csv">' +
          '<button class="browse-btn" data-inp="inp-betas" data-btype="file">...</button>' +
        '</div>' +
      '</div>' +
      '<div class="field full">' +
        '<label>Metadata CSV (geo_accession and disease.state columns)</label>' +
        '<div class="field-row">' +
          '<input type="text" id="inp-metadata-csv" placeholder="/data/metadata.csv">' +
          '<button class="browse-btn" data-inp="inp-metadata-csv" data-btype="file">...</button>' +
        '</div>' +
      '</div>';
  } else {
    if (samplesRow) samplesRow.style.display = "";
    div.innerHTML =
      '<div class="field full">' +
        '<label>Bismark genome directory</label>' +
        '<div class="field-row">' +
          '<input type="text" id="inp-bismark-genome" placeholder="'+EXTRA_INFO["bismark-genome"].ph+'">' +
          '<button class="browse-btn" data-inp="inp-bismark-genome" data-btype="dir">...</button>' +
        '</div>' +
      '</div>';
  }
  _wireBrowse(div);
  if (typeof _restoreInputs === "function") _restoreInputs();
}

// ── GEO download ──────────────────────────────────────────────────────────

var _geoRows = [];
var _geoCharKeys = [];
var _geoFilters = {};
var _geoSearch = "";
var _geoSortCol = null;
var _geoSortAsc = true;

function geoMsg(msg, color) {
  var el = document.getElementById("geo-msg");
  el.textContent = msg; el.style.color = color || "var(--red)";
}

function _geoSortedFiltered() {
  var q = _geoSearch;
  var rows = _geoRows.filter(function(row) {
    for (var col in _geoFilters) {
      var allowed = _geoFilters[col];
      if (allowed.size > 0 && !allowed.has(String(row[col] || ""))) return false;
    }
    if (q) {  // free-text search across GSM, title, and every characteristic value
      var hay = (row.gsm || "") + " " + (row.title || "");
      for (var ki = 0; ki < _geoCharKeys.length; ki++) hay += " " + (row[_geoCharKeys[ki]] || "");
      if (hay.toLowerCase().indexOf(q) === -1) return false;
    }
    return true;
  });
  if (_geoSortCol) {
    var sc = _geoSortCol, asc = _geoSortAsc;
    rows = rows.slice().sort(function(a, b) {
      var va = String(a[sc] || ""), vb = String(b[sc] || "");
      var na = parseFloat(va), nb = parseFloat(vb);
      if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
      return asc ? va.localeCompare(vb) : vb.localeCompare(va);
    });
  }
  return rows;
}

function geoSortBy(col) {
  if (_geoSortCol === col) { _geoSortAsc = !_geoSortAsc; }
  else { _geoSortCol = col; _geoSortAsc = true; }
  _buildGeoHeader();
  renderGeoTable();
}

function _buildGeoHeader() {
  var thead = document.getElementById("geo-thead");
  if (!thead) return;
  thead.innerHTML = "";
  var tr = document.createElement("tr");
  var cols = [{key:"gsm",label:"Sample (GSM)"},{key:"title",label:"Title"}];
  for (var ki = 0; ki < _geoCharKeys.length; ki++) {
    cols.push({key:_geoCharKeys[ki], label:_geoCharKeys[ki]});
  }
  cols.push({key:"_group",label:"Group"}, {key:"_srr",label:"SRR"});
  for (var ci = 0; ci < cols.length; ci++) {
    var col = cols[ci];
    var th = document.createElement("th");
    if (col.key !== "_group" && col.key !== "_srr") {
      th.className = "geo-sort";
      th.tabIndex = 0;
      th.setAttribute("role", "button");
      th.setAttribute("aria-label", "Sort by " + col.label);
      (function(k){
        th.onclick = function(){ geoSortBy(k); };
        th.onkeydown = function(e){
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); geoSortBy(k); }
        };
      })(col.key);
      var ind = (_geoSortCol === col.key)
        ? '<span class="sort-ind">'+(_geoSortAsc ? "&#9650;" : "&#9660;")+"</span>"
        : '<span class="sort-ind" style="opacity:0">&#9650;</span>';
      th.innerHTML = col.label + ind;
    } else {
      th.textContent = col.label;
    }
    tr.appendChild(th);
  }
  thead.appendChild(tr);
}

function _buildFacets() {
  var panel = document.getElementById("geo-facets");
  if (!panel) return;
  panel.innerHTML = "";
  var hasAny = false;
  for (var ki = 0; ki < _geoCharKeys.length; ki++) {
    var col = _geoCharKeys[ki];
    var seen = {}, unique = [];
    for (var ri = 0; ri < _geoRows.length; ri++) {
      var v = String(_geoRows[ri][col] || "");
      if (!seen[v]) { seen[v] = 0; }
      seen[v]++;
    }
    unique = Object.keys(seen).sort(function(a,b){
      var na=parseFloat(a),nb=parseFloat(b);
      if(!isNaN(na)&&!isNaN(nb)) return na-nb;
      return a.localeCompare(b);
    });
    if (unique.length < 2) continue;
    hasAny = true;

    var wrap = document.createElement("div");
    wrap.className = "facet-group";
    var btn = document.createElement("button");
    btn.className = "facet-btn"; btn.id = "facet-btn-"+col;
    btn.textContent = col;
    var popup = document.createElement("div");
    popup.className = "facet-popup"; popup.id = "facet-popup-"+col;
    popup.addEventListener("click", function(e){ e.stopPropagation(); });

    var ctrl = document.createElement("div");
    ctrl.className = "facet-ctrl";
    var selAllBtn = document.createElement("button");
    selAllBtn.className = "facet-ctrl-btn"; selAllBtn.textContent = "All";
    var desAllBtn = document.createElement("button");
    desAllBtn.className = "facet-ctrl-btn"; desAllBtn.textContent = "None";
    (function(c, p, u){
      selAllBtn.onclick = function(){
        p.querySelectorAll("input[type=checkbox]").forEach(function(cb){cb.checked=true;});
        delete _geoFilters[c];
        _updateFacetBtn(c, u.length, u.length);
        renderGeoTable();
      };
      desAllBtn.onclick = function(){
        p.querySelectorAll("input[type=checkbox]").forEach(function(cb){cb.checked=false;});
        _geoFilters[c] = new Set();
        _updateFacetBtn(c, 0, u.length);
        renderGeoTable();
      };
    })(col, popup, unique);
    ctrl.appendChild(selAllBtn); ctrl.appendChild(desAllBtn);
    popup.appendChild(ctrl);

    for (var ui = 0; ui < unique.length; ui++) {
      var val = unique[ui];
      var item = document.createElement("label");
      item.className = "facet-item";
      var cb = document.createElement("input");
      cb.type = "checkbox";
      cb.checked = !(_geoFilters[col] && !_geoFilters[col].has(val));
      cb.dataset.v = val;
      var span = document.createElement("span");
      span.textContent = val || "(empty)";
      var cnt = document.createElement("span");
      cnt.className = "facet-count"; cnt.textContent = seen[val];
      item.appendChild(cb); item.appendChild(span); item.appendChild(cnt);
      popup.appendChild(item);
    }

    (function(c, p, u){
      p.querySelectorAll("input[type=checkbox]").forEach(function(cb){
        cb.onchange = function(){
          var checked = Array.from(p.querySelectorAll("input[type=checkbox]:checked"))
            .map(function(x){ return x.dataset.v; });
          if (checked.length === u.length) { delete _geoFilters[c]; }
          else { _geoFilters[c] = new Set(checked); }
          _updateFacetBtn(c, checked.length, u.length);
          renderGeoTable();
        };
      });
    })(col, popup, unique);

    (function(p, b){
      b.onclick = function(e){ e.stopPropagation();
        document.querySelectorAll(".facet-popup.open").forEach(function(x){ if(x!==p) x.classList.remove("open"); });
        p.classList.toggle("open");
      };
    })(popup, btn);

    wrap.appendChild(btn); wrap.appendChild(popup);
    panel.appendChild(wrap);
  }
  panel.style.display = hasAny ? "flex" : "none";
}

function _updateFacetBtn(col, selected, total) {
  var btn = document.getElementById("facet-btn-"+col);
  if (!btn) return;
  if (selected === total) { btn.textContent = col; btn.classList.remove("active"); }
  else { btn.textContent = col+": "+selected+"/"+total+" selected"; btn.classList.add("active"); }
}

function renderGeoTable() {
  var tbody = document.getElementById("geo-tbody");
  if (!tbody) return;
  var rows = _geoSortedFiltered();
  var bc = document.getElementById("geo-bulk-count");
  if (bc) bc.textContent = rows.length + " of " + _geoRows.length + " shown";
  tbody.innerHTML = "";
  for (var i = 0; i < rows.length; i++) {
    var row = rows[i];
    var tr = document.createElement("tr");
    var tdG = document.createElement("td");
    tdG.style.cssText = "font-family:monospace;font-size:.72rem";
    tdG.textContent = row.gsm; tr.appendChild(tdG);
    var tdT = document.createElement("td");
    tdT.textContent = row.title; tr.appendChild(tdT);
    for (var ki = 0; ki < _geoCharKeys.length; ki++) {
      var td = document.createElement("td");
      td.textContent = row[_geoCharKeys[ki]] || ""; tr.appendChild(td);
    }
    var tdGrp = document.createElement("td");
    var sel = document.createElement("select");
    sel.className = "grp-sel"; sel.dataset.gsm = row.gsm;
    ["disease","control","skip"].forEach(function(v){
      var opt = document.createElement("option");
      opt.value = v; opt.textContent = v;
      if (v === row.group) opt.selected = true;
      sel.appendChild(opt);
    });
    // Capture both row ref AND this specific select element in the closure.
    // Using a bare `sel` variable (var-scoped) would read the LAST select
    // created in the loop once onchange fires — hence the second argument `s`.
    (function(r, s){ s.onchange = function(){ r.group = s.value; renderSummaryTable(); }; })(row, sel);
    tdGrp.appendChild(sel); tr.appendChild(tdGrp);
    var tdS = document.createElement("td");
    tdS.id = "srr-"+row.gsm;
    if (row.srr === null) {
      tdS.innerHTML = '<span class="spin" style="font-size:.8rem" aria-hidden="true">&#8635;</span>';
    } else if (!row.srr) {
      tdS.textContent = "not found"; tdS.className = "srr-err";
    } else {
      tdS.textContent = row.srr; tdS.className = "srr-cell";
    }
    tr.appendChild(tdS);
    tbody.appendChild(tr);
  }
  renderSummaryTable();
}

function renderSummaryTable() {
  var summaryDiv = document.getElementById("geo-summary");
  if (!summaryDiv) return;
  // Use all rows (not just filtered view) for the summary
  var selected = _geoRows.filter(function(r){ return r.group !== "skip"; });
  if (selected.length === 0) {
    summaryDiv.style.display = "none";
    var expDiv = document.getElementById("geo-csv-export");
    if (expDiv) expDiv.style.display = "none";
    return;
  }
  var disease = selected.filter(function(r){ return r.group === "disease"; }).length;
  var control = selected.filter(function(r){ return r.group === "control"; }).length;
  document.getElementById("geo-summary-title").textContent =
    selected.length + " sample" + (selected.length===1?"":"s") + " selected — " +
    disease + " disease, " + control + " control";
  var tbody = document.getElementById("geo-summary-tbody");
  tbody.innerHTML = "";
  for (var i = 0; i < selected.length; i++) {
    var row = selected[i];
    var tr = document.createElement("tr");
    var td1 = document.createElement("td");
    td1.style.cssText = "font-family:monospace;font-size:.72rem";
    td1.textContent = row.gsm; tr.appendChild(td1);
    var td2 = document.createElement("td");
    td2.textContent = row.title; tr.appendChild(td2);
    var td3 = document.createElement("td");
    td3.textContent = row.group;
    td3.className = row.group === "disease" ? "grp-disease" : "grp-control";
    tr.appendChild(td3);
    var td4 = document.createElement("td");
    if (row.srr === null) {
      td4.innerHTML = '<span class="spin" style="font-size:.8rem" aria-hidden="true">&#8635;</span>';
    } else {
      td4.textContent = row.srr || "not found";
      if (!row.srr) td4.className = "srr-err";
    }
    tr.appendChild(td4);
    tbody.appendChild(tr);
  }
  summaryDiv.style.display = "";
  _buildCSVSelector();
}

function _buildCSVSelector() {
  var exportDiv = document.getElementById("geo-csv-export");
  var container = document.getElementById("geo-csv-cols");
  if (!exportDiv || !container) return;

  // Fixed columns always available; characteristic columns added after metadata loads
  var cols = ["gsm", "title", "group", "srr"].concat(_geoCharKeys);
  var colLabels = {gsm:"Sample (GSM)", title:"Title", group:"Group", srr:"SRR"};

  // Preserve any existing checked state the user set
  var prevChecked = {};
  container.querySelectorAll("input[type=checkbox]").forEach(function(cb){
    prevChecked[cb.value] = cb.checked;
  });

  container.innerHTML = "";
  for (var i = 0; i < cols.length; i++) {
    var col = cols[i];
    var lbl = document.createElement("label");
    lbl.className = "csv-col-tag";
    var cb = document.createElement("input");
    cb.type = "checkbox";
    cb.value = col;
    // Default: checked; restore user's previous choice if it exists
    cb.checked = (col in prevChecked) ? prevChecked[col] : true;
    if (cb.checked) lbl.classList.add("checked");
    (function(c, l){ c.onchange = function(){ l.classList.toggle("checked", c.checked); }; })(cb, lbl);
    lbl.appendChild(cb);
    lbl.appendChild(document.createTextNode(" "+(colLabels[col] || col)));
    container.appendChild(lbl);
  }
  exportDiv.style.display = "";
}

function downloadSelectionCSV() {
  var selected = _geoRows.filter(function(r){ return r.group !== "skip"; });
  if (!selected.length) return;

  var cols = Array.from(
    document.querySelectorAll("#geo-csv-cols input[type=checkbox]:checked")
  ).map(function(cb){ return cb.value; });
  if (!cols.length) return;

  var colLabels = {gsm:"geo_accession", title:"title", group:"group", srr:"srr"};
  var header = cols.map(function(c){ return colLabels[c] || c; });

  var lines = [header.map(function(h){ return '"'+h.replace(/"/g,'""')+'"'; }).join(",")];
  for (var i = 0; i < selected.length; i++) {
    var row = selected[i];
    var vals = cols.map(function(c){
      var v = (c === "srr") ? (row.srr || "") : (row[c] || "");
      return '"'+String(v).replace(/"/g,'""')+'"';
    });
    lines.push(vals.join(","));
  }

  var blob = new Blob([lines.join("\\n")], {type:"text/csv"});
  var url = URL.createObjectURL(blob);
  var a = document.createElement("a");
  a.href = url; a.download = "geo_selection.csv";
  document.body.appendChild(a); a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function resetGeoGroups() {
  for (var i = 0; i < _geoRows.length; i++) { _geoRows[i].group = "skip"; }
  renderGeoTable();  // re-renders dropdowns from _geoRows (summary re-renders inside)
}

// Free-text filter over the sample table — composes with the facet filters.
function geoSearch(val) {
  _geoSearch = (val || "").trim().toLowerCase();
  renderGeoTable();
}

// Assign a group to every CURRENTLY-SHOWN row (after facets + search), so you can
// filter to a subset and label it all at once instead of row by row.
function bulkAssign(group) {
  var shown = _geoSortedFiltered();
  for (var i = 0; i < shown.length; i++) shown[i].group = group;
  renderGeoTable();
}

async function fetchGSE() {
  var inp = document.getElementById("inp-gse");
  var gse = inp.value.trim().toUpperCase();
  inp.value = gse;
  geoMsg(""); document.getElementById("geo-results").style.display = "none";
  if (!gse) { geoMsg("Enter a GSE accession."); return; }
  var btn = document.getElementById("btn-fetch");
  btn.disabled = true; btn.textContent = "Fetching…";
  _geoRows = []; _geoCharKeys = []; _geoFilters = {};
  _geoSearch = ""; _geoSortCol = null; _geoSortAsc = true;
  var searchInp = document.getElementById("geo-search");
  if (searchInp) searchInp.value = "";
  var panel = document.getElementById("geo-facets");
  if (panel) { panel.innerHTML = ""; panel.style.display = "none"; }
  var sumDiv = document.getElementById("geo-summary");
  if (sumDiv) sumDiv.style.display = "none";
  try {
    var res = await fetch("/api/geo/"+gse);
    var data = await res.json();
    if (!res.ok) { geoMsg(data.detail || "Series not found."); return; }

    document.getElementById("geo-series-info").textContent =
      data.gse+" — "+data.title+" — "+data.organism+" — "+data.n_samples+" samples";

    _geoRows = data.samples.map(function(s){
      return {gsm:s.gsm, title:s.title, srr:null, group:"skip"};
    });
    _buildGeoHeader();
    renderGeoTable();
    document.getElementById("geo-results").style.display = "";
    geoMsg("Loading metadata…", "var(--dim)");

    // Fetch characteristics in background — augments table when ready
    fetch("/api/geo/char/"+gse)
      .then(function(r){ return r.json(); })
      .then(function(chars){
        var keySet = {};
        for (var gsm in chars) {
          var c = (chars[gsm] && chars[gsm].characteristics) || {};
          for (var k in c) { keySet[k] = true; }
        }
        _geoCharKeys = Object.keys(keySet).sort();
        for (var i = 0; i < _geoRows.length; i++) {
          var row = _geoRows[i];
          var charInfo = (chars[row.gsm] && chars[row.gsm].characteristics) || {};
          for (var ki = 0; ki < _geoCharKeys.length; ki++) {
            row[_geoCharKeys[ki]] = charInfo[_geoCharKeys[ki]] || "";
          }
        }
        _buildGeoHeader();
        _buildFacets();
        renderGeoTable();  // renderGeoTable calls renderSummaryTable which calls _buildCSVSelector
      })
      .catch(function(){});  // characteristics unavailable — proceed without them

    // Fetch all SRRs for the series in one batch call (esearch → elink → esummary)
    // instead of one API call per sample — reduces N×3 NCBI calls to ~3 total.
    fetch("/api/geo/srr-all/"+gse)
      .then(function(r){ return r.json(); })
      .then(function(srrMap){
        for (var i = 0; i < _geoRows.length; i++) {
          var row = _geoRows[i];
          row.srr = srrMap[row.gsm] || "";
          var td = document.getElementById("srr-"+row.gsm);
          if (td) {
            if (row.srr) {
              td.textContent = row.srr; td.className = "srr-cell";
            } else {
              td.textContent = "not found"; td.className = "srr-err";
            }
          }
        }
        renderSummaryTable();
        geoMsg("");
      })
      .catch(function(e){ geoMsg("SRR lookup failed: "+e.message); });

  } catch(e) {
    geoMsg("Network error: "+e.message);
  } finally {
    btn.disabled = false; btn.textContent = "Fetch metadata";
  }
}

async function startGeoDownload() {
  document.getElementById("geo-dl-msg").textContent = "";
  var outdir = document.getElementById("inp-geo-outdir").value.trim() || "data/";
  var geoSamples = [];
  var visible = _geoSortedFiltered();
  for (var i = 0; i < visible.length; i++) {
    var row = visible[i];
    if (row.group === "skip") continue;
    if (row.srr === null) {
      document.getElementById("geo-dl-msg").textContent =
        "SRR not loaded for "+row.gsm+" yet — wait a moment and try again.";
      return;
    }
    if (!row.srr) {
      document.getElementById("geo-dl-msg").textContent =
        "No SRR found for "+row.gsm+". Mark it as ‘skip’ or use a different dataset.";
      return;
    }
    geoSamples.push({gsm:row.gsm, srr:row.srr, name:row.gsm, group:row.group});
  }
  var disease = geoSamples.filter(function(s){return s.group==="disease";});
  var control = geoSamples.filter(function(s){return s.group==="control";});
  if (!disease.length || !control.length) {
    document.getElementById("geo-dl-msg").textContent =
      "Assign at least one disease and one control sample.";
    return;
  }
  // Reads: "all" → full coverage (-1); a number → subsample size; blank/0 → server default.
  var rawReads = (document.getElementById("inp-geo-reads") || {value:""}).value.trim().toLowerCase();
  var reads = null;
  if (rawReads === "all" || rawReads === "full") reads = -1;
  else if (rawReads) { var n = parseInt(rawReads, 10); if (!isNaN(n)) reads = n; }
  await _kickoffRun(
    {track:"download_geo", outdir:outdir, geo_samples:geoSamples, reads:reads},
    "Download Data",
    1   // samples download in parallel — one timed batch step
  );
}

// ── GEO supplementary-matrix download (proteomics / methylation array) ──────

function updateDownloadMode() {
  var mode = (document.getElementById("inp-dl-mode") || {value:"fastq"}).value;
  var fastq = document.getElementById("dl-fastq-section");
  var matrix = document.getElementById("dl-matrix-section");
  if (fastq)  fastq.style.display  = mode === "matrix" ? "none" : "";
  if (matrix) matrix.style.display = mode === "matrix" ? "" : "none";
}

async function listSuppl() {
  var inp = document.getElementById("inp-suppl-gse");
  var gse = inp.value.trim().toUpperCase(); inp.value = gse;
  var msg = document.getElementById("suppl-msg");
  var listDiv = document.getElementById("suppl-list");
  msg.textContent = ""; msg.style.color = "var(--red)"; listDiv.style.display = "none";
  if (!gse) { msg.textContent = "Enter a GSE accession."; return; }
  var btn = document.getElementById("btn-suppl");
  btn.disabled = true; btn.textContent = "Listing…";
  try {
    var res = await fetch("/api/geo/suppl/" + encodeURIComponent(gse));
    var data = await res.json();
    if (!res.ok) { msg.textContent = data.detail || "No supplementary files found."; return; }
    var tb = document.getElementById("suppl-tbody");
    tb.innerHTML = "";
    data.files.forEach(function(f) {
      var tr = document.createElement("tr");
      var td1 = document.createElement("td");
      td1.textContent = f.name;
      td1.style.cssText = "font-family:monospace;font-size:.72rem";
      var td2 = document.createElement("td");
      var b = document.createElement("button");
      b.className = "btn ghost"; b.textContent = "↓ Download";
      b.style.cssText = "font-size:.72rem;padding:.2rem .6rem";
      b.setAttribute("aria-label", "Download " + f.name);
      (function(url, name){ b.onclick = function(){ downloadSuppl(url, name); }; })(f.url, f.name);
      td2.appendChild(b);
      tr.appendChild(td1); tr.appendChild(td2); tb.appendChild(tr);
    });
    listDiv.style.display = "";
    msg.style.color = "var(--dim)";
    msg.textContent = data.files.length + " file(s) — the matrix is usually the "
      + "*_betas / *_intensities / counts file.";
  } catch (e) {
    msg.textContent = "Network error: " + e.message;
  } finally {
    btn.disabled = false; btn.textContent = "List files";
  }
}

async function downloadSuppl(url, name) {
  var outdir = (document.getElementById("inp-suppl-outdir") || {value:"data/"}).value.trim() || "data/";
  await _kickoffRun(
    {track:"download_suppl", outdir:outdir, suppl_url:url, suppl_name:name},
    "Download matrix", 1
  );
}

async function downloadSupplUrl() {
  var url = (document.getElementById("inp-suppl-url") || {value:""}).value.trim();
  var dmsg = document.getElementById("suppl-dl-msg");
  dmsg.textContent = "";
  if (!url) { dmsg.textContent = "Paste a file URL, or use the GEO list above."; return; }
  await downloadSuppl(url, "");
}

// ── Pipeline run ──────────────────────────────────────────────────────────

function err(msg) {
  const el = document.getElementById("form-msg"); if (el) el.textContent = msg;
}

async function startRun() {
  err("");
  if (!track) { err("Select a track first."); return; }
  const outdir  = document.getElementById("inp-outdir").value.trim() || "results/";
  const threads = parseInt(document.getElementById("inp-threads").value) || 4;
  const explain = document.getElementById("inp-explain").checked;
  const dryRun  = (document.getElementById("inp-dryrun") || {checked:false}).checked;
  const body = {track:track, outdir:outdir, threads:threads, explain:explain, dry_run:dryRun};

  if (TRACKS[track].isMethylation) {
    const method = (document.getElementById("inp-meth-method") || {value:"wgbs"}).value;
    body.method = method;
    body.output_format = (document.getElementById("inp-format") || {value:"obama"}).value;
    if (method === "array") {
      const betas   = (document.getElementById("inp-betas")        || {value:""}).value.trim();
      const metaCsv = (document.getElementById("inp-metadata-csv") || {value:""}).value.trim();
      if (!betas)   { err("Beta matrix CSV path is required."); return; }
      if (!metaCsv) { err("Metadata CSV path is required."); return; }
      body.betas = betas;
      body.metadata_csv = metaCsv;
      await _kickoffRun(body, TRACKS[track].label+" (array)", 1);
    } else {
      const samples = document.getElementById("inp-samples").value.trim();
      const bg      = (document.getElementById("inp-bismark-genome") || {value:""}).value.trim();
      if (!samples) { err("Samples CSV path is required."); return; }
      if (!bg)      { err("Bismark genome directory is required."); return; }
      body.samples = samples;
      body.bismark_genome = bg;
      await _kickoffRun(body, TRACKS[track].label+" (WGBS)", 4);
    }
    return;
  }

  if (TRACKS[track].isProteomics) {
    const intensities = (document.getElementById("inp-intensities")  || {value:""}).value.trim();
    const metaCsv     = (document.getElementById("inp-metadata-csv") || {value:""}).value.trim();
    if (!intensities) { err("Intensity table path is required."); return; }
    if (!metaCsv)     { err("Metadata CSV path is required."); return; }
    body.intensities   = intensities;
    body.metadata_csv  = metaCsv;
    body.input_type    = (document.getElementById("inp-input-type")    || {value:"maxquant"}).value;
    body.intensity_col = (document.getElementById("inp-intensity-col") || {value:"LFQ"}).value;
    body.normalize     = (document.getElementById("inp-normalize")     || {value:"median"}).value;
    const mv = (document.getElementById("inp-min-valid") || {value:""}).value.trim();
    if (mv !== "") body.min_valid = parseFloat(mv);
    await _kickoffRun(body, TRACKS[track].label, 1);
    return;
  }

  const samples = document.getElementById("inp-samples").value.trim();
  if (!samples) { err("Samples CSV path is required."); return; }
  body.samples = samples;

  if (TRACKS[track].isRnaseq) {
    const aligner = (document.getElementById("inp-aligner") || {}).value || "salmon";
    body.aligner = aligner;
    const idxId = aligner==="salmon" ? "inp-salmon-index" : "inp-star-index";
    const idxEl = document.getElementById(idxId);
    const idx   = idxEl ? idxEl.value.trim() : "";
    if (!idx) { err((aligner==="salmon"?"Salmon":"STAR")+" index directory is required."); return; }
    body[aligner==="salmon" ? "salmon_index" : "star_index"] = idx;
    body.output_format = (document.getElementById("inp-format") || {value:"obama"}).value;
    if (aligner === "salmon") {
      const gtf = (document.getElementById("inp-gtf") || {value:""}).value.trim();
      if (gtf) body.gtf = gtf;
    }
  } else {
    for (const k of (TRACKS[track].extra || [])) {
      const el = document.getElementById("inp-"+k);
      const v  = el ? el.value.trim() : "";
      if (!v) { err(EXTRA_INFO[k].label+" is required."); return; }
      body[k.split("-").join("_")] = v;
    }
    if (track === "atacseq") {
      body.output_format = (document.getElementById("inp-format") || {value:"obama"}).value;
    }
    if (track === "chipseq") {
      body.output_format = (document.getElementById("inp-format") || {value:"matrix"}).value;
      body.peak_type = (document.getElementById("inp-peak-type") || {value:"narrow"}).value;
    }
  }

  await _kickoffRun(body, TRACKS[track].label, TRACKS[track].nsteps);
}

async function _kickoffRun(body, label, nsteps) {
  // Preflight: validate samplesheet, tools, and paths before starting (pipeline tracks;
  // the GEO download checks its own tools server-side). Surfaced in the setup form.
  if (body.track !== "download_geo" && body.track !== "download_suppl") {
    try {
      const cr = await (await fetch("/api/check", {
        method:"POST", headers:{"Content-Type":"application/json"}, body: JSON.stringify(body)
      })).json();
      if (cr.issues && cr.issues.length) {
        err("Preflight found " + cr.issues.length + " problem(s) — fix before running:  •  "
            + cr.issues.join("  •  "));
        return;
      }
    } catch (e) { /* if the check endpoint is unreachable, fall through and let the run report errors */ }
  }
  _runNsteps = nsteps;
  document.getElementById("setup-panel").style.display = "none";
  document.getElementById("log-panel").style.display = "flex";
  document.getElementById("log-title").textContent = label || "";
  document.getElementById("log").innerHTML = "";
  document.getElementById("log-status").innerHTML = '<span class="spin" aria-hidden="true">&#8635;</span>&nbsp;Running&hellip;';
  document.getElementById("btn-stop").style.display = "";
  document.getElementById("btn-newrun").style.display = "none";

  _clearStepTimer();
  stepIdx = 0;
  const sr = document.getElementById("steps-row");
  sr.innerHTML = "";
  for (let i = 0; i < nsteps; i++) {
    const d = document.createElement("div");
    d.className = "dot"+(i===0?" active":""); d.id = "dot-"+i;
    sr.appendChild(d);
  }
  document.getElementById("steps-legend").style.display = "";

  const res = await fetch("/api/run", {
    method:"POST", headers:{"Content-Type":"application/json"},
    body: JSON.stringify(body)
  });
  const data = await res.json();
  runId = data.run_id;

  sse = new EventSource("/api/stream/"+runId);
  sse.onmessage = function(e) {
    const d = JSON.parse(e.data);
    if (d.done) {
      sse.close(); runId = null;
      _clearStepTimer();
      const ok = d.exit_code === 0;
      document.getElementById("log-status").innerHTML = ok
        ? '<span style="color:var(--green)">✓ Complete</span>'
        : d.status==="cancelled"
          ? '<span style="color:var(--yellow)">Cancelled</span>'
          : '<span style="color:var(--red)">✗ Failed (exit '+d.exit_code+')</span> '
            + '<a href="#" onclick="reportRunError();return false;" '
            + 'style="color:var(--accent);text-decoration:none">&#9873; Report this error</a>';
      document.getElementById("btn-stop").style.display  = "none";
      document.getElementById("btn-newrun").style.display = "";
      if (ok) {
        document.querySelectorAll(".dot").forEach(function(dot) {
          dot.classList.remove("active"); dot.classList.add("done");
        });
      } else {
        document.querySelectorAll(".dot.active").forEach(function(dot) {
          dot.classList.remove("active"); dot.classList.add("err");
        });
      }
      return;
    }
    const lineDiv = appendLine(d.text, d.cls);
    if (d.cls==="step") {
      _startStepTimer(lineDiv, d.text);
      const cur = document.getElementById("dot-"+stepIdx);
      if (cur) { cur.classList.remove("active"); cur.classList.add("done"); }
      stepIdx = Math.min(stepIdx+1, _runNsteps-1);
      const nxt = document.getElementById("dot-"+stepIdx);
      if (nxt) nxt.classList.add("active");
    }
  };
  sse.onerror = function() {
    if (runId) {
      _clearStepTimer();
      appendLine("⚠ Connection lost. Check your outdir for results.", "error");
      document.getElementById("log-status").innerHTML =
        '<span style="color:var(--yellow)">Connection lost</span>';
      document.getElementById("btn-stop").style.display  = "none";
      document.getElementById("btn-newrun").style.display = "";
    }
    sse.close(); runId = null;
  };
}

// ── File browser ─────────────────────────────────────────────────────────

let _browserTarget = null;
let _browserType   = "any";
let _browserMode   = "path";   // "path" = fill an input; "sheet" = build samples.csv from a folder
let _browserOpener = null;     // element to restore focus to when the modal closes

// Make a non-button element behave like a button for keyboard users.
function _activatable(el, fn) {
  el.tabIndex = 0;
  el.setAttribute("role", "option");
  el.onclick = fn;
  el.onkeydown = function(e) {
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); fn(); }
  };
}

async function openBrowser(inputId, type) {
  _browserTarget = inputId;
  _browserType   = type;
  _browserOpener = document.activeElement;
  const el = document.getElementById(inputId);
  const start = el ? el.value.trim() : "";
  document.getElementById("browser-bar").style.display = type === "file" ? "none" : "";
  const modal = document.getElementById("browser-modal");
  modal.classList.add("open");
  await navigateBrowser(start || "~");
  // Move focus into the dialog (first folder/file row, else the Cancel button).
  const target = modal.querySelector(".browser-row") ||
                 modal.querySelector(".browser-header button");
  if (target) target.focus();
}

function buildSamplesheet() {
  _browserMode = "sheet";
  openBrowser("inp-samples", "dir");  // pick the FASTQ folder; "Select this folder" builds it
}

function closeBrowser() {
  document.getElementById("browser-modal").classList.remove("open");
  _browserMode = "path";
  if (_browserOpener && _browserOpener.focus) _browserOpener.focus();
  _browserOpener = null;
}

async function navigateBrowser(path) {
  try {
    const res  = await fetch("/api/browse?path="+encodeURIComponent(path));
    const data = await res.json();
    document.getElementById("browser-path").textContent = data.path;
    document.getElementById("browser-path").dataset.path = data.path;

    const list = document.getElementById("browser-list");
    list.innerHTML = "";

    if (data.parent) {
      const row = document.createElement("div");
      row.className = "browser-row";
      row.setAttribute("aria-label", "Parent folder");
      row.innerHTML = '<span class="browser-icon" aria-hidden="true">&#128193;</span><span style="color:var(--dim)">..</span>';
      _activatable(row, function() { navigateBrowser(data.parent); });
      list.appendChild(row);
    }

    for (const e of data.entries) {
      if (!e.is_dir && _browserType === "dir") continue;
      const row = document.createElement("div");
      row.className = "browser-row";
      const icon = e.is_dir ? "&#128193;" : "&#128196;";
      row.innerHTML = '<span class="browser-icon" aria-hidden="true">'+icon+'</span><span>'+e.name+'</span>';
      if (e.is_dir) {
        _activatable(row, function() { navigateBrowser(e.path); });
      } else {
        _activatable(row, function() { pickBrowserItem(e.path); });
      }
      list.appendChild(row);
    }
  } catch(ex) {
    document.getElementById("browser-list").innerHTML =
      '<div style="padding:.75rem;color:var(--red);font-size:.75rem">Error: '+ex.message+'</div>';
  }
}

async function selectBrowserDir() {
  const path = document.getElementById("browser-path").dataset.path;
  if (_browserMode === "sheet") {
    closeBrowser();
    if (!path) return;
    err("");
    try {
      const d = await (await fetch("/api/samplesheet?dir=" + encodeURIComponent(path))).json();
      if (d.ok) {
        const el = document.getElementById("inp-samples");
        if (el) el.value = d.path;
        err("Wrote samples.csv — " + d.n + " sample(s), " + d.paired
            + " paired. Now fill in the 'group' column (disease/control).");
      } else {
        err("No FASTQ files found in that folder.");
      }
    } catch (e) { err("Could not scan the folder."); }
    return;
  }
  if (path && _browserTarget) {
    const el = document.getElementById(_browserTarget);
    if (el) el.value = path;
  }
  closeBrowser();
}

function pickBrowserItem(path) {
  if (_browserTarget) {
    const el = document.getElementById(_browserTarget);
    if (el) el.value = path;
  }
  closeBrowser();
}

function appendLine(text, cls) {
  const log = document.getElementById("log");
  const div = document.createElement("div");
  div.className = "ll "+(cls||"out"); div.textContent = text;
  log.appendChild(div); log.scrollTop = log.scrollHeight;
  return div;
}

// ── Per-step timing ───────────────────────────────────────────────────────
// Each "step" SSE line starts a new timed segment; the previous one is frozen.
// The live clock is ground truth; "typical" is a rough keyword-based hint.
var _stepTimerInterval = null;
var _stepStartMs       = 0;
var _curStepTimeEl     = null;

function _fmtElapsed(ms) {
  var s = Math.floor(ms / 1000);
  var m = Math.floor(s / 60);
  var ss = s % 60;
  return (m < 10 ? "0" : "") + m + ":" + (ss < 10 ? "0" : "") + ss;
}

// Rough typical durations, keyed by what the step does. Assumes full-size
// input; on subsampled teaching data the live timer will read much lower.
function _typicalFor(text) {
  var t = text.toLowerCase();
  if (/downloading|prefetch|srr:/.test(t)) return "depends on run size & bandwidth";
  if (/multiqc/.test(t))              return "~10-40 s";
  if (/fastqc/.test(t))               return "~30 s-2 min";
  if (/bismark/.test(t))              return "~10-40 min";
  if (/\\bstar\\b/.test(t))           return "~3-10 min";
  if (/bowtie/.test(t))               return "~3-10 min";
  if (/salmon|quantif/.test(t))       return "~2-8 min";
  if (/macs2|peak/.test(t))           return "~1-5 min";
  if (/extract/.test(t))              return "~5-20 min";
  if (/filter/.test(t))               return "~1-4 min";
  if (/fastp|trim/.test(t))           return "~30 s-2 min";
  if (/obama|matrix/.test(t))         return "~10-60 s";
  return null;
}

function _clearStepTimer() {
  if (_stepTimerInterval) { clearInterval(_stepTimerInterval); _stepTimerInterval = null; }
  if (_curStepTimeEl) {
    // Freeze from the start timestamp (not the last tick) so it is exact.
    _curStepTimeEl.textContent = "\\u23f1 " + _fmtElapsed(Date.now() - _stepStartMs);
    _curStepTimeEl.classList.remove("running");
    _curStepTimeEl = null;
  }
}

function _startStepTimer(stepDiv, text) {
  _clearStepTimer();                       // freeze the previous step first
  _stepStartMs = Date.now();

  var timeEl = document.createElement("span");
  timeEl.className = "step-time running";
  timeEl.textContent = "\\u23f1 00:00";
  stepDiv.appendChild(timeEl);
  _curStepTimeEl = timeEl;

  var typical = _typicalFor(text);
  if (typical) {
    var t = document.createElement("div");
    t.className = "step-typical";
    t.textContent = "typical: " + typical;
    stepDiv.parentNode.insertBefore(t, stepDiv.nextSibling);
  }

  _stepTimerInterval = setInterval(function(){
    if (!_curStepTimeEl) return;
    _curStepTimeEl.textContent = "\\u23f1 " + _fmtElapsed(Date.now() - _stepStartMs);
  }, 1000);
}

async function stopRun() {
  if (!runId) return;
  await fetch("/api/cancel/"+runId, {method:"POST"});
  if (sse) sse.close();
  _clearStepTimer();
  appendLine("Run cancelled.", "error");
  document.getElementById("log-status").innerHTML =
    '<span style="color:var(--yellow)">Cancelled</span>';
  document.getElementById("btn-stop").style.display  = "none";
  document.getElementById("btn-newrun").style.display = "";
  runId = null;
}

function newRun() {
  _clearStepTimer();
  document.getElementById("log-panel").style.display = "none";
  document.getElementById("setup-panel").style.display = "";
}

// Close any open facet popup when clicking elsewhere
document.addEventListener("click", function(){
  document.querySelectorAll(".facet-popup.open").forEach(function(p){ p.classList.remove("open"); });
});

// Remember last-used paths across sessions (localStorage). Only empty fields are
// restored, so a value the user typed or the "from folder" builder set is never
// clobbered when fields are rebuilt.
function _restoreInputs() {
  document.querySelectorAll('#setup-panel input[type=text]').forEach(function(el){
    if (el.id && !el.value) {
      try { var v = localStorage.getItem("us:" + el.id); if (v) el.value = v; } catch (e) {}
    }
  });
}
document.addEventListener("change", function(e){
  var t = e.target;
  if (t && t.id && t.id.indexOf("inp-") === 0 && t.type === "text") {
    try { localStorage.setItem("us:" + t.id, t.value); } catch (e2) {}
  }
});

// Open a prefilled GitHub issue with environment diagnostics (and error output,
// when reporting a failed run).
async function reportBug(errorContext) {
  var diag = "";
  try { diag = (await (await fetch("/api/diagnostics")).json()).text; } catch (e) {}
  var body = "## What happened\\n\\n_Describe what you did and what went wrong._\\n\\n"
           + "## Environment\\n```\\n" + diag + "\\n```\\n";
  if (errorContext) body += "\\n## Error output\\n```\\n" + errorContext + "\\n```\\n";
  var url = "https://github.com/GustavoBertran/upstream/issues/new"
          + "?title=" + encodeURIComponent("[bug] ")
          + "&body=" + encodeURIComponent(body);
  window.open(url, "_blank");
}

function reportRunError() {
  var log = document.getElementById("log");
  reportBug(log ? log.textContent.slice(-2000) : "");
}

// Associate each field's <label> with its primary control so screen readers
// announce the label and clicking the label focuses the input. Fields are
// rebuilt per track, so re-run whenever the setup panel's DOM changes.
function _assocLabels() {
  document.querySelectorAll("#setup-panel .field").forEach(function(f) {
    var lbl = f.querySelector("label");
    if (!lbl || lbl.htmlFor || lbl.classList.contains("toggle")) return;
    var ctrl = f.querySelector("input, select, textarea");
    if (ctrl && ctrl.id) lbl.htmlFor = ctrl.id;
  });
}
(function() {
  var panel = document.getElementById("setup-panel");
  if (panel && typeof MutationObserver !== "undefined") {
    new MutationObserver(_assocLabels).observe(panel, {childList: true, subtree: true});
  }
})();

// File-browser modal: Escape to close, click backdrop to close, trap Tab inside.
(function() {
  var modal = document.getElementById("browser-modal");
  if (!modal) return;
  modal.addEventListener("mousedown", function(e) {
    if (e.target === modal) closeBrowser();   // click outside the box
  });
  modal.addEventListener("keydown", function(e) {
    if (!modal.classList.contains("open")) return;
    if (e.key === "Escape") { e.preventDefault(); closeBrowser(); return; }
    if (e.key !== "Tab") return;
    var f = Array.prototype.filter.call(
      modal.querySelectorAll('button, input, [tabindex]:not([tabindex="-1"])'),
      function(el) { return el.offsetParent !== null; });
    if (!f.length) return;
    var first = f[0], last = f[f.length - 1];
    if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
    else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
  });
})();

// Startup: restore the last-used lab group, else show the picker first.
(function() {
  var saved = null;
  try { saved = localStorage.getItem("us:profile"); } catch(e) {}
  if (saved && _profileById(saved) && resolveTracks(_profileById(saved)).valid.length) {
    enterProfile(saved);
  } else {
    showPicker();
  }
})();
_restoreInputs();
_assocLabels();
</script>
</body>
</html>
"""
