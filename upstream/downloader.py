"""
Data-location and download helpers for the upstream teaching tool.

GEO raw sequencing data is accessed via SRA: each GEO sample (GSM) maps to one
or more SRA run accessions (SRR).  Download uses the SRA Toolkit (fasterq-dump),
then optionally subsamples with seqtk to keep runtimes manageable in a lab session.

Search order when a student runs `upstream download --track <track>`:
  1. Shared server path (admin pre-staged via scripts/prepare_sample_data.sh)
  2. Local --outdir (student downloaded files themselves)
  3. Automatic download via fasterq-dump (if SRR accessions are configured)
  4. Manual instructions panel

Admin setup: fill in the 'srr' and 'geo' fields in CATALOG below, then run
scripts/prepare_sample_data.sh to stage the data for all students at once.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.table import Table

console = Console()

_SHARED_DEFAULT = Path("/data/upstream/shared")
_SEED = 42
_DEFAULT_NREADS = 1_000_000


# ── Catalog ───────────────────────────────────────────────────────────────
#
# HOW TO CONFIGURE:
#   For each sample, set 'srr' to the SRA run accession (e.g. "SRR1234567").
#   Find it on the GEO sample page (GSM…) → SRA Experiments → SRR link.
#   Set 'geo' to the GEO series accession (GSE…) for display purposes.
#
#   Leave srr=None for any sample not yet assigned; those samples will fall
#   through to the manual instructions panel.

CATALOG: dict[str, dict] = {
    "rnaseq": {
        "label": "RNA-seq — breast cancer",
        "geo": None,          # e.g. "GSE12345"
        "samples": [
            {"name": "disease", "group": "disease", "srr": None},
            {"name": "control", "group": "control", "srr": None},
        ],
    },
    "atacseq": {
        "label": "ATAC-seq — breast cancer",
        "geo": None,
        "samples": [
            {"name": "disease", "group": "disease", "srr": None},
            {"name": "control", "group": "control", "srr": None},
        ],
    },
    "methylation": {
        "label": "Methylation (WGBS/RRBS) — breast cancer",
        "geo": None,
        "samples": [
            {"name": "disease", "group": "disease", "srr": None},
            {"name": "control", "group": "control", "srr": None},
        ],
    },
    "qc": {
        "label": "QC — reuses rnaseq FASTQ",
        "geo": None,
        "samples": None,   # delegate to rnaseq at runtime
    },
}


# ── Public API ────────────────────────────────────────────────────────────


def print_catalog() -> None:
    """Print a summary table of configured datasets."""
    table = Table(
        title="Available Example Datasets",
        show_header=True,
        header_style="bold cyan",
        show_lines=True,
    )
    table.add_column("Track", style="bold", no_wrap=True)
    table.add_column("Dataset")
    table.add_column("GEO Series")
    table.add_column("SRR accessions")
    for track, info in CATALOG.items():
        if info["samples"] is None:
            srr_col = "→ uses rnaseq data"
        else:
            srrs = [s["srr"] or "[dim]not set[/dim]" for s in info["samples"]]
            srr_col = ", ".join(srrs)
        geo_col = info["geo"] or "[dim]not set[/dim]"
        table.add_row(track, info["label"], geo_col, srr_col)
    console.print()
    console.print(table)
    console.print(
        "\nRun [bold]upstream download --track <track> --outdir <dir>[/bold] to fetch data.\n"
        "Use [bold]--instructions[/bold] to see manual steps for any track.\n"
        "Admins: fill in the 'srr' and 'geo' fields in upstream/downloader.py CATALOG,\n"
        "then run scripts/prepare_sample_data.sh to stage shared data for all students.\n"
    )


def locate_or_download(
    track: str,
    outdir: Path,
    shared_data: Path = _SHARED_DEFAULT,
    instructions_only: bool = False,
    nreads: int = _DEFAULT_NREADS,
    subsample: bool = True,
    samples_csv: Optional[Path] = None,
) -> None:
    """
    Find or fetch example data for *track*, then write samples.csv into *outdir*.

    Search order:
    1. Shared server path (admin pre-staged).
    2. Local *outdir* (student already has files there).
    3. Download via fasterq-dump (if SRR accessions are configured in CATALOG).
    4. Print manual instructions and exit.
    """
    if track == "qc":
        console.print("[dim]QC track reuses rnaseq FASTQ data.[/dim]")
        track = "rnaseq"

    if track not in CATALOG:
        console.print(
            f"[bold red]Error:[/bold red] Unknown track '{track}'. "
            f"Choose from: {', '.join(CATALOG)}"
        )
        raise SystemExit(1)

    info = CATALOG[track]

    if instructions_only:
        _show_instructions(track, info)
        return

    if info["samples"] is None:
        _show_instructions(track, info)
        return

    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = samples_csv or (outdir / "samples.csv")

    # 1. Shared server data (pre-staged by admin)
    shared_track = shared_data / "fastq" / track
    expected = _expected_filenames(info["samples"])
    if _files_present(expected, shared_track):
        _write_csv(info["samples"], expected, shared_track, csv_path, track)
        return

    # 2. Local outdir (student put files there themselves)
    if _files_present(expected, outdir):
        _write_csv(info["samples"], expected, outdir, csv_path, track)
        return

    # 3. Download via SRA toolkit
    configured = [s for s in info["samples"] if s["srr"]]
    if len(configured) == len(info["samples"]):
        console.print(
            f"\n[bold]Downloading {info['label']}[/bold]"
            + (f"  (GEO: {info['geo']})" if info["geo"] else "")
            + "\n"
        )
        _download_all(info["samples"], outdir, nreads, subsample)
        _write_csv(info["samples"], expected, outdir, csv_path, track)
        return

    # 4. No data found, SRR not configured
    missing_srr = [s["name"] for s in info["samples"] if not s["srr"]]
    console.print(
        f"\n[yellow]No data found and SRR accessions not configured.[/yellow]\n"
        f"  Samples without SRR: {', '.join(missing_srr)}\n"
        f"  Checked shared path: [dim]{shared_track}[/dim]\n"
        f"  Checked local path:  [dim]{outdir}[/dim]\n"
    )
    _show_instructions(track, info)


# ── Download helpers ──────────────────────────────────────────────────────


def _require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        console.print(
            f"[bold red]Error:[/bold red] '{name}' not found on PATH. "
            f"Install the SRA Toolkit (conda install -c bioconda sra-tools) and retry."
        )
        raise SystemExit(1)
    return path


def _download_all(
    samples: list[dict],
    outdir: Path,
    nreads: int,
    subsample: bool,
) -> None:
    _require_tool("fasterq-dump")
    if subsample:
        _require_tool("seqtk")

    for s in samples:
        srr = s["srr"]
        name = s["name"]
        r1_final = outdir / f"{name}_R1.fastq.gz"
        r2_final = outdir / f"{name}_R2.fastq.gz"

        if r1_final.exists() and r2_final.exists():
            console.print(f"[dim]  {name}: FASTQ already present, skipping download.[/dim]")
            continue

        console.print(f"\n[bold cyan]{name}[/bold cyan] ({s['group']})  SRR: {srr}")
        _fasterq_dump(srr, name, outdir, r1_final, r2_final, nreads, subsample)
        console.print(f"  [green]✓[/green] {name} ready")


def _fasterq_dump(
    srr: str,
    name: str,
    outdir: Path,
    r1_final: Path,
    r2_final: Path,
    nreads: int,
    subsample: bool,
) -> None:
    """Run fasterq-dump for *srr*, producing gzipped paired FASTQ in *outdir*."""
    tmp_r1 = outdir / f"{srr}_1.fastq"
    tmp_r2 = outdir / f"{srr}_2.fastq"

    if not (tmp_r1.exists() and tmp_r2.exists()):
        console.print(f"  Fetching {srr} via fasterq-dump…")
        rc = subprocess.run([
            "fasterq-dump",
            "--split-files",
            "--outdir", str(outdir),
            "--progress",
            srr,
        ]).returncode
        if rc != 0:
            console.print(f"[bold red]Error:[/bold red] fasterq-dump failed for {srr}.")
            raise SystemExit(1)

    if subsample:
        console.print(f"  Subsampling to {nreads:,} reads…")
        _seqtk_sample(tmp_r1, r1_final, nreads)
        _seqtk_sample(tmp_r2, r2_final, nreads)
        tmp_r1.unlink(missing_ok=True)
        tmp_r2.unlink(missing_ok=True)
    else:
        console.print(f"  Compressing…")
        _gzip_file(tmp_r1, r1_final)
        _gzip_file(tmp_r2, r2_final)


def _seqtk_sample(src: Path, dest: Path, nreads: int) -> None:
    with dest.open("wb") as fh_out:
        p_seqtk = subprocess.Popen(
            ["seqtk", "sample", "-s", str(_SEED), str(src), str(nreads)],
            stdout=subprocess.PIPE,
        )
        p_gzip = subprocess.Popen(["gzip"], stdin=p_seqtk.stdout, stdout=fh_out)
        p_seqtk.stdout.close()
        p_gzip.wait()
        p_seqtk.wait()
    if p_seqtk.returncode != 0 or p_gzip.returncode != 0:
        dest.unlink(missing_ok=True)
        console.print(f"[bold red]Error:[/bold red] seqtk/gzip failed for {src.name}.")
        raise SystemExit(1)


def _gzip_file(src: Path, dest: Path) -> None:
    with dest.open("wb") as fh_out:
        subprocess.run(["gzip", "-c", str(src)], stdout=fh_out, check=True)
    src.unlink(missing_ok=True)


# ── CSV / detection helpers ───────────────────────────────────────────────


def _expected_filenames(samples: list[dict]) -> list[tuple[str, str]]:
    """Return (r1_filename, r2_filename) pairs matching what we write."""
    return [(f"{s['name']}_R1.fastq.gz", f"{s['name']}_R2.fastq.gz") for s in samples]


def _files_present(expected: list[tuple[str, str]], base: Path) -> bool:
    return all((base / r1).exists() and (base / r2).exists() for r1, r2 in expected)


def _write_csv(
    samples: list[dict],
    expected: list[tuple[str, str]],
    base: Path,
    csv_path: Path,
    track: str,
) -> None:
    rows = [
        {
            "name": s["name"],
            "group": s["group"],
            "r1": str((base / r1).resolve()),
            "r2": str((base / r2).resolve()),
        }
        for s, (r1, r2) in zip(samples, expected)
    ]
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["name", "group", "r1", "r2"])
        writer.writeheader()
        writer.writerows(rows)

    console.print(f"\n[bold green]✓ Samplesheet written:[/bold green] {csv_path}")
    console.print(f"  {len(rows)} sample(s) — data in [dim]{base}[/dim]\n")
    console.print(
        f"Run the pipeline with:\n"
        f"  [bold]upstream {track} --samples {csv_path} --outdir results/[/bold]\n"
    )


def _show_instructions(track: str, info: dict) -> None:
    srr_lines = ""
    if info.get("samples"):
        for s in info["samples"]:
            srr = s["srr"] or "SRR_NOT_CONFIGURED"
            srr_lines += f"  # {s['name']} ({s['group']}): {srr}\n"

    geo_ref = f"GEO series: {info['geo']}\n\n" if info.get("geo") else ""

    manual = (
        f"{geo_ref}"
        f"Download each SRR run from GEO via the SRA Toolkit:\n\n"
        f"```bash\n"
        f"# Install SRA Toolkit if needed:\n"
        f"#   conda install -c bioconda sra-tools\n\n"
        f"# Download paired FASTQ for each sample:\n"
        f"{srr_lines}"
        f"fasterq-dump --split-files SRR_ACCESSION --outdir ./fastq/ --progress\n\n"
        f"# Subsample to 1 M reads (keeps lab-session runtimes short):\n"
        f"seqtk sample -s 42 SRR_ACCESSION_1.fastq 1000000 | gzip > sample_R1.fastq.gz\n"
        f"seqtk sample -s 42 SRR_ACCESSION_2.fastq 1000000 | gzip > sample_R2.fastq.gz\n"
        f"```\n\n"
        f"Then place files in [dim]{_SHARED_DEFAULT}/fastq/{track}/[/dim] (admin) or run:\n"
        f"  upstream download --track {track} --outdir <dir>\n"
        f"to write samples.csv once the files are in place."
    )
    console.print(
        Panel(
            Markdown(manual),
            title=f"[bold]How to obtain {info['label']}[/bold]",
            border_style="yellow",
            padding=(1, 2),
        )
    )
