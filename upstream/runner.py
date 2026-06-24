import subprocess
import sys
from pathlib import Path

from rich.console import Console
from rich.rule import Rule

console = Console()

# When True, run()/pipe() print the command they WOULD run and return 0 without
# executing anything. Set by the CLI's --dry-run flag (and the web "show commands"
# toggle, via --dry-run). Lets students preview the exact tool invocations.
DRY_RUN = False


def run(cmd: list[str], cwd: Path | None = None) -> int:
    """Stream a command's stdout/stderr to the terminal. Returns the exit code."""
    display = " \\\n  ".join(str(c) for c in cmd)
    if DRY_RUN:
        console.print(f"\n[yellow][dry-run][/yellow] [dim]$ {display}[/dim]")
        return 0
    console.print(f"\n[dim]$ {display}[/dim]\n")
    result = subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd else None,
        stdout=sys.stdout,
        stderr=sys.stderr,
    )
    return result.returncode


def run_capture(cmd: list[str], out_path: Path, cwd: Path | None = None) -> int:
    """Run a command, writing its stdout to *out_path* (stderr still streams).

    Like run(), but for tools whose report goes to stdout (e.g. `bcftools stats`).
    Returns the exit code. Honors DRY_RUN.
    """
    display = " ".join(str(c) for c in cmd)
    if DRY_RUN:
        console.print(f"\n[yellow][dry-run][/yellow] [dim]$ {display} > {out_path}[/dim]")
        return 0
    console.print(f"\n[dim]$ {display} > {out_path}[/dim]\n")
    with open(out_path, "w") as fh:
        result = subprocess.run(
            [str(c) for c in cmd],
            cwd=str(cwd) if cwd else None,
            stdout=fh,
            stderr=sys.stderr,
        )
    return result.returncode


def pipe(cmd_a: list[str], cmd_b: list[str], cwd: Path | None = None) -> int:
    """Run cmd_a | cmd_b, streaming stderr of both. Returns cmd_b exit code."""
    display_a = " ".join(str(c) for c in cmd_a)
    display_b = " ".join(str(c) for c in cmd_b)
    if DRY_RUN:
        console.print(f"\n[yellow][dry-run][/yellow] [dim]$ {display_a} | {display_b}[/dim]")
        return 0
    console.print(f"\n[dim]$ {display_a} \\\n  | {display_b}[/dim]\n")
    p_a = subprocess.Popen(
        [str(c) for c in cmd_a],
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        cwd=str(cwd) if cwd else None,
    )
    p_b = subprocess.Popen(
        [str(c) for c in cmd_b],
        stdin=p_a.stdout,
        stdout=sys.stdout,
        stderr=sys.stderr,
        cwd=str(cwd) if cwd else None,
    )
    p_a.stdout.close()
    p_b.wait()
    p_a.wait()
    return p_b.returncode


def step_header(title: str, n: int, total: int) -> None:
    console.print(Rule(f"[bold cyan]{n}/{total} — {title}[/bold cyan]", style="cyan"))


def ok(message: str) -> None:
    console.print(f"[bold green]✓[/bold green] {message}")


def fail(message: str) -> None:
    console.print(f"[bold red]✗[/bold red] {message}")
