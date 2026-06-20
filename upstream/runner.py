import subprocess
import sys
from pathlib import Path

from rich.console import Console
from rich.rule import Rule

console = Console()


def run(cmd: list[str], cwd: Path | None = None) -> int:
    """Stream a command's stdout/stderr to the terminal. Returns the exit code."""
    display = " \\\n  ".join(str(c) for c in cmd)
    console.print(f"\n[dim]$ {display}[/dim]\n")
    result = subprocess.run(
        [str(c) for c in cmd],
        cwd=str(cwd) if cwd else None,
        stdout=sys.stdout,
        stderr=sys.stderr,
    )
    return result.returncode


def pipe(cmd_a: list[str], cmd_b: list[str], cwd: Path | None = None) -> int:
    """Run cmd_a | cmd_b, streaming stderr of both. Returns cmd_b exit code."""
    display_a = " ".join(str(c) for c in cmd_a)
    display_b = " ".join(str(c) for c in cmd_b)
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
