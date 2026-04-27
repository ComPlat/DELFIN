"""Programmatic Python API for DELFIN.

Wraps the CLI subcommands (``delfin/cli.py``) with typed Python functions
returning structured ``CommandResult`` objects.  Intended for embedding in
notebooks, agent tooling (MCP server), tests, and downstream automation.

The internal pattern is:

1. Build the argv list from typed kwargs.
2. Run ``cli_main(argv)`` with captured stdout/stderr.
3. Return a ``CommandResult`` containing exit code + captured streams.

This keeps the CLI as the single source of truth and avoids logic duplication.
"""
from __future__ import annotations

import contextlib
import io
import sys
from dataclasses import dataclass, field
from typing import Sequence


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class CommandResult:
    """Structured result of a DELFIN CLI subcommand.

    Attributes
    ----------
    returncode : int
        Exit code (0 = success).  ``cli_main`` returns this directly; on
        ``SystemExit`` we capture ``code`` from the exception.
    stdout : str
        Captured standard output.
    stderr : str
        Captured standard error.
    argv : list[str]
        The argv passed to ``cli_main`` (useful for debugging / logging).
    """

    returncode: int = 0
    stdout: str = ""
    stderr: str = ""
    argv: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.returncode == 0


# ---------------------------------------------------------------------------
# Internal CLI runner
# ---------------------------------------------------------------------------

def _has_flag(argv: Sequence[str], *flags: str) -> bool:
    for flag in flags:
        if flag in argv:
            return True
        prefix = f"{flag}="
        if any(arg.startswith(prefix) for arg in argv):
            return True
    return False


def _run_cli(argv: list[str], *, capture: bool = True) -> CommandResult:
    """Invoke ``delfin.cli.main`` with a built argv and capture the streams.

    Parameters
    ----------
    argv : list[str]
        The argv passed to ``cli_main`` (without the leading ``delfin``).
    capture : bool
        If True (default), capture stdout/stderr into the result.  If False,
        streams pass through to the caller's terminal.
    """
    from delfin.cli import main as cli_main

    if not capture:
        try:
            rc = cli_main(argv)
        except SystemExit as exc:
            rc = int(exc.code) if exc.code is not None else 0
        return CommandResult(returncode=rc, argv=list(argv))

    out_buf = io.StringIO()
    err_buf = io.StringIO()
    rc: int = 0
    with contextlib.redirect_stdout(out_buf), contextlib.redirect_stderr(err_buf):
        try:
            rc = cli_main(argv)
        except SystemExit as exc:
            rc = int(exc.code) if exc.code is not None else 0
    return CommandResult(
        returncode=int(rc or 0),
        stdout=out_buf.getvalue(),
        stderr=err_buf.getvalue(),
        argv=list(argv),
    )


# ---------------------------------------------------------------------------
# Pipeline (run / prepare) — kept rückwärtskompatibel
# ---------------------------------------------------------------------------

def run(
    control_file: str = "CONTROL.txt",
    *,
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str | None = None,
    extra_args: Sequence[str] | None = None,
) -> int:
    """Execute the DELFIN pipeline programmatically using CLI semantics.

    Returns the raw exit code (int) for backwards compatibility.
    Use :func:`pipeline_run` for a structured ``CommandResult``.
    """
    return pipeline_run(
        control_file=control_file,
        cleanup=cleanup,
        recalc=recalc,
        overwrite=overwrite,
        define=define,
        extra_args=extra_args,
        capture=False,
    ).returncode


def prepare(control_file: str = "CONTROL.txt", overwrite: bool = False) -> int:
    """Convenience wrapper for ``delfin --define`` to create CONTROL templates."""
    return pipeline_prepare(
        control_file=control_file, overwrite=overwrite, capture=False
    ).returncode


def pipeline_run(
    control_file: str = "CONTROL.txt",
    *,
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str | None = None,
    extra_args: Sequence[str] | None = None,
    capture: bool = True,
) -> CommandResult:
    """Execute the DELFIN pipeline and return a ``CommandResult``."""
    argv: list[str] = list(extra_args or [])

    if recalc and not _has_flag(argv, "--recalc"):
        argv.append("--recalc")
    if not cleanup and not _has_flag(argv, "--no-cleanup"):
        argv.append("--no-cleanup")
    if overwrite and not _has_flag(argv, "--overwrite"):
        argv.append("--overwrite")

    if not _has_flag(argv, "--control", "-F") and control_file != "CONTROL.txt":
        argv.extend(["--control", control_file])

    if define is not None and not _has_flag(argv, "--define", "-D"):
        argv.extend(["--define", define])

    return _run_cli(argv, capture=capture)


def pipeline_prepare(
    control_file: str = "CONTROL.txt",
    overwrite: bool = False,
    *,
    capture: bool = True,
) -> CommandResult:
    """Generate a CONTROL template via ``delfin --define <file>``."""
    argv = ["--define", control_file]
    if overwrite:
        argv.append("--overwrite")
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Tool / runtime checks (read-only, safe)
# ---------------------------------------------------------------------------

def qm_check(
    tools: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Inspect QM tool resolution (xtb, crest, etc.).  Read-only."""
    argv = ["qm_check"]
    if tools:
        argv.extend(list(tools))
    return _run_cli(argv, capture=capture)


def csp_check(*, capture: bool = True) -> CommandResult:
    """Check CSP (genarris) tool availability.  Read-only."""
    return _run_cli(["csp_check"], capture=capture)


def mlp_check(*, capture: bool = True) -> CommandResult:
    """Check MLP (torchani / AIMNet2 / MACE) backend availability.  Read-only."""
    return _run_cli(["mlp_check"], capture=capture)


def analysis_check(*, capture: bool = True) -> CommandResult:
    """Check analysis tools (Multiwfn, CENSO, ANMR, morfeus).  Read-only."""
    return _run_cli(["analysis_check"], capture=capture)


# ---------------------------------------------------------------------------
# QM tool runner
# ---------------------------------------------------------------------------

def qm_run(
    tool: str,
    tool_args: Sequence[str] | None = None,
    *,
    cwd: str = ".",
    capture_tool: bool = True,
    capture: bool = True,
) -> CommandResult:
    """Run a single QM tool via DELFIN's resolver.

    Parameters
    ----------
    tool : str
        Tool name (e.g. ``xtb``, ``crest``, ``xtb4stda``, ``stda``, ``dftb+``).
    tool_args : sequence of str
        Arguments forwarded to the tool.
    cwd : str
        Working directory for the tool process.
    capture_tool : bool
        Pass ``--capture`` so ORCA prints stdout/stderr after completion.
    capture : bool
        Capture DELFIN's own output into the ``CommandResult``.
    """
    argv = ["qm_run", tool, "--cwd", cwd]
    if capture_tool:
        argv.append("--capture")
    if tool_args:
        argv.append("--")
        argv.extend(list(tool_args))
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Process control (cleanup / stop)
# ---------------------------------------------------------------------------

def cleanup(
    *,
    orca: bool = False,
    dry_run: bool = False,
    workspace: str = ".",
    scratch: str | None = None,
    capture: bool = True,
) -> CommandResult:
    """Remove DELFIN scratch artifacts and optionally stop ORCA."""
    argv = ["cleanup", "--workspace", workspace]
    if orca:
        argv.append("--orca")
    if dry_run:
        argv.append("--dry-run")
    if scratch:
        argv.extend(["--scratch", scratch])
    return _run_cli(argv, capture=capture)


def stop(
    *,
    workspace: str = ".",
    signal_name: str = "INT",
    dry_run: bool = False,
    cleanup_after: bool = False,
    wait_seconds: float = 3.0,
    capture: bool = True,
) -> CommandResult:
    """Send a signal (INT/TERM/KILL) to running DELFIN processes.

    Parameters
    ----------
    signal_name : {"INT", "TERM", "KILL"}
        Signal to send.  Defaults to graceful SIGINT.
    cleanup_after : bool
        Run cleanup after signaling (only when not dry-run).
    """
    if signal_name not in {"INT", "TERM", "KILL"}:
        raise ValueError(
            f"signal_name must be one of INT/TERM/KILL, got {signal_name!r}"
        )
    argv = ["stop", "--workspace", workspace, "--signal", signal_name,
            "--wait-seconds", str(wait_seconds)]
    if dry_run:
        argv.append("--dry-run")
    if cleanup_after:
        argv.append("--cleanup")
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# ORCA invocation
# ---------------------------------------------------------------------------

def run_orca_input(
    input_file: str | None = None,
    output: str | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run ORCA on a ``.inp`` file (uses DELFIN's ORCA resolver)."""
    argv: list[str] = ["run_orca"]
    if input_file:
        argv.append(input_file)
    if output:
        argv.extend(["--output", output])
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# CO2 Coordinator
# ---------------------------------------------------------------------------

def co2(
    *,
    define: bool = False,
    force: bool = False,
    recalc: bool = False,
    charge: int | None = None,
    multiplicity: int | None = None,
    solvent: str | None = None,
    metal: str | None = None,
    broken_sym: str | None = None,
    capture: bool = True,
) -> CommandResult:
    """CO2 Coordinator workflow."""
    argv: list[str] = ["co2"]
    if define:
        argv.append("--define")
    if force:
        argv.append("--force")
    if recalc:
        argv.append("--recalc")
    if charge is not None:
        argv.extend(["--charge", str(charge)])
    if multiplicity is not None:
        argv.extend(["--multiplicity", str(multiplicity)])
    if solvent is not None:
        argv.extend(["--solvent", solvent])
    if metal is not None:
        argv.extend(["--metal", metal])
    if broken_sym is not None:
        argv.extend(["--broken_sym", broken_sym])
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Specialised workflows: TADF (xTB), hyperpolarisability
# ---------------------------------------------------------------------------

def tadf_xtb(
    extra_args: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run the TADF xTB workflow (``delfin tadf_xtb``)."""
    argv: list[str] = ["tadf_xtb"]
    if extra_args:
        argv.extend(list(extra_args))
    return _run_cli(argv, capture=capture)


def hyperpol(
    extra_args: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run the hyperpolarisability workflow (``delfin hyperpol``)."""
    argv: list[str] = ["hyperpol"]
    if extra_args:
        argv.extend(list(extra_args))
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Module surface
# ---------------------------------------------------------------------------

__all__ = [
    "CommandResult",
    "run",
    "prepare",
    "pipeline_run",
    "pipeline_prepare",
    "qm_check",
    "csp_check",
    "mlp_check",
    "analysis_check",
    "qm_run",
    "cleanup",
    "stop",
    "run_orca_input",
    "co2",
    "tadf_xtb",
    "hyperpol",
]
