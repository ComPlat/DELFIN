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
# On-demand operational-pattern lookup
# ---------------------------------------------------------------------------
#
# Concrete slash-chain recipes for common dashboard workflows.  These
# used to be baked into ``dashboard_agent.md`` but that paid token cost
# every turn — even for conversations that have nothing to do with
# batch jobs or recalc.  Now they live behind a typed lookup the agent
# only calls when it actually needs a workflow recipe.
#
# Adding a new pattern: drop an entry into ``_DASHBOARD_PATTERNS`` and
# add the name to the docstring for discoverability.

_DASHBOARD_PATTERNS: dict[str, str] = {
    "batch": (
        "## Batch jobs (Submit-Tab)\n\n"
        "- Build a batch from EVERY calc-folder's initial geometry:\n"
        "    `ACTION: /batch from-calc`\n"
        "- Build from a glob subset (e.g. all Casagrande folders):\n"
        "    `ACTION: /batch from-calc Casagrande*`\n"
        "- Add one SMILES line manually:\n"
        "    `ACTION: /batch add Name;SMILES;charge=…`\n"
        "- Show / clear the current batch text:\n"
        "    `ACTION: /batch show`  ·  `ACTION: /batch clear`\n\n"
        "**Never** assemble batch text by reading XYZ files yourself —\n"
        "`/batch from-calc` already collects every initial.xyz (or\n"
        "fallback input.txt / coords.xyz) across `calculations/` and\n"
        "writes a properly formatted block into the Submit-Tab textarea.\n\n"
        "Example:\n"
        "  > User: \"Bau einen Batch aus allen XYZ in calc/.\"\n"
        "  > Agent:\n"
        "  >\n"
        "  >     ACTION: /batch from-calc\n"
    ),
    "control_edit": (
        "## CONTROL.txt edits (Submit-Tab)\n\n"
        "- One key change → `ACTION: /control key <key> <value>`\n"
        "- Replace whole content (rare) → `ACTION: /control set <multi-line>`\n"
        "- Validate before submit → `ACTION: /control validate`\n\n"
        "Never paste the full CONTROL into chat — use `/control key` per change.\n"
    ),
    "smart_recalc": (
        "## Smart Recalc (Calc browser)\n\n"
        "Selecting Smart Recalc auto-loads CONTROL into the editor, so\n"
        "the right idiom is `/ui calc-editor replace <old> <new>`:\n\n"
        "  > User: \"smart recalc von Foo/bar mit PAL=1, 1 h\"\n"
        "  > Agent:\n"
        "  >\n"
        "  >     ACTION: /calc cd Foo/bar\n"
        "  >     ACTION: /calc select CONTROL.txt\n"
        "  >     ACTION: /ui calc-options value Smart Recalc\n"
        "  >     ACTION: /ui calc-editor replace PAL=40 PAL=1\n"
        "  >     ACTION: /ui calc-override-time value 01:00:00\n"
        "  >\n"
        "  > Then ASK before submitting:\n"
        "  >\n"
        "  >     ACTION: /ui calc-override-btn click\n\n"
        "`calc-override-btn` is the submit; `calc-override-time` is the\n"
        "time field. Don't confuse with `calc-recalc-btn` /\n"
        "`calc-submit-recalc-btn` — those are NOT the Smart-Recalc panel.\n"
    ),
    "submit_orca": (
        "## Submit a single ORCA job (ORCA Builder)\n\n"
        "- Set fields → `ACTION: /ui orca-method value PBE0`,\n"
        "  `/ui orca-basis value def2-TZVP`, `/ui orca-charge value 0`, …\n"
        "- Switch to the tab → `ACTION: /tab orca`\n"
        "- Submit (after the user explicitly OKs!) → `ACTION: /orca submit`\n"
    ),
    "analyze": (
        "## Analyze existing calculations\n\n"
        "- One folder, full → `ACTION: /analyze <dir>`\n"
        "- Just energies → `ACTION: /analyze energy <dir>`\n"
        "- SCF convergence → `ACTION: /analyze convergence <dir>`\n"
        "- Error scan → `ACTION: /analyze errors <dir>`\n"
        "- All folders overview → `ACTION: /analyze status`\n\n"
        "For multi-folder energy tables prefer\n"
        "`mcp__delfin-ops__extract_energy_table(folders=[...])` — it\n"
        "returns structured data the agent can format directly.\n"
    ),
    "recalc": (
        "## Recalc check / submit\n\n"
        "- Check one folder (safe, read-only) → `ACTION: /recalc check <dir>`\n"
        "- Scan everything (safe) → `ACTION: /recalc check-all`\n"
        "- Submit recalc (DESTRUCTIVE — needs explicit user OK):\n"
        "    `ACTION: /recalc <dir>` and confirm\n"
        "- Bulk auto-recalc only on explicit \"alle neuberechnen\":\n"
        "    `ACTION: /recalc auto`\n"
    ),
    "cancel": (
        "## Cancel jobs\n\n"
        "- One job → `ACTION: /cancel <job_id>` (after user OK)\n"
        "- All — only on explicit \"cancel all\":\n"
        "    `ACTION: /cancel all`\n"
    ),
}


def list_dashboard_patterns() -> list[str]:
    """Return the names of every operational pattern available for lookup.

    Use this when the agent isn't sure which pattern matches the user's
    request — a names-only list is cheap and lets the agent pick.
    """
    return sorted(_DASHBOARD_PATTERNS)


def get_dashboard_pattern(name: str) -> str:
    """Return the slash-chain recipe for a named dashboard workflow.

    Available names: ``batch``, ``control_edit``, ``smart_recalc``,
    ``submit_orca``, ``analyze``, ``recalc``, ``cancel``.

    Names are matched case-insensitively. An unknown name returns a
    short error string with the available choices so the agent can
    self-correct.
    """
    if not name:
        return (
            "No pattern requested. Call with one of: "
            + ", ".join(list_dashboard_patterns())
        )
    key = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    text = _DASHBOARD_PATTERNS.get(key)
    if text is None:
        return (
            f"Unknown pattern: {name!r}. Available: "
            + ", ".join(list_dashboard_patterns())
        )
    return text


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
    "list_dashboard_patterns",
    "get_dashboard_pattern",
]
