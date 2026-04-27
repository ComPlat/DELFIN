"""MCP server for DELFIN runtime actions.

Exposes ``delfin.api`` functions as MCP tools so an agent can invoke DELFIN
workflows in a structured, typed way (no fragile ACTION-string parsing).

Tools provided
--------------
Read-only (safe):
    - ``qm_check``         — inspect QM tool resolution
    - ``csp_check``        — CSP tool availability
    - ``mlp_check``        — MLP backend availability
    - ``analysis_check``   — analysis tool availability
    - ``stop_dry_run``     — list DELFIN processes that *would* be signaled

Mutating (require explicit ``allow_mutate=True``):
    - ``cleanup``          — remove scratch artifacts (supports dry_run)
    - ``stop``             — signal DELFIN processes
    - ``pipeline_run``     — full DELFIN pipeline
    - ``pipeline_prepare`` — generate CONTROL template
    - ``run_orca_input``   — run ORCA on a .inp file
    - ``co2``              — CO2 Coordinator workflow
    - ``tadf_xtb``         — TADF xTB workflow
    - ``hyperpol``         — hyperpolarisability workflow

Mutating tools refuse to execute unless the ``allow_mutate`` parameter is
passed explicitly.  This matches the dashboard's rule that destructive
agent actions require user confirmation.

The tool functions are defined at module level so they can be imported and
tested without the optional FastMCP dependency.  ``run_server`` only loads
FastMCP when actually starting the stdio server.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any

from delfin import api as delfin_api


# ---------------------------------------------------------------------------
# Result formatting
# ---------------------------------------------------------------------------

_STDOUT_TRUNC = 16_000  # chars; large enough for normal CLI output


def _format_result(rc: delfin_api.CommandResult, *, action: str, dry_run: bool = False) -> str:
    """Format a ``CommandResult`` as a JSON string for the MCP tool reply."""
    stdout = rc.stdout or ""
    stderr = rc.stderr or ""
    if len(stdout) > _STDOUT_TRUNC:
        stdout = stdout[:_STDOUT_TRUNC] + (
            f"\n... [truncated {len(rc.stdout) - _STDOUT_TRUNC} chars]"
        )
    if len(stderr) > _STDOUT_TRUNC:
        stderr = stderr[:_STDOUT_TRUNC] + (
            f"\n... [truncated {len(rc.stderr) - _STDOUT_TRUNC} chars]"
        )
    payload: dict[str, Any] = {
        "action": action,
        "dry_run": dry_run,
        "returncode": rc.returncode,
        "ok": rc.ok,
        "stdout": stdout,
        "stderr": stderr,
        "argv": rc.argv,
    }
    return json.dumps(payload, indent=2, ensure_ascii=False)


def _refuse_mutation(action: str) -> str:
    """Standard refusal payload when ``allow_mutate`` is not set."""
    return json.dumps({
        "action": action,
        "ok": False,
        "error": "mutation_blocked",
        "message": (
            f"The tool '{action}' modifies state. Pass allow_mutate=True "
            "to execute. The agent should ask the user before doing so."
        ),
    }, indent=2, ensure_ascii=False)


# ---------------------------------------------------------------------------
# Read-only tool implementations
# ---------------------------------------------------------------------------

def tool_qm_check(tools: str = "") -> str:
    """Check QM tool resolution (xtb, crest, xtb4stda, std2, stda, dftb+).

    Args:
        tools: optional comma-separated subset (e.g. "xtb,crest"). Empty = all.
    """
    names = [t.strip() for t in tools.split(",") if t.strip()] or None
    rc = delfin_api.qm_check(tools=names)
    return _format_result(rc, action="qm_check")


def tool_csp_check() -> str:
    """Check CSP (genarris) tool availability."""
    return _format_result(delfin_api.csp_check(), action="csp_check")


def tool_mlp_check() -> str:
    """Check MLP backend availability (torchani, AIMNet2, MACE)."""
    return _format_result(delfin_api.mlp_check(), action="mlp_check")


def tool_analysis_check() -> str:
    """Check analysis tools (Multiwfn, CENSO, ANMR, morfeus)."""
    return _format_result(delfin_api.analysis_check(), action="analysis_check")


def tool_stop_dry_run(workspace: str) -> str:
    """List DELFIN processes that would be signaled (no actual signal sent)."""
    rc = delfin_api.stop(workspace=workspace, dry_run=True)
    return _format_result(rc, action="stop_dry_run", dry_run=True)


# ---------------------------------------------------------------------------
# Mutating tool implementations (require allow_mutate=True)
# ---------------------------------------------------------------------------

def tool_cleanup(
    orca: bool = False,
    dry_run: bool = True,
    workspace: str = "",
    scratch: str = "",
    allow_mutate: bool = False,
) -> str:
    """Remove DELFIN scratch artifacts.

    Defaults to dry_run=True for safety.  To actually delete, pass
    dry_run=False AND allow_mutate=True.
    """
    if not dry_run and not allow_mutate:
        return _refuse_mutation("cleanup")
    rc = delfin_api.cleanup(
        orca=orca, dry_run=dry_run, workspace=workspace,
        scratch=scratch or None,
    )
    return _format_result(rc, action="cleanup", dry_run=dry_run)


def tool_stop(
    signal_name: str = "INT",
    workspace: str = "",
    dry_run: bool = True,
    cleanup_after: bool = False,
    wait_seconds: float = 3.0,
    allow_mutate: bool = False,
) -> str:
    """Signal running DELFIN processes."""
    if not dry_run and not allow_mutate:
        return _refuse_mutation("stop")
    rc = delfin_api.stop(
        workspace=workspace, signal_name=signal_name, dry_run=dry_run,
        cleanup_after=cleanup_after, wait_seconds=wait_seconds,
    )
    return _format_result(rc, action="stop", dry_run=dry_run)


def tool_pipeline_prepare(
    control_file: str = "CONTROL.txt",
    overwrite: bool = False,
    allow_mutate: bool = False,
) -> str:
    """Generate a CONTROL.txt template (``delfin --define``)."""
    if not allow_mutate:
        return _refuse_mutation("pipeline_prepare")
    rc = delfin_api.pipeline_prepare(control_file=control_file, overwrite=overwrite)
    return _format_result(rc, action="pipeline_prepare")


def tool_pipeline_run(
    control_file: str = "CONTROL.txt",
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str = "",
    extra_args: str = "",
    allow_mutate: bool = False,
) -> str:
    """Run the full DELFIN pipeline."""
    if not allow_mutate:
        return _refuse_mutation("pipeline_run")
    extras = extra_args.split() if extra_args else None
    rc = delfin_api.pipeline_run(
        control_file=control_file, cleanup=cleanup, recalc=recalc,
        overwrite=overwrite, define=define or None, extra_args=extras,
    )
    return _format_result(rc, action="pipeline_run")


def tool_run_orca_input(
    input_file: str = "",
    output: str = "",
    allow_mutate: bool = False,
) -> str:
    """Run ORCA on a .inp file via DELFIN's ORCA resolver."""
    if not allow_mutate:
        return _refuse_mutation("run_orca_input")
    rc = delfin_api.run_orca_input(
        input_file=input_file or None, output=output or None,
    )
    return _format_result(rc, action="run_orca_input")


def tool_co2(
    define: bool = False,
    force: bool = False,
    recalc: bool = False,
    charge: int = 0,
    multiplicity: int = 0,
    solvent: str = "",
    metal: str = "",
    broken_sym: str = "",
    allow_mutate: bool = False,
) -> str:
    """CO2 Coordinator workflow."""
    if not allow_mutate:
        return _refuse_mutation("co2")
    rc = delfin_api.co2(
        define=define, force=force, recalc=recalc,
        charge=charge if charge else None,
        multiplicity=multiplicity if multiplicity else None,
        solvent=solvent or None, metal=metal or None,
        broken_sym=broken_sym or None,
    )
    return _format_result(rc, action="co2")


def tool_tadf_xtb(extra_args: str = "", allow_mutate: bool = False) -> str:
    """Run the TADF xTB workflow."""
    if not allow_mutate:
        return _refuse_mutation("tadf_xtb")
    extras = extra_args.split() if extra_args else None
    rc = delfin_api.tadf_xtb(extra_args=extras)
    return _format_result(rc, action="tadf_xtb")


def tool_hyperpol(extra_args: str = "", allow_mutate: bool = False) -> str:
    """Run the hyperpolarisability workflow."""
    if not allow_mutate:
        return _refuse_mutation("hyperpol")
    extras = extra_args.split() if extra_args else None
    rc = delfin_api.hyperpol(extra_args=extras)
    return _format_result(rc, action="hyperpol")


# ---------------------------------------------------------------------------
# Server bootstrap
# ---------------------------------------------------------------------------

def run_server(argv: list[str] | None = None) -> None:
    """Start the MCP server on stdio."""
    parser = argparse.ArgumentParser(prog="delfin-ops-server")
    parser.add_argument(
        "--workspace",
        default=os.getcwd(),
        help="Default workspace directory used when tools take a workspace arg.",
    )
    args = parser.parse_args(argv)
    default_workspace = args.workspace

    from mcp.server.fastmcp import FastMCP

    mcp = FastMCP(
        "delfin-ops",
        instructions=(
            "DELFIN operations server. Run DELFIN workflows (pipeline, "
            "ORCA, TADF, hyperpol, CO2), inspect tool availability, and "
            "control running processes. Read-only tools are always safe; "
            "mutating tools require allow_mutate=True and should be "
            "confirmed with the user first."
        ),
    )

    # Read-only — register module functions directly
    mcp.tool()(tool_qm_check)
    mcp.tool()(tool_csp_check)
    mcp.tool()(tool_mlp_check)
    mcp.tool()(tool_analysis_check)

    # stop_dry_run needs the default workspace closed over
    @mcp.tool(name="stop_dry_run", description=tool_stop_dry_run.__doc__)
    def _stop_dry_run(workspace: str = "") -> str:
        return tool_stop_dry_run(workspace or default_workspace)

    # Mutating wrappers default workspace = server cwd
    @mcp.tool(name="cleanup", description=tool_cleanup.__doc__)
    def _cleanup(
        orca: bool = False,
        dry_run: bool = True,
        workspace: str = "",
        scratch: str = "",
        allow_mutate: bool = False,
    ) -> str:
        return tool_cleanup(
            orca=orca, dry_run=dry_run,
            workspace=workspace or default_workspace,
            scratch=scratch, allow_mutate=allow_mutate,
        )

    @mcp.tool(name="stop", description=tool_stop.__doc__)
    def _stop(
        signal_name: str = "INT",
        workspace: str = "",
        dry_run: bool = True,
        cleanup_after: bool = False,
        wait_seconds: float = 3.0,
        allow_mutate: bool = False,
    ) -> str:
        return tool_stop(
            signal_name=signal_name,
            workspace=workspace or default_workspace,
            dry_run=dry_run, cleanup_after=cleanup_after,
            wait_seconds=wait_seconds, allow_mutate=allow_mutate,
        )

    mcp.tool()(tool_pipeline_prepare)
    mcp.tool()(tool_pipeline_run)
    mcp.tool()(tool_run_orca_input)
    mcp.tool()(tool_co2)
    mcp.tool()(tool_tadf_xtb)
    mcp.tool()(tool_hyperpol)

    mcp.run(transport="stdio")
