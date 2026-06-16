"""MCP server exposing the DELFIN tools platform.

Lets the agent system, an IDE, or any MCP client bind to DELFIN's building
blocks and applications over a standard protocol — discover capabilities and
applications, read the central key vocabulary, validate and run workflows, and
probe / plan tool installation — without any Python coupling.

The request handlers are plain module-level functions returning JSON strings, so
they are testable without the optional ``mcp`` dependency; ``run_server`` (which
imports ``mcp``) only wires them to the stdio transport, mirroring
``delfin/doc_server/server.py``.
"""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, Optional

from delfin.tools import platform
from delfin.tools._application import get_application
from delfin.tools._keys import get_key, list_keys
from delfin.tools.manifest import manifest_json


def _dumps(obj: Any) -> str:
    return json.dumps(obj, indent=2, ensure_ascii=False, default=str)


# --- capabilities ---------------------------------------------------------


def h_get_manifest() -> str:
    return manifest_json()


def h_list_capabilities() -> str:
    return _dumps(platform.list_capabilities())


def h_describe_capability(name: str) -> str:
    from delfin.tools.manifest import contract_to_dict
    c = platform.describe_capability(name)
    if c is None:
        return _dumps({"error": f"unknown capability {name!r}"})
    return _dumps(contract_to_dict(c))


def h_catalog(by: str = "category") -> str:
    grouped = platform.catalog(by=by)
    return _dumps({k: [c.name for c in v] for k, v in grouped.items()})


# --- keys -----------------------------------------------------------------


def h_list_keys() -> str:
    from delfin.tools.manifest import _keyspec_to_dict
    return _dumps([_keyspec_to_dict(ks) for _, ks in sorted(list_keys().items())])


def h_describe_key(name: str) -> str:
    from delfin.tools.manifest import _keyspec_to_dict
    ks = get_key(name)
    if ks is None:
        return _dumps({"error": f"unknown key {name!r}"})
    return _dumps(_keyspec_to_dict(ks))


# --- applications ---------------------------------------------------------


def h_list_applications() -> str:
    return _dumps(platform.list_applications())


def h_describe_application(name: str) -> str:
    app = get_application(name)
    if app is None:
        return _dumps({"error": f"unknown application {name!r}"})
    return _dumps(app.to_dict())


def h_validate_application(
    name: str, inputs: Optional[Dict[str, Any]] = None, geometry: bool = False,
) -> str:
    inputs = inputs or {}
    rep = platform.validate_application(name, geometry=geometry, **inputs)
    if rep is None:
        return _dumps({"error": f"unknown application {name!r}"})
    return _dumps({
        "ok": rep.ok,
        "diagnostics": [
            {
                "step": d.step_name,
                "label": d.label,
                "level": d.level.value,
                "location": d.location,
                "missing_params": list(d.missing_params),
                "missing_inputs": list(d.missing_inputs),
                "messages": list(d.messages),
            }
            for d in rep.diagnostics
        ],
    })


def h_run_application(
    name: str, inputs: Optional[Dict[str, Any]] = None, cores: int = 1,
) -> str:
    inputs = inputs or {}
    res = platform.run_application(name, cores=cores, **inputs)
    return _dumps({
        "name": res.name,
        "ok": res.ok,
        "outputs": res.outputs,
        "error": res.error,
    })


# --- runs (async execution) -----------------------------------------------


def h_submit_application(
    name: str, inputs: Optional[Dict[str, Any]] = None, cores: int = 1,
) -> str:
    inputs = inputs or {}
    run_id = platform.submit_application(name, cores=cores, **inputs)
    return _dumps({"run_id": run_id})


def h_run_status(run_id: str) -> str:
    from dataclasses import asdict
    rec = platform.run_record(run_id)
    if rec is None:
        return _dumps({"error": f"unknown run {run_id!r}"})
    return _dumps(asdict(rec))


def h_list_runs() -> str:
    out = [
        {"id": r.id, "name": r.name, "status": r.status,
         "created_at": r.created_at, "outputs": r.outputs}
        for r in platform.list_runs()
    ]
    return _dumps(out)


def h_cancel_run(run_id: str) -> str:
    return _dumps({"cancelled": platform.cancel_run(run_id)})


# --- environment ----------------------------------------------------------


def h_probe() -> str:
    out = []
    for cr in platform.probe():
        out.append({
            "capability": cr.step_name,
            "ready": cr.ready,
            "missing": [
                {
                    "name": r.name, "kind": r.kind, "policy": r.policy,
                    "source": r.source, "hint": r.install_hint,
                }
                for r in cr.missing
            ],
        })
    return _dumps(out)


def h_install_plan() -> str:
    return _dumps(platform.install_plan())


# --- server wiring --------------------------------------------------------


def run_server(argv: Optional[list[str]] = None) -> None:
    """Start the MCP server on stdio."""
    argparse.ArgumentParser(prog="delfin-tools-server").parse_args(argv)

    from mcp.server.fastmcp import FastMCP

    mcp = FastMCP(
        "delfin-tools",
        instructions=(
            "DELFIN tools platform. Discover computational-chemistry building "
            "blocks (capabilities) and applications (named workflows with typed "
            "inputs/outputs), read the central key vocabulary (functional, basis, "
            "solvent, … with allowed values), validate and run workflows, and "
            "probe / plan tool installation. License-restricted tools (ORCA, "
            "Turbomole) are never auto-installed — guidance is returned instead."
        ),
    )

    @mcp.tool()
    def get_manifest() -> str:
        """Full machine-readable manifest: capabilities, applications, keys, schemas."""
        return h_get_manifest()

    @mcp.tool()
    def list_capabilities() -> str:
        """List all registered building blocks (tool names)."""
        return h_list_capabilities()

    @mcp.tool()
    def describe_capability(name: str) -> str:
        """Full contract (params, ports, requirements) for one capability."""
        return h_describe_capability(name)

    @mcp.tool()
    def catalog(by: str = "category") -> str:
        """Capabilities grouped by 'category', 'produces', or 'consumes'."""
        return h_catalog(by)

    @mcp.tool()
    def list_keys() -> str:
        """List the central key vocabulary (functional, basis, solvent, …) with allowed values."""
        return h_list_keys()

    @mcp.tool()
    def describe_key(name: str) -> str:
        """Type, default and allowed values for one well-known key (e.g. 'functional')."""
        return h_describe_key(name)

    @mcp.tool()
    def list_applications() -> str:
        """List all registered applications (named workflows)."""
        return h_list_applications()

    @mcp.tool()
    def describe_application(name: str) -> str:
        """An application's full definition: inputs, outputs, and pipeline spec."""
        return h_describe_application(name)

    @mcp.tool()
    def validate_application(name: str, inputs: dict | None = None,
                             geometry: bool = False) -> str:
        """Statically validate an application's inputs + pipeline (no execution)."""
        return h_validate_application(name, inputs, geometry)

    @mcp.tool()
    def run_application(name: str, inputs: dict | None = None, cores: int = 1) -> str:
        """Run a named application with named inputs; returns named outputs.

        Note: this can be long-running (it executes real QM jobs).
        """
        return h_run_application(name, inputs, cores)

    @mcp.tool()
    def submit_application(name: str, inputs: dict | None = None, cores: int = 1) -> str:
        """Submit an application run in the background; returns a run_id immediately."""
        return h_submit_application(name, inputs, cores)

    @mcp.tool()
    def run_status(run_id: str) -> str:
        """Status, outputs, events and metrics of a submitted run."""
        return h_run_status(run_id)

    @mcp.tool()
    def list_runs() -> str:
        """List submitted runs (newest first) with status and outputs."""
        return h_list_runs()

    @mcp.tool()
    def cancel_run(run_id: str) -> str:
        """Request cooperative cancellation of a running run."""
        return h_cancel_run(run_id)

    @mcp.tool()
    def probe() -> str:
        """Readiness of every capability against its required tools, with install hints."""
        return h_probe()

    @mcp.tool()
    def install_plan() -> str:
        """Which missing tools are auto-installable vs. manual (license-restricted)."""
        return h_install_plan()

    mcp.run(transport="stdio")


def main(argv: Optional[list[str]] = None) -> int:
    run_server(argv)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
