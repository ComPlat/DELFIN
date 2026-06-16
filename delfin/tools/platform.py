"""DELFIN tools platform — the stable facade other systems bind to.

One documented entry point over the building-block framework, so the agent
system, a UI, or another tool can drive it without reaching into internals:

* **capabilities** — discover and describe the registered tools (building blocks)
* **applications** — register / discover / validate / run named workflows that
  carry a typed input/output contract
* **environment** — probe which tools are installed and install the open-source
  ones (license-restricted tools get guidance, never auto-installed)

    from delfin.tools import platform

    platform.list_capabilities()
    platform.describe_application("opt_freq_energy")
    res = platform.run_application("opt_freq_energy", smiles="CCO", charge=0, cores=4)
    res.outputs                      # {"energy_Eh": ..., "gibbs_Eh": ...}

    platform.probe()                 # readiness of every capability
    platform.install_tools()         # plan (run=True to install open-source tools)

Execution is currently synchronous; this facade is the seam a future async
run-store / event stream plugs in behind without changing callers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.tools import _environment
from delfin.tools._application import (
    Application,
    extract_outputs,
    get_application,
)
from delfin.tools._application import list_applications as _list_applications
from delfin.tools._application import register_application as _register_application
from delfin.tools._catalog import catalog as _catalog
from delfin.tools._catalog import compatible_successors as _compatible_successors
from delfin.tools._catalog import describe as _describe
from delfin.tools._registry import list_steps


# ======================================================================
#  Capabilities (building blocks)
# ======================================================================


def list_capabilities() -> List[str]:
    """Names of all registered building blocks (tools)."""
    return sorted(list_steps().keys())


def describe_capability(name: str):
    """Full :class:`StepContract` for one capability, or ``None`` if unknown."""
    return _describe(name)


def catalog(*, by: str = "category") -> Dict[str, list]:
    """Grouped view of all capabilities (by category / produces / consumes)."""
    return _catalog(by=by)


def compatible_successors(name: str) -> List[str]:
    """Capabilities whose inputs *name* can satisfy (what can follow it)."""
    return _compatible_successors(name)


# ======================================================================
#  Applications (workflows with a contract)
# ======================================================================


@dataclass
class ApplicationResult:
    """Outcome of running an application through the platform."""

    name: str
    ok: bool
    outputs: Dict[str, Any] = field(default_factory=dict)
    validation: Any = None            # ValidationReport (or None)
    pipeline_result: Any = None       # PipelineResult (or None)
    error: Optional[str] = None


def register_application(app: Application) -> Application:
    """Register an application so it can be discovered and run by name."""
    return _register_application(app)


def list_applications() -> List[str]:
    """Names of all registered applications."""
    return sorted(_list_applications().keys())


def describe_application(name: str) -> Optional[Application]:
    """Return the :class:`Application` (its contract is the description)."""
    return get_application(name)


def validate_application(name: str, *, geometry: bool = False, **inputs: Any):
    """Validate an application's inputs and its materialised pipeline statically.

    Returns a :class:`ValidationReport`, or ``None`` if the application is unknown.
    Missing required inputs surface as unresolved-placeholder warnings/errors in
    the report (and via :meth:`Application.missing_inputs`).
    """
    app = get_application(name)
    if app is None:
        return None
    return app.validate(geometry=geometry, **inputs)


def run_application(
    name: str,
    *,
    cores: int = 1,
    geometry: Optional[str | Path] = None,
    work_dir: Optional[Path] = None,
    validate: bool = True,
    **inputs: Any,
) -> ApplicationResult:
    """Run a named application with named inputs and return its named outputs.

    Fails fast (without executing) on an unknown application, missing required
    inputs, or — when ``validate=True`` — a pipeline that does not statically
    validate.
    """
    app = get_application(name)
    if app is None:
        return ApplicationResult(name=name, ok=False, error=f"unknown application {name!r}")

    missing = app.missing_inputs(inputs)
    if missing:
        return ApplicationResult(
            name=name, ok=False,
            error=f"missing required input(s): {', '.join(missing)}",
        )

    pipeline = app.build(**inputs)
    report = pipeline.validate(geometry=bool(geometry))
    if validate and not report.ok:
        return ApplicationResult(
            name=name, ok=False, validation=report,
            error="application failed static validation; see validation report",
        )

    result = pipeline.run(cores=cores, geometry=geometry, work_dir=work_dir)
    outputs = extract_outputs(app, pipeline, result)
    return ApplicationResult(
        name=name, ok=result.ok, outputs=outputs,
        validation=report, pipeline_result=result,
    )


# ======================================================================
#  Environment (probe + install)
# ======================================================================


def probe() -> List["_environment.CapabilityReadiness"]:
    """Readiness of every capability against its declared tool requirements."""
    return _environment.probe()


def missing_tools() -> Dict[str, Dict[str, List[str]]]:
    """Unique missing requirements → the capabilities that need them."""
    return _environment.missing_tools()


def install_plan() -> Dict[str, object]:
    """Classify missing tools into auto-installable vs. manual (license-restricted)."""
    return _environment.install_plan()


def install_tools(*, run: bool = False, timeout: Optional[float] = None) -> Dict[str, object]:
    """Install the open-source tools DELFIN may legally install.

    License-restricted tools (ORCA, Turbomole) are never installed; their
    guidance is returned under ``plan.manual``.  ``run=False`` (default) only
    plans.
    """
    return _environment.install_tools(run=run, timeout=timeout)


__all__ = [
    # capabilities
    "list_capabilities",
    "describe_capability",
    "catalog",
    "compatible_successors",
    # applications
    "ApplicationResult",
    "register_application",
    "list_applications",
    "describe_application",
    "validate_application",
    "run_application",
    # environment
    "probe",
    "missing_tools",
    "install_plan",
    "install_tools",
]
