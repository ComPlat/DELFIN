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
from delfin.tools._keys import get_key as _get_key
from delfin.tools._keys import list_keys as _list_keys
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


def new_capability_template(name: str, *, category: str = "meta") -> str:
    """A ready-to-fill Python skeleton for a NEW building block (StepAdapter).

    For when no existing tool does what a workflow needs: save the result to
    ``~/.delfin/adapters/<name>.py`` (or ``$DELFIN_ADAPTERS_DIR``) and it is
    auto-discovered as a new capability, usable in pipelines like any built-in.
    """
    cls = "".join(w.capitalize() for w in name.replace("-", "_").split("_")) + "Adapter"
    return f'''"""Custom DELFIN building block: {name}."""

import time
from pathlib import Path

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec


class {cls}(StepAdapter):
    name = "{name}"
    description = "TODO: what this step does"
    produces_geometry = False          # True if it writes an output geometry
    category = "{category}"
    params = (
        # ParamSpec("charge", "int", required=True, description="Molecular charge"),
    )
    consumes = ()                      # e.g. ("geometry",) or ("qc_output",)
    produces = ()                      # e.g. ("qc_output",)
    data_keys = ()                     # e.g. (DataKeySpec("energy_Eh", "float", "Eh"),)
    requires_binaries = ()             # e.g. ("orca",)
    requires_python = ()               # e.g. ("rdkit",)

    def validate_params(self, **kwargs):
        # raise ValueError("'x' parameter is required") for any missing required key
        pass

    def execute(self, work_dir: Path, *, geometry=None, cores: int = 1, **kwargs) -> StepResult:
        start = time.monotonic()
        # TODO: do the work inside work_dir (never os.chdir). Call your tool, then:
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start,
                                 data={{}})


register({cls}())
'''


# ======================================================================
#  Keys (central well-known-parameter vocabulary)
# ======================================================================


def list_keys() -> List[str]:
    """Names of all well-known keys (functional, basis, solvent, …)."""
    return sorted(_list_keys().keys())


def describe_key(name: str):
    """The :class:`KeySpec` (type, default, allowed values) for one key, or ``None``."""
    return _get_key(name)


# ======================================================================
#  Manifest (machine-readable binding contract)
# ======================================================================


def manifest() -> Dict[str, Any]:
    """The full platform manifest (capabilities, applications, keys, schemas)."""
    from delfin.tools.manifest import build_manifest
    return build_manifest()


def manifest_json(*, indent: int = 2) -> str:
    """The full platform manifest as a JSON string."""
    from delfin.tools.manifest import manifest_json as _mj
    return _mj(indent=indent)


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


def validate_spec(spec: Dict[str, Any], *, inputs: Optional[Dict[str, Any]] = None,
                  geometry: bool = False):
    """Validate an ad-hoc PipelineSpec or Application dict *without* registering it.

    The build-loop feedback for an agent assembling a workflow: returns a
    :class:`ValidationReport` (``.ok`` / ``.errors`` / ``.summary()``).  Accepts
    either a raw pipeline spec (``{"name","steps",...}``) or an Application dict
    (``{"name","spec","inputs","outputs"}``); for the latter pass ``inputs`` to
    resolve template params.
    """
    inputs = inputs or {}
    if "spec" in spec and "steps" not in spec:
        from delfin.tools._application import Application
        return Application.from_dict(spec).validate(geometry=geometry, **inputs)
    from delfin.tools._serialize import from_dict
    from delfin.tools.pipeline import PipelineTemplate
    pipeline = from_dict(spec)
    if isinstance(pipeline, PipelineTemplate):
        pipeline = pipeline.build(**inputs)
    return pipeline.validate(geometry=geometry)


def save_application(app: Any, *, directory: Optional[str | Path] = None) -> str:
    """Register *app* and persist it so it shows up everywhere (Pipelines tab …).

    *app* is an :class:`Application` or its ``to_dict()`` dict.  Writes
    ``<name>.json`` into the user applications dir (``$DELFIN_APPLICATIONS_DIR``
    or ``~/.delfin/applications``) and returns the path.  This is how an agent
    makes a workflow it just built permanent.
    """
    import json

    from delfin.tools._application import (Application, _user_application_dirs,
                                           register_application)

    if isinstance(app, dict):
        app = Application.from_dict(app)
    register_application(app)

    target = Path(directory) if directory else _user_application_dirs()[0]
    target.mkdir(parents=True, exist_ok=True)
    path = target / f"{app.name}.json"
    path.write_text(json.dumps(app.to_dict(), indent=2, default=str), encoding="utf-8")
    return str(path)


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
#  Runs (async execution — submit / status / result / list / cancel)
# ======================================================================


def submit_application(
    name: str,
    *,
    cores: int = 1,
    maxcore: Optional[int] = None,
    geometry: Optional[str | Path] = None,
    work_dir: Optional[Path] = None,
    backend: str = "local",
    **inputs: Any,
) -> str:
    """Submit an application run in the background and return its run id.

    ``cores`` is the parallelism (ORCA PAL); ``maxcore`` the memory per core (MB).
    ``backend="local"`` runs on this machine; ``backend="slurm"`` submits an
    sbatch job that runs it on a compute node (results land in the shared run
    store). Non-blocking: poll :func:`run_status` / :func:`run_record`, or
    :func:`wait_run`, and cancel with :func:`cancel_run`.
    """
    from delfin.tools._runtime import get_runtime
    handle = get_runtime().submit_application(
        name, cores=cores, maxcore=maxcore, geometry=geometry, work_dir=work_dir,
        inputs=inputs, backend=backend,
    )
    return handle.id


def run_status(run_id: str) -> Optional[str]:
    """Current status of a run (``pending``/``running``/``success``/…), or ``None``."""
    from delfin.tools._runtime import get_runtime
    rec = get_runtime().get(run_id)
    return rec.status if rec else None


def run_record(run_id: str):
    """The full :class:`RunRecord` (inputs, status, outputs, events, metrics)."""
    from delfin.tools._runtime import get_runtime
    return get_runtime().get(run_id)


def list_runs() -> list:
    """All run records, newest first."""
    from delfin.tools._runtime import get_runtime
    return get_runtime().list_runs()


def run_diagnostics(run_id: str) -> Dict[str, Any]:
    """Failure diagnostics for a run: status, error, events, and its log files.

    The pipeline's output files (ORCA ``.out``, logs) live in the run's work dir
    under ``~/calc``; this lists them so an agent can study what went wrong and
    iterate.  Returns ``{"error": ...}`` if the run is unknown.
    """
    from delfin.tools._runtime import get_runtime

    rec = get_runtime().get(run_id)
    if rec is None:
        return {"error": f"unknown run {run_id!r}"}

    log_files: List[str] = []
    if rec.work_dir:
        wd = Path(rec.work_dir)
        if wd.is_dir():
            for f in sorted(wd.rglob("*")):
                if f.is_file() and (f.suffix in (".out", ".log", ".err", ".txt")):
                    log_files.append(str(f))
    return {
        "id": rec.id,
        "name": rec.name,
        "status": rec.status,
        "error": rec.error,
        "work_dir": rec.work_dir,
        "outputs": rec.outputs,
        "events": rec.events[-20:],
        "log_files": log_files[:50],
    }


def cancel_run(run_id: str) -> bool:
    """Request cooperative cancellation of a run (takes effect before the next step)."""
    from delfin.tools._runtime import get_runtime
    return get_runtime().cancel(run_id)


def wait_run(run_id: str, *, timeout: Optional[float] = None):
    """Block until a run finishes (or *timeout*); returns its :class:`RunRecord`."""
    from delfin.tools._runtime import RunHandle, get_runtime
    return RunHandle(run_id, get_runtime()).wait(timeout=timeout)


def run_metrics() -> Dict[str, Any]:
    """Aggregate the run history into timing / success-rate / cost signals."""
    from delfin.tools._metrics import aggregate_runs
    return aggregate_runs()


def application_keyfile(name: str) -> str:
    """Generate an editable CONTROL.txt-style keys file for an application.

    Each input key shows its default, type and (from the central key vocabulary)
    its allowed values.  Edit it and run via :func:`run_from_keyfile` or
    ``delfin-app run <file>``.
    """
    from delfin.tools._keyfile import application_keyfile as _kf
    return _kf(name)


def run_from_keyfile(path: str, *, cores: int = 1, submit: bool = False, **overrides: Any):
    """Run (or submit) an application described by a keys file."""
    from delfin.tools._keyfile import run_from_keyfile as _run
    return _run(path, cores=cores, submit=submit, **overrides)


# ======================================================================
#  Environment (probe + install)
# ======================================================================


def probe() -> List["_environment.CapabilityReadiness"]:
    """Readiness of every capability against its declared tool requirements."""
    return _environment.probe()


def missing_tools() -> Dict[str, Dict[str, List[str]]]:
    """Unique missing requirements → the capabilities that need them."""
    return _environment.missing_tools()


def tool_inventory() -> List[Dict[str, Any]]:
    """Per-tool install status (installed ✓ / missing, policy, hint, used_by)."""
    return _environment.tool_inventory()


def install_plan() -> Dict[str, object]:
    """Classify missing tools into auto-installable vs. manual (license-restricted)."""
    return _environment.install_plan()


def install_tools(
    *,
    run: bool = False,
    select: Optional[List[str]] = None,
    timeout: Optional[float] = None,
) -> Dict[str, object]:
    """Install the open-source tools DELFIN may legally install.

    ``select`` restricts which tools to install (subset of the plan's
    auto-installable tools); ``None`` installs all. License-restricted tools
    (ORCA, Turbomole, Multiwfn, …) are never installed — their guidance is
    returned under ``plan.manual``. ``run=False`` (default) only plans.
    """
    return _environment.install_tools(run=run, select=select, timeout=timeout)


__all__ = [
    # capabilities
    "list_capabilities",
    "describe_capability",
    "catalog",
    "compatible_successors",
    "new_capability_template",
    # keys + manifest
    "list_keys",
    "describe_key",
    "manifest",
    "manifest_json",
    # applications
    "ApplicationResult",
    "register_application",
    "list_applications",
    "describe_application",
    "validate_application",
    "validate_spec",
    "save_application",
    "run_application",
    # runs (async)
    "submit_application",
    "run_status",
    "run_record",
    "list_runs",
    "run_diagnostics",
    "cancel_run",
    "wait_run",
    "run_metrics",
    "application_keyfile",
    "run_from_keyfile",
    # environment
    "probe",
    "missing_tools",
    "tool_inventory",
    "install_plan",
    "install_tools",
]
