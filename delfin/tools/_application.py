"""Applications — workflows as named, versioned, callable functions.

An :class:`Application` binds a serialized pipeline (a PipelineSpec, see
:mod:`delfin.tools._serialize`) to a *contract*: declared **inputs** (reusing
:class:`~delfin.tools._spec.ParamSpec`) and declared **outputs** (named results
mapped from a step's ``data`` keys).  That turns a workflow into a function the
agent system / a UI / another tool can introspect and call by name with named
arguments, getting named results back — without reaching into pipeline internals.

Applications are fully serializable (``to_dict`` / ``from_dict``), so they can be
registered, listed, shipped over the wire, and bound to from outside Python.

    app = Application.from_pipeline(
        redox_template, name="redox", inputs=[...], outputs=[...])
    register_application(app)

    # later, via the platform facade:
    result = run_application("redox", smiles="CCO", charge=0)
    result.outputs  # {"e_oxidation": -123.4, ...}
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from delfin.tools._spec import ParamSpec


@dataclass(frozen=True)
class OutputSpec:
    """Declares one named application output and how to extract it.

    The value is read from ``result.data[key]`` of the step identified by
    ``step`` (its pipeline label, or its step name as a fallback) within the
    given ``branch`` (``""`` = the trunk).
    """

    name: str
    step: str                 # step label (preferred) or step name to read from
    key: str                  # data key within that step's StepResult.data
    branch: str = ""          # branch name; "" = trunk
    type: str = "float"
    unit: str = ""
    description: str = ""


@dataclass(frozen=True)
class Application:
    """A named, versioned workflow with a typed input/output contract."""

    name: str
    description: str = ""
    version: str = "1"
    category: str = ""
    inputs: Tuple[ParamSpec, ...] = ()
    outputs: Tuple[OutputSpec, ...] = ()
    spec: Dict[str, Any] = field(default_factory=dict)   # serialized PipelineSpec

    @property
    def required_inputs(self) -> Tuple[str, ...]:
        return tuple(p.name for p in self.inputs if p.required)

    def missing_inputs(self, values: Dict[str, Any]) -> List[str]:
        """Required input names not present in *values*."""
        return [p.name for p in self.inputs if p.required and p.name not in values]

    # --- construction -------------------------------------------------

    @classmethod
    def from_pipeline(
        cls,
        pipeline,
        *,
        name: str,
        description: str = "",
        version: str = "1",
        category: str = "",
        inputs=(),
        outputs=(),
    ) -> "Application":
        """Build an Application from a Pipeline / PipelineTemplate.

        The pipeline is serialized via :func:`delfin.tools._serialize.to_dict`
        (strict), so any non-serializable callable surfaces here rather than at
        run time.
        """
        from delfin.tools._serialize import to_dict

        return cls(
            name=name,
            description=description,
            version=version,
            category=category,
            inputs=tuple(inputs),
            outputs=tuple(outputs),
            spec=to_dict(pipeline, strict=True),
        )

    def with_input_defaults(self, values: Dict[str, Any]) -> Dict[str, Any]:
        """*values* augmented with each declared input's default (explicit wins).

        So a caller need only pass the required inputs: every optional input with
        a declared default is filled, which resolves the template's ``{...}``
        placeholders that would otherwise leak literally into a run.
        """
        merged = {p.name: p.default for p in self.inputs if p.default is not None}
        merged.update(values)
        return merged

    def build(self, **values: Any):
        """Materialise a concrete Pipeline, filling template params with *values*."""
        from delfin.tools._serialize import from_dict
        from delfin.tools.pipeline import PipelineTemplate

        obj = from_dict(self.spec)
        if isinstance(obj, PipelineTemplate):
            return obj.build(**self.with_input_defaults(values))
        return obj  # already a concrete Pipeline

    def validate(self, *, geometry: bool = False, **values: Any):
        """Statically validate the materialised pipeline (no execution)."""
        return self.build(**values).validate(geometry=geometry)

    # --- serialization ------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "delfin_application": self.version,
            "name": self.name,
            "description": self.description,
            "category": self.category,
            "inputs": [_paramspec_to_dict(p) for p in self.inputs],
            "outputs": [_outputspec_to_dict(o) for o in self.outputs],
            "spec": self.spec,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Application":
        return cls(
            name=data["name"],
            description=data.get("description", ""),
            version=str(data.get("delfin_application", data.get("version", "1"))),
            category=data.get("category", ""),
            inputs=tuple(_paramspec_from_dict(d) for d in data.get("inputs", [])),
            outputs=tuple(_outputspec_from_dict(d) for d in data.get("outputs", [])),
            spec=data.get("spec", {}),
        )


# --- ParamSpec / OutputSpec (de)serialization -----------------------------


def _paramspec_to_dict(p: ParamSpec) -> Dict[str, Any]:
    d: Dict[str, Any] = {"name": p.name, "type": p.type, "required": p.required}
    if p.default is not None:
        d["default"] = p.default
    if p.enum is not None:
        d["enum"] = list(p.enum)
    if p.description:
        d["description"] = p.description
    if p.unit:
        d["unit"] = p.unit
    return d


def _paramspec_from_dict(d: Dict[str, Any]) -> ParamSpec:
    return ParamSpec(
        name=d["name"],
        type=d.get("type", "str"),
        required=d.get("required", False),
        default=d.get("default"),
        enum=tuple(d["enum"]) if d.get("enum") else None,
        description=d.get("description", ""),
        unit=d.get("unit", ""),
    )


def _outputspec_to_dict(o: OutputSpec) -> Dict[str, Any]:
    d: Dict[str, Any] = {"name": o.name, "step": o.step, "key": o.key}
    if o.branch:
        d["branch"] = o.branch
    if o.type != "float":
        d["type"] = o.type
    if o.unit:
        d["unit"] = o.unit
    if o.description:
        d["description"] = o.description
    return d


def _outputspec_from_dict(d: Dict[str, Any]) -> OutputSpec:
    return OutputSpec(
        name=d["name"],
        step=d["step"],
        key=d["key"],
        branch=d.get("branch", ""),
        type=d.get("type", "float"),
        unit=d.get("unit", ""),
        description=d.get("description", ""),
    )


# --- output extraction ----------------------------------------------------


def _label_result_map(pipeline, result) -> Dict[tuple, Any]:
    """Best-effort ``{(branch, label_or_stepname): StepResult}`` for a run.

    Built by zipping the pipeline's declarative specs with the result list (1:1
    for ordinary/looped steps).  Dynamically-inserted steps (``add_reactive``)
    are not label-addressable; outputs referencing them resolve to ``None``.
    """
    out: Dict[tuple, Any] = {}

    def _index(branch: str, specs, results) -> None:
        for spec, res in zip(specs, results):
            out[(branch, spec.label)] = res
            out.setdefault((branch, res.step_name), res)

    _index("", pipeline._trunk, result.results)
    for bname, bres in result.branch_results.items():
        bpipe = pipeline._branches.get(bname)
        _index(bname, bpipe._trunk if bpipe else [], bres.results)
    return out


def extract_outputs(application: Application, pipeline, result) -> Dict[str, Any]:
    """Map a :class:`PipelineResult` to the application's named outputs."""
    label_map = _label_result_map(pipeline, result)
    out: Dict[str, Any] = {}
    for o in application.outputs:
        res = label_map.get((o.branch, o.step))
        out[o.name] = res.data.get(o.key) if (res is not None and res.data) else None
    return out


# --- registry -------------------------------------------------------------

_APP_REGISTRY: Dict[str, Application] = {}
_APP_DISCOVERED = False


def register_application(app: Application) -> Application:
    """Register an application instance under its name."""
    if not app.name:
        raise ValueError("Application has no name")
    _APP_REGISTRY[app.name] = app
    return app


def _user_application_dirs():
    import os
    from pathlib import Path

    dirs = []
    env = os.environ.get("DELFIN_APPLICATIONS_DIR")
    if env:
        dirs.append(Path(env))
    dirs.append(Path.home() / ".delfin" / "applications")
    return dirs


def _load_user_applications() -> None:
    """Load user applications from ~/.delfin/applications/ (``*.json`` / ``*.py``).

    Drop a serialized application (``delfin-app describe <name> > foo.json``,
    then edit) or a ``.py`` module that calls :func:`register_application`.
    Failures per file are skipped so one bad file never blocks discovery.
    """
    import json

    for directory in _user_application_dirs():
        try:
            if not directory.is_dir():
                continue
            for f in sorted(directory.glob("*.json")):
                try:
                    register_application(
                        Application.from_dict(json.loads(f.read_text(encoding="utf-8")))
                    )
                except Exception:  # noqa: BLE001
                    pass
            for f in sorted(directory.glob("*.py")):
                try:
                    import importlib.util
                    spec = importlib.util.spec_from_file_location(
                        f"delfin_user_app_{f.stem}", f)
                    if spec and spec.loader:
                        module = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(module)
                except Exception:  # noqa: BLE001
                    pass
        except Exception:  # noqa: BLE001
            continue


def _ensure_discovered() -> None:
    global _APP_DISCOVERED
    if _APP_DISCOVERED:
        return
    _APP_DISCOVERED = True
    try:
        __import__("delfin.tools.applications")
    except ImportError:
        pass
    _load_user_applications()


def get_application(name: str) -> Optional[Application]:
    """Look up a registered application by name (triggers lazy discovery)."""
    _ensure_discovered()
    return _APP_REGISTRY.get(name)


def list_applications() -> Dict[str, Application]:
    """All registered applications (triggers lazy discovery)."""
    _ensure_discovered()
    return dict(_APP_REGISTRY)


__all__ = [
    "OutputSpec",
    "Application",
    "extract_outputs",
    "register_application",
    "get_application",
    "list_applications",
]
