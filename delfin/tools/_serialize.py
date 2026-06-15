"""Serialize pipelines to/from plain data — the "building blocks as data" bridge.

A :class:`~delfin.tools.pipeline.Pipeline` / ``PipelineTemplate`` is normally
Python code.  This module gives the declarative skeleton (steps, kwargs,
branches, defaults, labels, and the data-expressible flow control) a round-trip
to a plain ``dict`` — hence JSON/YAML.  That lets a workflow be emitted by a UI
or an agent, validated with :meth:`Pipeline.validate` before running, saved,
diffed in git, or shared.

The callable boundary is explicit, never silent:

* ``if`` / ``loop`` conditions round-trip as the small data DSL lifted from the
  ``delfin-pipeline`` YAML loader (``last.data.key > 0`` etc.).  The built
  callable carries its source so :func:`to_dict` can recover it.
* ``compute`` / ``transform`` / ``reactive`` / ``fan_out`` / ``map`` callables
  round-trip via a **named-callable registry** (:func:`register_callable`) — the
  spec stores a ``ref`` resolved on load.
* An anonymous callable that is neither a DSL form nor registered cannot be
  serialized: :func:`to_dict` raises :class:`PipelineSerializationError` by
  default, or — with ``strict=False`` — emits a visible marker and lists the
  step under ``_non_serializable_steps``.

The JSON shape is a superset of the existing ``delfin-pipeline`` YAML format, so
every existing pipeline file still loads.
"""

from __future__ import annotations

import re
from typing import Any, Callable, Dict, List, Optional

SCHEMA_VERSION = "1"


class PipelineSerializationError(ValueError):
    """Raised when a pipeline cannot be losslessly serialized."""


# --- Named-callable registry ----------------------------------------------

_CALLABLE_REGISTRY: Dict[str, Callable] = {}


def register_callable(name: str) -> Callable[[Callable], Callable]:
    """Decorator: register *fn* under *name* so it can be referenced in a spec.

        @register_callable("calc_redox")
        def calc_redox(results, last, work_dir): ...

    A step that uses ``calc_redox`` then serializes as ``{"ref": "calc_redox"}``
    and :func:`from_dict` resolves it back.
    """

    def _decorator(fn: Callable) -> Callable:
        _CALLABLE_REGISTRY[name] = fn
        try:
            fn._delfin_ref = name  # type: ignore[attr-defined]
        except (AttributeError, TypeError):
            pass
        return fn

    return _decorator


def _ref_of(fn: Optional[Callable]) -> Optional[str]:
    return getattr(fn, "_delfin_ref", None) if fn is not None else None


# --- Condition / loop DSL (lifted from cli_pipeline, now round-trippable) ---

_OPS = {
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
    ">": lambda a, b: a > b,
    "<": lambda a, b: a < b,
    ">=": lambda a, b: a >= b,
    "<=": lambda a, b: a <= b,
}


def _coerce_val(raw: str) -> Any:
    raw = raw.strip()
    try:
        return int(raw)
    except ValueError:
        pass
    try:
        return float(raw)
    except ValueError:
        return raw.strip("'\"")


def build_condition(expr: str) -> Callable:
    """Build an ``add_if`` condition callable from a DSL string.

    Supports ``last.ok``, ``not last.ok`` and ``last.data.KEY OP VALUE``.  The
    returned callable carries ``_delfin_dsl`` so it can be serialized back.
    """
    expr = expr.strip()
    if expr == "not last.ok":
        fn = lambda results, last: last is not None and not last.ok  # noqa: E731
    elif expr == "last.ok":
        fn = lambda results, last: last is not None and last.ok  # noqa: E731
    else:
        m = re.match(r"last\.data\.(\w+)\s*(==|!=|>=|<=|>|<)\s*(.+)", expr)
        if not m:
            raise ValueError(f"Cannot parse condition expression: {expr!r}")
        key, op, val = m.group(1), m.group(2), _coerce_val(m.group(3))
        op_fn = _OPS[op]
        fn = lambda results, last, _k=key, _fn=op_fn, _v=val: (  # noqa: E731
            last is not None and _fn(last.data.get(_k), _v)
        )
    fn._delfin_dsl = expr  # type: ignore[attr-defined]
    fn._delfin_dsl_kind = "condition"  # type: ignore[attr-defined]
    return fn


def build_until(expr: str) -> Callable:
    """Build an ``add_loop`` until-callable from a DSL string.

    Supports ``result.data.KEY OP VALUE``.  Carries ``_delfin_dsl`` for round-trip.
    """
    expr = expr.strip()
    m = re.match(r"result\.data\.(\w+)\s*(==|!=|>=|<=|>|<)\s*(.+)", expr)
    if not m:
        raise ValueError(f"Cannot parse until expression: {expr!r}")
    key, op, val = m.group(1), m.group(2), _coerce_val(m.group(3))
    op_fn = _OPS[op]
    fn = lambda result, iteration, _k=key, _fn=op_fn, _v=val: (  # noqa: E731
        _fn(result.data.get(_k), _v)
    )
    fn._delfin_dsl = expr  # type: ignore[attr-defined]
    fn._delfin_dsl_kind = "until"  # type: ignore[attr-defined]
    return fn


def _dsl_of(fn: Optional[Callable]) -> Optional[str]:
    return getattr(fn, "_delfin_dsl", None) if fn is not None else None


# --- from_dict (data -> pipeline) -----------------------------------------

# Keys with structural meaning in a step dict; everything else is a flat kwarg.
_RESERVED = {
    "step", "type", "label", "condition", "until", "max_iter",
    "max_attempts", "delay", "ref", "module", "function", "kwargs",
    "geometry", "pipeline",
}


def _resolve_named_callable(spec: Dict[str, Any], kind: str) -> Callable:
    """Resolve a callable referenced by ``ref`` or legacy ``module``+``function``."""
    ref = spec.get("ref")
    if ref is not None:
        if ref not in _CALLABLE_REGISTRY:
            raise PipelineSerializationError(
                f"{kind} step references unknown callable {ref!r}; "
                f"register it with @register_callable({ref!r})"
            )
        return _CALLABLE_REGISTRY[ref]
    module, function = spec.get("module"), spec.get("function")
    if module and function:
        import importlib
        return getattr(importlib.import_module(module), function)
    raise ValueError(f"{kind} step requires 'ref' or 'module'+'function'")


def _step_kwargs(spec: Dict[str, Any]) -> Dict[str, Any]:
    """Merge nested ``kwargs`` with legacy flat sibling kwargs."""
    kwargs = dict(spec.get("kwargs", {}))
    for k, v in spec.items():
        if k not in _RESERVED:
            kwargs.setdefault(k, v)
    return kwargs


def _add_step(pipe, spec: Dict[str, Any], *, is_template: bool) -> None:
    step_type = spec.get("type", "normal")
    label = spec.get("label", "")
    kwargs = _step_kwargs(spec)

    # Templates only support plain add() (no flow control / meta steps).
    if is_template:
        if "step" not in spec:
            raise ValueError(f"Each template step needs a 'step' key, got: {spec}")
        pipe.add(spec["step"], label=label, **kwargs)
        return

    if step_type == "checkpoint":
        pipe.add_checkpoint(label=label or "checkpoint")
        return
    if step_type == "compute":
        pipe.add_compute(_resolve_named_callable(spec, "compute"), label=label or "compute")
        return
    if step_type == "reactive":
        pipe.add_reactive(_resolve_named_callable(spec, "reactive"), label=label or "reactive")
        return
    if step_type == "sub_pipeline":
        sub = from_dict(spec["pipeline"])
        pipe.add_sub_pipeline(sub, label=label)
        return

    if "step" not in spec:
        raise ValueError(f"Step must have a 'step' key, got: {spec}")
    step_name = spec["step"]
    geometry = spec.get("geometry")

    if step_type in ("normal", "step"):
        pipe.add(step_name, label=label, geometry=geometry, **kwargs)
    elif step_type == "if":
        cond = spec.get("condition")
        if not cond:
            raise ValueError("'if' step requires 'condition'")
        pipe.add_if(build_condition(cond), step_name, label=label, **kwargs)
    elif step_type == "loop":
        until = spec.get("until")
        if not until:
            raise ValueError("'loop' step requires 'until'")
        pipe.add_loop(step_name, until=build_until(until),
                      max_iter=spec.get("max_iter", 10), label=label, **kwargs)
    elif step_type == "retry":
        pipe.add_retry(step_name, max_attempts=spec.get("max_attempts", 3),
                       delay=spec.get("delay", 0.0), label=label, **kwargs)
    elif step_type == "transform":
        pipe.add_transform(step_name, _resolve_named_callable(spec, "transform"),
                           label=label, **kwargs)
    elif step_type == "fan_out":
        pipe.add_fan_out(step_name, _resolve_named_callable(spec, "fan_out"),
                         label=label, **kwargs)
    elif step_type == "map":
        pipe.add_map(step_name, _resolve_named_callable(spec, "map"),
                     label=label, **kwargs)
    else:
        raise ValueError(f"Unknown step type: {step_type!r}")


def from_dict(data: Dict[str, Any]):
    """Build a :class:`Pipeline` or ``PipelineTemplate`` from a plain dict.

    Accepts both the new nested-``kwargs`` form and the legacy flat form, so it
    is a drop-in for the existing ``delfin-pipeline`` YAML loader.
    """
    from delfin.tools.pipeline import Pipeline, PipelineTemplate

    if not isinstance(data, dict):
        raise ValueError(f"Pipeline spec must be a mapping, got {type(data).__name__}")
    if "name" not in data:
        raise ValueError("Pipeline spec must have a 'name' field")
    if "steps" not in data:
        raise ValueError("Pipeline spec must have a 'steps' list")

    is_template = bool(data.get("template", False) or data.get("kind") == "template")
    name = data["name"]
    defaults = data.get("defaults", {})
    pipe = (PipelineTemplate(name, defaults=defaults) if is_template
            else Pipeline(name, defaults=defaults))

    for spec in data["steps"]:
        if not isinstance(spec, dict):
            raise ValueError(f"Each step must be a dict, got: {spec}")
        _add_step(pipe, spec, is_template=is_template)

    for bname, bsteps in data.get("branches", {}).items():
        branch = pipe.branch(bname)
        for spec in bsteps:
            if not isinstance(spec, dict):
                raise ValueError(f"Branch step must be a dict, got: {spec}")
            _add_step(branch, spec, is_template=is_template)

    return pipe


# --- to_dict (pipeline -> data) -------------------------------------------


def _spec_to_dict(spec, *, strict: bool, unserializable: List[str]) -> Dict[str, Any]:
    """Serialize one _StepSpec to its canonical dict form."""
    label = spec.label or ""

    def _named(kind: str, fn) -> Dict[str, Any]:
        ref = _ref_of(fn)
        if ref is None:
            unserializable.append(spec.label or kind)
            if strict:
                raise PipelineSerializationError(
                    f"step {spec.label or kind!r} uses an unregistered {kind} "
                    f"callable; register it with @register_callable(...) to serialize"
                )
            return {"type": kind, "label": label, "ref": None,
                    "_unserializable": kind}
        return {"type": kind, "label": label, "ref": ref}

    def _dsl(kind_key: str, fn, extra: Dict[str, Any]) -> Dict[str, Any]:
        dsl = _dsl_of(fn)
        if dsl is None:
            unserializable.append(spec.label or spec.step_name)
            if strict:
                raise PipelineSerializationError(
                    f"step {spec.label or spec.step_name!r} uses a non-DSL "
                    f"{kind_key} callable that cannot be serialized"
                )
            out = {"step": spec.step_name, "label": label, kind_key: None,
                   "_unserializable": kind_key, "kwargs": dict(spec.kwargs)}
            out.update(extra)
            return out
        out = {"step": spec.step_name, "label": label, kind_key: dsl,
               "kwargs": dict(spec.kwargs)}
        out.update(extra)
        return out

    # Meta steps first.
    if spec.checkpoint:
        return {"type": "checkpoint", "label": label}
    if spec.sub_pipeline is not None:
        return {"type": "sub_pipeline", "label": label,
                "pipeline": to_dict(spec.sub_pipeline, strict=strict)}
    if spec.compute_fn is not None:
        return _named("compute", spec.compute_fn)
    if spec.reactive_fn is not None:
        return _named("reactive", spec.reactive_fn)

    # Flow-control variants on a real step.
    if spec.condition is not None:
        return {**_dsl("condition", spec.condition, {}), "type": "if"}
    if spec.loop_until is not None and spec.loop_max > 1:
        return {**_dsl("until", spec.loop_until, {"max_iter": spec.loop_max}),
                "type": "loop"}
    if spec.transform is not None:
        d = _named("transform", spec.transform)
        d.update({"step": spec.step_name, "kwargs": dict(spec.kwargs)})
        return d
    if spec.fan_out_fn is not None:
        d = _named("fan_out", spec.fan_out_fn)
        d.update({"step": spec.step_name, "kwargs": dict(spec.kwargs)})
        return d
    if spec.map_items is not None:
        d = _named("map", spec.map_items)
        d.update({"step": spec.step_name, "kwargs": dict(spec.kwargs)})
        return d
    if spec.retry_max and spec.retry_max > 1:
        return {"step": spec.step_name, "type": "retry", "label": label,
                "max_attempts": spec.retry_max, "delay": spec.retry_delay,
                "kwargs": dict(spec.kwargs)}

    # Plain step.
    out: Dict[str, Any] = {"step": spec.step_name, "label": label,
                           "kwargs": dict(spec.kwargs)}
    if spec.geometry_override is not None:
        out["geometry"] = str(spec.geometry_override)
    return out


def to_dict(pipeline, *, strict: bool = True) -> Dict[str, Any]:
    """Serialize a Pipeline / PipelineTemplate to a plain dict.

    ``strict=True`` (default) raises :class:`PipelineSerializationError` on any
    step whose callable cannot be represented as data; ``strict=False`` emits a
    visible marker instead and lists such steps under ``_non_serializable_steps``.
    """
    from delfin.tools.pipeline import PipelineTemplate

    kind = "template" if isinstance(pipeline, PipelineTemplate) else "pipeline"
    unserializable: List[str] = []
    steps = [_spec_to_dict(s, strict=strict, unserializable=unserializable)
             for s in pipeline._trunk]
    branches = {
        bname: [_spec_to_dict(s, strict=strict, unserializable=unserializable)
                for s in bpipe._trunk]
        for bname, bpipe in pipeline._branches.items()
    }
    out: Dict[str, Any] = {
        "delfin_pipeline_spec": SCHEMA_VERSION,
        "name": pipeline.name,
        "kind": kind,
        "defaults": dict(pipeline._defaults),
        "steps": steps,
    }
    if kind == "template":
        out["template"] = True
    if branches:
        out["branches"] = branches
    if unserializable and not strict:
        out["_non_serializable_steps"] = unserializable
    return out


# --- JSON helpers ----------------------------------------------------------


def to_json(pipeline, *, strict: bool = True, indent: int = 2) -> str:
    import json
    return json.dumps(to_dict(pipeline, strict=strict), indent=indent, default=str)


def from_json(text: str):
    import json
    return from_dict(json.loads(text))


__all__ = [
    "SCHEMA_VERSION",
    "PipelineSerializationError",
    "register_callable",
    "build_condition",
    "build_until",
    "from_dict",
    "to_dict",
    "to_json",
    "from_json",
]
