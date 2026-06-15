"""Build-time validation for pipelines — catch broken wiring before running.

Walks a :class:`~delfin.tools.pipeline.Pipeline` (trunk + branches + nested
sub-pipelines) *without executing anything*, threading an abstract set of
available capability tags.  For each step it checks that the adapter's declared
required parameters are supplied (or auto-wireable) and that its consumed
capability ports (``geometry`` etc.) are produced upstream.  External binaries
and Python modules are checked best-effort (the run host may differ), so those
surface as warnings rather than errors.

Soundness over completeness: steps the validator cannot see through — pure
``add_compute`` blocks, ``add_reactive`` dynamic insertion, ``add_transform``
kwargs rewrites — relax downstream checks to warnings instead of raising false
errors.  A clean report (``ok``) means "no statically detectable problem", not
"guaranteed to run".
"""

from __future__ import annotations

import enum
import importlib.util
import shutil
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

from delfin.tools._registry import get as _get_adapter
from delfin.tools._wiring import can_autowire


class DiagLevel(enum.Enum):
    OK = "ok"
    WARNING = "warning"
    ERROR = "error"


@dataclass
class StepDiagnostic:
    index: int
    label: str
    step_name: str
    level: DiagLevel
    location: str = "trunk"          # "trunk" | "branch:<name>" | "<loc>/sub:<name>"
    missing_params: Tuple[str, ...] = ()
    missing_inputs: Tuple[str, ...] = ()
    messages: Tuple[str, ...] = ()


@dataclass
class ValidationReport:
    pipeline_name: str
    diagnostics: List[StepDiagnostic] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(d.level is DiagLevel.ERROR for d in self.diagnostics)

    @property
    def errors(self) -> List[StepDiagnostic]:
        return [d for d in self.diagnostics if d.level is DiagLevel.ERROR]

    @property
    def warnings(self) -> List[StepDiagnostic]:
        return [d for d in self.diagnostics if d.level is DiagLevel.WARNING]

    def summary(self) -> str:
        head = "OK" if self.ok else "FAILED"
        lines = [f"Pipeline '{self.pipeline_name}' validation: {head}"]
        for d in self.diagnostics:
            tag = {DiagLevel.OK: "ok", DiagLevel.WARNING: "warn",
                   DiagLevel.ERROR: "ERR "}[d.level]
            lines.append(f"  [{tag}] {d.label} ({d.step_name}) @{d.location}")
            for m in d.messages:
                lines.append(f"         {m}")
        return "\n".join(lines)


# Capabilities a pure-python / dynamic block might conservatively be assumed to
# emit, so the validator does not raise false "missing geometry" past it.
_DYNAMIC_ASSUMED_CAPS = frozenset({"geometry"})


def _unresolved_placeholders(kwargs: Dict[str, Any]) -> List[str]:
    out = []
    for k, v in kwargs.items():
        if isinstance(v, str) and "{" in v and "}" in v:
            out.append(k)
    return out


def _check_binaries_and_python(contract) -> List[str]:
    msgs: List[str] = []
    for b in sorted(contract.requires_binaries):
        if shutil.which(b) is None:
            msgs.append(f"binary not found on PATH: {b!r} (may exist on run host)")
    for m in sorted(contract.requires_python):
        try:
            found = importlib.util.find_spec(m) is not None
        except (ImportError, ValueError):
            found = False
        if not found:
            msgs.append(f"python module not importable: {m!r}")
    return msgs


def _check_real_step(
    index: int,
    spec,
    avail_caps: Set[str],
    defaults: Dict[str, Any],
    location: str,
    relax: bool,
) -> Tuple[StepDiagnostic, Set[str]]:
    """Validate one real adapter step; return (diagnostic, produced capabilities)."""
    adapter = _get_adapter(spec.step_name)
    if adapter is None:
        return (
            StepDiagnostic(
                index, spec.label, spec.step_name, DiagLevel.ERROR, location,
                messages=(f"unknown step '{spec.step_name}' (not registered)",),
            ),
            set(),
        )

    contract = adapter.contract()
    eff: Dict[str, Any] = {**defaults, **spec.kwargs}
    messages: List[str] = []
    level = DiagLevel.OK

    # Flow-control that can inject params/geometry at runtime relaxes checks.
    relax_step = relax or spec.transform is not None or spec.map_items is not None
    geometry_supplied = (
        "geometry" in avail_caps
        or spec.geometry_override is not None
        or spec.fan_out_fn is not None
        or spec.map_items is not None
    )

    # 1. required parameters
    missing_params: List[str] = []
    for p in contract.required_params:
        if p in eff:
            continue
        if can_autowire(p, frozenset(avail_caps)):
            continue
        missing_params.append(p)
    if missing_params:
        if relax_step:
            level = DiagLevel.WARNING
            messages.append(
                "missing required param(s) "
                + ", ".join(missing_params)
                + " — may be injected at runtime"
            )
        else:
            level = DiagLevel.ERROR
            messages.append("missing required param(s): " + ", ".join(missing_params))

    # 2. consumed capability ports
    missing_inputs: List[str] = []
    for tag in sorted(contract.consumes):
        if tag == "geometry":
            satisfied = geometry_supplied
        else:
            satisfied = tag in avail_caps
        if not satisfied:
            missing_inputs.append(tag)
    if missing_inputs:
        if relax_step:
            if level is DiagLevel.OK:
                level = DiagLevel.WARNING
            messages.append(
                "input capability not statically available: "
                + ", ".join(missing_inputs)
            )
        else:
            level = DiagLevel.ERROR
            messages.append("missing input capability: " + ", ".join(missing_inputs))

    # 3. unresolved template placeholders (advisory)
    placeholders = _unresolved_placeholders(eff)
    if placeholders:
        if level is DiagLevel.OK:
            level = DiagLevel.WARNING
        messages.append("unresolved placeholder(s): " + ", ".join(sorted(placeholders)))

    # 4. binaries / python deps (advisory)
    env_msgs = _check_binaries_and_python(contract)
    if env_msgs:
        if level is DiagLevel.OK:
            level = DiagLevel.WARNING
        messages.extend(env_msgs)

    diag = StepDiagnostic(
        index, spec.label, spec.step_name, level, location,
        missing_params=tuple(missing_params),
        missing_inputs=tuple(missing_inputs),
        messages=tuple(messages),
    )
    return diag, set(contract.produces)


def _walk(
    trunk, branches, *,
    initial_caps: Set[str],
    defaults: Dict[str, Any],
    location: str,
) -> Tuple[List[StepDiagnostic], Set[str]]:
    avail: Set[str] = set(initial_caps)
    diags: List[StepDiagnostic] = []
    dynamic = False

    for i, spec in enumerate(trunk):
        # --- meta steps -------------------------------------------------
        if spec.sub_pipeline is not None:
            sub = spec.sub_pipeline
            diags.append(StepDiagnostic(
                i, spec.label, "_sub_pipeline", DiagLevel.OK, location,
                messages=(f"sub-pipeline '{sub.name}'",),
            ))
            sub_diags, sub_caps = _walk(
                sub._trunk, sub._branches,
                initial_caps=avail,
                defaults={**defaults, **getattr(sub, "_defaults", {})},
                location=f"{location}/sub:{sub.name}",
            )
            diags.extend(sub_diags)
            avail |= sub_caps
            continue
        if spec.checkpoint:
            diags.append(StepDiagnostic(
                i, spec.label, "_checkpoint", DiagLevel.OK, location,
            ))
            continue
        if spec.compute_fn is not None:
            diags.append(StepDiagnostic(
                i, spec.label, "_compute", DiagLevel.OK, location,
                messages=("pure-python compute step (not statically analysable)",),
            ))
            dynamic = True
            avail |= _DYNAMIC_ASSUMED_CAPS
            continue
        if spec.reactive_fn is not None:
            diags.append(StepDiagnostic(
                i, spec.label, "_reactive", DiagLevel.OK, location,
                messages=("dynamic step insertion (not statically analysable)",),
            ))
            dynamic = True
            avail |= _DYNAMIC_ASSUMED_CAPS
            continue

        # --- real adapter step -----------------------------------------
        diag, produces = _check_real_step(
            i, spec, avail, defaults, location, relax=dynamic,
        )
        diags.append(diag)
        avail |= produces

    final = set(avail)

    for bname, bpipe in branches.items():
        b_defaults = {**defaults, **getattr(bpipe, "_defaults", {})}
        b_diags, _ = _walk(
            bpipe._trunk, bpipe._branches,
            initial_caps=avail,
            defaults=b_defaults,
            location=f"branch:{bname}",
        )
        diags.extend(b_diags)

    return diags, final


def validate_pipeline(
    pipeline,
    *,
    geometry_available: bool = False,
    defaults: Optional[Dict[str, Any]] = None,
) -> ValidationReport:
    """Statically validate *pipeline* and return a :class:`ValidationReport`."""
    initial: Set[str] = {"geometry"} if geometry_available else set()
    base_defaults = defaults if defaults is not None else getattr(pipeline, "_defaults", {})
    diags, _ = _walk(
        pipeline._trunk, pipeline._branches,
        initial_caps=initial,
        defaults=dict(base_defaults or {}),
        location="trunk",
    )
    return ValidationReport(pipeline_name=getattr(pipeline, "name", "pipeline"),
                            diagnostics=diags)


__all__ = ["DiagLevel", "StepDiagnostic", "ValidationReport", "validate_pipeline"]
