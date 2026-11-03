"""Workflow descriptors — ONE contract for every DELFIN workflow.

A :class:`WorkflowDescriptor` extends :class:`~delfin.tools._application.Application`
(no parallel type) with the fields that make a workflow *introspectable and
extensible* by an agent at build time:

    kind            — "classic-control" | "pipeline" | "user"
    purpose         — one sentence on what the workflow is for
    entry           — how it is started (the dispatch gate / CLI)
    control_keys    — the CONTROL.txt keys it reads, with type, default,
                      allowed values and what each changes — DERIVED from
                      ``delfin.common.control_validator.CONTROL_FIELD_SPECS``
                      so a descriptor can never drift from the validator
    steps           — the building blocks (pipelines) or the ORCA jobs a
                      classic workflow writes (described as data)
    artifacts       — files / directories the workflow produces
    extension_points— where a step may be inserted (before/after/in place)
                      and what the inserted step receives
    resume          — resume / recovery behaviour
    backend         — local vs Slurm

The classic CONTROL.txt workflows (classic, OCCUPIER, ESD, IMAG, redox,
CO2 coordinator) are registered through this descriptor WITHOUT changing
how they run: the descriptor describes their jobs, CONTROL keys and entry
points; the real engines stay untouched.  ``register_classic_workflows``
only ever registers descriptions.

``list_workflows`` lists classic-control workflows, pipeline applications
and user workflows **together**, so an agent sees the whole repertoire as
one catalog.  ``add_step`` implements the extension rule "add step X after
step Y of workflow W": it produces a *new* checked pipeline spec that runs
W's steps unchanged and then X (or before / in place where the descriptor
allows), validated by the existing checks (``user_workflows.check_spec``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from delfin.tools._application import Application

__all__ = [
    "ControlKeySpec",
    "ExtensionPoint",
    "WorkflowDescriptor",
    "register_classic_workflows",
    "workflow_for",
    "list_workflows",
    "add_step",
    "control_key_spec",
    "CLASSIC_WORKFLOWS",
]


# ---------------------------------------------------------------------------
# Contract value types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ControlKeySpec:
    """Declarative description of one CONTROL.txt key a workflow reads.

    Derived from :class:`delfin.common.control_validator.FieldSpec` so the
    descriptor always matches control_validator.  ``changes`` is the "what
    this key does" help text from the validator.
    """

    name: str
    type: str = "str"          # "str"|"int"|"float"|"bool"|"list"|"path"
    default: Any = None
    allowed: Tuple[Any, ...] = ()   # allowed values, or empty = unconstrained
    changes: str = ""


@dataclass(frozen=True)
class ExtensionPoint:
    """Where a step may be added to a workflow, and what it receives."""

    position: str            # "before" | "after" | "in_place"
    anchor: str = ""         # step label / job name the step attaches to; "" = end
    receives: str = ""       # what the inserted step gets (geometry, gbw, state...)


# ---------------------------------------------------------------------------
# The descriptor — extends Application, no parallel type
# ---------------------------------------------------------------------------


class WorkflowDescriptor(Application):
    """A named, versioned workflow with a typed input/output contract PLUS the
    fields an agent needs to use, extend and integrate it.

    Every classic-workflow, pipeline and user workflow is described by this
    type; there is no second, parallel workflow type.  Construction from a
    serialized Application (no descriptor fields) is tolerated: the
    descriptor fields default and the workflow is still a valid Application.
    """

    def __init__(
        self,
        name: str,
        *,
        description: str = "",
        version: str = "1",
        category: str = "",
        inputs=(),
        outputs=(),
        spec: Optional[Dict[str, Any]] = None,
        kind: str = "pipeline",
        purpose: str = "",
        entry: str = "",
        control_keys=(),
        steps=(),
        artifacts=(),
        extension_points=(),
        resume: str = "",
        backend: str = "local",
    ):
        super().__init__(
            name=name,
            description=description,
            version=version,
            category=category,
            inputs=inputs,
            outputs=outputs,
            spec=spec if spec is not None else {},
        )
        self.kind = kind                      # "classic-control"|"pipeline"|"user"
        self.purpose = purpose or description  # one sentence
        self.entry = entry
        self.control_keys: Tuple[ControlKeySpec, ...] = tuple(control_keys)
        self.steps: Tuple[Dict[str, Any], ...] = tuple(steps)
        self.artifacts: Tuple[str, ...] = tuple(artifacts)
        self.extension_points: Tuple[ExtensionPoint, ...] = tuple(extension_points)
        self.resume = resume
        self.backend = backend

    # -- serialization (extends Application.to_dict/from_dict) --------------

    def to_dict(self) -> Dict[str, Any]:
        data = super().to_dict()
        data.update({
            "kind": self.kind,
            "purpose": self.purpose,
            "entry": self.entry,
            "control_keys": [_ck_to_dict(c) for c in self.control_keys],
            "steps": [dict(s) for s in self.steps],
            "artifacts": list(self.artifacts),
            "extension_points": [_ep_to_dict(e) for e in self.extension_points],
            "resume": self.resume,
            "backend": self.backend,
        })
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowDescriptor":
        # Tolerate a plain serialized Application (no descriptor fields).
        return cls(
            name=data["name"],
            description=data.get("description", ""),
            version=str(data.get("delfin_application", data.get("version", "1"))),
            category=data.get("category", ""),
            kind=data.get("kind", "pipeline"),
            purpose=data.get("purpose", ""),
            entry=data.get("entry", ""),
            control_keys=tuple(_ck_from_dict(d) for d in data.get("control_keys", [])),
            steps=tuple(dict(s) for s in data.get("steps", [])),
            artifacts=tuple(data.get("artifacts", [])),
            extension_points=tuple(
                _ep_from_dict(d) for d in data.get("extension_points", [])),
            resume=data.get("resume", ""),
            backend=data.get("backend", "local"),
            spec=data.get("spec", {}),
        )

    # -- introspection helpers -----------------------------------------------

    def control_key(self, name: str) -> Optional[ControlKeySpec]:
        return next((c for c in self.control_keys if c.name == name), None)

    def steps_as_spec(self) -> Dict[str, Any]:
        """The workflow's steps as a PipelineSpec-shaped dict, so the block list
        is something the pipeline server can render, validate or extend.

        For a classic-control workflow this is a *description* of the jobs it
        writes — the workflow itself is not reimplemented; nothing here runs.
        """
        return {
            "delfin_pipeline_spec": "1",
            "name": self.name,
            "kind": "pipeline",
            "defaults": {},
            "steps": [dict(s) for s in self.steps],
        }

    def allows(self, position: str, anchor: str = "") -> bool:
        """True if a step may be inserted at *position* after/before *anchor*."""
        for ep in self.extension_points:
            if ep.position == position and (not ep.anchor or ep.anchor == anchor):
                return True
        return False


# ---------------------------------------------------------------------------
# (de)serialization of the value types
# ---------------------------------------------------------------------------


def _ck_to_dict(c: ControlKeySpec) -> Dict[str, Any]:
    d: Dict[str, Any] = {"name": c.name, "type": c.type}
    if c.default is not None:
        d["default"] = c.default
    if c.allowed:
        d["allowed"] = list(c.allowed)
    if c.changes:
        d["changes"] = c.changes
    return d


def _ck_from_dict(d: Dict[str, Any]) -> ControlKeySpec:
    return ControlKeySpec(
        name=d["name"],
        type=d.get("type", "str"),
        default=d.get("default"),
        allowed=tuple(d["allowed"]) if d.get("allowed") else (),
        changes=d.get("changes", ""),
    )


def _ep_to_dict(e: ExtensionPoint) -> Dict[str, Any]:
    d: Dict[str, Any] = {"position": e.position}
    if e.anchor:
        d["anchor"] = e.anchor
    if e.receives:
        d["receives"] = e.receives
    return d


def _ep_from_dict(d: Dict[str, Any]) -> ExtensionPoint:
    return ExtensionPoint(
        position=d["position"],
        anchor=d.get("anchor", ""),
        receives=d.get("receives", ""),
    )


# ---------------------------------------------------------------------------
# CONTROL -key spec derivation (never drifts from control_validator)
# ---------------------------------------------------------------------------


def _infer_type(field) -> str:
    """Map a control_validator FieldSpec to a coarse type string."""
    fn = getattr(field.coerce, "__name__", "")
    if "int" in fn:
        return "int"
    if "float" in fn:
        return "float"
    if "bool" in fn or "yes_no" in fn or "on_off" in fn:
        return "bool"
    # Fall back to the default's Python type.
    if isinstance(field.default, bool):
        return "bool"
    if isinstance(field.default, int):
        return "int"
    if isinstance(field.default, float):
        return "float"
    if isinstance(field.default, (list, tuple)):
        return "list"
    return "str"


def control_key_spec(name: str) -> Optional[ControlKeySpec]:
    """The :class:`ControlKeySpec` for one CONTROL key, straight from the
    validator — None if control_validator does not know the key.

    ``changes`` comes from the validator's help text when present, else falls
    back to a neutral description so every declared key carries a reason."""
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    for f in CONTROL_FIELD_SPECS:
        if f.name == name:
            return ControlKeySpec(
                name=f.name,
                type=_infer_type(f),
                default=f.default,
                allowed=(),
                changes=f.help or (
                    f"CONTROL.txt key '{f.name}' that the workflow reads when "
                    "building its ORCA jobs"
                ),
            )
    return None


def _ck_specs(names: Iterable[str]) -> Tuple[ControlKeySpec, ...]:
    """Build ControlKeySpec tuple for the named CONTROL keys, refusing to
    invent a key control_validator does not know."""
    out: List[ControlKeySpec] = []
    for n in names:
        spec = control_key_spec(n)
        if spec is None:
            raise ValueError(
                f"workflow descriptor asks for CONTROL key {n!r} which "
                "control_validator does not know — descriptor may not invent keys")
        out.append(spec)
    return tuple(out)


# ---------------------------------------------------------------------------
# Classic workflow descriptors (describe; never reimplement)
# ---------------------------------------------------------------------------

CLASSIC_WORKFLOWS: Tuple[str, ...] = ("classic", "OCCUPIER", "ESD", "IMAG",
                                      "redox", "co2")

_CLASSIC_REGISTRY: Dict[str, WorkflowDescriptor] = {}
_CLASSIC_REGISTERED = False


def _classic_block(label: str, step: str, **kwargs: Any) -> Dict[str, Any]:
    return {"step": step, "label": label, "kwargs": kwargs}


def _build_classic_descriptors() -> List[WorkflowDescriptor]:
    """Truthful descriptors of the CONTROL.txt-driven classic workflows.

    ``steps`` names the building blocks / jobs each workflow writes,
    described as data.  Nothing here runs; the real engines (run_classic_phase,
    run_occuper_phase, run_esd_phase, the IMAG and CO2 coordinators) are
    untouched.  Extension-point anchors match the step labels below.
    """
    return [
        WorkflowDescriptor(
            "classic",
            kind="classic-control",
            purpose="Optimise, frequency-check and screen a structure through "
                    "the classic CONTROL.txt engine (method=classic).",
            entry="CONTROL.txt with method=classic; run_classic_phase dispatches it "
                  "(delfin/workflows/pipeline.py).",
            category="dft",
            control_keys=_ck_specs([
                "method", "charge", "functional", "main_basisset", "metal_basisset",
                "geom_opt", "freq_type", "PAL", "solvent", "initial_guess",
                "multiplicity_global_opt",
            ]),
            steps=[
                _classic_block("structure", "smiles_to_xyz",
                               smiles="{smiles}"),
                _classic_block("opt", "orca_opt"),
                _classic_block("freq", "orca_freq"),
            ],
            artifacts=("*.xyz", "*.gbw", "*.out", "DELFIN_Data.json",
                       "frequency_calculation_*"),
            extension_points=(
                ExtensionPoint("after", "opt", "geometry, gbw"),
                ExtensionPoint("after", "freq", "geometry, hessian"),
                ExtensionPoint("in_place", "opt", "geometry"),
            ),
            resume="Recovery via delfin/orca_recovery.py; ORCA restart files "
                   "let an interrupted opt/freq be resumed.",
            backend="local and Slurm",
        ),
        WorkflowDescriptor(
            "OCCUPIER",
            kind="classic-control",
            purpose="Optimise every plausible spin state of a metal complex and "
                    "select the preferred one (method=OCCUPIER).",
            entry="CONTROL.txt gate OCCUPIER or method=OCCUPIER; run_occuper_phase "
                  "dispatches it (delfin/workflows/pipeline.py).",
            category="dft",
            control_keys=_ck_specs([
                "method", "charge", "functional", "main_basisset", "metal_basisset",
                "multiplicity_global_opt", "OCCUPIER_method",
                "frequency_calculation_OCCUPIER", "geom_opt_OCCUPIER",
                "OCCUPIER_compare", "solvent", "PAL",
            ]),
            steps=[
                _classic_block("structure", "smiles_to_xyz",
                               smiles="{smiles}"),
                _classic_block("spin_states", "occupier"),
            ],
            artifacts=("input.xyz", "*.gbw", "*.out", "OCCUPIER.txt",
                       "input_*_OCCUPIER.xyz/.gbw", "DELFIN_Data.json"),
            extension_points=(
                ExtensionPoint("after", "spin_states", "chosen state, geometry"),
                ExtensionPoint("after", "structure", "geometry"),
            ),
            resume="Per-(frame,mult) folders let a run continue past a failed "
                   "candidate; OCCUPIER.txt records the preferred index.",
            backend="local and Slurm",
        ),
        WorkflowDescriptor(
            "ESD",
            kind="classic-control",
            purpose="Excited-state dynamics screening: TDDFT per state, then "
                    "emission/ISC/IC rates (ESD_modul=yes).",
            entry="CONTROL.txt with ESD_modul=yes; run_esd_phase dispatches it "
                  "(delfin/workflows/pipeline.py).",
            category="dft",
            control_keys=_ck_specs([
                "method", "charge", "functional", "main_basisset", "metal_basisset",
                "ESD_modus", "ESD_T1_opt", "TDDFT_nroots", "TDDFT_maxiter",
                "TDDFT_TDA", "TDDFT_followiroot", "TDDFT_SOC", "states",
                "ICs", "ISCs", "emission_rates", "solvent",
            ]),
            steps=[
                _classic_block("structure", "smiles_to_xyz",
                               smiles="{smiles}"),
                _classic_block("opt", "orca_opt"),
                _classic_block("excited", "orca_tddft"),
                _classic_block("esd", "esd"),
            ],
            artifacts=("*.out", "*.gbw", "emission_rates.json",
                       "profile_diagnosis.json", "DELFIN_Data.json"),
            extension_points=(
                ExtensionPoint("after", "excited", "gbw, excited states"),
                ExtensionPoint("after", "opt", "geometry, gbw"),
            ),
            resume="Per-state TDDFT folders resume independently on a partial "
                   "failure.",
            backend="local and Slurm",
        ),
        WorkflowDescriptor(
            "IMAG",
            kind="classic-control",
            purpose="Drive an imaginary-frequency saddle to a minimum by "
                    "reoptimising along the imaginary mode.",
            entry="CONTROL.txt gate IMAG_option / IMAG_scope; the IMAG engine "
                  "(delfin/cli_imag.py) dispatches it.",
            category="dft",
            control_keys=_ck_specs([
                "method", "charge", "functional", "main_basisset", "metal_basisset",
                "IMAG_scope", "IMAG_sp_energy_window", "IMAG_option",
                "IMAG_max_rounds", "allow_imaginary_freq", "solvent",
            ]),
            steps=[
                _classic_block("structure", "smiles_to_xyz",
                               smiles="{smiles}"),
                _classic_block("saddle_opt", "orca_opt"),
                _classic_block("freq", "orca_freq"),
                _classic_block("imag_fix", "imag_fix"),
            ],
            artifacts=("*.out", "*.gbw", "*.hess", "IMAG_*",
                       "DELFIN_Data.json"),
            extension_points=(
                ExtensionPoint("after", "imag_fix", "geometry, hessian"),
                ExtensionPoint("after", "freq", "hessian, imaginary modes"),
            ),
            resume="Reopt along the imaginary mode resumes from the recorded "
                   "saddle geometry; IMAG_max_rounds bounds the loop.",
            backend="local and Slurm",
        ),
        WorkflowDescriptor(
            "redox",
            kind="classic-control",
            purpose="Compute a redox potential: optimise neutral, then "
                    "single-point the oxidised/reduced states.",
            entry="CONTROL.txt with calc_potential_method; the classic redox "
                  "engine builds the ox/red steps.",
            category="dft",
            control_keys=_ck_specs([
                "method", "charge", "functional", "main_basisset", "metal_basisset",
                "solvent", "calc_potential_method",
            ]),
            steps=[
                _classic_block("structure", "smiles_to_xyz",
                               smiles="{smiles}"),
                _classic_block("opt_neutral", "orca_opt"),
                _classic_block("ox", "orca_sp"),
                _classic_block("red", "orca_sp"),
            ],
            artifacts=("*.out", "*.gbw", "e_oxidation.txt", "e_reduction.txt",
                       "DELFIN_Data.json"),
            extension_points=(
                ExtensionPoint("after", "opt_neutral", "geometry, gbw"),
                ExtensionPoint("after", "red", "geometry, gbw"),
            ),
            resume="Ox/red single points run independently and resume on "
                   "partial failure.",
            backend="local and Slurm",
        ),
        WorkflowDescriptor(
            "co2",
            kind="classic-control",
            purpose="Coordination workflow for CO2 adducts: chain several ORCA "
                    "runs across species (co2_coordination).",
            entry="CONTROL.txt with co2_coordination / co2_species_delta; the "
                  "CO2 coordinator (delfin/co2/CO2_Coordinator6.py) drives it.",
            category="dft",
            control_keys=_ck_specs([
                "method", "charge", "functional", "main_basisset", "metal_basisset",
                "solvent", "co2_coordination", "co2_species_delta",
            ]),
            steps=[
                _classic_block("structure", "smiles_to_xyz",
                               smiles="{smiles}"),
                _classic_block("opt", "orca_opt"),
                _classic_block("species_runs", "orca_sp"),
            ],
            artifacts=("*.out", "*.gbw", "chain_state.json", "DELFIN_Data.json"),
            extension_points=(
                ExtensionPoint("after", "opt", "geometry, gbw"),
                ExtensionPoint("after", "species_runs", "geometry, gbw"),
            ),
            resume="The coordinator chains runs with per-species checkpoints "
                   "and resumes from the last completed species.",
            backend="local and Slurm",
        ),
    ]


def register_classic_workflows() -> Dict[str, WorkflowDescriptor]:
    """Register (once, idempotently) a descriptor per classic workflow.

    Registration is purely descriptive: it builds the catalog entries; it
    never calls, replaces or otherwise perturbs the classic engines.  Returns
    the registry mapping workflow name -> descriptor.
    """
    global _CLASSIC_REGISTERED
    if not _CLASSIC_REGISTERED:
        for wf in _build_classic_descriptors():
            _CLASSIC_REGISTRY[wf.name] = wf
        _CLASSIC_REGISTERED = True
    return dict(_CLASSIC_REGISTRY)


def workflow_for(name: str) -> Optional[WorkflowDescriptor]:
    """The descriptor for a classic workflow, or None if unknown.

    Registry access does not trigger registration side effects beyond the
    classic descriptors (listed on demand) so a plain ``list_workflows`` call
    stays cheap and deterministic.
    """
    return _CLASSIC_REGISTRY.get(name)


# ---------------------------------------------------------------------------
# Catalog: classic + pipeline + user listed together
# ---------------------------------------------------------------------------


def _user_descriptor(entry: Dict[str, Any], root: Optional[Path]) -> WorkflowDescriptor:
    """Wrap one saved user workflow as a descriptor of kind 'user'."""
    spec = entry
    return WorkflowDescriptor(
        entry.get("name", "user_workflow"),
        kind="user",
        purpose=spec.get("description", ""),
        entry=f"saved user workflow ({spec.get('name', '?')})",
        spec={"name": spec.get("name", ""), "steps": spec.get("steps", []),
              "kind": "pipeline", "defaults": spec.get("defaults", {})},
        steps=spec.get("steps", []),
        backend="local and Slurm",
    )


def list_workflows(*, root: Optional[Path] = None) -> List[WorkflowDescriptor]:
    """Every workflow the platform knows, in the order classic-control,
    pipeline applications, then saved user workflows.

    One contract for the whole repertoire: each entry is a
    :class:`WorkflowDescriptor` (itself an Application) carrying kind,
    purpose, entry, control keys, steps, artifacts and extension points.
    """
    register_classic_workflows()
    out: List[WorkflowDescriptor] = []
    # 1. classic CONTROL.txt workflows.
    out.extend(_CLASSIC_REGISTRY.values())
    # 2. pipeline applications (kind "pipeline" by default).
    from delfin.tools._application import list_applications

    for app in list_applications().values():
        # A classic descriptor may share a short name, keep the pipeline one.
        if isinstance(app, WorkflowDescriptor):
            continue
        out.append(_as_descriptor(app))
    # 3. saved user workflows.
    from delfin.tools import user_workflows as uw

    for entry in uw.list_workflows(root=root):
        out.append(_user_descriptor(entry, root))
    return out


def _as_descriptor(app: Application) -> WorkflowDescriptor:
    """Promote a plain Application to a WorkflowDescriptor (kind pipeline).

    This is a view, not a new workflow type: the descriptor IS the
    Application, extended with the one-contract fields.
    """
    wf = WorkflowDescriptor(
        name=app.name,
        description=app.description,
        version=app.version,
        category=app.category,
        inputs=app.inputs,
        outputs=app.outputs,
        spec=app.spec,
        kind="pipeline",
        purpose=app.description,
        entry=f"run_application('{app.name}')",
        backend="local and Slurm",
    )
    return wf


# ---------------------------------------------------------------------------
# Extension: "add step X after step Y of workflow W"
# ---------------------------------------------------------------------------


def add_step(
    descriptor: WorkflowDescriptor,
    new_step: Dict[str, Any],
    *,
    after: Optional[str] = None,
    before: Optional[str] = None,
    replace: Optional[str] = None,
) -> Dict[str, Any]:
    """Produce a NEW checked pipeline spec by inserting *new_step* into
    workflow *descriptor*.

    ``after`` / ``before`` / ``replace`` select the anchor by step label
    (falling back to step name).  Exactly one of them must be given.  The
    original workflow is never mutated: a new ``PipelineSpec`` is returned
    whose steps are W's steps unchanged plus the inserted step, ready to be
    passed through the existing checks (``user_workflows.check_spec``).

    Raises
    ------
    ValueError
        when zero or more than one of after/before/replace is given.
    KeyError
        when the anchor step does not exist in the workflow.
    """
    modes = [m for m in (after, before, replace) if m is not None]
    if len(modes) != 1:
        raise ValueError("give exactly one of after= / before= / replace=")

    # Start from the workflow's steps (jobs for classic, real spec for others).
    base = descriptor.steps_as_spec()
    steps = [dict(s) for s in base["steps"]]

    def _key(s: Dict[str, Any]) -> str:
        return s.get("label") or s.get("step", "")

    labels = [_key(s) for s in steps]
    if replace is not None:
        if replace not in labels:
            raise KeyError(f"no step matched replace={replace!r}")
        out = []
        for s in steps:
            if _key(s) == replace:
                out.append(dict(new_step))
            else:
                out.append(s)
        steps = out
    elif after is not None or before is not None:
        anchor = after if after is not None else before
        if anchor not in labels:
            raise KeyError(f"no step matched anchor={anchor!r}")
        inserted = False
        out: List[Dict[str, Any]] = []
        for s in steps:
            if before is not None and _key(s) == before and not inserted:
                out.append(dict(new_step))
                inserted = True
            out.append(s)
            if after is not None and _key(s) == after:
                out.append(dict(new_step))
                inserted = True
        steps = out

    name = _sluggish(f"{descriptor.name}_{_new_label(new_step)}")
    spec = {
        "delfin_pipeline_spec": "1",
        "name": name,
        "kind": "pipeline",
        "defaults": dict(base.get("defaults", {})),
        "steps": steps,
    }
    return spec


def _new_label(step: Dict[str, Any]) -> str:
    return step.get("label") or step.get("step", "step")


def _sluggish(name: str) -> str:
    import re as _re

    s = _re.sub(r"[^a-zA-Z0-9_-]+", "_", name).strip("_").lower()
    return s or "workflow"
