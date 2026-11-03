"""Workflow builder — request text to a checked, savable workflow spec.

The agent-side path from "build me a workflow for X" to a draft spec that
has PASSED the mandatory checks, in this order, before anything is
submitted or saved:

1. resolve_spec   — fill contract defaults, show auto-wiring
2. validate_spec  — static wiring/contract validation
3. scientific_lint— physics plausibility (an ``error`` blocks)
4. dry run        — Pipeline.validate with the geometry flag

NOTE (R1, reviewer): the "dry run" stage is NOT an execution —
``Pipeline.dry_run`` (delfin/tools/pipeline.py:1648) is an alias of
``Pipeline.validate`` (:1629). Stage 4 therefore re-validates with the
geometry flag satisfied; the chain's distinct checks are
resolve -> validate -> lint -> cost. It is kept as an explicit stage so
a geometry-flagged validation always happens before submit, but it adds
no execution. The skill and docs say the same.

5. cost estimate  — cluster_efficiency.advise per step (DELFIN's own sizing)

A failed check returns the findings verbatim (the fix is named, never
applied silently); the result never claims a submit happened.

Sources: a built-in template by name (``delfin.tools.templates``), a
contributed workflow spec, a provided spec dict, or — with none given —
a small request-driven draft from the catalog. Extend/modify always
returns a NEW spec (see ``user_workflows.extend_template``); the source
is never mutated.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------
# Request vocabulary: which template a request maps to. The names come
# from delfin/tools/templates/__init__.py, never invented here.
# ---------------------------------------------------------------------

_SMILES_RE = re.compile(r"SMILES\s+([^\s()]+)", re.IGNORECASE)


def _known_templates() -> Dict[str, str]:
    from delfin.tools import templates

    names = list(getattr(templates, "__all__", []))
    hints = {
        "classic_opt_freq": "opt+freq",
        "imag_elimination": "imag",
        "redox_potential": "redox",
        "conformer_screening": "conformer",
        "multi_level_opt": "multi-level",
        "esd_states": "esd excited",
        "opi_opt_freq_dlpno": "dlpno",
    }
    return {n: hints.get(n, n.replace("_", " ")) for n in names}


def _pick_template(request: str) -> Optional[str]:
    text = request.lower()
    best = None
    # The rightmost matching hint wins: the head of a request names the
    # generic task ("opt+freq" -> classic), the tail names the
    # distinguishing step ("...then DLPNO" -> the OPI chain).
    for name, hint in _known_templates().items():
        pos = text.rfind(hint)
        if pos >= 0 and (best is None or pos > best[1]):
            best = (name, pos)
    return best[0] if best else None


def _smiles_of(request: str) -> Optional[str]:
    m = _SMILES_RE.search(request)
    return m.group(1).strip(",.;") if m else None


# ---------------------------------------------------------------------
# Step reasons: one line per step, grounded in what the step does.
# ---------------------------------------------------------------------

_STEP_REASONS: Dict[str, str] = {
    "smiles_to_xyz": "3D geometry from SMILES (RDKit embed + xTB clean-up)",
    "xtb_opt": "fast xTB pre-optimization brings the structure near the minimum",
    "orca_opt": "DFT geometry optimization to the converged minimum",
    "orca_freq": "frequencies verify the minimum and give thermochemistry",
    "orca_sp": "single-point energy at the final geometry",
    "orca_tddft": "excited states / absorption spectrum at the ground-state geometry",
    "xtb_sp": "fast xTB single point for a first energy estimate",
    "orca_opt_opi": "OPI-built DFT optimization; leaves the optimized geometry and .gbw for the next step",
    "orca_freq_opi": "OPI-built frequencies verify the minimum and give the thermal correction",
    "orca_dlpno_opi": "DLPNO-CCSD(T) single point on the chain's geometry, started from its orbitals",
    "orca_opt_freq_opi": "OPI-built combined optimization + frequency",
}


def _steps_reasons(spec: Dict[str, Any]) -> List[Dict[str, str]]:
    out = []
    for s in spec.get("steps", []):
        name = s.get("step", "?")
        out.append({
            "step": name,
            "label": s.get("label", name),
            "reason": _STEP_REASONS.get(
                name, f"step '{name}' from the DELFIN catalog (see describe_capability)"),
        })
    return out


# ---------------------------------------------------------------------
# The mandatory check chain.
# ---------------------------------------------------------------------

_CHECK_ORDER = ["resolve", "validate", "lint", "dry_run", "cost"]


def run_checks(spec: Dict[str, Any], *, system_size: Optional[int] = None,
               geometry: bool = True) -> Dict[str, Any]:
    """Run the five mandatory checks over a draft spec, in order.

    Returns ``{"ok", "order", "errors", "lint", "report", "cost"}``.
    ``ok`` is False if ANY check fails; nothing is executed — a dry run is
    static validation, and the cost estimate is DELFIN's own
    :func:`cluster_efficiency.advise` (no submission of any kind).
    """
    from delfin.tools import platform
    from delfin.tools import user_workflows as uw
    from delfin.agent import cluster_efficiency

    errors: List[str] = []
    lint: List[Dict[str, Any]] = []
    report = None
    resolved = None

    resolved = platform.resolve_spec(spec)
    report = platform.validate_spec(spec, geometry=geometry)
    if report is not None and not report.ok:
        errors.extend(report.errors)
    lint = platform.scientific_lint(spec)
    if not uw.lint_ok(lint):
        errors.extend(f"[lint:{f['rule']}] {f['step']}: {f['message']}"
                      for f in lint if f.get("level") == "error")

    # Check 4: a real dry run — re-build the pipeline from the spec and
    # validate it with the geometry flag (Pipeline.dry_run IS validate,
    # pipeline.py:1648; it executes nothing). Distinct from check 2 only
    # in that it goes through the pipeline object, catching spec shapes
    # validate_spec's application path may accept but the pipeline cannot.
    from delfin.tools._serialize import from_dict as _from_dict
    dry_run = _from_dict(spec).validate(geometry=geometry)
    if not dry_run.ok:
        errors.extend(f"[dry_run] {e}" for e in dry_run.errors)

    # Check 5: cost estimate. Provenance: without a caller-given size,
    # DELFIN's default (10 atoms) is ASSUMED and each entry says so —
    # a number without provenance is marked as an estimate.
    assumed = system_size is None
    cost: List[Dict[str, Any]] = []
    for entry in spec.get("steps", []):
        size = system_size if system_size else 10
        suggestion = cluster_efficiency.advise(entry.get("step", "?"), size)
        if assumed:
            suggestion["assumed_size"] = True
            suggestion["note"] = (
                "system_size not given; a default of 10 atoms was "
                "ASSUMED for sizing — give system_size for a measured "
                "estimate")
        cost.append(suggestion)

    return {
        "ok": not errors,
        "order": list(_CHECK_ORDER),
        "errors": errors,
        "lint": lint,
        "report": report,
        "resolved": resolved,
        "dry_run": dry_run,
        "cost": cost,
    }


# ---------------------------------------------------------------------
# Draft construction.
# ---------------------------------------------------------------------


def _draft_spec(request: str) -> Dict[str, Any]:
    import delfin.tools.templates as templates

    name = _pick_template(request) or "classic_opt_freq"
    tpl = getattr(templates, name, None)
    if tpl is None or not hasattr(tpl, "to_dict"):
        raise ValueError(
            f"unknown template {name!r}: pick one of "
            f"{sorted(getattr(templates, '__all__', []))}")
    spec = tpl.to_dict()
    spec.pop("delfin_pipeline_spec", None)
    spec = dict(spec)
    spec["name"] = f"draft_{name}"
    smiles = _smiles_of(request)
    if smiles and spec.get("steps"):
        first = dict(spec["steps"][0])
        if first.get("step") == "smiles_to_xyz":
            first["kwargs"] = dict(first.get("kwargs", {}), smiles=smiles)
            spec["steps"] = [first] + list(spec["steps"][1:])
    return spec


def build_workflow(request: str, *, base_template: Optional[str] = None,
                   base_spec: Optional[Dict[str, Any]] = None,
                   edits: Optional[Dict[str, Any]] = None,
                   system_size: Optional[int] = None,
                   geometry: bool = True) -> Dict[str, Any]:
    """Request text -> draft spec + mandatory checks + cost estimate.

    Precedence: ``base_spec`` (a provided spec dict) > ``base_template``
    (a named built-in template) > a request-driven draft from DELFIN's
    templates. ``edits`` (add/remove/replace/overrides) go through
    :func:`delfin.tools.user_workflows.extend_template`, which returns a
    NEW spec and leaves the source untouched.

    The result carries ``source`` (where the spec came from), ``spec``,
    ``steps_reasons`` (one line per step), ``checks`` and ``submitted``
    — always False here: submitting is a separate, explicit act that the
    caller only does after ``checks["ok"]`` is True.
    """
    from delfin.tools import user_workflows as uw

    if base_spec is not None:
        source = "provided spec"
        spec = dict(base_spec)
    elif base_template is not None:
        source = f"template:{base_template}"
        spec = uw.extend_template(base_template, new_name=f"draft_{base_template}",
                                  **(edits or {}))
    else:
        # A request that names a template maps to it explicitly; the
        # draft source records which template the request picked.
        picked = _pick_template(request)
        source = f"template:{picked}" if picked else "request draft"
        spec = _draft_spec(request)

    spec.setdefault("description", request.strip()[:200])
    checks = run_checks(spec, system_size=system_size, geometry=geometry)
    return {
        "request": request,
        "source": source,
        "spec": spec,
        "steps_reasons": _steps_reasons(spec),
        "checks": checks,
        "cost": checks["cost"],
        "submitted": False,
    }


__all__ = ["build_workflow", "run_checks"]
