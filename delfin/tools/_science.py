"""Scientific plausibility lint for pipeline specs — the domain layer.

Static contract validation (``_validate``) checks that params/ports are wired and
values respect declared enums.  This module adds *physics* sanity: rules a
chemist would catch on review but that are invisible to a structural validator —
a composite method paired with a redundant basis/dispersion, a relativistic
Hamiltonian without a relativistic basis, a solvation model with no solvent, an
impossible multiplicity.

Findings are advisory by default (``warning``); only genuinely invalid input
(e.g. multiplicity < 1) is an ``error``.  Each finding names the offending step
and the concrete fix, so an agent (or a person) can correct the spec before
spending compute.  Complements the agent framework's ``/check`` critic; this is
the platform-side, purely static counterpart.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

# Grimme "3c" composite methods: each ships with its own basis set AND its own
# dispersion/geometry corrections, so a separately specified basis/dispersion is
# redundant or silently ignored.
_COMPOSITE_METHODS = frozenset({
    "hf-3c", "pbeh-3c", "b97-3c", "r2scan-3c", "b3lyp-3c", "pbe-3c", "tpss-3c",
    "wb97x-3c",
})
# Relativistic Hamiltonians that require a relativistically recontracted basis.
_RELATIVISTIC = frozenset({"zora", "dkh", "dkh2", "x2c"})
# Markers of a relativistic basis set name.
_REL_BASIS_MARKERS = ("zora", "dkh", "x2c", "sarc")


@dataclass(frozen=True)
class ScienceFinding:
    """One scientific-plausibility finding for a step."""

    step: str
    level: str          # "warning" | "error"
    rule: str
    message: str


def _is_placeholder(v: Any) -> bool:
    return isinstance(v, str) and "{" in v and "}" in v


def _lint_step(label: str, kwargs: Dict[str, Any]) -> List[ScienceFinding]:
    out: List[ScienceFinding] = []
    method = str(kwargs.get("method", "")).lower()
    basis = str(kwargs.get("basis", "")).lower()
    relativity = str(kwargs.get("relativity", "")).lower()

    composite = method in _COMPOSITE_METHODS and not _is_placeholder(kwargs.get("method"))

    # 1. composite method + an explicit basis (the basis is part of the method)
    if composite and basis and not _is_placeholder(kwargs.get("basis")):
        out.append(ScienceFinding(
            label, "warning", "composite-with-basis",
            f"method '{method}' is a composite that defines its own basis; the "
            f"explicit basis '{kwargs.get('basis')}' is ignored — drop it"))

    # 2. composite method + dispersion (dispersion is already included)
    if composite and kwargs.get("dispersion") and not _is_placeholder(kwargs.get("dispersion")):
        out.append(ScienceFinding(
            label, "warning", "composite-with-dispersion",
            f"method '{method}' already includes dispersion; the explicit "
            f"'{kwargs.get('dispersion')}' is redundant — drop it"))

    # 3. relativistic Hamiltonian with a non-relativistic basis
    if (relativity in _RELATIVISTIC and basis and not _is_placeholder(kwargs.get("basis"))
            and not any(m in basis for m in _REL_BASIS_MARKERS)):
        out.append(ScienceFinding(
            label, "warning", "relativity-without-rel-basis",
            f"relativity '{relativity}' needs a relativistically recontracted "
            f"basis (e.g. ZORA-def2-…/DKH-def2-…/SARC-…), not '{kwargs.get('basis')}'"))

    # 4. a solvation model set but no solvent named (no effect)
    if (kwargs.get("solvent_model") and not kwargs.get("solvent")
            and not _is_placeholder(kwargs.get("solvent_model"))):
        out.append(ScienceFinding(
            label, "warning", "solvent-model-without-solvent",
            f"solvent_model '{kwargs.get('solvent_model')}' has no effect without "
            "a 'solvent' — set a solvent or drop the model"))

    # 5. an impossible multiplicity
    mult = kwargs.get("mult")
    if isinstance(mult, bool):
        mult = None
    if isinstance(mult, int) and mult < 1:
        out.append(ScienceFinding(
            label, "error", "invalid-multiplicity",
            f"multiplicity {mult} is invalid; spin multiplicity (2S+1) must be ≥ 1"))

    return out


def scientific_lint(spec: Dict[str, Any]) -> List[ScienceFinding]:
    """Run the scientific-plausibility rules over a pipeline/application spec.

    Accepts a serialized pipeline spec (``{"steps": [...]}``) or an application
    dict (``{"spec": {...}}``); walks trunk + branch steps and applies the rules
    to each step's concrete kwargs (placeholder values are skipped).
    """
    inner = spec.get("spec") if "spec" in spec and "steps" not in spec else spec
    findings: List[ScienceFinding] = []

    def _walk(steps):
        for st in steps or []:
            label = st.get("label") or st.get("step") or "?"
            findings.extend(_lint_step(label, st.get("kwargs", {}) or {}))

    _walk(inner.get("steps"))
    for bname, bsteps in (inner.get("branches") or {}).items():
        # branches serialize either as a list of steps or as a nested spec dict
        steps = bsteps.get("steps") if isinstance(bsteps, dict) else bsteps
        _walk(steps)
    return findings


__all__ = ["ScienceFinding", "scientific_lint"]
