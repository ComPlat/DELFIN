"""SV1 solvation realism rules for the scientific lint.

Wired into delfin/tools/_science.py by the operator's SV1 hook (the
WB1 extension pattern, one call per rule family). Two rules:

- ``lint_solvator_step`` (step level): SOLVATOR is an xTB feature —
  a ``method`` on an ``xtb_solvator`` step that is not GFN0/1/2-xTB is
  a physics ERROR (the ORCA input would be invalid, surfacing only
  after queue time). Source: ORCA 6 manual, %SOLVATOR block.
- ``lint_solvation_chain`` (spec level): explicit and implicit
  solvation on the same species double-count solvation. If a step
  AFTER the explicit shell (xtb_solvator / packmol_solvate) carries a
  DFT solvent model (CPCM/SMD/ALPB via solvent_model), the free-energy
  correction is counted twice. The finding is a WARNING because a
  later single-point at a different solvent treatment can be
  deliberate — but it must be justified, not silent.
"""

from __future__ import annotations

from typing import Any, Dict, List

from delfin.tools._science import ScienceFinding, _is_placeholder

_XTB_HAMILTONIANS = frozenset({"gfn0-xtb", "gfn1-xtb", "gfn2-xtb"})
_SHELL_STEPS = frozenset({"xtb_solvator", "packmol_solvate"})
_DFT_SOLVENT_MODELS = frozenset({"cpcm", "cpcm(w)", "smd", "smd(w)"})


def lint_solvator_step(label: str, step_name: str,
                       kwargs: Dict[str, Any]) -> List[ScienceFinding]:
    """Step-level rule: the SOLVATOR method must be an xTB Hamiltonian."""
    out: List[ScienceFinding] = []
    if step_name != "xtb_solvator":
        return out
    method = kwargs.get("method")
    if method is None or _is_placeholder(method):
        return out
    if str(method).strip().lower() not in _XTB_HAMILTONIANS:
        out.append(ScienceFinding(
            label, "error", "xtb-solvator-method-not-xtb",
            f"method {method!r} is not an xTB Hamiltonian; %SOLVATOR "
            "(ORCA manual, %SOLVATOR block) is an xTB feature and supports "
            "GFN0-xTB, GFN1-xTB or GFN2-xTB — a DFT method here produces an "
            "invalid ORCA input"))
    return out


def lint_solvation_chain(steps: List[Dict[str, Any]]) -> List[ScienceFinding]:
    """Spec-level rule: no implicit model on a step after the explicit shell."""
    out: List[ScienceFinding] = []
    shell_seen = ""
    for st in steps:
        if not isinstance(st, dict):
            continue
        name = str(st.get("step") or st.get("name") or "").strip().lower()
        label = str(st.get("label") or name or "?")
        kw = st.get("kwargs", {}) or {}
        model = kw.get("solvent_model")
        if (name not in _SHELL_STEPS and shell_seen
                and model is not None and not _is_placeholder(model)
                and str(model).strip().lower() in _DFT_SOLVENT_MODELS):
            out.append(ScienceFinding(
                label, "warning",
                "explicit-and-implicit-solvation-double-count",
                f"step '{label}' applies an implicit solvent model "
                f"({model}) after the explicit shell built by "
                f"'{shell_seen}' — the solvation free energy would be "
                "counted twice; either drop the implicit model on the "
                "explicitly solvated species or justify the combination"))
        if name in _SHELL_STEPS and not shell_seen:
            shell_seen = label
    return out


__all__ = ["lint_solvator_step", "lint_solvation_chain"]
