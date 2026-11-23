"""Excited-state decision rules for the scientific lint.

Called from ``delfin.tools._science._lint_step``, which needs the step
NAME as well as its kwargs -- these rules are about what a TD-DFT step
decided, and a step's identity is part of that. This is not a second
linter; it is where these rules live so the single entry point stays
``scientific_lint``.

Both rules say the same thing in two places: **an engine default is not
a decision.** A spec that leaves a choice to the program records no
choice, and a result whose level of theory was picked by a default
cannot be defended by the person who published it. The existing grader
turns that into pressure in the right direction: a ``warning`` must be
justified in the answer text (``agent/workflow_grader``), so a step that
decided deliberately and says why passes, and one that drifted into a
default does not.

Sources:
- Charge-transfer excitations and non-local exchange: Dreuw, Weisman and
  Head-Gordon, J. Chem. Phys. 2003, 119, 2943-2946,
  DOI 10.1063/1.1590951 -- the attraction between separated charges comes
  from non-local exchange, so a global hybrid cannot describe a
  long-range CT state at any basis.
- The Tamm-Dancoff approximation: Hirata and Head-Gordon,
  Chem. Phys. Lett. 1999, 314, 291-299,
  DOI 10.1016/S0009-2614(99)01149-5 -- TDA drops the de-excitation
  amplitudes, which changes the excitation energies it returns.
- That TDA is ORCA's own default for TD-DFT: DELFIN's knowledge
  dictionary, ``orca_knowledge/dictionary/
  properties_spectroscopy_excited.json`` under ``TDA`` ("TDA is the
  default for TD-DFT in ORCA and can be switched off to recover the
  full (RPA-like) response"). Read from the tree, not asserted here.
"""

from __future__ import annotations

from typing import Any, Dict, List, TYPE_CHECKING

if TYPE_CHECKING:  # annotation only -- real import is lazy (parent imports us)
    from delfin.tools._science import ScienceFinding

#: The steps these rules apply to: the ones that solve a TD-DFT response
#: equation. Named rather than pattern-matched, so a new excited-state
#: step is added here deliberately and not caught by accident.
_EXCITED_STEPS = ("orca_tddft", "esd", "esd_opi")

#: A declared excited state with charge-transfer character. The field is
#: ``excited_state``, the same one the P1 rules read for Rydberg and
#: dipole-bound states -- one vocabulary, not two.
_CT_MARKERS = ("charge_transfer", "charge-transfer", "ct", "push_pull",
               "push-pull", "donor_acceptor", "donor-acceptor",
               "intramolecular_ct", "mlct", "lmct")

#: Functionals whose exchange is long-range corrected or fully non-local.
#: A substring match on the method name, because ORCA accepts several
#: spellings of each (wB97X-D3, CAM-B3LYP, LC-BLYP, ...).
_LONG_RANGE_MARKERS = ("wb97", "cam-", "camb3", "lc-", "lrc-", "hse",
                       "m11", "mn15", "w b97", "ωb97", "tuned")

#: Range-separated functionals ORCA knows, offered by name when a global
#: hybrid is used on a CT state. Kept short on purpose: a list of twenty
#: is not advice.
_SUGGESTED = "wB97X-D4, CAM-B3LYP or LC-BLYP"


def _is_excited_step(step_name: str) -> bool:
    return str(step_name or "").strip().lower() in _EXCITED_STEPS


def _declares_charge_transfer(kwargs: Dict[str, Any]) -> bool:
    state = str(kwargs.get("excited_state") or "").strip().lower()
    if not state:
        return False
    # Whole-token match for the short marker: "ct" must not fire on
    # "octahedral" or "spectrum".
    tokens = {t for t in state.replace("-", "_").split("_") if t}
    if "ct" in tokens:
        return True
    return any(m in state for m in _CT_MARKERS if m != "ct")


def _is_long_range_corrected(method: str) -> bool:
    m = str(method or "").strip().lower()
    return any(marker in m for marker in _LONG_RANGE_MARKERS)


def lint_excited(step_name: str, label: str,
                 kwargs: Dict[str, Any]) -> List[ScienceFinding]:
    """The excited-state rungs: CT exchange, and TDA stated not inherited."""
    from delfin.tools._science import ScienceFinding, _is_placeholder
    out: List[ScienceFinding] = []
    if not _is_excited_step(step_name):
        return out

    # Rung 1: a long-range CT state needs non-local exchange. This is not
    # a basis-set question and cannot be fixed by one: the global hybrid
    # has no term that holds the separated charges together.
    method = kwargs.get("method")
    if (_declares_charge_transfer(kwargs) and isinstance(method, str)
            and method and not _is_placeholder(method)
            and not _is_long_range_corrected(method)):
        out.append(ScienceFinding(
            label, "warning", "ct-needs-long-range-exchange",
            f"excited_state '{kwargs.get('excited_state')}' has "
            f"charge-transfer character and '{method}' is a global hybrid: "
            "the attraction between the separated charges comes from "
            "non-local exchange, so this underestimates the excitation "
            f"energy at any basis set — use {_SUGGESTED}, or say in the "
            "answer why a global hybrid is the right choice here"))

    # Rung 2: whether TDA was used is a decision and belongs in the spec.
    # ORCA applies it by default, so a step that is silent about it has
    # not chosen the full response -- it has inherited an approximation
    # that moves the numbers it reports.
    if "tda" not in kwargs:
        out.append(ScienceFinding(
            label, "warning", "tda-not-stated",
            "this step does not say whether the Tamm-Dancoff "
            "approximation is used, and TD-DFT in ORCA applies it by "
            "default: TDA drops the de-excitation amplitudes and so "
            "changes the excitation energies reported. Set tda=true or "
            "tda=false, or name the choice in the answer — a default is "
            "not a decision"))
    return out


__all__ = ["lint_excited"]
