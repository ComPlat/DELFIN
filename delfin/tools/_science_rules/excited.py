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

from dataclasses import dataclass
from typing import Any, Dict, List, TYPE_CHECKING

if TYPE_CHECKING:  # annotation only -- real import is lazy (parent imports us)
    from delfin.tools._science import ScienceFinding


@dataclass(frozen=True)
class Rung:
    """One decision of the ladder, declared once.

    Two consumers read this: the lint rules below, which judge a spec
    against it, and ``ladder_text``, which renders it into the prompt the
    model reads BEFORE it builds a spec. One declaration, because the
    alternative was measured: the subagent roles were a registry and a
    hand-written table, and six of thirteen roles ended up callable and
    unnamed. A ladder the model is graded against and a ladder it is told
    about must be the same ladder.

    Why the prompt side matters at all, and not only the grading: on 558
    published vertical transition energies an unguided agent reached
    24.9 % coverage, and the same agent given the decision sequence up
    front reached 44.1 %. Grading alone does not produce that; it only
    says afterwards that the decision was missing.
    """
    code: str      #: the lint finding's rule name, so both sides match
    decide: str    #: what to decide, imperative, one line
    because: str   #: why it changes the answer
    closes: str    #: what in the spec closes it


#: The ladder for an excited-state calculation, in the order the
#: decisions actually bind: the functional governs whether the state can
#: be described at all, the response approximation moves the number it
#: returns.
RUNGS: tuple[Rung, ...] = (
    Rung(code="ct-needs-long-range-exchange",
         decide="choose the functional against the character of the state",
         because=("a global hybrid cannot hold a long-range "
                  "charge-transfer state together at any basis set, "
                  "because the attraction between the separated charges "
                  "comes from non-local exchange"),
         closes=("a range-separated functional (wB97X-D4, CAM-B3LYP, "
                 "LC-BLYP), or one sentence in the answer saying why a "
                 "global hybrid is right here")),
    Rung(code="tda-not-stated",
         decide="say whether the Tamm-Dancoff approximation is used",
         because=("ORCA applies it to TD-DFT by default, and it drops the "
                  "de-excitation amplitudes, which moves the excitation "
                  "energies you report"),
         closes="tda=true or tda=false in the step's kwargs"),
)


def _rung(code: str) -> Rung:
    """The declared rung for a lint code; a KeyError here means the two
    sides have drifted, which is the thing this file exists to prevent."""
    for rung in RUNGS:
        if rung.code == code:
            return rung
    raise KeyError(f"no declared rung for {code!r}")


def ladder_text() -> str:
    """The ladder as the model is told it, before it builds a spec."""
    lines = ["**Before an excited-state step, decide these and record "
             "them.** An engine default is not a decision, and a result "
             "whose level of theory nobody chose cannot be defended."]
    for i, rung in enumerate(RUNGS, 1):
        lines.append(f"{i}. **{rung.decide}** — {rung.because}. "
                     f"Closes it: {rung.closes}.")
    return "\n".join(lines)

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
        rung = _rung("ct-needs-long-range-exchange")
        out.append(ScienceFinding(
            label, "warning", rung.code,
            f"excited_state '{kwargs.get('excited_state')}' has "
            f"charge-transfer character and '{method}' is a global "
            f"hybrid: {rung.because}. Closes it: {rung.closes}"))

    # Rung 2: whether TDA was used is a decision and belongs in the spec.
    # ORCA applies it by default, so a step that is silent about it has
    # not chosen the full response -- it has inherited an approximation
    # that moves the numbers it reports.
    if "tda" not in kwargs:
        rung = _rung("tda-not-stated")
        out.append(ScienceFinding(
            label, "warning", rung.code,
            f"this step does not {rung.decide}: {rung.because}. "
            f"Closes it: {rung.closes} — a default is not a decision"))
    return out


__all__ = ["lint_excited"]
