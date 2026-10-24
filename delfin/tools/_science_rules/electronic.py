"""P1 electronic realism rules for the scientific lint.

Called from delfin.tools._science._lint_step — this is not a second linter,
it is where these rules live so the single entry point stays
``scientific_lint``.

Sources:
- Diffuse functions for anions: J. Phys. Chem. A 2013, 117, 5017–5027,
  DOI 10.1021/jp312755z; minimally augmented basis sets
  (ma-/maug-) J. Chem. Theory Comput. 2009, 5, 2948, DOI 10.1021/ct900566x;
  diffuse functions for Rydberg/dipole-bound states
  Theor. Chem. Acc. 2015, 134:30, DOI 10.1007/s00214-015-1647-1.
- Relativistic elements: DELFIN's own policy delfin/utils.py:203
  ``_should_use_rel`` (4d/5d/4f/5f/mixed) — reused by import, never copied.
- Parity: DELFIN's own rule delfin/config.py:709-764, quoting ORCA's stop
  "multiplicity (4) is even and number of electrons (242) is even ->
  impossible"; parity-derived multiplicity delfin/utils.py:414-428.
"""

from __future__ import annotations

from typing import Any, Dict, List, TYPE_CHECKING

if TYPE_CHECKING:  # annotation only — real import is lazy (parent imports us)
    from delfin.tools._science import ScienceFinding

# Diffuse basis-set markers. ma- and maug- are DELFIN's own vocabulary
# (delfin/common/control_validator.py:573-580, 734-741); aug- likewise.
_DIFFUSE_MARKERS = ("ma-def2", "maug-", "aug-", "pcH-", "def2-ma")

# A step-level hint that the electronic structure involves an anion-like,
# Rydberg or diffuse excited state.
_EXCITED_STATE_MARKERS = ("rydberg", "diffuse_excited", "dipole_bound")


def _needs_diffuse_basis(kwargs: Dict[str, Any]) -> bool:
    """An anion, or a step explicitly marked as a Rydberg/diffuse excited state."""
    from delfin.tools._science import _is_placeholder
    if _is_placeholder(kwargs.get("charge")) or _is_placeholder(kwargs.get("basis")):
        return False
    if isinstance(kwargs.get("charge"), int) and kwargs["charge"] < 0:
        return True
    es = str(kwargs.get("excited_state") or "").strip().lower()
    return any(m in es for m in _EXCITED_STATE_MARKERS)


def lint_electronic(label: str, kwargs: Dict[str, Any]) -> List[ScienceFinding]:
    """Rules 3, 4, 10 of the WB1 physical-realism package."""
    from delfin.tools._science import ScienceFinding, _is_placeholder
    out: List[ScienceFinding] = []

    # Rule 3: anions / Rydberg / diffuse excited states need diffuse functions.
    basis = kwargs.get("basis")
    if _needs_diffuse_basis(kwargs) and isinstance(basis, str) and basis:
        if not any(m in basis for m in _DIFFUSE_MARKERS):
            out.append(ScienceFinding(
                label, "warning", "diffuse-basis-required",
                f"charge {kwargs.get('charge')} / state '{kwargs.get('excited_state')}' "
                f"needs diffuse basis functions; '{basis}' has none — use a "
                "ma-def2-* / aug-cc-* / maug-cc-* basis"))

    # Rule 4: 4d/5d/4f/5f elements need a relativistic Hamiltonian + basis.
    # Decision reused from DELFIN's own policy, never a copy.
    metals = kwargs.get("metals")
    if isinstance(metals, list) and metals and not _is_placeholder(metals):
        from delfin.utils import _should_use_rel
        relativity = str(kwargs.get("relativity") or "").strip().lower()
        basis_s = str(basis or "")
        has_rel = any(m in basis_s for m in ("zora", "dkh", "x2c", "sarc"))
        if _should_use_rel([str(m) for m in metals]) and not (relativity or has_rel):
            out.append(ScienceFinding(
                label, "warning", "relativistic-required",
                f"metals {metals} are 4d/5d/4f/5f (or mixed): relativistic "
                "treatment required — set relativity (ZORA/X2C/DKH) and a "
                "recontracted basis"))

    # Rule 10: electron count compatible with multiplicity at every step.
    mult = kwargs.get("mult")
    if isinstance(mult, bool):
        mult = None
    nelec = kwargs.get("electrons")
    if isinstance(mult, int) and mult >= 1 and isinstance(nelec, int):
        # ORCA stops the run when parities disagree (delfin/config.py:709-764).
        if (mult % 2 == 0) == (nelec % 2 == 0):
            out.append(ScienceFinding(
                label, "error", "electron-multiplicity-parity",
                f"multiplicity {mult} is {'even' if mult % 2 == 0 else 'odd'} "
                f"and {nelec} electrons are "
                f"{'even' if nelec % 2 == 0 else 'odd'}: ORCA stops this run "
                "('multiplicity ... -> impossible'); Nα+Nβ = "
                f"{nelec} must equal electrons, Nα−Nβ = mult − 1"))
    return out


__all__ = ["lint_electronic"]
