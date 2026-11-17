"""CASSCF active-space proposal, validation and occupation checks.

Grounded in the verified knowledge cards delfin/orca_knowledge/cards/
casscf.json and nevpt2.json (ORCA 6.1 manual sections 3.13.x and 3.16),
not in a second rule list. The key facts used here:

- A CASSCF calculation always needs ``nel`` and ``norb`` in the %casscf
  block (card casscf.active-space.mandatory, manual 3.13.3 Basic
  Usage). There are no keywords ``Nactel`` or ``Nactorb``.
- Beyond roughly 14-16 active orbitals CASSCF becomes infeasible;
  RAS, ICE or DMRG are the alternatives (skill casscf-setup.md,
  warning rule).
- After convergence every active orbital's natural occupation should
  be clearly partial (about 0.02-1.98); an occupation pinned at 0 or 2
  means a wasted orbital and the space should be resized.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Sequence

__all__ = [
    "validate_casscf_block",
    "check_active_occupations",
    "OccupationReport",
    "ActiveSpaceProposal",
    "CasscfConfirmationGate",
    "MAX_ORBITALS",
    "parse_uno_occupations",
]

# CASSCF beyond ~16 orbitals is beyond what the module can set up.
MAX_ORBITALS = 16

# Allowed %casscf block keys. Each one is a real ORCA keyword from the
# verified cards: nel/norb/mult/nroots (mandatory set, manual 3.13.3),
# trafostep (casscf.trah.needs-aux-basis context; nevpt2 cards use
# TrafoStep RI), ptmethod (nevpt2.method-selection block form with the
# underscore), maxiter and weights are standard %casscf options.
_ALLOWED_KEYS = frozenset(
    {"nel", "norb", "mult", "nroots", "maxiter",
     "trafostep", "ptmethod", "weights"}
)
_MANDATORY_KEYS = ("nel", "norb")


def validate_casscf_block(block: Dict[str, str]) -> bool:
    """Validate a %casscf keyword block before it is written to an input.

    Unknown keys are refused with the offending name (a hallucinated
    keyword like ``Nactel`` must never reach an ORCA input file). The
    mandatory keys are ``nel`` and ``norb``.
    """
    for key in block:
        if key not in _ALLOWED_KEYS:
            raise ValueError(
                f"unknown %casscf keyword {key!r}: not in the ORCA 6 "
                f"manual (3.13.3 Basic Usage); real keys are "
                f"nel/norb/mult/nroots")
    missing = [k for k in _MANDATORY_KEYS if k not in block]
    if missing:
        raise ValueError(
            f"%casscf block is missing the mandatory key(s) {missing} "
            "(manual 3.13.3 Basic Usage)")
    nel, norb = int(block["nel"]), int(block["norb"])
    if nel <= 0:
        raise ValueError(
            f"nel={nel} is not a positive number of active electrons "
            "(manual 3.13.3: nel is the number of active electrons)")
    if norb <= 0:
        raise ValueError(
            f"norb={norb} is not a positive number of active orbitals "
            "(manual 3.13.3: norb is the number of active orbitals)")
    if norb > MAX_ORBITALS:
        raise ValueError(
            f"norb={norb} exceeds what CASSCF can do (max {MAX_ORBITALS} "
            "orbitals); use RAS, ICE or DMRG instead")
    if nel > 2 * norb:
        raise ValueError(
            f"nel={nel} exceeds 2*norb={2 * norb}: with at most 2 "
            "electrons per orbital, nel can be at most 2*norb")
    return True


# Pin thresholds for natural occupations (skill casscf-setup.md:
# clearly partial is about 0.02-1.98; a value at 0 or 2 is wasted).
_PIN_LOW, _PIN_HIGH = 0.02, 1.98


@dataclass
class OccupationReport:
    """Result of a natural-occupation check of the active orbitals."""

    occupations: List[float]
    pinned_indices: List[int] = field(default_factory=list)
    recommendation: str = ""

    @property
    def is_healthy(self) -> bool:
        # Empty occupations were never checked, so 'healthy' is not
        # claimed for them (an empty report must not read as healthy).
        return bool(self.occupations) and not self.pinned_indices


def parse_casscf_occupations(folder: str):
    """Read the converged CASSCF natural occupations from an ORCA 6
    output (the ``N(occ)=`` lines printed every CASSCF iteration,
    manual 3.13.x). Returns the LAST block (the converged iteration)
    in orbital order, or None when the folder has no .out or the
    output has no ``N(occ)=`` line (no data is not an empty list).
    """
    from pathlib import Path

    target = None
    for cand in sorted(Path(folder).glob("*.out"),
                       key=lambda p: p.stat().st_size, reverse=True):
        target = cand
        break
    if target is None:
        return None
    try:
        text = target.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    occ = None
    for line in text.split("\n"):
        if "N(occ)=" in line:
            rest = line.split("N(occ)=", 1)[1].split()
            occ = [float(x) for x in rest]
    return occ


def parse_uno_occupations(folder: str):
    """Read the UHF natural-orbital occupation numbers from an ORCA 6
    output produced with ``!UNO`` + ``Print[P_UNO_OccNum]`` (manual
    3.1.2: UNO saves the natural orbitals; 5.1.11: the occupation
    print). Returns the occupations in orbital-index order, or None
    when the folder has no .out or the output has no UNO block (no
    data is not an empty list).
    """
    from pathlib import Path

    target = None
    for cand in sorted(Path(folder).glob("*.out"),
                       key=lambda p: p.stat().st_size, reverse=True):
        target = cand
        break
    if target is None:
        return None
    try:
        text = target.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    m = re.search(
        r"UHF NATURAL ORBITALS\s*\n-+\s*\n\s*\n?Occupation Numbers:\n",
        text)
    if not m:
        return None
    occ: List[float] = []
    for line in text[m.end():].split("\n"):
        entries = re.findall(r"N\[\s*(\d+)\]=\s*(-?\d+\.\d+)", line)
        if not entries:
            if occ:
                break
            continue
        for idx_s, val_s in entries:
            idx = int(idx_s)
            if idx != len(occ):
                # Index gaps must not reorder the list silently.
                raise ValueError(
                    f"{target.name}: UNO occupation index {idx} out of "
                    f"order (expected {len(occ)})")
            occ.append(float(val_s))
    return occ or None


def check_active_occupations(
        occupations: Sequence[float]) -> OccupationReport:
    """Check CASSCF natural occupations of the active orbitals.

    An orbital pinned near 0 (empty) or 2 (doubly occupied) is wasted
    active space; the recommendation says which way to resize.

    Thresholds: an occupation is pinned when o < 0.02 or o > 1.98
    (STRICTLY outside); exactly 0.02 and 1.98 count as healthy
    (inclusive lower bound, inclusive upper bound). Empty input or
    None goes through check_active_occupations_safe, which defines
    that behavior; this function is undefined for an empty list.
    """
    occ = [float(x) for x in occupations]
    pinned = [i for i, o in enumerate(occ) if o < _PIN_LOW or o > _PIN_HIGH]
    if pinned:
        if any(o > _PIN_HIGH for i, o in enumerate(occ) if i in pinned):
            rec = (f"{len(pinned)} active orbital(s) pinned (indices "
                   f"{pinned}): shrink the space or rotate these orbitals "
                   "out of the active set")
        else:
            rec = (f"{len(pinned)} active orbital(s) empty (indices "
                   f"{pinned}): drop them from the active space")
    else:
        rec = ("all active occupations clearly partial; space is healthy")
    return OccupationReport(occupations=occ, pinned_indices=pinned,
                            recommendation=rec)


def check_active_occupations_safe(occupations) -> OccupationReport:
    """Defined behavior for degenerate inputs instead of an opaque crash.

    None is treated as 'not computed yet' (empty report, explicit
    note), an empty list as 'no active orbitals read' (also not
    healthy: nothing to check). Both carry a recommendation.
    """
    if occupations is None:
        return OccupationReport(
            occupations=[], pinned_indices=[], recommendation=(
                "no natural occupations available: the CASSCF run has "
                "not produced them yet (or the output was not parsed)"))
    occ = [float(x) for x in occupations]
    if not occ:
        return OccupationReport(
            occupations=[], pinned_indices=[], recommendation=(
                "no active orbitals were read: nothing was checked"))
    return check_active_occupations(occ)


# Chemistry vocabulary a justification must touch (skill casscf-setup.md
# section 1: propose the active space FROM the chemistry). The check is
# a floor, not a proof: it forces a chemical statement; coherence
# between the statement and the space is judged by the user/reviewer,
# not by keyword matching alone.
_JUSTIFICATION_WORDS = ("d", "orbital", "ligand", "pi", "bond", "metal",
                        "shell", "degener", "reference", "state",
                        "electron", "correlat")


@dataclass
class ActiveSpaceProposal:
    """A CAS(n,m) proposal with its written chemical justification.

    ``validate`` refuses an oversized space (beyond MAX_ORBITALS the
    alternatives RAS / ICE / DMRG are named in the error) and a
    justification without chemical content.

    The consistency rule (ONE rule, documented): the proposal must name
    the orbital class(es) that hold the active electrons
    (``named_orbital_classes``); the free-text justification must carry
    chemistry vocabulary. Whether the named classes are chemically
    COHERENT with the space is judged by the user, not by keyword
    matching.
    """

    n_electrons: int
    n_orbitals: int
    justification: str = ""
    orbital_rotations: str = ""
    named_orbital_classes: tuple = ()

    def validate(self) -> bool:
        if self.n_orbitals > MAX_ORBITALS:
            raise ValueError(
                f"CAS({self.n_electrons},{self.n_orbitals}) exceeds what "
                f"CASSCF can do (max {MAX_ORBITALS} orbitals); use RAS, "
                "ICE or DMRG instead")
        text = self.justification.strip().lower()
        if not text:
            raise ValueError(
                "a chemical justification is required: state WHY this "
                "active space (TM d shell, bonding ligand orbitals, "
                "full pi space, the bond being broken, ...)")
        if not any(re.search(rf"\b{re.escape(w)}\b", text)
                   for w in _JUSTIFICATION_WORDS):
            raise ValueError(
                f"the justification {self.justification!r} carries no "
                "chemistry: name the orbitals/electrons involved and why")
        if not self.named_orbital_classes:
            raise ValueError(
                "the proposal must name the orbital class(es) holding "
                "the active electrons (e.g. '5 metal d electrons + 1 "
                "pi bonding electron'); their coherence with the space "
                "is checked by the user")
        if self.n_electrons > 2 * self.n_orbitals:
            raise ValueError(
                f"nel={self.n_electrons} exceeds 2*norb="
                f"{2 * self.n_orbitals}")
        return True


class CasscfConfirmationGate:
    """The CASSCF step runs only after an explicit structured user
    confirmation.

    Confirmation is a typed flag (``confirm(user_confirmed=True)``),
    never free text: a text answer -- 'ok', 'yes', a unicode look-alike
    -- carries no authority and raises TypeError. ``start_casscf``
    raises PermissionError while the confirmation is missing.

    A gate is created fresh per workflow run: confirmation is never
    inherited across re-entries.
    """

    def __init__(self) -> None:
        self._confirmed = False

    def confirm(self, *, user_confirmed: bool) -> None:
        if not isinstance(user_confirmed, bool):
            raise TypeError(
                "confirmation is the structured flag user_confirmed=True "
                "from the user; text answers are not a confirmation")
        self._confirmed = bool(user_confirmed)

    def revoke(self) -> None:
        """The user withdraws the confirmation (e.g. after revising the
        proposal); the gate is closed again until re-confirmed."""
        self._confirmed = False

    def start_casscf(self) -> bool:
        if not self._confirmed:
            raise PermissionError(
                "the CASSCF step needs the user's explicit confirmation "
                "of the proposal before it starts")
        return True
