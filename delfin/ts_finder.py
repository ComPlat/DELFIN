"""Pure checks an expert uses to find and verify a transition state.

This module is the library behind the ``ts-finder`` and ``irc-checker``
sub-agent presets. It turns the method of ``orca-ts-workflow.md`` into
concrete, individually-checkable predicates:

* one level of theory across reactant / TS / product / ZPE,
* the single imaginary mode actually connects the intended bonds,
* the IRC ends are real minima that map to the stated endpoints.

None of these checks submit a job; they run on structures, Hessians and
OUTPUT files DELFIN has already produced (or that a test builds
synthetically). Every check returns ``CheckResult`` with an ``ok`` flag
and a list of human-readable problems, so the expert can act on them and
a test can assert on them precisely.

The verdict a whole TS search collapses to is ``summarize_ts_verdict``;
it is deliberately honest -- there is no "probably a TS".

Reuse, never re-implement:
* atom-order mapping uses ``delfin.atom_mapping.map_atoms``,
* mode-response uses the displaced-geometry pairs ``delfin.ts_goat``
  already scores (``_pair_response_fraction`` and its response window).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from delfin.atom_mapping import map_atoms
from delfin.ts_goat import (
    _MAX_RELEVANT_CONTACT_ANGSTROM,
    _RESPONSE_FRACTION,
    _DISTANCE_RESPONSE_THRESHOLD_ANGSTROM,
    _distance,
    _pair_response_fraction,
)

# Verdict labels -- honest and stable. A test pins them.
VERDICT_NOT_FOUND = "not found"
VERDICT_WRONG_MODE = "found but wrong mode"
VERDICT_VERIFIED = "found and verified"


@dataclass
class CheckResult:
    """Outcome of one pure check: passes with zero problems, else fails."""
    ok: bool
    problems: list[str] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.ok


# ---------------------------------------------------------------------------
# Level of theory
# ---------------------------------------------------------------------------

def check_single_level_of_theory(levels: dict) -> CheckResult:
    """A barrier is one level: reactant, TS, product AND ZPE must match.

    ``levels`` maps a stage name ("reactant", "ts", "product", "zpe") to
    the method line DELFIN recorded for it (functional/dispersion/basis/
    solvent, the ML_EFF charge·mult·method string the pipeline stores).
    A stage missing from the dict is a problem on its own -- the barrier
    must not be silent about an endpoint it did not run.

    Returns CheckResult; problems name every stage that differs from the
    first stage that carried a method.
    """
    named = sorted((s for s, m in (levels or {}).items() if m), key=list(levels))
    if not named:
        return CheckResult(False, ["no level of theory provided for any stage"])
    reference = named[0]
    problems: list[str] = []
    for stage, method in (levels or {}).items():
        if not method:
            problems.append(f"{stage}: no level of theory recorded")
        elif method != levels[reference]:
            problems.append(
                f"{stage}: level '{method}' differs from {reference} "
                f"'{levels[reference]}' -- a barrier assembled from "
                "different levels is not a barrier")
    # every stage the caller skipped (it is not in `levels` at all)
    expected = {"reactant", "ts", "product", "zpe"}
    for stage in sorted(expected - set(levels or {})):
        problems.append(f"{stage}: not provided -- cannot confirm one level")
    return CheckResult(not problems, problems)


# ---------------------------------------------------------------------------
# Atom-order endpoint mapping (reactant/product <-> IRC end / candidate)
# ---------------------------------------------------------------------------

def _read_xyz_lines(text: str):
    """Return (syms, xyz) from an xyz block or file text; None when no atoms.

    Accepts the multi-line xyz text a caller has already read from a file,
    so a check never reads paths the builder did not pass -- and tests can
    hand it synthetic blocks without touching the repo.
    """
    lines = [ln.strip() for ln in (text or "").splitlines() if ln.strip()]
    if len(lines) < 3:
        return None
    try:
        n_atoms = int(lines[0])
    except ValueError:
        return None
    syms: list[str] = []
    coords: list[list[float]] = []
    for raw in lines[2:2 + n_atoms]:
        parts = raw.split()
        if len(parts) < 4:
            continue
        try:
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            syms.append(parts[0])
        except ValueError:
            continue
    if len(syms) != n_atoms or len(syms) == 0:
        return None
    return syms, coords


def check_endpoint_mapping(reference_block, candidate_block) -> CheckResult:
    """Do ``reference_block`` and ``candidate_block`` map onto one another?

    Both are xyz text (a 2-line header + atom lines). The check verifies
    the two structures have the same molecular formula and that an atom
    mapping can be built between them (atom_mapping.map_atoms). It returns
    the RMSD in the problems text so a near-identical endpoint is told
    apart from a merely isomeric one. Returns ok=False on any
    unparseable block, a formula mismatch or a failed mapping.
    """
    ref = _read_xyz_lines(reference_block)
    cand = _read_xyz_lines(candidate_block)
    if ref is None:
        return CheckResult(False, ["reference block is not parseable xyz"])
    if cand is None:
        return CheckResult(False, ["candidate block is not parseable xyz"])
    try:
        result = map_atoms(ref[0], ref[1], cand[0], cand[1])
    except (ValueError, RuntimeError) as exc:
        return CheckResult(False, [f"atom mapping failed: {exc}"])
    if not result.get("verified"):
        return CheckResult(False, [
            f"structures do not map onto one another "
            f"(method={result.get('method')})"])
    rmsd = result.get("rmsd")
    problems: list[str] = []
    if rmsd is not None and rmsd > 0.5:
        problems.append(
            f"mapped with RMSD {rmsd:.3f} Angstrom -- the two structures "
            "share the atom order but are not the same geometry; is this "
            "the intended endpoint?")
    return CheckResult(not problems, problems)


# ---------------------------------------------------------------------------
# IRC endpoint match
# ---------------------------------------------------------------------------

def check_irc_endpoint_match(irc_end_block, reactant_block) -> CheckResult:
    """Is the given IRC end the stated reactant?

    An IRC end is only meaningful as a mechanism endpoint when it maps onto
    the molecule the caller claims it is (the reactant, or the product).
    This is ``check_endpoint_mapping`` with an IRC-flavoured message, kept
    as its own public function so the ``irc-checker`` preset and its tests
    have a stable name.
    """
    result = check_endpoint_mapping(reactant_block, irc_end_block)
    if result.ok:
        return result
    return CheckResult(False, [f"IRC end does not match the stated reactant",
                               *result.problems])


def check_irc_ends_are_minima(imag_counts_per_end: dict) -> CheckResult:
    """Each IRC end must be a real minimum: zero imaginary modes.

    ``imag_counts_per_end`` maps an end label ("forward", "reverse") to
    the number of imaginary modes its frequency run reported. A minimum
    with one or more imaginary modes is a saddle, not an endpoint, and the
    mechanism must not be claimed.
    """
    problems: list[str] = []
    for end, count in (imag_counts_per_end or {}).items():
        if count is None:
            problems.append(f"{end} IRC end: no frequency run to judge")
        elif count != 0:
            problems.append(
                f"{end} IRC end: {count} imaginary mode(s) -- not a minimum")
    if not imag_counts_per_end:
        problems.append("no IRC end frequency runs provided")
    return CheckResult(not problems, problems)


# ---------------------------------------------------------------------------
# Imaginary mode connects the intended bonds
# ---------------------------------------------------------------------------

def check_imaginary_mode_connects_bonds(
    hess_path,
    ts_xyz,
    expected_pairs,
    mode_index=None,
):
    """Does the (selected) imaginary mode act on the bond-forming/breaking
    pairs ``expected_pairs``?

    ``hess_path`` is a path to the frequency Hessian DELFIN wrote (the
    ``.hess`` from an ``orca_ts_freq_check``/``orca_ts_opt`` run).
    ``ts_xyz`` is xyz text of the TS geometry (same atom order as the
    Hessian). ``expected_pairs`` is the list of ``(i, j)`` 0-based atom
    indices whose bond this TS is supposed to form or break (the reactive
    centres, e.g. the H-transfer donor and acceptor).

    The check reuses the pair-response rule ts_goat already applies: a
    contact counts as "acted on by the mode" when displacing the geometry
    along the mode by 0.25 Angstrom changes its length by more than
    ``_RESPONSE_FRACTION`` of its equilibrium length and more than
    ``_DISTANCE_RESPONSE_THRESHOLD_ANGSTROM`` -- the same threshold that
    drives ts_goat's reactive-centre score. It reuses ``displaced_geometries``
    and ``_pair_response_fraction`` rather than re-deriving them.

    Returns CheckResult; when a Hessian cannot be read or the mode index
    is unknown, it fails with the reason rather than guessing.
    """
    try:
        from delfin.imag import displaced_geometries, imaginary_modes
    except Exception as exc:  # pragma: no cover -- import wiring
        return CheckResult(False, [f"imag-mode analysis unavailable: {exc}"])

    atoms = _read_xyz_lines(ts_xyz)
    if atoms is None:
        return CheckResult(False, ["TS geometry is not parseable xyz"])
    syms, coords = atoms
    n = len(syms)
    # ts_goat._distance expects 4-tuples (element, x, y, z) -- the shape
    # displaced_geometries returns. Reuse that shape so the same helper
    # measures the equilibrium and the displaced geometry identically.
    coords4 = [(syms[i], coords[i][0], coords[i][1], coords[i][2])
               for i in range(n)]
    if any(
        i >= n or j >= n or i < 0 or j < 0 or i == j
        for i, j in expected_pairs
    ):
        return CheckResult(False, [
            "expected bond pair lies outside the atom count "
            f"(n={n}, pairs={expected_pairs!r})"])

    hess = Path(hess_path)
    if not hess.is_file():
        return CheckResult(False, [f"no Hessian at {hess_path}"])

    if mode_index is None:
        try:
            # An empty config keeps the noise floor default: a clearly
            # imaginary mode below it is the reactive mode candidate.
            modes = imaginary_modes(hess, {})
        except Exception as exc:
            return CheckResult(
                False, [f"cannot read imaginary modes from Hessian: {exc}"])
        if not modes:
            return CheckResult(
                False, ["no imaginary mode found -- this is a minimum, "
                        "not a transition state"])
        # most negative first is how ts_goat pre-selects the reactive mode
        mode_index = modes[0][0]

    displaced = None
    try:
        displaced = displaced_geometries(hess, mode_index, 0.25)
    except (ValueError, OSError) as exc:
        return CheckResult(False, [
            f"displacement along mode {mode_index} unavailable: {exc}"])

    problems: list[str] = []
    for i, j in expected_pairs:
        d0 = _distance(coords4[i], coords4[j])
        if d0 > _MAX_RELEVANT_CONTACT_ANGSTROM:
            problems.append(
                f"pair ({i},{j}): {d0:.2f} A is beyond the bonding "
                f"window ({_MAX_RELEVANT_CONTACT_ANGSTROM:.1f} A)")
            continue
        fraction = _pair_response_fraction(i, j, d0, displaced)
        if d0 <= 0:
            problems.append(f"pair ({i},{j}): zero equilibrium distance")
            continue
        change = fraction * d0
        if change <= _DISTANCE_RESPONSE_THRESHOLD_ANGSTROM:
            problems.append(
                f"pair ({i},{j}): mode {mode_index} changes the bond by "
                f"{change * 1000:.1f} mA (< {_DISTANCE_RESPONSE_THRESHOLD_ANGSTROM * 1000:.0f} mA) "
                "-- the imaginary mode does not act on this bond")
    if not problems:
        problems.append(
            f"mode {mode_index} acts on all {len(expected_pairs)} "
            "expected bond-forming/breaking pair(s)")
    return CheckResult(not problems, problems)


def summarize_ts_verdict(found: bool, mode_ok: bool, irc_ok: bool) -> str:
    """The honest headline of a TS search, one of three labels.

    * not found          -- no candidate reached an OptTS (or none had a
                            single imaginary mode that was verified),
    * found but wrong mode -- a candidate exists but the imaginary mode
                            does not connect the intended bonds (and/or
                            the IRC ends do not match),
    * found and verified -- all checks passed: one level, one imaginary
                            mode acting on the right bonds, IRC to both
                            endpoints that are real minima.

    A candidate is only ``found and verified`` when every check is green;
    anything less is reported honestly, never rounded up.
    """
    if not found:
        return VERDICT_NOT_FOUND
    if not (mode_ok and irc_ok):
        return VERDICT_WRONG_MODE
    return VERDICT_VERIFIED
