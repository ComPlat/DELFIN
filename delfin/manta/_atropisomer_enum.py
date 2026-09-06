"""_atropisomer_enum.py — ADDITIVE completeness for AXIAL chirality (atropisomers).

THE GAP, MEASURED ON 16.08.2026 (not assumed).  The axis `ccdc_atropisomer_realized` stood
at **0 of 44** on `rows_spy5geom` (965 systems) -- not a single system realizes the axial
handedness of its crystal.  With random handedness ~22 hits would be expected; zero is
systematic.

THE CAUSE IS NOT THE GEOMETRY.  Measured over all 44 systems and **145 axes** (twist over
ALL frames, not only frame 0):

    mean twist   42.8 degrees      89 of 145 axes >= 30 degrees, only 22 < 10 degrees
    sign in the manifold:  both 48 (33%) | ONLY ONE 59 (41%) | never reliable 38 (26%)

The builder produces real, twisted axes -- `BIRVUW` reaches 85.2 degrees.  What is missing is
the OPPOSITE HANDEDNESS.  And the axis demands that EVERY crystal axis find its sign in the
manifold: with ~3.3 axes per system and only a third of the axes carrying both hands, a full
hit is arithmetically almost impossible.  **That explains 0/44 quantitatively, and it makes
the enumeration the answer -- not a geometry correction.**

MESOMERISM (user, 16.08.).  `BIRVUW` carries the axes N29-C30 and N18-C19 -- **aryl-amide**
axes.  Their stereogenicity arises only through the partial C-N double-bond character in the
first place: the free N lone pair conjugates into the carbonyl group, the bond acquires
double-bond character, and from that follow the barrier and the preferred twist.  Whoever
treats it as a single bond sees no axis there at all.  The case is still covered -- but
since 23.08. through the pi branch of `_atrop_is_sp2_center` (an amide N with C-N <= 1.37 A
falls into it by itself), no longer through a special case of its own.  The difference is
not cosmetic: the special case checked the same chemistry with a DIFFERENT threshold than
the eye, and so the builder saw axes there that the eye did not carry, and vice versa.

🔴 WHY THE MODULE HAD NO EFFECT UNTIL 23.08.2026 -- MEASURED, NOT ASSUMED.
`harness/atrop_schluessel_vergleich.py` holds both axis detections against each other on THE
SAME frames.  On 120 systems / 3350 frames:

                              before      after
    axes builder                 112          282
    signature different      103/103        0/282
    sign different                87            0
    axes only in builder           9            0

**On EVERY shared bond the signature was different, on 84 % the sign as well.**
Thus every appended mirror image was a frame in the WRONG bucket -- and a frame in the
wrong bucket costs exactly as much as a missing one.  Two independent causes:

  * the element table mapped everything outside its 16 entries to **0** (silver: 0 here,
    29 in the eye) -> different flank rank -> different reference atom -> different sign;
  * `_atrop_dihedral` formed the first vector reversed -> the dihedral was shifted by
    **180 degrees**.  The folded magnitude stays the same, which is why it never showed up.

Since 23.08. the axis definition is, stage for stage, that of the eye (adjacency, rings,
sp2 test, single-bond band, flank rank, dihedral, twist band 20-88 instead of 10-80).
`atrop_schluessel_vergleich.py` is the guard against them drifting apart again.

THE OPERATION IS EXACT, NOT AN OPTIMIZATION.  The mirror atropisomer has the same MAGNITUDE
of twist and the opposite SIGN.  A rotation of one side about the axis by `-2*theta` maps
`theta -> -theta`: sign flipped, magnitude preserved.  A rigid rotation about the BOND AXIS
changes exclusively the torsion -- all bond lengths and all bond angles stay exactly the
same.  Therefore this step needs NO relaxation and cannot BREAK a bond by construction.

⚠ BUT IT CAN CREATE A BOND, and that went unnoticed until 23.08.  The clash floor was at
1.70 A, the builder's organic bond threshold at r_i+r_j+0.25 -- for C-C **1.77 A**.  A
rotated atom was thus allowed to approach to 1.71 A, passed the clash check and got a
spurious bond; with that the flank profile changed, and so did the axis key.  Instead of
setting one number against another, `_atrop_topologie_gleich` has since checked the PROPERTY
itself: does the bond topology stay unchanged -- in BOTH adjacencies, the builder's and the
eye's?  Otherwise no frame.

ADDITIVE AND NEVER-WORSE.  Originals remain unchanged; only MISSING signs are appended.
Default OFF -> byte-identical.  Capped, and LOGGED at the cap (no silent truncation).

⚠ WHAT THIS MODULE DOES NOT DO: it does not judge WHICH sign is correct.  That would be a
crystal question, and the construction must do without the crystal.  It ensures that BOTH
are present in the manifold -- completeness rule, not a selection.

⚠ DELIMITATION (checked 16.08., before building): `_chirality_enumerator.py` handles the
Lambda/Delta HELICITY of the coordination polyhedron from chelate pairs -- metal-centered
chirality.  Its "axial" sites mean axial POSITIONS (axial versus equatorial), not axial
chirality.  There is no second atropisomer enumeration in the builder.

⚠ NAMES: the private helpers all carry the prefix `_atrop_`, and the entry point is called
`expand_atropisomers` instead of `expand_results`.  The watchdog `check_exists_first` had
reported the generic names (`_env_int`, `_is_enabled`, `_dihedral`, `expand_results`) as a
collision; instead of overriding it they were made unambiguous.
"""
from __future__ import annotations

import logging
import math
import os
from typing import Dict, List, Optional, Sequence, Set

import numpy as np

# The same primitives as `_stereocenter_enum` -- one source, no drift.
from delfin.manta._coord_angle_corrector import (
    _build_geometric_adjacency,
    _format_xyz,
    _is_metal_sym,
    _parse_xyz,
)

_LOG = logging.getLogger(__name__)

# ===== THE AXIS DEFINITION IS THAT OF THE EYE -- CHARACTER FOR CHARACTER (23.08.2026) =======
#
# MEASURED, not assumed (`harness/atrop_schluessel_vergleich.py` on 120 systems,
# 3350 frames): on the 103 bonds that BOTH sides held to be an axis, the signature was
# different in **103 of 103 cases** and the SIGN in **87 of 103**.
# The builder found 112 axes, the eye 2090.  Thus every mirror image from this module was a
# frame in the WRONG bucket -- and a frame in the wrong bucket costs exactly as much as a
# missing one.  That explains `atrop44` (16.08., 9x false in BOTH arms) and `atrop10k`
# (21.08., 35 -> 35 at +316 frames) without remainder.
#
# THE ONE EXAMPLE THAT SHOWS IT (ABOZIB C51-C52, the same bond, the same frame):
#     builder: ('C', ((6, 6, 6), 0), ...)  and  ('C', ((8, 0), 0), ...)
#     eye    : ('C', ((29, 8), 0), ...)    and  ('C', ((6, 6, 6), 0), ...)
# The `0` is SILVER.  The old table knew 16 elements and mapped everything else to 0; the
# eye falls back to `round(covalent radius * 20)`.  With that the flank rank flips, with it
# the reference atom of the dihedral angle, with it the sign.  "Deliberately kept SMALL,
# it is about an ORDERING" was exactly the fallacy: an ordering is only the same if both
# sides use the same ordering.
#
# NO IMPORT FROM THE EYE -- the construction must not depend on the measuring instrument
# ("eye = UPPER BOUND of the construction").  That is why the definition stands here a
# SECOND TIME, deliberately, and `harness/atrop_schluessel_vergleich.py` is the guard
# against the two drifting apart again: it must report 0 % SIG_VERSCHIEDEN.
#
# Origin of every single number: `weddell/detectors/atropisomer_sign.py`, lines 101-131.

# Twist band: below it planar (no sign), above it perpendicular (the two hands become
# indistinguishable).  FORMERLY 10/80 -- too wide at the bottom (the eye credits nothing
# below 20 degrees, every such frame was scrap), too narrow at the top (`BIRVUW` carries
# axes at 83.7 and 85.2 degrees, which the eye DEMANDS and the builder could not see).
_ATROP_TWIST_MIN_DEG = 20.0
_ATROP_TWIST_MAX_DEG = 88.0

# Atomic numbers for the CIP-like flank rank -- identical to `_Z` of the eye.
_ATROP_Z = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Al": 13, "Si": 14, "P": 15, "S": 16,
    "Cl": 17, "As": 33, "Se": 34, "Br": 35, "Te": 52, "I": 53,
}

# Covalent radii -- identical to `_COV_R` of the eye.  They carry THREE roles: the fallback
# for the atomic number, the adjacency and the axis length.  A deviation here shifts all
# three at once.
_ATROP_COV_R = {
    "H": 0.31, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "K": 2.03,
    "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "La": 2.07, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
}

_ATROP_SP2_ELEMS = {"C", "N", "O", "S", "Se", "B", "P"}   # sp2-capable ring/center elements
_ATROP_ARO_BOND_MAX = 1.46        # mean ring bond below this = aromatic/conjugated
_ATROP_RING_OOP_MAX = 0.35        # RMS deviation from the ring plane (A) below this = planar
_ATROP_SP2_SUM_MIN = 348.0        # sum of the three neighbor angles above this = planar sp2
# Axis single-bond band: above the double-bond shoulder (a C=C/C=N belongs to `ez_stereo`,
# not here) and below margin * sum of the covalent radii.
_ATROP_DOUBLE_CUT = {("C", "C"): 1.42, ("C", "N"): 1.37, ("N", "N"): 1.37, ("C", "O"): 1.36}
_ATROP_LEN_MARGIN = 1.18
_ATROP_ADJ_CUT_FRAC = 1.30        # adjacency of the eye: d <= 1.30 * (r_i + r_j)


def _atrop_znum(sym: str) -> int:
    """Identical to `_znum` of the eye -- table, otherwise covalent radius * 20."""
    return _ATROP_Z.get(sym, int(round(_ATROP_COV_R.get(sym, 0.9) * 20)))


def _atrop_env_int(name: str, default: int) -> int:
    try:
        return int(str(os.environ.get(name, default)).strip())
    except Exception:
        return default


def _atrop_env_float(name: str, default: float) -> float:
    try:
        return float(str(os.environ.get(name, default)).strip())
    except Exception:
        return default


def _atrop_enabled() -> bool:
    """THE ONE read site.  Default OFF -> byte-identical."""
    return (os.environ.get("DELFIN_FFFREE_ATROPISOMER_ENUM", "0") == "1"
            or os.environ.get("DELFIN_ATROPISOMER_ENUM", "0") == "1")


# ---------------------------------------------------------------------------------------------
# Axis detection -- geometric, without an element table and without the eye
# ---------------------------------------------------------------------------------------------

def _atrop_heavy_nbrs(syms: Sequence[str], nbrs: List[List[int]], i: int) -> List[int]:
    return [j for j in nbrs[i] if syms[j] != "H"]


def _atrop_adjacency(syms: Sequence[str], pts: np.ndarray,
                     cut_frac: float = _ATROP_ADJ_CUT_FRAC) -> List[List[int]]:
    """Adjacency OF THE EYE: d <= cut_frac * (r_i + r_j), heavy AND H.

    ⚠ NOT `_build_geometric_adjacency` (r_i+r_j+0.25 organic).  For C-C that is
    1.98 versus 1.77 A -- different neighborhoods, different flank profiles,
    different signatures.  The axis detection MUST see the neighbors that the eye
    sees; the clash check of the mirroring still uses the tighter build adjacency.
    """
    n = len(syms)
    r = np.array([_ATROP_COV_R.get(s, 1.5) for s in syms])
    P = np.asarray(pts, float)
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    cut = cut_frac * (r[:, None] + r[None, :])
    np.fill_diagonal(D, 1e9)
    nbr: List[List[int]] = [[] for _ in range(n)]
    ii, jj = np.where(D <= cut)
    for a, b in zip(ii.tolist(), jj.tolist()):
        if a < b:
            nbr[a].append(b)
            nbr[b].append(a)
    return nbr


def _atrop_aromatic_rings(nbrs: List[List[int]], syms: Sequence[str],
                          pts: np.ndarray) -> List[frozenset]:
    """Aromatic-like sp2 rings (size 5-6, sp2-capable elements, short ring bonds,
    planar) -- identical to `_aromatic_rings` of the eye.

    ⚠ REPLACES `_atrop_ring_of` (smallest ring from a breadth-first search, up to size 8).
    The ring determines which atoms are EXCLUDED in the flank profile; a different ring
    means a different profile.  A cyclohexane ring was a ring for the old version and
    not one for the eye.
    """
    P = np.asarray(pts, float)
    adjA = {i: [k for k in nbrs[i] if syms[k] in _ATROP_SP2_ELEMS]
            for i in range(len(syms)) if syms[i] in _ATROP_SP2_ELEMS}
    found = set()
    for start in adjA:
        stack = [(start, (start,))]
        while stack:
            node, path = stack.pop()
            for nb in adjA.get(node, ()):
                if nb == start and 5 <= len(path) <= 6:
                    found.add(frozenset(path))
                elif nb not in path and len(path) < 6:
                    stack.append((nb, path + (nb,)))
    rings = []
    for r in found:
        idx = list(r)
        bl = [float(np.linalg.norm(P[a] - P[b])) for ai, a in enumerate(idx)
              for b in idx[ai + 1:] if b in nbrs[a]]
        if not bl or (sum(bl) / len(bl)) > _ATROP_ARO_BOND_MAX:
            continue
        Q = P[idx]
        c = Q.mean(0)
        try:
            _u, _s, vt = np.linalg.svd(Q - c)
        except np.linalg.LinAlgError:
            continue
        oop = float(np.sqrt(np.mean(np.dot(Q - c, vt[2]) ** 2)))
        if oop > _ATROP_RING_OOP_MAX:
            continue
        rings.append(r)
    return rings


def _atrop_sum_angles(pts: np.ndarray, c: int, neigh: List[int]) -> float:
    if len(neigh) < 3:
        return 0.0
    s = 0.0
    for a in range(len(neigh)):
        for b in range(a + 1, len(neigh)):
            v1 = pts[neigh[a]] - pts[c]
            v2 = pts[neigh[b]] - pts[c]
            nn = float(np.linalg.norm(v1)) * float(np.linalg.norm(v2))
            if nn < 1e-9:
                continue
            s += math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2)) / nn))))
    return s


def _atrop_is_sp2_center(i: int, nbrs: List[List[int]], syms: Sequence[str],
                         pts: np.ndarray, ring_of: List[List[frozenset]]) -> bool:
    """sp2 unit -- identical to `_is_sp2_center` of the eye: ring member OR
    three-coordinate and planar OR with one neighbor at double-bond distance.

    ⚠ THE LAST BRANCH REPLACES `_atrop_is_mesomeric_amide`.  An amide N with C-N <= 1.37 A
    falls into it by itself -- without a special case, and above all: EXACTLY WHEN the
    eye does so too.  The special case checked the same chemistry with a different threshold.
    """
    if syms[i] == "H" or _is_metal_sym(syms[i]) or syms[i] not in _ATROP_SP2_ELEMS:
        return False
    if ring_of[i]:
        return True
    heavy = [k for k in nbrs[i] if syms[k] != "H" and not _is_metal_sym(syms[k])]
    if len(heavy) == 3 and _atrop_sum_angles(pts, i, heavy) >= _ATROP_SP2_SUM_MIN:
        return True
    for k in heavy:
        cut = (_ATROP_DOUBLE_CUT.get((syms[i], syms[k]))
               or _ATROP_DOUBLE_CUT.get((syms[k], syms[i])))
        d = float(np.linalg.norm(pts[i] - pts[k]))
        rr = _ATROP_COV_R.get(syms[i], 1.5) + _ATROP_COV_R.get(syms[k], 1.5)
        if (cut is not None and d <= cut) or d <= 0.93 * rr:
            return True
    return False


def _atrop_axis_is_single(i: int, j: int, syms: Sequence[str], pts: np.ndarray) -> bool:
    """Real single bond between two units -- identical to `_axis_is_single`."""
    d = float(np.linalg.norm(pts[i] - pts[j]))
    lo = (_ATROP_DOUBLE_CUT.get((syms[i], syms[j]))
          or _ATROP_DOUBLE_CUT.get((syms[j], syms[i])))
    rr = _ATROP_COV_R.get(syms[i], 1.5) + _ATROP_COV_R.get(syms[j], 1.5)
    if lo is not None and d <= lo:
        return False                 # clear double bond (E/Z range) -- no atrop axis
    return d <= _ATROP_LEN_MARGIN * rr


def _atrop_share_ring(i: int, j: int, ring_of: List[List[frozenset]]) -> bool:
    return any(r in ring_of[j] for r in ring_of[i])


def _atrop_side_atoms(syms: Sequence[str], nbrs: List[List[int]],
                      i: int, j: int) -> Optional[Set[int]]:
    """All atoms on the j side of the bond i-j (without i), metals as the boundary.

    None if the bond lies IN AN ORGANIC RING -- then there are no two sides, and a
    rotation would tear the structure apart.

    ⚠ METALS ARE NOT TRAVERSED (16.08.2026, from the self-test).  On `BEBGUL` this module
    found NO axis, while the eye carries `C26-C27` there with 81.2 degrees.  Reason: in a
    CHELATING biaryl BOTH aryl rings coordinate -- so there is a path from one side to the
    other, via the METAL.  The search thus took the axis for a ring bond.  That is a
    METALLACYCLE, not an organic ring; the eye does not see it, because it uses RDKit's
    ring detection on the organic graph.

    ⚠⚠ FROM THIS FOLLOWS A BUILD LIMIT THAT CANNOT BE PATCHED AWAY.  If BOTH sides hang on
    the metal, the opposite handedness can NOT be produced by rigid rotation -- it would
    tear the M-D bond apart.  Exactly the property this module relies on ("a rotation
    about the bond axis changes exclusively the torsion") no longer holds there.
    `_atrop_find_axes` rejects such axes EXPLICITLY instead of building them wrongly:
    their opposite handedness needs a RESEATING, not a rotation.
    """
    seen = {j}
    stack = [k for k in nbrs[j] if k != i]
    while stack:
        cur = stack.pop()
        if cur == i:
            return None                       # organic ring closure -> no rotation axis
        if cur in seen:
            continue
        seen.add(cur)
        if _is_metal_sym(syms[cur]):
            continue                          # metal: boundary of the organic side
        stack.extend(k for k in nbrs[cur] if k not in seen)
    return seen


def _atrop_dihedral(pts: np.ndarray, a: int, i: int, j: int, b: int) -> Optional[float]:
    """Signed dihedral a-i-j-b in (-180, 180] -- identical to `_signed_dihedral`
    of the eye.

    🔴 THE OLD VERSION WAS SHIFTED BY 180 DEGREES.  It formed the first vector as
    `pts[a] - pts[i]`, the eye as `P[i] - P[a]` -- opposite sign, hence `n1` and `m`
    reversed, hence `atan2(-y, -x) = atan2(y, x) +- pi`.  The FOLDED magnitude stays the
    same (which is why it never showed up: 1 of 103 bonds deviated in the twist), but
    the SIGN flips: +60 of the eye was -120 here.  Measured: 87 of 103 shared bonds
    carried opposite handedness.  An enumerator that appends the missing sign thereby
    systematically appended the one ALREADY PRESENT.
    """
    b1 = pts[i] - pts[a]
    b2 = pts[j] - pts[i]
    b3 = pts[b] - pts[j]
    nb2 = float(np.linalg.norm(b2))
    if nb2 < 1e-9:
        return None
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    if float(np.linalg.norm(n1)) < 1e-9 or float(np.linalg.norm(n2)) < 1e-9:
        return None
    m = np.cross(n1, b2 / nb2)
    return math.degrees(math.atan2(float(np.dot(m, n2)), float(np.dot(n1, n2))))


def _atrop_fold(dih: float) -> float:
    """|deviation from planar| in [0, 90] -- identical to `_fold` of the eye."""
    a = abs(dih)
    return min(a, 180.0 - a)


# ===== FLANKS INSTEAD OF THE WHOLE SIDE (16.08.2026, after the `atrop44` verdict) ===========
#
# WHY IT WAS CHANGED.  `atrop44` ran additively and without damage -- `never_worse_ok = true`,
# all eleven `ccdc_*_lost = 0`, 15 frames appended on 9 systems -- and the axis
# `ccdc_atropisomer_realized` still stayed at **9x false in BOTH arms**.  That is the
# pre-registered H0 case: I append the opposite handedness of MY axes, the eye checks
# ITS OWN, and the two signatures do not coincide.
#
# The old version keyed on the WHOLE side behind the axis (element count + sorted
# elements).  The eye keys on the TWO FLANKS at the axis atom -- and that is the quantity
# that defines the handedness in the first place: two axes with equal flanks, but
# differently sized residues behind them, are the same stereogenic situation.
#
# ⚠ AND THE REFERENCE POINT OF THE DIHEDRAL ANGLE MUST FOLLOW.  The old version took the
# heavy neighbor with the SMALLEST INDEX.  The eye measures via the HIGH-RANKING flank.  With
# a different reference atom the same frame can carry the opposite sign -- my "P" would
# then be the eye's "M", and the enumeration appended the handedness already present.
# That is the subtlest part of the changeover and shows up in no self-test that only
# compares its own output.
#
# NO IMPORT FROM THE EYE.  The construction must not depend on the measuring instrument
# ("eye = UPPER BOUND of the construction"), which is why the ring detection here is built
# from the adjacency and not taken over from RDKit.


def _atrop_sub_profile(o: int, exclude: frozenset,
                       nbrs: List[List[int]], syms: Sequence[str]) -> tuple:
    """CIP-like rank of the substituent at flank atom `o`, seen AWAY from the axis --
    identical to `_sub_profile` of the eye.

    Key = (heavy atomic numbers sorted descending up to two bonds away, -nH).
    Without atom index -- a purely chemical rank, comparable across frames.

    ⚠ `_atrop_znum` instead of `_ATROP_Z.get(..., 0)`: an element outside the table gets
    the eye's covalent-radius fallback, not zero.  Exactly that is what ABOZIB hung on
    (silver: 0 here, 29 there) -- and with it the flank rank, the reference atom and the sign.
    """
    heavy_z: List[int] = []
    n_h = 0
    seen = {o}
    shell1 = [k for k in nbrs[o] if k not in exclude]
    for k in shell1:
        if syms[k] == "H":
            n_h += 1
        else:
            heavy_z.append(_atrop_znum(syms[k]))
            seen.add(k)
    for k in [x for x in seen if x != o]:
        for q in nbrs[k]:
            if q in exclude or q in seen or q == o:
                continue
            if syms[q] == "H":
                n_h += 1
            else:
                heavy_z.append(_atrop_znum(syms[q]))
    return (tuple(sorted(heavy_z, reverse=True)), -n_h)


def _atrop_flanks(i: int, partner: int, nbrs: List[List[int]],
                  syms: Sequence[str], ring_of: List[List[frozenset]]):
    """(flankHigh, flankLow, keyHigh, keyLow) or None --
    identical to `_flanks` of the eye.

    Ring atom -> the two ORTHO ring neighbors of the FIRST ring; otherwise the two heavy
    non-partners.  None if there are not exactly two DISTINGUISHABLE flanks -- then the
    axis is not stereogenic (local mirror plane) and must not be enumerated at all.
    """
    rings_i = ring_of[i]
    if rings_i:
        ring = rings_i[0]
        flanks = [k for k in nbrs[i] if k in ring and k != partner]
        excl_base = frozenset(ring)
    else:
        flanks = [k for k in nbrs[i] if k != partner and syms[k] != "H"]
        excl_base = frozenset({i})
    if len(flanks) != 2:
        return None
    f1, f2 = flanks
    ex = excl_base | {i, partner}
    k1 = _atrop_sub_profile(f1, ex, nbrs, syms)
    k2 = _atrop_sub_profile(f2, ex, nbrs, syms)
    if k1 == (tuple(), 0) and k2 == (tuple(), 0):
        return None                      # no substituent -> no rotation barrier
    if k1 == k2:
        return None                      # indistinguishable -> achiral
    return (f1, f2, k1, k2) if k1 > k2 else (f2, f1, k2, k1)


def _atrop_axis_sig_flanks(syms: Sequence[str], i: int, j: int,
                           kHiA, kLoA, kHiB, kLoB) -> tuple:
    """Canonical axis signature from the FLANKS -- the same quantity the eye carries.
    The sides are interchangeable, hence sorted; the sign is invariant under side swap,
    because dihedral(a,i,j,b) == dihedral(b,j,i,a)."""
    sideA = (syms[i], kHiA, kLoA)
    sideB = (syms[j], kHiB, kLoB)
    return tuple(sorted((repr(sideA), repr(sideB))))


def _atrop_find_axes(syms: Sequence[str], pts: np.ndarray,
                     nbrs: Optional[List[List[int]]] = None) -> List[dict]:
    """All stereogenic axes of a frame -- stage for stage those of the eye (`_axes`),
    plus the ONE quantity the eye does not need and the builder cannot do without:
    `side`, the atoms that get rotated.

    ⚠ `nbrs` is ACCEPTED AND DISCARDED.  The callers pass through the build adjacency
    (r_i+r_j+0.25); the axis detection needs the eye's (1.30*(r_i+r_j)).  The parameter
    remains only so that no caller silently receives something other than it thinks --
    it is explicitly unused.

    ⚠⚠ THE BUILDER FINDS FEWER THAN THE EYE, AND THAT IS NOT A BUG: if both sides of an
    axis hang on the metal, there is no rotation axis (`_atrop_side_atoms` -> None).  Their
    opposite handedness needs a RESEATING.  That is COUNTED and reported, not concealed --
    a silent truncation looks like completeness from the outside.
    """
    del nbrs                                  # see docstring: deliberately discarded
    P = np.asarray(pts, float)
    n = len(syms)
    if n < 6:
        return []
    nbr = _atrop_adjacency(syms, P)
    rings = _atrop_aromatic_rings(nbr, syms, P)
    ring_of: List[List[frozenset]] = [[] for _ in range(n)]
    for r in rings:
        for a in r:
            ring_of[a].append(r)
    out: List[dict] = []
    n_metallbruecke = 0
    seen_bond = set()
    for i in range(n):
        if syms[i] == "H" or _is_metal_sym(syms[i]):
            continue
        for j in nbr[i]:
            if j <= i or syms[j] == "H" or _is_metal_sym(syms[j]):
                continue
            if (i, j) in seen_bond:
                continue
            seen_bond.add((i, j))
            # An axis is a SINGLE BOND BETWEEN TWO sp2 UNITS.
            if not (_atrop_is_sp2_center(i, nbr, syms, P, ring_of)
                    and _atrop_is_sp2_center(j, nbr, syms, P, ring_of)):
                continue
            if ring_of[i] and ring_of[j] and _atrop_share_ring(i, j, ring_of):
                continue                  # fused / same ring -> no axis between units
            if not (ring_of[i] or ring_of[j]):
                continue                  # at least one aryl unit required
            if not _atrop_axis_is_single(i, j, syms, P):
                continue
            fa = _atrop_flanks(i, j, nbr, syms, ring_of)
            fb = _atrop_flanks(j, i, nbr, syms, ring_of)
            if fa is None or fb is None:
                continue                  # no two distinguishable flanks -> not stereogenic
            fHiA, _fLoA, kHiA, kLoA = fa
            fHiB, _fLoB, kHiB, kLoB = fb
            # ⚠ The dihedral goes via the HIGH-RANKING flanks -- the same choice as in the eye.
            dih = _atrop_dihedral(P, fHiA, i, j, fHiB)
            if dih is None:
                continue
            fold = _atrop_fold(dih)
            if not (_ATROP_TWIST_MIN_DEG <= fold <= _ATROP_TWIST_MAX_DEG):
                continue                  # planar or perpendicular -> sign is noise
            # FROM HERE ON the axis is an axis for the eye.  Only now does the builder ask
            # whether it can ROTATE it at all -- the order matters, otherwise a build
            # limit disappears into the axis detection.
            #
            # ⚠ NOT ROTATABLE DOES NOT MEAN NOT BUILDABLE (23.08.2026).  Formerly these
            # axes dropped out here entirely, and with that the manifold was silently
            # incomplete on them -- measured 89 of them on 120 systems.  They now stay in
            # the list, with `rotatable=False`: the opposite handedness there needs a
            # RESEATING (whole-frame mirroring), not a rotation.  Whoever filters them out
            # measures their own tool limit and takes it for chemistry.
            side = _atrop_side_atoms(syms, nbr, i, j)      # None = ring / metal bridge
            _rot = bool(side) and len(side) >= 2
            if not _rot:
                n_metallbruecke += 1
            out.append({
                "i": i, "j": j, "dih": dih, "fold": fold,
                "sign": "P" if dih > 0 else "M",
                "side": side if _rot else None,
                "rotatable": _rot,
                "sig": _atrop_axis_sig_flanks(syms, i, j, kHiA, kLoA, kHiB, kLoB),
            })
    if n_metallbruecke:
        _LOG.debug("atropisomer-enum: %d Achse(n) sind fuer das Auge Achsen, aber nicht "
                   "drehbar (Ring oder Metallbruecke) -- sie brauchen eine Neusetzung",
                   n_metallbruecke)
    return sorted(out, key=lambda d: (d["i"], d["j"]))


# ---------------------------------------------------------------------------------------------
# The mirroring
# ---------------------------------------------------------------------------------------------

def _atrop_rotate_side(pts: np.ndarray, ax: dict, angle_deg: float) -> np.ndarray:
    """Rigid rotation of the j side about the axis i->j (Rodrigues).

    Changes exclusively the torsion: bond lengths and bond angles stay exact.
    """
    p = pts.copy()
    o = pts[ax["i"]]
    k = pts[ax["j"]] - o
    nk = float(np.linalg.norm(k))
    if nk < 1e-9:
        return p
    k = k / nk
    t = math.radians(angle_deg)
    ct, st = math.cos(t), math.sin(t)
    idx = sorted(ax["side"])
    v = p[idx] - o
    p[idx] = (v * ct
              + np.cross(np.broadcast_to(k, v.shape), v) * st
              + np.outer(v @ k, k) * (1.0 - ct)) + o
    return p


def _atrop_min_nonbonded(syms: Sequence[str], pts: np.ndarray, nbrs: List[List[int]],
                         moved: Set[int]) -> float:
    """Smallest distance between a moved and a resting atom; 1-2 and 1-3 are
    excluded, because those do not change by construction."""
    worst = 9.9
    still = [k for k in range(len(syms)) if k not in moved]
    if not still:
        return worst
    for m in sorted(moved):
        excl = set(nbrs[m]) | {m}
        for q in nbrs[m]:
            excl |= set(nbrs[q])
        for s in still:
            if s in excl:
                continue
            d = float(np.linalg.norm(pts[m] - pts[s]))
            if d < worst:
                worst = d
    return worst


def _atrop_topologie_gleich(syms: Sequence[str], alt: np.ndarray, neu: np.ndarray) -> bool:
    """Did the rotation leave the BOND TOPOLOGY unchanged?

    🔴 THE SECOND ROOT (23.08.2026).  The clash floor alone is NOT enough: it was at
    1.70 A, the builder's organic bond threshold at r_i+r_j+0.25 -- for C-C **1.77 A**,
    for C-N 1.72.  A rotated atom was thus allowed to approach to 1.71 A, passed the
    clash check and got a SPURIOUS BOND.  With that `_sub_profile` changes, with it the
    axis key -- and the appended frame lands in a NEW bucket, instead of filling the one
    that is missing the sign.  Exactly that is what PEHWEH shows: the ON arm has two
    signatures for the same bond, both still with the same sign.

    Setting one number against another number would be guesswork.  What is checked is
    the PROPERTY the module relies on: "a rigid rotation about the bond axis changes
    exclusively the torsion".  If a bond partner changes in the process, that is no longer
    true -- regardless of the distance at which it happens.  Checked in BOTH adjacencies:
    the builder's (which every repairer sees later) and the eye's (which decides the
    axis key).
    """
    for adj in (lambda P: _build_geometric_adjacency(list(syms), P)[0],
                lambda P: _atrop_adjacency(syms, P)):
        a = [frozenset(x) for x in adj(alt)]
        b = [frozenset(x) for x in adj(neu)]
        if a != b:
            return False
    return True


def _atrop_mirror_frame(xyz: str, ax: dict, clash_min: float) -> Optional[str]:
    """A frame with the sign of this ONE axis flipped.

    None if it clashes OR if the rotation changes the bond topology --
    a frame in the wrong signature bucket costs exactly as much as a missing one.
    """
    syms, pts, lines = _parse_xyz(xyz)
    nbrs, _bd = _build_geometric_adjacency(syms, pts)
    new = _atrop_rotate_side(pts, ax, -2.0 * ax["dih"])
    if _atrop_min_nonbonded(syms, new, nbrs, ax["side"]) < clash_min:
        return None
    if not _atrop_topologie_gleich(syms, pts, new):
        _LOG.debug("atropisomer-enum: Drehung um %d-%d aendert die Bindungstopologie "
                   "-- Frame verworfen", ax["i"], ax["j"])
        return None
    # Mind the order: _format_xyz(orig_lines, syms, positions) -- see
    # _coord_angle_corrector.py:135.  Swapped, it delivers silent nonsense.
    return _format_xyz(lines, syms, new)


def _atrop_reseat_frame(xyz: str) -> Optional[str]:
    """RESEATING instead of rotation: the whole-frame mirroring of the frame.

    WHAT FOR.  If BOTH sides of an axis hang on the metal -- a chelating biaryl --, the
    opposite handedness can NOT be produced by rigid rotation: it would tear the M-D
    bond apart.  Measured on 23.08. on 120 systems: **89 axes** that the eye carries and
    the builder cannot rotate.  Until now they silently dropped out.

    The whole-frame mirroring is the right operation, and for a reason, not out of
    convenience: it is an ISOMETRY.  Every distance is preserved exactly -- all bond
    lengths, all angles, every M-D bond too.  It therefore cannot tear anything apart by
    construction, and it flips the sign of EVERY stereogenic axis at once.

    ⚠ IT ALSO FLIPS EVERYTHING ELSE: Lambda/Delta at the metal, backbone centers.  That
    is not collateral damage but a DIFFERENT isomer -- and appended additively that means
    completeness, not replacement.  The original remains untouched.

    NO SECOND MIRRORING.  `_mirror_enum.mirror_frame` does exactly that, including the
    achirality check (if the frame is its own mirror image, there is nothing to gain and
    the function returns None).  Here it is used, not rebuilt.
    """
    try:
        from delfin.manta._mirror_enum import mirror_frame as _mf
    except Exception as exc:                       # pragma: no cover
        _LOG.debug("atropisomer-enum: Spiegelung nicht verfuegbar: %s", exc)
        return None
    try:
        return _mf(xyz)
    except Exception as exc:
        _LOG.debug("atropisomer-enum: Neusetzung fehlgeschlagen: %s", exc)
        return None


def _atrop_realized(xyz: str, sig: tuple, want: str) -> bool:
    """Does this frame carry the signature `sig` with the sign `want` -- REALLY?

    🔑 THE QUESTION THAT WAS MISSING FOR THE WHOLE OF 23.08.  `atrop44` and `atrop10k`
    appended frames and the axis stayed put; nobody ever checked whether the appended
    frame carries the sign it was supposed to carry at all.  For PEHWEH it did not --
    it landed in a NEW signature bucket, with the same sign.

    A frame that does not fulfil its purpose is no gain: it costs build and eye time,
    raises maxima and contributes nothing.  That is why here it is CHECKED after the
    operation, not hoped for beforehand.
    """
    try:
        A = _atrop_analyze(xyz)
    except Exception:
        return False
    if not A:
        return False
    return any(ax["sig"] == sig and ax["sign"] == want for ax in A["axes"])


def _atrop_analyze(xyz: str) -> Optional[dict]:
    syms, pts, _lines = _parse_xyz(xyz)
    if len(syms) < 4:
        return None
    axes = _atrop_find_axes(syms, pts)
    return {"axes": axes} if axes else None


# ---------------------------------------------------------------------------------------------
# The additive pass
# ---------------------------------------------------------------------------------------------

def expand_atropisomers(results):
    """ADDITIVE: append the MISSING sign for every stereogenic axis.

    Pass 1 collects over ALL frames which (axis signature, sign) occur -- exactly the
    point at which `_stereocenter_enum` failed on 10.08.: it reported completeness on a
    COARSER partition than the one it completed.  Here the key is the AXIS itself, i.e.
    the quantity at stake.

    Pass 2 builds, for every signature missing a sign, ONE frame from the
    representative.  Originals remain untouched -> never-worse by construction.
    """
    if not results:
        return results
    max_added = _atrop_env_int("DELFIN_ATROPISOMER_MAX_ADDED", 64)
    clash_min = _atrop_env_float("DELFIN_ATROPISOMER_CLASH_MIN", 1.70)

    present: Set[tuple] = set()                 # (sig, sign) already in the manifold
    reps: Dict[tuple, tuple] = {}               # sig -> (order, xyz, label, axis)
    for order, entry in enumerate(results):
        try:
            xyz = entry[0]
            lbl = entry[1] if len(entry) > 1 else ""
        except Exception:
            continue
        try:
            A = _atrop_analyze(xyz)
        except Exception as exc:
            _LOG.debug("atropisomer-enum: Frame %d nicht lesbar: %s", order, exc)
            continue
        if not A:
            continue
        for ax in A["axes"]:
            present.add((ax["sig"], ax["sign"]))
            reps.setdefault(ax["sig"], (order, xyz, lbl, ax))

    if not reps:
        return results

    # RESEATING for non-rotatable axes -- its own switch, default OFF, so that the rotation
    # path stays measurable unchanged and `atrop2` measures exactly what was pre-registered.
    # A second mechanism in the same run would be a second axis in the A/B.
    _reseat_on = _atrop_env_int("DELFIN_ATROPISOMER_RESEAT", 0) == 1
    added: List[tuple] = []
    capped = 0
    n_unverifiziert = 0
    n_nicht_drehbar = 0
    for sig, (order, xyz, lbl, ax) in sorted(reps.items(), key=lambda kv: kv[1][0]):
        want = "M" if ax["sign"] == "P" else "P"
        if (sig, want) in present:
            continue                             # both hands already there
        if len(added) >= max_added:
            capped += 1
            continue
        mx = None
        how = "atrop"
        if ax.get("rotatable", True):
            try:
                mx = _atrop_mirror_frame(xyz, ax, clash_min)
            except Exception as exc:
                _LOG.debug("atropisomer-enum: Drehung fehlgeschlagen (%s): %s", lbl, exc)
                mx = None
        else:
            # Both sides on the metal: a rotation would tear the M-D bond apart.
            # ⚠ NO FALLBACK FROM THE ROTATION TO THE MIRRORING.  A rotation that FAILS on
            # clash or topology is a justified refusal -- covering it up with an
            # operation that incidentally flips every other stereo element would be
            # cosmetics.  The reseating applies only where rotating is not POSSIBLE
            # at all.
            n_nicht_drehbar += 1
            if _reseat_on:
                mx = _atrop_reseat_frame(xyz)
                how = "atropR"
        if mx is None:
            continue                             # do not build, lose nothing
        # ===== VERIFICATION: DOES THE FRAME DO WHAT IT WAS BUILT FOR? ===================
        if not _atrop_realized(mx, sig, want):
            n_unverifiziert += 1
            _LOG.debug("atropisomer-enum: %s traegt %s auf %r NICHT -- verworfen",
                       lbl, want, sig)
            continue
        added.append((mx, f"{lbl}_{how}-{want}"))
        # The whole-frame mirroring flips ALL axes at once.  Whatever it realizes along
        # the way is recorded as present -- otherwise the next round builds the same
        # frame once more for a different signature.
        if how == "atropR":
            try:
                for _a in (_atrop_analyze(mx) or {}).get("axes", []):
                    present.add((_a["sig"], _a["sign"]))
            except Exception:
                pass
        present.add((sig, want))
    if n_unverifiziert:
        _LOG.warning("atropisomer-enum: %d Frame(s) gebaut und VERWORFEN, weil sie das "
                     "fehlende Vorzeichen nicht trugen", n_unverifiziert)
    if n_nicht_drehbar and not _reseat_on:
        _LOG.info("atropisomer-enum: %d Achse(n) nicht drehbar (Metallbruecke); "
                  "DELFIN_ATROPISOMER_RESEAT=1 setzt sie per Neusetzung", n_nicht_drehbar)

    if capped:
        # NO SILENT TRUNCATION: what the cap takes away is reported.
        _LOG.warning("atropisomer-enum: Deckel DELFIN_ATROPISOMER_MAX_ADDED=%d erreicht, "
                     "%d Achsen NICHT vervollstaendigt", max_added, capped)
    if not added:
        return results
    _LOG.info("atropisomer-enum: %d Gegenhaendigkeiten ergaenzt (%d Achsen gesamt)",
              len(added), len(reps))

    out = list(results)
    proto = results[0]
    for mx, lab in added:
        if isinstance(proto, tuple):
            out.append((mx, lab))
        elif isinstance(proto, list):
            out.append([mx, lab])
        else:
            out.append(mx)
    return out


def _atrop_why(syms, pts, nbrs, i, j) -> str:
    """At WHICH stage does the bond i-j fail?  A self-test that only reports 'nothing
    found' is worthless -- it must name the stage.

    ⚠ The stages stand here in THE SAME order as in `_atrop_find_axes`; if the two drift
    apart, the self-test names a stage at which it does not fail at all.
    `nbrs` is discarded here as well -- the check is against the adjacency of the eye.
    """
    del nbrs
    P = np.asarray(pts, float)
    if syms[i] == "H" or syms[j] == "H":
        return "H"
    if _is_metal_sym(syms[i]) or _is_metal_sym(syms[j]):
        return "Metall"
    nbr = _atrop_adjacency(syms, P)
    rings = _atrop_aromatic_rings(nbr, syms, P)
    ring_of: List[List[frozenset]] = [[] for _ in range(len(syms))]
    for r in rings:
        for a in r:
            ring_of[a].append(r)
    sp2i = _atrop_is_sp2_center(i, nbr, syms, P, ring_of)
    sp2j = _atrop_is_sp2_center(j, nbr, syms, P, ring_of)
    if not (sp2i and sp2j):
        return f"nicht sp2 (i={sp2i} deg={len(nbr[i])}, j={sp2j} deg={len(nbr[j])})"
    if ring_of[i] and ring_of[j] and _atrop_share_ring(i, j, ring_of):
        return "kondensiert / selber Ring"
    if not (ring_of[i] or ring_of[j]):
        return "keine Seite ist ein aromatischer Ring"
    if not _atrop_axis_is_single(i, j, syms, P):
        d = float(np.linalg.norm(P[i] - P[j]))
        return f"keine Einfachbindung (d={d:.2f})"
    fa = _atrop_flanks(i, j, nbr, syms, ring_of)
    fb = _atrop_flanks(j, i, nbr, syms, ring_of)
    if fa is None:
        return "Seite i: keine zwei unterscheidbaren Flanken"
    if fb is None:
        return "Seite j: keine zwei unterscheidbaren Flanken"
    d = _atrop_dihedral(P, fa[0], i, j, fb[0])
    if d is None:
        return "Dieder nicht berechenbar"
    f = _atrop_fold(d)
    if not (_ATROP_TWIST_MIN_DEG <= f <= _ATROP_TWIST_MAX_DEG):
        return f"Verdrillung {f:.1f} ausserhalb {_ATROP_TWIST_MIN_DEG}-{_ATROP_TWIST_MAX_DEG}"
    side = _atrop_side_atoms(syms, nbr, i, j)
    if not side:
        return "Bindung liegt im Ring oder ueberbrueckt das Metall -- nicht drehbar"
    if len(side) < 2:
        return "Seite zu klein"
    return "OK"


if __name__ == "__main__":  # pragma: no cover -- axis comparison against the eye
    # WHAT FOR.  `atrop44` measured on 16.08. that the enumeration appends frames and the
    # axis `ccdc_atropisomer_realized` nevertheless stays at zero: my axis signature was
    # not that of the eye.  This entry point prints the axes found in the same form as
    # `weddell/detectors/atropisomer_sign.py`, so that both lists can be read SIDE BY
    # SIDE -- atom pair, twist, sign.  If they do not coincide, a further A/B run is
    # pointless.
    #
    #     python delfin/manta/_atropisomer_enum.py <frame.xyz> [more.xyz ...]
    import sys as _sys
    for _p in _sys.argv[1:]:
        try:
            with open(_p) as _fh:
                _txt = _fh.read()
        except Exception as _e:
            print(f"  {_p}: nicht lesbar ({_e})")
            continue
        # Multi-frame XYZ: only the first frame, as the eye's command line does too.
        _lines = _txt.splitlines()
        try:
            _n = int(_lines[0].split()[0])
            _one = "\n".join(_lines[:_n + 2])
        except Exception:
            _one = _txt
        _A = _atrop_analyze(_one)
        print(f"\n=== {_p} ===")
        if not _A:
            print("  keine stereogene Achse gefunden -- Ablehnungsgruende je Kandidat:")
            _s, _P, _l = _parse_xyz(_one)
            _nb, _bd = _build_geometric_adjacency(_s, _P)
            _seen_reason = {}
            for _i in range(len(_s)):
                for _j in _nb[_i]:
                    if _j <= _i:
                        continue
                    _r = _atrop_why(_s, _P, _nb, _i, _j)
                    if _r in ("H", "Metall"):
                        continue
                    _seen_reason.setdefault(_r, []).append(f"{_s[_i]}{_i}-{_s[_j]}{_j}")
            for _r, _ex in sorted(_seen_reason.items(), key=lambda kv: -len(kv[1])):
                print(f"    {len(_ex):4d}x  {_r}   z.B. {' '.join(_ex[:3])}")
            continue
        for _ax in _A["axes"]:
            print(f"  {_ax['i']}-{_ax['j']}  twist {_ax['fold']:5.1f}  "
                  f"sign {_ax['sign']}  sig {_ax['sig']}")


