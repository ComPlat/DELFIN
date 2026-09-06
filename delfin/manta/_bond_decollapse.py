"""Iter-25 (2026-05-20) — final bond-decollapse corrector.

Calibration breakthrough (project_calibration_breakthrough_2026_05_20): ~79% of
hapto structures emit COLLAPSED ligands — heavy-heavy pairs at 0.24-1.2 Å
(atoms overlapping / superimposed phenyls / substituents fused onto their
parent C).  The collapse originates in several builder paths (scaffold
substituent placement, BFS ring building, post-UFF on unparametrized metals);
rather than chase every origin, this is a FINAL guaranteed pass that runs on
the emitted XYZ after all builders/UFF/bandages and pushes any collapsed
non-metal geometry back to physical bond lengths.

Mechanism (per frame): a short spring relaxation on the heavy-atom graph —
  * bonded heavy-heavy pairs are pulled toward their element-pair ideal length
  * non-bonded heavy pairs closer than a floor are pushed apart (resolves
    superimposed phenyl rings)
  * metals AND every atom within M-D bonding distance are FROZEN, so the
    coordination sphere / η-ring / M-D invariant is never disturbed
H atoms are dragged rigidly with their heavy parent.

Per-frame rollback: only keep the relaxed frame if it strictly REDUCES the
collapsed-bond count and introduces no new metal-donor break.  Bit-exact when
no collapsed bond is present (early return).  XRD-validated target: drop the
collapse rate from ~79%.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Tuple

import numpy as np

# ===== THE METAL LIST WAS A TRUNCATED HAND COPY (2026-08-14) =====
# Canonical is smiles_converter._METALS with 68 elements.  Here there were 37, and the
# 33 missing ones are no exotics but whole blocks:
#     s-block      Li Na K Rb Cs Be Mg Ca Sr Ba          (complete)
#     lanthanides  Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb   (12 of 15; only La/Ce/Lu were there)
#     actinides    Ac Th Pa U Np Pu                      (COMPLETELY missing)
#     p-block      Al Ga In Tl Po
#
# WHAT DAMAGE THIS DOES: _is_metal("U") returned False.  A U-O coordination bond thus counted
# as an ordinary heavy-heavy bond, got an invented radius via _COV.get(s, 0.9) and was
# "repaired" by the decollapse corrector -- those are the 32 invented bonds, 22 of them
# uranyl.  A metal that is not recognised as a metal is also not frozen (the frozen set
# further down hangs on _is_metal).
#
# ⚠ NOT imported, because smiles_converter in turn pulls manta modules (cycle).
# The list is therefore kept here AND checked by harness/polyhedra_audit.py against the
# canonical one, so that the copy cannot drift apart unnoticed again.
_METALS_HIST = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
                   "Hf Ta W Re Os Ir Pt Au Hg La Ce Lu Sn Pb Ge Sb Bi".split())
_METALS_ADDED = set("Li Na K Rb Cs Be Mg Ca Sr Ba "
                    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
                    "Ac Th Pa U Np Pu Al Ga In Tl Po".split())
_METALS = _METALS_HIST | (
    _METALS_ADDED if os.environ.get("DELFIN_FFFREE_DECOLLAPSE_METALS", "0") == "1" else set())
_COV = {'C':0.76,'N':0.71,'O':0.66,'H':0.31,'S':1.05,'Cl':1.02,'P':1.07,'F':0.57,
        'Br':1.20,'I':1.39,'B':0.84,'Si':1.11,'Se':1.20,'As':1.19,'Te':1.38}
_MD_TOL = 0.05          # M-D invariant tolerance (Å)
_MD_FACTOR = 1.30       # M-D bond detection factor
_NONBOND_FLOOR = 0.78   # non-bonded heavy pair < this·Σcov → push apart
# ===== THE COLLAPSE FLOOR, SINCE 2026-08-11 IN ONE PLACE AND SWITCHABLE =====
# It stood as the literal 0.82 in THREE places (here in _count_collapsed, in
# converter_backend._coll_floor, as the default of assemble_complex._collapsed_heavy_bonds_strict)
# and was thus neither findable nor measurable.  Default unchanged 0.82 -> byte-identical.
#
# WHY IT NEEDS TO BE MEASURED.  Measured on 11.08. on 8974 CCDC crystals with 287043
# heavy-atom bonds (radii and bond detection exactly as here):
#     distribution d/Sigma_cov   p0.1 = 0.768   p1 = 0.808   median 0.919   minimum 0.511
# The self-gate rejects an ENTIRE config as soon as ONE bond lies below it.  Counted per
# crystal that means:
#     floor 0.82 -> 1558 of 8974 crystals affected = 17,4 %
#           0.78 ->  549 = 6,1 %      0.75 -> 77 = 0,9 %      0.72 -> 34 = 0,4 %
# The floor would thus reject every sixth EXPERIMENTALLY DETERMINED crystal as collapsed
# -- against the standing rule that the eye stays silent on clean CCDC structures.  Two
# independent measurements from the same week fit this: 94,5 % of ALL self-gate
# rejections are COLLAPSED_BOND, and 193 of 996 systems die on this one number.
# 17,4 % against 19,4 % -- that matches.
#
# ⚠ AND THE CODE CONTRADICTED ITSELF: refine.py:93 builds against _COLLAPSE = 0.70, the
# non-bond floor here stands at 0.78, judgement was passed with 0.82.  The builder was
# measured more strictly than it built.
#
# ⚠ THE 0.82 IS THUS NOT REFUTED.  The crystal measurement says how many REAL bonds the
# floor hits -- not how many FALSE ones it still catches.  Lowering it on that basis
# alone would be exactly the Goodhart error this repo otherwise hunts.  The decision is
# made by an A/B 0.82 against 0.75 under UNION, criterion ccdc_isomer_realized and
# defect-free manifolds -- not the frame count.
#
# TWO SWITCHES, and the reason is the harness mechanics: `loop.py --on FLAG` always sets
# a switch to "1" and therefore cannot switch a VALUE.  An A/B over a numeric value
# therefore needs a BOOLEAN axis -- otherwise the value would have to stand in BOTH arms,
# and that is exactly where cap10swingunion failed ("UNDECLARED AXIS: both arms carry
# it, the A/B cannot see it").
#   DELFIN_FFFREE_COLLAPSE_FLOOR_CAL=1    -> the crystal-calibrated floor 0.75  (A/B axis)
#   DELFIN_FFFREE_COLLAPSE_FLOOR=<number> -> overrides both  (manual sweep)
# 0.75 is not chosen freely: it is the threshold with 0,9 % false positives on the 8974
# crystals AND at the same time the lower bound of the distortion band in refine.py:139
# (dev < -0.25).  At this value the code says the same thing in three places.
_CAL = os.environ.get("DELFIN_FFFREE_COLLAPSE_FLOOR_CAL", "0") == "1"
try:
    COLLAPSE_FLOOR = float(os.environ.get("DELFIN_FFFREE_COLLAPSE_FLOOR",
                                          "0.75" if _CAL else "0.82"))
except Exception:
    COLLAPSE_FLOOR = 0.75 if _CAL else 0.82
_MAX_STEP = 0.30        # max per-atom displacement per relaxation pass (Å) — keeps
                        # the vdw-aware repulsion from blowing coordinates up to NaN

# Iter-26 (Task #32): vdw/angle-awareness so the decollapse pass cannot trade a
# collapse fix for vdw-clash / H-planarity / angle regressions (the Iter-25
# concentrated trade vs iter23: vdw +20%, F3_angle +33%, F20 +3.6pp).
# Bondi vdW radii — verbatim from quality_framework/scripts/nonbonded_vdw_check.py
# so the in-corrector clash proxy matches the detector bit-for-bit on radii.
_VDW = {
    "H": 1.20, "He": 1.40, "Li": 1.81, "Be": 1.52, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.15, "Ti": 2.11, "V": 2.07, "Cr": 2.06,
    "Mn": 2.05, "Fe": 2.04, "Co": 2.00, "Ni": 1.97, "Cu": 1.96, "Zn": 2.01,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.32, "Zr": 2.23, "Nb": 2.18, "Mo": 2.17,
    "Tc": 2.16, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11, "Cd": 2.18,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43, "Hf": 2.23, "Ta": 2.22, "W": 2.18,
    "Re": 2.16, "Os": 2.16, "Ir": 2.13, "Pt": 2.13, "Au": 2.14, "Hg": 2.23,
    "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
}
_VDW_DEFAULT = 1.7      # detector default for missing element
_VDW_FACTOR = 0.85      # clash if d < factor·(vdW_a + vdW_b)  (matches detector)
_ANGLE_MIN_DEG = 70.0   # optional angle proxy: heavy-heavy-heavy angle < this = bad
_F20_OOP_TOL = 0.20     # ring-H out-of-plane tolerance (Å)   (matches detector)
_F20_RING_TOL = 0.10    # ring planarity tolerance (Å)        (matches detector)
_F20_AROMATIC = {"C", "N", "O", "S"}
# ===== THE THRESHOLD STOOD THREE TIMES IN THE TREE (19.08.2026) =======================
# Here stood a literal 1.46 of its own with the comment "matches
# _arom_planarize._AROMATIC_BOND_MAX" -- i.e. a hand-maintained alignment, and exactly
# such a thing drifts apart: on 18.08. the sp3 repair hit two of the three copies, this
# third one not.  The number now comes from ``_arom_criterion``, and that is also where
# the switch ``DELFIN_FFFREE_AROM_CRITERION_RADII`` hangs (default 0 -> byte-identical).
#
# ⚠ THE IMPORT DOES NOT BREAK THIS MODULE'S CYCLE-FREENESS.  ``_arom_criterion`` pulls
# ONLY ``os`` and ``typing`` at module level; it fetches the radii table only at the
# first normalised ring evaluation.  This import therefore also shifts no import-time
# env read (``polyhedra`` reads DELFIN_FFFREE_COV_COMPLETE at import, and
# _bond_decollapse is pulled early -- which is why the lazy construction is mandatory there).
# ``_AROMATIC_BOND_MAX`` is only passed through (the name is preserved),
# the decision is made in ``ring_rejected_by_length``.
from delfin.manta._arom_criterion import (   # noqa: E402,F401
    _AROMATIC_BOND_MAX,
    ring_rejected_by_length,
)
_AROM_M_COORD_DIST = 2.6    # atom within this of a metal counts as coordinated;
                            # a ring TOUCHING the coordination sphere is never
                            # aromatic-retargeted (LUMTAP never-worse scope)

# Iter 31 (aromatic bond-length ROOT seat, DELFIN_FFFREE_AROM_SEAT): when on, the
# per-bond SPRING TARGET of this final decollapse pass uses the delocalised CCDC
# aromatic length for bonds KNOWN to be aromatic ring bonds, so the spring cannot
# re-elongate an aromatic ring bond that refine/UFF already seated back toward the
# single-bond covalent sum (~1.52 for C–C).  Default OFF → byte-identical.  ONLY
# the spring target consults this (via aromatic=True); the bond-graph, collapse
# and M–D primitives keep the covalent-sum ideal (aromatic=False) so connectivity
# perception is unchanged.
_AROM_SEAT = os.environ.get("DELFIN_FFFREE_AROM_SEAT", "0") == "1"


def _ideal_bond(a: str, b: str, aromatic: bool = False) -> float:
    if aromatic and _AROM_SEAT:
        try:
            from delfin.fffree.aromatic_bond_targets import aromatic_target
            t = aromatic_target(a, b)
        except Exception:
            t = None
        if t is not None:
            return t
    # THE METALS THAT NEVER STOOD HERE (switch and table in `_elements`, default
    # OFF -> byte-identical; the read site there remains the ONLY one).  `_COV` below
    # carries 15 elements and NO metal: for Fe-N it yields 0.90 + 0.71 = 1.61 A --
    # an invented number that carries the SELF-GATE and that the seven eye detectors
    # compensate for through inflated factors.
    # ⚠ WITH THE SWITCH THESE FACTORS ARE DOUBLY COMPENSATED: metric_md_direction
    # (1.40), metric_donor_collapse / metric_h_axis / metric_md_angle_realism (1.45)
    # and metric_coord_geom / metric_coord_shape (1.65) read the same function and
    # then reach TOO FAR.  This is NOT a null test -- the first run measures the COUPLING,
    # not the benefit, and the factors belong back at 1.30 in the same step.
    from delfin.manta import _elements as _EL
    if _EL.cov_radii_enabled():
        return _EL.covalent_radius(a) + _EL.covalent_radius(b)
    return _COV.get(a, 0.9) + _COV.get(b, 0.9)


def _is_metal(s: str) -> bool:
    # ONE SOURCE (14.08.2026): delfin/manta/_elements.py, 68 elements.  Default OFF ->
    # byte-identical.  This spot has the greatest leverage in the whole tree: around 25
    # eye predicates point here WITHOUT a fallback of their own (eye.py:61 and the
    # metric_*), and by default it knows only 37 elements -- in the standard run Na, K,
    # Ca, Mg, U, Th and Al are therefore not metals in the EYE.  The import sits inside
    # the function because a module import built a cycle here; after the first call it
    # costs nothing.
    from delfin.manta import _elements as _EL
    if _EL.unified_enabled():
        return _EL.is_metal(s)
    return s in _METALS


def _norm_iso(sym: str) -> str:
    """Map D/T onto H -- ONLY for the geometry logic (2026-08-14).

    Deuterium was nowhere recognised as hydrogen: `_COV` knows no 'D', so
    `_COV.get('D', 0.9)` gave a radius of 0.9 instead of 0.31 -- almost three times --,
    and all queries `syms[i] == 'H'` (H parent assignment, vdW counting, angle gate,
    collapse counting) missed it.  A D atom thus counted as a HEAVY ATOM with an
    invented radius: it attracted phantom bonds and triggered collapse findings that
    are not real.  Affected: 219 of 307370 crystals.

    ⚠ The OUTPUT stays untouched: the reassembly below writes `p[0]`, i.e. the
    ORIGINAL SYMBOL from the input line.  A D does not become an H here, it is only
    finally COMPUTED as one.
    """
    if sym in ("D", "T") and os.environ.get("DELFIN_FFFREE_DECOLLAPSE_ISOTOPES", "0") == "1":
        return "H"
    return sym


def _parse(xyz: str):
    lines = xyz.splitlines()
    syms: List[str] = []
    pts: List[List[float]] = []
    keep: List[str] = []
    for ln in lines:
        p = ln.split()
        if len(p) == 4:
            try:
                xyz_v = [float(p[1]), float(p[2]), float(p[3])]
            except ValueError:
                keep.append(ln)
                continue
            syms.append(_norm_iso(p[0])); pts.append(xyz_v)
        else:
            keep.append(ln)
    return syms, np.array(pts, dtype=float), lines


def _count_collapsed(syms, P, bonds) -> int:
    n = 0
    for i, j in bonds:
        if _is_metal(syms[i]) or _is_metal(syms[j]):
            continue
        if float(np.linalg.norm(P[i] - P[j])) < COLLAPSE_FLOOR * _ideal_bond(syms[i], syms[j]):
            n += 1
    return n


def _md_snapshot(syms, P) -> Dict[Tuple[int, int], float]:
    snap = {}
    n = len(syms)
    for i in range(n):
        if not _is_metal(syms[i]):
            continue
        for j in range(n):
            if i == j or _is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d <= _MD_FACTOR * _ideal_bond(syms[i], syms[j]):
                snap[(i, j)] = d
    return snap


def _md_broken(snap, P) -> bool:
    for (i, j), d0 in snap.items():
        if abs(float(np.linalg.norm(P[i] - P[j])) - d0) > _MD_TOL:
            return True
    return False


def _geometric_bonds(syms, P) -> List[Tuple[int, int]]:
    """Heavy-heavy + X-H pairs within 1.3·Σcov are treated as bonds.  Atom-order
    independent (the emitted XYZ may not match the mol's atom indexing), and
    collapse-robust (a 0.24 Å pair is still < 1.3·Σcov so it is captured)."""
    n = len(syms)
    bonds: List[Tuple[int, int]] = []
    for i in range(n):
        if _is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if _is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * _ideal_bond(syms[i], syms[j]):
                bonds.append((i, j))
    return _cut_geminal(syms, P, bonds)


def _cut_geminal(syms, P, bonds):
    """Throw geminal 1-3 pairs out of the bond list -- via the ANGLE, not the distance.

    THE PROBLEM.  Above, `d < 1.30 * Sigma_cov` alone decides.  On an IDEAL geometry
    that is harmless: a C-C-C pair only drops below the threshold under ~80 degrees.
    But this module works precisely on COLLAPSED frames -- and when the bonds are
    short, the 1-3 distance shrinks with them.  With 1.2-A bonds an 85-degree pair
    already lies at 1.62 A and is counted as a bond.  Then _count_collapsed counts this
    PHANTOM BOND as a collapse, and the corrector "repairs" something that is not a
    bond at all.

    WHY THE ANGLE AND NOT THE THRESHOLD.  On 11.08. the scalar 1.30 was swept: NO value
    reaches 0 false findings at 0 losses, because real bonds extend up to
    1.290 * Sigma_cov -- there is no gap.  The error is STRUCTURAL, not a matter of the
    threshold.  The angle, by contrast, SURVIVES THE COLLAPSE: if all bonds shrink by
    the same factor, A-X-B stays unchanged.  It is thus the only quantity here that
    still says something when the lengths no longer say anything.

    THE SEPARATION.  A REAL three-membered ring (cyclopropane) has ~60 degrees at the
    apex and its A-B bond MUST stay.  A geminal pair has >= 85 degrees.  Measured, there
    is NOTHING between 82.0 and 84.4 degrees -- the separation is empty, so the cut is
    not calibrated but read off.  Effect: `smiles_topology` 6/509 -> 0.

    Default OFF -> byte-identical.  Threshold sweepable via DELFIN_FFFREE_GEMINAL_CUT_DEG.
    """
    if os.environ.get("DELFIN_FFFREE_GEMINAL_CUT", "0") != "1":
        return bonds
    cut = float(os.environ.get("DELFIN_FFFREE_GEMINAL_CUT_DEG", "85"))
    adj: Dict[int, set] = {}
    for i, j in bonds:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    drop = set()
    for i, j in bonds:
        # Common neighbours = triangles in the candidate graph.  ONLY there can a
        # 1-3 pair appear as a bond at all; 1-4 and beyond share none.
        for k in adj.get(i, ()) & adj.get(j, ()):
            va, vb = P[i] - P[k], P[j] - P[k]
            na = float(np.linalg.norm(va)); nb = float(np.linalg.norm(vb))
            if na < 1e-9 or nb < 1e-9:
                continue
            c = float(np.dot(va, vb)) / (na * nb)
            ang = math.degrees(math.acos(max(-1.0, min(1.0, c))))
            if ang >= cut:
                drop.add((i, j))
                break
    return [b for b in bonds if b not in drop]


def _aromatic_ring_bonds(syms, P, bonds) -> set:
    """Set of (min,max) heavy–heavy bonds lying in a geometric 5/6-ring of
    aromatic-eligible atoms (C/N/O/S) whose mean intra-ring bond is in the
    aromatic band; used only to select the aromatic spring target; returns the
    empty set unless the seat flag is on (caller guards).

    ⚠ THE RING FINDER REMAINS SEPARATE, ONLY THE GATE IS SHARED (19.08.2026).
    This copy sees a DIFFERENT ring set than ``_arom_planarize``, for two
    demonstrable reasons:
      * the bond perception is a different one -- here ``d < 1.30 * _ideal_bond``
        (for C-C thus 1.98 A), there ``d < Sigma_r_cov + 0.25`` (1.77 A).  This
        copy therefore systematically sees MORE edges and more rings.
      * every ring that TOUCHES the coordination sphere drops out entirely here
        (LUMTAP never-worse scope) -- there it is anchored and treated as well.
    The ring finder therefore stays where it is; shared from now on is solely the
    question "is this ring-length signature aromatic", and that lives in
    ``_arom_criterion.ring_rejected_by_length``."""
    n = len(syms)
    arom_nbr: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        if (syms[i] in _F20_AROMATIC and syms[j] in _F20_AROMATIC
                and not _is_metal(syms[i]) and not _is_metal(syms[j])):
            arom_nbr[i].append(j); arom_nbr[j].append(i)
    # coordination sphere (atoms within _AROM_M_COORD_DIST of a metal) + full heavy
    # +H adjacency, so a ring that TOUCHES it (a ring atom is, or is bonded to, a
    # coordinated atom) is left byte-identical: reshaping a coordinated ring distorts
    # the coordination geometry the eye perceives (the LUMTAP never-worse scope that
    # the refine seat also uses).
    coord: set = set()
    metals = [i for i in range(n) if _is_metal(syms[i])]
    for m in metals:
        for a in range(n):
            if a != m and float(np.linalg.norm(P[a] - P[m])) < _AROM_M_COORD_DIST:
                coord.add(a)
    full_adj: List[List[int]] = [[] for _ in range(n)]
    if coord:
        for i, j in bonds:
            full_adj[i].append(j); full_adj[j].append(i)
    rings: set = set()
    for start in range(n):
        if syms[start] not in _F20_AROMATIC:
            continue
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 6:
                continue
            for nx in arom_nbr[cur]:
                if nx == path[0] and len(path) >= 5:
                    rings.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 6:
                    stack.append((nx, path + [nx]))
    out: set = set()
    for ring in rings:
        rset = set(ring)
        if coord and any((a in coord) or any(nb in coord for nb in full_adj[a])
                         for a in ring):
            continue                       # coordinated ring — leave byte-identical
        lens: List[float] = []
        edges: List[Tuple[int, int]] = []
        for i in ring:
            for j in arom_nbr[i]:
                if j in rset and j > i:
                    lens.append(float(np.linalg.norm(P[i] - P[j])))
                    edges.append((i, j))
        # ONE SOURCE (see import above): switch OFF = old comparison
        # ``mittel >= 1.46``, switch ON = mean of d/(r_i+r_j) >= 0.939.
        if not lens or ring_rejected_by_length(syms, edges, lens):
            continue                       # saturated / non-aromatic ring
        for (i, j) in edges:
            out.add((min(i, j), max(i, j)))
    return out


def _neighbor_exclusion(n: int, bonds: List[Tuple[int, int]]) -> List[set]:
    """1-2 + 1-3 neighbor exclusion sets derived purely from the geometric
    bonds (no SMILES graph).  Mirrors the detector's 1-3 closure
    (nonbonded_vdw_check._smiles_to_atom_pairs): excl[i] = adj[i] ∪
    (⋃_{j∈adj[i]} adj[j]) − {i}.  Used by both the vdw proxy and the
    pair-aware repulsion floor so true steric clashes are treated
    identically to the detector while 1-3 contacts are left alone."""
    adj: List[set] = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    excl: List[set] = []
    for i in range(n):
        s = set(adj[i])
        for j in adj[i]:
            s.update(adj[j])
        s.discard(i)
        excl.append(s)
    return excl


def _count_vdw_clashes(syms, P, excl) -> int:
    """Geometry-only reproduction of nonbonded_vdw_check.vdw_clash_report:
    non-bonded (not 1-2/1-3) heavy pair is a clash if d < 0.85·(vdW_a+vdW_b).
    Heavy-only, metal-metal pairs skipped."""
    n = len(syms)
    cnt = 0
    for i in range(n):
        if syms[i] == 'H':
            continue
        mi = _is_metal(syms[i])
        ri = _VDW.get(syms[i], _VDW_DEFAULT)
        ex = excl[i]
        for j in range(i + 1, n):
            if syms[j] == 'H' or j in ex:
                continue
            if mi and _is_metal(syms[j]):
                continue
            min_d = _VDW_FACTOR * (ri + _VDW.get(syms[j], _VDW_DEFAULT))
            if float(np.linalg.norm(P[i] - P[j])) < min_d:
                cnt += 1
    return cnt


def _count_h_planar_viol(syms, P) -> int:
    """F20 proxy — port of compute_ligand_angles.sp2_h_planarity_check onto
    (syms, P).  Find 5/6-rings of aromatic-eligible heavy atoms (C/N/O/S),
    confirm ring planar (≤ _F20_RING_TOL), count ring-bonded H whose OOP
    distance from the ring plane exceeds _F20_OOP_TOL.  Self-contained
    distance-bond adjacency (includes H) so it matches the detector."""
    n = len(syms)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        ri = _COV.get(syms[i], 0.76)
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * (ri + _COV.get(syms[j], 0.76)):
                nbrs[i].append(j)
                nbrs[j].append(i)
    heavy_only: List[List[int]] = [
        [j for j in nbrs[i] if syms[j] != 'H' and syms[i] in _F20_AROMATIC
         and syms[j] in _F20_AROMATIC]
        for i in range(n)
    ]
    rings: set = set()
    for start in range(n):
        if syms[start] not in _F20_AROMATIC:
            continue
        stack = [(start, [start])]
        while stack:
            cur, path = stack.pop()
            if len(path) > 6:
                continue
            for nx in heavy_only[cur]:
                if nx == path[0] and len(path) >= 5:
                    rings.add(tuple(sorted(path)))
                    continue
                if nx in path:
                    continue
                if len(path) < 6:
                    stack.append((nx, path + [nx]))
    cnt = 0
    for ring in rings:
        ring_pts = np.array([P[i] for i in ring])
        centroid = ring_pts.mean(axis=0)
        centered = ring_pts - centroid
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        normal = vh[-1] / max(float(np.linalg.norm(vh[-1])), 1e-9)
        if float(np.max(np.abs(centered @ normal))) > _F20_RING_TOL:
            continue
        for ri_idx in ring:
            for hi in nbrs[ri_idx]:
                if syms[hi] != 'H':
                    continue
                if float(abs(np.dot(P[hi] - centroid, normal))) > _F20_OOP_TOL:
                    cnt += 1
    return cnt


def _count_bond_distort(syms, P, bonds) -> int:
    """F3_bond proxy — heavy-heavy bonds whose length deviates outside the CSD
    band [-0.25, +0.085] of the element-pair ideal (matches the _violations
    distort band).  The decollapse moves non-frozen atoms to fix collapse but
    can stretch/compress other bonds out of band; this lets the acceptance gate
    reject such trades (the F3_bond regression seen vs iter23/iter25)."""
    cnt = 0
    for i, j in bonds:
        if syms[i] == "H" or syms[j] == "H":
            continue
        d = float(np.linalg.norm(P[i] - P[j]))
        ideal = _ideal_bond(syms[i], syms[j])
        if ideal <= 1e-6:
            continue
        dev = (d - ideal) / ideal
        if dev < -0.25 or dev > 0.085:
            cnt += 1
    return cnt


def _count_bad_angles(syms, P, heavy_nbr) -> int:
    """Optional angle proxy — count heavy-heavy-heavy bond angles below
    _ANGLE_MIN_DEG (physically implausible for any hybridization, exactly
    what collapse + repulsion-overshoot produce).  Monotone proxy for F3,
    not the exact metric (F3 needs RDKit hybridization)."""
    cnt = 0
    for c in range(len(syms)):
        if _is_metal(syms[c]) or syms[c] == 'H':
            continue
        nb = [k for k in heavy_nbr[c] if syms[k] != 'H']
        for a_i in range(len(nb)):
            va = P[nb[a_i]] - P[c]
            na = float(np.linalg.norm(va))
            if na < 1e-6:
                continue
            for b_i in range(a_i + 1, len(nb)):
                vb = P[nb[b_i]] - P[c]
                nbn = float(np.linalg.norm(vb))
                if nbn < 1e-6:
                    continue
                cosang = float(np.dot(va, vb) / (na * nbn))
                cosang = max(-1.0, min(1.0, cosang))
                if np.degrees(np.arccos(cosang)) < _ANGLE_MIN_DEG:
                    cnt += 1
    return cnt


def correct_xyz(mol, xyz: str) -> str:
    """Return a decollapsed copy of ``xyz`` (or the original if no gain).

    ``mol`` is accepted for signature parity but NOT used — connectivity is
    derived geometrically so the corrector is robust to mol/XYZ atom-order
    drift (the actual failure mode that made a mol-bond version detect 0
    collapses on a clearly-collapsed OKOKUU frame).
    """
    try:
        syms, P, _ = _parse(xyz)
    except Exception:
        return xyz
    n = len(syms)
    if n < 3:
        return xyz

    bonds = _geometric_bonds(syms, P)
    n_collapsed = _count_collapsed(syms, P, bonds)
    if n_collapsed == 0:
        return xyz  # bit-exact early return

    # Freeze metals + coordination sphere (atoms within M-D distance).
    frozen = set(i for i in range(n) if _is_metal(syms[i]))
    md = _md_snapshot(syms, P)
    for (mi, dj) in md:
        frozen.add(dj)
    md_snap = dict(md)

    # H → heavy parent (for rigid drag)
    h_parent: Dict[int, int] = {}
    heavy_nbr: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i, j in bonds:
        if syms[i] == 'H' and syms[j] != 'H':
            h_parent[i] = j
        elif syms[j] == 'H' and syms[i] != 'H':
            h_parent[j] = i
        if syms[i] != 'H' and syms[j] != 'H':
            heavy_nbr[i].append(j); heavy_nbr[j].append(i)

    # Iter-26: 1-2/1-3 exclusion + pre-relaxation baselines for the vdw/F20/
    # angle firewall (reject any relaxed frame that worsens a measured axis).
    excl = _neighbor_exclusion(n, bonds)
    anglegate = os.environ.get("DELFIN_BOND_DECOLLAPSE_ANGLEGATE", "0") == "1"
    clashes0 = _count_vdw_clashes(syms, P, excl)
    h_planar0 = _count_h_planar_viol(syms, P)
    bond0 = _count_bond_distort(syms, P, bonds)
    angle0 = _count_bad_angles(syms, P, heavy_nbr) if anglegate else 0

    work = P.copy()
    bonded_set = set()
    for i, j in bonds:
        bonded_set.add((min(i, j), max(i, j)))

    # ---- Perf hoist (byte-identical): precompute the per-pass-invariant
    # structure of both force loops ONCE.  The 250-pass relaxation re-derived
    # the same bond list, the same O(n^2) non-bonded pair set and the same
    # static floors on every pass; only the distances change.  Pulling the
    # invariants out and vectorizing the distance/force math (matmul-based norm
    # = bit-identical to the per-row np.linalg.norm via the same BLAS ddot
    # kernel; interleaved np.add.at scatter = bit-identical to the per-pair
    # `forces[a]+=f; forces[b]-=f` accumulation order) leaves the emitted
    # coordinates byte-for-byte unchanged while removing the Python per-pair
    # overhead and the millions of scalar np.linalg.norm calls.

    # Bonded springs: heavy-heavy bond index arrays + ideal target lengths.
    _bi: List[int] = []
    _bj: List[int] = []
    _btgt: List[float] = []
    # Aromatic ring bonds get the delocalised spring target when the seat flag is
    # on; empty set otherwise → every _ideal_bond call below is aromatic=False →
    # byte-identical spring targets.
    _arom_set = _aromatic_ring_bonds(syms, P, bonds) if _AROM_SEAT else frozenset()
    for i, j in bonds:
        if syms[i] == 'H' or syms[j] == 'H':
            continue
        _bi.append(i); _bj.append(j)
        _btgt.append(_ideal_bond(
            syms[i], syms[j], (min(i, j), max(i, j)) in _arom_set))
    bond_i = np.asarray(_bi, dtype=np.intp)
    bond_j = np.asarray(_bj, dtype=np.intp)
    bond_tgt = np.asarray(_btgt, dtype=float)
    n_bonds_hh = bond_i.shape[0]
    # interleaved scatter index for bonded force accumulation [i0,j0,i1,j1,...]
    bond_scatter = np.empty(2 * n_bonds_hh, dtype=np.intp)
    if n_bonds_hh:
        bond_scatter[0::2] = bond_i
        bond_scatter[1::2] = bond_j

    # Non-bonded repulsion: enumerate heavy non-bonded pairs in the exact
    # original iteration order (a outer 0..n, b inner a+1..n), precompute the
    # static per-pair floor (depends only on syms/excl, never on positions).
    _na: List[int] = []
    _nb: List[int] = []
    _nfloor: List[float] = []
    for a in range(n):
        if syms[a] == 'H':
            continue
        sa = syms[a]
        excl_a = excl[a]
        metal_a = _is_metal(sa)
        vdw_a = _VDW.get(sa, _VDW_DEFAULT)
        for b in range(a + 1, n):
            if syms[b] == 'H':
                continue
            if (a, b) in bonded_set:
                continue
            sb = syms[b]
            floor = _NONBOND_FLOOR * _ideal_bond(sa, sb)
            # Iter-26: for TRUE non-bonded pairs (not 1-3 neighbors) lift the
            # floor to the vdw-clash threshold so the relaxation stops
            # settling atoms into the clash band; 1-3 contacts keep the
            # gentler covalent floor (else rings blow apart).
            if b not in excl_a and not (metal_a and _is_metal(sb)):
                floor = max(floor, _VDW_FACTOR * (vdw_a + _VDW.get(sb, _VDW_DEFAULT)))
            _na.append(a); _nb.append(b); _nfloor.append(floor)
    nb_a = np.asarray(_na, dtype=np.intp)
    nb_b = np.asarray(_nb, dtype=np.intp)
    nb_floor = np.asarray(_nfloor, dtype=float)
    n_pairs = nb_a.shape[0]

    def _row_norm(deltas: np.ndarray) -> np.ndarray:
        """Per-row Euclidean norm, BIT-IDENTICAL to a loop of
        ``float(np.linalg.norm(row))``.  np.linalg.norm on a 1-D vector reduces
        to ``sqrt(row.dot(row))`` (BLAS ddot); batched matmul dispatches the
        same per-row ddot kernel, so the FMA/summation order matches exactly
        (verified array_equal across scales).  axis-reductions / (x**2).sum
        do NOT match (1-ULP drift on ~10% of elements)."""
        if deltas.shape[0] == 0:
            return np.zeros(0, dtype=float)
        sq = np.matmul(deltas[:, None, :], deltas[:, :, None]).reshape(-1)
        return np.sqrt(sq)

    for _pass in range(250):
        forces = np.zeros_like(work)
        moved = False
        # bonded springs → ideal length (vectorized, byte-identical)
        if n_bonds_hh:
            diff_b = work[bond_j] - work[bond_i]
            d_b = _row_norm(diff_b)
            # rare collapse fallback: a superimposed bonded pair (d<1e-6) is
            # re-seeded from a per-pass RNG, exactly as the scalar loop did.
            # This is order-sensitive (one RNG advanced per triggering bond in
            # bond order) so it is replayed scalar; it almost never fires.
            zero_mask = d_b < 1e-6
            if zero_mask.any():
                diff_b = diff_b.copy()
                for _k in np.nonzero(zero_mask)[0]:
                    rdiff = np.random.default_rng(_pass).standard_normal(3)
                    diff_b[_k] = rdiff
                    d_b[_k] = float(np.linalg.norm(rdiff))
            active = np.abs(d_b - bond_tgt) > 0.02
            if active.any():
                moved = True
                coef = 0.5 * (d_b - bond_tgt)
                fb = coef[:, None] * (diff_b / d_b[:, None])
                fb = np.where(active[:, None], fb, 0.0)
                vals = np.empty((2 * n_bonds_hh, 3), dtype=float)
                vals[0::2] = fb
                vals[1::2] = -fb
                np.add.at(forces, bond_scatter, vals)
        # non-bonded repulsion (resolve superimposed atoms) — vectorized
        if n_pairs:
            diff_p = work[nb_b] - work[nb_a]
            d_p = _row_norm(diff_p)
            active = (d_p < nb_floor) & (d_p > 1e-6)
            if active.any():
                moved = True
                coef = 0.5 * (d_p - nb_floor)
                # guard the division for inactive rows (d may be 0); masked out
                d_safe = np.where(active, d_p, 1.0)
                fp = coef[:, None] * (diff_p / d_safe[:, None])
                fp = np.where(active[:, None], fp, 0.0)
                idx = np.empty(2 * n_pairs, dtype=np.intp)
                idx[0::2] = nb_a
                idx[1::2] = nb_b
                vals = np.empty((2 * n_pairs, 3), dtype=float)
                vals[0::2] = fp
                vals[1::2] = -fp
                np.add.at(forces, idx, vals)
        if not moved:
            break
        # apply, but never move frozen atoms.  Clamp the per-pass displacement:
        # the vdw-aware repulsion floor can otherwise accumulate unbounded force
        # over the passes on a heavily-clashing frame and blow coordinates up to
        # inf -> NaN (Iter-26 regression).  A bounded step keeps the relaxation
        # stable and still converges.
        prev = work.copy()
        for i in range(n):
            if i in frozen or syms[i] == 'H':
                continue
            step = forces[i]
            sn = float(np.linalg.norm(step))
            if sn > _MAX_STEP:
                step = step * (_MAX_STEP / sn)
            work[i] = work[i] + step
        if not np.all(np.isfinite(work)):     # numeric blow-up -> abort cleanly
            work = prev
            break
    # drag H rigidly to follow their parent's displacement
    for h, par in h_parent.items():
        work[h] = work[h] + (work[par] - P[par])

    # acceptance (Iter-26 multi-axis firewall): keep the relaxed frame ONLY if
    # collapse strictly drops AND no measured axis worsens AND M-D intact.
    # Cheapest checks first so most rejections short-circuit before the SVD.
    if not np.all(np.isfinite(work)):     # never emit a non-finite frame
        return xyz
    if _count_collapsed(syms, work, bonds) >= n_collapsed:
        return xyz
    if _count_vdw_clashes(syms, work, excl) > clashes0:          # vdw (mandatory)
        return xyz
    if _count_h_planar_viol(syms, work) > h_planar0:            # F20 firewall
        return xyz
    if _count_bond_distort(syms, work, bonds) > bond0:          # F3_bond firewall
        return xyz
    if anglegate and _count_bad_angles(syms, work, heavy_nbr) > angle0:  # opt-in
        return xyz
    if _md_broken(md_snap, work):
        return xyz

    # rebuild xyz preserving original non-atom lines
    out_lines = xyz.splitlines()
    ai = 0
    res = []
    for ln in out_lines:
        p = ln.split()
        if len(p) == 4:
            try:
                float(p[1])
            except ValueError:
                res.append(ln); continue
            x, y, z = work[ai]
            res.append(f"{p[0]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            ai += 1
        else:
            res.append(ln)
    return "\n".join(res) + ("\n" if xyz.endswith("\n") else "")


def correct_results(mol, results):
    """Apply correct_xyz to every (xyz, label) frame in results."""
    if not results:
        return results
    out = []
    for item in results:
        try:
            xyz, label = item
        except Exception:
            out.append(item); continue
        try:
            out.append((correct_xyz(mol, xyz), label))
        except Exception:
            out.append((xyz, label))
    return out
