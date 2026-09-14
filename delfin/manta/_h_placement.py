"""_h_placement -- ONE hydrogen repair, built for the FF-free path.

THE MEASUREMENT THAT JUSTIFIES THIS MODULE (19.08.2026, 9 pools, 248 099 frames).
Four detector families together carry around 45 000 hard frames = 28 % of the
hardness mass, and none of them had an active mechanism:

    xh_hh_clash + core_hhclash   16 881   10.7 %   H...H below the crystal floor
    xh_stretch  + xh_orphan      14 165    9.0 %   X-H length / no parent
    methyl_broken                10 049    6.4 %   H-C-H more than 15 degrees off
    h_axis_H_proximal_via_donor   7 205    4.6 %   donor H points AT the metal

It is ONE physical error under six names: the heavy-atom skeleton is seated,
rotated, radially corrected -- and the hydrogen stays lying where the
embedding had put it.  The evidence is the guaranteed frame overlap
(pigeonhole argument, lower bound): methyl_broken & smiles_hyb-angle 68.4 %,
smiles_topology & xh_hh_clash 68.1 %, smiles_hyb-angle & xh_hh_clash 67.8 %.

WHY THIS MOTION CLASS AND NO OTHER.  The cost law is measured at four
points: ordering ~0 · isometry +0.98 pp · rigid rotation with a new
conformation +6.57 pp · re-embedding +11.9 pp hard frames.  All three stages
here touch EXCLUSIVELY hydrogen atoms; every heavy atom and every metal
stays put byte-exactly (checked at the end, not asserted).  That is the
cheapest class of all -- a re-embedding would be the most expensive and would
have to fail by the same law.

THE THREE STAGES, in this order:

  A  UMBRELLA (methyl_broken).  Delegates to ``_vsepr_repair.repair_terminal_groups``
     -- the existing, XYZ-only, never-refuted EX3 repairer -- and afterwards
     TAKES OVER exclusively the hydrogens.  Every heavy-atom
     displacement (CF3, SO3) is rejected.  Stage A is therefore NOT a new build,
     but the hooking-up of an existing mechanism.

  B  LENGTH (xh_stretch, xh_orphan).  Sets the X-H distance to the sum of the
     covalent radii, RADIALLY along the existing direction.  Not a single
     angle changes, hence no signed volume either: stage B mathematically
     cannot touch stereochemistry.

  C  ROTOR (xh_hh_clash, h_axis_H_proximal_via_donor).  Rotates the H group of a
     rotor (heavy atom with EXACTLY ONE heavy neighbour: CH3, NH2, OH, SH)
     rigidly about the X-Y axis onto the grid angle that maximises the closest
     contact.  Bond lengths and bond angles are preserved exactly.

STEREOCHEMISTRY IS PROVABLY INTACT HERE, not hopefully.  On 19.08. a plane
projection flattened 61 sp3 centres and flipped two signs; exactly that
cannot happen here:
  * Stage A fires only on centres with THREE terminal neighbours of the same
    element.  Three constitutionally identical substituents -> the centre is by
    definition not a stereocentre.
  * Stage B changes only radii, no directions -> the sign of the
    determinant is invariant.
  * Stage C fires only on centres with EXACTLY ONE heavy neighbour; all
    remaining substituents are hydrogens, hence identical -> no stereocentre.
It is measured anyway: ``stereo_signature`` / ``stereo_delta`` compare
the signed volume before and after, and the census prints it alongside.

THRESHOLDS.  The FIRING thresholds deliberately lie INSIDE the detector bands
of the eye (0.85/1.15 against the detector limits 0.75/1.25), so that a repair
never happens on the edge.  The TARGET values are physical -- the Cordero
covalent-radius sum from ``_elements.COV_R``, the same convention as in the eye
(``find_xh_integrity``), and the reference-free measurement over 296 696 X-H bonds
of the champion archive confirms it as the MEDIAN (aromatic C-H 1.080 against sum
1.07; N-H 1.034 against 1.02; O-H 0.990 against 0.97).  What is wrong is not the
rule, what is wrong is the TAIL -- and that is exactly what stage B corrects.  No
detector value is rebuilt and no threshold is aimed at.

NAMING.  All module-internal helpers carry the prefix ``_hp_``.  That is
not style but a lesson: the generic names (`_adjacency`, `_is_metal`,
`_parse_xyz`) exist a dozen times in the twin tree with SLIGHTLY
different semantics, and that is exactly where the radius census got stuck on
16.08.  An unambiguous name makes the copy visible instead of disguising it.

WIRING.  Default OFF (``DELFIN_FFFREE_H_PLACEMENT``), and this module has
ZERO call sites today -- the build is therefore byte-identical, because not a
single line in the build path is different.  The required call site is in
``smiles_converter.py`` immediately before the FF-free ``return _ff, None``, next to
``_apply_mirror_enum_if_enabled``; it is named in the report and deliberately
NOT built because of file territory.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from delfin.manta import _elements as _el

# --- Switch ---------------------------------------------------------------
# THE ONE READ SITE.  A second os.environ.get on the same name is exactly
# the construction on which the fire census failed on 14.08.
FLAG = "DELFIN_FFFREE_H_PLACEMENT"
# THE THIRD GUARD (2026-09-10, register #416): a whole-frame contact census
# before/after all three stages.  Default OFF -> byte-identical.  Same
# construction as the topology gate and the stereo gate: measured, not
# assumed, and the WHOLE frame rolls back.  One read site, like FLAG.
FLAG_CONTACT_GATE = "DELFIN_FFFREE_H_CONTACT_GATE"
# THE PARENT-REGAIN RULE OF THE TOPOLOGY GATE (2026-09-10, register #420).
# The topology gate reads bonds with the eye's perception, for which a
# stretched (1.60 A) or orphaned X-H is UNBONDED; repairing it to 1.07 A then
# reads as "topology changed" and the gate rolls the frame back.  Stage B was
# dead under it (hplace6k2: xh_orphan 121 -> 121).  With this rule an H may
# gain exactly ONE thing: its own parent.  Default OFF -> byte-identical.
FLAG_PARENT_REGAIN = "DELFIN_FFFREE_H_PARENT_REGAIN"
# THE EYE GATE (2026-09-14, register #444).  hplace6k3b: 0 capabilities lost,
# 5 gained, hard frames -3.2 points -- and still 5 systems with a WORSE ligand
# (broken_regressed) and 3 with a new donor lone-pair clash.  Read in the
# frames: only H moved; the core H...H clashes were gone, and in their place
# stood inter-ligand contacts the eye's `intclash_pair` (0.70 x vdW sum,
# ligand-wise) and `donor_lone_pair_clash` reject.  The contact gate compares
# COUNTS on the module's OWN floors, so a frame may trade a mild contact for a
# harder one and pass.  The eye gate compares the SEVERITY PROFILE under the
# eye's own criteria: no more violating pairs than before, no pair of the
# sorted profile worse than its counterpart, no more H in free lone-pair
# cones.  Same lesson as the pair-gate factor 0.70 (register #3xx): the
# builder's gate must be the eye's gate.  Default OFF -> byte-identical.
FLAG_EYE_GATE = "DELFIN_FFFREE_H_EYE_GATE"
# ROTOR ONLY AT sp3 CENTRES (2026-09-14, register #444).  The rotor turns a
# terminal H group about the centre-neighbour axis.  At a centre with fewer
# than four substituents that is not a group-16 hydroxyl/thiol, the H position
# is fixed by the plane (=N-H, =CH2, aryl-NH2, formyl H): a turn takes the H
# OUT of the conjugation plane (FIDSOD, WUFMUF conj_planar 0 -> 51/81) or
# pyramidalises an sp2 N-H (ODOGOF, 11 deg).  Rule: a centre may be turned iff
# it is tetrahedral (four substituents), a two-substituent group-16 atom, or a
# three-substituent centre that is already pyramidal (angle sum < 350 deg).
# Read from the geometry and the graph, no element list beyond group 16.
# Default OFF -> byte-identical.
FLAG_SP3_ONLY = "DELFIN_FFFREE_H_SP3_ONLY"

# --- Thresholds -----------------------------------------------------------
# Firing INSIDE the detector bands (eye: collision < 0.75 · stretch > 1.25).
_XH_FIRE_LO = 0.85
_XH_FIRE_HI = 1.15
# Crystal floor for H...H; identical to find_xh_integrity (1.50 A).
_HH_FLOOR = 1.50
# Real overlap of H against a heavy non-metal: 0.85 x vdW sum, but ONLY
# for genuinely non-bonded pairs (graph distance >= 4).  Without this condition
# every normal 1,3 contact reports a violation.
_H_HEAVY_FRAC = 0.85
_NONBONDED_MIN_HOPS = 4
# Parent detection: identical to metric_h_axis._HEAVY_PARENT_BOND.
_PARENT_FACTOR = 1.45
# Bond graph (non-metals only), factor as in _bond_criterion / _isolated_reseat.
_BOND_FACTOR = 1.30
# An H farther away than this is not "left lying" but lost;
# pulling it 1.5 A would no longer be an isometry.
_MAX_PULL = 2.20
# The eye's criterion for "donor H points at the metal" (metric_h_axis).
_PROX_H_M_MAX = 2.00
_PROX_DELTA = 0.05
# Hydride detection: NOT "some metal within reach", but "the NEAREST
# heavy atom IS a metal" -- exactly the rule by which find_xh_integrity files
# an M-H as a hydride and takes it out of the X-H axis.
#
# WHY THE FIRST VERSION WAS WRONG.  It excluded every H that ANY metal
# came closer to than 2.00 A -- and the finding ``h_axis_H_proximal_via_donor``
# IS by definition an H with d(M,H) < 2.00.  So the rule excluded exactly the
# class that stage C is built for; the module-internal census could never
# report it as anything but 0, and that looked like "does not occur".
#
# ⚠ A REMAINDER STAYS, and it belongs here and not in a footnote.  On the
# 2144 measured frames ``xh_n_m_hydride`` drops 21 -> 18, and the cause is
# NOT this rule but a METAL DEFINITION that diverges:
# ``_elements.METALS`` files Te and Sb as METALLOIDS, ``find_xh_integrity``
# files them as METALS.  In ALEKOS an H sits at 1.249 A from Te and 1.281 A
# from C -- 0.03 A decide who is the parent.  The module reads Te-H, sets
# it to the covalent sum 1.690 (textbook Te-H 1.65-1.70), and after that C at
# 1.150 A is the nearest heavy atom, hence no longer a hydride.  ``n_m_hydride_bad``
# stays 13 in BOTH arms, stretch/collision/orphan drop -- nothing was
# destroyed.  But "nearest heavy atom" is FRAGILE on squashed frames, and
# that is an open weakness, not a settled question.
# Rotor grid: 24 steps of 15 degrees.  Fixed, deterministic, no optimisation.
_ROTOR_STEPS = 24
# An improvement below this value is noise and is rejected.
_EPS_GAIN = 1e-3
# Detector thresholds for the module-internal census (see _hp_frame_census).
_DET_ORPHAN_MAX = 1.6
_DET_COLL_FRAC = 0.75
_DET_COLL_ABS = 0.90
_DET_STRETCH_FRAC = 1.25
_DET_STRETCH_ABS = 1.30
_DET_METHYL_DEG = 15.0
_TETRA_DEG = 109.471

_VDW: Dict[str, float] = {}
try:                                     # pragma: no cover - import path
    from delfin.manta._vdw_radii import VDW_RADII as _VDW  # type: ignore
except Exception:                        # pragma: no cover
    _VDW = {}
_VDW_DEFAULT = 1.80


# ---------------------------------------------------------------------------
# XYZ input/output -- header lines and atom order stay untouched.
# ---------------------------------------------------------------------------
def _hp_atom_line(line: str) -> bool:
    parts = line.split()
    if len(parts) < 4 or not parts[0][:1].isalpha():
        return False
    try:
        float(parts[1]); float(parts[2]); float(parts[3])
    except (ValueError, IndexError):
        return False
    return True


def _hp_read(xyz_str: str) -> Tuple[List[str], np.ndarray, List[str]]:
    """Symbols, coordinates and the ORIGINAL LINES.

    Nothing is reordered and nothing is thrown away: ``_hp_write`` rewrites
    exclusively the numbers of the atom lines; every other line (count,
    comment, blank line) passes through unchanged.
    """
    syms: List[str] = []
    pts: List[List[float]] = []
    lines = xyz_str.splitlines()
    for line in lines:
        if not _hp_atom_line(line):
            continue
        parts = line.split()
        syms.append(parts[0])
        pts.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if not pts:
        return syms, np.zeros((0, 3), dtype=float), lines
    return syms, np.asarray(pts, dtype=float), lines


def _hp_write(orig_lines: Sequence[str], syms: Sequence[str],
              P: np.ndarray) -> str:
    out: List[str] = []
    k = 0
    for line in orig_lines:
        if _hp_atom_line(line) and k < len(syms):
            x, y, z = P[k]
            out.append(f"{syms[k]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            k += 1
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if orig_lines else "")


# ---------------------------------------------------------------------------
# Element knowledge -- NO further radius table.  _elements is the canon.
# ---------------------------------------------------------------------------
def _hp_cov(sym: str) -> float:
    return _el.covalent_radius(sym)


def _hp_metal(sym: str) -> bool:
    """Delegates to the canon ``_elements.is_metal`` -- no predicate of its own."""
    return _el.is_metal(sym)


def _hp_vdw(sym: str) -> float:
    return float(_VDW.get(_el.normalise(sym), _VDW_DEFAULT))


def _hp_xh_target(parent_sym: str) -> float:
    """Target length of an X-H bond = sum of covalent radii (Cordero)."""
    return _hp_cov("H") + _hp_cov(parent_sym)


# ---------------------------------------------------------------------------
# Graph.
# ---------------------------------------------------------------------------
def _hp_h_contacts(syms: Sequence[str], P: np.ndarray) -> Dict[int, frozenset]:
    """For every H: the set of atoms it is bonded to under the covalent-radius rule,
    metals INCLUDED.  This is the topology the eye perceives for a hydrogen; the
    repair may change positions, never this set (register #357: when it did, every
    frame of the system lost `topo_correct_frame`)."""
    out: Dict[int, frozenset] = {}
    hs = set(_hp_hydrogens(syms))
    for h in hs:
        out[h] = frozenset()
    # ONE TRUTH, NOT TWO: the eye's topology detector perceives bonds with
    # `_bond_decollapse._geometric_bonds` (weddell find_smiles_topology_match:73).  The first
    # version of this gate used this module's own rule (_BOND_FACTOR * COV_R) and caught 4 of
    # the 5 hplace6k losers but not TAFROI: rotor-moved H land exactly in the margin where the
    # two rules disagree.  So the gate reads the same perception the verdict reads.
    try:
        from delfin.manta._bond_decollapse import _geometric_bonds
        pairs = _geometric_bonds(list(syms), np.asarray(P, float))
    except Exception:
        pairs = None
    if pairs is None:                      # fallback: this module's own rule
        n = len(syms)
        rh = _hp_cov("H")
        pairs = [(h, j) for h in hs for j in range(n) if j != h
                 and float(np.linalg.norm(P[h] - P[j]))
                 < _BOND_FACTOR * (rh + _hp_cov(_el.normalise(syms[j])))]
    acc: Dict[int, set] = {h: set() for h in hs}
    for i, j in pairs:
        if i in acc:
            acc[i].add(j)
        if j in acc:
            acc[j].add(i)
    for h in hs:
        out[h] = frozenset(acc[h])
    return out


def _hp_topo_change_is_parent_regain(before: Dict[int, frozenset],
                                     after: Dict[int, frozenset],
                                     parents: Dict[int, int]) -> bool:
    """True iff EVERY H whose bond set changed went from EMPTY to exactly
    {its own parent} -- the module's parent (``parents_of_h``: nearest heavy
    non-metal within 1.45 x target or 2.20 A).  Anything else -- a foreign
    heavy atom, a metal, a lost parent, a second bond -- is NOT a regain and
    the gate must fire as before.  Register #420: the rule is universal (no
    element, no system), and it is the ONLY transition the gate of #357 had
    forbidden that a repair legitimately needs."""
    for h, s_after in after.items():
        s_before = before.get(h, frozenset())
        if s_after == s_before:
            continue
        p = parents.get(h, -1)
        if p < 0 or s_before or s_after != frozenset({int(p)}):
            return False
    return True


def _hp_hydrogens(syms: Sequence[str]) -> List[int]:
    return [i for i, s in enumerate(syms) if _el.normalise(s) == "H"]


def _hp_metals(syms: Sequence[str]) -> List[int]:
    return [i for i, s in enumerate(syms) if _hp_metal(s)]


def _hp_graph(syms: Sequence[str], P: np.ndarray) -> List[List[int]]:
    """Bond graph via covalent radii; metals stay outside.

    Metals are deliberately NOT linked: the M-D bond is not a covalent bond in
    the sense of this distance rule, and taking it along would turn half the
    coordination sphere into 1,3 neighbours.  Unlike
    ``_isolated_reseat._adjacency`` (metals included) and
    ``_vsepr_repair._adjacency`` (15 elements, fallback 0.90) this graph reads
    the canon ``_elements.COV_R``.
    """
    n = len(syms)
    adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        si = _el.normalise(syms[i])
        if _hp_metal(si):
            continue
        for j in range(i + 1, n):
            sj = _el.normalise(syms[j])
            if _hp_metal(sj):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * (_hp_cov(si) + _hp_cov(sj)):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _hp_link_parents(adj: List[List[int]], parents: Dict[int, int]) -> None:
    """Enter the parent bond INTO the graph -- and why that is necessary.

    MEASURED in the self-test, and it is not a special case but the regular
    case of this module: a stretched C-H at 1.60 A lies ABOVE the bond limit
    of the graph (1.30 x 1.07 = 1.39), so the H counts as unbonded there.  Then
    its 1,3 neighbour is not excluded, the neighbour counts as a genuine
    non-bonded contact -- and the rollback rejects exactly the correction
    that would have fixed the error.  The repairer would have blocked itself
    on its own defect.  ``parents_of_h`` knows the bond (it draws the wider
    limit), so it is entered here after the fact.
    """
    for h, p in parents.items():
        if p not in adj[h]:
            adj[h].append(p)
        if h not in adj[p]:
            adj[p].append(h)


# ---------------------------------------------------------------------------
# THE HOLE IN THE MODULE'S OWN PROOF -- and the path that was NOT taken.
#
# Measured 19.08. on 2144 frames: ONE genuine stereo candidate among 1076 flipped
# its sign, with unchanged neighbourhood.  The proof "stage B changes only
# radii, so the sign cannot flip" holds only if the H has exactly ONE heavy
# atom as neighbour.  A squashed H between two heavy atoms is bonded to BOTH
# in the distance graph; radial to its parent is then NOT radial to the other
# centre, and that centre's determinant may flip.
#
# ⛔ THE OBVIOUS PATH WAS MEASURED AND REJECTED.  "Simply skip such H"
# (a filter on exactly one heavy graph neighbour) makes the module safe and
# at the same time almost ineffective -- because the ambiguous H ARE the
# collisions.  On the same 2144 frames, eye ``find_xh_integrity``:
#     without filter   xh_n_collision 80 -> 48,  hh_n_hard 225 -> 142
#     with    filter   xh_n_collision 80 -> 77,  hh_n_hard 225 -> 212
# So the filter costs around 90 % of the effect on the hardest stage in order
# to protect ONE frame out of 2144.  That is the wrong trade.
#
# ✅ INSTEAD, A GATE RATHER THAN A PROHIBITION: the signed volume is measured
# before and after the repair, and a frame in which a GENUINE stereo candidate
# flips or is flattened is rejected as a WHOLE.  Reach stays, the damage
# becomes impossible instead of unlikely -- the same construction as the
# rollbacks in _isolated_reseat and _fix_sp3_h_tetrahedrality.
# ---------------------------------------------------------------------------


def _hp_within(adj: Sequence[Sequence[int]], start: int, hops: int) -> Set[int]:
    """All atoms up to ``hops`` bonds away (including ``start``)."""
    seen = {int(start)}
    front = [int(start)]
    for _ in range(hops):
        nxt: List[int] = []
        for a in front:
            for b in adj[a]:
                if b not in seen:
                    seen.add(b)
                    nxt.append(b)
        front = nxt
        if not front:
            break
    return seen


def parents_of_h(syms: Sequence[str], P: np.ndarray) -> Dict[int, int]:
    """H index -> index of the heavy parent atom (non-metal).

    Three rules, all three taken over from the eye instead of newly invented:
      * If the NEAREST heavy atom is a METAL, the H is a hydride and gets
        no parent at all -- exactly the classification of
        ``find_xh_integrity``.  Hydrides are not our error.
      * Otherwise the parent is the nearest heavy non-metal.
      * It counts if d < 1.45 x target length (normal bond) OR
        d <= 2.20 A (H left lying that is still reachable).
    """
    out: Dict[int, int] = {}
    n = len(syms)
    if n == 0 or P.shape[0] != n:
        return out
    for h in _hp_hydrogens(syms):
        best_j = -1
        best_d = 1e9
        near_j = -1
        near_d = 1e9
        for j in range(n):
            if j == h:
                continue
            sj = _el.normalise(syms[j])
            if sj == "H":
                continue
            d = float(np.linalg.norm(P[h] - P[j]))
            if d < near_d:                 # nearest heavy atom, metals included
                near_d = d
                near_j = j
            if _hp_metal(sj):
                continue
            if d < best_d:                 # nearest heavy NON-metal
                best_d = d
                best_j = j
        if near_j >= 0 and _hp_metal(syms[near_j]):
            continue                       # hydride -> untouched
        if best_j < 0:
            continue
        if (best_d < _PARENT_FACTOR * _hp_xh_target(syms[best_j])
                or best_d <= _MAX_PULL):
            out[h] = best_j
    return out


# ---------------------------------------------------------------------------
# Clearance -- the ONE objective of all three stages.
# ---------------------------------------------------------------------------
def _hp_clearance(syms: Sequence[str], P: np.ndarray, h: int, parent: int,
                  skip: Set[int]) -> float:
    """Smallest NORMALISED distance of the H to everything that concerns it.

    Normalised means d / floor: >= 1 is clean, < 1 is a violation.
    Three floors, each from the detector it serves:
      H...H            1.50 A                  (find_xh_integrity)
      H...heavy atom   0.85 x vdW sum          (only at graph distance >= 4)
      H...metal        min(2.00, d(M,X)-0.05)  (metric_h_axis: the H may not be
                       closer to the metal than its own parent)
    A far-away metal does not constrain the objective by itself -- the
    quotient then becomes large without needing a special rule for it.
    """
    worst = 1e9
    for j in range(len(syms)):
        if j == h or j in skip:
            continue
        sj = _el.normalise(syms[j])
        d = float(np.linalg.norm(P[h] - P[j]))
        if _hp_metal(sj):
            d_mx = (float(np.linalg.norm(P[parent] - P[j]))
                    if parent >= 0 else 1e9)
            thr = min(_PROX_H_M_MAX, d_mx - _PROX_DELTA)
            if thr <= 1e-6:
                continue
        elif sj == "H":
            thr = _HH_FLOOR
        else:
            thr = _H_HEAVY_FRAC * (_hp_vdw("H") + _hp_vdw(sj))
        worst = min(worst, d / thr)
    return worst if worst < 1e9 else 1e9


def _hp_metal_clear(syms: Sequence[str], P: np.ndarray, h: int,
                    parent: int) -> float:
    """ONLY the metal axis, as a SEPARATE number -- and why it has to be one.

    MEASURED 19.08. on 2144 frames: with the metal axis INSIDE the minimum of
    ``_hp_clearance``, ``h_prox_donor`` rose from 48 to 51.  The reason is
    not a chemistry error but arithmetic: a minimum hides its
    sub-axes.  If the closest contact of an H stood at 0.60 (H...H), a
    correction was allowed to push the metal axis from 1.20 to 0.80 -- the
    minimum still rose from 0.60 to 0.80, the rollback saw an improvement, and
    the donor-H-at-metal finding was newly created.  The same shape as
    "a mean cannot see a dead class": the axis needs its own gate,
    not a place in a sum.
    """
    worst = 1e9
    for j in range(len(syms)):
        if j == h or not _hp_metal(syms[j]):
            continue
        d = float(np.linalg.norm(P[h] - P[j]))
        d_mx = (float(np.linalg.norm(P[parent] - P[j]))
                if parent >= 0 else 1e9)
        thr = min(_PROX_H_M_MAX, d_mx - _PROX_DELTA)
        if thr <= 1e-6:
            continue
        worst = min(worst, d / thr)
    return worst if worst < 1e9 else 1e9


def _hp_skip(syms: Sequence[str], adj: Sequence[Sequence[int]], h: int,
             parent: int, extra: Sequence[int] = ()) -> Set[int]:
    """What does not count as a contact for this H.

    Its own parent, the group members passed in -- and all
    HEAVY ATOMS within three bonds (non-bonded counts from four onwards).
    Hydrogens always remain in the test: 1,3 and 1,4 H...H are exactly the
    contacts the rotor is meant to rotate away, and their floor (1.50 A) lies
    well below any healthy geminal distance (~1.78 A).
    """
    out: Set[int] = set(int(x) for x in extra)
    if parent >= 0:
        out.add(int(parent))
    for j in _hp_within(adj, h, _NONBONDED_MIN_HOPS - 1):
        if _el.normalise(syms[j]) != "H":
            out.add(int(j))
    return out


def _hp_contact_census(syms: Sequence[str], P: np.ndarray,
                       adj: Sequence[Sequence[int]], parents: Dict[int, int]
                       ) -> Tuple[int, int, int]:
    """(H...H, H...heavy, H...metal): pairs below the module's OWN floors, whole frame.

    The same three floors and the same skip as ``_hp_clearance`` /
    ``_hp_metal_clear`` -- but COUNTED per pair instead of taken as a minimum,
    and over ALL hydrogens instead of the one being moved.  Both differences
    are the point (register #416):

    * Stage A has no rollback at all; stages B and C roll back on a MINIMUM,
      so a second contact may get worse as long as it stays above the new
      minimum (the arithmetic ``_hp_metal_clear`` names for the metal axis,
      here on every axis); and no stage can see what the stages do TOGETHER.
      Measured on hplace6k2: the repair created H...H contacts BELOW its own
      1.50 A floor on 11 of 24 regressed systems (+45 pairs), all of them
      genuine inter-ligand pairs, and RIMKON lost its last two valid frames
      that way (four blocker terms from one cause).
    * Three separate numbers, no sum: a sum hides a dead axis exactly as a
      minimum hides a sub-axis.

    ``adj`` / ``parents`` are passed in so that before and after are measured
    on the SAME graph (the one of the input geometry); each H...H pair counts
    once.  This is a census by the module's floors, not the eye -- families
    the module does not model (donor lone-pair clash, H-anomaly triangle)
    are invisible here and stay with the verdict.
    """
    n_hh = n_hv = n_hm = 0
    r_h = _hp_vdw("H")
    for h in _hp_hydrogens(syms):
        p = parents.get(h, -1)
        skip = _hp_skip(syms, adj, h, p)
        for j in range(len(syms)):
            if j == h or j in skip:
                continue
            sj = _el.normalise(syms[j])
            d = float(np.linalg.norm(P[h] - P[j]))
            if _hp_metal(sj):
                d_mx = (float(np.linalg.norm(P[p] - P[j]))
                        if p >= 0 else 1e9)
                thr = min(_PROX_H_M_MAX, d_mx - _PROX_DELTA)
                if thr > 1e-6 and d < thr:
                    n_hm += 1
            elif sj == "H":
                if j > h and d < _HH_FLOOR:
                    n_hh += 1
            elif d < _H_HEAVY_FRAC * (r_h + _hp_vdw(sj)):
                n_hv += 1
    return n_hh, n_hv, n_hm


# ---------------------------------------------------------------------------
# STAGE A -- the umbrella (methyl_broken).  Existing mechanism, made H-only.
# ---------------------------------------------------------------------------
def _hp_stage_umbrella(syms: List[str], P: np.ndarray,
                       tol_deg: float = _DET_METHYL_DEG) -> int:
    """Apply ``_vsepr_repair.repair_terminal_groups`` -- keep ONLY the H.

    The existing repairer sets a distorted terminal EX3 group onto its ideal
    VSEPR position and leaves centre and anchor standing.  So only group
    members are moved, and those are either all H (CH3, NH3)
    or all heavy (CF3, SO3).  By rejecting every non-H here,
    exactly the H share remains; the heavy-atom skeleton is untouched.

    ``tol_deg`` = 15 is the threshold of the detector ``methyl_broken`` itself;
    the repairer's module default (20) would leave part of the reported
    cases untouched.
    """
    try:
        from delfin.manta import _vsepr_repair as _vr
    except Exception:
        return 0
    try:
        block = "".join(
            f"{s:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n"
            for s, p in zip(syms, P)
        )
        fixed = _vr.repair_terminal_groups(block, tol=float(tol_deg))
        if fixed == block:
            return 0
        f_syms, f_P, _ = _hp_read(fixed)
        if len(f_syms) != len(syms) or f_P.shape[0] != P.shape[0]:
            return 0
        moved = 0
        for i, s in enumerate(syms):
            if _el.normalise(s) != "H":
                continue                   # heavy atoms stay where they are
            if float(np.linalg.norm(f_P[i] - P[i])) <= 1e-9:
                continue
            P[i] = f_P[i]
            moved += 1
        return moved
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# STAGE B -- the length (xh_stretch, xh_orphan).
# ---------------------------------------------------------------------------
def _hp_stage_length(syms: List[str], P: np.ndarray, parents: Dict[int, int],
                     adj: Sequence[Sequence[int]]) -> int:
    """Set X-H radially to the target length.  Direction stays, angles stay."""
    moved = 0
    for h in sorted(parents):
        p = parents[h]
        v = P[h] - P[p]
        d = float(np.linalg.norm(v))
        if d < 1e-6:
            continue                      # degenerate: there is no direction
        target = _hp_xh_target(syms[p])
        if target <= 1e-6:
            continue
        if _XH_FIRE_LO <= d / target <= _XH_FIRE_HI:
            continue
        skip = _hp_skip(syms, adj, h, p)
        before = _hp_clearance(syms, P, h, p, skip)
        before_m = _hp_metal_clear(syms, P, h, p)
        old = P[h].copy()
        P[h] = P[p] + v / d * target
        after = _hp_clearance(syms, P, h, p, skip)
        after_m = _hp_metal_clear(syms, P, h, p)
        # TWO rollbacks, not one.  The first protects the closest contact
        # of all, the second the metal axis as a quantity of its own --
        # otherwise the minimum hides it (see _hp_metal_clear).
        if (after < min(1.0, before) - _EPS_GAIN
                or after_m < min(1.0, before_m) - _EPS_GAIN):
            P[h] = old
            continue
        moved += 1
    return moved


# ---------------------------------------------------------------------------
# STAGE C -- the rotor (xh_hh_clash, h_axis_H_proximal_via_donor).
# ---------------------------------------------------------------------------
_SP3_ANGLE_SUM_MAX = 350.0   # deg; three substituents summing to >= this are planar (sp2)


def _hp_angle_sum(P: np.ndarray, centre: int, subs: Sequence[int]) -> float:
    """Sum of the three angles at ``centre`` over three substituents (deg); 0 if degenerate."""
    vecs = []
    for j in subs:
        v = P[j] - P[centre]
        nv = float(np.linalg.norm(v))
        if nv < 1e-6:
            return 0.0
        vecs.append(v / nv)
    total = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            cosang = max(-1.0, min(1.0, float(np.dot(vecs[i], vecs[j]))))
            total += float(np.degrees(np.arccos(cosang)))
    return total


def _hp_centre_turnable(syms: Sequence[str], P: np.ndarray, c: int,
                        nb: int, hs: Sequence[int],
                        adj: Optional[Sequence[Sequence[int]]] = None) -> bool:
    """May the H group at centre ``c`` be turned about the c-nb axis?

    Tetrahedral (four substituents): yes.  Two substituents at a group-16
    centre (hydroxyl, thiol): yes -- the lone pairs are the missing
    substituents.  Three substituents: only if the centre is pyramidal
    (sum of the three angles < 350 deg); a planar centre is sp2 and its H is
    fixed by the plane.  Two substituents elsewhere (=N-H, formyl, acetylenic
    H): no.

    THE GRAPH RULE (hplace6k4p, 14.09., register #451): a HETEROATOM centre
    (N, O, S, ...) with a lone pair next to an sp2 neighbour is conjugated
    (amide, aniline, amidine, enol, phenol) and its H are fixed IN the plane
    -- whatever the input geometry says.  The acceptance run showed why the
    geometric test is not enough: the assembler hands over amide NH2 slightly
    pyramidal (angle sum 335-345 deg), the rotor takes that as a licence and
    turns them 50-80 deg out of the conjugation plane (FIDSOD, WUFMUF).  The
    neighbour is sp2-like when it carries three heavy neighbours in a plane
    (angle sum >= 345 deg) or exactly two heavy neighbours at an angle above
    150 deg (sp centre).  A carbon centre (methyl on an arene) stays turnable.
    """
    n_sub = 1 + len(hs)
    sym_c = _el.normalise(syms[c])
    # THE METAL IS A SUBSTITUENT (hplace6k4p3, 14.09., register #451).  The
    # graph keeps metals outside (see _hp_graph), so a coordinated alcohol
    # O(H)(C)-M looks like a hydroxyl with ONE heavy neighbour and is turned
    # about C-O; the H lands in the M-O-C plane and the eye reads the sp3
    # donor as TRIGONAL-PLANAR (lone_pair_bv_sev: TAFROI, six Zn-O(H) donors
    # flattened in 6 of 9 frames; MEBRET Cu-O19).  Counted with the metal the
    # centre has two heavy substituents -- its H is no free rotor.  Universal:
    # any centre (also a metal-bound CH2) within a coordinate bond of a metal.
    for m in _hp_metals(syms):
        d_mc = float(np.linalg.norm(P[m] - P[c]))
        if d_mc <= (_hp_cov(_el.normalise(syms[m])) + _hp_cov(sym_c)) * _LP_MD_RATIO:
            return False
    if sym_c != "C" and adj is not None and n_sub <= 3:
        heavy_nb = [j for j in adj[nb] if _el.normalise(syms[j]) != "H"]
        if len(heavy_nb) == 3 and _hp_angle_sum(P, nb, heavy_nb) >= _SP3_ANGLE_SUM_MAX - 5.0:
            return False                   # conjugated lone pair: H stay in the plane
        if len(heavy_nb) == 2:
            v1 = P[heavy_nb[0]] - P[nb]
            v2 = P[heavy_nb[1]] - P[nb]
            n1, n2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
            if n1 > 1e-6 and n2 > 1e-6:
                cosang = max(-1.0, min(1.0, float(np.dot(v1, v2)) / (n1 * n2)))
                if float(np.degrees(np.arccos(cosang))) > 150.0:
                    return False           # sp neighbour (nitrile-like), H fixed
    if n_sub >= 4:
        return True
    if n_sub == 2:
        return sym_c in _LP_GROUP16
    return _hp_angle_sum(P, c, [nb] + list(hs)) < _SP3_ANGLE_SUM_MAX


def rotor_groups(syms: Sequence[str], parents: Dict[int, int],
                 adj: Sequence[Sequence[int]],
                 P: Optional[np.ndarray] = None, sp3_only: bool = False
                 ) -> List[Tuple[int, int, List[int]]]:
    """(centre, heavy neighbour, H list) for centres with EXACTLY ONE heavy
    neighbour and at least one H.

    This condition is at the same time the stereo proof: all remaining
    substituents of the centre are hydrogens, hence constitutionally identical
    -- such a centre cannot be a stereocentre.

    ``sp3_only`` (register #444): additionally require that the centre may be
    turned at all -- see :func:`_hp_centre_turnable`; needs ``P``.
    """
    h_of: Dict[int, List[int]] = {}
    for h, p in parents.items():
        h_of.setdefault(p, []).append(h)
    out: List[Tuple[int, int, List[int]]] = []
    for c in sorted(h_of):
        heavy_nb = [j for j in adj[c] if _el.normalise(syms[j]) != "H"]
        if len(heavy_nb) != 1:
            continue
        hs = sorted(h_of[c])
        if sp3_only and P is not None and not _hp_centre_turnable(
                syms, P, c, heavy_nb[0], hs, adj=adj):
            continue
        out.append((c, heavy_nb[0], hs))
    return out


def _hp_rotate(pts: np.ndarray, origin: np.ndarray, axis: np.ndarray,
               ang: float) -> np.ndarray:
    """Rodrigues -- rigid rotation; lengths and angles are preserved exactly."""
    k = axis / float(np.linalg.norm(axis))
    v = pts - origin
    c = float(np.cos(ang))
    s = float(np.sin(ang))
    return (origin + v * c + np.cross(k, v) * s
            + np.outer(v.dot(k), k) * (1.0 - c))


# ---------------------------------------------------------------------------
# The FREE lone-pair cone of group-16 donors -- vendored from the eye (2026-09-06, #360).
# ---------------------------------------------------------------------------
# hplace6k lost TAFROI although no H changed its bonding set: the rotor stage removed all
# seven H...H clashes of the system and turned the H into the free lone-pair cone of an sp3
# O/S donor -- `donor_lone_pair_clash` (hard) then fired on all 9 frames.  The rule and its
# numbers are the eye's (weddell/detectors/find_sp3_carbon_donor_geometry.py:296-347), written
# here as constants so that builder and verdict cannot drift apart (pattern
# `fffree/aromatic_bond_targets.py`).  A group-16 donor bound to a metal through ONE lone pair
# keeps its SECOND lone pair at the fourth tetrahedral vertex; nothing may be rotated into it.
_LP_GROUP16 = frozenset({"O", "S", "Se", "Te"})
_LP_CONE_COS = 0.82      # cos(~35 deg): inside this half-angle of the free-lp axis = "in the cone"
_LP_CLASH_MAX = 2.9      # A: a non-bonded atom this close along the lp axis clashes the lone pair
_LP_MD_RATIO = 1.20      # M-D contact counts as a coordinate bond up to this multiple of the radii sum


def _hp_free_lone_pairs(syms: Sequence[str], P: np.ndarray):
    """(donor, unit vector of the free lone pair, atoms exempt) for every group-16 donor that is
    bound to exactly one metal and carries exactly two non-metal substituents (sp3, one free lp)."""
    out = []
    n = len(syms)
    norm = [_el.normalise(s) for s in syms]
    metals = [i for i in range(n) if _hp_metal(norm[i])]
    for d in range(n):
        if norm[d] not in _LP_GROUP16:
            continue
        dp = P[d]
        m_bound = [m for m in metals
                   if float(np.linalg.norm(P[m] - dp)) <= (_hp_cov(norm[m]) + _hp_cov(norm[d])) * _LP_MD_RATIO]
        if len(m_bound) != 1:
            continue
        m = m_bound[0]
        subs = [k for k in range(n) if k not in (d, m) and not _hp_metal(norm[k])
                and float(np.linalg.norm(P[k] - dp)) < (_hp_cov(norm[d]) + _hp_cov(norm[k])) * _BOND_FACTOR]
        if len(subs) != 2:
            continue

        def _u(i):
            v = P[i] - dp
            return v / (float(np.linalg.norm(v)) + 1e-12)
        lp = -(_u(subs[0]) + _u(subs[1]) + _u(m))
        nn = float(np.linalg.norm(lp))
        if nn < 0.30:                      # near-planar donor -> no clear lp direction
            continue
        out.append((d, lp / nn, set(subs) | {m, d}))
    return out


def _hp_in_lone_pair_cone(P: np.ndarray, h: int, lps) -> bool:
    """Does H sit inside the free lone-pair cone of any donor in `lps`?  Same test as the eye."""
    for d, lp, exempt in lps:
        if h in exempt:
            continue
        xv = P[h] - P[d]
        dist = float(np.linalg.norm(xv))
        if dist < 0.4 or dist > _LP_CLASH_MAX:
            continue
        if float(np.dot(xv / (dist + 1e-12), lp)) > _LP_CONE_COS:
            return True
    return False


def _hp_lp_count(syms: Sequence[str], P: np.ndarray) -> int:
    """Number of hydrogens sitting in a free lone-pair cone -- recomputed from scratch, because an
    O-H / S-H rotor moves a SUBSTITUENT of the donor and with it the lone-pair axis, so other H
    can fall into the cone without moving (measured: TAFROI 2 -> 12 with the group-only check)."""
    lps = _hp_free_lone_pairs(syms, P)
    if not lps:
        return 0
    return sum(1 for h in _hp_hydrogens(syms) if _hp_in_lone_pair_cone(P, h, lps))


def _hp_stage_rotor(syms: List[str], P: np.ndarray, parents: Dict[int, int],
                    adj: Sequence[Sequence[int]]) -> int:
    """Rotor H rigidly about the centre-neighbour axis onto the best grid angle."""
    moved = 0
    for centre, nb, hs in rotor_groups(syms, parents, adj, P=P,
                                       sp3_only=h_sp3_only_enabled()):
        axis = P[centre] - P[nb]
        if float(np.linalg.norm(axis)) < 1e-6:
            continue
        skips = {h: _hp_skip(syms, adj, h, centre,
                             extra=list(hs) + [centre, nb]) for h in hs}

        def _group_clear() -> float:
            return min(_hp_clearance(syms, P, h, centre, skips[h]) for h in hs)

        def _group_metal() -> float:
            return min(_hp_metal_clear(syms, P, h, centre) for h in hs)

        base = _group_clear()
        base_m = _group_metal()
        if base >= 1.0:
            continue                       # nothing violated -> touch nothing
        orig = P[hs].copy()
        # Lone-pair cones: counted over ALL H and recomputed per angle, see _hp_lp_count.
        base_lp = _hp_lp_count(syms, P)
        best_ang = 0.0
        best_val = base
        for step in range(1, _ROTOR_STEPS):
            ang = 2.0 * np.pi * step / _ROTOR_STEPS
            P[hs] = _hp_rotate(orig, P[centre], axis, ang)
            if _group_metal() < min(1.0, base_m) - _EPS_GAIN:
                continue                   # metal axis protected independently
            if base_lp is not None and _hp_lp_count(syms, P) > base_lp:
                continue                   # no angle may put MORE H into free lone-pair cones (#360)
            val = _group_clear()
            if val > best_val + _EPS_GAIN:
                best_val = val
                best_ang = ang
        if best_ang == 0.0:
            P[hs] = orig
            continue
        P[hs] = _hp_rotate(orig, P[centre], axis, best_ang)
        moved += len(hs)
    return moved


# ---------------------------------------------------------------------------
# Stereochemistry -- signed volume, so that the proof has a measurement.
# ---------------------------------------------------------------------------
def stereo_signature(syms: Sequence[str], P: np.ndarray
                     ) -> List[Tuple[int, int, float, int, Tuple[int, ...]]]:
    """(centre, sign, |volume|, number of H, neighbours) per 4-neighbour centre.

    The three neighbours with the smallest indices span the determinant.
    Since this module never changes the atom order, before and after are
    directly comparable.

    ⚠ THE NUMBER OF H IS NOT DECORATION, IT IS HALF THE MEASUREMENT.  A centre
    with TWO or more hydrogens has two constitutionally identical
    substituents and is therefore NOT a stereocentre -- its determinant
    sign means nothing, and flipping a methyl over inevitably flips it.
    First measurement on 50 archive systems: 12 sign changes
    over 2204 centres, and the question "at which ones" is answered only by this field.
    Whoever reads only the total reads an artefact.
    """
    adj = _hp_graph(syms, P)
    # The same parent back-filling as in the repairer -- otherwise a centre
    # with a stretched X-H would have THREE neighbours before and FOUR after, and
    # the before/after comparison would concern different sets.
    _hp_link_parents(adj, parents_of_h(syms, P))
    out: List[Tuple[int, int, float]] = []
    for c in range(len(syms)):
        sc = _el.normalise(syms[c])
        if sc == "H" or _hp_metal(sc):
            continue
        nbrs = sorted(adj[c])
        if len(nbrs) != 4:
            continue
        a = P[nbrs[0]] - P[c]
        b = P[nbrs[1]] - P[c]
        d3 = P[nbrs[2]] - P[c]
        vol = float(np.dot(np.cross(a, b), d3))
        n_h = sum(1 for j in nbrs if _el.normalise(syms[j]) == "H")
        out.append((c, (1 if vol > 0 else (-1 if vol < 0 else 0)),
                    abs(vol), n_h, tuple(nbrs)))
    return out


def stereo_delta(before, after) -> Tuple[int, int, int, int, int, int]:
    """(common, flips, flat, candidates, candidate_flips, neighbour_changes).

    "Candidate" = centre with AT MOST ONE hydrogen.  Only there can a sign
    mean chemistry at all; from two H onwards two substituents are
    identical, and the determinant is a numbering question.
    "flattened" = |volume| falls below 10 % of its initial value -- the
    signature of the plane projection that destroyed 61 centres on 19.08.

    ⚠ A CENTRE WITH A CHANGED NEIGHBOURHOOD IS NOT COUNTED AS A FLIP,
    but reported separately.  Measured 19.08.: a single apparent flip
    among 1077 candidates, and it came from a squashed H changing its
    parent (Te 1.249 A against C 1.281 A -- the rule "nearest
    heavy atom" decides there by 0.03 A).  Thus before and after span
    DIFFERENT tripods; comparing the sign would be meaningless.
    The number does not disappear because of this, it just stands under the right name.
    """
    b = {t[0]: t for t in before}
    common = flips = flat = cand = cand_flips = nb_changed = 0
    for c, s, v, nh, nbrs in after:
        if c not in b:
            continue
        common += 1
        _, s0, v0, nh0, nbrs0 = b[c]
        if nbrs0 != nbrs:
            nb_changed += 1
            continue
        flipped = (s0 != 0 and s != 0 and s0 != s)
        if flipped:
            flips += 1
        if v0 > 1e-6 and v < 0.10 * v0:
            flat += 1
        if nh0 <= 1 and nh <= 1:
            cand += 1
            if flipped:
                cand_flips += 1
    return common, flips, flat, cand, cand_flips, nb_changed


# ---------------------------------------------------------------------------
# Public interface.
# ---------------------------------------------------------------------------
def h_placement_enabled() -> bool:
    """THE ONE read site of DELFIN_FFFREE_H_PLACEMENT (default 0)."""
    return os.environ.get(FLAG, "0") == "1"


def h_contact_gate_enabled() -> bool:
    """THE ONE read site of DELFIN_FFFREE_H_CONTACT_GATE (default 0)."""
    return os.environ.get(FLAG_CONTACT_GATE, "0") == "1"


def h_parent_regain_enabled() -> bool:
    """THE ONE read site of DELFIN_FFFREE_H_PARENT_REGAIN (default 0)."""
    return os.environ.get(FLAG_PARENT_REGAIN, "0") == "1"


def h_eye_gate_enabled() -> bool:
    """THE ONE read site of DELFIN_FFFREE_H_EYE_GATE (default 0)."""
    return os.environ.get(FLAG_EYE_GATE, "0") == "1"


def h_sp3_only_enabled() -> bool:
    """THE ONE read site of DELFIN_FFFREE_H_SP3_ONLY (default 0)."""
    return os.environ.get(FLAG_SP3_ONLY, "0") == "1"


# ---------------------------------------------------------------------------
# THE EYE GATE (register #444) -- the eye's own inter-ligand criterion.
# ---------------------------------------------------------------------------
# The factor the eye runs `intclash_pair` with (adapt_inter_ligand_clash,
# threshold 0.70: silent on 60/60 clean crystals, the real inter-ligand packing
# tail bottoms at 0.725).  H...H pairs keep the core floor 1.50 A (the eye
# filters intclash H-H above it, register #417/#422).  Same vdW table (Bondi,
# delfin.manta._vdw_radii == weddell find_inter_ligand_clash._VDW_RADII).
_EYE_PAIR_FACTOR = 0.70


def _hp_components(syms: Sequence[str], adj: Sequence[Sequence[int]]
                   ) -> List[int]:
    """Ligand id per atom: connected component of the metal-free bond graph
    (the eye's ligand split); metals get -1."""
    n = len(syms)
    comp = [-1] * n
    cid = 0
    for i in range(n):
        if comp[i] != -1 or _hp_metal(_el.normalise(syms[i])):
            continue
        comp[i] = cid
        stack = [i]
        while stack:
            a = stack.pop()
            for b in adj[a]:
                if comp[b] == -1 and not _hp_metal(_el.normalise(syms[b])):
                    comp[b] = cid
                    stack.append(b)
        cid += 1
    return comp


def _hp_eye_profile(syms: Sequence[str], P: np.ndarray, comp: Sequence[int],
                    parents: Dict[int, int]) -> List[float]:
    """Sorted (worst first) ratios d / floor of every H-involving pair of two
    DIFFERENT ligands that the eye would report: H...H below 1.50 A, H...heavy
    below 0.70 x vdW sum.  Metals are their own axis (contact census).

    An H reads the ligand of its parent (a stretched or orphaned H is
    otherwise its own component and would count against its own ligand).
    Index-free on purpose: a rigid turn of a CH3 permutes H labels without
    changing the geometry, and the profile must not change with it.
    """
    n = len(syms)
    r_h = _hp_vdw("H")
    out: List[float] = []
    hs = _hp_hydrogens(syms)
    hset = set(hs)

    def _lig(i: int) -> int:
        if i in hset:
            p = parents.get(i, -1)
            if p >= 0:
                return comp[p]
        return comp[i]

    for h in hs:
        lh = _lig(h)
        for j in range(n):
            if j == h:
                continue
            sj = _el.normalise(syms[j])
            if _hp_metal(sj):
                continue
            if sj == "H" and j < h:
                continue                   # each H...H pair once
            lj = _lig(j)
            if lj < 0 or lj == lh:
                continue
            d = float(np.linalg.norm(P[h] - P[j]))
            thr = _HH_FLOOR if sj == "H" else _EYE_PAIR_FACTOR * (r_h + _hp_vdw(sj))
            if d < thr:
                out.append(d / thr)
    out.sort()
    return out


def _hp_lp_pairs(syms: Sequence[str], P: np.ndarray) -> Set[Tuple[int, int]]:
    """{(occupant, donor)}: every NON-METAL atom sitting in the free lone-pair
    cone of a donor -- exactly the eye's test (find_sp3_carbon_donor_geometry
    ._lone_pair_occupied scans every atom, not only H).

    A SET, not a count (hplace6k4p, register #451): an O-H rotor that leaves
    cone A and enters cone B keeps the count, and the eye then reports the
    clash on the new pair.  And the occupant is NOT only H (hplace6k4p2):
    the free-lone-pair AXIS of a coordinated O-H donor is computed from its
    substituents, one of which is the H -- turning that H swings the axis,
    and a cis donor oxygen 2.3-2.4 A away that was outside the cone is now
    inside (TAFROI O16/O44, MEBQUI, MEBRET).  Only H moved, a heavy occupant
    appeared.  So every non-metal atom counts as an occupant, before and
    after; no (occupant, donor) pair may be new.
    """
    lps = _hp_free_lone_pairs(syms, P)
    out: Set[Tuple[int, int]] = set()
    if not lps:
        return out
    n = len(syms)
    for x in range(n):
        if _hp_metal(_el.normalise(syms[x])):
            continue
        for d, lp, exempt in lps:
            if x in exempt:
                continue
            if _hp_in_lone_pair_cone(P, x, [(d, lp, exempt)]):
                out.add((x, d))
    return out


def _hp_eye_worse(before: Sequence[float], after: Sequence[float],
                  lp_before: Set[Tuple[int, int]],
                  lp_after: Set[Tuple[int, int]]) -> bool:
    """Is the after-profile worse under the eye's criteria?  More violating
    pairs, or any pair of the sorted profile below its counterpart, or an H in
    a free lone-pair cone it was not in before."""
    if len(after) > len(before) or not lp_after <= lp_before:
        return True
    for a, b in zip(after, before):
        if a < b - _EPS_GAIN:
            return True
    return False


def repair_xyz(xyz: str, *, stats: Optional[dict] = None) -> str:
    """Repair ungated -- for self-test and measurement, NOT in the build path.

    Returns the INPUT OBJECT unchanged if no H was moved, and likewise
    if the final check finds a moved heavy atom: then the run is
    invalid, and an invalid run must not cost anything.
    """
    if not xyz:
        return xyz
    try:
        syms, P, lines = _hp_read(xyz)
        if P.shape[0] == 0 or not _hp_hydrogens(syms):
            return xyz
        P = P.astype(float).copy()
        frozen = P.copy()
        sig_before = stereo_signature(syms, frozen)
        contacts_before = _hp_h_contacts(syms, frozen)
        # Third guard: the census graph is built ONCE on the input geometry and
        # reused after the repair, so before and after are the same instrument.
        # Only computed when the gate is on -- OFF costs nothing and changes nothing.
        gate_on = h_contact_gate_enabled()
        eye_on = h_eye_gate_enabled()
        if gate_on or eye_on:
            adj0 = _hp_graph(syms, frozen)
            par0 = parents_of_h(syms, frozen)
            _hp_link_parents(adj0, par0)
        if gate_on:
            census_before = _hp_contact_census(syms, frozen, adj0, par0)
        if eye_on:
            comp0 = _hp_components(syms, adj0)
            eye_before = _hp_eye_profile(syms, frozen, comp0, par0)
            lp_before_n = _hp_lp_pairs(syms, frozen)
        n_a = _hp_stage_umbrella(syms, P)
        adj = _hp_graph(syms, P)
        parents = parents_of_h(syms, P)
        _hp_link_parents(adj, parents)
        n_b = _hp_stage_length(syms, P, parents, adj)
        n_c = _hp_stage_rotor(syms, P, parents, adj)
        # This module's claim, checked instead of asserted.
        for i, s in enumerate(syms):
            if _el.normalise(s) != "H" and float(
                    np.linalg.norm(P[i] - frozen[i])) > 1e-9:
                if stats is not None:
                    stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                                  "moved": 0, "aborted_heavy_moved": 1,
                                  "aborted_stereo": 0, "aborted_topology": 0})
                return xyz
        # THE TOPOLOGY GATE (2026-09-06, register #357).  hplace6k lost 5 capabilities
        # and 4 crystal isomers (GIYBOH JOCCUC MEBRET RIMKON TAFROI): no H detector
        # moved, but `topo_correct_frame` went true -> false on EVERY frame of those
        # systems.  A moved H had landed where the eye's bond perception reads a
        # different graph (bonded to a foreign heavy atom or metal, or detached from
        # its parent).  Same class as the cyclam pucker sibling (#353): a stage
        # without an assurance at its exit.  Rule from now on: each stage checks the
        # topology at its exit and rolls back.  Here: the set of atoms each H is
        # bonded to (covalent-radius rule, metals INCLUDED because an H on a metal is
        # exactly the case the eye punishes) must be IDENTICAL before and after.
        contacts_after = _hp_h_contacts(syms, P)
        if contacts_after != contacts_before and not (
                h_parent_regain_enabled()
                and _hp_topo_change_is_parent_regain(contacts_before,
                                                     contacts_after, parents)):
            # The parent-regain rule (#420, default OFF) lets exactly one
            # transition through: an H whose bond set was EMPTY (stretched or
            # orphaned under the eye's perception) and is now {its own parent}.
            # That is stage B's repair, and nothing else.
            if stats is not None:
                stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                              "moved": 0, "aborted_heavy_moved": 0,
                              "aborted_stereo": 0, "aborted_topology": 1})
            return xyz
        moved = n_a + n_b + n_c
        if moved == 0:
            if stats is not None:
                stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                              "moved": 0, "aborted_heavy_moved": 0,
                              "aborted_stereo": 0})
            return xyz
        # THE STEREO GATE.  No prohibition up front, a verdict afterwards: if in
        # this frame a GENUINE stereo candidate (at most one H at the centre)
        # flips its sign, or a centre is flattened, the WHOLE frame falls
        # back.  See the block above _hp_within: the cheaper path
        # (not touching such H in the first place) costs ~90 % of the effect on
        # the hardest stage and was therefore rejected.
        _, _, _flat, _, _cand_flips, _ = stereo_delta(
            sig_before, stereo_signature(syms, P))
        if _cand_flips or _flat:
            if stats is not None:
                stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                              "moved": 0, "aborted_heavy_moved": 0,
                              "aborted_stereo": 1})
            return xyz
        # THE CONTACT GATE (2026-09-10, register #416).  hplace6k2 improved every
        # verdict term over hplace6k and still did not land: 20 systems with a
        # frame that got WORSE, and RIMKON lost its last two valid frames (four
        # blocker terms, one cause: core H-H clash 0 -> 2).  The module had
        # created contacts BELOW ITS OWN 1.50 A floor -- because stage A has no
        # rollback, B and C roll back on a minimum, and nothing looked at the
        # frame as a whole.  Rule: the three contact counts of the module's own
        # floors may not rise on ANY axis; if one does, the whole frame falls
        # back.  Sized on 510 judged systems / 7171 frame pairs: 8.7 % of the
        # frames roll back, 38 % stay improved, 20 of 24 regressed systems are
        # covered, RIMKON completely (13 of 16 frames).
        if gate_on:
            census_after = _hp_contact_census(syms, P, adj0, par0)
            if any(a > b for a, b in zip(census_after, census_before)):
                if stats is not None:
                    stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                                  "moved": 0, "aborted_heavy_moved": 0,
                                  "aborted_stereo": 0, "aborted_contact": 1})
                return xyz
        # THE EYE GATE (2026-09-14, register #444): the severity profile under
        # the eye's own inter-ligand criterion and the lone-pair cone count may
        # not get worse -- measured on the input graph, like the census.
        if eye_on:
            eye_after = _hp_eye_profile(syms, P, comp0, par0)
            if _hp_eye_worse(eye_before, eye_after, lp_before_n,
                             _hp_lp_pairs(syms, P)):
                if stats is not None:
                    stats.update({"umbrella": 0, "length": 0, "rotor": 0,
                                  "moved": 0, "aborted_heavy_moved": 0,
                                  "aborted_stereo": 0, "aborted_contact": 0,
                                  "aborted_eye": 1})
                return xyz
        if stats is not None:
            stats.update({"umbrella": n_a, "length": n_b, "rotor": n_c,
                          "moved": moved, "aborted_heavy_moved": 0,
                          "aborted_stereo": 0, "aborted_contact": 0,
                          "aborted_eye": 0})
        return _hp_write(lines, syms, P)
    except Exception:
        return xyz


def apply_xyz(xyz: str) -> str:
    """Repair gated.  OFF -> the input object, unchanged."""
    if not h_placement_enabled():
        return xyz
    return repair_xyz(xyz)


def apply_to_results(results):
    """Frame list (xyz, label, ...) -> frame list.  Gated, fail-safe.

    Same shape as ``_h_vsepr_realism.correct_results`` and ``_me_bond_snap``, so
    that the call site at the FF-free exit looks like its neighbours.  The
    frame count never changes: this is a corrector, not an enumerator.
    """
    if not results or not h_placement_enabled():
        return results
    out = []
    for entry in results:
        try:
            out.append((repair_xyz(entry[0]),) + tuple(entry[1:]))
        except Exception:
            out.append(entry)
    return out


# ---------------------------------------------------------------------------
# Census inside the module -- direction indicator, NO substitute for the eye.
# ---------------------------------------------------------------------------
def _hp_frame_census(syms: Sequence[str], P: np.ndarray) -> Dict[str, int]:
    """The defect families by the PUBLISHED thresholds of the eye.

    ⚠ This is a REPLICA of the thresholds from ``find_xh_integrity``,
    ``full_verdict.adapt_methyl_quality`` and ``metric_h_axis`` -- not a second
    opinion of the eye.  It shows the DIRECTION of the effect inside the module
    itself; the verdict is passed by the real detector on the written frame pairs.
    """
    n = len(syms)
    res = {"xh_stretch": 0, "xh_collision": 0, "xh_orphan": 0,
           "hh_clash": 0, "methyl_broken": 0, "h_prox_donor": 0}
    hs = _hp_hydrogens(syms)
    metals = _hp_metals(syms)
    for h in hs:
        best_j, best_d = -1, 1e9
        for j in range(n):
            if j == h or _el.normalise(syms[j]) == "H":
                continue
            d = float(np.linalg.norm(P[h] - P[j]))
            if d < best_d:
                best_d, best_j = d, j
        if best_j < 0 or _hp_metal(syms[best_j]):
            continue                       # M-H is a finding of its own
        if best_d > _DET_ORPHAN_MAX:
            res["xh_orphan"] += 1
            continue
        target = _hp_xh_target(syms[best_j])
        if best_d < _DET_COLL_FRAC * target or best_d < _DET_COLL_ABS:
            res["xh_collision"] += 1
        elif best_d > _DET_STRETCH_FRAC * target or best_d > _DET_STRETCH_ABS:
            res["xh_stretch"] += 1
    for a in range(len(hs)):
        for b in range(a + 1, len(hs)):
            if float(np.linalg.norm(P[hs[a]] - P[hs[b]])) < _HH_FLOOR:
                res["hh_clash"] += 1
    parents = parents_of_h(syms, P)
    h_of: Dict[int, List[int]] = {}
    for h, p in parents.items():
        h_of.setdefault(p, []).append(h)
    for c, hlist in sorted(h_of.items()):
        if len(hlist) != 3:
            continue
        worst = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                u = P[hlist[i]] - P[c]
                v = P[hlist[j]] - P[c]
                nn = float(np.linalg.norm(u) * np.linalg.norm(v))
                if nn < 1e-9:
                    continue
                ang = float(np.degrees(np.arccos(
                    max(-1.0, min(1.0, float(np.dot(u, v)) / nn)))))
                worst = max(worst, abs(ang - _TETRA_DEG))
        if worst > _DET_METHYL_DEG:
            res["methyl_broken"] += 1
    for h, p in sorted(parents.items()):
        for m in metals:
            d_mh = float(np.linalg.norm(P[h] - P[m]))
            if d_mh >= _PROX_H_M_MAX:
                continue
            if d_mh <= float(np.linalg.norm(P[p] - P[m])) - _PROX_DELTA:
                res["h_prox_donor"] += 1
                break
    return res


def split_frames(text: str) -> List[str]:
    """Split a multi-frame XYZ (concatenated) into single frames."""
    lines = text.splitlines()
    frames: List[str] = []
    i = 0
    while i < len(lines):
        head = lines[i].split()
        if len(head) == 1 and head[0].isdigit():
            k = int(head[0])
            block = lines[i:i + 2 + k]
            if len(block) == 2 + k:
                frames.append("\n".join(block) + "\n")
            i += 2 + k
        else:
            i += 1
    return frames


# ---------------------------------------------------------------------------
# Self-test.
# ---------------------------------------------------------------------------
def _hp_selftest() -> int:
    fails = 0

    def _expect(name: str, ok: bool, detail: str = "") -> None:
        nonlocal fails
        if not ok:
            fails += 1
        mark = "ok  " if ok else "FAIL"
        print(f"  [{mark}] {name}" + (f" -- {detail}" if detail else ""))

    print("== Schalter ==")
    os.environ.pop(FLAG, None)
    _expect("Vorgabe AUS", h_placement_enabled() is False)
    stretched = ("5\ntest\n"
                 "C       0.000000     0.000000     0.000000\n"
                 "C       1.520000     0.000000     0.000000\n"
                 "H      -0.500000     1.520000     0.000000\n"
                 "H      -0.363000    -0.520000     0.901000\n"
                 "H      -0.363000    -0.520000    -0.901000\n")
    _expect("AUS ist byte-identisch", apply_xyz(stretched) == stretched)
    os.environ[FLAG] = "1"
    _expect("AN wird gelesen", h_placement_enabled() is True)

    print("== Stufe B unter dem Topologietor (register #420) ==")
    os.environ.pop(FLAG_PARENT_REGAIN, None)
    _expect("Elternteil-Erlass: Vorgabe AUS", h_parent_regain_enabled() is False)
    _st_topo: Dict[str, int] = {}
    _out_topo = repair_xyz(stretched, stats=_st_topo)
    _expect("ohne Erlass nimmt das Topologietor die Laengenreparatur zurueck (#420)",
            _out_topo is stretched and _st_topo.get("aborted_topology", 0) == 1,
            str(_st_topo))
    os.environ[FLAG_PARENT_REGAIN] = "1"
    _expect("Elternteil-Erlass: AN wird gelesen", h_parent_regain_enabled() is True)

    print("== Stufe B: Laenge (mit Erlass) ==")
    st: Dict[str, int] = {}
    out = repair_xyz(stretched, stats=st)
    _s0, P0, _ = _hp_read(stretched)
    _s1, P1, _ = _hp_read(out)
    d_before = float(np.linalg.norm(P0[2] - P0[0]))
    d_after = float(np.linalg.norm(P1[2] - P1[0]))
    _expect("gestrecktes C-H wird auf die Sollaenge gesetzt",
            d_before > 1.5 and abs(d_after - 1.07) < 1e-3,
            f"{d_before:.3f} -> {d_after:.3f}")
    _expect("Schweratome unbewegt",
            float(np.linalg.norm(P1[0] - P0[0])) < 1e-9
            and float(np.linalg.norm(P1[1] - P0[1])) < 1e-9)
    v0 = P0[2] - P0[0]
    v1 = P1[2] - P1[0]
    cosang = float(np.dot(v0, v1) / (np.linalg.norm(v0) * np.linalg.norm(v1)))
    _expect("Richtung unveraendert (kein Winkel bewegt sich)",
            abs(cosang - 1.0) < 1e-9)
    _expect("Zaehlung meldet Stufe B", st.get("length", 0) >= 1, str(st))

    print("== Stufe B: Waise ==")
    orphan = ("4\ntest\n"
              "N       0.000000     0.000000     0.000000\n"
              "C       1.470000     0.000000     0.000000\n"
              "H       0.000000     1.850000     0.000000\n"
              "H      -0.340000    -0.480000     0.830000\n")
    so, Po, _ = _hp_read(orphan)
    c0 = _hp_frame_census(so, Po)
    sr, Pr, _ = _hp_read(repair_xyz(orphan))
    c1 = _hp_frame_census(sr, Pr)
    _expect("verwaistes H bekommt seinen Elternteil zurueck",
            c0["xh_orphan"] == 1 and c1["xh_orphan"] == 0,
            f"{c0['xh_orphan']} -> {c1['xh_orphan']}")
    os.environ.pop(FLAG_PARENT_REGAIN, None)
    _expect("Elternteil-Erlass wieder AUS", h_parent_regain_enabled() is False)

    print("== Stufe A: Dach ==")
    # C2 hangs on C1 so that C1 is NOT a terminal group -- otherwise
    # repair_terminal_groups considers the centre anchorless and does not touch it.
    broken = ("6\ntest\n"
              "C       0.000000     0.000000     0.000000\n"
              "C       1.520000     0.000000     0.000000\n"
              "H      -0.400000     0.980000     0.000000\n"
              "H      -0.400000    -0.490000     0.849000\n"
              "H      -1.070000     0.000000    -0.100000\n"
              "C       2.040000     1.430000     0.000000\n")
    sb, Pb, _ = _hp_read(broken)
    cb0 = _hp_frame_census(sb, Pb)
    sb1, Pb1, _ = _hp_read(repair_xyz(broken))
    cb1 = _hp_frame_census(sb1, Pb1)
    _expect("gebrochenes Methyl wird repariert",
            cb0["methyl_broken"] == 1 and cb1["methyl_broken"] == 0,
            f"{cb0['methyl_broken']} -> {cb1['methyl_broken']}")
    _expect("Methylreparatur laesst die Schweratome stehen",
            float(np.linalg.norm(Pb1[0] - Pb[0])) < 1e-9
            and float(np.linalg.norm(Pb1[1] - Pb[1])) < 1e-9)

    print("== Stufe C: Rotor ==")
    # Methyl on C0-C1-C2, plus an N-H whose H sits exactly on a methyl H
    # (1.09 A, below the crystal floor 1.50).  All X-H stand exactly at their
    # target length, so that stages A and B demonstrably contribute nothing.
    import math as _m
    _rows = ["9", "test",
             f"C    {0.0:12.6f} {0.0:12.6f} {0.0:12.6f}",
             f"C    {1.53:12.6f} {0.0:12.6f} {0.0:12.6f}",
             f"C    {2.05:12.6f} {1.45:12.6f} {0.0:12.6f}"]
    _r = 1.07 * _m.sin(_m.radians(109.471))
    _x = -1.07 * _m.cos(_m.radians(109.471)) * -1.0
    for k in range(3):
        a = 2 * _m.pi * k / 3
        _rows.append(f"H    {_x:12.6f} {_r * _m.cos(a):12.6f} "
                     f"{_r * _m.sin(a):12.6f}")
    _rows.append(f"N    {_x:12.6f} {_r + 2.11:12.6f} {0.0:12.6f}")
    _rows.append(f"H    {_x:12.6f} {_r + 1.09:12.6f} {0.0:12.6f}")
    _rows.append(f"C    {_x + 1.47:12.6f} {_r + 2.11:12.6f} {0.0:12.6f}")
    clashing = "\n".join(_rows) + "\n"
    se, Pe, _ = _hp_read(clashing)
    ce0 = _hp_frame_census(se, Pe)
    _st_c: Dict[str, int] = {}
    se1, Pe1, _ = _hp_read(repair_xyz(clashing, stats=_st_c))
    ce1 = _hp_frame_census(se1, Pe1)
    _expect("H...H-Konflikt wird entdreht",
            ce0["hh_clash"] > 0 and ce1["hh_clash"] < ce0["hh_clash"],
            f"{ce0['hh_clash']} -> {ce1['hh_clash']}")
    _expect("nur der Rotor hat gearbeitet",
            _st_c.get("umbrella", -1) == 0 and _st_c.get("length", -1) == 0
            and _st_c.get("rotor", 0) > 0, str(_st_c))
    _expect("Rotor laesst ALLE Schweratome stehen",
            all(float(np.linalg.norm(Pe1[i] - Pe[i])) < 1e-9
                for i, s in enumerate(se) if s != "H"))
    dh0 = float(np.linalg.norm(Pe[3] - Pe[0]))
    dh1 = float(np.linalg.norm(Pe1[3] - Pe1[0]))
    _expect("Rotor erhaelt die C-H-Laenge exakt", abs(dh0 - dh1) < 1e-6,
            f"{dh0:.6f} -> {dh1:.6f}")

    print("== Stereochemie ==")
    chiral = ("5\ntest\n"
              "C       0.000000     0.000000     0.000000\n"
              "F       1.350000     0.000000     0.000000\n"
              "Cl     -0.560000     1.680000     0.000000\n"
              "Br     -0.640000    -0.900000     1.640000\n"
              "H      -0.300000    -0.420000    -0.760000\n")
    sc0, Pc0, _ = _hp_read(chiral)
    sig0 = stereo_signature(sc0, Pc0)
    sc1, Pc1, _ = _hp_read(repair_xyz(chiral))
    common, flips, flat, cand, cand_flips, _nbch = stereo_delta(
        sig0, stereo_signature(sc1, Pc1))
    _expect("Stereozentrum erkannt", cand >= 1,
            f"gemeinsam={common} Kandidaten={cand}")
    _expect("kein Vorzeichenwechsel", flips == 0 and cand_flips == 0,
            f"Wechsel={flips}/{cand_flips}")
    _expect("nicht plattgedrueckt", flat == 0, f"flach={flat}")
    # A methyl is NOT a stereocentre -- it must count as a non-candidate,
    # otherwise the census counts every rotated-away methyl as stereo damage.
    _mc = [t for t in stereo_signature(sb, Pb) if t[3] >= 2]
    _expect("Methylzentrum ist kein Stereokandidat", len(_mc) >= 1,
            f"Zentren mit >=2 H: {len(_mc)}")

    print("== Dritter Waechter: Kontakttor (register #416) ==")
    os.environ.pop(FLAG_CONTACT_GATE, None)
    _expect("Kontakttor: Vorgabe AUS", h_contact_gate_enabled() is False)
    # (1) The census counts PAIRS per axis by the module's own floors.  Two
    # ethane-like stubs far apart, then one H of each brought to 1.40 A (below
    # the 1.50 A floor) and to 1.60 A (above it).
    def _two_stubs(d_hh: float) -> str:
        return ("6\ntest\n"
                "C       0.000000     0.000000     0.000000\n"
                "C       1.520000     0.000000     0.000000\n"
                "H      -0.357000     1.009000     0.000000\n"
                f"H      -0.357000    {1.009 + d_hh:12.6f}     0.000000\n"
                f"C      -0.357000    {1.009 + d_hh + 1.070:12.6f}     0.000000\n"
                f"C      -0.357000    {1.009 + d_hh + 2.590:12.6f}     0.000000\n")

    def _census_of(frame: str) -> Tuple[int, int, int]:
        s_, P_, _ = _hp_read(frame)
        a_ = _hp_graph(s_, P_)
        p_ = parents_of_h(s_, P_)
        _hp_link_parents(a_, p_)
        return _hp_contact_census(s_, P_, a_, p_)
    _expect("Zensus zaehlt ein H...H-Paar unter 1,50 A genau einmal",
            _census_of(_two_stubs(1.40))[0] == 1, str(_census_of(_two_stubs(1.40))))
    _expect("Zensus ist still bei 1,60 A",
            _census_of(_two_stubs(1.60)) == (0, 0, 0), str(_census_of(_two_stubs(1.60))))
    # (2) THE MECHANISM, isolated from any stage's geometry: a stand-in stage of
    # stage A's shape (moves H, no rollback of its own) swings the methyl H of
    # the FIRST stub onto the H of the second, keeping its C-H length and its
    # bond -- so the topology gate stays silent -- while the rotor is a no-op
    # (the real rotor would swing it away again on a fixture this empty; on
    # RIMKON it could not, 13 of 16 frames).  Without the gate the frame goes
    # out with a new contact; with the gate the WHOLE frame falls back.
    # Fixture: an ethane stub (C0-C1, H2 on C0 pointing +y) and, below it, a
    # free-standing O-H (H3 at y=-2.459, O at -3.429, C at -4.859).  Nothing is
    # in contact.  The stand-in swings H2 to (-0.357, -1.009, 0): still 1.07 A
    # from C0, 2.13 A from C1 (no bond under any perception), and 1.45 A from
    # H3 -- below the 1.50 A floor.  O instead of C as H3's parent keeps the
    # H...heavy axis clean (0.85 x vdW: 2.31 A for O against 2.42 measured).
    stub = ("6\ntest\n"
            "C       0.000000     0.000000     0.000000\n"
            "C       1.520000     0.000000     0.000000\n"
            "H      -0.357000     1.009000     0.000000\n"
            "H      -0.357000    -2.459000     0.000000\n"
            "O      -0.357000    -3.429000     0.000000\n"
            "C      -0.357000    -4.859000     0.000000\n")
    _target = np.array([-0.357, -1.009, 0.0])
    _real_umbrella, _real_rotor = _hp_stage_umbrella, _hp_stage_rotor

    def _bad_umbrella(syms, P, tol_deg=_DET_METHYL_DEG):
        P[2] = _target.copy()
        return 1

    def _no_rotor(syms, P, parents, adj):
        return 0
    globals()["_hp_stage_umbrella"] = _bad_umbrella
    globals()["_hp_stage_rotor"] = _no_rotor
    try:
        _st_off: Dict[str, int] = {}
        _out_off = repair_xyz(stub, stats=_st_off)
        _expect("ohne Tor: die Stellvertreterstufe erzeugt einen neuen H...H-Kontakt",
                _out_off != stub and _census_of(stub)[0] == 0
                and _census_of(_out_off)[0] == 1
                and _st_off.get("aborted_topology", 0) == 0,
                f"{_census_of(stub)} -> {_census_of(_out_off)}  {_st_off}")
        os.environ[FLAG_CONTACT_GATE] = "1"
        _expect("Kontakttor: AN wird gelesen", h_contact_gate_enabled() is True)
        _st_on: Dict[str, int] = {}
        _out_on = repair_xyz(stub, stats=_st_on)
        _expect("mit Tor: der GANZE Frame faellt zurueck (Eingabeobjekt)",
                _out_on is stub and _st_on.get("aborted_contact", 0) == 1,
                str(_st_on))
        # The parent-regain rule (#420) must NOT widen the topology gate for
        # anything but a regained parent: swing H2 onto the foreign O (1.00 A,
        # a bond under any perception) -- with the rule ON the frame still
        # falls back, and it is the TOPOLOGY gate that says so.
        os.environ.pop(FLAG_CONTACT_GATE, None)
        os.environ[FLAG_PARENT_REGAIN] = "1"
        _s_st, _P_st, _ = _hp_read(stub)
        _foreign = _P_st[4] + np.array([0.0, 1.0, 0.0])   # 1.00 A above the O

        def _bad_umbrella_foreign(syms, P, tol_deg=_DET_METHYL_DEG):
            P[2] = _foreign.copy()
            return 1
        globals()["_hp_stage_umbrella"] = _bad_umbrella_foreign
        _st_f: Dict[str, int] = {}
        _out_f = repair_xyz(stub, stats=_st_f)
        _expect("Erlass AN: ein H an ein FREMDES Atom faellt weiter zurueck (Topologietor)",
                _out_f is stub and _st_f.get("aborted_topology", 0) == 1, str(_st_f))
        os.environ.pop(FLAG_PARENT_REGAIN, None)
    finally:
        globals()["_hp_stage_umbrella"] = _real_umbrella
        globals()["_hp_stage_rotor"] = _real_rotor
    # (3) A legitimate repair must survive the gate: the real stages on the
    # rotor fixture remove an H...H clash and create nothing on any axis.
    _c0 = _census_of(clashing)
    _st_r: Dict[str, int] = {}
    _out_r = repair_xyz(clashing, stats=_st_r)
    _c1 = _census_of(_out_r)
    _expect("mit Tor: eine echte Reparatur bleibt (Rotor entdreht H...H)",
            _out_r != clashing and _c1[0] < _c0[0]
            and all(x <= y for x, y in zip(_c1, _c0))
            and _st_r.get("aborted_contact", -1) == 0,
            f"{_c0} -> {_c1}  {_st_r}")
    _expect("Zaehler meldet 'gelaufen' getrennt von 'getroffen'",
            "aborted_contact" in _st_on and "aborted_contact" in _st_r)
    os.environ.pop(FLAG_CONTACT_GATE, None)
    _expect("Kontakttor wieder AUS", h_contact_gate_enabled() is False)

    print("== Rotor nur an sp3-Zentren (register #444) ==")
    os.environ.pop(FLAG_SP3_ONLY, None)
    _expect("sp3-Regel: Vorgabe AUS", h_sp3_only_enabled() is False)

    def _groups(block: str, sp3: bool):
        s_, P_, _ = _hp_read(block)
        a_ = _hp_graph(s_, P_)
        p_ = parents_of_h(s_, P_)
        _hp_link_parents(a_, p_)
        return rotor_groups(s_, p_, a_, P=P_, sp3_only=sp3)

    imine_nh = ("3\ntest\n"
                "C       0.000000     0.000000     0.000000\n"
                "N       1.300000     0.000000     0.000000\n"
                "H       0.795000     0.875000     0.000000\n")
    hydroxyl = ("3\ntest\n"
                "C       0.000000     0.000000     0.000000\n"
                "O       1.430000     0.000000     0.000000\n"
                "H       1.750000     0.900000     0.000000\n")
    methyl = ("5\ntest\n"
              "C       0.000000     0.000000     0.000000\n"
              "C       1.520000     0.000000     0.000000\n"
              "H       1.883000     1.028000     0.000000\n"
              "H       1.883000    -0.514000     0.890000\n"
              "H       1.883000    -0.514000    -0.890000\n")
    amine_nh2 = ("4\ntest\n"
                 "C       0.000000     0.000000     0.000000\n"
                 "N       1.470000     0.000000     0.000000\n"
                 "H       1.820000     0.950000     0.000000\n"
                 "H       1.820000    -0.480000     0.830000\n")
    planar_nh2 = ("4\ntest\n"
                  "C       0.000000     0.000000     0.000000\n"
                  "N       1.470000     0.000000     0.000000\n"
                  "H       1.970000     0.866000     0.000000\n"
                  "H       1.970000    -0.866000     0.000000\n")
    _expect("ohne Regel ist =N-H ein Rotor (die Stufe, die ODOGOF pyramidalisierte)",
            len(_groups(imine_nh, False)) == 1)
    _expect("mit Regel: =N-H (zwei Substituenten, kein Gruppe-16) wird nicht gedreht",
            len(_groups(imine_nh, True)) == 0)
    _expect("mit Regel: O-H (Gruppe 16, zwei Substituenten) bleibt Rotor",
            len(_groups(hydroxyl, True)) == 1)
    _expect("mit Regel: CH3 (tetraedrisch) bleibt Rotor",
            len(_groups(methyl, True)) == 1)
    _expect("mit Regel: pyramidales NH2 (Winkelsumme ~329) bleibt Rotor",
            len(_groups(amine_nh2, True)) == 1)
    _expect("mit Regel: planares NH2 (Winkelsumme 360) wird nicht gedreht",
            len(_groups(planar_nh2, True)) == 0)
    # THE METAL IS A SUBSTITUENT (#451, hplace6k4p3): the same hydroxyl
    # coordinated to Zn (O...Zn 1.95 A, pyramidal) is no rotor -- turning it
    # about C-O flattened the sp3 donor (TAFROI, MEBRET).  A far metal (4 A)
    # is no substituent, the hydroxyl stays a rotor.
    hydroxyl_on_metal = ("4\ntest\n"
                         "C       0.000000     0.000000     0.000000\n"
                         "O       1.430000     0.000000     0.000000\n"
                         "H       1.750000     0.900000     0.000000\n"
                         "Zn      1.900000    -0.600000     1.800000\n")
    hydroxyl_far_metal = ("4\ntest\n"
                          "C       0.000000     0.000000     0.000000\n"
                          "O       1.430000     0.000000     0.000000\n"
                          "H       1.750000     0.900000     0.000000\n"
                          "Zn      1.900000    -0.600000     4.000000\n")
    _expect("ohne Regel ist O-H am Metall ein Rotor (die Stufe, die TAFROI/MEBRET planarisierte)",
            len(_groups(hydroxyl_on_metal, False)) == 1)
    _expect("mit Regel: O-H am Metall (Zn 1,95 A) ist kein Rotor -- das Metall ist Substituent",
            len(_groups(hydroxyl_on_metal, True)) == 0)
    _expect("mit Regel: O-H mit fernem Metall (4 A) bleibt Rotor",
            len(_groups(hydroxyl_far_metal, True)) == 1)
    # THE GRAPH RULE (#451): an amide NH2 handed over slightly pyramidal (angle
    # sum ~335) next to a planar carbonyl carbon is conjugated -- not a rotor.
    # The same slightly pyramidal NH2 on an sp3 carbon (an amine) stays one,
    # and a methyl on an sp2 carbon (toluene) stays one too.
    amide_nh2 = ("6\ntest\n"
                 "C       0.000000     0.000000     0.000000\n"
                 "O       0.620000     1.060000     0.000000\n"
                 "C      -1.520000     0.000000     0.000000\n"
                 "N       0.700000    -1.220000     0.000000\n"
                 "H       0.250000    -2.100000     0.300000\n"
                 "H       1.700000    -1.250000    -0.300000\n")
    amine_nh2_sp3 = ("7\ntest\n"
                     "C       0.000000     0.000000     0.000000\n"
                     "H      -0.363000     1.028000     0.000000\n"
                     "H      -0.363000    -0.514000     0.890000\n"
                     "C      -0.507000    -0.717000    -1.242000\n"
                     "N       1.470000     0.000000     0.000000\n"
                     "H       1.820000     0.950000     0.000000\n"
                     "H       1.820000    -0.480000     0.830000\n")
    toluene_ch3 = ("6\ntest\n"
                   "C       0.000000     0.000000     0.000000\n"
                   "C       0.700000     1.212000     0.000000\n"
                   "C       0.700000    -1.212000     0.000000\n"
                   "C      -1.510000     0.000000     0.000000\n"
                   "H      -1.873000     1.028000     0.000000\n"
                   "H      -1.873000    -0.514000     0.890000\n")
    _expect("Graphregel: leicht pyramidales Amid-NH2 an planarem C=O-Kohlenstoff wird nicht gedreht",
            len(_groups(amide_nh2, True)) == 0 and len(_groups(amide_nh2, False)) == 1)
    _expect("Graphregel: dasselbe NH2 an sp3-Kohlenstoff (Amin) bleibt Rotor",
            len(_groups(amine_nh2_sp3, True)) == 1)
    _expect("Graphregel: CH2/CH3 an sp2-Kohlenstoff (Toluol) bleibt Rotor",
            len(_groups(toluene_ch3, True)) == 1)
    os.environ[FLAG_SP3_ONLY] = "1"
    _expect("sp3-Regel: AN wird gelesen", h_sp3_only_enabled() is True)
    os.environ.pop(FLAG_SP3_ONLY, None)

    print("== Augen-Tor: Schweregrad statt Zaehler (register #444) ==")
    os.environ.pop(FLAG_EYE_GATE, None)
    _expect("Augen-Tor: Vorgabe AUS", h_eye_gate_enabled() is False)
    # Methane next to a foreign O (own ligand): H1...O 1.81 A is an intclash
    # pair of the eye (floor 0.70 x 2.72 = 1.904) and a census pair of the
    # module (0.85 x 2.72 = 2.31).  A stand-in stage pushes H1 to 1.70 A:
    # the census COUNT stays 1 (the contact gate passes), the eye's severity
    # profile gets worse (the eye gate must roll back).
    trade = ("6\ntest\n"
             "C       0.000000     0.000000     0.000000\n"
             "H       1.090000     0.000000     0.000000\n"
             "H      -0.363000     1.028000     0.000000\n"
             "H      -0.363000    -0.514000     0.890000\n"
             "H      -0.363000    -0.514000    -0.890000\n"
             "O       2.900000     0.000000     0.000000\n")
    _s_t, _P_t, _ = _hp_read(trade)
    _a_t = _hp_graph(_s_t, _P_t)
    _p_t = parents_of_h(_s_t, _P_t)
    _hp_link_parents(_a_t, _p_t)
    _c_t = _hp_components(_s_t, _a_t)
    _prof0 = _hp_eye_profile(_s_t, _P_t, _c_t, _p_t)
    _expect("Profil vorher: genau ein Paar H1...O unter dem Augenboden",
            len(_prof0) == 1 and abs(_prof0[0] - 1.81 / 1.904) < 1e-3, str(_prof0))
    _P_perm = _P_t.copy()
    _P_perm[[1, 2, 3, 4]] = _P_t[[2, 3, 4, 1]]
    _expect("Profil ist indexfrei (permutierte H, gleiche Geometrie)",
            _hp_eye_profile(_s_t, _P_perm, _c_t, _p_t) == _prof0)

    _real_umbrella_e = globals()["_hp_stage_umbrella"]
    _real_rotor_e = globals()["_hp_stage_rotor"]
    _real_length_e = globals()["_hp_stage_length"]

    def _no_umbrella_e(syms_, P_, tol_deg=0.0):
        return 0

    def _no_length_e(syms_, P_, parents_, adj_):
        return 0

    def _worse_rotor_e(syms_, P_, parents_, adj_):
        P_[1] = np.array([1.20, 0.0, 0.0])      # H1...O 1.81 -> 1.70 A
        return 1

    def _cure_rotor_e(syms_, P_, parents_, adj_):
        P_[1] = np.array([0.545, 0.944, 0.0])   # H1...O 1.81 -> 2.54 A
        return 1

    globals()["_hp_stage_umbrella"] = _no_umbrella_e
    globals()["_hp_stage_length"] = _no_length_e
    globals()["_hp_stage_rotor"] = _worse_rotor_e
    try:
        os.environ[FLAG_CONTACT_GATE] = "1"
        _st_c: Dict[str, int] = {}
        _out_c = repair_xyz(trade, stats=_st_c)
        _expect("Kontakttor allein laesst den Tausch durch (Zaehler 1 -> 1)",
                _out_c != trade and _st_c.get("aborted_contact", -1) == 0, str(_st_c))
        os.environ[FLAG_EYE_GATE] = "1"
        _expect("Augen-Tor: AN wird gelesen", h_eye_gate_enabled() is True)
        _st_e: Dict[str, int] = {}
        _out_e = repair_xyz(trade, stats=_st_e)
        _expect("Augen-Tor nimmt das schlechtere Profil zurueck (0,951 -> 0,893)",
                _out_e is trade and _st_e.get("aborted_eye", 0) == 1, str(_st_e))
        os.environ.pop(FLAG_CONTACT_GATE, None)
        globals()["_hp_stage_rotor"] = _cure_rotor_e
        _st_g: Dict[str, int] = {}
        _out_g = repair_xyz(trade, stats=_st_g)
        _expect("Augen-Tor laesst eine echte Heilung durch (Paar verschwindet)",
                _out_g != trade and _st_g.get("aborted_eye", -1) == 0, str(_st_g))
    finally:
        globals()["_hp_stage_umbrella"] = _real_umbrella_e
        globals()["_hp_stage_length"] = _real_length_e
        globals()["_hp_stage_rotor"] = _real_rotor_e
        os.environ.pop(FLAG_EYE_GATE, None)
        os.environ.pop(FLAG_CONTACT_GATE, None)
    _expect("Augen-Tor wieder AUS", h_eye_gate_enabled() is False)

    print("== Idempotenz ==")
    once = repair_xyz(stretched)
    _expect("zweiter Lauf aendert nichts mehr", once == repair_xyz(once))

    print("== Metallnaehe ==")
    # nearest heavy atom IS the metal -> hydride, off limits.
    hydride = ("3\ntest\n"
               "Fe      0.000000     0.000000     0.000000\n"
               "H       1.600000     0.000000     0.000000\n"
               "C       3.500000     0.000000     0.000000\n")
    _expect("terminales Hydrid wird nicht angefasst",
            repair_xyz(hydride) == hydride)
    # Donor H that points AT the metal: nearest heavy atom is the donor,
    # hence NOT a hydride -- the repairer must be able to see this class.
    prox = ("5\ntest\n"
            "Fe      0.000000     0.000000     0.000000\n"
            "N       2.150000     0.000000     0.000000\n"
            "C       3.620000     0.000000     0.000000\n"
            "H       1.900000     0.480000     0.000000\n"
            "H       2.500000    -0.900000     0.400000\n")
    sp_, Pp_, _ = _hp_read(prox)
    _expect("Donor-H am Metall bekommt seinen Elternteil (kein Hydrid)",
            3 in parents_of_h(sp_, Pp_) and parents_of_h(sp_, Pp_)[3] == 1,
            str(parents_of_h(sp_, Pp_)))
    _expect("und der Zensus meldet ihn",
            _hp_frame_census(sp_, Pp_)["h_prox_donor"] == 1,
            str(_hp_frame_census(sp_, Pp_)))

    os.environ.pop(FLAG, None)
    _expect("Schalter wieder AUS", h_placement_enabled() is False)
    return fails


def _hp_run_census(paths: List[str], limit: int = 0) -> None:
    import glob as _glob
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(_glob.glob(os.path.join(p, "*.xyz"))))
        else:
            files.append(p)
    if limit > 0:
        files = files[:limit]
    keys = ["xh_stretch", "xh_collision", "xh_orphan", "hh_clash",
            "methyl_broken", "h_prox_donor"]
    tot_b = {k: 0 for k in keys}
    tot_a = {k: 0 for k in keys}
    fr_b = {k: 0 for k in keys}
    fr_a = {k: 0 for k in keys}
    n_frames = n_changed = n_ident = 0
    st_sum = {"umbrella": 0, "length": 0, "rotor": 0,
              "aborted_heavy_moved": 0, "aborted_stereo": 0,
              "aborted_topology": 0, "aborted_contact": 0,
              "aborted_eye": 0}
    s_common = s_flips = s_flat = s_cand = s_cand_flips = s_nbch = 0
    for fp in files:
        try:
            with open(fp, "r") as fh:
                text = fh.read()
        except Exception:
            continue
        for frame in split_frames(text):
            syms, P, _ = _hp_read(frame)
            if P.shape[0] == 0:
                continue
            n_frames += 1
            cb = _hp_frame_census(syms, P)
            for k in keys:
                tot_b[k] += cb[k]
                fr_b[k] += 1 if cb[k] else 0
            sig0 = stereo_signature(syms, P)
            st: Dict[str, int] = {}
            out = repair_xyz(frame, stats=st)
            for k in st_sum:
                st_sum[k] += st.get(k, 0)
            if out == frame:
                n_ident += 1
                for k in keys:
                    tot_a[k] += cb[k]
                    fr_a[k] += 1 if cb[k] else 0
                continue
            n_changed += 1
            syms2, P2, _ = _hp_read(out)
            ca = _hp_frame_census(syms2, P2)
            for k in keys:
                tot_a[k] += ca[k]
                fr_a[k] += 1 if ca[k] else 0
            c, f, fl, cd, cf, nbc = stereo_delta(
                sig0, stereo_signature(syms2, P2))
            s_common += c
            s_flips += f
            s_flat += fl
            s_cand += cd
            s_cand_flips += cf
            s_nbch += nbc
    print(f"Dateien {len(files)}  Frames {n_frames}  "
          f"veraendert {n_changed}  unveraendert {n_ident}")
    print(f"bewegte H: Dach {st_sum['umbrella']}  Laenge {st_sum['length']}  "
          f"Rotor {st_sum['rotor']}  "
          f"Abbruch-Schweratom {st_sum['aborted_heavy_moved']}  "
          f"Abbruch-Stereotor {st_sum['aborted_stereo']}  "
          f"Abbruch-Topologie {st_sum['aborted_topology']}  "
          f"Abbruch-Kontakttor {st_sum['aborted_contact']}  "
          f"Abbruch-Augentor {st_sum['aborted_eye']}")
    print(f"{'Befund':<16}{'Treffer vor':>13}{'nach':>9}"
          f"{'Frames vor':>13}{'nach':>9}")
    for k in keys:
        print(f"{k:<16}{tot_b[k]:>13}{tot_a[k]:>9}{fr_b[k]:>13}{fr_a[k]:>9}")
    print(f"Zentren mit 4 Nachbarn verglichen {s_common}  "
          f"Vorzeichenwechsel {s_flips}  plattgedrueckt {s_flat}")
    print(f"davon ECHTE Stereokandidaten (<=1 H) {s_cand}  "
          f"Vorzeichenwechsel {s_cand_flips}")
    print(f"Zentren mit GEAENDERTER Nachbarschaft (nicht vergleichbar) {s_nbch}")


def _hp_byteid(paths: List[str], limit: int = 0) -> int:
    """OFF proof on REAL frames: identical OBJECT, not just equal text.

    The switch is NOT set here.  What is checked is the gated
    interface (``apply_xyz`` / ``apply_to_results``) -- i.e. exactly what
    a future call site would use.  ``is`` comparison instead of ``==``:
    an equal string would already be good, the same object is better, because
    it proves that not even a reformat took place.
    """
    import glob as _glob
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(_glob.glob(os.path.join(p, "*.xyz"))))
        else:
            files.append(p)
    if limit > 0:
        files = files[:limit]
    n = bad = 0
    for fp in files:
        try:
            with open(fp, "r") as fh:
                text = fh.read()
        except Exception:
            continue
        frames = split_frames(text)
        for frame in frames:
            n += 1
            if apply_xyz(frame) is not frame:
                bad += 1
        got = apply_to_results([(f, "x") for f in frames])
        for k, f in enumerate(frames):
            if got[k][0] is not f:
                bad += 1
    print(f"Schalter AUS: {n} Frames, ueber apply_xyz UND apply_to_results, "
          f"nicht-identische Rueckgaben: {bad}")
    print("BYTE-IDENTISCH" if bad == 0 else "VERLETZT")
    return bad


def _hp_write_pairs(paths: List[str], out_before: str, out_after: str,
                    limit: int = 0) -> None:
    """Before/after as two directories -- so that the REAL detector judges."""
    import glob as _glob
    files: List[str] = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(_glob.glob(os.path.join(p, "*.xyz"))))
        else:
            files.append(p)
    if limit > 0:
        files = files[:limit]
    os.makedirs(out_before, exist_ok=True)
    os.makedirs(out_after, exist_ok=True)
    n = 0
    for fp in files:
        try:
            with open(fp, "r") as fh:
                text = fh.read()
        except Exception:
            continue
        stem = os.path.splitext(os.path.basename(fp))[0]
        for i, frame in enumerate(split_frames(text)):
            syms, P, _ = _hp_read(frame)
            if P.shape[0] == 0:
                continue
            with open(os.path.join(out_before, f"{stem}_{i:04d}.xyz"), "w") as fh:
                fh.write(frame)
            with open(os.path.join(out_after, f"{stem}_{i:04d}.xyz"), "w") as fh:
                fh.write(repair_xyz(frame))
            n += 1
    print(f"{n} Framepaare geschrieben: {out_before} / {out_after}")


if __name__ == "__main__":
    import sys

    _argv = sys.argv[1:]

    def _hp_pop_limit(rest: List[str]) -> Tuple[List[str], int]:
        if "--limit" in rest:
            k = rest.index("--limit")
            return rest[:k] + rest[k + 2:], int(rest[k + 1])
        return rest, 0

    if _argv and _argv[0] == "--census":
        _rest, _lim = _hp_pop_limit(_argv[1:])
        os.environ[FLAG] = "1"
        _hp_run_census(_rest, limit=_lim)
        raise SystemExit(0)
    if _argv and _argv[0] == "--byteid":
        _rest, _lim = _hp_pop_limit(_argv[1:])
        os.environ.pop(FLAG, None)
        raise SystemExit(1 if _hp_byteid(_rest, limit=_lim) else 0)
    if _argv and _argv[0] == "--pairs":
        _rest, _lim = _hp_pop_limit(_argv[1:])
        os.environ[FLAG] = "1"
        _hp_write_pairs(_rest[:-2], _rest[-2], _rest[-1], limit=_lim)
        raise SystemExit(0)
    print("== _h_placement Selbsttest ==")
    _n_fail = _hp_selftest()
    print("== " + ("ALLE BESTANDEN" if _n_fail == 0
                   else f"{_n_fail} FEHLER") + " ==")
    raise SystemExit(1 if _n_fail else 0)
