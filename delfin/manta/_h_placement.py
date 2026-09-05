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
def rotor_groups(syms: Sequence[str], parents: Dict[int, int],
                 adj: Sequence[Sequence[int]]
                 ) -> List[Tuple[int, int, List[int]]]:
    """(centre, heavy neighbour, H list) for centres with EXACTLY ONE heavy
    neighbour and at least one H.

    This condition is at the same time the stereo proof: all remaining
    substituents of the centre are hydrogens, hence constitutionally identical
    -- such a centre cannot be a stereocentre.
    """
    h_of: Dict[int, List[int]] = {}
    for h, p in parents.items():
        h_of.setdefault(p, []).append(h)
    out: List[Tuple[int, int, List[int]]] = []
    for c in sorted(h_of):
        heavy_nb = [j for j in adj[c] if _el.normalise(syms[j]) != "H"]
        if len(heavy_nb) != 1:
            continue
        out.append((c, heavy_nb[0], sorted(h_of[c])))
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


def _hp_stage_rotor(syms: List[str], P: np.ndarray, parents: Dict[int, int],
                    adj: Sequence[Sequence[int]]) -> int:
    """Rotor H rigidly about the centre-neighbour axis onto the best grid angle."""
    moved = 0
    for centre, nb, hs in rotor_groups(syms, parents, adj):
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
        best_ang = 0.0
        best_val = base
        for step in range(1, _ROTOR_STEPS):
            ang = 2.0 * np.pi * step / _ROTOR_STEPS
            P[hs] = _hp_rotate(orig, P[centre], axis, ang)
            if _group_metal() < min(1.0, base_m) - _EPS_GAIN:
                continue                   # metal axis protected independently
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
                                  "aborted_stereo": 0})
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
        if stats is not None:
            stats.update({"umbrella": n_a, "length": n_b, "rotor": n_c,
                          "moved": moved, "aborted_heavy_moved": 0,
                          "aborted_stereo": 0})
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

    print("== Stufe B: Laenge ==")
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
              "aborted_heavy_moved": 0, "aborted_stereo": 0}
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
          f"Abbruch-Stereotor {st_sum['aborted_stereo']}")
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
