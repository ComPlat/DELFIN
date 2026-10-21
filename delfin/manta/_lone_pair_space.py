"""delfin.manta._lone_pair_space -- the lone pairs of a coordinating sp3 donor take up space.

WHY (08.10.2026).  The FF-free seating puts the metal on "minus the sum of the unit vectors
from the donor to its own neighbours" (``assemble_orient._donor_and_lp``, the centroid rule in
``backbone_reembed._graft_ligand``).  That direction IS the lone pair only for a donor that has
exactly ONE lone pair (NR3, PR3, pyridine N).  A group-16 donor with two substituents (aqua,
alcohol, thioether) keeps TWO lone pairs; they sit on the two free vertices of its tetrahedron,
above and below the substituent plane -- not on the in-plane bisector.  Seating the metal on
the bisector puts it INTO the plane of its substituents: the donor is built trigonal planar
(angle sum exactly 360 deg).  ``_vsepr_reconstruct`` then sets every substituent to 109.47 deg
from that in-plane axis, keeping their in-plane azimuths, which opens H-O-H of a coordinated
water to 360 - 2 x 109.47 = 141.1 deg.  A donor with ONE substituent (hydroxo, thiolate) has
three free vertices; the centroid rule seats the metal antiperiplanar (M-O-H 180 deg).

THE RULE.  Lone pairs are pseudo-particles on the free vertices of the donor's ideal
tetrahedron, fitted onto its own substituents.  The metal takes ONE free vertex, the remaining
lone pair(s) keep the others: the substituents keep their tetrahedral distance from the metal
AND from the remaining lone pair, and the metal-donor axis runs along a lone pair.  The
pseudo-particles are pure geometry; they are never written to a frame.

SCOPE (graph-only, no element-pair table, no crystal statistics):
  * donor element of group 16 (O, S, Se, Te), not aromatic, not in a ring, only single bonds,
    and not sp/sp2 by RDKit's conjugation-aware hybridisation (phenolate, carboxylate, carbonyl
    O stay trigonal);
  * one or two substituents (hydrogens counted), all of them explicit atoms of the fragment;
  * an oxygen only when it carries hydrogen (aqua, hydroxo, alcohol): ether and alkoxide O are
    commonly trigonal at a metal and are left to the old rule.

Opt-in via ``DELFIN_FFFREE_LONE_PAIR_SPACE=1``; default OFF -> no call site reaches the new
code, byte-identical.
"""
from __future__ import annotations

import math
import os
from typing import List, Optional

import numpy as np

FLAG = "DELFIN_FFFREE_LONE_PAIR_SPACE"            # the seat pass on the finished frame (#525 II)
# The first version (#525 I) changed the donor direction INSIDE the seating.  lonepair_mini:
# the tilted H changed which chelate conformer won the clash selection (IJAVAT config 3 then
# collapsed onto the water O) -> configs, isomers and frames lost.  Kept dark as *_EARLY.
EARLY_FLAG = "DELFIN_FFFREE_LONE_PAIR_SPACE_EARLY"

_GROUP16 = frozenset(("O", "S", "Se", "Te"))
_TET = math.degrees(math.acos(-1.0 / 3.0))          # 109.4712 deg
_HALF_TET = math.radians(_TET / 2.0)                 # free vertices +- this from -bisector


def lone_pair_space_enabled() -> bool:
    return os.environ.get(EARLY_FLAG, "0") == "1"


def is_lone_pair_donor(atom) -> bool:
    """True for a coordinating group-16 donor whose lone pairs occupy tetrahedral vertices
    (see SCOPE in the module docstring).  ``atom`` is an RDKit atom of the FREE ligand."""
    try:
        sym = atom.GetSymbol()
        if sym not in _GROUP16 or atom.GetIsAromatic() or atom.IsInRing():
            return False
        for b in atom.GetBonds():
            if b.GetBondTypeAsDouble() != 1.0 or b.GetIsAromatic():
                return False
        if str(atom.GetHybridization()) in ("SP", "SP2"):
            return False
        n_h = atom.GetTotalNumHs(includeNeighbors=True)
        n_heavy = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)
        if n_h + n_heavy not in (1, 2):
            return False
        if sym == "O" and n_h == 0:
            return False
        return True
    except Exception:
        return False


def _perp(s: np.ndarray) -> Optional[np.ndarray]:
    """Deterministic unit vector perpendicular to ``s`` (same convention as DONOR_BEND)."""
    tmp = np.array([1.0, 0.0, 0.0]) if abs(float(s[0])) < 0.9 else np.array([0.0, 1.0, 0.0])
    p = tmp - s * float(np.dot(tmp, s))
    n = float(np.linalg.norm(p))
    return p / n if n > 1e-6 else None


def free_vertices(subs: List[np.ndarray]) -> Optional[List[np.ndarray]]:
    """Unit directions of the free vertices of the ideal tetrahedron fitted onto the donor's
    substituent vectors ``subs`` (one or two).  Deterministic order: the vertex the metal takes
    first, the remaining lone pair(s) after it.  None when undefined."""
    try:
        u = [np.asarray(v, float) / float(np.linalg.norm(v)) for v in subs]
    except Exception:
        return None
    if any(not np.all(np.isfinite(v)) for v in u):
        return None
    if len(u) == 2:
        b = u[0] + u[1]
        nb = float(np.linalg.norm(b))
        n = np.cross(u[0], u[1])
        nn = float(np.linalg.norm(n))
        if nb < 1e-3 or nn < 1e-3:
            return None                         # linear or collapsed substituents
        b = b / nb
        n = n / nn
        c, s = math.cos(_HALF_TET), math.sin(_HALF_TET)
        return [-b * c + n * s, -b * c - n * s]
    if len(u) == 1:
        s0 = u[0]
        p = _perp(s0)
        if p is None:
            return None
        q = np.cross(s0, p)
        ct, st = math.cos(math.radians(_TET)), math.sin(math.radians(_TET))
        out = []
        for k in range(3):
            phi = 2.0 * math.pi * k / 3.0
            out.append(ct * s0 + st * (math.cos(phi) * p + math.sin(phi) * q))
        return out
    return None


# ---- THE LENGTH OF A NEUTRAL O-H DONOR BOND (DELFIN_FFFREE_LONE_PAIR_MD, default OFF) ----------
# Every M-D distance of the FF-free seating is the Cordero covalent-radius sum (Zn-O 1.22 + 0.66
# = 1.88 A).  That is a COVALENT single bond; water or an alcohol donates a lone pair to a
# cation and sits at the IONIC contact distance -- Zn(II)-OH2 at CN 5 is ~2.0-2.1 A.  The rule
# used here is the classical one for that contact: Shannon effective ionic radius of the metal
# ion at its oxidation state and coordination number plus that of three-coordinate oxide
# (R. D. Shannon, Acta Cryst. A32 (1976) 751; derived from oxide/fluoride structures, no
# molecular-crystal statistics).  It applies ONLY to a neutral O donor that carries H and has two
# substituents (aqua, alcohol); unknown metal, oxidation state or CN -> None -> the caller keeps
# its covalent sum.  High-spin radii where the spin state matters.
LONE_PAIR_MD_FLAG = "DELFIN_FFFREE_LONE_PAIR_MD"
_R_O_CN3 = 1.36
_SHANNON = {
    ("Li", 1): {4: 0.59, 6: 0.76}, ("Na", 1): {4: 0.99, 5: 1.00, 6: 1.02, 7: 1.12, 8: 1.18},
    ("K", 1): {6: 1.38, 7: 1.46, 8: 1.51}, ("Mg", 2): {4: 0.57, 5: 0.66, 6: 0.72},
    ("Ca", 2): {6: 1.00, 7: 1.06, 8: 1.12}, ("Sr", 2): {6: 1.18, 8: 1.26},
    ("Ba", 2): {6: 1.35, 8: 1.42}, ("Al", 3): {4: 0.39, 5: 0.48, 6: 0.535},
    ("Sc", 3): {6: 0.745}, ("Ti", 4): {6: 0.605}, ("V", 3): {6: 0.64},
    ("V", 4): {5: 0.53, 6: 0.58}, ("Cr", 2): {6: 0.80}, ("Cr", 3): {6: 0.615},
    ("Mn", 2): {4: 0.66, 5: 0.75, 6: 0.83, 7: 0.90, 8: 0.96}, ("Mn", 3): {5: 0.58, 6: 0.645},
    ("Fe", 2): {4: 0.63, 6: 0.78}, ("Fe", 3): {4: 0.49, 5: 0.58, 6: 0.645},
    ("Co", 2): {4: 0.58, 5: 0.67, 6: 0.745}, ("Co", 3): {6: 0.545},
    ("Ni", 2): {4: 0.55, 5: 0.63, 6: 0.69}, ("Cu", 1): {4: 0.60},
    ("Cu", 2): {4: 0.57, 5: 0.65, 6: 0.73}, ("Zn", 2): {4: 0.60, 5: 0.68, 6: 0.74, 8: 0.90},
    ("Ga", 3): {6: 0.62}, ("Y", 3): {6: 0.90, 8: 1.019}, ("Zr", 4): {6: 0.72, 8: 0.84},
    ("Ru", 3): {6: 0.68}, ("Pd", 2): {4: 0.64}, ("Ag", 1): {4: 1.00, 6: 1.15},
    ("Cd", 2): {4: 0.78, 5: 0.87, 6: 0.95, 7: 1.03, 8: 1.10}, ("In", 3): {6: 0.80},
    ("Pt", 2): {4: 0.60}, ("Hg", 2): {4: 0.96, 6: 1.02}, ("Pb", 2): {6: 1.19},
    ("Bi", 3): {6: 1.03}, ("La", 3): {6: 1.032, 8: 1.16, 9: 1.216},
    ("Ce", 3): {6: 1.01, 8: 1.143}, ("Pr", 3): {6: 0.99, 8: 1.126},
    ("Nd", 3): {6: 0.983, 8: 1.109, 9: 1.163}, ("Sm", 3): {6: 0.958, 8: 1.079},
    ("Eu", 3): {6: 0.947, 8: 1.066}, ("Gd", 3): {6: 0.938, 8: 1.053},
    ("Tb", 3): {6: 0.923, 8: 1.04}, ("Dy", 3): {6: 0.912, 8: 1.027},
    ("Ho", 3): {6: 0.901, 8: 1.015}, ("Er", 3): {6: 0.89, 8: 1.004},
    ("Tm", 3): {6: 0.88, 8: 0.994}, ("Yb", 3): {6: 0.868, 8: 0.985},
    ("Lu", 3): {6: 0.861, 8: 0.977}, ("U", 6): {6: 0.73, 7: 0.81, 8: 0.86},
}
_CURRENT_OS = None


def lone_pair_md_enabled() -> bool:
    return os.environ.get(LONE_PAIR_MD_FLAG + "_EARLY", "0") == "1"


def set_current_os(ox) -> None:
    """Oxidation state of the complex being built (one per complex, like polyhedra's CN)."""
    global _CURRENT_OS
    try:
        _CURRENT_OS = int(ox) if ox is not None else None
    except Exception:
        _CURRENT_OS = None


def _shannon(metal: str, ox, cn) -> Optional[float]:
    tbl = _SHANNON.get((metal, ox))
    if not tbl or not cn:
        return None
    if cn in tbl:
        return tbl[cn]
    best = min(tbl, key=lambda c: (abs(c - int(cn)), c))
    return tbl[best]


def dative_oh_md(metal: str, atom, cn) -> Optional[float]:
    """Ionic contact distance metal-O for a neutral O donor carrying H with two substituents
    (aqua, alcohol); None when the rule does not apply or a number is unknown."""
    try:
        if atom is None or atom.GetSymbol() != "O" or not is_lone_pair_donor(atom):
            return None
        n_h = atom.GetTotalNumHs(includeNeighbors=True)
        n_heavy = sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)
        if n_h < 1 or n_h + n_heavy != 2:
            return None
        r = _shannon(metal, _CURRENT_OS, cn)
        if r is None:
            return None
        return float(r + _R_O_CN3)
    except Exception:
        return None


def metal_direction(mol, P: np.ndarray, di: int) -> Optional[np.ndarray]:
    """Unit vector from donor ``di`` towards the metal: a free tetrahedral vertex of the donor
    (its other free vertices stay lone pairs).  ``mol``/``P`` are the free ligand and its
    coordinates in the same atom order.  None when the rule does not apply -> the caller keeps
    its old direction."""
    try:
        atom = mol.GetAtomWithIdx(int(di))
    except Exception:
        return None
    if not is_lone_pair_donor(atom):
        return None
    if atom.GetNumImplicitHs() + atom.GetNumExplicitHs() != 0:
        return None                              # an H without coordinates: shape not fittable
    nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
    if len(nbrs) not in (1, 2):
        return None
    P = np.asarray(P, float)
    if P.ndim != 2 or P.shape[0] <= max(nbrs + [int(di)]):
        return None
    subs = [P[j] - P[int(di)] for j in nbrs]
    if any(float(np.linalg.norm(v)) < 1e-6 for v in subs):
        return None
    fv = free_vertices(subs)
    if not fv:
        return None
    a = fv[0]
    na = float(np.linalg.norm(a))
    if na < 1e-6 or not np.all(np.isfinite(a)):
        return None
    return a / na


# =====================================================================================
# THE SEAT ON THE FINISHED FRAME (#525 II, DELFIN_FFFREE_LONE_PAIR_SPACE / _LONE_PAIR_MD)
# =====================================================================================
# The direction rule above is right, but WHERE it acts decides what else changes.  Inside the
# seating it moved the H of a water before the chelate was placed; the conformer selection of
# the chelate reads clashes against everything already seated, picked another conformer and
# the self-gate rejected configurations (lonepair_mini: IJAVAT / LIWLAF lost the crystal isomer,
# 33 frames fewer).  And a pyramidal water whose SECOND lone pair points at a neighbour is a new
# defect (donor_lone_pair_clash) -- the pseudo-particle has to be repelled, which needs the
# neighbours to exist.
#
# So the lone pairs are seated when every atom of the frame stands: per monodentate group-16
# donor that was built planar (angle sum >= 355) or linear (M-O-H >= 165), the ligand is turned
# rigidly about its donor so that a free tetrahedral vertex carries the metal, and among the
# vertices x spins about the M-D axis the candidate is taken whose free lone pair(s) are EMPTY
# (no atom within 35 deg and 2.9 A, the eye's cone) and whose atoms keep the largest distance to
# everything else.  If no candidate has an empty lone pair and at least the old clearance, the
# donor stays as it was.  Heavy atoms of other ligands, the metal and the M-D distance are never
# moved (unless _LONE_PAIR_MD lengthens a neutral O-H donor bond to its ionic contact, again
# only when that keeps the clearance).
_SEAT_PLANAR_SUM = 355.0
_SEAT_LINEAR = 165.0
_LP_CONE_COS = 0.82
_LP_CONE_DIST = 2.9
_N_SPIN = 12
_MAX_LIG_ATOMS = 40
_VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
        "Cl": 1.75, "Br": 1.85, "I": 1.98, "Se": 1.90, "Te": 2.06, "B": 1.92, "Si": 2.10}


def seat_enabled() -> bool:
    return os.environ.get(FLAG, "0") == "1"


def seat_md_enabled() -> bool:
    return os.environ.get(LONE_PAIR_MD_FLAG, "0") == "1"


def _rot_a_to_b(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); c = float(np.dot(a, b)); s = float(np.linalg.norm(v))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        p = _perp(a)
        return 2.0 * np.outer(p, p) - np.eye(3)        # 180 deg about a perpendicular axis
    k = v / s
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    ang = math.atan2(s, c)
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)


def _rot_axis(axis: np.ndarray, ang: float) -> np.ndarray:
    k = axis / np.linalg.norm(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * (K @ K)


def _graph(syms, P):
    from delfin.manta import _elements as _el
    n = len(syms)
    met = [_el.is_metal(s) for s in syms]
    r = [_el.covalent_radius(s) for s in syms]
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    nb = [[] for _ in range(n)]
    mb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if met[i] and met[j]:
                continue
            if met[i] or met[j]:
                if D[i, j] < 1.30 * (r[i] + r[j]):
                    (mb[j] if met[i] else mb[i]).append(i if met[i] else j)
            elif D[i, j] < 1.20 * (r[i] + r[j]):
                nb[i].append(j); nb[j].append(i)
    return met, nb, mb, D


def _component(nb, d):
    seen = {d}; st = [d]
    while st:
        a = st.pop()
        for b in nb[a]:
            if b not in seen:
                seen.add(b); st.append(b)
    return seen


def _in_cycle(nb, d, comp):
    """True if two neighbours of d are connected without passing through d."""
    ns = nb[d]
    if len(ns) < 2:
        return False
    seen = {d, ns[0]}; st = [ns[0]]
    while st:
        a = st.pop()
        for b in nb[a]:
            if b == d or b in seen:
                continue
            if b in ns:
                return True
            seen.add(b); st.append(b)
    return False


_HB_ACCEPTORS = frozenset(("N", "O", "F"))
_HB_MIN = 1.5


def _clearance(P, syms, moved, others, polar_h=()):
    """Smallest contact ratio d / (vdW_a + vdW_b) between moved atoms and the rest, WITHOUT the
    polar-H...acceptor pairs (donor H against N/O/F): those are hydrogen bonds and are guarded on
    their own (_polar_min).  One mixed number let the seat trade a short O-H...O contact for a
    real H...H / H...C clash (LIWLAF, lonepair_mini2)."""
    best = 9.9
    for a in moved:
        ra = _VDW.get(syms[a], 1.8)
        dv = np.linalg.norm(P[others] - P[a], axis=1)
        for b, d in zip(others, dv):
            if a in polar_h and syms[b] in _HB_ACCEPTORS:
                continue
            q = float(d) / (ra + _VDW.get(syms[b], 1.8))
            if q < best:
                best = q
    return best


def _polar_min(P, syms, polar_h, others):
    """Shortest polar-H...acceptor distance (A)."""
    best = 9.9
    acc = [b for b in others if syms[b] in _HB_ACCEPTORS]
    if not acc:
        return best
    for a in polar_h:
        dv = np.linalg.norm(P[acc] - P[a], axis=1)
        best = min(best, float(dv.min()))
    return best


# A contact ratio d / (vdW_a + vdW_b) at or above this is "clear" for the seat (the acceptance
# below already used it).  #525 VI (EZIBOF, lonepair_mini4): the spin guard compared ONE number,
# the smallest ratio over all moved atoms.  A water sitting in a pocket between two chelate
# rings (H...C 1.68 A, ratio 0.58) took the orientation with the best minimum (0.60) -- and
# with it put its other H 1.74-1.87 A from three ring carbons, one of them 2.04 -> 1.87 A (0.70 ->
# 0.65), i.e. it traded one short contact for several.  A minimum cannot see a trade; every
# PAIR can: no contact below the contact floor may appear or get shorter.  The floor is the
# hand-set 0.70 x vdW-sum contact threshold that the inter-ligand clash term and the H-placement
# eye gate use as well (checked to stay silent on clean crystal structures); 0.75 as floor (first
# try) also blocked trades that end above it and gave back seat gains without a new contact.
_CONTACT_FLOOR = 0.70
_CLEAR = 0.75


def _contact_worsened(P0, P1, syms, moved, others, polar_h=()):
    """True if some pair (moved atom, other atom) is below _CONTACT_FLOOR in P1 and shorter than in P0.
    Polar-H...acceptor pairs are hydrogen bonds and guarded by _polar_min, as in _clearance."""
    oth = np.asarray(others, int)
    rb = np.array([_VDW.get(syms[b], 1.8) for b in others])
    acc = np.array([syms[b] in _HB_ACCEPTORS for b in others])
    for a in moved:
        s = _VDW.get(syms[a], 1.8) + rb
        q1 = np.linalg.norm(P1[oth] - P1[a], axis=1) / s
        q0 = np.linalg.norm(P0[oth] - P0[a], axis=1) / s
        bad = (q1 < _CONTACT_FLOOR) & (q1 < q0 - 1e-9)
        if a in polar_h:
            bad &= ~acc
        if bool(bad.any()):
            return True
    return False


def _lp_occupied(P, d, subs, m, others):
    """Number of atoms in the free lone-pair cone(s) of donor d (eye criterion)."""
    dp = P[d]
    u = [(P[j] - dp) / np.linalg.norm(P[j] - dp) for j in subs]
    um = (P[m] - dp) / np.linalg.norm(P[m] - dp)
    if len(u) == 2:
        lp = -(u[0] + u[1] + um)
        if float(np.linalg.norm(lp)) < 0.30:
            lps = []
        else:
            lps = [lp / np.linalg.norm(lp)]
    else:
        fv = free_vertices([u[0], um])
        lps = fv or []
    occ = 0
    for lp in lps:
        for x in others:
            v = P[x] - dp
            dist = float(np.linalg.norm(v))
            if 0.4 < dist <= _LP_CONE_DIST and float(np.dot(v / dist, lp)) > _LP_CONE_COS:
                occ += 1
    return occ


def _angle(a, b):
    c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def seat_frame(syms, P, ox=None):
    """Seat the lone pairs of every eligible monodentate group-16 donor of ONE frame.
    Returns (new_P, n_seated, n_kept); new_P is P itself when nothing changed."""
    P = np.asarray(P, float)
    n = len(syms)
    if n == 0 or not np.all(np.isfinite(P)):
        return P, 0, 0
    met, nb, mb, D = _graph(syms, P)
    metals = [i for i in range(n) if met[i]]
    if not metals:
        return P, 0, 0
    cn = {m: sum(1 for i in range(n) if m in mb[i]) for m in metals}
    Q = P.copy()
    n_seat = n_keep = 0
    for d in range(n):
        s = syms[d]
        if s not in _GROUP16 or len(mb[d]) != 1:
            continue
        m = mb[d][0]
        subs = nb[d]
        if len(subs) not in (1, 2):
            continue
        hs = [j for j in subs if syms[j] == "H"]
        heavy = [j for j in subs if syms[j] != "H"]
        if not hs:
            # a donor without H (thioether, thiolate): turning it turns a backbone whose rotamer
            # the crystal fixes (FIYDIB, lonepair_mini2: ccdc_backbone_lost) -- out of scope
            continue
        if s == "O":
            if any(len(nb[c]) != 4 for c in heavy):      # phenol / acid O-H: conjugated, left alone
                continue
        comp = _component(nb, d)
        if len(comp) > _MAX_LIG_ATOMS or any(mb[a] for a in comp if a != d):
            continue                                       # not a monodentate ligand
        if _in_cycle(nb, d, comp):
            continue
        # built planar / linear?
        um = Q[m] - Q[d]
        if len(subs) == 2:
            ssum = (_angle(um, Q[subs[0]] - Q[d]) + _angle(um, Q[subs[1]] - Q[d])
                    + _angle(Q[subs[0]] - Q[d], Q[subs[1]] - Q[d]))
            if ssum < _SEAT_PLANAR_SUM:
                continue
        else:
            if _angle(um, Q[subs[0]] - Q[d]) < _SEAT_LINEAR:
                continue
        moved = sorted(a for a in comp if a != d)
        others = [i for i in range(n) if i not in comp and not met[i]]
        if not others:
            continue
        clear0 = _clearance(Q, syms, moved + [d], others, set(hs))
        pol0 = _polar_min(Q, syms, hs, others)
        # target M-D length (ionic contact for a neutral O-H donor), else unchanged
        md0 = float(np.linalg.norm(um))
        md_t = md0
        if seat_md_enabled() and s == "O" and len(subs) == 2 and ox is not None:
            r = _shannon(syms[m], ox, cn.get(m))
            if r is not None and r + _R_O_CN3 > md0:
                md_t = r + _R_O_CN3
        uM = um / md0
        fv = free_vertices([Q[j] - Q[d] for j in subs])
        if not fv:
            continue
        best = None
        for k, v in enumerate(fv):
            R0 = _rot_a_to_b(v, uM)
            for t in range(_N_SPIN):
                R = _rot_axis(uM, 2.0 * math.pi * t / _N_SPIN) @ R0
                C = Q.copy()
                shift = -uM * (md_t - md0)                 # donor moves away from the metal
                for a in comp:
                    C[a] = (Q[a] - Q[d]) @ R.T + Q[d] + shift
                if len(subs) == 2 and hs:
                    # the H of the donor takes its own tetrahedral vertex: the in-plane seat
                    # opened H-O-H / C-O-H to 141 deg (_vsepr_reconstruct), a rigid turn keeps that
                    _anc = heavy[0] if heavy else subs[0]
                    _oth = [j for j in subs if j != _anc][0]
                    if syms[_oth] == "H":
                        _fv2 = free_vertices([uM, C[_anc] - C[d]])
                        if _fv2:
                            _cur = C[_oth] - C[d]
                            _L = float(np.linalg.norm(_cur))
                            _v = max(_fv2, key=lambda w: float(np.dot(w, _cur)))
                            C[_oth] = C[d] + _v * _L
                occ = _lp_occupied(C, d, subs, m, others)
                if occ:
                    continue
                if _polar_min(C, syms, hs, others) < min(pol0, _HB_MIN) - 1e-9:
                    continue                               # an H-bond may not get shorter than 1.5 A
                if _contact_worsened(Q, C, syms, moved + [d], others, set(hs)):
                    continue                               # no short contact traded for another
                cl = _clearance(C, syms, moved + [d], others, set(hs))
                if best is None or cl > best[0] + 1e-9:
                    best = (cl, C)
        if os.environ.get("DELFIN_FFFREE_LONE_PAIR_TRACE", "0") == "1":
            import sys as _sys
            _sys.stderr.write("[lp-seat] donor %s%d clear0=%.3f best=%s\n" % (s, d, clear0, None if best is None else round(best[0], 3)))
        if best is not None and best[0] >= min(clear0, _CLEAR) - 1e-9:
            Q = best[1]
            n_seat += 1
        else:
            n_keep += 1
    if n_seat == 0:
        return P, 0, n_keep
    return Q, n_seat, n_keep


def seat_results(results, smiles=None):
    """Apply seat_frame to every (xyz, label) of an FF-free result list.  Off -> results."""
    if not seat_enabled() or not results:
        return results
    from delfin.manta._fix_sp2n_planarize import _parse_xyz, _format_xyz
    ox = None
    if seat_md_enabled() and smiles:
        try:
            from rdkit import Chem
            from delfin.manta._oxidation_state import _ox_cut as oxidation_state  # #528: the cut, not the donor table
            from delfin.manta import _elements as _el
            om = Chem.MolFromSmiles(smiles)
            if om is None:
                om = Chem.MolFromSmiles(smiles, sanitize=False)
            if om is not None:
                mi = [a.GetIdx() for a in om.GetAtoms() if _el.is_metal(a.GetSymbol())]
                if len(mi) == 1:
                    r = oxidation_state(om, mi[0])
                    ox = r[0] if r else None
        except Exception:
            ox = None
    out = []
    for e in results:
        try:
            xyz, lbl = e[0], e[1]
            syms, pts, lines = _parse_xyz(xyz)
            newP, k, _ = seat_frame(syms, pts, ox)
            if k and np.all(np.isfinite(newP)):
                out.append((_format_xyz(lines, syms, newP), lbl) + tuple(e[2:]))
            else:
                out.append(e)
        except Exception:
            out.append(e)
    return out
