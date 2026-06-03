"""delfin.fffree.amide_vsepr_template — Build-time amide planarity template.

Construction-side fix #1 (User mandate 2026-06-03): F24_amide voll-pool +99.5%
remains after GRIP polish because the construction step leaves amide nitrogen
atoms pyramidal.  Polish (sp2-flatten + post_grip_corrector) only re-projects
sp2 candidates whose mean angle already falls into [115, 130]; an amide-N at
~104 deg slips out of the geometric window but it IS an amide.

This module implements a GRAPH-BASED amide motif detector and a strict
planarity projector:

    1. Find every C(=O) carbonyl unit (C with neighbour O at d < 1.30 A).
    2. Find every N bonded to that C with 3 total non-H neighbours -> amide-N.
       Also captures imide, carbamate, urea (any C-N where C is a carbonyl).
    3. For each amide-N: project the N atom onto the plane defined by the
       (carbonyl-C, alpha-substituents) plane via SVD of the three N-neighbour
       displacement vectors.  After projection the 3rd-substituent OOP becomes
       0 by construction.
    4. M-D INVARIANT: donor atoms (and the metal) are NEVER moved.  If the
       amide-N IS a donor it is skipped (rare but possible -- an
       N(C=O)-coord deprotonated amidate to a metal).

Universal, deterministic, geometry-only.  Env-gated default-OFF byte-identical.

API:
    enforce_amide_planarity(syms, P, metal_idx=-1, donor_idxs=()) -> (P_new, n_fixed)

Env flag:
    DELFIN_FFFREE_AMIDE_VSEPR=1    (per-fix activation)
    DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1   (master activation)
"""
from __future__ import annotations

import os
from typing import List, Sequence, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Detector calibration:
#   Carbonyl C=O bond cutoff (covers C=O ~1.21 A in aldehyde/ketone/amide).
_CARBONYL_CO_MAX = 1.30
#   Bond-perception cutoff (× covalent-sum) for the amide-N neighbour graph.
_BOND_FACTOR = 1.30
#   OOP threshold above which we project (matches F24 tol 0.20 A; we use a
#   slightly more sensitive 0.10 A to leave headroom for post-polish drift).
_OOP_THRESH = 0.10
#   M-D invariant verification tolerance (mirror post_grip_corrector).
_MD_INVARIANT_TOL = 0.05


def _flag_active() -> bool:
    """True iff DELFIN_FFFREE_AMIDE_VSEPR=1 or master flag set."""
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_AMIDE_VSEPR", "0") == "1"


def _heavy_adj(syms: Sequence[str], P: np.ndarray) -> List[List[int]]:
    """Heavy (non-metal, non-H) covalent neighbour list."""
    n = len(syms)
    adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        if syms[i] == "H" or _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * _bd._ideal_bond(syms[i], syms[j]):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _is_carbonyl_carbon(c: int, syms: Sequence[str], P: np.ndarray,
                        adj: List[List[int]]) -> bool:
    """C is a carbonyl carbon iff it has any O-neighbour at distance < 1.30 A."""
    if syms[c] != "C":
        return False
    for j in adj[c]:
        if syms[j] != "O":
            continue
        d = float(np.linalg.norm(P[c] - P[j]))
        if d < _CARBONYL_CO_MAX:
            return True
    return False


def find_amide_nitrogens(
    syms: Sequence[str], P: np.ndarray,
    exclude_idxs: Sequence[int] = (),
) -> List[Tuple[int, int]]:
    """Return list of (n_idx, carbonyl_c_idx) for each amide-N motif.

    An amide-N is any N atom that:
      - has 3 non-metal heavy neighbours (so substituent template applies),
      - has at least one neighbour C that is a carbonyl (C=O within 1.30 A).

    This captures amide, imide, carbamate, urea, sulfonamide-N (excluded by
    requiring carbonyl-C explicitly; sulfonamide skipped — sp2-N adjacent to S=O,
    different geometry).
    """
    n = len(syms)
    adj = _heavy_adj(syms, P)
    excl = set(int(i) for i in exclude_idxs)
    out: List[Tuple[int, int]] = []
    for ni in range(n):
        if syms[ni] != "N" or ni in excl:
            continue
        if len(adj[ni]) != 3:           # require sp2-amide-N (3 substituents)
            continue
        c_idx = None
        for j in adj[ni]:
            if _is_carbonyl_carbon(j, syms, P, adj):
                c_idx = j
                break
        if c_idx is not None:
            out.append((ni, c_idx))
    return out


def _project_onto_plane(point: np.ndarray, plane_pt: np.ndarray,
                        normal: np.ndarray) -> np.ndarray:
    """Orthogonal projection of `point` onto the plane through `plane_pt`
    with unit `normal`.  Returns the projected coordinate.  Deterministic.
    """
    v = point - plane_pt
    off = float(np.dot(v, normal))
    return point - off * normal


def enforce_amide_planarity(
    syms: Sequence[str], P: np.ndarray,
    metal_idx: int = -1, donor_idxs: Sequence[int] = (),
) -> Tuple[np.ndarray, int]:
    """Project amide nitrogens onto their carbonyl-substituent plane.

    The plane is fitted by SVD of the 3 N-neighbour displacement vectors:
    when the 3 vectors are coplanar (correct sp2-amide) the smallest singular
    value -> 0 and SVD returns the normal of that plane.  Projecting N onto
    that plane sets its 3rd-substituent OOP to 0.

    Donor atoms and the metal are NEVER moved (skipped when ni in donor_set).
    A rigid-H drag re-applies the same shift to any H whose only heavy
    neighbour is the moved N.

    Returns (P_new, n_fixed). n_fixed counts amide-N actually projected
    (motif found AND OOP > _OOP_THRESH).

    Env-gated default-OFF: when no flag set, returns (P_copy, 0) byte-identical.
    """
    P_old = np.asarray(P, dtype=float)
    P_new = P_old.copy()
    if not _flag_active():
        return P_new, 0
    if len(syms) == 0:
        return P_new, 0

    donor_set = set(int(d) for d in donor_idxs)
    if metal_idx >= 0:
        donor_set.add(int(metal_idx))

    n = len(syms)
    # build H-inclusive adjacency for the H-drag (geometric bonds incl. H)
    h_adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if _bd._is_metal(syms[i]) or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P_new[i] - P_new[j]))
            if d < _BOND_FACTOR * _bd._ideal_bond(syms[i], syms[j]):
                h_adj[i].append(j)
                h_adj[j].append(i)

    n_fixed = 0
    amides = find_amide_nitrogens(syms, P_new, exclude_idxs=tuple(donor_set))
    for ni, _c_idx in amides:
        # heavy neighbours of N (excluding metal)
        heavies = [k for k in h_adj[ni]
                   if syms[k] != "H" and not _bd._is_metal(syms[k])]
        if len(heavies) != 3:
            continue
        # fit plane through the 3 substituent positions (NOT through N).
        # The plane is defined by the 3 neighbour-points; we then project N
        # onto that plane.  This sets the OOP of the 3rd neighbour to 0
        # because the plane IS the 3-neighbour plane.
        Q = np.array([P_new[k] for k in heavies], dtype=float)
        centroid = Q.mean(axis=0)
        try:
            _, _, Vt = np.linalg.svd(Q - centroid, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        normal = Vt[-1]
        nn = float(np.linalg.norm(normal))
        if nn < 1e-9:
            continue
        normal = normal / nn

        off = float(np.dot(P_new[ni] - centroid, normal))
        if abs(off) < _OOP_THRESH:
            continue                   # already planar; nothing to do
        move = -off * normal           # move N into the plane
        P_candidate = P_new.copy()
        P_candidate[ni] = P_new[ni] + move
        # rigid-H drag: H whose only heavy neighbour is `ni`
        for h in range(n):
            if syms[h] != "H":
                continue
            heavies_of_h = [k for k in h_adj[h] if syms[k] != "H"]
            if heavies_of_h == [ni]:
                P_candidate[h] = P_new[h] + move

        # M-D invariant check (defence-in-depth; we don't touch donors so it
        # should hold by construction, but verify).
        if metal_idx >= 0:
            ok = True
            for d in donor_idxs:
                if d < 0 or d >= n:
                    continue
                d_old = float(np.linalg.norm(P_old[d] - P_old[metal_idx]))
                d_new = float(np.linalg.norm(P_candidate[d] - P_candidate[metal_idx]))
                if abs(d_old - d_new) > _MD_INVARIANT_TOL:
                    ok = False
                    break
            if not ok:
                continue
        P_new = P_candidate
        n_fixed += 1
    return P_new, n_fixed


__all__ = ["enforce_amide_planarity", "find_amide_nitrogens", "_flag_active"]
