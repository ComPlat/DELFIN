"""delfin.fffree.post_grip_corrector — Post-GRIP geometry correctors.

Phase 5 (User 2026-06-02, post-GRIP-hapto-VOLLPOOL forensik).  Three
chained geometric correctors that run AFTER GRIP polish and BEFORE the
final ff-free clash-relief:

  Fix 1 (POST_GRIP_FLATTEN): re-apply sp2 in-plane projection on every
    sp2 3-coordinate non-metal heavy atom that has drifted >0.02 A out
    of its 3-neighbour plane during GRIP optimisation.  Re-uses the
    same projection algorithm as the pre-GRIP sp2 corrector.  Affects
    F20 (aromatic-H planarity), F24 (amide-N planarity) and F23 (carbonyl-C / sp2 functional groups).

  Fix 2 (POST_GRIP_HAXIS_FIX): for each donor with a bonded H or methyl-H
    subgroup, check whether any bonded H sits on the M-D axis (within
    AXIS_PERP_MAX of the M-D line AND between M and D).  If yes, rotate
    the donor's substituent subtree 180 deg around the M-D axis to flip
    H away.  Pure axis rotation preserves the M-D distance exactly.
    Affects haxis (mode 2 = H_on_MD_axis, mode 3 = H_proximal_via_donor).

  Fix 3 (POST_GRIP_MDSHORT_FIX): scan all non-donor heavy atoms; if a
    secondary contact d(M, X) < 0.80 x ideal_bond(M, X) AND d < 1.45 x
    ideal (so the detector still considers it relevant), translate X
    along the (X - M) direction to exactly 0.85 x ideal.  Verify the
    parent-X bonds remain within 1.5 x ideal.  Affects mdshort
    (secondary contact mode).

All three correctors:
  - Verify M-D invariant (donor distances unchanged within 0.05 A) AFTER
    each step; roll back the step on failure.
  - Preserve topology: each step is purely geometric and only moves
    selected atoms minimally.
  - Are deterministic (sorted iteration order, no random sampling).
  - Universal: keyed on graph patterns (sp2-fingerprint, donor-H
    bonding, secondary-contact distance), not on SMILES strings.
  - Env-gated, default-OFF byte-identical when env unset.

Env flags:
  DELFIN_FFFREE_POST_GRIP_FLATTEN  (Fix 1, default 0)
  DELFIN_FFFREE_HAXIS_FIX           (Fix 2, default 0)
  DELFIN_FFFREE_MDSHORT_FIX         (Fix 3, default 0)
  DELFIN_FFFREE_POST_GRIP_ALL       (master, enables all three, default 0)

The master `post_grip_corrections(...)` entrypoint applies the enabled
fixes in order and returns the corrected coordinates.  It is the only
function the caller should use.
"""
from __future__ import annotations

import math
import os
from typing import List, Optional, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Reuse same thresholds as the matching detectors so the corrector
# operates on EXACTLY the cases the detector flags.
_F20_OOP_TOL = 0.20            # match sp2_h_planarity_check
_F24_OOP_TOL = 0.20            # match amide_imine_planarity_check
_SP2_FLATTEN_THRESH = 0.05     # only project if the atom is >0.05 A off
_HAXIS_PERP_MAX = 0.50         # match metric_h_axis.AXIS_PERP_MAX
_HAXIS_H_M_MAX = 1.50          # match metric_h_axis.AXIS_H_M_MAX
_HAXIS_PROX_DELTA = 0.05       # match metric_h_axis.PROX_DELTA
_HAXIS_PROX_H_M_MAX = 2.20     # match metric_h_axis.PROX_H_M_MAX
_MDSHORT_FACTOR = 0.80         # match metric_md_short_collapse.SHORT_FACTOR
_MDSHORT_TARGET_FACTOR = 0.85  # push back to this fraction of ideal
_MD_NEIGHBOUR_FACTOR = 1.45    # match metric_md_short_collapse
_MD_INVARIANT_TOL = 0.05       # M-D donor distance preservation tolerance
_DONORS = {"N", "O", "S", "P", "C", "F", "Cl", "Br", "I", "Se", "As", "Te"}


# ---------------------------------------------------------------------------
# Env-flag gates
# ---------------------------------------------------------------------------
def _flag(name: str, default: str = "0") -> bool:
    if os.environ.get("DELFIN_FFFREE_POST_GRIP_ALL", "0") == "1":
        return True
    return os.environ.get(name, default) == "1"


def flatten_active() -> bool:
    return _flag("DELFIN_FFFREE_POST_GRIP_FLATTEN")


def haxis_active() -> bool:
    return _flag("DELFIN_FFFREE_HAXIS_FIX")


def mdshort_active() -> bool:
    return _flag("DELFIN_FFFREE_MDSHORT_FIX")


def any_active() -> bool:
    return flatten_active() or haxis_active() or mdshort_active()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_adj(
    syms: Sequence[str], P: np.ndarray
) -> Tuple[List[List[int]], List[List[int]]]:
    """Geometric bond graph (1.30 x covalent-sum cutoff).

    Returns (adj_with_h, heavy_adj_no_metal).  Skips metal-X bonds for
    heavy_adj so the sp2 fingerprint doesn't trip on M-donor coordination.
    """
    n = len(syms)
    adj_with_h: List[List[int]] = [[] for _ in range(n)]
    heavy_adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * _bd._ideal_bond(syms[i], syms[j]):
                adj_with_h[i].append(j)
                adj_with_h[j].append(i)
                if (
                    syms[i] != "H"
                    and syms[j] != "H"
                    and not _bd._is_metal(syms[i])
                    and not _bd._is_metal(syms[j])
                ):
                    heavy_adj[i].append(j)
                    heavy_adj[j].append(i)
    return adj_with_h, heavy_adj


def _md_invariant_ok(
    P_old: np.ndarray, P_new: np.ndarray, donor_idxs: Sequence[int],
    metal_idx: int, tol: float = _MD_INVARIANT_TOL,
) -> bool:
    """Check that every donor-metal distance is preserved within tol."""
    if metal_idx < 0 or metal_idx >= len(P_old):
        return True
    for d in donor_idxs:
        if d < 0 or d >= len(P_old):
            continue
        d_old = float(np.linalg.norm(P_old[d] - P_old[metal_idx]))
        d_new = float(np.linalg.norm(P_new[d] - P_new[metal_idx]))
        if abs(d_old - d_new) > tol:
            return False
    return True


def _angle_deg(P: np.ndarray, i: int, j: int, k: int) -> float:
    v1 = P[i] - P[j]
    v2 = P[k] - P[j]
    n1 = float(np.linalg.norm(v1))
    n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return float("nan")
    c = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(math.degrees(math.acos(c)))


def _axis_rot(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rotation matrix for `theta` radians around `axis` (right-hand)."""
    a = axis / max(float(np.linalg.norm(axis)), 1e-12)
    c = math.cos(theta)
    s = math.sin(theta)
    x, y, z = float(a[0]), float(a[1]), float(a[2])
    return np.array([
        [c + x * x * (1 - c),   x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c),   y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ], dtype=float)


def _subtree_from(
    adj: List[List[int]], start: int, blocked: int, n: int
) -> Set[int]:
    """Atoms reachable from `start` over `adj` without crossing `blocked`."""
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        for j in adj[a]:
            if j == blocked or j in seen:
                continue
            seen.add(j)
            stack.append(j)
    return seen


# ---------------------------------------------------------------------------
# Fix 1: post-GRIP sp2 flatten
# ---------------------------------------------------------------------------
def apply_sp2_flatten(
    syms: Sequence[str], P: np.ndarray,
    metal_idx: int = -1, donor_idxs: Sequence[int] = (),
) -> Tuple[np.ndarray, int]:
    """Re-project sp2 3-coordinate non-metal heavy atoms onto their
    neighbour plane.

    Returns (P_new, n_projected).  Atoms in the donor set or metal are
    NOT moved (preserves coordination invariants).
    """
    P_new = np.asarray(P, dtype=float).copy()
    n = len(syms)
    if n == 0:
        return P_new, 0

    adj_with_h, heavy_adj = _build_adj(syms, P_new)
    donor_set = set(int(d) for d in donor_idxs) | ({int(metal_idx)} if metal_idx >= 0 else set())
    moved = 0

    # sp2 candidate detection: graph-based + geometric.  An atom is sp2 iff:
    #   (a) C/N/O with exactly 3 heavy neighbours
    #   (b) at least ONE bond to a double-bond partner (graph evidence via
    #       d(a, nbr) < 1.30 A for C=O, < 1.32 A for C=N, etc.) -- this rules
    #       out sp3 N (amines) which lack double-bond evidence; OR
    #   (c) part of a planar 5/6 ring (aromatic) -- already in_plane.
    # If neither (b) nor (c) holds, we additionally require the geometric
    # sp2 fingerprint (mean angle in [115, 130]); pure-graph sp2 (b) gets the
    # benefit of the doubt and the wider [100, 135] window so amide-N that has
    # already been pyramidalised by GRIP is recognised.
    SP2_LOW_GRAPH, SP2_HIGH_GRAPH = 100.0, 135.0    # for graph-sp2
    SP2_LOW_GEOM, SP2_HIGH_GEOM = 115.0, 130.0      # for geom-only sp2

    def _has_double_bond_partner(a):
        """Heuristic: does atom `a` have ANY heavy neighbour at a distance
        characteristic of a double bond?  Returns True for C=O, C=N, C=C,
        N=O etc.  Used to recognise sp2 atoms from the geometry graph alone.
        """
        for nb in heavy_adj[a]:
            d = float(np.linalg.norm(P_new[a] - P_new[nb]))
            # Double bond cutoffs (Å): C=O 1.21, C=N 1.28, C=C 1.34, N=O 1.21
            # Use a conservative 1.32 floor (matches RDKit's bond-order heuristic).
            pair = (syms[a], syms[nb])
            if pair in {("C","O"), ("O","C")} and d < 1.30:
                return True
            if pair in {("C","N"), ("N","C")} and d < 1.32:
                return True
            if pair in {("C","C")} and d < 1.36:
                return True
            if pair in {("N","O"), ("O","N")} and d < 1.30:
                return True
            if pair in {("N","N")} and d < 1.30:
                return True
        return False

    for a in range(n):
        if syms[a] == "H":
            continue
        if _bd._is_metal(syms[a]):
            continue
        if a in donor_set:
            continue                     # never move a coordinated donor
        if syms[a] not in {"C", "N", "O"}:
            continue
        nbrs = heavy_adj[a]
        if len(nbrs) != 3:
            continue
        # Skip if any neighbour is a metal (defence; heavy_adj already excludes)
        if any(_bd._is_metal(syms[k]) for k in nbrs):
            continue
        # Skip if any neighbour is in the donor set (preserve M-coord substituents)
        if any(k in donor_set for k in nbrs):
            continue
        # Geometric sp2 fingerprint
        a01 = _angle_deg(P_new, nbrs[0], a, nbrs[1])
        a02 = _angle_deg(P_new, nbrs[0], a, nbrs[2])
        a12 = _angle_deg(P_new, nbrs[1], a, nbrs[2])
        if any(math.isnan(x) for x in (a01, a02, a12)):
            continue
        mean_ang = (a01 + a02 + a12) / 3.0

        # If atom has graph-evidence of double bond (carbonyl C, amide-N
        # adjacent to C=O, etc.) we use the wider [100, 135] window so
        # pyramidalised sp2 atoms get caught.  Otherwise (sp3 candidates
        # like amine N at ~109 deg) we use the tighter [115, 130] window
        # so they are NOT flattened.
        if _has_double_bond_partner(a):
            lo, hi = SP2_LOW_GRAPH, SP2_HIGH_GRAPH
        else:
            # Special case: amide-N is N bonded to a carbonyl-C (C with C=O).
            # The N itself has no double-bond partner but is sp2.  Detect
            # via the carbonyl-C neighbour.
            is_amide_n = (syms[a] == "N" and any(
                syms[nb] == "C" and _has_double_bond_partner(nb)
                and any(syms[k] == "O" and float(np.linalg.norm(P_new[nb] - P_new[k])) < 1.30
                        for k in heavy_adj[nb])
                for nb in nbrs))
            if is_amide_n:
                lo, hi = SP2_LOW_GRAPH, SP2_HIGH_GRAPH
            else:
                lo, hi = SP2_LOW_GEOM, SP2_HIGH_GEOM
        if not (lo <= mean_ang <= hi):
            continue

        Q = np.array([P_new[k] for k in nbrs])
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

        off = float(np.dot(P_new[a] - centroid, normal))
        if abs(off) < _SP2_FLATTEN_THRESH:
            continue
        move = -off * normal
        P_new[a] = P_new[a] + move

        # Rigid-H drag for the single H whose only heavy neighbour is `a`
        for h in range(n):
            if syms[h] != "H":
                continue
            heavies = [k for k in adj_with_h[h] if syms[k] != "H"]
            if heavies == [a]:
                P_new[h] = P_new[h] + move
        moved += 1
    return P_new, moved


# ---------------------------------------------------------------------------
# Fix 2: post-GRIP H-on-M-D-axis rotation
# ---------------------------------------------------------------------------
def _find_haxis_violations(
    syms: Sequence[str], P: np.ndarray, metal_idx: int,
    donor_idxs: Sequence[int],
) -> List[Tuple[int, int, str]]:
    """Replicate metric_h_axis mode 2 (axis) and mode 3 (proximal-via-donor)
    detection, return list of (h_idx, donor_idx, mode).  Only flags Hs whose
    parent is in the donor set so we know who to rotate."""
    n = len(syms)
    if metal_idx < 0 or metal_idx >= n:
        return []
    out: List[Tuple[int, int, str]] = []
    seen_h: Set[int] = set()
    P = np.asarray(P, dtype=float)

    adj_with_h, _ = _build_adj(syms, P)

    # Mode 2: H on M-D axis (between M and D, perp < 0.50, d_HM < 1.50)
    for di in donor_idxs:
        if di < 0 or di >= n:
            continue
        v_md = P[di] - P[metal_idx]
        d_md = float(np.linalg.norm(v_md))
        if d_md < 1e-6:
            continue
        axis = v_md / d_md
        for h in range(n):
            if syms[h] != "H":
                continue
            if h in seen_h:
                continue
            w = P[h] - P[metal_idx]
            t = float(np.dot(w, axis))
            perp = float(np.linalg.norm(w - t * axis))
            d_hm = float(np.linalg.norm(P[h] - P[metal_idx]))
            if perp < _HAXIS_PERP_MAX and 0.2 < t < d_md - 0.2 and d_hm < _HAXIS_H_M_MAX:
                # find this H's parent
                heavies = [k for k in adj_with_h[h] if syms[k] != "H"
                           and not _bd._is_metal(syms[k])]
                if not heavies:
                    continue
                parent = heavies[0]
                if parent == di or parent in donor_idxs:
                    out.append((h, parent, "H_on_MD_axis"))
                    seen_h.add(h)

    # Mode 3: H proximal via donor (d_HM < 2.20 AND d_HM <= d_MX - 0.05)
    for h in range(n):
        if syms[h] != "H":
            continue
        if h in seen_h:
            continue
        # parent: closest bonded non-metal heavy atom
        parent_idx = -1
        parent_d = 1e9
        for k in adj_with_h[h]:
            if syms[k] == "H" or _bd._is_metal(syms[k]):
                continue
            d_hk = float(np.linalg.norm(P[h] - P[k]))
            if d_hk < parent_d:
                parent_d = d_hk
                parent_idx = k
        if parent_idx < 0:
            continue
        if syms[parent_idx] == "B" or syms[parent_idx] not in _DONORS:
            continue
        if parent_idx not in donor_idxs:
            continue
        d_mh = float(np.linalg.norm(P[h] - P[metal_idx]))
        if d_mh >= _HAXIS_PROX_H_M_MAX:
            continue
        d_mx = float(np.linalg.norm(P[parent_idx] - P[metal_idx]))
        if d_mh <= d_mx - _HAXIS_PROX_DELTA:
            out.append((h, parent_idx, "H_proximal_via_donor"))
            seen_h.add(h)
    return out


def _detect_hapto_donors(
    syms: Sequence[str], P: np.ndarray, metal_idx: int,
    donor_idxs: Sequence[int],
) -> Set[int]:
    """Return the subset of donor_idxs that belong to a hapto-pi cluster.

    A donor qualifies as hapto when there are >=2 other donors of the SAME
    element bonded to the same metal AND within 3.0 A of each other (i.e.
    a Cp / arene / allyl / diene ring).  These donors are excluded from
    haxis rotation -- their geometry is part of a rigid pi-frame.
    """
    donor_list = [int(d) for d in donor_idxs]
    hapto: Set[int] = set()
    for di in donor_list:
        same_elem_close = 0
        for dj in donor_list:
            if dj == di:
                continue
            if syms[di] != syms[dj]:
                continue
            if float(np.linalg.norm(P[di] - P[dj])) < 3.0:
                same_elem_close += 1
        if same_elem_close >= 2:
            hapto.add(di)
    return hapto


def _carries_stereocenter(
    adj: List[List[int]], rotate_set: Set[int], syms: Sequence[str],
) -> bool:
    """A 180 deg rotation around the M-D axis flips the chirality of any
    sp3 center inside the rotate set.  Detect sp3 C/Si/P/N atoms in
    rotate_set with 4 distinct heavy neighbours -- those are stereocentres.
    """
    for a in rotate_set:
        if syms[a] not in {"C", "Si", "N", "P"}:
            continue
        heavies = [k for k in adj[a] if syms[k] != "H" and not _bd._is_metal(syms[k])]
        if len(heavies) >= 4:
            return True
    return False


def apply_haxis_rotation(
    syms: Sequence[str], P: np.ndarray, metal_idx: int,
    donor_idxs: Sequence[int],
) -> Tuple[np.ndarray, int]:
    """For each detected haxis violation, rotate the donor's substituent
    subtree by 180 deg around the M-D axis.  M-D invariant is preserved
    exactly because the donor itself is ON the rotation axis.

    Skips:
      - chelate donors (rotation would move another donor)
      - hapto-pi donors (rigid pi-frame)
      - donors whose rotate set carries a stereocenter (180 deg flip
        inverts chirality)

    Returns (P_new, n_rotated).
    """
    P_new = np.asarray(P, dtype=float).copy()
    n = len(syms)
    if metal_idx < 0 or n == 0 or not donor_idxs:
        return P_new, 0

    adj_with_h, _ = _build_adj(syms, P_new)
    metal_set = {int(metal_idx)}
    donor_set = set(int(d) for d in donor_idxs)
    hapto_set = _detect_hapto_donors(syms, P_new, metal_idx, donor_idxs)

    rotated_donors: Set[int] = set()
    violations = _find_haxis_violations(syms, P_new, metal_idx, donor_idxs)

    for h_idx, donor_idx, mode in violations:
        if donor_idx in rotated_donors:
            continue
        if donor_idx in hapto_set:
            continue                     # hapto-pi: skip
        # Build the subtree from donor, NOT crossing the metal.
        # We want to rotate ONLY the H-side substituents, not the whole ligand.
        # Strategy: find the donor's heavy-atom neighbours EXCLUDING the metal;
        # for each of those parent atoms, take their subtree (not crossing the
        # donor) and rotate them.

        v_md = P_new[donor_idx] - P_new[metal_idx]
        d_md = float(np.linalg.norm(v_md))
        if d_md < 1e-6:
            continue
        axis = v_md / d_md

        # Find heavy non-metal neighbours of the donor.
        donor_atom_nbrs = [k for k in adj_with_h[donor_idx]
                           if syms[k] != "H" and not _bd._is_metal(syms[k])]
        if not donor_atom_nbrs:
            continue
        # No neighbour-subtree can contain another donor (chelate protection)
        any_other_donor = False
        nbr_subtrees = []
        for nbr in donor_atom_nbrs:
            st = _subtree_from(adj_with_h, nbr, donor_idx, n)
            if (donor_set - {donor_idx}) & st:
                any_other_donor = True
                break
            nbr_subtrees.append(st)
        if any_other_donor:
            continue

        # Build the union of all subtrees + their H atoms
        rotate_set: Set[int] = set()
        for st in nbr_subtrees:
            rotate_set |= st
        for h in range(n):
            if syms[h] != "H":
                continue
            heavies = [k for k in adj_with_h[h] if syms[k] != "H" and not _bd._is_metal(syms[k])]
            if heavies and set(heavies).issubset(rotate_set):
                rotate_set.add(h)
        rotate_set -= metal_set
        rotate_set.discard(donor_idx)
        if not rotate_set:
            continue

        # Stereo protection: skip if rotate_set carries an sp3 chiral atom.
        if _carries_stereocenter(adj_with_h, rotate_set, syms):
            continue

        # Save snapshot for rollback
        P_snap = P_new.copy()

        # Rotate 180 deg around (M, donor) axis: axis passes through donor.
        R = _axis_rot(axis, math.pi)
        pivot = P_new[donor_idx]
        for k in sorted(rotate_set):
            P_new[k] = (P_new[k] - pivot) @ R.T + pivot

        # Verify M-D invariant
        if not _md_invariant_ok(P_snap, P_new, donor_idxs, metal_idx):
            P_new = P_snap
            continue

        # Verify NO new inter-fragment clash: any rotated heavy atom must not
        # be within 1.5 A of any non-rotated, non-bonded heavy atom of another
        # ligand.  Approximation: pairwise distance check, exclude pairs that
        # already had a 1-2 bond.
        clash = False
        rotated_heavy = [k for k in rotate_set if syms[k] != "H"]
        for k in rotated_heavy:
            for j in range(n):
                if j in rotate_set:
                    continue
                if syms[j] == "H" or _bd._is_metal(syms[j]):
                    continue
                d_new = float(np.linalg.norm(P_new[k] - P_new[j]))
                if d_new < 1.5 and d_new < float(np.linalg.norm(P_snap[k] - P_snap[j])):
                    clash = True
                    break
            if clash:
                break
        if clash:
            P_new = P_snap
            continue

        # Verify the haxis violation for this h_idx is gone (or improved)
        new_vios = _find_haxis_violations(syms, P_new, metal_idx, donor_idxs)
        still_h = {(hi, di, m) for hi, di, m in new_vios if hi == h_idx}
        if still_h:
            P_new = P_snap
            continue
        rotated_donors.add(donor_idx)

    return P_new, len(rotated_donors)


# ---------------------------------------------------------------------------
# Fix 3: post-GRIP mdshort secondary-contact repulsion
# ---------------------------------------------------------------------------
def apply_mdshort_repulsion(
    syms: Sequence[str], P: np.ndarray, metal_idx: int,
    donor_idxs: Sequence[int],
) -> Tuple[np.ndarray, int]:
    """Push back any NON-donor heavy atom X that is within
    SHORT_FACTOR x ideal_bond(M, X) of the metal.  Push along (X - M) to
    exactly 0.85 x ideal so the detector no longer flags.  Donors are
    skipped (the M-D invariant guard already protects them).

    Returns (P_new, n_pushed).
    """
    P_new = np.asarray(P, dtype=float).copy()
    n = len(syms)
    if metal_idx < 0 or n == 0:
        return P_new, 0

    adj_with_h, _ = _build_adj(syms, P_new)
    donor_set = set(int(d) for d in donor_idxs)
    pushed = 0

    # Pre-compute the set of "donor-adjacent" atoms: any heavy non-metal atom
    # within 2 graph steps of a declared donor (1-2 or 1-3 to donor) is part of
    # the donor's first / second coordination shell and MUST NOT be pushed --
    # it carries the donor's natural geometry.  Pushing such an atom would
    # break the donor's local angle structure and induce funcgrp / interlig
    # regressions (smoke500 voll-pool 2026-06-02).
    donor_adjacent: Set[int] = set()
    for di in donor_set:
        if 0 <= di < n:
            donor_adjacent.add(di)
            for nbr in adj_with_h[di]:
                if syms[nbr] == "H" or _bd._is_metal(syms[nbr]):
                    continue
                donor_adjacent.add(nbr)
                for nbr2 in adj_with_h[nbr]:
                    if syms[nbr2] == "H" or _bd._is_metal(syms[nbr2]):
                        continue
                    donor_adjacent.add(nbr2)

    for x in range(n):
        if x == metal_idx:
            continue
        if syms[x] == "H":
            continue
        if _bd._is_metal(syms[x]):
            continue
        if syms[x] not in _DONORS:
            continue
        if x in donor_set:
            continue                     # M-D invariant protects this
        if x in donor_adjacent:
            continue                     # 1-2 / 1-3 to a donor -- preserve donor geometry
        ideal = _bd._ideal_bond(syms[metal_idx], syms[x])
        d = float(np.linalg.norm(P_new[x] - P_new[metal_idx]))
        if d < 1e-6:
            # zero-overlap; need a direction to push -- use a random axis
            push_dir = np.array([1.0, 0.0, 0.0])
        else:
            push_dir = (P_new[x] - P_new[metal_idx]) / d
        if d >= _MDSHORT_FACTOR * ideal:
            continue
        if d >= _MD_NEIGHBOUR_FACTOR * ideal:
            continue                     # outside detector's window already
        target_d = _MDSHORT_TARGET_FACTOR * ideal
        new_pos = P_new[metal_idx] + push_dir * target_d

        # Move x AND its H children (so X-H bonds preserved).  Heavy-parent
        # bonds may stretch but only if X had multiple heavy parents.  Verify.
        shift = new_pos - P_new[x]
        P_snap = P_new.copy()
        P_new[x] = new_pos
        for h in range(n):
            if syms[h] != "H":
                continue
            heavies = [k for k in adj_with_h[h] if syms[k] != "H"
                       and not _bd._is_metal(syms[k])]
            if heavies == [x]:
                P_new[h] = P_new[h] + shift

        # Topology check: every old bond of x to heavy non-metal must still be
        # within 1.5 x ideal (no bond-break) AND no new inter-atom near-contact
        # introduced.
        x_heavies = [k for k in adj_with_h[x] if syms[k] != "H"
                     and not _bd._is_metal(syms[k])]
        topo_ok = True
        for k in x_heavies:
            d_new = float(np.linalg.norm(P_new[x] - P_new[k]))
            if d_new > 1.5 * _bd._ideal_bond(syms[x], syms[k]):
                topo_ok = False
                break
        if not topo_ok:
            P_new = P_snap
            continue

        # Inter-atom clash check: no new heavy-heavy pair closer than 1.4 A.
        clash = False
        for k in range(n):
            if k == x or k == metal_idx or syms[k] == "H" or _bd._is_metal(syms[k]):
                continue
            if k in x_heavies:
                continue                 # bonded, expected to be close
            d_new = float(np.linalg.norm(P_new[x] - P_new[k]))
            d_old = float(np.linalg.norm(P_snap[x] - P_snap[k]))
            if d_new < 1.4 and d_new < d_old:
                clash = True
                break
        if clash:
            P_new = P_snap
            continue

        # M-D invariant
        if not _md_invariant_ok(P_snap, P_new, donor_idxs, metal_idx):
            P_new = P_snap
            continue

        pushed += 1
    return P_new, pushed


# ---------------------------------------------------------------------------
# Master entrypoint
# ---------------------------------------------------------------------------
def post_grip_corrections(
    syms: Sequence[str], P: np.ndarray, metal_idx: int,
    donor_idxs: Sequence[int],
) -> Tuple[np.ndarray, dict]:
    """Run all enabled post-GRIP correctors in fixed order.

    Returns (P_new, diagnostics_dict).  diagnostics has keys
    {'flatten_n', 'haxis_n', 'mdshort_n'} and 'any_enabled'.

    Order: flatten -> haxis -> mdshort.  This order matters because the
    sp2-flatten can introduce small M-D drift that the M-D guard catches;
    haxis rotation preserves M-D exactly; mdshort runs last so it sees the
    final geometry.
    """
    diag = {"flatten_n": 0, "haxis_n": 0, "mdshort_n": 0,
            "any_enabled": any_active()}
    if not any_active():
        return np.asarray(P, dtype=float), diag

    P_cur = np.asarray(P, dtype=float).copy()

    if flatten_active():
        P_try, n_moved = apply_sp2_flatten(
            syms, P_cur, metal_idx=metal_idx, donor_idxs=donor_idxs
        )
        if _md_invariant_ok(P_cur, P_try, donor_idxs, metal_idx):
            P_cur = P_try
            diag["flatten_n"] = n_moved

    if haxis_active():
        P_try, n_rot = apply_haxis_rotation(
            syms, P_cur, metal_idx=metal_idx, donor_idxs=donor_idxs
        )
        if _md_invariant_ok(P_cur, P_try, donor_idxs, metal_idx):
            P_cur = P_try
            diag["haxis_n"] = n_rot

    if mdshort_active():
        P_try, n_push = apply_mdshort_repulsion(
            syms, P_cur, metal_idx=metal_idx, donor_idxs=donor_idxs
        )
        if _md_invariant_ok(P_cur, P_try, donor_idxs, metal_idx):
            P_cur = P_try
            diag["mdshort_n"] = n_push

    return P_cur, diag
