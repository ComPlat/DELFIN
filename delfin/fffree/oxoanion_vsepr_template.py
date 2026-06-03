"""delfin.fffree.oxoanion_vsepr_template -- Build-time oxoanion VSEPR template.

Construction-side fix (User mandate 2026-06-03 task #34 re-activation):
voll-pool f8c9905 verdict showed nitrate_pct_files_with_violation +139 % and
nitrate_pct_no3_broken +93 % vs pgcorr-v3.  The original iter-32a-3 oxoanion
template was either not activated in the race-stack or never wired in.

This module implements a UNIVERSAL graph-based oxoanion detector and a strict
rigid-VSEPR projector for the common nitrate/perchlorate/sulfate/phosphate
geometries:

    NO3-   trigonal planar D3h   X-O ideal 1.245 A  O-X-O 120 deg
    ClO4-  tetrahedral    Td     X-O ideal 1.450 A  O-X-O 109.47 deg
    SO4^2- tetrahedral    Td     X-O ideal 1.490 A  O-X-O 109.47 deg
    PO4^3- tetrahedral    Td     X-O ideal 1.540 A  O-X-O 109.47 deg

Algorithm (deterministic, geometry-only):

  1. Identify the central atom X (N/Cl/S/P) by looking for a heavy atom whose
     non-metal, non-H heavy-neighbour set is *only* oxygens AND whose oxygen
     count matches the VSEPR template (3 -> D3h, 4 -> Td).
  2. Take the centroid of the current 3 (or 4) oxygens and the centroid->X
     direction to anchor a local frame.
  3. Rebuild the oxygen positions on the ideal-template directions, rotated
     into the local frame so that the position of the most M-coordinating
     oxygen is preserved (when a metal is present and one O is closer to it
     than the others) -- this keeps the M-O coordination invariant.
  4. M-D INVARIANT: the central atom X is NEVER moved.  Donor oxygens (any
     O whose index is in ``donor_idxs``) are NEVER moved; the template is
     applied around them so the metal-donor distance is preserved by
     construction.  If MORE than one O is a donor (kappa^2-bidentate
     nitrate/carboxylate-like) the template is skipped (too constrained
     -- the geometry is already fixed by two donors).

The fix is universal across all metals and ligand environments (the
detector is purely graph-geometric, no atom-specific lookup tables).
FF-free.  Default-OFF byte-identical.

Env flags:

    DELFIN_FFFREE_OXOANION_VSEPR=1         (per-fix activation)
    DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1   (global construction master)

The module exports:

    enforce_oxoanion_vsepr(syms, P, metal_idx=-1, donor_idxs=()) -> (P_new, n_fixed)
    find_oxoanions(syms, P)               -> list[ (x_idx, [o_idxs]) ]
"""
from __future__ import annotations

import os
from typing import List, Sequence, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Ideal X-O bond lengths and target O-X-O angles per oxoanion family.
# Sources: typical CSD/COD medians, also the textbook VSEPR ideal angles.
_OXOANION_TEMPLATES = {
    # central_element : (n_oxygens, X-O ideal A, O-X-O ideal deg, geometry tag)
    "N":  (3, 1.245, 120.000, "D3h"),
    "Cl": (4, 1.450, 109.471, "Td"),
    "S":  (4, 1.490, 109.471, "Td"),
    "P":  (4, 1.540, 109.471, "Td"),
}

# Bond cutoff factor (multiplier of covalent-sum) for the heavy-neighbour
# perception graph.  Matches the convention used in amide_vsepr_template.
_BOND_FACTOR = 1.30

# OOP / angle threshold above which we project.  We pick 5 deg on the
# O-X-O angle  -- well below the +139 % "broken" cutoff in the nitrate
# detector  -- so any meaningfully-distorted oxoanion is repaired but
# perfectly-built ones are left byte-identical.
_ANGLE_THRESH_DEG = 5.0

# M-D invariant verification tolerance (mirror amide_vsepr_template).
_MD_INVARIANT_TOL = 0.05


def _flag_active() -> bool:
    """True iff DELFIN_FFFREE_OXOANION_VSEPR=1 or a master flag is set."""
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_OXOANION_VSEPR", "0") == "1"


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


def find_oxoanions(
    syms: Sequence[str], P: np.ndarray,
) -> List[Tuple[int, List[int], str]]:
    """Return list of (x_idx, [o_idxs sorted], geometry_tag).

    A central atom X is reported iff:
      - X is in {N, Cl, S, P},
      - X has EXACTLY n_oxygens heavy-neighbours (n_oxygens per template),
      - ALL of X's heavy-neighbours are oxygen (no C-N-O etc. -- pure
        oxoanion only; this excludes nitro/sulfonyl-attached-to-C groups
        which have their own templates),
      - none of those oxygens has another heavy-neighbour besides X
        (terminal oxygens; protonated O-H is allowed because H is not in
        the heavy graph).

    Universal, deterministic, geometry-only.
    """
    n = len(syms)
    adj = _heavy_adj(syms, P)
    out: List[Tuple[int, List[int], str]] = []
    for xi in range(n):
        s = syms[xi]
        tpl = _OXOANION_TEMPLATES.get(s)
        if tpl is None:
            continue
        n_o_ideal, _bl, _ang, geom = tpl
        nbrs = adj[xi]
        if len(nbrs) != n_o_ideal:
            continue
        if any(syms[j] != "O" for j in nbrs):
            continue
        # require terminal oxygens (no further heavy-neighbours besides xi)
        ok = True
        for j in nbrs:
            others = [k for k in adj[j] if k != xi]
            if others:
                ok = False
                break
        if not ok:
            continue
        out.append((xi, sorted(int(j) for j in nbrs), geom))
    return out


def _td_directions() -> np.ndarray:
    """Four ideal Td unit vectors (alternating cube corners)."""
    a = 1.0 / np.sqrt(3.0)
    return np.array([
        [ a,  a,  a],
        [ a, -a, -a],
        [-a,  a, -a],
        [-a, -a,  a],
    ], dtype=float)


def _d3h_directions() -> np.ndarray:
    """Three ideal D3h unit vectors (in xy-plane, 120 deg apart)."""
    return np.array([
        [1.0, 0.0, 0.0],
        [np.cos(np.deg2rad(120.0)), np.sin(np.deg2rad(120.0)), 0.0],
        [np.cos(np.deg2rad(240.0)), np.sin(np.deg2rad(240.0)), 0.0],
    ], dtype=float)


def _kabsch_rotation(P_src: np.ndarray, P_tgt: np.ndarray) -> np.ndarray:
    """Optimal rotation matrix that maps rows of P_src onto rows of P_tgt
    (both already centred).  Returns 3x3 R with det(R) > 0.
    """
    H = P_src.T @ P_tgt
    U, _S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def _o_o_angles_deg(P_center: np.ndarray, P_oxygens: np.ndarray) -> np.ndarray:
    """All pairwise O-X-O angles in degrees (upper triangle)."""
    v = P_oxygens - P_center
    nrm = np.linalg.norm(v, axis=1)
    bad = (nrm < 1e-9)
    if bool(bad.any()):
        nrm = np.where(bad, 1e-9, nrm)
    u = v / nrm[:, None]
    ang = []
    n = len(u)
    for i in range(n):
        for j in range(i + 1, n):
            c = float(np.clip(np.dot(u[i], u[j]), -1.0, 1.0))
            ang.append(np.degrees(np.arccos(c)))
    return np.asarray(ang, dtype=float)


def _max_angle_deviation(
    P_center: np.ndarray, P_oxygens: np.ndarray, ideal_deg: float,
) -> float:
    """Max |angle - ideal_deg| over all O-X-O pairs.  0 if no oxygens."""
    a = _o_o_angles_deg(P_center, P_oxygens)
    if a.size == 0:
        return 0.0
    return float(np.max(np.abs(a - ideal_deg)))


def enforce_oxoanion_vsepr(
    syms: Sequence[str], P: np.ndarray,
    metal_idx: int = -1, donor_idxs: Sequence[int] = (),
) -> Tuple[np.ndarray, int]:
    """Project oxoanion oxygens onto their ideal-VSEPR geometry.

    Algorithm:
      1. Find oxoanion motifs (X = N/Cl/S/P with all-O neighbours).
      2. For each motif, build the ideal template (Td / D3h) at the
         observed X-O average bond length (so we don't fight the rest of
         the structure's bond-length context).
      3. Rotate the template onto the current oxygen-set via Kabsch +
         keep the central atom X fixed.
      4. If a single oxygen is a donor, anchor that oxygen first (Kabsch
         with extra weight) so the M-O coordination stays invariant.
         If two or more oxygens are donors, SKIP the projection (the
         geometry is already constrained by two donors -- forcing the
         template would destroy the chelate).
      5. M-D invariant check before commit.

    Donor atoms and the metal are NEVER moved (skipped explicitly).
    Universal, deterministic, byte-identical to HEAD when no env flag set.

    Returns (P_new, n_fixed). n_fixed = number of oxoanions actually
    projected (motif found AND max O-X-O angle deviation > _ANGLE_THRESH_DEG).
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

    n_fixed = 0
    motifs = find_oxoanions(syms, P_new)
    for xi, o_idxs, geom in motifs:
        # Skip if multiple oxygens are donors (over-constrained chelate).
        n_donor_o = sum(1 for o in o_idxs if o in donor_set)
        if n_donor_o > 1:
            continue
        # Skip if the central X is a donor (rare -- M-N(=O)... eta1 nitrate
        # via N).  In that case we don't reorient the oxygens because the
        # M-N direction is the constraint.
        if xi in donor_set:
            continue

        # Compute current angle deviation; skip if already ideal.
        n_o, bl_ideal, ang_ideal, _ = _OXOANION_TEMPLATES[syms[xi]]
        if len(o_idxs) != n_o:
            continue
        P_o = P_new[o_idxs]
        dev = _max_angle_deviation(P_new[xi], P_o, ang_ideal)
        if dev < _ANGLE_THRESH_DEG:
            continue

        # Use the observed average X-O bond length (preserves the rest of
        # the structure's scale; the cost of a 0.02 A bond drift is much
        # smaller than the O-X-O angle restoration).
        bls = np.linalg.norm(P_o - P_new[xi], axis=1)
        bl_obs = float(np.mean(bls))
        # Clamp away from pathological collapse / explosion.
        bl_use = float(np.clip(bl_obs, 0.85 * bl_ideal, 1.20 * bl_ideal))

        # Build ideal template at X.
        if geom == "Td":
            dirs = _td_directions()
        else:  # D3h
            dirs = _d3h_directions()
        tpl = dirs * bl_use  # shape (n_o, 3)

        # Source vectors (current, centred at X).
        src = P_o - P_new[xi]

        # If a single oxygen is a donor, the geometry is over-constrained
        # if we move X (X and O_donor both anchored).  We therefore keep
        # X fixed, keep O_donor fixed, and rotate the remaining (n_o - 1)
        # template oxygens around the X->O_donor axis at the ideal angle
        # (120 deg in D3h, 109.47 deg in Td) at the average X-O bond length,
        # choosing the azimuthal phase that maximises overlap with the
        # current other oxygens.  Deterministic.
        donor_o = next((o for o in o_idxs if o in donor_set), None)
        P_candidate = P_new.copy()
        if donor_o is not None:
            v_d = P_new[donor_o] - P_new[xi]
            d_donor = float(np.linalg.norm(v_d))
            if d_donor < 1e-9:
                continue
            axis = v_d / d_donor
            # Build an orthonormal frame (axis, e1, e2).
            # Pick a reference vector not parallel to axis.
            ref = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 \
                else np.array([0.0, 1.0, 0.0])
            e1 = ref - np.dot(ref, axis) * axis
            e1n = float(np.linalg.norm(e1))
            if e1n < 1e-9:
                continue
            e1 = e1 / e1n
            e2 = np.cross(axis, e1)
            # Ideal template direction for a non-donor O: it sits at the
            # ideal O-X-O angle from the donor direction, at d_use radius.
            # We keep the donor bond length intact (matches MO invariant)
            # and use the observed mean for the other oxygens.
            theta = np.deg2rad(ang_ideal)
            # Each non-donor O at: r * (cos(theta) * axis + sin(theta) *
            # (cos(phi_k) * e1 + sin(phi_k) * e2)).
            n_other = n_o - 1
            # phi spacing: equispaced around the axis, deterministic.
            # For D3h with one donor on +axis: the other 2 sit at theta=120
            # from axis, spaced 180 deg in phi -> back-to-back (planar).
            # For Td with one donor on +axis: the other 3 sit at theta
            # =109.47 from axis, spaced 120 deg in phi.
            phis_base = np.array(
                [2.0 * np.pi * k / n_other for k in range(n_other)],
                dtype=float,
            )
            # Identify the indices of the non-donor oxygens (in o_idxs order).
            other_idx_in_olist = [k for k in range(n_o) if o_idxs[k] != donor_o]
            other_src = np.array(
                [src[k] for k in other_idx_in_olist], dtype=float,
            )
            # Find the phase phi0 (and assignment) that minimises sum of
            # squared distances between template-placed oxygens and the
            # current oxygens.  Deterministic search over phi0 grid +
            # nearest-neighbour assignment.
            best_rmsd = float("inf")
            best_phi0 = 0.0
            best_assign: List[int] = list(range(n_other))
            n_grid = 360 * 2  # 0.5 deg resolution -- deterministic
            for gi in range(n_grid):
                phi0 = 2.0 * np.pi * gi / n_grid
                cand = np.zeros((n_other, 3), dtype=float)
                for k in range(n_other):
                    phi = phis_base[k] + phi0
                    cand[k] = bl_use * (
                        np.cos(theta) * axis
                        + np.sin(theta) * (np.cos(phi) * e1 + np.sin(phi) * e2)
                    )
                # nearest-neighbour assignment current->template
                # (greedy, deterministic by index)
                used = np.zeros(n_other, dtype=bool)
                assign = [0] * n_other
                rmsd = 0.0
                for ci in range(n_other):
                    best_t = -1
                    best_d2 = float("inf")
                    for ti in range(n_other):
                        if used[ti]:
                            continue
                        d2 = float(np.sum((other_src[ci] - cand[ti]) ** 2))
                        if d2 < best_d2 - 1e-12:
                            best_d2 = d2
                            best_t = ti
                    if best_t < 0:
                        best_t = next(ti for ti in range(n_other)
                                      if not used[ti])
                    used[best_t] = True
                    assign[ci] = best_t
                    rmsd += best_d2
                if rmsd < best_rmsd - 1e-12:
                    best_rmsd = rmsd
                    best_phi0 = phi0
                    best_assign = list(assign)
            # Apply the best phase.
            phi0 = best_phi0
            for ci, ki in enumerate(other_idx_in_olist):
                phi = phis_base[best_assign[ci]] + phi0
                disp = bl_use * (
                    np.cos(theta) * axis
                    + np.sin(theta) * (np.cos(phi) * e1 + np.sin(phi) * e2)
                )
                P_candidate[o_idxs[ki]] = P_new[xi] + disp
            # Donor oxygen never moved.
        else:
            # No donor among oxygens: free Kabsch alignment of the ideal
            # template onto the current oxygens.  Greedy permutation +
            # Kabsch rotation centred on X.
            template_used = np.zeros(n_o, dtype=bool)
            tpl_perm = np.zeros_like(tpl)
            for k in range(n_o):
                best = -1
                best_d = float("inf")
                for ti in range(n_o):
                    if template_used[ti]:
                        continue
                    d = float(np.linalg.norm(src[k] - tpl[ti]))
                    if d < best_d - 1e-12:
                        best_d = d
                        best = ti
                if best < 0:
                    best = next(ti for ti in range(n_o)
                                if not template_used[ti])
                tpl_perm[k] = tpl[best]
                template_used[best] = True
            try:
                R = _kabsch_rotation(tpl_perm, src)
            except np.linalg.LinAlgError:
                continue
            rotated = tpl_perm @ R.T
            for k in range(n_o):
                P_candidate[o_idxs[k]] = P_new[xi] + rotated[k]

        # NaN / finite guard.
        if not np.all(np.isfinite(P_candidate)):
            continue

        # M-D invariant check (defence-in-depth; donor not touched but
        # verify nonetheless, and verify the central X is not moved).
        if float(np.linalg.norm(P_candidate[xi] - P_new[xi])) > 1e-9:
            continue
        if metal_idx >= 0:
            ok = True
            for d in donor_idxs:
                if d < 0 or d >= len(syms):
                    continue
                d_old = float(np.linalg.norm(P_old[d] - P_old[metal_idx]))
                d_new = float(np.linalg.norm(P_candidate[d] - P_candidate[metal_idx]))
                if abs(d_old - d_new) > _MD_INVARIANT_TOL:
                    ok = False
                    break
            if not ok:
                continue

        # Net-improvement gate: the new max-angle-deviation must be
        # strictly smaller than the old one (accept-if-better).
        new_dev = _max_angle_deviation(
            P_candidate[xi], P_candidate[o_idxs], ang_ideal,
        )
        if new_dev >= dev - 1e-6:
            continue

        P_new = P_candidate
        n_fixed += 1
    return P_new, n_fixed


__all__ = [
    "enforce_oxoanion_vsepr",
    "find_oxoanions",
    "_flag_active",
]
