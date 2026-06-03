"""delfin.fffree.mddir_alignment — Build-time M-D direction alignment.

Construction-side fix #2 (User mandate 2026-06-03): mddir voll-pool 17.55%
file-rate.  3 modes from the metric_md_direction detector:
  - aromatic_face : M over donor pi-face (M-D vector ~|| ring-normal)
  - donor_linearized : max(M-D-X) > 168 deg (donor sp instead of sp3)
  - sp3_anti : M on substituent side (M-D...substituent-centroid angle < 90)

ALL three modes are correctable by a pure rotation of the donor's substituent
subtree around the M-D axis (preserves the M-D distance exactly).

  Aromatic donor (planar, 2-3 in-plane neighbours):
    The lone pair lies IN the donor plane perpendicular to the substituent
    direction.  When M is over the pi-face (theta = angle(M->D, ring-normal)
    is small), we rotate the substituent subtree around the M-D axis by 90
    degrees so the ring normal becomes perpendicular to M-D = M in-plane.

  sp3 donor (2-4 substituents non-coplanar):
    Anti-substituent-centroid mode: rotate the substituent subtree by 180 deg
    around the M-D axis so the substituent centroid swings to the opposite
    side of D (M side -> substituent side becomes wrong-side -> after 180
    rotation the substituent centroid is now anti to M).

  donor_linearized:
    If max(M-D-X) > 168 deg the donor is sp not sp3; this is a substituent-
    placement defect that the build cannot fix by pure rotation around M-D
    (the linear D-X axis is the axis itself).  We DETECT and skip — leave
    these to post-processing.  Logged as `n_skipped_linearized`.

M-D INVARIANT: axis rotation preserves donor-metal distance exactly.

Universal, deterministic, geometry-only.  Env-gated default-OFF byte-identical.

API:
    align_donor_lonepairs(syms, P, metal_idx, donor_idxs) -> (P_new, n_fixed, n_skipped)

Env flag:
    DELFIN_FFFREE_MDDIR_ALIGN=1    (per-fix activation)
    DELFIN_FFFREE_CONSTRUCTION_FIX_ALL=1   (master activation)
"""
from __future__ import annotations

import math
import os
from typing import List, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Match metric_md_direction tolerances so we correct exactly the cases the
# detector flags.  Slightly more aggressive thresholds (we want zero post-build
# flags); the actual rotation respects the M-D axis exactly.
_FACE_TOL_DEG = 35.0           # match metric_md_direction.FACE_TOL_DEG
_SP3_ANTI_TOL_DEG = 90.0       # match metric_md_direction.SP3_ANTI_TOL_DEG
_LINEAR_TOL_DEG = 168.0        # match metric_md_direction.LINEAR_TOL_DEG
_COPLANAR_TOL = 0.25           # planarity SVD ratio cutoff
_MD_BOND_FACTOR = 1.40         # M-D bond detection cutoff (× cov sum)
_DONOR_ELEMENTS = {"N", "O", "P", "S"}
_BOND_FACTOR = 1.30            # general bond perception (× cov sum)
_MD_INVARIANT_TOL = 0.05       # post-rotation invariant check


def _flag_active() -> bool:
    if os.environ.get("DELFIN_FFFREE_CONSTRUCTION_FIX_ALL", "0") == "1":
        return True
    return os.environ.get("DELFIN_FFFREE_MDDIR_ALIGN", "0") == "1"


def _heavy_h_neighbors(syms: Sequence[str], P: np.ndarray, di: int,
                       exclude_metals: bool = True) -> List[int]:
    """Non-metal covalent neighbours of donor `di` (incl. H)."""
    n = len(syms)
    nbrs: List[int] = []
    for j in range(n):
        if j == di:
            continue
        if exclude_metals and _bd._is_metal(syms[j]):
            continue
        d = float(np.linalg.norm(P[di] - P[j]))
        if d < _BOND_FACTOR * _bd._ideal_bond(syms[di], syms[j]):
            nbrs.append(j)
    return nbrs


def _full_adj(syms: Sequence[str], P: np.ndarray) -> List[List[int]]:
    """Full bond graph (incl. H, excl. metals)."""
    n = len(syms)
    adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        if _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < _BOND_FACTOR * _bd._ideal_bond(syms[i], syms[j]):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _subtree_atoms(adj: List[List[int]], start: int, blocked: int,
                   n: int) -> Set[int]:
    """All atoms reachable from `start` without crossing `blocked`.

    Used to define which atoms ROTATE with the donor's substituent block:
    everything on the donor's side of the M-D axis (NOT the metal, NOT the
    other ligands).
    """
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


def _axis_rot(axis: np.ndarray, theta: float) -> np.ndarray:
    """Right-hand rotation matrix for `theta` radians around `axis`."""
    a = axis / max(float(np.linalg.norm(axis)), 1e-12)
    c = math.cos(theta)
    s = math.sin(theta)
    x, y, z = float(a[0]), float(a[1]), float(a[2])
    return np.array([
        [c + x * x * (1 - c),     x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
        [y * x * (1 - c) + z * s, c + y * y * (1 - c),     y * z * (1 - c) - x * s],
        [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
    ], dtype=float)


def _is_donor_planar(syms: Sequence[str], P: np.ndarray, di: int,
                     nbrs: Sequence[int]) -> Tuple[bool, np.ndarray]:
    """Returns (planar, n_hat) for donor `di`."""
    if len(nbrs) not in (2, 3):
        return False, np.zeros(3)
    pts = set(j for j in nbrs if syms[j] != "H")
    if len(pts) < 3:
        # extend to second-shell heavy neighbours (aromatic ring case)
        for j in nbrs:
            for k in _heavy_h_neighbors(syms, P, j):
                if syms[k] != "H" and k != di:
                    pts.add(k)
    pts.discard(di)
    if len(pts) < 3:
        return False, np.zeros(3)
    V = np.array([P[j] - P[di] for j in pts], dtype=float)
    try:
        _, s, vt = np.linalg.svd(V)
    except np.linalg.LinAlgError:
        return False, np.zeros(3)
    if len(s) < 3 or s[0] < 1e-9:
        return False, np.zeros(3)
    planar = (s[2] / s[0]) < _COPLANAR_TOL
    n_hat = vt[2]
    nn = float(np.linalg.norm(n_hat))
    if nn < 1e-9:
        return False, np.zeros(3)
    return planar, n_hat / nn


def _angle_between_unit(u: np.ndarray, v: np.ndarray) -> float:
    """Angle in degrees between two non-zero vectors (folded to [0,180])."""
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu < 1e-9 or nv < 1e-9:
        return float("nan")
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return math.degrees(math.acos(c))


def align_donor_lonepairs(
    syms: Sequence[str], P: np.ndarray,
    metal_idx: int, donor_idxs: Sequence[int],
) -> Tuple[np.ndarray, int, int]:
    """Per-donor rotate the substituent subtree around the M-D axis to put the
    donor lone pair on the M side.

    Modes handled:
      aromatic_face   -> rotate subtree 90 deg around M-D axis (ring flips
                          from face-on to edge-on -> M in plane).
      sp3_anti        -> rotate subtree 180 deg around M-D axis (substituent
                          centroid swings to opposite side -> M anti).
      donor_linearized -> SKIPPED (cannot be fixed by axis rotation; the linear
                          D-X axis lies on the rotation axis).

    Returns (P_new, n_fixed, n_skipped).

    Env-gated default-OFF: byte-identical when no flag set.
    """
    P_old = np.asarray(P, dtype=float)
    P_new = P_old.copy()
    if not _flag_active():
        return P_new, 0, 0
    n = len(syms)
    if metal_idx < 0 or metal_idx >= n:
        return P_new, 0, 0

    adj = _full_adj(syms, P_new)
    n_fixed = 0
    n_skipped = 0
    # Deterministic iteration: sorted donor indices
    for di in sorted(donor_idxs):
        if di < 0 or di >= n:
            continue
        if syms[di] not in _DONOR_ELEMENTS:
            continue
        v_md_raw = P_new[metal_idx] - P_new[di]
        nv = float(np.linalg.norm(v_md_raw))
        if nv < 1e-6:
            continue
        v_md = v_md_raw / nv         # unit vector D -> M (matches detector convention)
        nbrs = _heavy_h_neighbors(syms, P_new, di)
        if len(nbrs) == 0:
            continue
        # Build subtree (rotates with donor) — start at each non-metal neighbour
        # and grow without crossing di.  Donor itself stays fixed.
        subtree: Set[int] = set()
        for nb in nbrs:
            subtree |= _subtree_atoms(adj, nb, di, n)
        # NEVER move the metal or other donors
        subtree.discard(metal_idx)
        for od in donor_idxs:
            subtree.discard(int(od))
        if not subtree:
            continue

        planar, n_hat = _is_donor_planar(syms, P_new, di, nbrs)
        # Check linearization first (max M-D-X angle).
        max_mdx = 0.0
        for j in nbrs:
            u = P_new[j] - P_new[di]
            a = _angle_between_unit(P_new[metal_idx] - P_new[di], u)
            if not math.isnan(a) and a > max_mdx:
                max_mdx = a
        if max_mdx > _LINEAR_TOL_DEG:
            n_skipped += 1
            continue

        do_rotation = False
        theta = 0.0
        # Rotation axis through donor (preserves M-D distance exactly).
        # For 180-deg flips we can use any axis perpendicular to M-D; for the
        # aromatic-face 90-deg fix we need an axis perpendicular to BOTH M-D
        # and the ring normal so the rotation tips the ring from face-on to
        # edge-on.  The M-D axis itself is NOT useful when the ring normal IS
        # M-D (the rotation is a no-op).
        axis_rot_vec: np.ndarray = np.zeros(3)
        if planar:
            # aromatic_face: theta = angle(D->M, ring-normal) folded [0,90]
            # Detector uses absolute dot -> symmetric in sign anyway.
            cosv = float(np.clip(abs(np.dot(v_md, n_hat)), 0.0, 1.0))
            face_theta = math.degrees(math.acos(cosv))
            if face_theta < _FACE_TOL_DEG:
                # rotate 90 deg around (M->D) x ring-normal so the ring tips
                # from face-on to edge-on -> M moves into the ring plane
                axis_candidate = np.cross(v_md, n_hat)
                if float(np.linalg.norm(axis_candidate)) < 1e-6:
                    # M-D and ring-normal are exactly parallel: pick any
                    # perpendicular axis (use the donor->first-substituent
                    # direction projected perpendicular to M-D).
                    nb0 = nbrs[0]
                    u = P_new[nb0] - P_new[di]
                    u = u - float(np.dot(u, v_md)) * v_md
                    nu = float(np.linalg.norm(u))
                    if nu < 1e-6:
                        # extreme degenerate case: try +x perpendicular to M-D
                        axis_candidate = np.array([v_md[1], -v_md[0], 0.0])
                        if float(np.linalg.norm(axis_candidate)) < 1e-6:
                            axis_candidate = np.array([0.0, v_md[2], -v_md[1]])
                    else:
                        axis_candidate = u / nu
                axis_rot_vec = axis_candidate
                theta = math.radians(90.0)
                do_rotation = True
        elif 2 <= len(nbrs) <= 4:
            c = np.mean([P_new[j] for j in nbrs], axis=0)
            v_sub = c - P_new[di]
            ns = float(np.linalg.norm(v_sub))
            if ns >= 1e-6:
                v_sub_u = v_sub / ns
                # phi = angle(D->M, D->centroid) -- detector convention.
                # Flag when phi < 90 (M and substituents on SAME side of D);
                # correct ~180 means M and centroid on OPPOSITE sides (anti).
                phi = math.degrees(math.acos(float(np.clip(
                    np.dot(v_md, v_sub_u), -1.0, 1.0))))
                if phi < _SP3_ANTI_TOL_DEG:
                    # 180-deg flip around an axis perpendicular to M-D so the
                    # substituent centroid swings to the opposite side of D.
                    # Pick axis = (M-D) x (D->centroid) -> perpendicular to
                    # both -> 180 swings centroid through M-D axis to anti.
                    axis_candidate = np.cross(v_md, v_sub_u)
                    if float(np.linalg.norm(axis_candidate)) < 1e-6:
                        # parallel case: any axis perpendicular to M-D works
                        axis_candidate = np.array([v_md[1], -v_md[0], 0.0])
                        if float(np.linalg.norm(axis_candidate)) < 1e-6:
                            axis_candidate = np.array([0.0, v_md[2], -v_md[1]])
                    axis_rot_vec = axis_candidate
                    theta = math.radians(180.0)
                    do_rotation = True

        if not do_rotation:
            continue

        # Rotation around `axis_rot_vec` through D (donor stays fixed -> M-D
        # distance preserved EXACTLY by construction).
        R = _axis_rot(axis_rot_vec, theta)
        P_candidate = P_new.copy()
        center = P_new[di]
        for a in sorted(subtree):
            v = P_new[a] - center
            P_candidate[a] = center + R @ v

        # M-D invariant check: donor IS the rotation center -> distance to
        # metal preserved exactly.  But verify all donors stayed put
        # numerically.
        ok = True
        for od in donor_idxs:
            if od < 0 or od >= n:
                continue
            d_old = float(np.linalg.norm(P_old[od] - P_old[metal_idx]))
            d_new = float(np.linalg.norm(P_candidate[od] - P_candidate[metal_idx]))
            if abs(d_old - d_new) > _MD_INVARIANT_TOL:
                ok = False
                break
        if not ok:
            continue
        P_new = P_candidate
        n_fixed += 1
    return P_new, n_fixed, n_skipped


__all__ = ["align_donor_lonepairs", "_flag_active"]
