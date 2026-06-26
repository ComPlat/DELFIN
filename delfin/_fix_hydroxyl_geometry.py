"""Targeted #329 fixer — hydroxyl (C-O-H) geometry correction.

Data-driven root finding (2026-06-23) over 2022 V2R structures (155 hydroxyl
groups): the C-O-H angle is systematically too wide (43% > 115 deg, median
111.9; ideal ~108.5) and the C-O bond too short (26% < 1.36 A; alcohol ideal
~1.43, phenol ~1.36).  The existing sp3-H tetrahedrality fixer (F19) targets
109.47 deg with a 10 deg tolerance, so it only fires above ~119.5 deg (the
moderately-wide bulk is untouched) and it NEVER corrects the C-O bond length.

This fixer repairs BOTH, only on pendant (non-coordinating) hydroxyls:
  * C-O bond length -> hybridisation-aware ideal (carboxyl 1.31 / aromatic
    phenol 1.36 / sp3 alcohol 1.43), by translating the rigid {O,H} unit along
    the C->O axis (O stays bonded to the same C; O-H + C-O-H preserved by the
    move, the angle is then set separately).
  * C-O-H angle -> 108.5 deg, by rotating H about O in the C-O-H plane
    (O-H length preserved).

Doctrine (same family as _fix_sp3_h_tetrahedrality):
- Metals NEVER moved; a hydroxyl O bonded to any metal is SKIPPED (coordinating
  O is a donor — its position is a coordination property, not ours to touch).
- Only the hydroxyl O and its single H move; the parent C and everything else
  stay put -> topology preserved by construction (C-O stays in bonding range).
- Per-hydroxyl rollback: if the move introduces a NEW non-bonded clash, revert
  that hydroxyl only.
- Pure geometry, deterministic, no force field, no env flag inside the module
  (the pipeline wrapper decides).
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from delfin._fix_sp3_h_tetrahedrality import (
    _COV_RADII,
    _VDW_RADII,
    _build_geometric_adjacency,
    _compute_angle_deg,
    _format_xyz,
    _is_metal_sym,
    _parse_xyz,
)

_IDEAL_COH_DEG = 108.5          # alcohol/phenol C-O-H
_CO_ALCOHOL = 1.43              # sp3 C-OH
_CO_PHENOL = 1.36               # aromatic C-OH
_CO_CARBOXYL = 1.31            # HO-C(=O)
_OH_DEFAULT = 0.97             # used only if an embedded O-H is degenerate


def _rings_up_to_6(adj: List[List[int]], syms: List[str]) -> List[List[int]]:
    """All simple rings of size 5-6 over heavy atoms (for aromaticity test)."""
    n = len(syms)
    R: List[List[int]] = []
    seen = set()

    def dfs(start: int, cur: List[int], visited: set):
        last = cur[-1]
        for nb in adj[last]:
            if syms[nb] == "H":
                continue
            if nb == start and 5 <= len(cur) <= 6:
                k = frozenset(cur)
                if k not in seen:
                    seen.add(k); R.append(list(cur))
            elif nb not in visited and len(cur) < 6:
                dfs(start, cur + [nb], visited | {nb})

    for s in range(n):
        if syms[s] != "H" and not _is_metal_sym(syms[s]):
            dfs(s, [s], {s})
    return R


def _carbon_is_aromatic(c_idx: int, pts: np.ndarray, rings: List[List[int]]) -> bool:
    """True if C sits in a planar 5-6 ring (max out-of-plane < 0.20 A)."""
    for rg in rings:
        if c_idx not in rg:
            continue
        ring_pts = pts[rg]
        centroid = ring_pts.mean(axis=0)
        _, _, vt = np.linalg.svd(ring_pts - centroid)
        oop = float(np.abs((ring_pts - centroid) @ vt[2]).max())
        if oop < 0.20:
            return True
    return False


def _carbon_is_carboxyl(c_idx: int, o_idx: int, adj: List[List[int]],
                        syms: List[str], pts: np.ndarray) -> bool:
    """True if the parent C carries a second O at carbonyl distance (C=O <1.30)."""
    for nb in adj[c_idx]:
        if nb == o_idx or syms[nb] != "O":
            continue
        if float(np.linalg.norm(pts[c_idx] - pts[nb])) < 1.30:
            return True
    return False


def _move_introduces_clash(new_pos: np.ndarray, moved_idx: int,
                           keep: set, syms: List[str], pts: np.ndarray,
                           factor: float = 0.70) -> bool:
    """True iff ``moved_idx`` at ``new_pos`` overlaps any atom not in ``keep``
    (bonded partners) below ``factor``*Sum_vdw, ignoring normal bonded range."""
    ri_v = _VDW_RADII.get(syms[moved_idx], 1.7)
    ri_c = _COV_RADII.get(syms[moved_idx], 1.5)
    for j in range(len(syms)):
        if j == moved_idx or j in keep:
            continue
        d = float(np.linalg.norm(new_pos - pts[j]))
        rj_c = _COV_RADII.get(syms[j], 1.5)
        if 0.85 * (ri_c + rj_c) <= d <= 1.15 * (ri_c + rj_c):
            continue  # bonded range
        if d < factor * (ri_v + _VDW_RADII.get(syms[j], 1.7)):
            return True
    return False


def fix_hydroxyl_geometry(xyz_str: str, mol=None, *,
                          angle_max_deg: float = 115.0,
                          angle_min_deg: float = 100.0,
                          length_tol: float = 0.05,
                          do_length: bool = False,
                          ) -> Tuple[str, Dict]:
    """Correct pendant-hydroxyl C-O-H angle (and, if ``do_length``, the C-O bond
    length).  Returns ``(new_xyz, report)``; report has ``angle_fixed``,
    ``length_fixed``, ``topology_preserved`` (always True by construction) and
    ``n_hydroxyl``.  No-op (returns input) on parse failure / nothing qualifies.

    ``do_length`` defaults OFF: data (2026-06-23, V2R) showed the ANGLE
    correction is a clean win (C-O-H >115 deg 46%->2%, median 113->108.5, no new
    real clashes), but the C-O-length correction was marginal/slightly negative
    because geometric phenol/carboxyl classification is unreliable (and much of
    the "C-O too short" census signal was phenols at 1.36, which are correct).
    The length leg stays available for a future mol-based hybridisation call."""
    report = {"angle_fixed": 0, "length_fixed": 0, "n_hydroxyl": 0,
              "topology_preserved": True}
    try:
        syms, pts, lines = _parse_xyz(xyz_str)
    except Exception:
        return xyz_str, report
    n = len(syms)
    if n < 3 or pts.shape[0] != n:
        return xyz_str, report
    pts = pts.copy()
    adj = _build_geometric_adjacency(syms, pts)
    rings = _rings_up_to_6(adj, syms)
    changed = False
    for o in range(n):
        if syms[o] != "O":
            continue
        nbrs = adj[o]
        hs = [k for k in nbrs if syms[k] == "H"]
        cs = [k for k in nbrs if syms[k] == "C"]
        heavy = [k for k in nbrs if syms[k] != "H"]
        metal_nb = any(_is_metal_sym(syms[k]) for k in nbrs)
        # pendant hydroxyl: exactly 1 H, exactly 1 heavy (a C), no metal contact
        if metal_nb or len(hs) != 1 or len(heavy) != 1 or len(cs) != 1:
            continue
        report["n_hydroxyl"] += 1
        h = hs[0]; c = cs[0]
        # --- target C-O length (hybridisation-aware) ---
        if _carbon_is_carboxyl(c, o, adj, syms, pts):
            co_ideal = _CO_CARBOXYL
        elif _carbon_is_aromatic(c, pts, rings):
            co_ideal = _CO_PHENOL
        else:
            co_ideal = _CO_ALCOHOL
        co_vec = pts[o] - pts[c]
        co_len = float(np.linalg.norm(co_vec))
        if co_len < 1e-6:
            continue
        # 1-3 exclusion sets: a moved atom legitimately sits ~1.9 A from its
        # 1-3 partner (e.g. H···C across O), which is NOT a clash.  Exclude the
        # moved atom's bonded + 1-3 neighbours from the clash test (mirrors the
        # sibling-exclusion in _fix_sp3_h_tetrahedrality._h_introduces_clash).
        keep_h = {o, c}                       # H is 1-2 to O, 1-3 to C
        keep_o = {c, h} | set(adj[c])         # O is 1-2 to C/H, 1-3 to C's nbrs
        # --- 1) fix C-O length: rigid {O,H} translation along C->O axis ---
        if do_length and abs(co_len - co_ideal) > length_tol:
            new_o = pts[c] + co_vec / co_len * co_ideal
            trans = new_o - pts[o]
            new_h = pts[h] + trans
            if (not _move_introduces_clash(new_o, o, keep_o, syms, pts)
                    and not _move_introduces_clash(new_h, h, keep_h, syms, pts)):
                pts[o] = new_o; pts[h] = new_h; changed = True
                report["length_fixed"] += 1
        # --- 2) fix C-O-H angle -> 108.5 deg (rotate H about O, keep O-H len) ---
        ang = _compute_angle_deg(pts[c], pts[o], pts[h])
        if math.isnan(ang) or not (ang > angle_max_deg or ang < angle_min_deg):
            continue
        oh = pts[h] - pts[o]
        oh_len = float(np.linalg.norm(oh))
        if oh_len < 1e-6:
            oh_len = _OH_DEFAULT
        u = pts[c] - pts[o]; un = float(np.linalg.norm(u))
        if un < 1e-6:
            continue
        u = u / un
        # perpendicular component of O->H w.r.t. u (in the C-O-H plane)
        w = oh - np.dot(oh, u) * u
        wn = float(np.linalg.norm(w))
        if wn < 1e-6:
            continue
        w = w / wn
        th = math.radians(_IDEAL_COH_DEG)
        new_dir = math.cos(th) * u + math.sin(th) * w
        new_h = pts[o] + oh_len * new_dir
        if not _move_introduces_clash(new_h, h, keep_h, syms, pts):
            pts[h] = new_h; changed = True
            report["angle_fixed"] += 1
    if not changed:
        return xyz_str, report
    try:
        return _format_xyz(lines, syms, pts), report
    except Exception:
        return xyz_str, report
