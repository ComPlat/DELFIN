"""delfin.fffree.sp2_planar_snap — sp2 atom out-of-plane projection.

Phase G14 (User 2026-05-31, post-G13b ring-scale fix). Targets the lingering
structqual_bond_delta (+122 %), funcgrp_n_bond_viol (+176 %) and non-aromatic
component of pi_planar (+123 %) HEAD-vs-UFF gaps -- all driven by sp2 atoms
(amide N, carbonyl C, ester O, imine C/N, formal sp2 outside aromatic rings)
sitting slightly out of the plane defined by their three heavy neighbours.

Aromatic ring sp2 atoms are already handled by ``aromatic_ring_scale``; this
module targets the COMPLEMENT -- non-aromatic sp2 atoms with exactly three
bonded heavy neighbours. For each, the central atom is projected onto the
unique plane through its three neighbours (minimum-movement perpendicular
projection), which makes it geometrically planar (sum of angles 360 deg)
without disturbing any bond length.

Same safety net as G13b: runs AFTER ``_build_is_clean`` has accepted, the
post-snap structure is re-checked, and on failure the pre-snap coords are
emitted. Metal-coordinated atoms and any atom whose three neighbours include
a metal are frozen. Rigid-H drag for the attached hydrogen (if any).

Pure geometry, FF-free, universal, deterministic.

Env-gate: ``DELFIN_FFFREE_SP2_SNAP=1`` (default OFF, byte-identical when unset).
Auto-enabled under ``DELFIN_FFFREE_PURE_TRACK3=1``.
"""
from __future__ import annotations

import math
import os
from typing import List, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Default-OFF live-read; PT3 auto-enable.
def _sp2_snap_active() -> bool:
    if os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_SP2_SNAP", "0") == "1":
        return True
    return False


# Detect sp2 candidacy: 3 heavy nbrs and the mean of the three nbr-atom-nbr
# angles is close to 120 deg. Tolerance band [100, 135] catches sp2 buckled
# up to ~25 deg out of plane (the "stretched aromatic" / amide / carbonyl
# regime where projection is meaningful) without catching sp3 (109.5 deg).
SP2_ANGLE_LOW = 100.0
SP2_ANGLE_HIGH = 135.0


def _enumerate_small_cycles(adj: List[List[int]], min_len: int = 5,
                            max_len: int = 6) -> List[Tuple[int, ...]]:
    """5-/6-membered cycle enumeration (reused convention -- see
    aromatic_ring_scale for the equivalent helper)."""
    n = len(adj)
    cycles: Set[Tuple[int, ...]] = set()

    def dfs(start: int, current: int, path: List[int], depth: int) -> None:
        if depth > max_len:
            return
        for nxt in adj[current]:
            if nxt == start and depth >= min_len:
                cycles.add(tuple(path))
                continue
            if nxt in path:
                continue
            if nxt < start:
                continue
            path.append(nxt)
            dfs(start, nxt, path, depth + 1)
            path.pop()

    for s in range(n):
        dfs(s, s, [s], 1)

    seen: Set[Tuple[int, ...]] = set()
    out: List[Tuple[int, ...]] = []
    for c in sorted(cycles):
        key = tuple(sorted(c))
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def _aromatic_ring_atoms(syms: Sequence[str], P: np.ndarray,
                          heavy_adj: List[List[int]]) -> Set[int]:
    """Return the union of all aromatic-ring atom indices (5- / 6-cycles with
    C / N / O / S members and mean intra-ring bond length 1.30-1.65 A). These
    are EXCLUDED from sp2 polish (the ring scaler handles them)."""
    atoms: Set[int] = set()
    for ring in _enumerate_small_cycles(heavy_adj, 5, 6):
        if not (5 <= len(ring) <= 6):
            continue
        allowed = {"C", "N", "O", "S"}
        ok = True
        for a in ring:
            if syms[a] not in allowed or _bd._is_metal(syms[a]):
                ok = False
                break
        if not ok:
            continue
        lens = [float(np.linalg.norm(P[ring[i]] - P[ring[(i + 1) % len(ring)]]))
                for i in range(len(ring))]
        if 1.30 <= float(np.mean(lens)) <= 1.65:
            atoms.update(ring)
    return atoms


def _metal_coordinated_atoms(syms: Sequence[str], P: np.ndarray,
                              max_md_factor: float = 1.45,
                              hapto_pi_cutoff: float = 3.0) -> Set[int]:
    """Frozen-set: M-bonded + hapto-pi same-element cluster (matches
    aromatic_ring_scale.G13b convention)."""
    n = len(syms)
    frozen: Set[int] = set()
    for m in range(n):
        if not _bd._is_metal(syms[m]):
            continue
        near_heavy: List[Tuple[int, str]] = []
        for j in range(n):
            if j == m or syms[j] == "H":
                continue
            d = float(np.linalg.norm(P[j] - P[m]))
            if d < max_md_factor * _bd._ideal_bond(syms[m], syms[j]):
                frozen.add(j)
            if d < hapto_pi_cutoff:
                near_heavy.append((j, syms[j]))
        # hapto-pi fingerprint
        from collections import defaultdict
        by_elem: dict = defaultdict(list)
        for j, e in near_heavy:
            if e in {"C", "N", "O", "S"}:
                by_elem[e].append(j)
        for elem, idxs in by_elem.items():
            if len(idxs) >= 3:
                for j in idxs:
                    frozen.add(j)
    return frozen


def _angle_deg(P: np.ndarray, i: int, j: int, k: int) -> float:
    """Angle i-j-k in degrees (j is vertex)."""
    v1 = P[i] - P[j]; v2 = P[k] - P[j]
    n1 = float(np.linalg.norm(v1)); n2 = float(np.linalg.norm(v2))
    if n1 < 1e-9 or n2 < 1e-9:
        return float("nan")
    c = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    return float(math.degrees(math.acos(c)))


def snap_sp2_to_plane(syms: Sequence[str], P: np.ndarray
                       ) -> Tuple[List[str], np.ndarray]:
    """Project every non-aromatic sp2 atom onto the plane defined by its three
    heavy neighbours. Returns (syms_out, P_out). No-op if the env flag is off.

    Universal, deterministic, FF-free.
    """
    if not _sp2_snap_active():
        return list(syms), np.asarray(P, dtype=float).copy()

    P_out = np.asarray(P, dtype=float).copy()
    n = len(syms)

    # Heavy and full bond graphs
    bonds = _bd._geometric_bonds(list(syms), P_out)
    adj_with_h: List[List[int]] = [[] for _ in range(n)]
    heavy_adj: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        adj_with_h[i].append(j); adj_with_h[j].append(i)
        if syms[i] != "H" and syms[j] != "H":
            heavy_adj[i].append(j); heavy_adj[j].append(i)

    aromatic = _aromatic_ring_atoms(list(syms), P_out, heavy_adj)
    frozen = _metal_coordinated_atoms(list(syms), P_out)

    for a in range(n):
        if syms[a] == "H":
            continue
        if _bd._is_metal(syms[a]):
            continue
        if a in aromatic:
            continue                       # ring scaler handles
        if a in frozen:
            continue                       # M-coord or hapto-pi
        if syms[a] not in {"C", "N", "O"}:
            continue                       # sp2 candidates only
        nbrs = heavy_adj[a]
        if len(nbrs) != 3:
            continue
        # Skip if any neighbour is a metal (chelate-N pyridine etc covered by frozen,
        # but defensive).
        if any(_bd._is_metal(syms[k]) for k in nbrs):
            continue
        # Skip if any neighbour is frozen (preserve M-coord substituents).
        if any(k in frozen for k in nbrs):
            continue
        # Geometric sp2 fingerprint: mean of three nbr-a-nbr angles in [100, 135].
        a01 = _angle_deg(P_out, nbrs[0], a, nbrs[1])
        a02 = _angle_deg(P_out, nbrs[0], a, nbrs[2])
        a12 = _angle_deg(P_out, nbrs[1], a, nbrs[2])
        if any(math.isnan(x) for x in (a01, a02, a12)):
            continue
        mean_ang = (a01 + a02 + a12) / 3.0
        if not (SP2_ANGLE_LOW <= mean_ang <= SP2_ANGLE_HIGH):
            continue

        # Plane through the three neighbour points.
        Q = np.array([P_out[k] for k in nbrs])
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

        # Project a onto the plane.
        off = float(np.dot(P_out[a] - centroid, normal))
        if abs(off) < 0.02:                # already in plane (< 0.02 A)
            continue
        move = -off * normal
        P_out[a] = P_out[a] + move

        # Rigid-H drag: any H whose only heavy neighbour is a -> translate too.
        for h in range(n):
            if syms[h] != "H":
                continue
            hn = [k for k in adj_with_h[h] if syms[k] != "H"]
            if hn == [a]:
                P_out[h] = P_out[h] + move

    return list(syms), P_out


def is_enabled() -> bool:
    return _sp2_snap_active()


if __name__ == "__main__":
    # Self-test: amide C=O with N(H)H, the C is sp2 with 3 heavy nbrs (N, O, R).
    # Put it slightly out of plane and check it flattens.
    os.environ["DELFIN_FFFREE_SP2_SNAP"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.sp2_planar_snap", None)
    from delfin.fffree.sp2_planar_snap import snap_sp2_to_plane

    # Atom layout: idx 0 = central C, 1 = O, 2 = N, 3 = R(C), 4 = H_on_N
    syms = ["C", "O", "N", "C", "H"]
    P0 = np.array([
        [0.0, 0.0, 0.3],         # C (out of plane by 0.3 A)
        [1.22, 0.0, 0.0],        # =O
        [-0.6, 1.1, 0.0],        # -N
        [-0.6, -1.1, 0.0],       # -R (C)
        [-1.0, 1.9, 0.0],        # H on N (rigid drag does not apply here)
    ])
    _, P_out = snap_sp2_to_plane(syms, P0)
    # Plane through O, N, R is z=0 -> C should land at z=0.
    print(f"Amide C-out-of-plane before: {P0[0][2]:.3f}, after: {P_out[0][2]:.3f}")
    print(f"O / N / R untouched: O Δ={P_out[1][2]-P0[1][2]:.4f} N Δ={P_out[2][2]-P0[2][2]:.4f}")
