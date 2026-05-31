"""delfin.fffree.aromatic_ring_scale — rigid uniform ring scaling.

Phase G13 (User 2026-05-31, post-G12 negative result + Mogul ringfix forensik).

After the autonomous Mogul L1 cyclopropane blindspot fix
([[feedback_mogul_ringfix_clean_ranking_2026_05_31]]) the true TOP defect in
BEST-3 sample 500 became C sp2 bond 1.81 vs COD ideal 1.40 (+29 %, 32 files
in the PPh3 / PR3-aryl class). Direct experiment showed that ``refine.refine``
under PURE_TRACK3 cannot fix this (KIHVIH: 17 stretched > 1.50 Å, mean
1.559 Å -> after full refine pass, 15 stretched, mean 1.544 Å; only 2 bonds
compressed). Root cause: pairwise accept-if-better gate; compressing one ring
bond distorts the rest of the ring, the trial loss increases, refine rejects
([[feedback_refine_aromatic_pairwise_plateau]]).

Fix: detect aromatic ring geometrically, scale it uniformly about its centroid
to bring the mean intra-ring bond to 1.40 A. Rigid as a single body, so the
ring's PLANARITY + RELATIVE ORIENTATION + RELATIVE ANGLES are preserved (only
the size changes). Ring-attached H ride rigidly with their parent. Rings
containing a metal-coordinated atom are FROZEN -- preserves M-D invariant +
chelate-attachment geometry.

The G12 lesson [[feedback_g12_aromsnap_failed_isomer_selection]]: a per-atom
geometry modifier wired BEFORE `_build_is_clean` perturbs the acceptance
basket enough that the primary isomer gets rejected and the FF-free path falls
through to alternative isomers with up to 17 A atom shifts. G13 avoids this
by running AFTER `_build_is_clean` has accepted and falling back to pre-scale
coords if the post-scale structure would now fail the clean check.

Pure geometry, FF-free, universal across all aromatic ring sizes 5-6.
Deterministic. Env-gate: ``DELFIN_FFFREE_RING_SCALE=1`` (default OFF,
byte-identical when unset). Auto-enabled under PURE_TRACK3.
"""
from __future__ import annotations

import os
from typing import Iterable, List, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


# Target ring C-C / heteroaromatic bond (aromatic centre). Empirical COD median
# for in-ring sp2 C-C: 1.397 A; for sp2 C-N: 1.34 A; for sp2 C-O: 1.36 A. We
# scale toward 1.40 A as the canonical aromatic-class target -- universal,
# element-agnostic, and within 3 % of all common aromatic heteroatom bonds.
TARGET_BOND = 1.40
SCALE_CLAMP_MIN = 0.85   # max compression (refuse if needs > 15 % shrink)
SCALE_CLAMP_MAX = 1.05   # max expansion (refuse if needs > 5 % grow)
# Skip scaling when mean ring bond is already within this band of target.
TOLERANCE = 0.05         # +/- 0.05 A = within ~3.5 % of 1.40


def _ring_scale_active() -> bool:
    """Live-read env flags (call-time, never import-time -- pytest fixtures
    can flip them between processes)."""
    if os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_RING_SCALE", "0") == "1":
        return True
    return False


def _enumerate_small_cycles(adj: List[List[int]], min_len: int = 5,
                            max_len: int = 6) -> List[Tuple[int, ...]]:
    """Find all simple cycles of length [min_len, max_len] in an undirected
    heavy-atom bond graph. Universal DFS, deterministic ordering."""
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
            if nxt < start:        # avoid duplicate traversal
                continue
            path.append(nxt)
            dfs(start, nxt, path, depth + 1)
            path.pop()

    for s in range(n):
        dfs(s, s, [s], 1)

    # De-dup by sorted-tuple identity (rotational duplicates collapse).
    seen: Set[Tuple[int, ...]] = set()
    out: List[Tuple[int, ...]] = []
    for c in sorted(cycles):
        key = tuple(sorted(c))
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def _is_aromatic_ring(ring: Sequence[int], syms: Sequence[str],
                      P: np.ndarray) -> bool:
    """Aromatic test: 5-/6-membered, all C / N / O / S, no metal, AND mean
    closed-loop nearest-neighbour bond length in the aromatic range
    1.30-1.65 A (admits both kekulised single bonds and aromatic doubles,
    plus stretched aromatic up to 1.65 A which is exactly what we want to
    rescue)."""
    n = len(ring)
    if n < 5 or n > 6:
        return False
    allowed = {"C", "N", "O", "S"}
    for a in ring:
        if syms[a] not in allowed:
            return False
        if _bd._is_metal(syms[a]):
            return False
    lens = []
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        lens.append(float(np.linalg.norm(P[a] - P[b])))
    mean_len = float(np.mean(lens))
    return 1.30 <= mean_len <= 1.65


def _metal_coordinated_atoms(syms: Sequence[str], P: np.ndarray,
                              max_md_factor: float = 1.45) -> Set[int]:
    """Indices of atoms within max_md_factor * ideal_bond of any metal.
    These rings get frozen (preserves M-D invariant + chelate attachment)."""
    n = len(syms)
    frozen: Set[int] = set()
    for m in range(n):
        if not _bd._is_metal(syms[m]):
            continue
        for j in range(n):
            if j == m or syms[j] == "H":
                continue
            d = float(np.linalg.norm(P[j] - P[m]))
            if d < max_md_factor * _bd._ideal_bond(syms[m], syms[j]):
                frozen.add(j)
    return frozen


def _ring_attached_hs(ring_set: Set[int], syms: Sequence[str],
                      adj_with_h: List[List[int]]) -> List[Tuple[int, int]]:
    """List of (h_idx, parent_ring_atom) for H whose only heavy neighbour is
    a ring atom. These ride rigidly with their parent's net displacement
    (preserves X-H bond length to first order)."""
    pairs: List[Tuple[int, int]] = []
    for h in range(len(syms)):
        if syms[h] != "H":
            continue
        heavy_nbrs = [k for k in adj_with_h[h] if syms[k] != "H"]
        if len(heavy_nbrs) != 1:
            continue
        parent = heavy_nbrs[0]
        if parent in ring_set:
            pairs.append((h, parent))
    return pairs


def scale_aromatic_rings(syms: Sequence[str], P: np.ndarray,
                          target_bond: float = TARGET_BOND,
                          clamp_min: float = SCALE_CLAMP_MIN,
                          clamp_max: float = SCALE_CLAMP_MAX,
                          tolerance: float = TOLERANCE
                          ) -> Tuple[List[str], np.ndarray]:
    """Rigid uniform scaling of every detected aromatic ring about its
    centroid, with rigid-H drag and metal-coordination freeze.

    Universal across all TMC classes. Deterministic. FF-free. No-op if the
    env flag is unset (returns the input unchanged).

    Returns (syms_out, P_out) where P_out is a fresh copy.
    """
    if not _ring_scale_active():
        return list(syms), np.asarray(P, dtype=float).copy()

    P_out = np.asarray(P, dtype=float).copy()
    n = len(syms)

    # Build heavy-atom bond graph for ring detection; include H bonds
    # separately for the rigid-drag map.
    bonds = _bd._geometric_bonds(list(syms), P_out)
    adj_with_h: List[List[int]] = [[] for _ in range(n)]
    heavy_adj: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        adj_with_h[i].append(j); adj_with_h[j].append(i)
        if syms[i] != "H" and syms[j] != "H":
            heavy_adj[i].append(j); heavy_adj[j].append(i)

    cycles = _enumerate_small_cycles(heavy_adj, min_len=5, max_len=6)
    if not cycles:
        return list(syms), P_out

    frozen = _metal_coordinated_atoms(list(syms), P_out)

    for ring in cycles:
        if not _is_aromatic_ring(ring, syms, P_out):
            continue
        ring_set = set(ring)
        # Freeze: any ring atom that is coordinated to a metal -> the whole
        # ring is FROZEN (no scale). This protects pi-bound rings and rings
        # that carry a coordinating heteroatom (e.g. pyridine N coordinated
        # to M); their bond lengths come from the chelate / hapto template.
        if ring_set & frozen:
            continue
        ring_indices = list(ring)
        pts = P_out[ring_indices]
        centroid = pts.mean(axis=0)

        # Compute current mean intra-ring bond length.
        m = len(ring_indices)
        lens = []
        for i in range(m):
            a = ring_indices[i]; b = ring_indices[(i + 1) % m]
            lens.append(float(np.linalg.norm(P_out[a] - P_out[b])))
        mean_len = float(np.mean(lens))
        # Already on target?
        if abs(mean_len - target_bond) <= tolerance:
            continue
        # Scale factor; clamp.
        s = target_bond / mean_len
        s = max(clamp_min, min(clamp_max, s))
        if abs(s - 1.0) < 1e-4:
            continue

        # Uniform scaling about centroid. Frozen-set protection applied
        # below per atom (any one frozen blocks the ring above; so here all
        # ring atoms are free).
        moves = {}
        for idx in ring_indices:
            new_pos = centroid + s * (P_out[idx] - centroid)
            moves[idx] = new_pos - P_out[idx]
            P_out[idx] = new_pos
        # Rigid-H drag.
        h_pairs = _ring_attached_hs(ring_set, syms, adj_with_h)
        for h, parent in h_pairs:
            P_out[h] = P_out[h] + moves[parent]

    return list(syms), P_out


def is_enabled() -> bool:
    """Whether the env gate is active in this process."""
    return _ring_scale_active()


if __name__ == "__main__":
    # Self-test: stretched benzene 1.65 A -> rescale to 1.40 A.
    os.environ["DELFIN_FFFREE_RING_SCALE"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.aromatic_ring_scale", None)
    from delfin.fffree.aromatic_ring_scale import scale_aromatic_rings

    syms = ["C"] * 6 + ["H"] * 6
    # Hexagon with C-C 1.65 (stretched) in xy, H 1.08 outward
    R_ring = 1.65  # circumradius = bond for regular hexagon
    R_h = R_ring + 1.08
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    P0 = np.array(
        [[R_ring * np.cos(t), R_ring * np.sin(t), 0.0] for t in th]
        + [[R_h * np.cos(t), R_h * np.sin(t), 0.0] for t in th]
    )
    # Measure
    def _mean_cc(P):
        m = 0
        ls = []
        for i in range(6):
            ls.append(float(np.linalg.norm(P[i] - P[(i + 1) % 6])))
        return float(np.mean(ls))

    print(f"Stretched benzene: mean C-C before = {_mean_cc(P0):.3f}")
    _, P_out = scale_aromatic_rings(syms, P0)
    print(f"After ring_scale:      mean C-C after  = {_mean_cc(P_out):.3f}  (target 1.40)")
    # Check C-H bond preservation
    for k in range(6):
        d_b = float(np.linalg.norm(P0[k + 6] - P0[k]))
        d_a = float(np.linalg.norm(P_out[k + 6] - P_out[k]))
        print(f"  C{k}-H{k}: before={d_b:.3f}  after={d_a:.3f}")
