"""delfin.fffree.aromatic_snap — geometric aromatic-ring planarity snap.

Phase G12 (User 2026-05-31): Mogul aggregator on fb1ae9a-PT3 voll-pool flagged
the TOP fragment-severity blindspot as ``C sp2`` bonds stretched to 1.835 Å vs
COD median 1.397 Å (+31 %, 41/500 files, total severity 26 307). Forensic
inspection of one example (D-KIHVIH_4d_Pd_CN2.xyz) confirmed that aromatic ring
carbons in PPh3/PR3-aryl ligands sit ~1.96 Å apart on the long axis after
chelate placement + refine — the rings are slightly buckled, and the refine
loop's bond_distort band rejects the corrective moves because they break a
balancing constraint.

The legacy/UFF pipeline already has ``_snap_aromatic_rings_to_plane`` in
``smiles_converter.py:21750`` (called 5× in the UFF path) but it was never
wired into the FF-free path — a clear missed link. Rather than depending on a
mol_template (which the FF-free path does not always have at this stage), this
module performs the same SVD plane-fit purely from (syms, P) using geometric
ring detection.

  Algorithm:
    1. Geometric bond graph: bd._geometric_bonds on heavy atoms
    2. Cycle detection (DFS), keep cycles of size 5-6
    3. Aromatic-ring filter: all members C / N / O / S and the mean intra-ring
       bond length ~1.32-1.50 Å (aromatic / kekulised). Skip metallacycles.
    4. SVD best-fit plane → project each ring atom onto the plane (minimal-
       movement perpendicular component).
    5. Ring-attached H rides rigidly (preserve X-H length): a hydrogen whose
       only heavy neighbour is a ring atom is translated by the same vector
       its parent moved.
    6. M-D invariant guard: ring atoms that are also coordinated to a metal
       (donor atoms) are FROZEN (no move) — preserves M-D bond length.

Pure geometry, FF-free, universal across all TMC classes, no force field, no
SMILES-specific shortcuts. Deterministic (no RNG).

Env-gate: ``DELFIN_FFFREE_AROMSNAP=1`` (default OFF, byte-identical when unset).
Auto-enabled under ``DELFIN_FFFREE_PURE_TRACK3=1`` (race-strategy: bring the
aromatic-ring planarity to legacy parity inside the FF-free path).
"""
from __future__ import annotations

import os
from typing import List, Sequence, Set, Tuple

import numpy as np

import delfin._bond_decollapse as _bd


def _aromsnap_active() -> bool:
    """Live-read the env flags so toggling between processes/tests works
    without re-importing.  Universal pattern: env-gates always read at call
    time, never at import time (otherwise pytest fixtures can't flip them).
    """
    if os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_AROMSNAP", "0") == "1":
        return True
    return False


def _enumerate_small_cycles(adj: List[List[int]], max_len: int = 6,
                            min_len: int = 5) -> List[Tuple[int, ...]]:
    """Find all simple cycles of length [min_len, max_len] in undirected graph.

    Universal DFS-based cycle enumeration. Deterministic order (atom-index sort).
    Returns each cycle once, as a tuple of atom indices sorted-rotation-normalised.
    """
    n = len(adj)
    cycles: Set[Tuple[int, ...]] = set()

    def dfs(start: int, current: int, path: List[int], depth: int) -> None:
        if depth > max_len:
            return
        for nxt in adj[current]:
            if nxt == start and depth >= min_len:
                # found cycle
                normalised = tuple(sorted(path))
                if normalised not in cycles:
                    # store ordered path (rotation canonicalised)
                    cycles.add(tuple(path))
                continue
            if nxt in path:
                continue
            if nxt < start:  # avoid duplicate traversal
                continue
            path.append(nxt)
            dfs(start, nxt, path, depth + 1)
            path.pop()

    for s in range(n):
        dfs(s, s, [s], 1)

    # De-duplicate by sorted-tuple identity (handles rotational duplicates)
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
    """Geometric aromatic-ring test.

    Conditions (all must hold):
    - All ring members are C, N, O or S (heteroaromatic admissible).
    - Mean intra-ring nearest-neighbour bond length ∈ [1.30, 1.50] Å (aromatic /
      kekulised C-C / C-N / C-O / C-S range).
    - No metal among ring atoms (skip metallacycles — they have their own
      planarity treatment in assemble_complex).
    """
    n = len(ring)
    if n < 5 or n > 6:
        return False
    aromatic_elements = {"C", "N", "O", "S"}
    for a in ring:
        if syms[a] not in aromatic_elements:
            return False
        if _bd._is_metal(syms[a]):
            return False
    # Intra-ring nearest-neighbour bond lengths (closed cycle).
    lens = []
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        lens.append(float(np.linalg.norm(P[a] - P[b])))
    mean_len = float(np.mean(lens))
    return 1.30 <= mean_len <= 1.55


def _ring_attached_hs(ring: Sequence[int], syms: Sequence[str],
                       P: np.ndarray, adj: List[List[int]]) -> List[Tuple[int, int]]:
    """Return list of (h_idx, parent_ring_atom) for each H whose only heavy
    neighbour is a ring atom. These ride rigidly with their parent's move.
    """
    ring_set = set(ring)
    pairs: List[Tuple[int, int]] = []
    for h in range(len(syms)):
        if syms[h] != "H":
            continue
        heavy_nbrs = [k for k in adj[h] if syms[k] != "H"]
        if len(heavy_nbrs) != 1:
            continue
        parent = heavy_nbrs[0]
        if parent in ring_set:
            pairs.append((h, parent))
    return pairs


def _metal_coordinated_atoms(syms: Sequence[str], P: np.ndarray,
                              max_md_factor: float = 1.45) -> Set[int]:
    """Return indices of atoms within max_md_factor * ideal_bond of any metal.

    These are FROZEN during the aromatic snap to preserve M-D invariant
    (per [[feedback_md_invariant]]: M-D rotations within 0.05 Å are mandatory).
    """
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


def snap_aromatic_rings(syms: Sequence[str], P: np.ndarray,
                         rms_threshold: float = 0.05) -> Tuple[List[str], np.ndarray]:
    """Snap all detected aromatic rings to their best-fit plane.

    Returns (syms_out, P_out). Operates on a copy of P; original is untouched.
    Ring-attached H atoms ride rigidly with their parent's projection vector
    (preserves X-H bond length). Metal-coordinated atoms are frozen (preserves
    M-D invariant ≤ 0.05 Å). No-op if no aromatic ring is found.

    Universal across all TMC classes. Deterministic. FF-free.
    """
    if not _aromsnap_active():
        return list(syms), np.asarray(P, dtype=float).copy()
    P_out = np.asarray(P, dtype=float).copy()
    n = len(syms)

    # Build heavy-atom bond graph (ignore H for cycle search; H attached via adj
    # for ride-along H tracking).
    bonds = _bd._geometric_bonds(list(syms), P_out)
    adj: List[List[int]] = [[] for _ in range(n)]
    heavy_adj: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        adj[i].append(j); adj[j].append(i)
        if syms[i] != "H" and syms[j] != "H":
            heavy_adj[i].append(j); heavy_adj[j].append(i)

    cycles = _enumerate_small_cycles(heavy_adj, max_len=6, min_len=5)
    if not cycles:
        return list(syms), P_out

    frozen = _metal_coordinated_atoms(list(syms), P_out)

    any_changed = False
    for ring in cycles:
        if not _is_aromatic_ring(ring, syms, P_out):
            continue
        ring_indices = list(ring)
        pts = P_out[ring_indices]
        centroid = pts.mean(axis=0)
        centered = pts - centroid
        try:
            _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        normal = Vt[-1]
        nn = float(np.linalg.norm(normal))
        if nn < 1e-9:
            continue
        normal = normal / nn

        # Out-of-plane component per ring atom
        offsets = centered @ normal     # (m,)
        rms = float(np.sqrt(float((offsets * offsets).mean())))
        if rms < rms_threshold:
            continue  # already planar enough

        # Project: subtract out-of-plane component. Frozen atoms move 0.
        h_pairs = _ring_attached_hs(ring_indices, syms, P_out, adj)
        moves = {}
        for k, idx in enumerate(ring_indices):
            if idx in frozen:
                moves[idx] = np.zeros(3)
                continue
            moves[idx] = -offsets[k] * normal
        # Apply ring projection
        for idx, mv in moves.items():
            P_out[idx] = P_out[idx] + mv
        # Ride-along H: copy parent's move (rigid drag preserves X-H length)
        for h, parent in h_pairs:
            if parent in moves:
                P_out[h] = P_out[h] + moves[parent]
        any_changed = True

    return list(syms), P_out


def is_enabled() -> bool:
    """True iff AROMSNAP is active in the current process environment."""
    return _aromsnap_active()


if __name__ == "__main__":
    # Self-test: build a buckled benzene (6 C + 6 H), snap to plane.
    os.environ["DELFIN_FFFREE_AROMSNAP"] = "1"
    import importlib, sys
    sys.modules.pop("delfin.fffree.aromatic_snap", None)
    from delfin.fffree.aromatic_snap import snap_aromatic_rings

    syms = ["C"] * 6 + ["H"] * 6
    # Hexagon in xy + perturbation along z (buckle)
    th = np.linspace(0, 2 * np.pi, 7)[:-1]
    R_ring = 1.4
    R_h = 2.5
    P0 = np.array(
        [[R_ring * np.cos(t), R_ring * np.sin(t),
          0.3 * (i % 2)] for i, t in enumerate(th)] +
        [[R_h * np.cos(t), R_h * np.sin(t), 0.0] for t in th]
    )
    syms_out, P_out = snap_aromatic_rings(syms, P0)
    # Check planarity
    pts_ring = P_out[:6]
    centroid = pts_ring.mean(axis=0)
    centered = pts_ring - centroid
    _, _, Vt = np.linalg.svd(centered, full_matrices=False)
    normal = Vt[-1]
    offsets = centered @ normal
    rms_after = float(np.sqrt(float((offsets * offsets).mean())))
    pts_ring_before = P0[:6]
    centered_b = pts_ring_before - pts_ring_before.mean(axis=0)
    _, _, Vt_b = np.linalg.svd(centered_b, full_matrices=False)
    offsets_b = centered_b @ Vt_b[-1]
    rms_before = float(np.sqrt(float((offsets_b * offsets_b).mean())))
    print(f"Benzene snap: rms_before={rms_before:.4f} -> rms_after={rms_after:.4f}")
    # Check H followed parent
    for k in range(6):
        d = float(np.linalg.norm(P_out[k + 6] - P_out[k]))
        print(f"  C{k}-H{k}: d={d:.3f} (was {float(np.linalg.norm(P0[k+6]-P0[k])):.3f})")


def snap_aromatic_rings_regular_polygon(syms: Sequence[str], P: np.ndarray,
                                        bond_target: float = 1.40,
                                        rms_threshold: float = 0.05
                                        ) -> Tuple[List[str], np.ndarray]:
    """G12b (Phase G12b, follow-up): plane snap PLUS in-plane regular-polygon snap.

    Phase G12 fixes ring planarity via SVD plane projection but does NOT fix
    in-plane C-C bond stretch (KIHVIH 1.96 → 1.95 = only −0.6 % improvement on
    the stretched bond). The Mogul TOP-blindspot ``C sp2 +31 %`` is exactly that
    stretch.  This routine adds a second step:

      1. Plane snap (out-of-plane → 0), same as snap_aromatic_rings.
      2. In-plane regular polygon snap: rescale each ring atom to the ideal
         radius r = bond_target / (2·sin(π/N)) while preserving the original
         angular position (no atom-swap / no chirality flip).

    M-D invariant guard: ring atoms coordinated to a metal are FROZEN (no move,
    plane or in-plane). Ring-attached H rigidly tracks parent net displacement
    (preserves X-H length).

    Pure geometry, FF-free, universal across all aromatic ring sizes (5-7).
    Deterministic. Default OFF (returns identity when AROMSNAP not active).
    Currently NOT wired into _maybe_relax; enable in next iter when G12 smoke
    confirms planarity-only is insufficient (will be exposed via a dedicated
    env flag DELFIN_FFFREE_AROMSNAP_REGULAR=1 once smoke evidence supports it).
    """
    if not _aromsnap_active():
        return list(syms), np.asarray(P, dtype=float).copy()
    P_out = np.asarray(P, dtype=float).copy()
    n = len(syms)

    bonds = _bd._geometric_bonds(list(syms), P_out)
    adj: List[List[int]] = [[] for _ in range(n)]
    heavy_adj: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        adj[i].append(j); adj[j].append(i)
        if syms[i] != "H" and syms[j] != "H":
            heavy_adj[i].append(j); heavy_adj[j].append(i)

    cycles = _enumerate_small_cycles(heavy_adj, max_len=7, min_len=5)
    if not cycles:
        return list(syms), P_out

    frozen = _metal_coordinated_atoms(list(syms), P_out)

    for ring in cycles:
        if not _is_aromatic_ring(ring, syms, P_out):
            continue
        ring_indices = list(ring)
        N = len(ring_indices)
        r_ideal = bond_target / (2.0 * float(np.sin(np.pi / N)))

        pts = P_out[ring_indices]
        centroid = pts.mean(axis=0)
        centered = pts - centroid
        try:
            _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            continue

        # Plane basis: u, v in-plane; w = ring normal.
        w = Vt[-1]
        nn = float(np.linalg.norm(w))
        if nn < 1e-9:
            continue
        w = w / nn
        u = Vt[0]
        u = u / max(1e-9, float(np.linalg.norm(u)))
        v = np.cross(w, u)

        # Project each atom into the (u, v, w) frame
        coords_uvw = np.array(
            [[float(np.dot(p, u)), float(np.dot(p, v)), float(np.dot(p, w))]
             for p in centered]
        )

        # Compute current radial position + angle in plane
        radii = np.sqrt(coords_uvw[:, 0] ** 2 + coords_uvw[:, 1] ** 2)
        angles = np.arctan2(coords_uvw[:, 1], coords_uvw[:, 0])

        # Targets: rescale to r_ideal, keep angular position, project w=0
        new_coords_uvw = np.array(
            [[r_ideal * float(np.cos(a)), r_ideal * float(np.sin(a)), 0.0]
             for a in angles]
        )

        # Back to xyz frame
        h_pairs = _ring_attached_hs(ring_indices, syms, P_out, adj)
        moves = {}
        for k, idx in enumerate(ring_indices):
            if idx in frozen:
                moves[idx] = np.zeros(3)
                continue
            new_xyz = (centroid +
                       new_coords_uvw[k, 0] * u +
                       new_coords_uvw[k, 1] * v +
                       new_coords_uvw[k, 2] * w)
            moves[idx] = new_xyz - P_out[idx]
        for idx, mv in moves.items():
            P_out[idx] = P_out[idx] + mv
        for h, parent in h_pairs:
            if parent in moves:
                P_out[h] = P_out[h] + moves[parent]

    return list(syms), P_out
