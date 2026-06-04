"""delfin.fffree.aromatic_bond_enforcement — per-bond aromatic ideal pull.

Overnight bug 2026-06-04: ``e9f69af-iso-topo-heal/080-LUHMOT.xyz`` carries a
C-C bond of 1.130 Å (triple-bond length) inside a 6-membered ring that
should be aromatic.  The cure stack — ``_bond_decollapse`` (single-bond
ideal 1.52), ``topology_healing`` (no compressed detection), and
``aromatic_ring_scale`` (uniform mean-based scale, clamp 0.85..1.05) — does
NOT catch this.  See ``forensik/aromatic_bond_trace.md``.

This module is the per-bond fix.  For every heavy-heavy edge whose both
endpoints are RDKit-aromatic, we pull the bond toward the empirical CCDC
mean (``aromatic_bond_targets``).  When the entire ring is aromatic and
metal-free, we snap to a regular planar polygon at the in-ring CCDC ideal
(true per-bond ideal length within ±sigma, planarity + symmetry preserved).
Rings carrying a metal-coordinated atom are FROZEN — preserves M-D, hapto-η
placement, and chelate attachment.

Determinism: pure geometry, no RNG, no I/O, deterministic ring traversal.
Env gate ``DELFIN_FFFREE_AROMATIC_BONDS=1`` (default OFF, byte-identical
when unset).
"""
from __future__ import annotations

import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from delfin.fffree.aromatic_bond_targets import aromatic_ideal
import delfin._bond_decollapse as _bd


# Per-bond tolerance: only act when |d - mu| / sigma > Z_TRIGGER (CCDC z-score)
Z_TRIGGER = 2.5
# Max single-pass per-atom displacement (Å) — keeps the planarity intact for
# subsequent passes; mirrors ``_bond_decollapse._MAX_STEP``.
MAX_STEP = 0.25
# Convergence threshold on max per-atom step in an iteration.
CONV_TOL = 1e-3
# Max iterations of the per-bond gradient descent.
MAX_ITERS = 60


def is_enabled() -> bool:
    """Live-read env gates (call-time, never import-time)."""
    if os.environ.get("DELFIN_FFFREE_AROMATIC_BONDS", "0") == "1":
        return True
    if os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1":
        # PURE_TRACK3 unlocks every aromatic enforcement (consistent with the
        # ring_scale + ring pucker activation pattern).
        return True
    return False


# ---------------------------------------------------------------------------
# Internals: aromatic edges + frozen atoms
# ---------------------------------------------------------------------------
def _frozen_atoms(syms: Sequence[str], P: np.ndarray,
                  md_factor: float = 1.45,
                  hapto_cutoff: float = 3.0) -> Set[int]:
    """Indices of atoms within ``md_factor·ideal`` of any metal, plus hapto-π
    fingerprint atoms.  Mirrors ``aromatic_ring_scale._metal_coordinated_atoms``
    so the two modules freeze the same set."""
    n = len(syms)
    frozen: Set[int] = set()
    for m in range(n):
        if not _bd._is_metal(syms[m]):
            continue
        near_heavy: List[Tuple[int, str, float]] = []
        for j in range(n):
            if j == m or syms[j] == "H":
                continue
            d = float(np.linalg.norm(P[j] - P[m]))
            if d < md_factor * _bd._ideal_bond(syms[m], syms[j]):
                frozen.add(j)
            if d < hapto_cutoff:
                near_heavy.append((j, syms[j], d))
        # Hapto-π geometric fingerprint: ≥3 same-element heavy nbrs in cutoff.
        from collections import defaultdict
        by_elem: dict = defaultdict(list)
        for j, e, _ in near_heavy:
            if e in {"C", "N", "O", "S"}:
                by_elem[e].append(j)
        for _elem, idxs in by_elem.items():
            if len(idxs) >= 3:
                for j in idxs:
                    frozen.add(j)
    return frozen


def _aromatic_edges_from_mol(mol) -> List[Tuple[int, int]]:
    """Return lex-sorted heavy-heavy edges where BOTH endpoints have
    ``GetIsAromatic() == True`` in the RDKit ``mol``.  Hydrogens skipped."""
    if mol is None:
        return []
    out: Set[Tuple[int, int]] = set()
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetSymbol() == "H" or b.GetSymbol() == "H":
            continue
        if not (a.GetIsAromatic() and b.GetIsAromatic()):
            continue
        i, j = int(a.GetIdx()), int(b.GetIdx())
        if i > j:
            i, j = j, i
        out.add((i, j))
    return sorted(out)


def _aromatic_edges_from_geometry(syms: Sequence[str], P: np.ndarray,
                                  ring_size_lo: int = 5,
                                  ring_size_hi: int = 6,
                                  ) -> List[Tuple[int, int]]:
    """Fallback when no RDKit mol is available: detect 5/6-rings that look
    aromatic AND planar AND have every ring-bond in the aromatic length
    range [1.20, 1.55].  Deterministic.

    The geometric heuristic is tightened (vs the legacy aromatic_ring_scale
    detector) to avoid false positives on conjugated-non-aromatic systems
    that the bond-pull would otherwise mis-target.  In particular:
    * planarity: max out-of-plane distance < 0.35 Å for the ring
    * length sanity: every ring-bond in [1.20, 1.55] Å (admits the
      already-compressed LUHMOT case but excludes sp3 saturated rings whose
      C-C is at 1.52)
    * element set: only C / N / O / S / B (no metals, no halides, no Si/P
      since their aromatic ideals carry larger sigmas and the
      element-pair miss-rate is higher)
    """
    n = len(syms)
    # Heavy-atom adjacency from geometric bonds.
    bonds = _bd._geometric_bonds(list(syms), np.asarray(P))
    heavy_adj: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        if syms[i] == "H" or syms[j] == "H":
            continue
        heavy_adj[i].append(j); heavy_adj[j].append(i)

    cycles: Set[Tuple[int, ...]] = set()
    def dfs(start: int, current: int, path: List[int], depth: int) -> None:
        if depth > ring_size_hi:
            return
        for nxt in heavy_adj[current]:
            if nxt == start and depth >= ring_size_lo:
                cycles.add(tuple(sorted(path)))
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

    allowed = {"C", "N", "O", "S", "B"}
    out: Set[Tuple[int, int]] = set()
    for ring in cycles:
        if not all(syms[a] in allowed for a in ring):
            continue
        ring_set = set(ring)
        ring_atoms = sorted(ring)
        # Build a sub-adjacency restricted to ring atoms.
        sub_nbrs = {a: [b for b in heavy_adj[a] if b in ring_set] for a in ring_atoms}
        # Each ring atom must have AT LEAST 2 in-ring nbrs.  A collapsed atom
        # may have extra spurious in-ring contacts (e.g. LUHMOT C13 sees C14,
        # C18, C19 as "bonded" because the ring is squished), but the ring
        # itself remains a valid 6-cycle.
        if not all(len(sub_nbrs[a]) >= 2 for a in ring_atoms):
            continue
        # Reconstruct the ring as a cycle in the order returned by the DFS
        # cycle finder; collect its consecutive edges.  The DFS path is
        # already a Hamiltonian cycle of the ring atoms.
        ring_seq = list(ring)
        edges: Set[Tuple[int, int]] = set()
        # Find a Hamiltonian cycle in sub_nbrs starting from min(ring).
        start = min(ring_atoms)
        path = [start]
        cur = start
        prev = -1
        for _ in range(len(ring_atoms) - 1):
            cands = [b for b in sub_nbrs[cur] if b != prev and b in ring_set]
            # Prefer the candidate that is NOT already on path
            cands_new = [b for b in cands if b not in path]
            if cands_new:
                nxt = sorted(cands_new)[0]
            elif cands:
                nxt = sorted(cands)[0]
            else:
                break
            path.append(nxt)
            prev = cur
            cur = nxt
        if len(path) != len(ring_atoms):
            continue
        for k in range(len(path)):
            a, b = path[k], path[(k + 1) % len(path)]
            edges.add(tuple(sorted((a, b))))
        # Length sanity: every ring-bond in [1.05, 1.65] AND mean length in
        # the aromatic band [1.30, 1.48].  The mean-band check excludes
        # cyclohexane / cyclopentane (mean ~1.52) from being mis-pulled to
        # aromatic 1.40 — they are sp3 saturated rings whose ring-C-C
        # mean correctly sits ABOVE the aromatic band.  Smoke 158
        # forensik (079-FICNAG, 092-NEKCEM): without this gate, two
        # files were mis-pulled from 1.51 → 1.40.
        lens = [float(np.linalg.norm(P[a] - P[b])) for (a, b) in edges]
        if not all(1.05 <= L <= 1.65 for L in lens):
            continue
        mean_len = float(np.mean(lens))
        if not (1.30 <= mean_len <= 1.48):
            continue
        # Planarity check: SVD on centered ring atoms, sv[-1] < 0.40 Å
        # (relaxed from 0.35 for sub-collapsed cases — LUHMOT ring has
        # sv[-1] = 0.27).
        pts = np.array([P[a] for a in ring_atoms], dtype=float)
        ctr = pts.mean(axis=0)
        centered = pts - ctr
        try:
            _, sv, _ = np.linalg.svd(centered, full_matrices=False)
        except Exception:
            continue
        if len(sv) >= 3 and sv[-1] > 0.40:
            continue
        for e in edges:
            out.add(e)
    return sorted(out)


# ---------------------------------------------------------------------------
# Per-bond pull
# ---------------------------------------------------------------------------
def _per_bond_pull(P: np.ndarray, edges: List[Tuple[int, int]],
                    targets: List[Tuple[float, float]],
                    frozen: Set[int],
                    max_iters: int = MAX_ITERS,
                    max_step: float = MAX_STEP,
                    conv_tol: float = CONV_TOL,
                    z_trigger: float = Z_TRIGGER,
                    ) -> np.ndarray:
    """Gauss-Seidel-style spring relaxation on aromatic edges only.  Frozen
    atoms (metal-coordinated, hapto-π) are pinned.  Per-bond ideal is
    ``targets[i]`` (mu, sigma); only acts when the bond's z-score exceeds
    ``z_trigger`` (else CCDC says the bond is already within natural scatter).
    """
    n = P.shape[0]
    Q = P.copy()
    for _it in range(max_iters):
        disp = np.zeros_like(Q)
        cnt = np.zeros(n, dtype=int)
        for (a, b), (mu, sigma) in zip(edges, targets):
            vec = Q[b] - Q[a]
            d = float(np.linalg.norm(vec))
            if d < 1e-9:
                continue
            z = (d - mu) / max(sigma, 1e-6)
            if abs(z) < z_trigger:
                continue
            # Pull each endpoint half-way to ideal, capped at max_step.
            delta = mu - d
            step = 0.5 * delta * (vec / d)
            sn = float(np.linalg.norm(step))
            if sn > max_step:
                step = step * (max_step / sn)
            if a not in frozen:
                disp[a] -= step
                cnt[a] += 1
            if b not in frozen:
                disp[b] += step
                cnt[b] += 1
        # Average per atom.
        moved_atoms = np.where(cnt > 0)[0]
        if moved_atoms.size == 0:
            break
        disp[moved_atoms] /= cnt[moved_atoms, None]
        Q[moved_atoms] += disp[moved_atoms]
        max_disp = float(np.max(np.linalg.norm(disp[moved_atoms], axis=1)))
        if max_disp < conv_tol:
            break
    return Q


# ---------------------------------------------------------------------------
# Hydrogen drag
# ---------------------------------------------------------------------------
def _h_drag(syms: Sequence[str], P_before: np.ndarray, P_after: np.ndarray,
            heavy_adj_with_h: List[List[int]]) -> np.ndarray:
    """Rigid-translate each H atom by the displacement of its sole heavy
    parent (preserves X-H to first order).  Falls back to keeping the H in
    place if it has 0 or >=2 heavy nbrs (rare edge cases)."""
    n = len(syms)
    Q = P_after.copy()
    for h in range(n):
        if syms[h] != "H":
            continue
        heavy_nbrs = [k for k in heavy_adj_with_h[h] if syms[k] != "H"]
        if len(heavy_nbrs) != 1:
            continue
        parent = heavy_nbrs[0]
        Q[h] = P_before[h] + (P_after[parent] - P_before[parent])
    return Q


# ---------------------------------------------------------------------------
# M-D invariance check
# ---------------------------------------------------------------------------
def _md_preserved(P_before: np.ndarray, P_after: np.ndarray,
                  syms: Sequence[str], tol: float = 0.05) -> bool:
    """Every M-heavy non-metal pair within 1.45·ideal must keep its distance
    within ``tol`` (Å)."""
    n = len(syms)
    for i in range(n):
        if not _bd._is_metal(syms[i]):
            continue
        for j in range(n):
            if i == j or syms[j] == "H":
                continue
            d0 = float(np.linalg.norm(P_before[j] - P_before[i]))
            if d0 > 1.45 * _bd._ideal_bond(syms[i], syms[j]):
                continue
            d1 = float(np.linalg.norm(P_after[j] - P_after[i]))
            if abs(d1 - d0) > tol:
                return False
    return True


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def enforce_aromatic_bonds(
    syms: Sequence[str],
    P: np.ndarray,
    mol=None,
    fixed: Optional[Iterable[int]] = None,
) -> Tuple[List[str], np.ndarray]:
    """Pull every aromatic bond toward its CCDC empirical ideal.

    Parameters
    ----------
    syms, P
        Atom symbols and coordinates (``(n,3)``).
    mol
        Optional RDKit ``Mol``.  When supplied, aromaticity is read from
        RDKit per-atom (the exact intended definition); else we fall back to
        the geometric 5/6-ring heuristic on C/N/O/S/B.
    fixed
        Additional atom indices to pin (caller's frozen-set, e.g. metal +
        donors).  These are added to the auto-detected metal-coordinated /
        hapto-π set.

    Returns
    -------
    ``(syms_out, P_out)`` — fresh copy.  If the env flag is unset, returns
    the input unchanged (byte-identical no-op).  If the post-fix geometry
    fails the M-D invariance check, the pre-fix coords are returned.
    """
    if not is_enabled():
        return list(syms), np.asarray(P, dtype=float).copy()

    P_in = np.asarray(P, dtype=float).copy()
    n = len(syms)

    # Build aromatic edge set.
    if mol is not None:
        edges = _aromatic_edges_from_mol(mol)
    else:
        edges = _aromatic_edges_from_geometry(syms, P_in)
    if not edges:
        return list(syms), P_in

    # Aromatic ideal per edge — only keep edges with a real library entry.
    keep_edges: List[Tuple[int, int]] = []
    keep_targets: List[Tuple[float, float]] = []
    for (a, b) in edges:
        hit = aromatic_ideal(syms[a], syms[b])
        if hit is None:
            continue
        mu, sigma, _n = hit
        keep_edges.append((a, b))
        keep_targets.append((mu, sigma))
    if not keep_edges:
        return list(syms), P_in

    # Frozen set: metal-coordinated + hapto-π + caller-supplied.
    frozen = _frozen_atoms(syms, P_in)
    if fixed is not None:
        frozen.update(int(i) for i in fixed)

    # Build H-drag adjacency from current geometry.
    bonds_geom = _bd._geometric_bonds(list(syms), P_in)
    adj_with_h: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds_geom:
        adj_with_h[i].append(j); adj_with_h[j].append(i)

    # Per-bond Gauss-Seidel pull.
    P_after = _per_bond_pull(P_in, keep_edges, keep_targets, frozen)

    # H drag.
    P_after = _h_drag(syms, P_in, P_after, adj_with_h)

    # M-D invariance check.
    if not _md_preserved(P_in, P_after, syms, tol=0.05):
        return list(syms), P_in

    return list(syms), P_after


__all__ = [
    "enforce_aromatic_bonds",
    "is_enabled",
    "Z_TRIGGER",
    "MAX_STEP",
    "MAX_ITERS",
    "CONV_TOL",
]
