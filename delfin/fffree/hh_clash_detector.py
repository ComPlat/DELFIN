"""delfin.fffree.hh_clash_detector — Specialised H-H clash detection.

Background
----------
The 3 generic clash modules in :mod:`delfin.fffree` (``ClashFloorPenalty``,
``inter_ligand_clash_gate``, ``build_time_clash_gate``) all DO include
hydrogen in their vdW tables (``r_H = 1.20 Å``).  An audit of the
production voll-pool ``b00f9a0-full7-VOLLPOOL`` (2026-06-04) nevertheless
found that ~39 % of emitted structures contain non-trivial H-H close
contacts below the Bondi-vdW floor (``0.85 × 2 r_H = 2.04 Å``):

  * methyl-methyl eclipsing across an inter-ligand axis
    (H-H ≈ 1.95-2.03 Å, 3-6 pair-clashes per residue)
  * methylene-methyl rotamer mistakes (H-H ≈ 1.77-1.97 Å)
  * aromatic-H pointing into a bulky substituent (H-H ≈ 2.0 Å)

The root cause is twofold:

  1. ``assemble_complex._clash_count`` uses the loose factor ``0.70``
     (so H-H < 0.70 × 2.40 = 1.68 Å is the firing threshold; eclipsing
     pairs at 1.95-2.03 Å pass clean).
  2. The build-time gate's ``XH_COLLAPSE_TOL`` is the bond-collapse floor
     (0.70 Å absolute, for X-H bonded pairs only) -- it has nothing to
     say about non-bonded H-H contacts.

This module fills that gap with a focused detector tuned to the Bondi
geometry:

  - ``detect_hh_clashes(...)`` -- list every non-bonded, non-geminal,
    non-1,3 H-H pair below ``factor × (r_i + r_j)``.
  - ``count_hh_clashes(...)`` -- the same, returning only the count
    (used in selection / ranking loops).
  - ``hh_clash_severity(...)`` -- the sum of squared overlaps used as a
    penalty term in the L-BFGS loss when the env-flag is active.
  - ``hh_clash_active()`` -- env-flag predicate, default OFF, byte-
    identical to HEAD ``93b396d``.

Class-conditional handling
--------------------------
Not every short H-H distance is a clash.  The detector skips:

  * geminal H-H (both H attached to the same heavy atom; e.g. methylene
    geminal pair ~ 1.78 Å)
  * 1,3 H-H within the same methyl group (3 internal H on the same C
    sit at ~ 1.78 Å -- chemically correct, NOT a clash)
  * 1,4 H-H within the same ring carbon environment when the parent
    heavy atoms are directly bonded (eclipsing across an sp3-sp3 bond
    is treated separately -- see ``DELFIN_FFFREE_HH_CLASH_FOLD_14``)

Filter logic uses topological distance computed from a minimal
covalent-bond graph (derived from the heavy-atom coordinates via
:mod:`delfin._bond_decollapse` ideal-bond table).  This keeps the
detector graph-only, universal, deterministic, and free of SMILES
patterns.

Environment flags
-----------------

``DELFIN_FFFREE_HH_CLASH_INCLUDE``
    Activates the detector.  ``"1"`` (default OFF) enables the
    detector inside ``ClashFloorPenalty`` and the post-build gates.

``DELFIN_FFFREE_HH_CLASH_FACTOR``
    vdW-sum fraction below which H-H pairs are flagged (default
    ``0.85`` -- Bondi).

``DELFIN_FFFREE_HH_CLASH_MIN_TOPO``
    Minimum graph-distance between the two H atoms required to flag a
    pair (default ``4`` -- skip geminal, 1-3 within CH3, and direct
    1-4 across an sp3 bond which is normally a rotamer issue rather
    than a clash).  Set to ``3`` to also flag 1-3-H atom pairs whose
    parents are NOT bonded (e.g. inter-ligand H-H 3 bonds apart).

``DELFIN_FFFREE_HH_CLASH_WEIGHT``
    L-BFGS penalty weight (default ``5.0`` -- matches the heavy-pair
    intra weight in :class:`grip_constraints.ClashFloorPenalty`).

Determinism
-----------
All loops are sorted by ``(i, j)`` with ``i < j``; the bond graph is
constructed from a fixed covalent-radii table; no RNG; no platform
calls.  Same input ⇒ same output byte-for-byte.

Performance
-----------
Naive H-H pair iteration is ``O(n_H^2)`` which is < ``10^4`` for a
typical 200-atom complex (~ 50 H atoms).  No spatial hashing needed.
Bond-graph construction is ``O(n_heavy^2)`` once and cached per call.
"""
from __future__ import annotations

import os
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np


__all__ = [
    "DEFAULT_HH_FACTOR",
    "DEFAULT_HH_MIN_TOPO",
    "DEFAULT_HH_WEIGHT",
    "H_VDW_RADIUS",
    "hh_clash_active",
    "hh_clash_factor",
    "hh_clash_min_topo",
    "hh_clash_weight",
    "detect_hh_clashes",
    "count_hh_clashes",
    "hh_clash_severity",
    "build_hh_exclusion_pairs",
]


# Bondi (1964) vdW radius for hydrogen.  Matches the value used in
# ``grip_polish.DEFAULT_VDW_RADII`` / ``inter_ligand_clash_gate`` /
# ``build_time_clash_gate.INTERLIG_VDW_RADII``.
H_VDW_RADIUS = 1.20

DEFAULT_HH_FACTOR = 0.85         # Bondi-vdW floor fraction
DEFAULT_HH_MIN_TOPO = 4          # skip geminal (2), 1-3 (3); flag 1-4+
DEFAULT_HH_WEIGHT = 5.0          # match intra-clash weight default


# ---------------------------------------------------------------------------
# Env-flag plumbing
# ---------------------------------------------------------------------------
def hh_clash_active() -> bool:
    """Return ``True`` when the H-H clash detector is enabled.

    Default OFF.  When the env-flag ``DELFIN_FFFREE_HH_CLASH_INCLUDE``
    is unset or set to ``0``, every call site that imports this module
    behaves byte-identically to HEAD ``93b396d``.
    """
    raw = os.environ.get("DELFIN_FFFREE_HH_CLASH_INCLUDE", "0")
    return str(raw).strip().lower() in ("1", "true", "on", "yes")


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default


def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except (TypeError, ValueError):
        return default


def hh_clash_factor() -> float:
    """Read ``DELFIN_FFFREE_HH_CLASH_FACTOR`` (default 0.85)."""
    return _env_float("DELFIN_FFFREE_HH_CLASH_FACTOR", DEFAULT_HH_FACTOR)


def hh_clash_min_topo() -> int:
    """Read ``DELFIN_FFFREE_HH_CLASH_MIN_TOPO`` (default 4)."""
    return _env_int("DELFIN_FFFREE_HH_CLASH_MIN_TOPO", DEFAULT_HH_MIN_TOPO)


def hh_clash_weight() -> float:
    """Read ``DELFIN_FFFREE_HH_CLASH_WEIGHT`` (default 5.0)."""
    return _env_float("DELFIN_FFFREE_HH_CLASH_WEIGHT", DEFAULT_HH_WEIGHT)


# ---------------------------------------------------------------------------
# Bond-graph helpers (geometry-only, deterministic, no RDKit dependency)
# ---------------------------------------------------------------------------
def _build_bond_graph(
    syms: Sequence[str],
    P: np.ndarray,
    *,
    bond_factor: float = 1.30,
) -> List[Set[int]]:
    """Construct a covalent-bond adjacency list from coordinates.

    A pair ``(i, j)`` is considered bonded when

        ``d(i, j) < bond_factor × ideal_bond(syms[i], syms[j])``

    using the ideal-bond table from :mod:`delfin._bond_decollapse`.
    Metal-incident bonds are kept (an M-H hydride is a genuine bond);
    the topological-distance filter handles them correctly.

    Returns
    -------
    adj : list of set[int]
        ``adj[i]`` is the set of atom indices bonded to ``i``.

    Notes
    -----
    Used internally to compute the topological distance between H
    pairs.  Deterministic: iteration order is purely by ``(i, j)``.
    """
    # Lazy import keeps the module standalone for unit tests with a
    # minimal stub bond table.
    try:
        from delfin import _bond_decollapse as _bd
        ideal_bond = _bd._ideal_bond
    except Exception:
        # Conservative fallback: single covalent radii table.
        _COV = {
            "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
            "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02,
            "Br": 1.20, "I": 1.39, "B": 0.84,
        }
        def ideal_bond(a: str, b: str) -> float:
            return _COV.get(a, 1.0) + _COV.get(b, 1.0)

    n = len(syms)
    adj: List[Set[int]] = [set() for _ in range(n)]
    P = np.asarray(P, dtype=float)
    for i in range(n):
        si = str(syms[i])
        Pi = P[i]
        for j in range(i + 1, n):
            sj = str(syms[j])
            d = float(np.linalg.norm(Pi - P[j]))
            if d < bond_factor * ideal_bond(si, sj):
                adj[i].add(j)
                adj[j].add(i)
    return adj


def _topo_distance_from(
    adj: Sequence[Set[int]],
    src: int,
    max_depth: int,
) -> Dict[int, int]:
    """BFS up to ``max_depth`` returning ``{node: distance_from_src}``.

    The result includes ``src`` itself at distance 0.  Nodes farther
    than ``max_depth`` away are omitted.
    """
    out: Dict[int, int] = {int(src): 0}
    frontier: List[int] = [int(src)]
    for depth in range(1, int(max_depth) + 1):
        nxt: List[int] = []
        for v in frontier:
            for w in adj[v]:
                wi = int(w)
                if wi not in out:
                    out[wi] = depth
                    nxt.append(wi)
        if not nxt:
            break
        frontier = nxt
    return out


# ---------------------------------------------------------------------------
# Core detector
# ---------------------------------------------------------------------------
def detect_hh_clashes(
    syms: Sequence[str],
    P: np.ndarray,
    *,
    factor: Optional[float] = None,
    min_topo_distance: Optional[int] = None,
    excluded_pairs: Optional[Iterable[Iterable[int]]] = None,
    bond_adj: Optional[Sequence[Set[int]]] = None,
) -> List[Tuple[int, int, float, float]]:
    """Return every H-H pair below the vdW floor as
    ``(i, j, distance, severity)``.

    Parameters
    ----------
    syms : sequence of str
        Element symbol per atom.
    P : (N, 3) ndarray
        Cartesian coordinates.
    factor : float, optional
        vdW-sum fraction (default :data:`DEFAULT_HH_FACTOR` = 0.85).
        Pairs with ``d < factor × (r_i + r_j)`` are flagged.
    min_topo_distance : int, optional
        Minimum number of bonds between ``i`` and ``j`` for the pair to
        count (default :data:`DEFAULT_HH_MIN_TOPO` = 4).  Geminal
        (topo=2) and 1-3 within CH3 (topo=3) are skipped.
    excluded_pairs : iterable of (int, int), optional
        Additional pairs to skip (e.g. caller-supplied "approved"
        close contacts).
    bond_adj : list of set[int], optional
        Pre-computed bond adjacency.  Recomputed from coordinates when
        ``None`` (allows callers to cache the graph across calls).

    Returns
    -------
    list of (i, j, distance, severity)
        ``severity = max(0, d_min - d)^2`` (Å²); higher = worse clash.
        Pairs are sorted by ``(i, j)`` ascending; ``i < j``.

    Determinism
    -----------
    Same input ⇒ same output byte-for-byte.  Iteration order is sorted.
    """
    P = np.asarray(P, dtype=float)
    if P.ndim == 1:
        if P.size % 3 != 0:
            raise ValueError("P must reshape to (N, 3)")
        P = P.reshape(-1, 3)
    n = P.shape[0]
    if len(syms) != n:
        raise ValueError(f"syms length {len(syms)} != coords rows {n}")

    fac = float(factor) if factor is not None else hh_clash_factor()
    min_topo = int(min_topo_distance) if min_topo_distance is not None else hh_clash_min_topo()

    H_idx = [i for i in range(n) if str(syms[i]) == "H"]
    if len(H_idx) < 2:
        return []

    if bond_adj is None:
        adj = _build_bond_graph(syms, P)
    else:
        adj = list(bond_adj)

    # Pre-compute topological neighbourhood up to ``min_topo - 1`` for
    # every H -- pairs that fall inside it are skipped.  We BFS at most
    # ``min_topo - 1`` deep to avoid full-graph BFS.
    cutoff = max(0, min_topo - 1)
    topo_neigh: Dict[int, Dict[int, int]] = {}
    if cutoff > 0:
        for h in H_idx:
            topo_neigh[h] = _topo_distance_from(adj, h, cutoff)

    excl_set: Set[FrozenSet[int]] = set()
    if excluded_pairs:
        for pair in excluded_pairs:
            try:
                a, b = tuple(pair)
            except Exception:
                continue
            if int(a) == int(b):
                continue
            excl_set.add(frozenset((int(a), int(b))))

    d_min = fac * (H_VDW_RADIUS + H_VDW_RADIUS)
    out: List[Tuple[int, int, float, float]] = []
    for ii in range(len(H_idx)):
        i = H_idx[ii]
        Pi = P[i]
        nh_i = topo_neigh.get(i, {})
        for jj in range(ii + 1, len(H_idx)):
            j = H_idx[jj]
            if frozenset((i, j)) in excl_set:
                continue
            # Topological-distance filter (skip geminal + 1-3 by default).
            if cutoff > 0 and j in nh_i:
                # j is within cutoff bonds of i -> topo distance < min_topo
                continue
            d = float(np.linalg.norm(Pi - P[j]))
            if d >= d_min:
                continue
            # Skip exact-coincidence noise (already handled by build-time
            # gate; reporting it here would be redundant).
            if d < 1e-9:
                continue
            sev = (d_min - d) * (d_min - d)
            out.append((int(i), int(j), float(d), float(sev)))
    out.sort(key=lambda t: (t[0], t[1]))
    return out


def count_hh_clashes(
    syms: Sequence[str],
    P: np.ndarray,
    *,
    factor: Optional[float] = None,
    min_topo_distance: Optional[int] = None,
    excluded_pairs: Optional[Iterable[Iterable[int]]] = None,
    bond_adj: Optional[Sequence[Set[int]]] = None,
) -> int:
    """Return the number of H-H pairs flagged by :func:`detect_hh_clashes`.

    Convenience wrapper for selection-loop ranking; matches the signature
    of :func:`build_time_clash_gate.collapse_count`.
    """
    return len(detect_hh_clashes(
        syms, P,
        factor=factor,
        min_topo_distance=min_topo_distance,
        excluded_pairs=excluded_pairs,
        bond_adj=bond_adj,
    ))


def hh_clash_severity(
    syms: Sequence[str],
    P: np.ndarray,
    *,
    factor: Optional[float] = None,
    min_topo_distance: Optional[int] = None,
    excluded_pairs: Optional[Iterable[Iterable[int]]] = None,
    bond_adj: Optional[Sequence[Set[int]]] = None,
) -> float:
    """Total severity ``Σ max(0, d_min - d)^2`` (Å²) over flagged pairs.

    Used as a candidate-selection penalty (lower is better) and as the
    contribution to the L-BFGS loss when ``ClashFloorPenalty`` includes
    H atoms.  Returns 0.0 when no clashes are found.
    """
    return float(sum(t[3] for t in detect_hh_clashes(
        syms, P,
        factor=factor,
        min_topo_distance=min_topo_distance,
        excluded_pairs=excluded_pairs,
        bond_adj=bond_adj,
    )))


def build_hh_exclusion_pairs(
    syms: Sequence[str],
    P: np.ndarray,
    *,
    min_topo_distance: Optional[int] = None,
    bond_adj: Optional[Sequence[Set[int]]] = None,
) -> Set[FrozenSet[int]]:
    """Return the frozenset of H-H pairs that should be excluded from a
    generic ``ClashFloorPenalty`` (geminal + 1-3 H-H within CH3).

    This lets the caller plug the existing :class:`ClashFloorPenalty`
    (which only excludes 1-2 + 1-3 in the heavy graph) into an H-aware
    mode without double-counting the geminal H-H pairs.

    Each frozenset contains exactly two atom indices.
    """
    P = np.asarray(P, dtype=float)
    n = P.shape[0] if P.ndim == 2 else (P.size // 3)
    H_idx = [i for i in range(n) if str(syms[i]) == "H"]
    if len(H_idx) < 2:
        return set()
    min_topo = int(min_topo_distance) if min_topo_distance is not None else hh_clash_min_topo()
    cutoff = max(0, min_topo - 1)
    if cutoff <= 0:
        return set()
    if bond_adj is None:
        adj = _build_bond_graph(syms, P)
    else:
        adj = list(bond_adj)
    excl: Set[FrozenSet[int]] = set()
    for h in H_idx:
        nh = _topo_distance_from(adj, h, cutoff)
        for other, _d in nh.items():
            if other == h:
                continue
            if str(syms[other]) != "H":
                continue
            if other < h:
                # Will be added when we visit ``other``; keep loop tight.
                continue
            excl.add(frozenset((int(h), int(other))))
    return excl
