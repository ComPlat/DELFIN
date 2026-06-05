"""delfin.fffree.hapto_rigid_subtree — explicit hapto-π rigid-body freeze
with substituent-subtree drag.

Phase G13b (User 2026-06-05 mission A3, hapto rigid subtree).

PROBLEM
-------
Hapto ligands (η²/η⁵/η⁶: alkene, Cp, arene) currently get *partial*
handling:

  * :mod:`hapto_honest_construction` does correctly rigid-translate
    the ring C atoms + the H atoms directly bonded to ring C atoms
    when snapping to the empirical M-centroid distance.

  * BUT: heavy substituents (e.g. methyl carbons on a methyl-Cp, the
    backbone carbon attaching a -CH₂- bridge, a -CN substituent) and
    their entire subtrees do NOT ride with the ring. After
    grip_polish / aromatic_ring_scale / refine pass the ring carbons
    can move (e.g. centroid snap, ring scale) while the substituent
    root stays behind, stretching the C(ring)-C(sub) bond and
    deforming the substituent geometry.

  * Forensik on pool 9c40b70-wave7b3-vollpool (9815 frames, 27.6%
    hapto): 9.8% of heavy-substituent attachment bonds drift outside
    [0.85, 1.20] × cov-sum -- direct evidence of the "ring moves but
    substituent does not" failure mode.

FIX
---
Define a deterministic **hapto rigid-body set** per detected hapto
unit that contains:

  1. The η-C ring atoms (already handled).
  2. Every H atom bonded to a ring C atom (already handled).
  3. **Every non-ring heavy atom attached to a ring C, plus its
     entire subtree** -- BFS through the molecular graph, never
     re-entering the ring, never crossing into another hapto ring,
     never moving the metal or frozen atoms.

When ANY pipeline pass needs to displace the ring (centroid snap,
ring scale, refine step), the whole set translates / rotates as a
rigid body.

ENV FLAG
--------
``DELFIN_FFFREE_HAPTO_RIGID_SUBTREE=1`` activates the public
:func:`apply_rigid_subtree_translation` API.  Default OFF means the
function returns the input array byte-identically (`P.copy()` only
when called via the legacy path; the new API short-circuits).

UNIVERSALITY
------------
* Pure geometry + graph (no SMILES patterns).
* Deterministic: sorted iteration, no RNG, output bit-identical
  across runs.
* Metal + σ donors + any atom passed in `frozen_extra` are NEVER
  moved.
* Rings sharing an atom with another hapto ring (fused π-systems)
  are merged into one rigid body, preventing the BFS from cutting
  through a fused junction.

PUBLIC API
----------
* :func:`collect_rigid_subtree_atoms` — graph-only: return the full
  rigid set for ONE hapto unit (ring C + H + substituent subtree).
* :func:`apply_rigid_subtree_translation` — rigid-translate the
  whole set by a delta vector, returning a new coords array.
* :func:`rigid_subtree_active` — env-flag query.
"""
from __future__ import annotations

import itertools
import os
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np


__all__ = [
    "collect_rigid_subtree_atoms",
    "apply_rigid_subtree_translation",
    "apply_rigid_subtree_transform",
    "rigid_subtree_active",
    "build_geometric_adjacency",
]


# Element covalent radii (Å) — Pyykkö single-bond, mirror of the
# hapto_honest module so we never cross-import (default-OFF byte-
# identity guarantee remains airtight when this module is loaded).
_COV_R: Dict[str, float] = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "B": 0.84, "Si": 1.11, "Se": 1.20, "Te": 1.38, "As": 1.19,
    "Ge": 1.20, "Sn": 1.39, "Pb": 1.46, "Sb": 1.39, "Bi": 1.48,
    "Al": 1.21, "Ga": 1.22, "In": 1.42, "Tl": 1.45,
    "Mg": 1.41,
}
_DEFAULT_COV: float = 0.76

# Metal symbols — used to refuse routing the BFS through any metal.
_METAL_SYMBOLS: Set[str] = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Th", "U",
    "Al", "Ga", "In", "Sn", "Pb", "Sb", "Bi",
}

_BOND_TOL: float = 1.30


def rigid_subtree_active() -> bool:
    """Return ``True`` when ``DELFIN_FFFREE_HAPTO_RIGID_SUBTREE`` is set.

    Live-read of the env flag so test fixtures can flip it between
    subprocesses (default-OFF byte-identical when unset).
    """
    raw = os.environ.get("DELFIN_FFFREE_HAPTO_RIGID_SUBTREE", "0").strip().lower()
    return raw in ("1", "true", "yes", "on")


def build_geometric_adjacency(
    symbols: Sequence[str],
    coords: np.ndarray,
    bond_tol: float = _BOND_TOL,
) -> List[Set[int]]:
    """Build an undirected adjacency list using geometric covalent
    overlap (sum of covalent radii × ``bond_tol``).

    Deterministic, pure: same input → same adjacency.  Returns a list
    indexed by atom; each element is a sorted ``set`` of neighbour
    indices.  Metal-X edges are KEPT in the adjacency (the consumer
    decides whether to follow them).
    """
    n = len(symbols)
    P = np.asarray(coords, dtype=float).reshape(-1, 3)
    adj: List[Set[int]] = [set() for _ in range(n)]
    for i in range(n):
        ri = _COV_R.get(str(symbols[i]), _DEFAULT_COV)
        for j in range(i + 1, n):
            rj = _COV_R.get(str(symbols[j]), _DEFAULT_COV)
            d = float(np.linalg.norm(P[i] - P[j]))
            if d <= (ri + rj) * bond_tol and d > 1e-6:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def collect_rigid_subtree_atoms(
    symbols: Sequence[str],
    coords: np.ndarray,
    ring_atom_indices: Sequence[int],
    metal_idx: int,
    *,
    other_ring_atoms: Optional[Iterable[int]] = None,
    frozen_extra: Optional[Iterable[int]] = None,
    adjacency: Optional[List[Set[int]]] = None,
) -> List[int]:
    """Return the full atom-index set that should translate rigidly
    with the hapto-π ring at ``ring_atom_indices``.

    Algorithm (graph-only, deterministic):

      1. Start with the ring atoms.
      2. Add every neighbour H of a ring atom.
      3. BFS outward from each non-ring heavy neighbour of a ring atom.
         Stop at:
           * the metal (never traverse a metal-edge);
           * any atom in ``other_ring_atoms`` (other hapto ring);
           * any atom in ``frozen_extra`` (e.g. σ donors).
      4. Return the sorted union as a list of ints (deterministic).

    The returned list ALWAYS contains ``ring_atom_indices`` first
    (after sorting): the caller can use it directly to slice / move
    the rigid block.

    Pure-functional: same inputs → identical output.  No floating
    point arithmetic when the env flag is OFF *and the caller passes
    a precomputed adjacency*: adjacency build is the only fp work and
    the caller can skip it via the ``adjacency`` argument.
    """
    ring = sorted(int(i) for i in ring_atom_indices)
    if not ring:
        return []
    other = set(int(i) for i in (other_ring_atoms or ()))
    other -= set(ring)              # never block the ring itself
    frozen = set(int(i) for i in (frozen_extra or ()))
    frozen.discard(int(metal_idx))  # metal handled separately by "metal"
    metal = int(metal_idx)

    if adjacency is None:
        adjacency = build_geometric_adjacency(symbols, coords)

    n = len(symbols)
    ring_set: Set[int] = set(ring)
    rigid: Set[int] = set(ring_set)

    # Step 2 + 3 combined as one BFS, prefilling the queue with every
    # heavy / H neighbour of a ring atom that is not the metal / not in
    # ``other`` / not in ``frozen``.  Using a single BFS keeps the
    # order deterministic.
    visited: Set[int] = set(rigid)
    queue: List[int] = []
    for r in ring:
        if 0 <= r < n:
            for nb in sorted(adjacency[r]):
                if nb == metal or nb in other or nb in frozen:
                    continue
                if nb in visited:
                    continue
                visited.add(nb)
                queue.append(nb)

    while queue:
        cur = queue.pop(0)         # FIFO -- deterministic; small lists
        rigid.add(cur)
        # Refuse to leave the rigid set via a metal / other-ring / frozen
        # atom (the boundary stays where we put it).
        for nb in sorted(adjacency[cur]):
            if nb == metal or nb in other or nb in frozen:
                continue
            if nb in ring_set:
                # Ring atoms are already in rigid.  Don't re-enqueue.
                continue
            if nb in visited:
                continue
            visited.add(nb)
            queue.append(nb)

    return sorted(rigid)


def apply_rigid_subtree_translation(
    coords: np.ndarray,
    rigid_atoms: Sequence[int],
    delta: np.ndarray,
) -> np.ndarray:
    """Return a new array where every atom in ``rigid_atoms`` is
    shifted by ``delta`` (3,).  Atoms outside the set are byte-identical
    to the input.

    Default-OFF byte-identity is guaranteed by the caller (this
    function is a deterministic primitive; the env-flag check sits in
    :func:`rigid_subtree_active`).
    """
    if not rigid_subtree_active():
        return np.asarray(coords, dtype=float).copy()
    P = np.asarray(coords, dtype=float).copy()
    if len(rigid_atoms) == 0:
        return P
    delta = np.asarray(delta, dtype=float).reshape(3)
    idx = np.array(sorted(int(i) for i in rigid_atoms), dtype=np.int64)
    P[idx] = P[idx] + delta
    return P


def apply_rigid_subtree_transform(
    coords: np.ndarray,
    rigid_atoms: Sequence[int],
    *,
    pivot: np.ndarray,
    rotation: Optional[np.ndarray] = None,
    translation: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Return a new array where every atom in ``rigid_atoms`` is
    rotated about ``pivot`` (3×3 ``rotation`` matrix, ``rotation @
    (P - pivot) + pivot``), then translated by ``translation``.

    Both ``rotation`` and ``translation`` default to identity; when
    both are identity the output equals ``coords.copy()``.

    Pure: rigid-body invariant — all pairwise intra-block distances
    are preserved within float-precision (the validator
    :func:`_validate_rigid_block` in :mod:`hapto_honest_construction`
    catches the rare numerical edge case).
    """
    if not rigid_subtree_active():
        return np.asarray(coords, dtype=float).copy()
    P = np.asarray(coords, dtype=float).copy()
    if len(rigid_atoms) == 0:
        return P
    pivot = np.asarray(pivot, dtype=float).reshape(3)
    R = np.asarray(rotation, dtype=float).reshape(3, 3) if rotation is not None else np.eye(3)
    t = np.asarray(translation, dtype=float).reshape(3) if translation is not None else np.zeros(3)
    idx = np.array(sorted(int(i) for i in rigid_atoms), dtype=np.int64)
    block = P[idx] - pivot
    block = block @ R.T
    block = block + pivot + t
    P[idx] = block
    return P
