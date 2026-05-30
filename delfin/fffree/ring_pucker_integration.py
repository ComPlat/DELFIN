"""delfin.fffree.ring_pucker_integration — integrate Cremer-Pople ring pucker
enumeration with the fffree pipeline (post-assembly variant generator).

Universal: works for any TMC or organic molecule with rings.
FF-free: pure geometric Cremer-Pople synthesis.
Deterministic: stable enumeration order.

Env-gate: DELFIN_FFFREE_RING_PUCKER=1 (auto-on under PURE_TRACK3).
"""
from __future__ import annotations
import os
from typing import Iterator, List, Optional, Sequence, Tuple

import numpy as np

from .ring_pucker import (
    canonical_pucker_states, set_pucker, _restore_ring_bonds,
    compute_pucker, _RING_PUCKER
)


def _is_aromatic_ring(mol, ring_atoms: Sequence[int]) -> bool:
    """Check if all atoms in the ring are aromatic (RDKit IsAromatic).
    Aromatic rings stay planar — no pucker enumeration."""
    try:
        return all(mol.GetAtomWithIdx(int(i)).GetIsAromatic() for i in ring_atoms)
    except Exception:
        return False


def _contains_metal(mol, ring_atoms: Sequence[int]) -> bool:
    """Detect chelate ring (contains a metal atom)."""
    try:
        from delfin._bond_decollapse import _is_metal
        return any(_is_metal(mol.GetAtomWithIdx(int(i)).GetSymbol()) for i in ring_atoms)
    except Exception:
        return False


def enumerate_mol_ring_conformers(
    mol,
    coords: np.ndarray,
    max_per_ring: int = 4,
    skip_aromatic: bool = True,
    chelate_only: bool = False,
) -> Iterator[np.ndarray]:
    """Enumerate ring-pucker conformer variants of a molecule.

    Parameters
    ----------
    mol : RDKit Mol (for ring topology + aromaticity detection)
    coords : (M, 3) atom positions (must align with mol atom order)
    max_per_ring : cap on canonical pucker states per ring
    skip_aromatic : if True, planar aromatic rings keep their flat geometry
    chelate_only : if True, only enumerate puckers for rings containing a metal
        (useful for chelate δ/λ enumeration)

    Yields
    ------
    coords_variant : (M, 3) array with target pucker applied to puckerable rings.
        Bond lengths restored via iterative axis-pull.

    Notes
    -----
    Returns NOTHING if env DELFIN_FFFREE_RING_PUCKER is unset (default OFF,
    byte-identical to no integration).

    The first yielded variant is the original (no pucker change) so callers
    can treat the original as the leader of the conformer family.
    """
    if not _RING_PUCKER:
        return
    try:
        rings = list(mol.GetRingInfo().AtomRings())
    except Exception:
        rings = []
    # Filter: only puckerable rings
    puckerable: List[Sequence[int]] = []
    for r in rings:
        if len(r) < 4 or len(r) > 12:
            continue
        if skip_aromatic and _is_aromatic_ring(mol, r):
            continue
        if chelate_only and not _contains_metal(mol, r):
            continue
        puckerable.append(r)
    if not puckerable:
        return
    # Yield original first
    yield coords.copy()
    # Cartesian product over per-ring states
    per_ring_states = [canonical_pucker_states(len(r))[:max_per_ring] for r in puckerable]
    from itertools import product as _prod
    for combo in _prod(*per_ring_states):
        P = coords.copy()
        for ring_atoms, (Q_t, phi_t, label) in zip(puckerable, combo):
            ring_arr = np.array([P[int(i)] for i in ring_atoms])
            new_ring = set_pucker(ring_arr, Q_t, phi_t)
            new_ring = _restore_ring_bonds(new_ring)
            for k, i in enumerate(ring_atoms):
                P[int(i)] = new_ring[k]
        yield P


def enumerate_chelate_delta_lambda(
    mol, coords: np.ndarray, max_per_ring: int = 2,
) -> Iterator[np.ndarray]:
    """Convenience: enumerate Δ/Λ (delta/lambda) chelate ring conformers only.

    Chelates 5-rings (M-N-C-C-N) and 6-rings (M-N-C-C-C-N) have two energy-minima
    fold directions. Cremer-Pople with max_per_ring=2 captures the dominant
    chair/envelope pair.

    Universal: works for any chelate ring size, any metal, any donor atom.
    """
    yield from enumerate_mol_ring_conformers(
        mol, coords, max_per_ring=max_per_ring, skip_aromatic=True, chelate_only=True
    )
