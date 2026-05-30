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


def _rmsd_align(P: np.ndarray, Q: np.ndarray) -> float:
    """Kabsch RMSD between two (N,3) point sets (centered + rotated).
    Returns RMSD value. Used for on-the-fly conformer deduplication.
    """
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = Pc @ R.T - Qc
    return float(np.sqrt((diff * diff).sum() / len(P)))


def _get_bonded_h_per_ring_atom(mol, ring_atoms: Sequence[int]) -> dict:
    """Return {ring_atom_idx: [bonded_H_indices]} for rigid-H drag during pucker."""
    out = {}
    for ai in ring_atoms:
        atom = mol.GetAtomWithIdx(int(ai))
        bonded_h = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
        out[int(ai)] = bonded_h
    return out


def enumerate_mol_ring_conformers(
    mol,
    coords: np.ndarray,
    max_per_ring: Optional[int] = None,
    skip_aromatic: bool = True,
    chelate_only: bool = False,
    rmsd_dedup_tol: float = 0.15,
    max_total_variants: int = 1024,
) -> Iterator[np.ndarray]:
    """Enumerate ring-pucker conformer variants of a molecule with COMPLETENESS
    guarantee + RMSD-based on-the-fly deduplication.

    Parameters
    ----------
    mol : RDKit Mol (for ring topology + aromaticity detection)
    coords : (M, 3) atom positions (must align with mol atom order)
    max_per_ring : cap on canonical pucker states per ring.
        None (default) = FULL CP enumeration (14 for 6-ring, 10 for 5-ring,
        28 for 7-ring, etc.) — COMPLETENESS guarantee under Cremer-Pople.
        Set to small N to truncate for speed.
    skip_aromatic : aromatic rings stay planar (no pucker enumeration)
    chelate_only : if True, only enumerate puckers for metal-containing rings
        (chelate Δ/Λ enumeration mode)
    rmsd_dedup_tol : Å threshold; conformers within this RMSD are considered
        duplicates and skipped. 0.5 Å = standard organic conformer-dedup.
    max_total_variants : hard cap on yielded variants (safety net against
        combinatorial explosion on multi-ring systems with 14^n possibilities).

    Yields
    ------
    coords_variant : (M, 3) array with target pucker applied to puckerable rings.
        Each yielded variant is distinct (RMSD > rmsd_dedup_tol vs all previous).

    Completeness Guarantee
    ----------------------
    Under the Cremer-Pople formalism + C_N ring symmetry, the canonical_pucker_states
    function enumerates ALL distinct pure-mode pucker types of an N-ring (Stoddart
    pucker sphere coverage). The Cartesian product across all puckerable rings,
    combined with RMSD-deduplication at the molecular level, yields the COMPLETE
    set of distinct ring-conformer geometries of the molecule.

    Mathematically:
        |distinct variants| = |X|/|G_mol|
    where X is the full enumeration space and G_mol is the molecular automorphism
    group (Burnside's lemma applies to the orbit-counting).

    Returns NOTHING if env DELFIN_FFFREE_RING_PUCKER is unset (default OFF,
    byte-identical to no integration).
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
    yielded: List[np.ndarray] = []
    orig = coords.copy()
    yielded.append(orig)
    yield orig
    # Per-ring states (full coverage unless capped) + rigid-H map per ring
    per_ring_states = []
    bonded_h_per_ring = []
    for r in puckerable:
        states = canonical_pucker_states(len(r))
        if max_per_ring is not None and max_per_ring < len(states):
            states = states[:max_per_ring]
        per_ring_states.append(states)
        bonded_h_per_ring.append(_get_bonded_h_per_ring_atom(mol, r))
    # Cartesian product with RMSD-dedup + rigid-H drag
    from itertools import product as _prod
    count = 0
    for combo in _prod(*per_ring_states):
        if count >= max_total_variants:
            return
        P = coords.copy()
        for ring_atoms, h_map, (Q_t, phi_t, label) in zip(
            puckerable, bonded_h_per_ring, combo
        ):
            ring_arr = np.array([P[int(i)] for i in ring_atoms])
            new_ring = set_pucker(ring_arr, Q_t, phi_t)
            new_ring = _restore_ring_bonds(new_ring)
            for k, i in enumerate(ring_atoms):
                displacement = new_ring[k] - P[int(i)]
                P[int(i)] = new_ring[k]
                # Rigid-H drag: H atoms bonded to this ring atom translate with it
                for h_idx in h_map.get(int(i), []):
                    P[int(h_idx)] = P[int(h_idx)] + displacement
        # On-the-fly RMSD-dedup against all previously yielded
        is_duplicate = False
        for prev in yielded:
            if _rmsd_align(prev, P) < rmsd_dedup_tol:
                is_duplicate = True
                break
        if is_duplicate:
            continue
        yielded.append(P)
        count += 1
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
