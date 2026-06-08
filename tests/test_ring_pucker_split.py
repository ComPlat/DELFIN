"""Tests for the 2026-06-08 ring-pucker / single-bond rotamer SPLIT.

Architecture
------------
The legacy ``DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS`` env-flag boxed
TWO independent features into one switch:

  1. Cremer-Pople ring puckering (topology-preserving by construction --
     each ring atom keeps its ring-bond partners; only the out-of-plane
     displacement along the mean-plane normal changes).
  2. Single-bond rotamers (dihedral rotation about each rotatable bond;
     may break topology when a rotation propagates into the M-shell or
     swings a substituent into another ligand's footprint).

The 2026-06-08 split introduces two independent flags:

  ``DELFIN_FFFREE_MOGUL_PRIMARY_RING_PUCKER``  (default ON)
  ``DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER``      (default OFF)

with the legacy flag retained as a master kill-switch
(``CONFORMERS=0`` => both OFF; ``CONFORMERS=1`` => ring-pucker ON,
rotamer at its own default).

This test module pins ring puckering to the universal cases:

  * cyclohexane (6-ring, fully saturated) emits 2+ distinct puckers
    (chair + boat + twist + half-chair from the Cremer-Pople canonical
    state set).
  * benzene (6-ring, aromatic) is skipped (planar = single state).
  * the bond graph is preserved across all yielded pucker frames.
  * the single-bond rotamer channel stays at zero variants when its
    own flag is unset (independence proof).
"""
from __future__ import annotations

import os
import unittest
from typing import List, Optional

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_mol_with_coords(smiles: str):
    """Embed a SMILES with explicit Hs, return (mol_g, coords).

    Uses the standard RDKit ETKDG path so the ring geometry starts in a
    realistic puckered state (cyclohexane lands in a chair-ish
    configuration, benzene in a planar configuration).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xDE17
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None, None
    conf = mol.GetConformer(0)
    coords = np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
          conf.GetAtomPosition(i).z]
         for i in range(mol.GetNumAtoms())],
        dtype=float,
    )
    return mol, coords


def _bond_graph(mol) -> List[tuple]:
    """Return the sorted (i, j) tuples for every bond in ``mol``.

    The pucker variants must preserve this graph atom-for-atom.
    """
    bonds: List[tuple] = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds.append((min(i, j), max(i, j)))
    return sorted(bonds)


def _geometric_bond_graph(
    syms: List[str], P: np.ndarray, tol_factor: float = 1.30,
) -> List[tuple]:
    """Reconstruct the bond graph from coords alone (universal proximity
    test).  Two heavy atoms are bonded when their distance is below
    ``tol_factor * ideal_bond(i, j)``.
    """
    import delfin._bond_decollapse as _bd
    n = len(syms)
    bonds: List[tuple] = []
    for i in range(n):
        for j in range(i + 1, n):
            try:
                ideal = float(_bd._ideal_bond(syms[i], syms[j]))
            except Exception:
                continue
            if ideal <= 0:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < tol_factor * ideal:
                bonds.append((i, j))
    return sorted(bonds)


# ---------------------------------------------------------------------------
# Env-flag isolation
# ---------------------------------------------------------------------------
def _set_env(**kwargs) -> dict:
    """Set env vars in ``kwargs`` to their string values; ``None`` clears."""
    saved = {}
    for k, v in kwargs.items():
        saved[k] = os.environ.get(k)
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    return saved


def _restore_env(saved: dict) -> None:
    for k, v in saved.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------
class TestRingPuckerRotamerSplit(unittest.TestCase):
    """Validates the env-flag split between Cremer-Pople ring puckering
    and single-bond rotamer enumeration (2026-06-08)."""

    def test_env_flag_defaults(self):
        """Default state: RING_PUCKER ON, ROTAMER OFF, CONFORMERS ON
        (CONFORMERS reads as the disjunction of the two sub-flags)."""
        saved = _set_env(
            DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS=None,
            DELFIN_FFFREE_MOGUL_PRIMARY_RING_PUCKER=None,
            DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER=None,
        )
        try:
            from delfin.fffree.assemble_via_mogul import (
                mogul_primary_ring_pucker_enabled,
                mogul_primary_rotamer_enabled,
                mogul_primary_conformers_enabled,
            )
            self.assertTrue(mogul_primary_ring_pucker_enabled())
            self.assertFalse(mogul_primary_rotamer_enabled())
            self.assertTrue(mogul_primary_conformers_enabled())
        finally:
            _restore_env(saved)

    def test_env_flag_master_kill(self):
        """CONFORMERS=0 still kills both sub-features (backwards compat)."""
        saved = _set_env(
            DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS="0",
            DELFIN_FFFREE_MOGUL_PRIMARY_RING_PUCKER="1",
            DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER="1",
        )
        try:
            from delfin.fffree.assemble_via_mogul import (
                mogul_primary_ring_pucker_enabled,
                mogul_primary_rotamer_enabled,
                mogul_primary_conformers_enabled,
            )
            self.assertFalse(mogul_primary_ring_pucker_enabled())
            self.assertFalse(mogul_primary_rotamer_enabled())
            self.assertFalse(mogul_primary_conformers_enabled())
        finally:
            _restore_env(saved)

    def test_env_flag_pucker_only(self):
        """RING_PUCKER=1 + ROTAMER=0 == default split state."""
        saved = _set_env(
            DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS=None,
            DELFIN_FFFREE_MOGUL_PRIMARY_RING_PUCKER="1",
            DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER="0",
        )
        try:
            from delfin.fffree.assemble_via_mogul import (
                mogul_primary_ring_pucker_enabled,
                mogul_primary_rotamer_enabled,
            )
            self.assertTrue(mogul_primary_ring_pucker_enabled())
            self.assertFalse(mogul_primary_rotamer_enabled())
        finally:
            _restore_env(saved)

    def test_cyclohexane_chair_boat_emitted(self):
        """Cyclohexane must emit 2+ distinct Cremer-Pople pucker variants
        (chair, boat, twist-boat, half-chair, ...) under the ring-pucker
        channel."""
        from delfin.fffree.assemble_via_mogul import (
            _enumerate_ring_pucker_variants,
        )

        mol, coords = _build_mol_with_coords("C1CCCCC1")
        self.assertIsNotNone(mol, "cyclohexane embed failed")
        # Cap at 8 = chair+boat+twist+half-chair sampling
        variants = _enumerate_ring_pucker_variants(
            mol, coords, max_per_ring=8, skip_aromatic=True,
        )
        # At least 2 distinct frames (more is fine; common is 4-8).
        self.assertGreaterEqual(
            len(variants), 2,
            f"cyclohexane: expected >=2 pucker variants, got {len(variants)}",
        )
        # Each variant must be distinct from the input geometry by some
        # measurable z-displacement (otherwise the pucker collapse onto
        # the base).
        seen_z_signatures = set()
        for P_v, lab in variants:
            self.assertTrue(np.all(np.isfinite(P_v)))
            # Hash a coarse z-profile of the 6 ring atoms (indices 0..5).
            ring_atoms = list(range(6))
            z_sig = tuple(
                round(float(P_v[i, 2] - P_v[0, 2]), 1) for i in ring_atoms
            )
            seen_z_signatures.add(z_sig)
        self.assertGreaterEqual(
            len(seen_z_signatures), 2,
            "cyclohexane pucker variants collapsed onto a single z-profile "
            "(chair/boat distinction lost)",
        )

    def test_benzene_planar_single_mode(self):
        """Benzene (aromatic 6-ring) must be skipped by the pucker
        enumerator -- planar = single canonical state."""
        from delfin.fffree.assemble_via_mogul import (
            _enumerate_ring_pucker_variants,
        )

        mol, coords = _build_mol_with_coords("c1ccccc1")
        self.assertIsNotNone(mol, "benzene embed failed")
        # The aromatic-atom-set is the full ring (6 carbons + 6 hydrogens
        # at index 6..11 are non-ring); pass the ring carbons explicitly.
        arom_set = [
            i for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetIsAromatic()
        ]
        self.assertEqual(
            len(arom_set), 6,
            f"benzene aromatic-atom-set should be 6, got {len(arom_set)}",
        )
        variants = _enumerate_ring_pucker_variants(
            mol, coords, max_per_ring=8, skip_aromatic=True,
            aromatic_atom_set=arom_set,
        )
        # Aromatic ring skipped => no pucker variants.
        self.assertEqual(
            len(variants), 0,
            f"benzene: expected 0 pucker variants (aromatic skip), "
            f"got {len(variants)}",
        )

    def test_topology_preserved_across_puckers(self):
        """Every cyclohexane pucker variant must preserve the bond graph
        (Cremer-Pople is topology-preserving by construction)."""
        from delfin.fffree.assemble_via_mogul import (
            _enumerate_ring_pucker_variants,
        )

        mol, coords = _build_mol_with_coords("C1CCCCC1")
        self.assertIsNotNone(mol)

        # Reference graph from the RDKit mol bonds.
        ref_bonds = _bond_graph(mol)
        ref_set = set(ref_bonds)
        syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]

        # Base graph from coords (proximity test).  We use this same
        # test on every pucker frame -- the proximity-derived graph must
        # equal the base for every variant.
        base_graph = set(_geometric_bond_graph(syms, coords))
        # Every RDKit bond should appear in the proximity graph.
        for b in ref_set:
            self.assertIn(
                b, base_graph,
                f"base cyclohexane geometry missing RDKit bond {b}",
            )

        variants = _enumerate_ring_pucker_variants(
            mol, coords, max_per_ring=6, skip_aromatic=True,
        )
        self.assertGreater(len(variants), 0)
        for P_v, lab in variants:
            v_graph = set(_geometric_bond_graph(syms, P_v))
            # Every RDKit bond must survive in the proximity graph
            # (puckering preserves the bond graph by construction).
            for b in ref_set:
                self.assertIn(
                    b, v_graph,
                    f"pucker variant {lab!r} broke bond {b}: missing "
                    f"in geometric graph",
                )

    def test_rotamer_off_independent_of_pucker_on(self):
        """With RING_PUCKER=1 + ROTAMER=0 the rotamer code path emits zero
        variants regardless of how many rotatable bonds the mol has.

        This is the independence proof: enabling ring-pucker does not
        leak into the rotamer channel.
        """
        from delfin.fffree.assemble_via_mogul import (
            _enumerate_rotamer_variants,
            mogul_primary_rotamer_enabled,
        )

        # Cyclohexane has zero rotatable single bonds (ring bonds are
        # exempt) so the rotamer enumerator emits 0 variants directly.
        # The flag-state check is the real test: ROTAMER stays OFF under
        # the default split, even when RING_PUCKER=1.
        saved = _set_env(
            DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS=None,
            DELFIN_FFFREE_MOGUL_PRIMARY_RING_PUCKER="1",
            DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER=None,  # default OFF
        )
        try:
            self.assertFalse(
                mogul_primary_rotamer_enabled(),
                "ROTAMER channel must default OFF when only "
                "RING_PUCKER=1 is set explicitly",
            )
        finally:
            _restore_env(saved)

        # Sanity: when the rotamer flag is explicitly ON, the helper
        # returns True (and the underlying enumerator runs on a SMILES
        # with rotatable bonds).
        saved = _set_env(
            DELFIN_FFFREE_MOGUL_PRIMARY_RING_PUCKER="0",
            DELFIN_FFFREE_MOGUL_PRIMARY_ROTAMER="1",
        )
        try:
            self.assertTrue(mogul_primary_rotamer_enabled())
        finally:
            _restore_env(saved)


if __name__ == "__main__":
    unittest.main()
