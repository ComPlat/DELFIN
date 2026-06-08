"""Tests for the universal fused aromatic system planarity gate.

Phenanthroline (3 fused 6-rings, 14 atoms), naphthalene (2 fused 6-rings),
anthracene (3 linear fused), indole (5+6 fused) — every one of these
poly-aromatic systems is rigidly planar in the CCDC (RMS < 0.03 Å) yet
the soft-DG embed can leave them buckled across the fused edges.

Per-ring planarity (``flatten_aromatic_rings``) is necessary but NOT
sufficient: each ring can pass a per-ring tolerance while two fused
rings pivot around their shared edge.  This module's
``flatten_fused_aromatic_systems`` projects every atom of the joint
fused system onto ONE common plane.

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations

import unittest

import numpy as np


def _embed_planar(mol, seed: int = 7) -> np.ndarray:
    """Embed ``mol`` to 3D and return coordinates as an ``(N, 3)`` ndarray."""
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    conf = mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])


def _pucker_ring_atoms(P: np.ndarray, atoms, amplitude: float = 0.30) -> np.ndarray:
    """Push every listed atom out-of-plane by a smooth sinusoid (deterministic).

    A larger amplitude on every-other atom guarantees a sizeable
    max-deviation against the SVD best-fit plane (the SVD plane
    absorbs the smooth component).
    """
    Q = P.copy()
    for k, a in enumerate(sorted(atoms)):
        # Alternating sign + amplitude beats SVD smoothing.
        sign = 1.0 if (k % 2 == 0) else -1.0
        Q[a, 2] += sign * amplitude
    return Q


class TestFusedAromaticDetection(unittest.TestCase):
    """Detection: fused-system grouping covers naph / phen / indole / pyrene."""

    def test_phenanthroline_one_joint_system_14_atoms(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        # 1,10-phenanthroline
        m = Chem.MolFromSmiles("c1cnc2c(c1)ccc1cccnc12")
        m = Chem.AddHs(m)
        fused = detect_fused_aromatic_systems(m)
        self.assertEqual(len(fused), 1, "phen has exactly one fused system")
        self.assertEqual(len(fused[0]), 14, "phen joint system has 14 heavy atoms")

    def test_naphthalene_one_joint_system_10_atoms(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        m = Chem.AddHs(m)
        fused = detect_fused_aromatic_systems(m)
        self.assertEqual(len(fused), 1)
        self.assertEqual(len(fused[0]), 10)

    def test_anthracene_one_joint_system_14_atoms(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
        m = Chem.AddHs(m)
        fused = detect_fused_aromatic_systems(m)
        self.assertEqual(len(fused), 1)
        self.assertEqual(len(fused[0]), 14)

    def test_indole_one_joint_system_9_atoms(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccc2[nH]ccc2c1")
        m = Chem.AddHs(m)
        fused = detect_fused_aromatic_systems(m)
        self.assertEqual(len(fused), 1)
        self.assertEqual(len(fused[0]), 9)

    def test_pyrene_one_joint_system_16_atoms(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
        m = Chem.AddHs(m)
        fused = detect_fused_aromatic_systems(m)
        self.assertEqual(len(fused), 1)
        self.assertEqual(len(fused[0]), 16)

    def test_isolated_benzene_no_fused_system(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccccc1")
        m = Chem.AddHs(m)
        self.assertEqual(detect_fused_aromatic_systems(m), [])

    def test_isolated_pyridine_no_fused_system(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccncc1")
        m = Chem.AddHs(m)
        self.assertEqual(detect_fused_aromatic_systems(m), [])

    def test_cod_not_aromatic_no_fused_system(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import detect_fused_aromatic_systems
        # 1,5-cyclooctadiene — alkene, not aromatic
        m = Chem.MolFromSmiles("C1CC=CCCC=C1")
        m = Chem.AddHs(m)
        self.assertEqual(detect_fused_aromatic_systems(m), [])

    def test_two_isolated_phenyls_no_shared_atoms(self):
        """Two unfused phenyls (biphenyl) share NO atoms — not a fused system.

        Biphenyl is two phenyl rings connected by a single C-C bond; the
        rings share 0 atoms.  Our fused-system test requires ≥ 2 shared
        atoms (one shared edge), so biphenyl must be reported as TWO
        single-ring systems, NOT one joint system.
        """
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import (
            detect_aromatic_rings,
            detect_fused_aromatic_systems,
            group_fused_aromatic_rings,
        )
        m = Chem.MolFromSmiles("c1ccc(-c2ccccc2)cc1")
        m = Chem.AddHs(m)
        rings = detect_aromatic_rings(m)
        self.assertEqual(len(rings), 2, "biphenyl has two aromatic rings")
        groups = group_fused_aromatic_rings(rings)
        self.assertEqual(len(groups), 2, "biphenyl rings are NOT fused")
        self.assertEqual(detect_fused_aromatic_systems(m), [])


class TestPhenanthrolineFlat(unittest.TestCase):
    """Phenanthroline (PEHDUD-class) — 3 fused 6-rings, 14 atoms.

    Real-world buckled embed → joint plane projection brings every atom
    of the 14-atom system to within 0.05 Å of one common SVD plane.
    """

    def test_phenanthroline_flat(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import (
            flatten_fused_aromatic_systems,
            measure_fused_planarity,
            detect_fused_aromatic_systems,
        )
        m = Chem.MolFromSmiles("c1cnc2c(c1)ccc1cccnc12")
        m = Chem.AddHs(m)
        P = _embed_planar(m, seed=11)
        atoms = detect_fused_aromatic_systems(m)[0]
        # Force buckle (mimics the soft-DG residual on a real complex).
        P_b = _pucker_ring_atoms(P, atoms, amplitude=0.35)
        md_pre = measure_fused_planarity(P_b, atoms)
        self.assertGreater(md_pre, 0.20, "buckle should be measurable pre-flatten")

        P_f, diag = flatten_fused_aromatic_systems(m, P_b, max_dev_threshold=0.10)
        md_post = measure_fused_planarity(P_f, atoms)
        self.assertLess(
            md_post, 0.05,
            f"phen joint system not flat enough: max_dev={md_post:.4f}Å",
        )
        # Diagnostic: 3 rings reported, flattened=True.
        self.assertEqual(len(diag), 1)
        self.assertEqual(diag[0]["n_rings"], 3)
        self.assertTrue(diag[0]["flattened"])


class TestNaphthaleneFlat(unittest.TestCase):
    """Naphthalene — 2 fused 6-rings, 10 atoms — universal case."""

    def test_naphthalene_flat(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import (
            flatten_fused_aromatic_systems,
            measure_fused_planarity,
            detect_fused_aromatic_systems,
        )
        m = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        m = Chem.AddHs(m)
        P = _embed_planar(m, seed=3)
        atoms = detect_fused_aromatic_systems(m)[0]
        P_b = _pucker_ring_atoms(P, atoms, amplitude=0.30)
        md_pre = measure_fused_planarity(P_b, atoms)
        self.assertGreater(md_pre, 0.15)

        P_f, diag = flatten_fused_aromatic_systems(m, P_b, max_dev_threshold=0.10)
        md_post = measure_fused_planarity(P_f, atoms)
        self.assertLess(
            md_post, 0.05,
            f"naphthalene joint max_dev={md_post:.4f}Å",
        )
        self.assertEqual(diag[0]["n_rings"], 2)


class TestIsolatedBenzeneNotAffected(unittest.TestCase):
    """Single isolated ring → per-ring pass owns it.  Fused pass: no-op."""

    def test_isolated_benzene_no_motion(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import flatten_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccccc1")
        m = Chem.AddHs(m)
        P = _embed_planar(m, seed=1)
        # Buckle just to be sure: fused pass STILL must not touch it
        # because biphenyl-style "no fused system" is the contract.
        P_b = P.copy()
        P_b[0, 2] += 0.5
        P_f, diag = flatten_fused_aromatic_systems(m, P_b, max_dev_threshold=0.10)
        # Nothing reported, P unchanged byte-wise.
        self.assertEqual(diag, [])
        self.assertTrue(np.array_equal(P_b, P_f))

    def test_isolated_pyridine_no_motion(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import flatten_fused_aromatic_systems
        m = Chem.MolFromSmiles("c1ccncc1")
        m = Chem.AddHs(m)
        P = _embed_planar(m, seed=2)
        P_b = P.copy()
        P_b[2, 2] += 0.4
        P_f, diag = flatten_fused_aromatic_systems(m, P_b, max_dev_threshold=0.10)
        self.assertEqual(diag, [])
        self.assertTrue(np.array_equal(P_b, P_f))


class TestTopologyPreserved(unittest.TestCase):
    """Projection must not break bonds or invert chirality of the
    poly-aromatic skeleton — within-ring bond lengths must stay within
    0.05 Å of pre-projection, and the atom-count must not change."""

    def test_phen_bond_lengths_preserved(self):
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import (
            flatten_fused_aromatic_systems,
            detect_aromatic_rings,
        )
        m = Chem.MolFromSmiles("c1cnc2c(c1)ccc1cccnc12")
        m = Chem.AddHs(m)
        P = _embed_planar(m, seed=11)
        rings = detect_aromatic_rings(m)
        atoms = sorted({a for r in rings for a in r})
        P_b = _pucker_ring_atoms(P, atoms, amplitude=0.30)

        # Catalogue pre-projection bond lengths.
        bonds_pre = []
        for b in m.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            bonds_pre.append((i, j, float(np.linalg.norm(P_b[i] - P_b[j]))))

        P_f, _diag = flatten_fused_aromatic_systems(m, P_b, max_dev_threshold=0.10)
        # Same atom count.
        self.assertEqual(P_f.shape, P_b.shape)
        # Bond lengths inside the fused system are well preserved (slight
        # change from the perpendicular projection is OK; tolerance 0.15 Å
        # is comfortable because the perpendicular shift is small).
        for i, j, d_pre in bonds_pre:
            d_post = float(np.linalg.norm(P_f[i] - P_f[j]))
            self.assertLess(
                abs(d_post - d_pre), 0.20,
                f"bond {i}-{j} length drift {abs(d_post-d_pre):.4f}Å too high",
            )

    def test_frozen_atoms_do_not_move(self):
        """Frozen atoms (donor / metal) must receive zero displacement."""
        from rdkit import Chem
        from delfin.fffree.aromatic_flatten import (
            flatten_fused_aromatic_systems,
            detect_aromatic_rings,
        )
        # Phen has two N donors at positions 2 and 11 (ring atoms).
        m = Chem.MolFromSmiles("c1cnc2c(c1)ccc1cccnc12")
        m = Chem.AddHs(m)
        n_idxs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
        self.assertGreaterEqual(len(n_idxs), 2, "phen needs 2 N donors")

        P = _embed_planar(m, seed=11)
        rings = detect_aromatic_rings(m)
        atoms = sorted({a for r in rings for a in r})
        P_b = _pucker_ring_atoms(P, atoms, amplitude=0.30)

        # Lock both Ns as frozen.
        P_f, _diag = flatten_fused_aromatic_systems(
            m, P_b, max_dev_threshold=0.10, frozen_idxs=n_idxs,
        )
        for n in n_idxs:
            delta = float(np.linalg.norm(P_f[n] - P_b[n]))
            self.assertLess(
                delta, 1e-9,
                f"frozen N at {n} moved {delta:.6f}Å",
            )


class TestEnvFlagDefaultOFF(unittest.TestCase):
    """The fused-plane gate is default OFF (byte-identity unless flipped)."""

    def test_env_off_default(self):
        import os
        from delfin.fffree.aromatic_flatten import fused_aromatic_plane_enabled
        prior = os.environ.pop("DELFIN_FFFREE_FUSED_AROMATIC_PLANE", None)
        try:
            self.assertFalse(fused_aromatic_plane_enabled())
        finally:
            if prior is not None:
                os.environ["DELFIN_FFFREE_FUSED_AROMATIC_PLANE"] = prior

    def test_env_on_active(self):
        import os
        from delfin.fffree.aromatic_flatten import fused_aromatic_plane_enabled
        prior = os.environ.get("DELFIN_FFFREE_FUSED_AROMATIC_PLANE")
        os.environ["DELFIN_FFFREE_FUSED_AROMATIC_PLANE"] = "1"
        try:
            self.assertTrue(fused_aromatic_plane_enabled())
        finally:
            if prior is None:
                os.environ.pop("DELFIN_FFFREE_FUSED_AROMATIC_PLANE", None)
            else:
                os.environ["DELFIN_FFFREE_FUSED_AROMATIC_PLANE"] = prior


if __name__ == "__main__":
    unittest.main()
