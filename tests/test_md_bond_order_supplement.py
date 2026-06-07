"""Tests for the M-D bond-order supplement library.

The supplement disaggregates M-D bonds by CCDC bond_type (single, double,
triple, aromatic) so chemically distinct populations (e.g. W=N imido vs
W-pyridine) are not pooled into a single mean.  Audit verdict 2026-06-07:
the main library v5_TM had healthy n>>1000 coverage for W/Mo/Tc/Re donor
elements but used a hyb-only key that collapsed bond-order distinctions.

The supplement is additive (env-gated) and byte-identical when the env-flag
``DELFIN_GRIP_MD_BO_SUPPLEMENT`` is unset.

Author: hmaximilian
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path

import numpy as np

SUP_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/"
    "grip_lib_md_bo_supplement_5d_mid.npz"
)
MAIN_LIB_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_TM.npz"
)


@unittest.skipUnless(SUP_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "5d-mid supplement / main library not present")
class TestMDSupplementLookup(unittest.TestCase):
    """Lookup-level: the supplement returns chemistry-distinct distances."""

    def setUp(self):
        self._old_env = os.environ.get("DELFIN_GRIP_MD_BO_SUPPLEMENT")
        os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = str(SUP_PATH)
        # Force module re-import so the supplement loads at GripLibrary init.
        import importlib
        import delfin.fffree.grip_mogul_lookup as _m
        importlib.reload(_m)
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        self.lib = GripLibrary(MAIN_LIB_PATH)

    def tearDown(self):
        if self._old_env is None:
            del os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"]
        else:
            os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = self._old_env

    def test_supplement_loaded(self):
        self.assertGreater(self.lib.n_md_bo, 100,
                           "supplement should hold hundreds of M-D keys")
        self.assertEqual(self.lib._md_bo_metals,
                         frozenset({"W", "Mo", "Tc", "Re"}))

    def test_w_n_amine_vs_imido_disaggregation(self):
        """W-N(sp3, single bond) is amine ~2.22 Å; W-N(sp2, double) is imido ~1.75 Å."""
        amine = self.lib._md_bo_lookup("W", "N", "sp3", "1", min_n=10)
        imido = self.lib._md_bo_lookup("W", "N", "sp2", "2", min_n=10)
        self.assertIsNotNone(amine, "W-N(sp3, bo=1) should resolve")
        self.assertIsNotNone(imido, "W-N(sp2, bo=2) should resolve")
        # The chemistry distinction MUST be preserved.
        self.assertGreater(amine[0] - imido[0], 0.30,
                           "W-N amine should be ≥0.30 Å longer than W=N imido")
        # Specific empirical means (with healthy n)
        self.assertAlmostEqual(amine[0], 2.22, delta=0.10)
        self.assertAlmostEqual(imido[0], 1.76, delta=0.10)

    def test_w_o_alkoxide_vs_oxo_disaggregation(self):
        """W-O(sp3,1) ≈ 1.93 Å (alkoxide); W=O(sp2,2) ≈ 1.71 Å (oxo)."""
        alkoxide = self.lib._md_bo_lookup("W", "O", "sp3", "1", min_n=10)
        oxo = self.lib._md_bo_lookup("W", "O", "sp2", "2", min_n=10)
        self.assertIsNotNone(alkoxide)
        self.assertIsNotNone(oxo)
        self.assertGreater(alkoxide[0] - oxo[0], 0.15)

    def test_mo_re_tc_coverage(self):
        """All four 5d-mid metals should expose the common donor element bonds."""
        for metal in ("Mo", "Re", "Tc"):
            hit_n = self.lib._md_bo_lookup(metal, "N", "sp3", "1", min_n=10)
            self.assertIsNotNone(hit_n, f"{metal}-N(sp3,1) should resolve")
            # Should be in the κ¹-amine range (2.05-2.30) -- realistic empirical
            self.assertGreater(hit_n[0], 2.0)
            self.assertLess(hit_n[0], 2.35)

    def test_non_supplemented_metal_returns_none(self):
        """A 3d metal (e.g. Fe) is not in this 5d-mid supplement."""
        hit = self.lib._md_bo_lookup("Fe", "N", "sp3", "1", min_n=5)
        self.assertIsNone(hit, "Fe is not in 5d-mid supplement")


@unittest.skipUnless(SUP_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "5d-mid supplement / main library not present")
class TestMDSupplementOffByteIdentical(unittest.TestCase):
    """When the env-flag is unset, the supplement must be invisible."""

    def setUp(self):
        # Snapshot/clear
        self._old_env = os.environ.pop("DELFIN_GRIP_MD_BO_SUPPLEMENT", None)
        import importlib
        import delfin.fffree.grip_mogul_lookup as _m
        importlib.reload(_m)
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        self.lib = GripLibrary(MAIN_LIB_PATH)

    def tearDown(self):
        if self._old_env is not None:
            os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = self._old_env

    def test_off_means_zero(self):
        self.assertEqual(self.lib.n_md_bo, 0)
        # Lookup with bo tag should return None when supplement is off.
        hit = self.lib._md_bo_lookup("W", "N", "sp3", "1", min_n=5)
        self.assertIsNone(hit)


@unittest.skipUnless(SUP_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "5d-mid supplement / main library not present")
class TestMDBondOrderTagDetection(unittest.TestCase):
    """RDKit-level: ``_detect_md_bond_order_tag`` returns the right tag."""

    def test_imido_double_bond_tag(self):
        from rdkit import Chem
        # F[W]=N (imido) — simplified
        smi = "F[W](F)(F)(F)=NC"
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        self.assertIsNotNone(mol)
        from delfin.fffree.mogul_bounds import _detect_md_bond_order_tag
        # Locate W and the doubly-bonded N
        W_idx = next(i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "W")
        N_idx = None
        for nb in mol.GetAtomWithIdx(W_idx).GetNeighbors():
            if nb.GetSymbol() == "N":
                N_idx = nb.GetIdx()
                break
        self.assertIsNotNone(N_idx)
        tag = _detect_md_bond_order_tag(mol, W_idx, N_idx)
        self.assertEqual(tag, "2", "W=N should be tagged '2' (double)")

    def test_aromatic_donor_tag(self):
        from rdkit import Chem
        smi = "F[W](F)(F)(F)(F)(F)n1ccccc1"
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        # Use the assemble_via_mogul helper to sanitize properly
        from delfin.fffree.assemble_via_mogul import _full_complex_mol
        mol = _full_complex_mol(smi)
        self.assertIsNotNone(mol)
        from delfin.fffree.mogul_bounds import _detect_md_bond_order_tag
        W_idx = next(i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "W")
        N_idx = None
        for nb in mol.GetAtomWithIdx(W_idx).GetNeighbors():
            if nb.GetSymbol() == "N":
                N_idx = nb.GetIdx()
                break
        self.assertIsNotNone(N_idx)
        tag = _detect_md_bond_order_tag(mol, W_idx, N_idx)
        # Either 'a' (aromatic donor) -- the chemistry-correct tag for pyridine
        self.assertEqual(tag, "a")


if __name__ == "__main__":
    unittest.main()
