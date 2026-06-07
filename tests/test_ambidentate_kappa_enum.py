"""Tests for the universal ambidentate κⁿ binding-mode enumeration.

Validates:

1.  ``test_carboxylate_emits_kappa2`` — acetate (CC(=O)[O]M) emits both
    κ¹-O (each oxygen as linkage isomer) AND κ²-OO bidentate variants.
2.  ``test_nitrate_emits_kappa3`` — nitrate (M-ONO₂) emits κ¹-O variants
    (one per oxygen, deduped by canonical rank) plus κ²-OO bidentate
    plus κ³-OOO tridentate variants.
3.  ``test_thiocyanate_linkage_isomers`` — SCN⁻ emits both κ-S (seed)
    and κ-N (linkage isomer) plus κ²-NS chelate.
4.  ``test_pyridine_no_extras`` — monodentate pyridine has only one donor
    so no κⁿ variants beyond the seed.
5.  ``test_off_byte_identical`` — when the env-flag is unset, the
    enumerator emits only the seed κ¹ structure (byte-identical with
    pre-κⁿ code path).
6.  ``test_co_linkage_only`` — CO ligand emits κ¹-C (seed) + κ¹-O
    (isocarbonyl linkage isomer), no κ² (C and O too close).
7.  ``test_universal_no_smarts`` — verifies the module never matches via
    SMARTS patterns (read the source).
"""
from __future__ import annotations

import os
import unittest

from rdkit import Chem

import delfin._bond_decollapse as _bd


def _make_mol(smi: str):
    """Prepare a mol for κⁿ enumeration (mirrors `assemble_via_mogul`)."""
    try:
        from delfin.smiles_converter import _prepare_mol_for_embedding
        mol = _prepare_mol_for_embedding(smi)
    except Exception:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return None
    if not any(a.GetSymbol() == "H" for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol)
    try:
        Chem.FastFindRings(mol)
    except Exception:
        pass
    return mol


def _metal_and_donors(mol):
    metals = [a.GetIdx() for a in mol.GetAtoms()
              if _bd._is_metal(a.GetSymbol())]
    if not metals:
        return None, []
    m = metals[0]
    donors = [int(n.GetIdx()) for n in mol.GetAtomWithIdx(m).GetNeighbors()]
    return int(m), donors


class TestUniversalKappaEnumeration(unittest.TestCase):
    """Universal κⁿ binding mode detection (graph + Lewis octet only)."""

    def setUp(self):
        # κⁿ enumeration is gated by DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM;
        # tests run with the flag set so the module emits beyond the seed.
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM", None)

    def _modes_for(self, smi: str):
        from delfin.fffree.ambidentate_kappa_enum import (
            enumerate_kappa_modes_for_ligand,
            enumerate_complex_kappa_variants,
        )
        mol = _make_mol(smi)
        self.assertIsNotNone(mol, f"failed to parse {smi!r}")
        m, donors = _metal_and_donors(mol)
        self.assertIsNotNone(m, f"no metal found in {smi!r}")
        self.assertGreater(len(donors), 0, f"no donors in {smi!r}")
        per_ligand = {}
        for seed in donors:
            per_ligand[seed] = enumerate_kappa_modes_for_ligand(mol, m, seed)
        full = enumerate_complex_kappa_variants(mol, m, donors)
        return mol, m, donors, per_ligand, full

    def test_carboxylate_emits_kappa2(self):
        """SIYMEU-class: Ag-acetate emits κ²-OO bidentate option."""
        smi = "CC(=O)[O][Ag][NH3]"
        mol, m, donors, per_ligand, full = self._modes_for(smi)
        # Find the O-seed ligand (acetate, not NH3).
        o_seed = next(d for d in donors
                      if mol.GetAtomWithIdx(d).GetSymbol() == "O")
        modes = per_ligand[o_seed]
        kappas = sorted({mode["kappa"] for mode in modes})
        self.assertEqual(kappas, [1, 2],
                         f"acetate should expose κ1 + κ2; got {kappas}")
        # κ² mode must have both donor elements = O.
        k2 = next(m_ for m_ in modes if m_["kappa"] == 2)
        self.assertEqual(k2["donor_elements"], ("O", "O"),
                         f"κ² donors should be (O, O); got {k2['donor_elements']}")
        # Full-complex variants should include κ²-OO variant.
        k2_full = [v for v in full
                   if not v["is_seed_default"] and "k2-OO" in v["label_suffix"]]
        self.assertGreater(len(k2_full), 0,
                           "no κ²-OO full-complex variant emitted")

    def test_nitrate_emits_kappa3(self):
        """Nitrate (NO3⁻) emits κ¹-O / κ²-OO / κ³-OOO denticities."""
        smi = "[NH3][Ni][O][N+](=O)[O-]"
        mol, m, donors, per_ligand, full = self._modes_for(smi)
        # Find the O-seed (nitrate-O, the one bonded to Ni in the input).
        o_seed = next(d for d in donors
                      if mol.GetAtomWithIdx(d).GetSymbol() == "O")
        modes = per_ligand[o_seed]
        kappas = sorted({mode["kappa"] for mode in modes})
        # Must contain κ1, κ2, and κ3.
        self.assertIn(1, kappas)
        self.assertIn(2, kappas)
        self.assertIn(3, kappas)

    def test_thiocyanate_linkage_isomers(self):
        """SCN⁻ on Ni: emits κ¹-S (seed) + κ¹-N (linkage isomer) + κ²-NS."""
        smi = "[NH3][Ni][S]C#N"
        mol, m, donors, per_ligand, full = self._modes_for(smi)
        s_seed = next(d for d in donors
                      if mol.GetAtomWithIdx(d).GetSymbol() == "S")
        modes = per_ligand[s_seed]
        # Distinct donor element sets across the κⁿ modes.
        elem_sets = {tuple(sorted(m_["donor_elements"])) for m_ in modes}
        self.assertIn(("S",), elem_sets,
                      f"thiocyanate κ¹-S missing; got {elem_sets}")
        self.assertIn(("N",), elem_sets,
                      f"thiocyanate κ¹-N linkage isomer missing; got {elem_sets}")
        self.assertIn(("N", "S"), elem_sets,
                      f"thiocyanate κ²-NS missing; got {elem_sets}")

    def test_co_linkage_only(self):
        """CO ligand: κ¹-C (seed) + κ¹-O (isocarbonyl); NO κ² (too close)."""
        smi = "[O+]#[C-][Ni][NH3]"
        mol, m, donors, per_ligand, full = self._modes_for(smi)
        c_seed = next(d for d in donors
                      if mol.GetAtomWithIdx(d).GetSymbol() == "C")
        modes = per_ligand[c_seed]
        kappas = sorted({mode["kappa"] for mode in modes})
        # CO is 2 atoms only -> κ¹-C and κ¹-O linkage isomers.  No κ²
        # because the donors are bonded directly (graph distance 1)
        # which fails the bite-compatibility test.
        self.assertEqual(kappas, [1],
                         f"CO should only expose κ¹; got {kappas}")
        elem_sets = {tuple(sorted(m_["donor_elements"])) for m_ in modes}
        self.assertIn(("C",), elem_sets)
        self.assertIn(("O",), elem_sets)

    def test_pyridine_no_extras(self):
        """Pyridine has one N donor; no κⁿ alternates."""
        smi = "[NH3][Ni][N]1=CC=CC=C1"
        mol, m, donors, per_ligand, full = self._modes_for(smi)
        # Find the pyridine N seed.
        py_seed = None
        for d in donors:
            atom = mol.GetAtomWithIdx(d)
            if atom.GetSymbol() != "N":
                continue
            # Pyridine N: aromatic ring neighbour.
            if any(nb.GetIsAromatic() or nb.GetSymbol() == "C"
                   for nb in atom.GetNeighbors()):
                # The NH3 has only H neighbours plus Ni; pyridine N has C
                # neighbours.  Skip NH3 (H neighbours).
                if any(nb.GetSymbol() == "C" for nb in atom.GetNeighbors()):
                    py_seed = d
                    break
        self.assertIsNotNone(py_seed, "pyridine N seed not found")
        modes = per_ligand[py_seed]
        # Only the seed κ¹-N -- no other potential donors in the pyridine ring
        # are coordination-capable (the other ring atoms are C with no lone
        # pair, plus a remote N if any -- but plain pyridine has only one N).
        self.assertEqual(len(modes), 1,
                         f"pyridine should expose exactly 1 mode (seed κ¹); "
                         f"got {len(modes)}")
        self.assertTrue(modes[0]["is_seed_only"])

    def test_full_complex_default_always_first(self):
        """The κⁿ enumerator emits the seed-κ¹ default variant first."""
        smi = "CC(=O)[O][Ag][NH3]"
        _, _, _, _, full = self._modes_for(smi)
        self.assertGreaterEqual(len(full), 1)
        # First variant must be the canonical κ¹-seed (no rewiring).
        self.assertTrue(full[0]["is_seed_default"],
                        f"first variant should be seed-default; got {full[0]}")


class TestKappaEnumOffPath(unittest.TestCase):
    """When the κⁿ env-flag is OFF, the enumerator emits only the seed
    κ¹ structure (byte-identical with the pre-κⁿ code path).
    """

    def setUp(self):
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM", None)

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM", None)

    def test_off_seed_only(self):
        from delfin.fffree.ambidentate_kappa_enum import (
            enumerate_kappa_modes_for_ligand,
        )
        smi = "CC(=O)[O][Ag][NH3]"
        mol = _make_mol(smi)
        m, donors = _metal_and_donors(mol)
        o_seed = next(d for d in donors
                      if mol.GetAtomWithIdx(d).GetSymbol() == "O")
        modes = enumerate_kappa_modes_for_ligand(mol, m, o_seed)
        # OFF path: only the seed κ¹.
        self.assertEqual(len(modes), 1,
                         f"OFF path should emit only seed κ¹; got {len(modes)}")
        self.assertTrue(modes[0]["is_seed_only"])


class TestUniversalityChecks(unittest.TestCase):
    """The module must not match SMARTS patterns or per-anion templates."""

    def test_no_smarts_in_source(self):
        """Read the module source: ZERO SMARTS strings allowed."""
        import delfin.fffree.ambidentate_kappa_enum as _M
        src = open(_M.__file__).read()
        # SMARTS patterns are RDKit `MolFromSmarts` calls.
        self.assertNotIn("MolFromSmarts", src,
                         "ambidentate_kappa_enum must not use SMARTS")
        # Per-anion SMILES patterns are also forbidden (`"[S]C#N"` in code).
        # The module CAN contain SMILES in docstrings + tests; we only ban
        # the function-body matchers.  Heuristic: a raw `[S]C#N` string outside
        # the docstring would signal templating.  Skip if it appears only in
        # comments / docstrings.
        # We just assert MolFromSmarts is absent.

    def test_no_per_anion_dict(self):
        """No AMBIDENTATE_LIGANDS-style per-anion lookup dict in module."""
        import delfin.fffree.ambidentate_kappa_enum as _M
        # The universal module exposes find_potential_donors_in_ligand and
        # enumerate_kappa_modes_for_ligand; it must not expose a per-anion
        # template dictionary like the legacy linkage_isomers module does.
        forbidden = {"AMBIDENTATE_LIGANDS", "AMBIDENTATE_PATTERNS"}
        for name in forbidden:
            self.assertFalse(hasattr(_M, name),
                             f"universal module must not expose {name}")


if __name__ == "__main__":
    unittest.main()
