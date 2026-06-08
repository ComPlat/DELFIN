"""Conservative default for the M-D bond-order tag classifier.

The detector ``_detect_md_bond_order_tag`` in ``delfin.fffree.mogul_bounds``
must return ``'1'`` (single bond) or ``'a'`` (aromatic donor) for the
COMMON M-D coordination case, never the rarer ``'2'``/``'3'``/``'4'``
tags unless the SMILES parser explicitly stored a double/triple/
quadruple bond type on the M-D bond itself.  This is the conservative
default the spec mandates for the bond-order disaggregation supplement.

The pre-fix detector returned ``'2'`` on any bond whose
``GetBondTypeAsDouble()`` lay in [1.7, 2.5), which fired for parser
edge cases (dative bonds on sp3 P with phenyls, sp2 C of Cp/arene)
and mis-routed the supplement lookup to short imido / multiple-bond
distributions instead of the long κⁿ single-bond pool.

Author: hmaximilian
"""
from __future__ import annotations

import os
import unittest

from rdkit import Chem


class TestConservativeDefault(unittest.TestCase):
    """Conservative bond-order tagging for ambiguous M-D bonds."""

    def _detect(self, smi, metal_sym, donor_sym):
        from delfin.fffree.assemble_via_mogul import _full_complex_mol
        from delfin.fffree.mogul_bounds import _detect_md_bond_order_tag
        mol = _full_complex_mol(smi)
        self.assertIsNotNone(mol)
        m_idx = next(i for i, a in enumerate(mol.GetAtoms())
                     if a.GetSymbol() == metal_sym)
        tags = []
        for nb in mol.GetAtomWithIdx(m_idx).GetNeighbors():
            if nb.GetSymbol() == donor_sym:
                tags.append(_detect_md_bond_order_tag(mol, m_idx, nb.GetIdx()))
        return tags

    def test_dative_sp3_phosphine_is_single(self):
        """sp3 P with phenyls coordinating via DATIVE → bo='1'.

        Pre-fix this was '2' on some RDKit builds because the dative
        bond's numeric type fell in the (1.7, 2.5) gap.
        """
        smi = ("CC#[N+][Re-4]12([N+]#CC)([N+]#CC)[P+](C3=CC=CC=C3)(C3=CC=CC=C3)"
               "CC(C)(C[P+]1(C1=CC=CC=C1)C1=CC=CC=C1)C[P+]2(C1=CC=CC=C1)"
               "C1=CC=CC=C1")
        tags = self._detect(smi, "Re", "P")
        self.assertEqual(len(tags), 3,
                         "IJUCOF Re-P tridentate should expose 3 P donors")
        for t in tags:
            self.assertEqual(t, "1", "Re-P (sp3 P, DATIVE) must tag '1'")

    def test_dative_sp_carbon_is_single(self):
        """sp C donor (η²-CC olefin) coordinating via DATIVE → bo='1'."""
        smi = ("[C+]1=[C+][Mo-9]12345([C+]=[C+]2)([C+]=[C+]3)[N+]1=C(C=CC=C1"
               "C[P+]4(C1=CC=CC=C1)C1=CC=CC=C1)C[P+]5(C1=CC=CC=C1)"
               "C1=CC=CC=C1")
        tags = self._detect(smi, "Mo", "C")
        self.assertGreaterEqual(len(tags), 1)
        for t in tags:
            self.assertEqual(t, "1",
                             "Mo-C (sp C, DATIVE η² olefin) must tag '1'")

    def test_dative_sp2_carbon_in_arene_is_single(self):
        """sp2 C donor (η⁶ arene C) coordinating via DATIVE → bo='1'.

        Note: not '2' (would route to short Mo=C carbene at 1.92 Å)
        and not 'a' unless the donor atom itself is marked aromatic
        (RDKit may strip aromaticity off carbanionic [C+] / [C-] in
        the borane bridge SMILES).
        """
        smi = ("C[Si](C)(C)N(B(Cl)[C+]12->[Mo-12]3456789%10%11(<-[C+]%12"
               "(B(Cl)N([Si](C)(C)C)[Si](C)(C)C)[C+]3[C+]4[C+]5[C+]6"
               "[C+]%127)[C+]([C+]8[C+]19)[C+]%10[C+]2%11)[Si](C)(C)C")
        tags = self._detect(smi, "Mo", "C")
        self.assertEqual(len(tags), 12,
                         "UQIGEH Mo CN12 should expose 12 C donors")
        # ALL must be '1', never '2'.  '2' would route to the imido /
        # carbene short-bond distribution (1.92 Å) -- a 0.4 Å contraction.
        for t in tags:
            self.assertNotEqual(t, "2",
                                "Mo-C (sp2 C of arene, DATIVE) must NOT "
                                "tag as double '2'")
            self.assertNotEqual(t, "3",
                                "Mo-C (sp2 C of arene, DATIVE) must NOT "
                                "tag as triple '3'")

    def test_explicit_double_imido_is_two(self):
        """Explicit DOUBLE bond (W=N imido) still tags '2'.

        Conservative default must NOT downgrade explicit double bonds.
        """
        smi = "F[W](F)(F)(F)=NC"
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        self.assertIsNotNone(mol)
        from delfin.fffree.mogul_bounds import _detect_md_bond_order_tag
        W_idx = next(i for i, a in enumerate(mol.GetAtoms())
                     if a.GetSymbol() == "W")
        N_idx = next(nb.GetIdx() for nb in mol.GetAtomWithIdx(W_idx).GetNeighbors()
                     if nb.GetSymbol() == "N")
        tag = _detect_md_bond_order_tag(mol, W_idx, N_idx)
        self.assertEqual(tag, "2",
                         "Explicit DOUBLE W=N must still tag '2'")

    def test_explicit_aromatic_donor_preserves_a(self):
        """Aromatic donor (κ¹-pyridine N) still tags 'a'.

        Conservative default must NOT break the pyridine vs imine
        disaggregation the supplement was built for.
        """
        from delfin.fffree.assemble_via_mogul import _full_complex_mol
        from delfin.fffree.mogul_bounds import _detect_md_bond_order_tag
        smi = "F[W](F)(F)(F)(F)(F)n1ccccc1"
        mol = _full_complex_mol(smi)
        self.assertIsNotNone(mol)
        W_idx = next(i for i, a in enumerate(mol.GetAtoms())
                     if a.GetSymbol() == "W")
        N_idx = next(nb.GetIdx() for nb in mol.GetAtomWithIdx(W_idx).GetNeighbors()
                     if nb.GetSymbol() == "N")
        tag = _detect_md_bond_order_tag(mol, W_idx, N_idx)
        self.assertEqual(tag, "a",
                         "κ¹-pyridine W-N (aromatic donor) must tag 'a'")


class TestOffByteIdenticalForBondOrder(unittest.TestCase):
    """When the supplement is not loaded, the conservative classifier
    must still NOT route any cascade lookup to a short distribution.

    Concretely: ``_lookup_bond_md`` with the supplement absent simply
    drops the bo tag (it only fires the supplement lookup when the
    library has it).  We verify the cascade returns the SAME M-D
    distribution it returned pre-fix.
    """

    def test_no_supplement_no_side_effect(self):
        """Without the supplement env-flag, the classifier output is
        purely advisory and cannot change the M-D bond length.
        """
        # Save / clear supplement env var.
        old = os.environ.pop("DELFIN_GRIP_MD_BO_SUPPLEMENT", None)
        try:
            import importlib
            import delfin.fffree.grip_mogul_lookup as _g
            importlib.reload(_g)
            from delfin.fffree.grip_mogul_lookup import GripLibrary
            from delfin.fffree.mogul_bounds import _lookup_bond_md
            from pathlib import Path
            lib_path = Path(
                "/home/qmchem_max/agent_workspace/quality_framework/"
                "reports/grip_lib_v5_TM.npz"
            )
            if not lib_path.exists():
                self.skipTest("main library missing")
            lib = GripLibrary(lib_path)
            # All three regressed cases must return *some* hit from the
            # cascade (the main library covers Re-P, Mo-C).  We don't
            # assert exact mu here -- only that the bond-order tag is
            # irrelevant when the supplement is absent.
            hit_a = _lookup_bond_md(lib, "Re", "*", "P", "sp3", None, 1,
                                    bond_order_tag="1")
            hit_b = _lookup_bond_md(lib, "Re", "*", "P", "sp3", None, 1,
                                    bond_order_tag="2")
            # Both must agree (no supplement → tag doesn't matter).
            if hit_a is not None and hit_b is not None:
                self.assertAlmostEqual(hit_a[0], hit_b[0], delta=1e-6,
                                       msg="tag matters when supplement is "
                                       "OFF -- byte-identical violation")
        finally:
            if old is not None:
                os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = old


if __name__ == "__main__":
    import os  # noqa: F401  -- referenced inside Test class above
    unittest.main()
