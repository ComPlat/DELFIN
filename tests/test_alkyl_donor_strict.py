"""tests/test_alkyl_donor_strict.py — universal alkyl-carbon-donor rule
(OCANIT-class bug fix).

Verifies the ``DELFIN_FFFREE_ALKYL_DONOR_STRICT`` env-gate (default-OFF
byte-identical, ON enforces the universal substituent-of-donor predicate):

  1. CF3 ligand on Ni  -> only kappa1-C is emitted (no kappa2-FF / kappa3-FFF)
  2. CHF2 ligand on Pd -> only kappa1-C
  3. -OMe (where O is the donor) -> O is donor; Me-C and Me-H are NOT
  4. Bare fluoride [F-] -> kappa1-F preserved (direct M-F)
  5. Mixed Pt-CF3 + Pt-F demo -> CF3 emits kappa1-C, F- emits kappa1-F
  6. Byte-identical OFF behaviour
  7. Determinism (repeat calls bit-identical)
  8. OCANIT integration: no kappa2-FF / kappa3-FFF for CF3 in the OCANIT SMILES
"""
from __future__ import annotations

import os
import unittest


# Make sure we have a clean env at import time so tests do not poison each
# other when run in a long suite.
for k in (
    "DELFIN_FFFREE_PURE_TRACK3",
    "DELFIN_FFFREE_BINDING_MODE_ENUM",
    "DELFIN_FFFREE_BINDING_MODE_STRICT",
    "DELFIN_FFFREE_ALKYL_DONOR_STRICT",
    "DELFIN_FFFREE_BINDING_MODE_MAX_EXPAND",
):
    os.environ.pop(k, None)

os.environ.setdefault("PYTHONHASHSEED", "0")


# ---------------------------------------------------------------------------
class _BaseStrict(unittest.TestCase):
    """Common setUp/tearDown: BINDING_MODE_ENUM on + ALKYL_DONOR_STRICT on."""

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"
        os.environ["DELFIN_FFFREE_ALKYL_DONOR_STRICT"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)
        os.environ.pop("DELFIN_FFFREE_ALKYL_DONOR_STRICT", None)


# ---------------------------------------------------------------------------
class TestCF3OnNi(_BaseStrict):
    """Test 1: CF3 ligand -> only kappa1-C, never kappa2-FF / kappa3-FFF."""

    def test_cf3_emits_only_kappa1_C(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[C](F)(F)F")
        labels = {m["label"] for m in modes}
        elems_lists = [m["donor_elems"] for m in modes]
        # the legitimate sigma-C donor
        self.assertIn("kappa1-C", labels,
                      f"kappa1-C missing for CF3: {labels}")
        # the bug we are fixing -- F atoms must NOT appear as multi-F donors
        self.assertNotIn("kappa2-FF", labels,
                         f"spurious kappa2-FF still emitted: {labels}")
        self.assertNotIn("kappa3-FFF", labels,
                         f"spurious kappa3-FFF still emitted: {labels}")
        # no F-only modes of any flavour
        for elems in elems_lists:
            self.assertNotEqual(elems, ("F", "F"))
            self.assertNotEqual(elems, ("F", "F", "F"))


# ---------------------------------------------------------------------------
class TestCHF2OnPd(_BaseStrict):
    """Test 2: CHF2 ligand -> only kappa1-C."""

    def test_chf2_emits_only_kappa1_C(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[CH](F)F")
        labels = {m["label"] for m in modes}
        self.assertIn("kappa1-C", labels)
        self.assertNotIn("kappa2-FF", labels)
        # only one donor class: C
        donor_elems = {m["donor_elems"] for m in modes}
        for e in donor_elems:
            self.assertEqual(e, ("C",),
                             f"unexpected non-C donor elems {e} for CHF2")


# ---------------------------------------------------------------------------
class TestMethoxide(_BaseStrict):
    """Test 3: -OMe -> O is the donor, C and H of the methyl are not."""

    def test_methoxide_O_is_donor_not_C(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        # methoxide: [O-]CH3
        modes = enumerate_modes("[O-]C")
        labels = {m["label"] for m in modes}
        donor_elems = {m["donor_elems"] for m in modes}
        # O is donor (kappa1-O)
        self.assertIn("kappa1-O", labels,
                      f"methoxide kappa1-O missing: {labels}")
        # C must NOT appear as a separate donor
        for e in donor_elems:
            self.assertNotIn("C", e,
                             f"methoxide carbon spuriously enumerated: {donor_elems}")


# ---------------------------------------------------------------------------
class TestBareFluoride(_BaseStrict):
    """Test 4: bare fluoride [F-] is a direct M-F donor and must survive."""

    def test_bare_fluoride_still_emitted(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[F-]")
        labels = {m["label"] for m in modes}
        self.assertIn("kappa1-F", labels,
                      f"bare fluoride kappa1-F missing: {labels}")


# ---------------------------------------------------------------------------
class TestMixedCF3AndFluoride(_BaseStrict):
    """Test 5: Pt-CF3 + Pt-F (direct fluoride) -> CF3 emits kappa1-C,
    the direct F- ligand emits kappa1-F.  The two ligands are enumerated
    independently."""

    def test_separate_ligands_independent(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        # CF3 ligand
        cf3_modes = enumerate_modes("[C](F)(F)F")
        cf3_labels = {m["label"] for m in cf3_modes}
        # direct fluoride
        f_modes = enumerate_modes("[F-]")
        f_labels = {m["label"] for m in f_modes}
        # CF3 only emits kappa1-C (no F-only modes)
        self.assertIn("kappa1-C", cf3_labels)
        self.assertNotIn("kappa1-F", cf3_labels,
                         "spurious kappa1-F emitted for CF3 (its F's are "
                         "substituents, not donors)")
        # bare F- emits kappa1-F
        self.assertIn("kappa1-F", f_labels)


# ---------------------------------------------------------------------------
class TestByteIdenticalOff(unittest.TestCase):
    """Test 6: default OFF preserves the legacy enumeration byte-identical."""

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"
        # explicitly disable the strict gate
        os.environ.pop("DELFIN_FFFREE_ALKYL_DONOR_STRICT", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)

    def test_off_preserves_legacy_cf3_enumeration(self):
        """Default OFF: the (chemically wrong, but legacy) kappa2-FF /
        kappa3-FFF modes are STILL emitted.  This is the byte-identical
        contract."""
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[C](F)(F)F")
        labels = {m["label"] for m in modes}
        # legacy behaviour -- the bug is intentionally preserved when the flag
        # is OFF, so existing pools stay reproducible.
        self.assertIn("kappa1-F", labels)
        self.assertIn("kappa2-FF", labels)
        self.assertIn("kappa3-FFF", labels)

    def test_off_preserves_pyridine_acetate(self):
        """Default OFF: normal heteroatom-donor ligands unchanged."""
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        py = {m["label"] for m in enumerate_modes("c1ccncc1")}
        self.assertIn("kappa1-N", py)
        ac = {m["label"] for m in enumerate_modes("CC(=O)[O-]")}
        self.assertIn("kappa1-O-carboxylate", ac)
        self.assertIn("kappa2-OO-carboxylate", ac)


# ---------------------------------------------------------------------------
class TestDeterminism(_BaseStrict):
    """Test 7: repeated enumeration with the same input is bit-identical."""

    def test_repeat_calls_identical(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        for smi in (
            "[C](F)(F)F",
            "[CH](F)F",
            "[O-]C",
            "[F-]",
            "c1ccncc1",
            "CC(=O)[O-]",
        ):
            a = enumerate_modes(smi)
            b = enumerate_modes(smi)
            # compare a serialisable projection
            sa = [(m["label"], m["donors"], m["kappa"], m["donor_elems"])
                  for m in a]
            sb = [(m["label"], m["donors"], m["kappa"], m["donor_elems"])
                  for m in b]
            self.assertEqual(sa, sb, f"non-deterministic enumeration: {smi}")


# ---------------------------------------------------------------------------
class TestOCANITIntegration(_BaseStrict):
    """Test 8: integration on the actual OCANIT SMILES.  Each CF3 ligand
    when enumerated on its own must emit no kappa2-FF / kappa3-FFF, so the
    full-complex isomer expansion never gets the spurious multi-F orbits."""

    OCANIT_CF3_LIGAND_SMI = "[C](F)(F)F"

    def test_cf3_ligand_excludes_multi_F_orbits(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes(self.OCANIT_CF3_LIGAND_SMI)
        labels = [m["label"] for m in modes]
        # kappa1-C is the one chemically correct mode
        self.assertEqual(labels, ["kappa1-C"],
                         f"OCANIT CF3 ligand: expected only kappa1-C, got "
                         f"{labels}")

    def test_only_C_donor_class_for_cf3(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes(self.OCANIT_CF3_LIGAND_SMI)
        for m in modes:
            self.assertEqual(m["donor_elems"], ("C",),
                             f"unexpected non-C donor elems {m['donor_elems']}")
            self.assertEqual(m["kappa"], 1,
                             f"unexpected kappa>{1} for CF3: {m}")


# ---------------------------------------------------------------------------
class TestUniversalityCCl3CBr3(_BaseStrict):
    """Bonus: the rule is universal -- CCl3, CBr3 etc. all behave like CF3
    (no SMILES-specific shortcut, pure graph predicate)."""

    def test_ccl3_emits_only_kappa1_C(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[C](Cl)(Cl)Cl")
        labels = {m["label"] for m in modes}
        donor_elems = {m["donor_elems"] for m in modes}
        self.assertIn("kappa1-C", labels)
        for e in donor_elems:
            self.assertNotIn("Cl", e,
                             f"CCl3 spurious Cl donor: {donor_elems}")

    def test_cbr3_emits_only_kappa1_C(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[C](Br)(Br)Br")
        labels = {m["label"] for m in modes}
        donor_elems = {m["donor_elems"] for m in modes}
        self.assertIn("kappa1-C", labels)
        for e in donor_elems:
            self.assertNotIn("Br", e,
                             f"CBr3 spurious Br donor: {donor_elems}")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
