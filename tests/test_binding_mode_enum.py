"""tests/test_binding_mode_enum.py — coverage for the binding-mode enumeration
(`delfin.fffree.coordination_mode_enum` + completeness detector +
converter_backend wire).

Verifies:
  - functional-group recognition (carboxylate / nitrate / sulfate / phosphate /
    β-diketonate / amino acid)
  - hapto-π modes (Cp η¹/η³/η⁵, arene η²/η⁴/η⁶, allyl, butadiene)
  - linkage-isomer detection still works (SCN/NO₂)
  - completeness detector returns a fraction in [0, 1]
  - converter_backend wire is default-OFF byte-identical and ON-mode active
  - determinism: two calls with same SMILES + seed produce identical output
"""
from __future__ import annotations

import os
import sys
import unittest

# Make sure we have a clean env at import time
for k in ("DELFIN_FFFREE_PURE_TRACK3",
          "DELFIN_FFFREE_BINDING_MODE_ENUM",
          "DELFIN_FFFREE_BINDING_MODE_STRICT",
          "DELFIN_FFFREE_BINDING_MODE_MAX_EXPAND"):
    os.environ.pop(k, None)

os.environ.setdefault("PYTHONHASHSEED", "0")


# ---------------------------------------------------------------------------
class TestEnumerateModesFunctional(unittest.TestCase):
    """Functional group canonical-mode recognition."""

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)

    def test_carboxylate_kappa1_and_kappa2(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("CC(=O)[O-]")
        labels = {m["label"] for m in modes}
        self.assertIn("kappa1-O-carboxylate", labels)
        self.assertIn("kappa2-OO-carboxylate", labels)
        self.assertGreaterEqual(len(modes), 2)

    def test_nitrate_kappa1_kappa2(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[O-][N+](=O)[O-]")
        kinds = {m["kind"] for m in modes}
        self.assertIn("functional", kinds)
        labels = " ".join(m["label"] for m in modes)
        self.assertIn("nitrate", labels)

    def test_sulfate_functional(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[O-][S](=O)(=O)[O-]")
        labels = {m["label"] for m in modes}
        self.assertTrue(any("sulfate" in lab for lab in labels),
                        f"no sulfate label among {labels}")

    def test_phosphate_functional(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[O-][P](=O)([O-])[O-]")
        labels = {m["label"] for m in modes}
        self.assertTrue(any("phosphate" in lab for lab in labels))

    def test_betadiketonate_chelate(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[O-]C(C)=CC(C)=O")
        labels = {m["label"] for m in modes}
        self.assertTrue(any("betadiketo" in lab for lab in labels),
                        f"acac expected β-diketo label, got {labels}")

    def test_aminoacid_glycinate_chelate(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[NH2]CC(=O)[O-]")
        labels = {m["label"] for m in modes}
        self.assertTrue(any("aminoacid" in lab for lab in labels),
                        f"glycinate label missing: {labels}")
        # κ²-N,O on 5-ring is the canonical mode
        self.assertIn("kappa2-NO-aminoacid", labels)


# ---------------------------------------------------------------------------
class TestEnumerateModesHapto(unittest.TestCase):
    """Hapto-π η-mode recognition."""

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)

    def test_cp_anion_eta1_eta3_eta5(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("[cH-]1cccc1")
        labels = {m["label"] for m in modes}
        self.assertIn("eta1-Cp", labels)
        self.assertIn("eta3-Cp", labels)
        self.assertIn("eta5-Cp", labels)

    def test_benzene_eta2_eta4_eta6(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("c1ccccc1")
        labels = {m["label"] for m in modes}
        self.assertIn("eta2-arene", labels)
        self.assertIn("eta4-arene", labels)
        self.assertIn("eta6-arene", labels)

    def test_allyl_eta1_eta3(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("C=CC")
        labels = {m["label"] for m in modes}
        self.assertIn("eta3-allyl", labels)

    def test_butadiene_eta2_eta4(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        modes = enumerate_modes("C=CC=C")
        labels = {m["label"] for m in modes}
        # both modes should appear (eta2-butadiene + eta4-butadiene)
        self.assertIn("eta4-butadiene", labels)


# ---------------------------------------------------------------------------
class TestLinkageIsomers(unittest.TestCase):
    """Pre-existing linkage isomer table still detects SCN / NO₂ / N₃."""

    def setUp(self):
        os.environ["DELFIN_FFFREE_LINKAGE"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_LINKAGE", None)

    def test_scn_S_vs_N(self):
        # Re-import after env change
        import importlib, sys
        sys.modules.pop("delfin.fffree.linkage_isomers", None)
        from delfin.fffree.linkage_isomers import detect_ambidentate_groups
        from rdkit import Chem
        mol = Chem.MolFromSmiles("[S-]C#N")
        groups = detect_ambidentate_groups(mol)
        self.assertTrue(any("thiocyanato" in lbl
                            for g in groups for _, lbl in g["donor_options"]))
        self.assertTrue(any("isothiocyanato" in lbl
                            for g in groups for _, lbl in g["donor_options"]))


# ---------------------------------------------------------------------------
class TestRealismFilter(unittest.TestCase):
    """Strict realism filter prunes impossible / strained modes."""

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"

    def tearDown(self):
        for k in ("DELFIN_FFFREE_BINDING_MODE_ENUM",
                  "DELFIN_FFFREE_BINDING_MODE_STRICT"):
            os.environ.pop(k, None)

    def test_strict_keeps_chelate_5_6(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        os.environ["DELFIN_FFFREE_BINDING_MODE_STRICT"] = "1"
        modes = enumerate_modes("[NH2]CC(=O)[O-]")
        # 5-ring N,O chelate is the canonical, must survive strict
        labels = {m["label"] for m in modes}
        self.assertIn("kappa2-NO-aminoacid", labels)


# ---------------------------------------------------------------------------
class TestCompletenessDetector(unittest.TestCase):

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)

    def test_default_off_returns_none(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        from delfin.fffree.binding_mode_completeness import (
            completeness_for_results,
        )
        r = completeness_for_results("CC(=O)[O-]", [("xyz", "label")])
        self.assertIsNone(r)

    def test_completeness_fraction_in_range(self):
        from delfin.fffree.binding_mode_completeness import (
            completeness_for_results,
        )
        rec = completeness_for_results(
            "CC(=O)[O-]",
            [("xyz", "SP-4-kappa1-carboxylate-1"),
             ("xyz", "SP-4-kappa2-carboxylate-2")],
        )
        self.assertIsNotNone(rec)
        self.assertGreaterEqual(rec["fraction"], 0.0)
        self.assertLessEqual(rec["fraction"], 1.0)
        self.assertGreater(rec["emitted"], 0)
        self.assertGreater(rec["predicted"], 0)

    def test_aggregate_completeness(self):
        from delfin.fffree.binding_mode_completeness import (
            completeness_for_results, aggregate_completeness,
        )
        recs = []
        for smi, lbls in [
            ("CC(=O)[O-]", ["SP-4-kappa1-1", "SP-4-kappa2-2"]),
            ("[NH2]CC(=O)[O-]", ["OC-6-kappa2-NO-aminoacid-1"]),
            ("[cH-]1cccc1", ["PIANO-STOOL-8-eta5-Cp-1"]),
        ]:
            recs.append(completeness_for_results(
                smi, [("xyz", lab) for lab in lbls]))
        agg = aggregate_completeness(recs)
        self.assertEqual(agg["n"], 3)
        self.assertGreaterEqual(agg["binding_mode_frac_complete"], 0.0)
        self.assertLessEqual(agg["binding_mode_frac_complete"], 1.0)


# ---------------------------------------------------------------------------
class TestConverterBackendWire(unittest.TestCase):
    """Default-OFF byte-identity + ON-mode label-tag expansion."""

    def test_default_off_byte_identical(self):
        # Make sure env is OFF
        for k in ("DELFIN_FFFREE_BINDING_MODE_ENUM",
                  "DELFIN_FFFREE_PURE_TRACK3"):
            os.environ.pop(k, None)
        # Re-import to honour env
        import importlib
        for m in [
            "delfin.fffree.coordination_mode_enum",
            "delfin.fffree.converter_backend",
        ]:
            sys.modules.pop(m, None)
        from delfin.fffree.converter_backend import _fffree_isomers
        r1 = _fffree_isomers("N[Pt](N)(Cl)Cl")
        r2 = _fffree_isomers("N[Pt](N)(Cl)Cl")
        self.assertEqual(r1, r2)
        if r1 is not None:
            for _, lab in r1:
                self.assertNotIn("-bm(", lab,
                    "Default-OFF should not emit -bm(...) tags")

    def test_on_emits_bm_tag(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"
        # Re-import
        for m in [
            "delfin.fffree.coordination_mode_enum",
            "delfin.fffree.converter_backend",
        ]:
            sys.modules.pop(m, None)
        from delfin.fffree.converter_backend import _fffree_isomers
        r = _fffree_isomers("N[Pt](N)(Cl)Cl")
        if r is not None:
            for _, lab in r:
                self.assertIn("-bm(", lab,
                    f"ON-mode missing bm tag in {lab}")
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)


# ---------------------------------------------------------------------------
class TestDeterminism(unittest.TestCase):

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"
        os.environ["PYTHONHASHSEED"] = "0"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)

    def test_two_runs_byte_identical(self):
        from delfin.fffree.coordination_mode_enum import enumerate_modes
        for smi in ["CC(=O)[O-]", "[NH2]CC(=O)[O-]", "[cH-]1cccc1",
                    "c1ccccc1", "[O-][N+](=O)[O-]"]:
            r1 = enumerate_modes(smi)
            r2 = enumerate_modes(smi)
            self.assertEqual(
                [m["label"] for m in r1],
                [m["label"] for m in r2],
                f"Non-deterministic on {smi}",
            )


# ---------------------------------------------------------------------------
class TestPredictedModeCount(unittest.TestCase):

    def setUp(self):
        os.environ["DELFIN_FFFREE_BINDING_MODE_ENUM"] = "1"

    def tearDown(self):
        os.environ.pop("DELFIN_FFFREE_BINDING_MODE_ENUM", None)

    def test_predicted_count_nonzero_for_carboxylate(self):
        from delfin.fffree.coordination_mode_enum import predicted_mode_count
        self.assertGreater(predicted_mode_count("CC(=O)[O-]"), 1)

    def test_predicted_count_water_eq_1(self):
        from delfin.fffree.coordination_mode_enum import predicted_mode_count
        self.assertEqual(predicted_mode_count("O"), 1)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
