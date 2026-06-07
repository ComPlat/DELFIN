"""Tests for the f-block-aware polyhedron picker (2026-06-07).

Validates the contract described in
``delfin.fffree.f_block_picker``:

  * Y / Cd / Mo / Sc at CN 7-12 receive an expanded candidate list (SAP-8,
    DD-8, TTP-9, CSAP-9, BICAP-10, ...) when DELFIN_FFFREE_F_BLOCK_PICK=1.
  * Lanthanides (Z 57-71) and actinides (Z 89-103) at CN 7-12 are routed
    through the f-block picker.
  * Env-flag OFF => byte-identical with legacy first-candidate rule.
  * Picker is deterministic and lex-stable across runs.
  * f-block metal classification (Pyykkö covalent radius table).

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import os
import unittest

import numpy as np


class _EnvSnapshotMixin:
    """Save/restore env state around each test."""

    _SAVED_KEYS = (
        "DELFIN_FFFREE_F_BLOCK_PICK",
        "DELFIN_FFFREE_FBLOCK_CN8_12",
        "DELFIN_FFFREE_MOGUL_PRIMARY",
        "DELFIN_FFFREE_PURE_TRACK3",
    )

    def setUp(self):
        super().setUp()
        self._saved = {k: os.environ.get(k) for k in self._SAVED_KEYS}

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        super().tearDown()


class TestEnvFlag(_EnvSnapshotMixin, unittest.TestCase):
    """Env-flag behaviour: default-OFF, MOGUL_PRIMARY auto-ON."""

    def test_off_when_unset_and_no_mogul_primary(self):
        from delfin.fffree.f_block_picker import f_block_picker_enabled
        os.environ.pop("DELFIN_FFFREE_F_BLOCK_PICK", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        self.assertFalse(f_block_picker_enabled())

    def test_auto_on_under_mogul_primary(self):
        from delfin.fffree.f_block_picker import f_block_picker_enabled
        os.environ.pop("DELFIN_FFFREE_F_BLOCK_PICK", None)
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        self.assertTrue(f_block_picker_enabled())

    def test_explicit_off_overrides_mogul_primary(self):
        from delfin.fffree.f_block_picker import f_block_picker_enabled
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "0"
        self.assertFalse(f_block_picker_enabled())

    def test_explicit_on(self):
        from delfin.fffree.f_block_picker import f_block_picker_enabled
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        self.assertTrue(f_block_picker_enabled())


class TestPyykkoTable(unittest.TestCase):
    """The Pyykkö covalent-radius table covers all f-block + key metals."""

    def test_all_lanthanides_present(self):
        from delfin.fffree.f_block_picker import PYYKKO_COV_RADIUS
        for s in ("La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
                  "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"):
            self.assertIn(s, PYYKKO_COV_RADIUS, f"missing Ln: {s}")
            self.assertGreater(PYYKKO_COV_RADIUS[s], 1.5)
            self.assertLess(PYYKKO_COV_RADIUS[s], 1.9)

    def test_actinides_present(self):
        from delfin.fffree.f_block_picker import PYYKKO_COV_RADIUS
        for s in ("Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm"):
            self.assertIn(s, PYYKKO_COV_RADIUS, f"missing An: {s}")

    def test_mission_test_metals_present(self):
        from delfin.fffree.f_block_picker import PYYKKO_COV_RADIUS
        for s in ("Y", "Sc", "Cd", "Mo"):
            self.assertIn(s, PYYKKO_COV_RADIUS, f"missing target metal: {s}")

    def test_pyykko_delta_la_o_close_to_2_43(self):
        # r(La) ≈ 1.80, r(O) ≈ 0.63 → ideal La-O = 2.43 Å.
        from delfin.fffree.f_block_picker import pyykko_delta
        delta = pyykko_delta("La", ["O"] * 8, [2.43] * 8)
        self.assertIsNotNone(delta)
        self.assertLess(abs(delta), 0.005)

    def test_pyykko_delta_unknown_metal_returns_none(self):
        from delfin.fffree.f_block_picker import pyykko_delta
        self.assertIsNone(pyykko_delta("Xx", ["O"], [2.5]))


class TestHighCnTarget(unittest.TestCase):
    """Atomic-number gate for the high-CN picker scope."""

    def test_y_cn8_in_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertTrue(is_high_cn_non_fblock_target("Y", 8))

    def test_sc_cn7_in_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertTrue(is_high_cn_non_fblock_target("Sc", 7))

    def test_cd_cn8_in_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertTrue(is_high_cn_non_fblock_target("Cd", 8))

    def test_mo_cn8_in_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertTrue(is_high_cn_non_fblock_target("Mo", 8))

    def test_la_cn8_in_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertTrue(is_high_cn_non_fblock_target("La", 8))

    def test_cu_cn4_out_of_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertFalse(is_high_cn_non_fblock_target("Cu", 4))

    def test_fe_cn6_out_of_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertFalse(is_high_cn_non_fblock_target("Fe", 6))

    def test_unknown_metal_out_of_scope(self):
        from delfin.fffree.f_block_picker import is_high_cn_non_fblock_target
        self.assertFalse(is_high_cn_non_fblock_target("Xx", 8))


class TestCandidateList(_EnvSnapshotMixin, unittest.TestCase):
    """f_block_candidates_for_cn returns the canonical merged list."""

    def test_cn8_includes_sap8_and_dd8(self):
        from delfin.fffree.f_block_picker import f_block_candidates_for_cn
        cands = f_block_candidates_for_cn("La", 8)
        # SAP-8 must be present
        self.assertTrue(any("SAP-8" in c for c in cands),
                        f"SAP-8 missing from {cands}")
        # DD-8 must be present
        self.assertTrue(any("DD-8" in c for c in cands),
                        f"DD-8 missing from {cands}")

    def test_cn9_includes_ttp9_and_csap9(self):
        from delfin.fffree.f_block_picker import f_block_candidates_for_cn
        cands = f_block_candidates_for_cn("La", 9)
        self.assertTrue(any("TTP-9" in c for c in cands))
        self.assertTrue(any("CSAP-9" in c for c in cands))

    def test_cn12_includes_icosahedron(self):
        from delfin.fffree.f_block_picker import f_block_candidates_for_cn
        cands = f_block_candidates_for_cn("Ce", 12)
        self.assertTrue(any("IH-12" in c for c in cands))


class TestPicker(_EnvSnapshotMixin, unittest.TestCase):
    """End-to-end picker behaviour for the three smoke-100 test cases."""

    def test_y_cn8_sap8_picked(self):
        """Y CN=8 with 8 oxygen donors → SAP-8 should win (D4d high-sym)."""
        from delfin.fffree.f_block_picker import f_block_polyhedron_picker
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        picked = f_block_polyhedron_picker(
            metal_sym="Y",
            cn=8,
            donors=["O"] * 8,
            force=True,
        )
        # SAP-8 (D4d) is the highest-symmetry CN8 polyhedron; for an
        # all-oxygen donor set with ionic Y-O bonding the picker should
        # land on SAP-8 (either the canonical name or the SQAP-8 alias).
        self.assertIsNotNone(picked)
        self.assertTrue(
            "SAP-8" in picked or "SQAP-8" in picked,
            f"expected SAP-8 family for Y(O)8, got {picked!r}",
        )

    def test_lanthanide_high_cn_route_to_fblock(self):
        """La / Ce / U at CN 8-12 should select an f-block polyhedron."""
        from delfin.fffree.f_block_picker import f_block_polyhedron_picker
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        for metal, cn, expect_family in [
            ("La", 8, ("SAP-8", "DD-8", "SQAP-8")),
            ("Ce", 9, ("TTP-9", "CSAP-9")),
            ("Nd", 10, ("BICAP-10",)),
            ("U", 12, ("IH-12",)),
        ]:
            picked = f_block_polyhedron_picker(
                metal_sym=metal,
                cn=cn,
                donors=["O"] * cn,
                force=True,
            )
            self.assertIsNotNone(picked, f"None returned for {metal} CN={cn}")
            self.assertTrue(
                any(fam in picked for fam in expect_family),
                f"{metal} CN={cn}: expected one of {expect_family}, got {picked!r}",
            )

    def test_cd_cn8_picks_high_symmetry(self):
        """Cd CN=8 (e.g. Cd-uracil mixed-donor) → SAP-8 / DD-8 family."""
        from delfin.fffree.f_block_picker import f_block_polyhedron_picker
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        picked = f_block_polyhedron_picker(
            metal_sym="Cd",
            cn=8,
            donors=["O", "O", "O", "O", "N", "N", "N", "N"],
            force=True,
        )
        self.assertIsNotNone(picked)
        self.assertTrue(
            "SAP-8" in picked or "SQAP-8" in picked or "DD-8" in picked,
            f"expected SAP-8/DD-8 family for Cd(O4N4), got {picked!r}",
        )

    def test_sc_cn7_picks_pb7(self):
        """Sc CN=7 → pentagonal bipyramid (only CN7 polyhedron in the table)."""
        from delfin.fffree.f_block_picker import f_block_polyhedron_picker
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        picked = f_block_polyhedron_picker(
            metal_sym="Sc",
            cn=7,
            donors=["O"] * 7,
            force=True,
        )
        self.assertEqual(picked, "PB-7 pentagonal bipyramid")


class TestOffByteIdentical(_EnvSnapshotMixin, unittest.TestCase):
    """Env-OFF behaviour: picker returns first candidate (= legacy)."""

    def test_off_byte_identical(self):
        """When DELFIN_FFFREE_F_BLOCK_PICK is unset AND MOGUL_PRIMARY is unset,
        the picker MUST return the first candidate from
        ``polyhedra.geometries_for_cn`` (the legacy first-rule), so the
        whole-stack output is byte-identical to HEAD.
        """
        from delfin.fffree import polyhedra as _polyhedra
        from delfin.fffree.f_block_picker import f_block_polyhedron_picker
        os.environ.pop("DELFIN_FFFREE_F_BLOCK_PICK", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        os.environ.pop("DELFIN_FFFREE_FBLOCK_CN8_12", None)
        # Y CN=8 with no f-block gate: legacy is SQAP-8
        legacy_cands = _polyhedra.geometries_for_cn(8, "Y")
        self.assertEqual(legacy_cands, ["SQAP-8 square antiprism"])
        picked = f_block_polyhedron_picker(
            metal_sym="Y", cn=8, donors=["O"] * 8,
        )
        self.assertEqual(picked, "SQAP-8 square antiprism")

    def test_polyhedra_off_byte_identical_for_y_cn8(self):
        """When the F_BLOCK_PICK gate is off, ``geometries_for_cn(8, 'Y')``
        returns the legacy single-candidate list — the picker extension
        is fully invisible."""
        from delfin.fffree import polyhedra as _polyhedra
        os.environ.pop("DELFIN_FFFREE_F_BLOCK_PICK", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        os.environ.pop("DELFIN_FFFREE_FBLOCK_CN8_12", None)
        cands = _polyhedra.geometries_for_cn(8, "Y")
        self.assertEqual(cands, ["SQAP-8 square antiprism"])

    def test_polyhedra_on_extends_y_cn8(self):
        """When F_BLOCK_PICK=1, ``geometries_for_cn(8, 'Y')`` extends to
        SAP-8 + DD-8 even for non-f-block Y."""
        from delfin.fffree import polyhedra as _polyhedra
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        os.environ.pop("DELFIN_FFFREE_FBLOCK_CN8_12", None)
        cands = _polyhedra.geometries_for_cn(8, "Y")
        # Must include SQAP-8 (legacy) plus SAP-8 + DD-8 (extension)
        self.assertIn("SQAP-8 square antiprism", cands)
        self.assertIn("SAP-8 square antiprism", cands)
        self.assertIn("DD-8 dodecahedron", cands)


class TestDeterminism(_EnvSnapshotMixin, unittest.TestCase):
    """The picker must be deterministic across repeated calls."""

    def test_repeated_calls_byte_identical(self):
        from delfin.fffree.f_block_picker import f_block_polyhedron_picker
        os.environ["DELFIN_FFFREE_F_BLOCK_PICK"] = "1"
        outs = set()
        for _ in range(20):
            p = f_block_polyhedron_picker(
                metal_sym="Eu", cn=9, donors=["O"] * 9, force=True,
            )
            outs.add(p)
        self.assertEqual(len(outs), 1,
                         f"non-deterministic picker output: {outs}")


class TestPyykkoDelta(_EnvSnapshotMixin, unittest.TestCase):
    """Pyykkö-Δ metric matches the smoke 100 report expectation."""

    def test_pyykko_delta_y_cn8_o(self):
        """Y CN=8 with M-D = (2.30, 2.32, 2.34) etc. gives Δ ~ 0.10-0.25.

        The smoke 100 report quoted ``Y CN=8: Δ-Pyykkö 0.217`` —
        we don't reproduce that here (we don't have the per-file MD
        list), we just verify that the metric is in the expected order
        of magnitude for representative Y-O distances.
        """
        from delfin.fffree.f_block_picker import pyykko_delta
        md = [2.30, 2.31, 2.33, 2.34, 2.36, 2.38, 2.40, 2.45]
        delta = pyykko_delta("Y", ["O"] * 8, md)
        self.assertIsNotNone(delta)
        # ideal r(Y) + r(O) = 1.63 + 0.63 = 2.26
        # mean(md) = 2.359 → mean |d - 2.26| ~ 0.099
        self.assertGreater(delta, 0.05)
        self.assertLess(delta, 0.30)


if __name__ == "__main__":
    unittest.main()
