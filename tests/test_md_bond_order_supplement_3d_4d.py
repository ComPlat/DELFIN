"""Tests for the 3d/4d M-D bond-order supplement library.

The 3d/4d supplement extends the 5d-mid (W/Mo/Tc/Re) library to cover the
remaining transition metals with multi-modal M-D distance distributions:
3d (Ti/V/Cr/Mn/Fe/Co/Ni/Cu), 4d-early (Zr/Nb), 4d-late (Pd/Rh), 5d-early
(Hf/Ta) and 5d-late (Pt/Ir/Au).

User-Befund GIBSOA: a Cr complex shows M-O distances 1.507, 1.519, 1.573,
1.864, 2.053, 2.131 Å -- the multi-modal collapse of Cr=O oxo (~1.65 Å,
BO=2) pooled with Cr-O alkoxide (~1.95 Å, BO=1).  Disaggregation by
bond-order tag separates the populations.

The supplement is additive (env-gated) and byte-identical when the env-flag
``DELFIN_GRIP_MD_BO_SUPPLEMENT`` is unset.  When the env-var is set with an
``os.pathsep``-separated list of supplement paths, ALL listed supplements
load and the lookup picks the first hit covering the requested metal.

Author: hmaximilian
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path

import numpy as np

SUP_3D_4D_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/"
    "grip_lib_md_bo_supplement_3d_4d.npz"
)
SUP_5D_MID_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/"
    "grip_lib_md_bo_supplement_5d_mid.npz"
)
MAIN_LIB_PATH = Path(
    "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_TM.npz"
)


@unittest.skipUnless(SUP_3D_4D_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "3d/4d supplement / main library not present")
class Test3D4DSupplementLoads(unittest.TestCase):
    """The 3d/4d supplement loads and exposes the expected target metals."""

    def setUp(self):
        self._old_env = os.environ.get("DELFIN_GRIP_MD_BO_SUPPLEMENT")
        os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = str(SUP_3D_4D_PATH)
        import importlib
        import delfin.fffree.grip_mogul_lookup as _m
        importlib.reload(_m)
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        self.lib = GripLibrary(MAIN_LIB_PATH)

    def tearDown(self):
        if self._old_env is None:
            os.environ.pop("DELFIN_GRIP_MD_BO_SUPPLEMENT", None)
        else:
            os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = self._old_env

    def test_3d_4d_supplement_loads(self):
        self.assertGreater(self.lib.n_md_bo, 100,
                           "3d/4d supplement should hold hundreds of M-D keys")
        # 17 target metals: 3d-early, 3d-late, 4d-early, 4d-late, 5d-early, 5d-late
        expected_metals = {
            "Ti", "V", "Cr", "Mn",        # 3d-early
            "Fe", "Co", "Ni", "Cu",       # 3d-late
            "Zr", "Nb",                   # 4d-early
            "Pd", "Rh",                   # 4d-late
            "Hf", "Ta",                   # 5d-early
            "Pt", "Ir", "Au",             # 5d-late
        }
        self.assertEqual(self.lib._md_bo_metals, frozenset(expected_metals),
                         "3d/4d supplement must claim all 17 target metals")


@unittest.skipUnless(SUP_3D_4D_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "3d/4d supplement / main library not present")
class TestCrOOxoVsAlkoxideDisaggregated(unittest.TestCase):
    """The Cr-O bond is bi-modal (oxo vs alkoxide); the supplement must
    separate the two populations.

    GIBSOA evidence: Cr-O 1.507, 1.519, 1.573, 1.864, 2.053, 2.131 Å.
    Pre-supplement (hyb-only key): single μ around 1.8 Å -- wrong for both
    populations.  Post-supplement (bo-tag key):
      * Cr=O(sp2, '@bo2') = oxo,         μ ≈ 1.62-1.75 Å
      * Cr-O(sp3, '@bo1') = alkoxide,    μ ≈ 1.85-2.05 Å
    """

    def setUp(self):
        self._old_env = os.environ.get("DELFIN_GRIP_MD_BO_SUPPLEMENT")
        os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = str(SUP_3D_4D_PATH)
        import importlib
        import delfin.fffree.grip_mogul_lookup as _m
        importlib.reload(_m)
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        self.lib = GripLibrary(MAIN_LIB_PATH)

    def tearDown(self):
        if self._old_env is None:
            os.environ.pop("DELFIN_GRIP_MD_BO_SUPPLEMENT", None)
        else:
            os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = self._old_env

    def test_cr_o_oxo_vs_alkoxide_disaggregated(self):
        alkoxide = self.lib._md_bo_lookup("Cr", "O", "sp3", "1", min_n=10)
        oxo = self.lib._md_bo_lookup("Cr", "O", "sp2", "2", min_n=10)
        self.assertIsNotNone(alkoxide, "Cr-O(sp3, bo=1) should resolve")
        self.assertIsNotNone(oxo, "Cr=O(sp2, bo=2) should resolve")
        # Chemistry: oxo (M=O double bond) is strictly shorter than alkoxide.
        self.assertGreater(alkoxide[0] - oxo[0], 0.10,
                           "Cr-O alkoxide should be ≥0.10 Å longer than Cr=O oxo")
        # Realistic empirical means: oxo ≈ 1.55-1.75 Å; alkoxide ≈ 1.85-2.10 Å.
        self.assertGreater(oxo[0], 1.50)
        self.assertLess(oxo[0], 1.80)
        self.assertGreater(alkoxide[0], 1.75)
        self.assertLess(alkoxide[0], 2.20)

    def test_fe_o_disaggregated(self):
        """Fe-O double-bonded ferryl (~1.62 Å) vs alkoxide / hydroxide (~1.95 Å)."""
        alkoxide = self.lib._md_bo_lookup("Fe", "O", "sp3", "1", min_n=10)
        ferryl = self.lib._md_bo_lookup("Fe", "O", "sp2", "2", min_n=10)
        self.assertIsNotNone(alkoxide, "Fe-O(sp3, bo=1) should resolve")
        if ferryl is not None:
            # When the ferryl population exists, the disaggregation must hold.
            self.assertGreater(alkoxide[0] - ferryl[0], 0.10,
                               "Fe-O single bond should be longer than Fe=O double")

    def test_v_o_disaggregated(self):
        """V=O vanadyl (~1.59 Å) vs V-O alkoxide/carboxylate (~1.95 Å)."""
        alkoxide = self.lib._md_bo_lookup("V", "O", "sp3", "1", min_n=10)
        vanadyl = self.lib._md_bo_lookup("V", "O", "sp2", "2", min_n=10)
        self.assertIsNotNone(alkoxide, "V-O(sp3, bo=1) should resolve")
        if vanadyl is not None:
            self.assertGreater(alkoxide[0] - vanadyl[0], 0.10,
                               "V-O single bond should be longer than V=O double")


@unittest.skipUnless(SUP_3D_4D_PATH.exists() and MAIN_LIB_PATH.exists(),
                     "3d/4d supplement / main library not present")
class TestOffByteIdentical(unittest.TestCase):
    """When the env-flag is unset, the 3d/4d supplement is invisible."""

    def setUp(self):
        self._old_env = os.environ.pop("DELFIN_GRIP_MD_BO_SUPPLEMENT", None)
        import importlib
        import delfin.fffree.grip_mogul_lookup as _m
        importlib.reload(_m)
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        self.lib = GripLibrary(MAIN_LIB_PATH)

    def tearDown(self):
        if self._old_env is not None:
            os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = self._old_env

    def test_off_byte_identical(self):
        """env-flag unset => n_md_bo == 0 and lookups return None."""
        self.assertEqual(self.lib.n_md_bo, 0,
                         "env-flag unset must zero out the supplement")
        self.assertEqual(self.lib._md_bo_supplements, [],
                         "no supplements should be loaded")
        # Sample lookups must all return None.
        self.assertIsNone(self.lib._md_bo_lookup("Cr", "O", "sp2", "2", min_n=5))
        self.assertIsNone(self.lib._md_bo_lookup("Fe", "N", "sp3", "1", min_n=5))


@unittest.skipUnless(
    SUP_3D_4D_PATH.exists() and SUP_5D_MID_PATH.exists() and MAIN_LIB_PATH.exists(),
    "3d/4d + 5d-mid supplement / main library not present",
)
class TestDualSupplementCascade(unittest.TestCase):
    """Loading BOTH supplements (3d/4d + 5d-mid) covers the full transition
    series.  The lookup must pick the right supplement for each metal.
    """

    def setUp(self):
        self._old_env = os.environ.get("DELFIN_GRIP_MD_BO_SUPPLEMENT")
        # Colon-separated: 5d-mid first, 3d/4d second (order respected by lookup).
        os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = (
            f"{SUP_5D_MID_PATH}{os.pathsep}{SUP_3D_4D_PATH}"
        )
        import importlib
        import delfin.fffree.grip_mogul_lookup as _m
        importlib.reload(_m)
        from delfin.fffree.grip_mogul_lookup import GripLibrary
        self.lib = GripLibrary(MAIN_LIB_PATH)

    def tearDown(self):
        if self._old_env is None:
            os.environ.pop("DELFIN_GRIP_MD_BO_SUPPLEMENT", None)
        else:
            os.environ["DELFIN_GRIP_MD_BO_SUPPLEMENT"] = self._old_env

    def test_both_supplements_loaded(self):
        self.assertEqual(len(self.lib._md_bo_supplements), 2,
                         "exactly two supplements should be loaded")
        # Metals union must contain BOTH 5d-mid and 3d/4d targets.
        for m in ("W", "Mo", "Tc", "Re"):
            self.assertIn(m, self.lib._md_bo_metals,
                          f"5d-mid metal {m} must appear in union")
        for m in ("Ti", "V", "Cr", "Fe", "Pd", "Ir", "Pt"):
            self.assertIn(m, self.lib._md_bo_metals,
                          f"3d/4d metal {m} must appear in union")

    def test_w_via_5d_mid_supplement(self):
        """W-N(sp3, bo=1) is amine, should resolve via the 5d-mid supplement."""
        hit = self.lib._md_bo_lookup("W", "N", "sp3", "1", min_n=10)
        self.assertIsNotNone(hit, "W-N(sp3, bo=1) must resolve when both loaded")
        self.assertGreater(hit[0], 2.0)
        self.assertLess(hit[0], 2.4)

    def test_cr_via_3d_4d_supplement(self):
        """Cr=O via the 3d/4d supplement -- not covered by 5d-mid."""
        oxo = self.lib._md_bo_lookup("Cr", "O", "sp2", "2", min_n=10)
        self.assertIsNotNone(oxo, "Cr=O(sp2, bo=2) must resolve when both loaded")
        self.assertGreater(oxo[0], 1.50)
        self.assertLess(oxo[0], 1.80)


if __name__ == "__main__":
    unittest.main()
