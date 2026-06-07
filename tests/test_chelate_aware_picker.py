"""Tests for the chelate-aware polyhedron picker (2026-06-07).

Validates the per-class behaviour described in the project mission
``DELFIN universal mission: chelate-aware polyhedron picker``:

  * 2× bidentate-NN at CN=4 → picks the geometry whose vertex-vertex
    angles best match the chelate's intrinsic 5-ring bite (~78°).
  * 1× bidentate-NN + 2 monodentate at CN=4 → picks the geometry where
    the bidentate-NN pair lands on the cis-vertex.
  * No chelate ligand → returns the first candidate (legacy first-rule).
  * Env-flag OFF → byte-identical with legacy first-candidate rule.

The contract is geometric (no SMILES, no per-class branches), so the
tests use synthetic ``chelate_info`` dicts directly — bypassing the
SMILES → ``decompose`` path so the unit-level contract is isolated from
the heavier integration.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import os
import unittest
from typing import Dict, List

import numpy as np

from delfin.fffree.chelate_aware_picker import (
    RING_SIZE_TO_BITE_DEG,
    chelate_aware_picker_enabled,
    chelate_info_from_decompose,
    pick_polyhedron_chelate_aware,
    polyhedron_vertex_angles,
    ring_size_to_bite_angle,
    score_polyhedron_for_chelates,
)


class _EnvSnapshotMixin:
    """Save/restore env state around each test."""

    _SAVED_KEYS = (
        "DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK",
        "DELFIN_FFFREE_MOGUL_PRIMARY",
    )

    def setUp(self):  # noqa: D401
        super().setUp()
        self._saved = {k: os.environ.get(k) for k in self._SAVED_KEYS}

    def tearDown(self):  # noqa: D401
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        super().tearDown()


class TestRingSizeToBite(unittest.TestCase):
    """The empirical chelate-bite table is monotonic and covers 4-8 rings."""

    def test_table_monotonic_in_ring_size(self):
        sizes = sorted(RING_SIZE_TO_BITE_DEG.keys())
        bites = [RING_SIZE_TO_BITE_DEG[s] for s in sizes]
        self.assertEqual(sizes, sorted(sizes))
        self.assertEqual(bites, sorted(bites))

    def test_table_covers_4_to_8(self):
        for s in (4, 5, 6, 7, 8):
            self.assertIn(s, RING_SIZE_TO_BITE_DEG)

    def test_lookup_out_of_range_returns_none(self):
        self.assertIsNone(ring_size_to_bite_angle(2))
        self.assertIsNone(ring_size_to_bite_angle(20))

    def test_lookup_in_range(self):
        self.assertAlmostEqual(ring_size_to_bite_angle(5), 78.0, places=2)


class TestPolyhedronVertexAngles(unittest.TestCase):
    """Geometric sanity checks on the angle-matrix builder."""

    def test_sp4_has_only_90_and_180(self):
        angles = polyhedron_vertex_angles("SP-4 square planar")
        self.assertIsNotNone(angles)
        nondiag = angles[~np.eye(angles.shape[0], dtype=bool)]
        uniq = sorted(set(np.round(nondiag, 0)))
        self.assertEqual(uniq, [90.0, 180.0])

    def test_t4_has_only_109_5(self):
        angles = polyhedron_vertex_angles("T-4 tetrahedron")
        self.assertIsNotNone(angles)
        nondiag = angles[~np.eye(angles.shape[0], dtype=bool)]
        uniq = sorted(set(np.round(nondiag, 1)))
        self.assertEqual(uniq, [109.5])

    def test_tbp5_has_90_120_180(self):
        angles = polyhedron_vertex_angles("TBP-5 trigonal bipyramid")
        self.assertIsNotNone(angles)
        nondiag = angles[~np.eye(angles.shape[0], dtype=bool)]
        uniq = sorted(set(np.round(nondiag, 0)))
        self.assertEqual(uniq, [90.0, 120.0, 180.0])

    def test_spy5_has_no_120(self):
        """SPY-5 vertex angles are ~82°, ~89°, ~164° — no 120° pair.

        Critical for the AXOKAZ-style picker decision: a tridentate-mer
        chelate (two 5-ring bites of ~78°) lands closer to SPY-5's 82°
        basal-edges than to TBP-5's 120° equatorial-equatorial pair.
        """
        angles = polyhedron_vertex_angles("SPY-5 square pyramid")
        self.assertIsNotNone(angles)
        nondiag = angles[~np.eye(angles.shape[0], dtype=bool)]
        self.assertFalse(np.any(np.isclose(nondiag, 120.0, atol=1.0)))

    def test_unknown_geometry_returns_none(self):
        self.assertIsNone(polyhedron_vertex_angles("BOGUS-99 nonsense"))


class TestScoreFunction(unittest.TestCase):
    """Score = MSE between picked vertex pair angles and chelate bites."""

    def test_empty_bites_zero_score(self):
        for g in ("T-4 tetrahedron", "SP-4 square planar", "OC-6 octahedron"):
            self.assertEqual(score_polyhedron_for_chelates(g, 4, []), 0.0)

    def test_sp4_beats_t4_for_5ring_chelate(self):
        # Single bidentate 5-ring chelate (78°).  SP-4 has 90° (dev 12)
        # while T-4 has 109.5° (dev 31.5).  SP-4 wins.
        s_sp4 = score_polyhedron_for_chelates("SP-4 square planar", 4, [78.0])
        s_t4 = score_polyhedron_for_chelates("T-4 tetrahedron", 4, [78.0])
        self.assertLess(s_sp4, s_t4)

    def test_spy5_beats_tbp5_for_2_5ring_chelates(self):
        # AXOKAZ-style: two 5-ring chelate bites (terpy mer).  SPY-5's
        # basal edges are ~82° (dev 4) vs TBP-5's 120° equatorial pair
        # (dev 42).  SPY-5 must score lower.
        s_spy5 = score_polyhedron_for_chelates(
            "SPY-5 square pyramid", 5, [78.0, 78.0],
        )
        s_tbp5 = score_polyhedron_for_chelates(
            "TBP-5 trigonal bipyramid", 5, [78.0, 78.0],
        )
        self.assertIsNotNone(s_spy5)
        self.assertIsNotNone(s_tbp5)
        self.assertLess(s_spy5, s_tbp5)


class TestPickPolyhedron(_EnvSnapshotMixin, unittest.TestCase):
    """End-to-end picker on the four user-eye CN cases."""

    def setUp(self):
        super().setUp()
        os.environ["DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK"] = "1"

    def test_pick_polyhedron_2_bidentate_picks_cis_td(self):
        """2 bidentate-NN with 5-ring bite at CN=4.

        Both T-4 and SP-4 host 2 disjoint vertex pairs.
        * T-4: all pairs 109.5° → MSE = 2 * (31.5)² = 1984.5 deg²
        * SP-4: 4 cis-pairs at 90° + 2 trans-pairs at 180°.  Greedy picks
          the two 90° pairs → MSE = 2 * (12)² = 288 deg² ⇒ wins.

        The picker must therefore return SP-4 over T-4 even when the
        first-candidate (legacy) rule would have returned T-4.  Pass
        explicit candidates so the d8 default doesn't override.
        """
        chelate_info = [
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
        ]
        picked = pick_polyhedron_chelate_aware(
            cn=4,
            chelate_info=chelate_info,
            metal_sym="",
            # Explicit candidates: T-4 first, so legacy rule would pick it.
            candidates=["T-4 tetrahedron", "SP-4 square planar"],
        )
        self.assertEqual(picked, "SP-4 square planar")

    def test_pick_polyhedron_1_bidentate_2_monodentate(self):
        """1 bidentate-NN + 2 monodentate ligands at CN=4.

        Only ONE chelate-bite-constrained vertex pair.  SP-4 with one
        90° pair (dev 12) wins over T-4 with one 109.5° pair (dev 31.5).
        """
        chelate_info = [
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
        ]
        picked = pick_polyhedron_chelate_aware(
            cn=4,
            chelate_info=chelate_info,
            metal_sym="",
            candidates=["T-4 tetrahedron", "SP-4 square planar"],
        )
        self.assertEqual(picked, "SP-4 square planar")

    def test_pick_polyhedron_axokaz_tridentate_mer_picks_spy5(self):
        """AXOKAZ-style: tridentate-mer + 2 monodentate at CN=5.

        Two adjacent 5-ring chelate bites (~78°).  SPY-5 fits both at
        ~82° (dev 4); TBP-5 forces one to 120° (dev 42).  SPY-5 wins.
        """
        chelate_info = [
            {"bite_deg": 78.0, "denticity": 3, "ring_size": 5},
            {"bite_deg": 78.0, "denticity": 3, "ring_size": 5},
        ]
        picked = pick_polyhedron_chelate_aware(
            cn=5,
            chelate_info=chelate_info,
            metal_sym="",
            candidates=["TBP-5 trigonal bipyramid", "SPY-5 square pyramid"],
        )
        self.assertEqual(picked, "SPY-5 square pyramid")

    def test_pick_polyhedron_no_chelate_unchanged(self):
        """No chelate ligand → picker returns the first candidate."""
        for chelate_info in (None, [], [{"bite_deg": None}]):
            picked = pick_polyhedron_chelate_aware(
                cn=4,
                chelate_info=chelate_info,
                metal_sym="",
                candidates=["T-4 tetrahedron", "SP-4 square planar"],
            )
            self.assertEqual(picked, "T-4 tetrahedron")


class TestEnvFlagOff(_EnvSnapshotMixin, unittest.TestCase):
    """Env-flag OFF: picker returns first-candidate (legacy behaviour)."""

    def test_off_returns_first_candidate(self):
        os.environ["DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK"] = "0"
        # 2 strong chelate bites — would otherwise force SP-4 over T-4.
        chelate_info = [
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
        ]
        picked = pick_polyhedron_chelate_aware(
            cn=4,
            chelate_info=chelate_info,
            metal_sym="",
            candidates=["T-4 tetrahedron", "SP-4 square planar"],
        )
        self.assertEqual(picked, "T-4 tetrahedron")

    def test_default_off_without_mogul_primary(self):
        # No env-flags set at all → picker disabled, first candidate wins.
        os.environ.pop("DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        self.assertFalse(chelate_aware_picker_enabled())
        chelate_info = [
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
            {"bite_deg": 78.0, "denticity": 2, "ring_size": 5},
        ]
        picked = pick_polyhedron_chelate_aware(
            cn=4,
            chelate_info=chelate_info,
            metal_sym="",
            candidates=["T-4 tetrahedron", "SP-4 square planar"],
        )
        self.assertEqual(picked, "T-4 tetrahedron")

    def test_default_on_under_mogul_primary(self):
        # MOGUL_PRIMARY on + picker unset → picker auto-enabled.
        os.environ.pop("DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK", None)
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        self.assertTrue(chelate_aware_picker_enabled())

    def test_explicit_off_overrides_mogul_primary(self):
        os.environ["DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK"] = "0"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        self.assertFalse(chelate_aware_picker_enabled())


class TestChelateInfoFromDecompose(_EnvSnapshotMixin, unittest.TestCase):
    """Integration with ``decompose``: ring-size detection from graph.

    Uses simple synthetic SMILES (ethylenediamine-Ni, etc.) so the test
    is independent of the heavy CCDC SMILES.  RDKit is the only external
    dependency.
    """

    def setUp(self):
        super().setUp()
        self._pt3_saved = os.environ.get("DELFIN_FFFREE_PURE_TRACK3")
        os.environ["DELFIN_FFFREE_PURE_TRACK3"] = "1"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        os.environ["DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK"] = "1"

    def tearDown(self):
        if self._pt3_saved is None:
            os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        else:
            os.environ["DELFIN_FFFREE_PURE_TRACK3"] = self._pt3_saved
        super().tearDown()

    def test_en2_ni_picks_sp4_via_decompose(self):
        """Ni(en)2 (2× 5-ring chelate at CN=4, d8 → SP-4 default).

        The default geometry for d8-Ni CN=4 is already SP-4 from the
        ``_default_geometry`` lookup, so the picker should KEEP SP-4
        (the strong-bite chelate-pair matches the SP-4 90° cis-edge).
        Tests the chelate_info_from_decompose extraction on a real
        decompose result.
        """
        from delfin.fffree import decompose as _DEC
        smi = "[Ni]12(NCCN1)NCCN2"
        d = _DEC.decompose(smi)
        if d is None:
            self.skipTest("decompose returned None for synthetic en2-Ni SMILES")
        self.assertEqual(d.get("metal"), "Ni")
        self.assertTrue(d.get("has_chelate"))
        info = chelate_info_from_decompose(d)
        # Each chelate emits ONE bite (denticity-1).  Two chelates → 2 bites.
        self.assertEqual(len(info), 2)
        for c in info:
            self.assertEqual(c["ring_size"], 5)
            self.assertAlmostEqual(c["bite_deg"], 78.0, places=2)
        picked = pick_polyhedron_chelate_aware(
            cn=int(d["cn"]),
            chelate_info=info,
            metal_sym=d["metal"],
        )
        # d8-Ni CN=4 default is SP-4 — picker should keep it (5-ring bite
        # matches SP-4's 90° cis edge better than T-4's 109.5°).
        self.assertEqual(picked, "SP-4 square planar")


if __name__ == "__main__":
    unittest.main()
