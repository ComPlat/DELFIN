"""Tests for ``delfin.fffree.main_group_polyhedron`` — VSEPR-aware
polyhedron picker for main-group / p-block metals with stereo-active
lone pair (2026-06-07).

Validates:
  * LP detection rule per element (Sn(II) yes, Sn(IV) no, etc.).
  * Picker selects the correct LP-aware polyhedron for the affected
    (metal, CN, ox state) cases:
      - In(I) CN=1 → LP1-1, CN=2 → LP2-2 bent.
      - Sn(II) CN=3 → LP3-3 pyramidal, CN=4 → LP4-4 see-saw.
      - Sb(III)/Bi(III) CN=3 → LP3-3 pyramidal.
      - Po(IV) CN=5 → LP5-5 square-pyramid.
  * Vertex-vertex angle invariants for each LP-aware polyhedron:
      - LP3-3 D-D ≈ 109.5° (tetrahedral-pyramid, not 120° planar).
      - LP4-4 axial-axial ≈ 180°, equatorial-equatorial ≈ 120°,
        axial-equatorial ≈ 90°.
      - LP5-5 base-base ≈ 90° or 180°, apex-base ≈ 90°.
  * Sn(IV), Pb(IV), Al(III), Ga(III) bypass the picker entirely (no LP).
  * Env-flag OFF → byte-identical with legacy ``polyhedra.geometries_for_cn``
    (no LP-aware names injected, decompose returns the legacy geometry).
  * 5 validation cases from the user-supplied DG-CLEAN smoke 100 outputs
    (D-IHIVOL, D-ILEDAH, D-TENMIL, D-XEXGOZ + a synthetic Sb(III) case).

Universal contract: no SMILES patterns, no per-class branches.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import importlib
import math
import os
import unittest
from typing import Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Env-snapshot mixin — save/restore env so tests can flip flags freely.
# ---------------------------------------------------------------------------


class _EnvSnapshotMixin:
    _SAVED_KEYS = (
        "DELFIN_FFFREE_MAIN_GROUP_LP",
        "DELFIN_FFFREE_MOGUL_PRIMARY",
        "DELFIN_FFFREE_PURE_TRACK3",
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
        # Reset cached modules so reload picks up the new env.
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            mod = importlib.import_module(m)
            importlib.reload(mod)
        super().tearDown()


def _ang(a, b) -> float:
    """Angle (deg) between two 3-vectors."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    c = float(np.dot(a, b) / (na * nb))
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


# ---------------------------------------------------------------------------
# LP-detection rule
# ---------------------------------------------------------------------------


class TestHasStereoActiveLP(_EnvSnapshotMixin, unittest.TestCase):
    """Periodic-table LP-presence rule."""

    def test_sn_ii_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertTrue(has_stereo_active_lp("Sn", 2))

    def test_sn_iv_has_no_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertFalse(has_stereo_active_lp("Sn", 4))

    def test_pb_ii_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertTrue(has_stereo_active_lp("Pb", 2))

    def test_pb_iv_has_no_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertFalse(has_stereo_active_lp("Pb", 4))

    def test_sb_iii_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertTrue(has_stereo_active_lp("Sb", 3))

    def test_sb_v_has_no_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertFalse(has_stereo_active_lp("Sb", 5))

    def test_bi_iii_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertTrue(has_stereo_active_lp("Bi", 3))

    def test_in_i_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertTrue(has_stereo_active_lp("In", 1))

    def test_in_iii_has_no_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertFalse(has_stereo_active_lp("In", 3))

    def test_tl_i_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertTrue(has_stereo_active_lp("Tl", 1))

    def test_al_never_has_lp(self):
        """Al has no LP at any oxidation state (light p-block)."""
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        for ox in (1, 2, 3, 4, 5):
            self.assertFalse(has_stereo_active_lp("Al", ox))

    def test_ga_never_has_lp(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        for ox in (1, 2, 3, 4, 5):
            self.assertFalse(has_stereo_active_lp("Ga", ox))

    def test_fe_never_has_lp(self):
        """d-block: never carries a stereo-active LP via the main-group rule."""
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        for ox in (0, 1, 2, 3, 4, 5, 6):
            self.assertFalse(has_stereo_active_lp("Fe", ox))

    def test_no_ox_state_returns_false(self):
        from delfin.fffree.main_group_polyhedron import has_stereo_active_lp
        self.assertFalse(has_stereo_active_lp("Sn", None))


# ---------------------------------------------------------------------------
# Picker — geometry name for (metal, CN, ox state)
# ---------------------------------------------------------------------------


class TestMainGroupPolyhedronPicker(_EnvSnapshotMixin, unittest.TestCase):
    """The picker chooses the right LP-aware polyhedron."""

    def setUp(self):
        super().setUp()
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        # Reload modules so the gate is observed.
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))

    def test_sn_ii_cn3_picks_pyramidal(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Sn", 3, 2)
        self.assertEqual(g, "LP3-3 pyramidal lone-pair")

    def test_sn_ii_cn4_picks_seesaw(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Sn", 4, 2)
        self.assertEqual(g, "LP4-4 see-saw lone-pair")

    def test_sn_iv_cn4_picker_returns_none(self):
        """Sn(IV) has no LP → picker returns None → legacy T-4 path."""
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Sn", 4, 4)
        self.assertIsNone(g)

    def test_in_i_cn2_picks_bent(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("In", 2, 1)
        self.assertEqual(g, "LP2-2 bent lone-pair")

    def test_in_i_cn1_picks_linear_mono(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("In", 1, 1)
        self.assertEqual(g, "LP1-1 linear-mono lone-pair")

    def test_sb_iii_cn3_picks_pyramidal(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Sb", 3, 3)
        self.assertEqual(g, "LP3-3 pyramidal lone-pair")

    def test_bi_iii_cn3_picks_pyramidal(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Bi", 3, 3)
        self.assertEqual(g, "LP3-3 pyramidal lone-pair")

    def test_pb_ii_cn5_picks_square_pyramid(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Pb", 5, 2)
        self.assertEqual(g, "LP5-5 square-pyramid lone-pair")

    def test_al_iii_returns_none(self):
        """Al(III) — never an LP, picker returns None."""
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Al", 4, 3)
        self.assertIsNone(g)

    def test_fe_returns_none(self):
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        for cn in (3, 4, 5, 6):
            g = main_group_polyhedron_for("Fe", cn, 2)
            self.assertIsNone(g, f"Fe CN={cn} picked an LP geom: {g}")

    def test_cn6_with_lp_active_returns_none(self):
        """CN ≥ 6 — LP delocalised, fall through to legacy octahedron."""
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Sn", 6, 2)
        self.assertIsNone(g)

    def test_picker_with_unknown_ox_state_falls_back_to_single_lp(self):
        """When ox state is None and the metal has exactly one LP-active
        ox state, the picker assumes LP at low CN."""
        from delfin.fffree.main_group_polyhedron import main_group_polyhedron_for
        g = main_group_polyhedron_for("Sn", 3, None)
        self.assertEqual(g, "LP3-3 pyramidal lone-pair")
        # High CN without ox state → assume LP-inactive (high ox).
        g = main_group_polyhedron_for("Sn", 6, None)
        self.assertIsNone(g)


# ---------------------------------------------------------------------------
# Vertex-vertex angle invariants
# ---------------------------------------------------------------------------


class TestLPPolyhedronGeometry(_EnvSnapshotMixin, unittest.TestCase):
    """The LP-aware vertex arrays have the correct VSEPR angles."""

    def test_lp3_pyramidal_dd_angle_is_tetrahedral(self):
        """LP3-3 D-D angle ≈ 109.5° (NH3-like, NOT 120° trigonal-planar)."""
        from delfin.fffree.main_group_polyhedron import (
            effective_ref_vectors_main_group,
        )
        V = effective_ref_vectors_main_group("LP3-3 pyramidal lone-pair")
        self.assertEqual(V.shape, (3, 3))
        for i in range(3):
            for j in range(i + 1, 3):
                self.assertAlmostEqual(_ang(V[i], V[j]), 109.47, delta=0.5)

    def test_lp4_seesaw_axial_axial_180(self):
        from delfin.fffree.main_group_polyhedron import (
            effective_ref_vectors_main_group,
        )
        V = effective_ref_vectors_main_group("LP4-4 see-saw lone-pair")
        self.assertEqual(V.shape, (4, 3))
        # v0, v1 = axial (TBP parent vertices 0, 1).
        self.assertAlmostEqual(_ang(V[0], V[1]), 180.0, delta=1.0)
        # v2, v3 = equatorial (TBP parent vertices 2, 3) — at 120° spacing
        # because the LP sits at TBP vertex 4.
        self.assertAlmostEqual(_ang(V[2], V[3]), 120.0, delta=1.0)
        # ax–eq = 90°.
        self.assertAlmostEqual(_ang(V[0], V[2]), 90.0, delta=1.0)
        self.assertAlmostEqual(_ang(V[1], V[3]), 90.0, delta=1.0)

    def test_lp5_square_pyramid_apex_base_90(self):
        from delfin.fffree.main_group_polyhedron import (
            effective_ref_vectors_main_group,
        )
        V = effective_ref_vectors_main_group("LP5-5 square-pyramid lone-pair")
        self.assertEqual(V.shape, (5, 3))
        # Donors at parent OC-6 vertices [+x, -x, +y, -y, +z].  LP at -z.
        # apex (v4) to base (v0..v3) = 90°.
        for i in range(4):
            self.assertAlmostEqual(_ang(V[4], V[i]), 90.0, delta=1.0)
        # opposite base donors = 180° (+x ↔ -x, +y ↔ -y).
        self.assertAlmostEqual(_ang(V[0], V[1]), 180.0, delta=1.0)
        self.assertAlmostEqual(_ang(V[2], V[3]), 180.0, delta=1.0)
        # adjacent base donors = 90°.
        self.assertAlmostEqual(_ang(V[0], V[2]), 90.0, delta=1.0)

    def test_lp2_bent_dd_angle_is_120(self):
        """LP2-2 bent: parent SP-3 trigonal planar with LP at one vertex →
        donors at 0 and 1 give a 120° angle (SO2-like bent)."""
        from delfin.fffree.main_group_polyhedron import (
            effective_ref_vectors_main_group,
        )
        V = effective_ref_vectors_main_group("LP2-2 bent lone-pair")
        self.assertEqual(V.shape, (2, 3))
        self.assertAlmostEqual(_ang(V[0], V[1]), 120.0, delta=1.0)

    def test_all_lp_vectors_are_unit_norm(self):
        from delfin.fffree.main_group_polyhedron import (
            MAIN_GROUP_GEOM_BY_CN,
            effective_ref_vectors_main_group,
        )
        for cn, names in MAIN_GROUP_GEOM_BY_CN.items():
            for name in names:
                V = effective_ref_vectors_main_group(name)
                norms = np.linalg.norm(V, axis=1)
                np.testing.assert_allclose(
                    norms, np.ones(V.shape[0]), atol=1e-9,
                    err_msg=f"{name} norms {norms}"
                )

    def test_unknown_geometry_raises(self):
        from delfin.fffree.main_group_polyhedron import (
            effective_ref_vectors_main_group,
        )
        with self.assertRaises(KeyError):
            effective_ref_vectors_main_group("not-a-polyhedron")


# ---------------------------------------------------------------------------
# Integration: polyhedra.geometries_for_cn injects LP-aware names
# ---------------------------------------------------------------------------


class TestPolyhedraIntegration(_EnvSnapshotMixin, unittest.TestCase):
    """``polyhedra.geometries_for_cn`` exposes LP-aware names under ON."""

    def test_off_byte_identical_for_sn(self):
        """OFF → Sn CN=3 returns ['SP-3 trigonal planar', 'T-3 T-shape']."""
        os.environ.pop("DELFIN_FFFREE_MAIN_GROUP_LP", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.polyhedra import geometries_for_cn
        g = geometries_for_cn(3, "Sn")
        self.assertEqual(g, ["SP-3 trigonal planar", "T-3 T-shape"])
        g = geometries_for_cn(2, "In")
        self.assertEqual(g, ["L-2 linear"])

    def test_off_byte_identical_for_d_block(self):
        os.environ.pop("DELFIN_FFFREE_MAIN_GROUP_LP", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.polyhedra import geometries_for_cn
        # d-block untouched.
        self.assertEqual(
            geometries_for_cn(6, "Fe"),
            ["OC-6 octahedron", "TPR-6 trigonal prism"],
        )
        self.assertEqual(
            geometries_for_cn(4, "Pt"),
            ["T-4 tetrahedron", "SP-4 square planar"],
        )

    def test_on_state_lp_prepended_for_sn_cn3(self):
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.polyhedra import geometries_for_cn
        g = geometries_for_cn(3, "Sn")
        self.assertEqual(g[0], "LP3-3 pyramidal lone-pair")

    def test_on_state_d_block_unchanged(self):
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.polyhedra import geometries_for_cn
        # Fe / Pt — the per-CN list is metal-agnostic (the d8-vs-Td choice
        # for Pt happens in ``decompose._default_geometry``, not here).
        # The picker MUST NOT inject any LP-aware names for d-block input.
        for m in ("Fe", "Pt", "Cu", "Pd"):
            for cn in (3, 4, 5, 6):
                cands = geometries_for_cn(cn, m)
                for c in cands:
                    self.assertFalse(
                        str(c).startswith("LP"),
                        f"d-block metal {m} CN={cn} got LP candidate {c!r}",
                    )

    def test_ref_vectors_dispatches_lp_names(self):
        """polyhedra.ref_vectors resolves LP-aware names to donor-only vectors."""
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.polyhedra import ref_vectors
        V = ref_vectors("LP3-3 pyramidal lone-pair")
        self.assertEqual(V.shape, (3, 3))
        # All unit norm.
        np.testing.assert_allclose(
            np.linalg.norm(V, axis=1), np.ones(3), atol=1e-9
        )


# ---------------------------------------------------------------------------
# Integration: decompose._default_geometry uses LP-aware names
# ---------------------------------------------------------------------------


class TestDecomposeIntegration(_EnvSnapshotMixin, unittest.TestCase):
    """``decompose._default_geometry`` returns the LP-aware polyhedron."""

    def test_off_legacy_unchanged_for_sn(self):
        os.environ.pop("DELFIN_FFFREE_MAIN_GROUP_LP", None)
        os.environ.pop("DELFIN_FFFREE_MOGUL_PRIMARY", None)
        os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.decompose import _default_geometry
        self.assertEqual(_default_geometry("Sn", 3), "SP-3 trigonal planar")
        self.assertEqual(_default_geometry("In", 2), "L-2 linear")

    def test_on_lp_replaces_legacy_for_sn(self):
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.decompose import _default_geometry
        self.assertEqual(
            _default_geometry("Sn", 3), "LP3-3 pyramidal lone-pair"
        )
        self.assertEqual(
            _default_geometry("In", 2), "LP2-2 bent lone-pair"
        )
        self.assertEqual(
            _default_geometry("Bi", 3), "LP3-3 pyramidal lone-pair"
        )
        self.assertEqual(
            _default_geometry("Sn", 4), "LP4-4 see-saw lone-pair"
        )

    def test_on_state_high_cn_falls_through_to_legacy(self):
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.decompose import _default_geometry
        # Sn CN=6 → no LP-aware entry → legacy OC-6.
        self.assertEqual(_default_geometry("Sn", 6), "OC-6 octahedron")


# ---------------------------------------------------------------------------
# Default ON under MOGUL_PRIMARY gate
# ---------------------------------------------------------------------------


class TestMogulPrimaryAutoOn(_EnvSnapshotMixin, unittest.TestCase):
    """When MOGUL_PRIMARY=1 and the LP gate is unset, the LP picker is auto-on."""

    def test_auto_on_under_mogul_primary(self):
        os.environ.pop("DELFIN_FFFREE_MAIN_GROUP_LP", None)
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.main_group_polyhedron import main_group_lp_enabled
        self.assertTrue(main_group_lp_enabled())

    def test_explicit_off_overrides_mogul_primary(self):
        """The explicit env-flag wins over the MOGUL_PRIMARY auto-on rule."""
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "0"
        os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))
        from delfin.fffree.main_group_polyhedron import main_group_lp_enabled
        self.assertFalse(main_group_lp_enabled())


# ---------------------------------------------------------------------------
# 5-case validation from the smoke 100 main-group fails
# ---------------------------------------------------------------------------


class TestSmoke100ValidationCases(_EnvSnapshotMixin, unittest.TestCase):
    """Validate the 5 main-group cases from DG-CLEAN smoke 100.

    Each case names the refcode + metal + observed CN that the smoke run
    showed as over-coordinated.  The LP-aware picker MUST return a geometry
    name whose donor-vertex angles differ from the legacy choice (so the
    bounds matrix will pull donors to the correct distorted positions).
    """

    def setUp(self):
        super().setUp()
        os.environ["DELFIN_FFFREE_MAIN_GROUP_LP"] = "1"
        for m in (
            "delfin.fffree.main_group_polyhedron",
            "delfin.fffree.polyhedra",
            "delfin.fffree.decompose",
        ):
            importlib.reload(importlib.import_module(m))

    def test_d_ihivol_in_cn2(self):
        """D-IHIVOL: In(I)-Ru complex, In has CN=2 σ-donors.

        Legacy: L-2 linear (180° between donors) → embed over-fills the
        M-shell because the LP space is left empty and the optimiser
        pulls a non-donor close.  Fix: LP2-2 bent (120° between donors)
        leaves the LP vertex at the trigonal-planar third site and the
        donors stay distinct.
        """
        from delfin.fffree.decompose import _default_geometry
        g = _default_geometry("In", 2)
        self.assertEqual(g, "LP2-2 bent lone-pair")

    def test_d_iledah_sn_cn1(self):
        """D-ILEDAH: Sn(II)-Ru complex, Sn has CN=1 σ-donor.

        Legacy: no CN=1 polyhedron at all → fall back triggers M-shell
        over-fill (gets 4).  Fix: LP1-1 linear-mono lone-pair gives an
        axial Sn-D pair on a trigonal-planar frame with 2 LPs.
        """
        from delfin.fffree.main_group_polyhedron import (
            main_group_polyhedron_for,
        )
        g = main_group_polyhedron_for("Sn", 1, 2)
        self.assertEqual(g, "LP1-1 linear-mono lone-pair")

    def test_d_tenmil_sn_cn3(self):
        """D-TENMIL: Sn(II), CN=3.

        Legacy: SP-3 trigonal planar (120° flat) → embed places all 3
        donors in a plane that doesn't exist physically.  Fix: LP3-3
        pyramidal (109.5° at metal, NH3-like).
        """
        from delfin.fffree.decompose import _default_geometry
        g = _default_geometry("Sn", 3)
        self.assertEqual(g, "LP3-3 pyramidal lone-pair")

    def test_d_xexgoz_in_cn1(self):
        """D-XEXGOZ: In(I)-Ru, In CN=1 σ-donor.  Same fix as ILEDAH."""
        from delfin.fffree.main_group_polyhedron import (
            main_group_polyhedron_for,
        )
        g = main_group_polyhedron_for("In", 1, 1)
        self.assertEqual(g, "LP1-1 linear-mono lone-pair")

    def test_synthetic_sb_iii_cn3(self):
        """Synthetic Sb(III) CN=3 case.

        Sb(III) is the classic VSEPR pyramidal example (SbX3 like SbF3).
        Legacy: SP-3 trigonal planar.  Fix: LP3-3 pyramidal.
        """
        from delfin.fffree.decompose import _default_geometry
        g = _default_geometry("Sb", 3)
        self.assertEqual(g, "LP3-3 pyramidal lone-pair")


# ---------------------------------------------------------------------------
# Helpers exported by the module
# ---------------------------------------------------------------------------


class TestModuleSurface(_EnvSnapshotMixin, unittest.TestCase):
    def test_parent_and_donor_vertex_lookup(self):
        from delfin.fffree.main_group_polyhedron import (
            parent_polyhedron,
            donor_vertex_indices,
            lp_vertex_index,
        )
        self.assertEqual(
            parent_polyhedron("LP3-3 pyramidal lone-pair"), "T-4 tetrahedron"
        )
        self.assertEqual(
            donor_vertex_indices("LP3-3 pyramidal lone-pair"), (0, 1, 2)
        )
        self.assertEqual(
            lp_vertex_index("LP3-3 pyramidal lone-pair"), 3
        )
        self.assertEqual(
            parent_polyhedron("LP4-4 see-saw lone-pair"),
            "TBP-5 trigonal bipyramid",
        )
        self.assertEqual(
            lp_vertex_index("LP4-4 see-saw lone-pair"), 4
        )

    def test_aliases_resolve(self):
        from delfin.fffree.main_group_polyhedron import (
            effective_ref_vectors_main_group,
        )
        V1 = effective_ref_vectors_main_group("LP3-3")
        V2 = effective_ref_vectors_main_group("LP3-3 pyramidal lone-pair")
        np.testing.assert_array_equal(V1, V2)

    def test_main_group_lp_metal_classification(self):
        from delfin.fffree.main_group_polyhedron import is_main_group_lp_metal
        self.assertTrue(is_main_group_lp_metal("Sn"))
        self.assertTrue(is_main_group_lp_metal("Pb"))
        self.assertTrue(is_main_group_lp_metal("Sb"))
        self.assertTrue(is_main_group_lp_metal("Bi"))
        self.assertTrue(is_main_group_lp_metal("In"))
        self.assertTrue(is_main_group_lp_metal("Tl"))
        # NOT in table: Al, Ga, Fe (no stereo-active LP at relevant ox).
        self.assertFalse(is_main_group_lp_metal("Al"))
        self.assertFalse(is_main_group_lp_metal("Ga"))
        self.assertFalse(is_main_group_lp_metal("Fe"))
        self.assertFalse(is_main_group_lp_metal("Cu"))


if __name__ == "__main__":
    unittest.main()
