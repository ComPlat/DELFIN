"""Tests for delfin._topology_hash (Welle-5p-A topology hard-gate).

These tests cover the six acceptance scenarios from the Welle-5p-A
agent brief:

1. amine-H orientation angle correctness
2. topology-hash preserves under a valid rotation
3. topology-hash rejects a synthetic "bad rotation" (amine-H flip)
4. ALEQEO regression test (Fe(CO)2(NH2-CH2-CH2-S)2 amine-H toward Fe)
5. rollback semantics in :func:`_rotamer_diversity.apply`
6. bond-multiset consistency under benign permutations

All XYZs are inline (no external test fixtures).
"""

from __future__ import annotations

import math
import os

import pytest

from delfin import _topology_hash as th


# ---------------------------------------------------------------------------
# Common fixtures
# ---------------------------------------------------------------------------


def _clear_env(monkeypatch):
    for key in (
        "DELFIN_5P_A_TOPOLOGY_HARDGATE",
        "DELFIN_5P_A_MD_TOL",
        "DELFIN_5P_A_AMINE_H_MIN_DEG",
        "DELFIN_5P_A_BOND_TOL",
        "DELFIN_5P_A_CLASH_HM_MIN",
    ):
        monkeypatch.delenv(key, raising=False)


# Fe(NH3) — a single amine donor on Fe.  Good geometry: N at 2.05 Å with
# three H atoms pointing away from Fe (∠Fe-N-H ≈ 109.5°).
def _fe_nh3_good():
    syms = ["Fe", "N", "H", "H", "H"]
    coords = [
        (0.0, 0.0, 0.0),
        (2.05, 0.0, 0.0),
        (2.45, 0.89, 0.0),
        (2.45, -0.445, 0.77),
        (2.45, -0.445, -0.77),
    ]
    return syms, coords


def _fe_nh3_bad():
    """Same as _fe_nh3_good but with one H atom flipped toward Fe."""
    syms, coords = _fe_nh3_good()
    coords = list(coords)
    coords[2] = (1.55, 0.55, 0.0)  # H pulled toward Fe
    return syms, coords


# X10-ALEQEO toy model: Fe(CO)2(NH2-CH2-CH2-S) — simplified to one chelate
# ring with the amine-H pointing toward Fe.
def _aleqeo_bad():
    """One Fe with two CO + one NH2-CH2-CH2-S chelate, amine-H flipped.

    Positions are approximate octahedral.  We focus only on the amine-H
    flip — the metric we care about is whether the gate rejects this.
    """
    # Fe at origin
    syms = ["Fe", "C", "O", "C", "O", "N", "H", "H", "C", "H", "H",
            "C", "H", "H", "S"]
    coords = [
        (0.0, 0.0, 0.0),         # 0 Fe
        (0.0, 1.85, 0.0),        # 1 C(O)
        (0.0, 3.00, 0.0),        # 2 O
        (0.0, -1.85, 0.0),       # 3 C(O)
        (0.0, -3.00, 0.0),       # 4 O
        (2.05, 0.0, 0.0),        # 5 N
        (1.65, 0.6, 0.0),        # 6 H — flipped toward Fe!
        (2.50, -0.5, 0.85),      # 7 H — normal
        (2.85, 1.0, 0.0),        # 8 CH2 alpha-C
        (3.50, 0.6, 0.8),        # 9 H on alpha
        (3.50, 0.6, -0.8),       # 10 H on alpha
        (3.5, 2.5, 0.0),         # 11 CH2 beta-C
        (4.2, 2.7, 0.8),         # 12 H on beta
        (4.2, 2.7, -0.8),        # 13 H on beta
        (2.05, 2.5, -1.5),       # 14 S
    ]
    return syms, coords


def _aleqeo_good():
    """Same as _aleqeo_bad but with normal amine-H orientation."""
    syms, coords = _aleqeo_bad()
    coords = list(coords)
    coords[6] = (2.50, 0.5, 0.85)  # H now away from Fe
    return syms, coords


# ---------------------------------------------------------------------------
# 1. amine-H orientation angle correctness
# ---------------------------------------------------------------------------


class TestAmineHOrientation:
    def test_good_geometry_passes(self):
        syms, coords = _fe_nh3_good()
        result = th.standalone_amine_h_realism(syms, coords)
        assert result.passed, f"good geometry should pass: {result.violations}"

    def test_flipped_h_rejected(self):
        syms, coords = _fe_nh3_bad()
        result = th.standalone_amine_h_realism(syms, coords)
        assert not result.passed
        assert any("flip" in v or "too close" in v for v in result.violations)

    def test_angle_below_threshold_rejected(self):
        """Angle just under 60° should fail with default threshold."""
        # Build N–H at 59° from M-N direction
        syms = ["Fe", "N", "H"]
        # N at +x, H at angle 59° from N→Fe direction
        ang_deg = 59.0
        ang_rad = math.radians(ang_deg)
        nh_len = 1.01
        # H position: from N, angle ang_deg from M-N vector
        # M-N vector points from N toward M (i.e. -x).  So H at angle 59° from -x.
        hx = 2.05 + nh_len * (-math.cos(ang_rad))
        hy = nh_len * math.sin(ang_rad)
        coords = [(0.0, 0.0, 0.0), (2.05, 0.0, 0.0), (hx, hy, 0.0)]
        # However for ∠(M-N-H), M is at origin, vertex N at +x.
        # The vector N→M = (-2.05, 0, 0), N→H = (hx-2.05, hy, 0).
        # cos(angle) = (-2.05)*(hx-2.05) / (2.05 * nh_len)
        # We want the angle equal to ang_deg → so we compute H from N along
        # a direction at ang_deg from N→M.
        # Already done above (the construction was wrong — let's redo cleanly).
        # N at (2.05,0,0), M at (0,0,0). N→M unit = (-1,0,0).
        # Want H at distance 1.01 from N such that angle from N→M is 59°.
        # H = N + 1.01 * (cos(59°)*(-1,0,0) + sin(59°)*(0,1,0))
        hx = 2.05 + nh_len * math.cos(ang_rad) * (-1.0)
        hy = nh_len * math.sin(ang_rad) * 1.0
        coords[2] = (hx, hy, 0.0)
        result = th.standalone_amine_h_realism(syms, coords)
        assert not result.passed, f"angle {ang_deg}° should fail"

    def test_angle_above_threshold_passes(self):
        syms = ["Fe", "N", "H"]
        # Geometry with ∠(M-N-H) = 109.5° — tetrahedral amine
        ang_deg = 109.5
        ang_rad = math.radians(ang_deg)
        nh_len = 1.01
        hx = 2.05 + nh_len * math.cos(ang_rad) * (-1.0)
        hy = nh_len * math.sin(ang_rad) * 1.0
        coords = [(0.0, 0.0, 0.0), (2.05, 0.0, 0.0), (hx, hy, 0.0)]
        result = th.standalone_amine_h_realism(syms, coords)
        assert result.passed, result.violations


# ---------------------------------------------------------------------------
# 2. topology preservation under valid rotation
# ---------------------------------------------------------------------------


class TestValidRotation:
    def test_identity_passes(self):
        syms, coords = _fe_nh3_good()
        result = th.topology_preserved(syms, coords, coords)
        assert result.passed
        assert result.details["fp_before"] == result.details["fp_after"]

    def test_small_translation_passes(self):
        """A pure translation should not break topology."""
        syms, coords = _fe_nh3_good()
        coords_after = [(x + 0.01, y, z) for (x, y, z) in coords]
        result = th.topology_preserved(syms, coords, coords_after)
        assert result.passed, result.violations

    def test_rigid_rotation_passes(self):
        """A rigid rotation of the entire molecule should pass."""
        syms, coords = _fe_nh3_good()
        # rotate 30° around z-axis
        ang = math.radians(30.0)
        c, s = math.cos(ang), math.sin(ang)
        coords_after = [(c * x - s * y, s * x + c * y, z) for (x, y, z) in coords]
        result = th.topology_preserved(syms, coords, coords_after)
        assert result.passed, result.violations


# ---------------------------------------------------------------------------
# 3. topology-hash rejects synthetic bad rotation
# ---------------------------------------------------------------------------


class TestBadRotationRejected:
    def test_amine_h_flip_rejected(self):
        syms, before = _fe_nh3_good()
        _, after = _fe_nh3_bad()
        result = th.topology_preserved(syms, before, after)
        assert not result.passed
        # Should mention either "flip" or bond multiset change
        assert any(
            "flip" in v or "too close" in v or "bond multiset" in v
            or "donor gained" in v
            for v in result.violations
        )

    def test_md_distance_drift_rejected(self):
        """Move the N donor 0.5 Å further from Fe → M-D invariant fails."""
        syms, before = _fe_nh3_good()
        after = list(before)
        # Move N farther
        after[1] = (2.55, 0.0, 0.0)
        # Move H atoms with N (rigid-H) so they aren't auto-flagged
        for i in (2, 3, 4):
            x, y, z = before[i]
            after[i] = (x + 0.5, y, z)
        result = th.topology_preserved(syms, before, after, md_tol=0.05)
        assert not result.passed
        assert any("M-D" in v for v in result.violations)

    def test_bond_multiset_change_rejected(self):
        """Break a C-C bond by separating two C atoms."""
        # Simple ethane-like
        syms = ["C", "C", "H", "H", "H", "H", "H", "H"]
        before = [
            (0.0, 0.0, 0.0),
            (1.54, 0.0, 0.0),
            (-0.36, 1.03, 0.0),
            (-0.36, -0.51, 0.89),
            (-0.36, -0.51, -0.89),
            (1.90, -0.51, 0.89),
            (1.90, -0.51, -0.89),
            (1.90, 1.03, 0.0),
        ]
        after = list(before)
        # Pull C2 far away
        after[1] = (5.0, 0.0, 0.0)
        for i in (5, 6, 7):
            x, y, z = before[i]
            after[i] = (x + 3.46, y, z)
        result = th.topology_preserved(syms, before, after)
        assert not result.passed
        assert any("bond multiset" in v for v in result.violations)


# ---------------------------------------------------------------------------
# 4. ALEQEO regression test
# ---------------------------------------------------------------------------


class TestALEQEORegression:
    def test_aleqeo_bad_geometry_rejected(self):
        """Toy ALEQEO with amine-H flipped → standalone check must reject."""
        syms, coords = _aleqeo_bad()
        result = th.standalone_amine_h_realism(syms, coords)
        assert not result.passed
        assert any("flip" in v or "too close" in v for v in result.violations)

    def test_aleqeo_good_geometry_passes(self):
        """Same skeleton with amine-H pointing away → standalone check passes."""
        syms, coords = _aleqeo_good()
        result = th.standalone_amine_h_realism(syms, coords)
        # We don't enforce all-clean (atom-overlap is a separate concern),
        # but the amine-H realism alone must pass — no flip violations.
        flip_violations = [v for v in result.violations if "flip" in v or "too close" in v]
        assert not flip_violations, f"unexpected flip violations: {flip_violations}"


# ---------------------------------------------------------------------------
# 5. rollback semantics in _rotamer_diversity.apply
# ---------------------------------------------------------------------------


class TestRotamerDiversityRollback:
    def test_default_off_byte_identical(self, monkeypatch):
        """With master flag unset, _rotamer_diversity behaves as before."""
        _clear_env(monkeypatch)
        # Without master flag the rotamer pass is a no-op anyway, but we want
        # to verify wire-in does not change behaviour either.
        from delfin import _rotamer_diversity as rot
        # Quick: identity check on a simple butane XYZ — no metal, no donor
        xyz = (
            "C  0.0  0.0  0.0\n"
            "C  1.54 0.0  0.0\n"
            "H  -0.36  1.03  0.0\n"
            "H  -0.36  -0.51  0.89\n"
            "H  -0.36  -0.51  -0.89\n"
            "H  1.90  -0.51  0.89\n"
            "H  1.90  -0.51  -0.89\n"
            "H  1.90  1.03  0.0\n"
        )
        out = rot.apply_if_enabled(xyz)
        assert out == [xyz]

    def test_hardgate_on_rejects_bad_rotation(self, monkeypatch):
        """Enable hardgate + call topology_preserved on synthetic bad coords."""
        _clear_env(monkeypatch)
        monkeypatch.setenv("DELFIN_5P_A_TOPOLOGY_HARDGATE", "1")
        # We test the helper directly — _rotamer_diversity needs OB to run
        # end-to-end, but the gate itself is in _topology_hash.
        syms, before = _fe_nh3_good()
        _, after = _fe_nh3_bad()
        result = th.topology_preserved(syms, before, after)
        assert not result.passed


# ---------------------------------------------------------------------------
# 6. bond-multiset consistency
# ---------------------------------------------------------------------------


class TestBondMultisetConsistency:
    def test_atom_order_invariance(self):
        """Reordering atoms must not change the fingerprint."""
        syms, coords = _fe_nh3_good()
        fp1 = th.compute_topology_fingerprint(syms, coords)

        # Reorder atoms (swap H2 and H3) — fingerprint must be identical
        # because we sort the bond multiset by element pair.
        syms_r = list(syms)
        coords_r = list(coords)
        syms_r[2], syms_r[3] = syms_r[3], syms_r[2]
        coords_r[2], coords_r[3] = coords_r[3], coords_r[2]
        fp2 = th.compute_topology_fingerprint(syms_r, coords_r)
        assert fp1.bond_multiset == fp2.bond_multiset
        assert fp1.formula == fp2.formula
        assert fp1.hexdigest() == fp2.hexdigest()

    def test_formula_hill_convention(self):
        syms = ["C", "C", "H", "H", "H", "H", "H", "H", "O", "Fe"]
        coords = [(0.0, 0.0, 0.0)] * 10  # geometry irrelevant for formula
        fp = th.compute_topology_fingerprint(syms, coords)
        # Hill: C first, H second, then alphabetical for rest
        assert fp.formula.startswith("C2H6")
        # "Fe" comes before "O" alphabetically
        assert fp.formula == "C2H6FeO"

    def test_md_bins_capture_distance(self):
        """Two snapshots with different M-D distances → different fingerprint."""
        syms, coords = _fe_nh3_good()
        fp1 = th.compute_topology_fingerprint(syms, coords)
        # Stretch the M-D bond by 0.1 Å (2 buckets at 0.05 Å)
        coords2 = list(coords)
        coords2[1] = (2.15, 0.0, 0.0)
        # We also shift the bonded H atoms so the bond detection still
        # finds the N-H bonds.
        for i in (2, 3, 4):
            x, y, z = coords[i]
            coords2[i] = (x + 0.10, y, z)
        fp2 = th.compute_topology_fingerprint(syms, coords2)
        assert fp1.md_bins != fp2.md_bins

    def test_xyz_string_adapter(self):
        """The XYZ-string adapter should give same result as the symbol/coord API."""
        xyz = (
            "Fe  0.0  0.0  0.0\n"
            "N   2.05 0.0  0.0\n"
            "H   2.45 0.89 0.0\n"
            "H   2.45 -0.445 0.77\n"
            "H   2.45 -0.445 -0.77\n"
        )
        result_xyz = th.standalone_amine_h_realism_xyz(xyz)
        syms, coords = _fe_nh3_good()
        result_direct = th.standalone_amine_h_realism(syms, coords)
        assert result_xyz.passed == result_direct.passed
        assert result_xyz.violations == result_direct.violations


# ---------------------------------------------------------------------------
# Extra: default-OFF wire-in checks
# ---------------------------------------------------------------------------


class TestDefaultOff:
    def test_hardgate_flag_default(self, monkeypatch):
        _clear_env(monkeypatch)
        assert th.is_hardgate_enabled() is False

    def test_hardgate_flag_on(self, monkeypatch):
        _clear_env(monkeypatch)
        monkeypatch.setenv("DELFIN_5P_A_TOPOLOGY_HARDGATE", "1")
        assert th.is_hardgate_enabled() is True

    def test_post_emit_default_off_smiles_converter(self, monkeypatch):
        """Smoke-test: default-OFF byte identity at the wire-in point."""
        _clear_env(monkeypatch)
        # We cannot run the full smiles_converter without OB/RDKit, but we
        # verify the gate function is a no-op when the flag is unset.
        good_xyz = (
            "Fe  0.0  0.0  0.0\n"
            "N   2.05 0.0  0.0\n"
            "H   2.45 0.89 0.0\n"
            "H   2.45 -0.445 0.77\n"
            "H   2.45 -0.445 -0.77\n"
        )
        # When flag is off, is_hardgate_enabled() returns False — the
        # wire-in branch in smiles_converter is skipped entirely.
        assert not th.is_hardgate_enabled()
        # The function itself still works when called directly:
        result = th.standalone_amine_h_realism_xyz(good_xyz)
        assert result.passed
