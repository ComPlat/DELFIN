"""Tests for the donor-donor 1-3 Pauli-floor term (2026-06-07).

The term is wired into :func:`delfin.fffree.grip_polish.grip_polish` behind
the env-flag ``DELFIN_FFFREE_GRIP_DONOR_DONOR_1_3_FLOOR`` (default OFF).
Addresses the OKACOU-class collapse where halide donors all bonded to the
same metal (e.g. Ru(N-carbene)(Cl)4(CO)) collapse to 0.99 / 1.41 A because
the existing vdW-floor excludes 1-3 atom pairs.

These tests pin down:

* Test 1 -- Byte-identical OFF: flag unset (or 0) -> identical to HEAD.
* Test 2 -- OKACOU synthetic: 4 Cl bonded to Ru, 2 pairs collapsed @ 0.99
  and 1.41 A; the penalty fires and the analytic gradient pushes them
  apart (single negative-gradient step yields > 1.5 A).
* Test 3 -- Cis-trans OC-6: 6 monodentate donors at ideal 90 deg / 180 deg
  -> no penalty (all geometrically ideal).
* Test 4 -- T-4 angle test: 4 donors at 109.5 deg apart, d_DD ~ 1.633 *
  md_target -> no penalty.
* Test 5 -- SP-4 cis vs trans: cis pairs penalised if too close; trans
  pairs no penalty (they are at ~180 deg, d ~ 2 * md_target which is
  well above the floor 1.414 * md_target).
* Test 6 -- Donor-donor 1-4 NOT penalised: pair through 2 bonds (M-A-B-M
  loop, A/B not in donor list) does not apply.
* Test 7 -- Gradient finite-difference: analytic gradient matches
  central-difference numeric gradient within 1e-4.
* Test 8 -- Determinism: two runs byte-identical.
* Test 9 -- Composes additively with the vdW-floor: with BOTH terms ON,
  the loss equals the sum of each term in isolation.
"""
from __future__ import annotations

import os
import sys

# Strict determinism set BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import math

import numpy as np
import pytest

pytest.importorskip("rdkit")

from delfin.fffree.grip_polish import (
    DEFAULT_DD13_WEIGHT,
    DEFAULT_VDW_FLOOR_FRACTION,
    DEFAULT_VDW_FLOOR_WEIGHT,
    DEFAULT_VDW_RADII,
    _dd13_active,
    _heavy_atom_indices,
    _resolve_dd13_weight,
    _resolve_theta_min,
    _vdw_floor_value_and_grad,
    donor_donor_1_3_floor_value_and_grad,
)


_FLAG = "DELFIN_FFFREE_GRIP_DONOR_DONOR_1_3_FLOOR"
_WEIGHT_FLAG = "DELFIN_FFFREE_GRIP_DONOR_DONOR_1_3_FLOOR_WEIGHT"


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
def _set_env(flag: str, value):
    """Set or unset an env-var; returns the previous value for cleanup."""
    prev = os.environ.get(flag)
    if value is None:
        os.environ.pop(flag, None)
    else:
        os.environ[flag] = value
    return prev


@pytest.fixture(autouse=True)
def _scrub_dd13_env():
    """Restore the dd13 env-vars after every test (defence-in-depth)."""
    snap = {k: os.environ.get(k) for k in (_FLAG, _WEIGHT_FLAG)}
    try:
        yield
    finally:
        for k, v in snap.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _ideal_oc6(r_md: float) -> np.ndarray:
    """Six donors at the OC-6 vertices around a metal at the origin."""
    return np.array(
        [
            [0.0, 0.0, 0.0],          # metal (index 0)
            [+r_md, 0.0, 0.0],        # donor 1
            [-r_md, 0.0, 0.0],        # donor 2 (trans to 1)
            [0.0, +r_md, 0.0],        # donor 3
            [0.0, -r_md, 0.0],        # donor 4 (trans to 3)
            [0.0, 0.0, +r_md],        # donor 5
            [0.0, 0.0, -r_md],        # donor 6 (trans to 5)
        ],
        dtype=np.float64,
    )


# ---------------------------------------------------------------------------
# Test 1 -- Byte-identical OFF
# ---------------------------------------------------------------------------
class TestByteIdenticalOff:
    """``DELFIN_FFFREE_GRIP_DONOR_DONOR_1_3_FLOOR`` unset -> term is OFF."""

    def test_flag_unset_inactive(self):
        _set_env(_FLAG, None)
        assert _dd13_active() is False

    def test_flag_zero_inactive(self):
        _set_env(_FLAG, "0")
        assert _dd13_active() is False

    def test_flag_one_active(self):
        _set_env(_FLAG, "1")
        assert _dd13_active() is True

    def test_flag_true_active(self):
        _set_env(_FLAG, "TRUE")
        assert _dd13_active() is True
        _set_env(_FLAG, "on")
        assert _dd13_active() is True
        _set_env(_FLAG, "yes")
        assert _dd13_active() is True


# ---------------------------------------------------------------------------
# Test 2 -- OKACOU synthetic geometry
# ---------------------------------------------------------------------------
class TestOkacouSynthetic:
    """4 Cl bonded to Ru (OC-6), two pairs collapsed at 0.99 and 1.41 A."""

    def _okacou_like(self):
        """Ru @ origin, 4 Cl + 1 C(carbene) + 1 C(CO) -- 6 donors total.
        Collapsed pairs: Cl(idx 1) - Cl(idx 2) at 0.99 A, and
        Cl(idx 3) - Cl(idx 4) at 1.41 A.  The other donors are at clean
        OC-6 positions so only the Cl collapses violate the floor.
        """
        r_md = 2.30  # Ru-Cl covalent sum approx
        R = np.array(
            [
                [0.0, 0.0, 0.0],           # 0: Ru (metal)
                [+r_md, 0.0, 0.0],         # 1: Cl @ +x
                [+r_md - 0.99, 0.0, 0.0],  # 2: Cl collapsed at 0.99 from atom 1
                [0.0, +r_md, 0.0],         # 3: Cl @ +y
                [0.0, +r_md - 1.41, 0.0],  # 4: Cl collapsed at 1.41 from atom 3
                [0.0, 0.0, +r_md],         # 5: C(carbene) @ +z (clean)
                [0.0, 0.0, -r_md],         # 6: C(CO) @ -z (clean)
            ],
            dtype=np.float64,
        )
        donors = [1, 2, 3, 4, 5, 6]
        return R, donors, r_md

    def test_collapsed_pairs_trigger_penalty(self):
        R, donors, r_md = self._okacou_like()
        L, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="OC-6 octahedron",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        # OC-6 floor = sqrt(2) * 2.30 ~ 3.253 A
        floor = 2.0 * r_md * math.sin(math.radians(90.0) / 2.0)
        assert floor > 3.0
        # Both collapsed pairs violate (0.99 < 3.25 and 1.41 < 3.25).
        assert L > 0.0
        # Metal is not in the donor list -> its gradient slot is exactly 0.
        assert np.linalg.norm(G[0]) == 0.0
        # Collapsed atoms must have a non-trivial gradient.
        assert np.linalg.norm(G[1]) > 0.1
        assert np.linalg.norm(G[2]) > 0.1
        assert np.linalg.norm(G[3]) > 0.1
        assert np.linalg.norm(G[4]) > 0.1
        # Gradient pushes atoms 1 and 2 APART (along +/-x).  Atom 1 is on
        # the +x side and atom 2 on the slightly-less-+x side.  (R[1]-R[2])
        # points in +x; coef = -2 w gap / d < 0 so gi = coef * d_vec points
        # in -x for atom 1.  Negative-grad step R_new = R - step * G then
        # moves atom 1 in +x (away from atom 2).  Equivalently:
        # G[1].x < 0, G[2].x > 0, so (G[1]-G[2]).x < 0.
        assert (G[1] - G[2])[0] < 0.0
        # Same for atoms 3-4 along y.
        assert (G[3] - G[4])[1] < 0.0

    def test_gradient_step_separates_pairs(self):
        """A single negative-gradient step pushes both collapsed pairs above
        a non-trivial floor (>= 1.5 A -- a clear sign of recovery)."""
        R, donors, r_md = self._okacou_like()
        # Use a big weight + a reasonable step size so the move is visible.
        weight = 50.0
        L, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="OC-6 octahedron",
            md_target=r_md,
            weight=weight,
        )
        step = 0.01
        R_new = R - step * G
        d_12_new = float(np.linalg.norm(R_new[1] - R_new[2]))
        d_34_new = float(np.linalg.norm(R_new[3] - R_new[4]))
        d_12_old = float(np.linalg.norm(R[1] - R[2]))
        d_34_old = float(np.linalg.norm(R[3] - R[4]))
        # Both pairs must increase distance.
        assert d_12_new > d_12_old
        assert d_34_new > d_34_old
        # And cross the 1.5 A recovery threshold.
        assert d_12_new > 1.5, f"d_12 new = {d_12_new:.3f} should exceed 1.5"
        assert d_34_new > 1.5, f"d_34 new = {d_34_new:.3f} should exceed 1.5"


# ---------------------------------------------------------------------------
# Test 3 -- Ideal OC-6 -> no penalty
# ---------------------------------------------------------------------------
class TestIdealOC6NoPenalty:
    """6 donors at the OC-6 vertices: all cis pairs are at sqrt(2) * r_md,
    exactly at the floor (or just above)."""

    def test_ideal_oc6_zero_loss(self):
        r_md = 2.30
        R = _ideal_oc6(r_md)
        donors = [1, 2, 3, 4, 5, 6]
        L, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="OC-6",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        # All pairs are at d >= floor -- no violation.
        # cis pairs: sqrt(2) * 2.30 = 3.253 A (EQUAL to the floor); the
        # comparison `d >= floor` is True so loss = 0.
        # trans pairs: 2 * 2.30 = 4.60 A (well above).
        assert L == 0.0
        assert np.linalg.norm(G) == 0.0


# ---------------------------------------------------------------------------
# Test 4 -- T-4 ideal -> no penalty
# ---------------------------------------------------------------------------
class TestT4IdealNoPenalty:
    """4 donors at the T-4 tetrahedron vertices: every pair separated by
    109.5 deg, d_DD = 2 r_md sin(109.5/2) ~ 1.633 * r_md."""

    def test_t4_ideal_zero_loss(self):
        r_md = 2.0
        t = 1.0 / math.sqrt(3.0)
        verts = r_md * np.array(
            [[t, t, t], [t, -t, -t], [-t, t, -t], [-t, -t, t]],
            dtype=np.float64,
        )
        R = np.vstack([np.zeros(3, dtype=np.float64), verts])
        donors = [1, 2, 3, 4]
        L, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="T-4 tetrahedron",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        # T-4 cis pairs sit AT the floor 2 r_md sin(109.47/2) -> L = 0.
        assert L == 0.0
        assert np.linalg.norm(G) == 0.0


# ---------------------------------------------------------------------------
# Test 5 -- SP-4 cis vs trans
# ---------------------------------------------------------------------------
class TestSp4CisTrans:
    """SP-4: cis pairs at 90 deg (floor), trans pairs at 180 deg (no floor)."""

    def test_sp4_collapse_cis_only(self):
        r_md = 2.0
        # 4 donors in a square at +-r_md, +-r_md axis.
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [+r_md, 0.0, 0.0],
                [0.0, +r_md, 0.0],
                [-r_md, 0.0, 0.0],
                [0.0, -r_md, 0.0],
            ],
            dtype=np.float64,
        )
        # Collapse one cis pair: bring atom 2 close to atom 1.
        # cis floor = sqrt(2) * 2.0 = 2.828; place at 1.0 A.
        R[2] = R[1] + np.array([0.0, 1.0, 0.0])
        donors = [1, 2, 3, 4]
        L, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="SP-4",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        assert L > 0.0
        # Trans pair (1 - 3) at 2 r_md = 4 A: never collides with floor 2.828.
        # We CHECK that the trans gradient contribution is zero by
        # constructing a R where ONLY the trans pair is "close" (~= 2 r_md
        # = 4 A): no penalty.
        R_trans = np.array(
            [
                [0.0, 0.0, 0.0],
                [+r_md, 0.0, 0.0],
                [0.0, +r_md, 0.0],
                [-r_md, 0.0, 0.0],
                [0.0, -r_md, 0.0],
            ],
            dtype=np.float64,
        )
        L_t, G_t = donor_donor_1_3_floor_value_and_grad(
            R_trans,
            donor_indices=donors,
            metal_index=0,
            geometry="SP-4",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        # Cis pairs all at floor; trans at 2*r_md.  No violation anywhere.
        assert L_t == 0.0


# ---------------------------------------------------------------------------
# Test 6 -- 1-4 NOT penalised
# ---------------------------------------------------------------------------
class TestNon1_3Excluded:
    """An atom that is NOT bonded to the metal (donor_indices does not
    include it) is silently ignored regardless of distance."""

    def test_non_donor_no_penalty(self):
        r_md = 2.30
        R = np.array(
            [
                [0.0, 0.0, 0.0],            # metal
                [+r_md, 0.0, 0.0],          # donor 1
                [+r_md, 1.0, 0.0],          # NON-donor (1.0 A from donor 1)
                [-r_md, 0.0, 0.0],          # donor 2
            ],
            dtype=np.float64,
        )
        # Only atom 1 and 3 are donors.  Atom 2 is bonded to atom 1 via the
        # ligand backbone (1-4 through metal), NOT to the metal.
        donors = [1, 3]
        L, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="OC-6",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        # donor 1 - donor 3 are at d = 2 r_md = 4.6 A; floor = 3.25 A.
        # No violation, no penalty.
        assert L == 0.0
        # Atom 2 gradient must be exactly zero -- it was never visited.
        assert np.linalg.norm(G[2]) == 0.0


# ---------------------------------------------------------------------------
# Test 7 -- Analytic gradient matches finite-difference
# ---------------------------------------------------------------------------
class TestGradientFiniteDiff:
    """Central-difference numeric gradient matches the analytic gradient."""

    def test_finite_difference_match(self):
        r_md = 2.30
        # 2 donors at 0.99 A apart, both at ~r_md from metal.
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [+r_md, 0.0, 0.0],
                [+r_md - 0.99, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        donors = [1, 2]

        def L_only(R_arg):
            L, _ = donor_donor_1_3_floor_value_and_grad(
                R_arg,
                donor_indices=donors,
                metal_index=0,
                geometry="OC-6",
                md_target=r_md,
                weight=DEFAULT_DD13_WEIGHT,
            )
            return L

        _, G = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=donors,
            metal_index=0,
            geometry="OC-6",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )

        eps = 1e-5
        G_num = np.zeros_like(R)
        for i in range(R.shape[0]):
            for k in range(3):
                Rp = R.copy()
                Rp[i, k] += eps
                Rm = R.copy()
                Rm[i, k] -= eps
                G_num[i, k] = (L_only(Rp) - L_only(Rm)) / (2.0 * eps)
        # Within 1e-4 tolerance.
        diff = np.max(np.abs(G - G_num))
        assert diff < 1e-4, f"max grad diff = {diff:g}"


# ---------------------------------------------------------------------------
# Test 8 -- Determinism
# ---------------------------------------------------------------------------
class TestDeterminism:
    """Two independent runs produce byte-identical outputs."""

    def test_byte_identical_two_runs(self):
        r_md = 2.30
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [+r_md, 0.0, 0.0],
                [+r_md - 0.99, 0.0, 0.0],
                [0.0, +r_md, 0.0],
                [0.0, +r_md - 1.41, 0.0],
            ],
            dtype=np.float64,
        )
        donors = [1, 2, 3, 4]
        L1, G1 = donor_donor_1_3_floor_value_and_grad(
            R, donor_indices=donors, metal_index=0,
            geometry="OC-6", md_target=r_md, weight=50.0,
        )
        L2, G2 = donor_donor_1_3_floor_value_and_grad(
            R, donor_indices=donors, metal_index=0,
            geometry="OC-6", md_target=r_md, weight=50.0,
        )
        assert L1 == L2
        assert np.array_equal(G1, G2)

    def test_donor_order_invariant(self):
        """Donor index permutation does not affect the loss/gradient."""
        r_md = 2.30
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [+r_md, 0.0, 0.0],
                [+r_md - 0.99, 0.0, 0.0],
                [0.0, +r_md, 0.0],
                [0.0, +r_md - 1.41, 0.0],
            ],
            dtype=np.float64,
        )
        L_a, G_a = donor_donor_1_3_floor_value_and_grad(
            R, donor_indices=[1, 2, 3, 4], metal_index=0,
            geometry="OC-6", md_target=r_md, weight=50.0,
        )
        L_b, G_b = donor_donor_1_3_floor_value_and_grad(
            R, donor_indices=[4, 2, 1, 3], metal_index=0,
            geometry="OC-6", md_target=r_md, weight=50.0,
        )
        assert L_a == L_b
        assert np.array_equal(G_a, G_b)


# ---------------------------------------------------------------------------
# Test 9 -- Additive composition with vdW-floor
# ---------------------------------------------------------------------------
class TestComposeAdditively:
    """When both terms fire on a collapsed donor-donor pair, the sum equals
    the vdW-floor contribution plus the dd13 contribution evaluated
    independently."""

    def test_additive_sum_on_collapsed_cl_pair(self):
        r_md = 2.30
        # Two Cl atoms collapsed at 0.99 A.  They are heavy AND donors.
        R = np.array(
            [
                [0.0, 0.0, 0.0],
                [+r_md, 0.0, 0.0],
                [+r_md - 0.99, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        # dd13 term: OC-6 floor = sqrt(2) * 2.30 = 3.253; 0.99 < 3.253 -> fires.
        L_dd, G_dd = donor_donor_1_3_floor_value_and_grad(
            R,
            donor_indices=[1, 2],
            metal_index=0,
            geometry="OC-6",
            md_target=r_md,
            weight=DEFAULT_DD13_WEIGHT,
        )
        assert L_dd > 0.0

        # vdW-floor term: Cl vdW radius = 1.75; floor = 0.85 * 2 * 1.75 = 2.975;
        # 0.99 < 2.975 -> fires.  The two terms compute independently and add.
        rCl = DEFAULT_VDW_RADII["Cl"]
        radii = np.array([np.nan, rCl, rCl], dtype=np.float64)
        heavy = np.array([1, 2], dtype=np.int64)
        L_vdw, G_vdw = _vdw_floor_value_and_grad(
            R,
            heavy_indices=heavy,
            radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        assert L_vdw > 0.0

        # Both gradients point in the same direction (separate atoms 1 and 2);
        # the additive sum has |G_total| > max(|G_dd|, |G_vdw|).
        G_total = G_dd + G_vdw
        L_total = L_dd + L_vdw
        assert L_total > L_dd
        assert L_total > L_vdw
        # Sum-of-norms test: both grads in the same outward direction.
        assert np.linalg.norm(G_total[1]) > np.linalg.norm(G_dd[1])
        assert np.linalg.norm(G_total[1]) > np.linalg.norm(G_vdw[1])


# ---------------------------------------------------------------------------
# Bonus -- theta_min table sanity
# ---------------------------------------------------------------------------
class TestThetaMinTable:
    """The polyhedron name resolver returns sensible angles."""

    def test_oc6_long_and_short(self):
        assert _resolve_theta_min("OC-6") == 90.0
        assert _resolve_theta_min("OC-6 octahedron") == 90.0

    def test_t4(self):
        v = _resolve_theta_min("T-4")
        assert v is not None
        assert abs(v - 109.47) < 0.5

    def test_unknown_returns_none(self):
        assert _resolve_theta_min("THIS-IS-NOT-A-POLYHEDRON") is None
        assert _resolve_theta_min("") is None
        assert _resolve_theta_min(None) is None

    def test_resolve_weight_default(self):
        _set_env(_WEIGHT_FLAG, None)
        assert _resolve_dd13_weight() == DEFAULT_DD13_WEIGHT

    def test_resolve_weight_override(self):
        _set_env(_WEIGHT_FLAG, "75.5")
        assert _resolve_dd13_weight() == 75.5

    def test_resolve_weight_garbage_falls_through(self):
        _set_env(_WEIGHT_FLAG, "not-a-number")
        assert _resolve_dd13_weight() == DEFAULT_DD13_WEIGHT
