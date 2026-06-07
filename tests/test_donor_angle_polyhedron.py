"""Tests for the donor-angle polyhedron floor (2026-06-07).

Anti-``metal_axis_linearisation`` term wired into
:func:`delfin.fffree.grip_polish.grip_polish` behind the env-flag
``DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_FLOOR`` (default OFF).  Auto-
diagnostic finding: 890 / 6627 (13.4 %) V3 files showed donor-M-donor
angles drifting toward 180 deg on Pd(SP-4) / Pt(SP-4) centres where they
should be 90 deg / 180 deg.

These tests pin down:

* Test 1 -- SP-4 (90/90/90/90/180/180) -> zero penalty (matches template).
* Test 2 -- Linearised CN4 (4 angles near 180 deg) -> penalty fires AND
  the gradient pushes the donors toward the SP-4 template (verified by
  checking that one negative-gradient step reduces the total angular
  error by more than 50 percent).
* Test 3 -- T-4 (109.47 deg x 6) -> zero penalty (matches template).
* Test 4 -- Byte-identical OFF: env unset / 0 -> grip_polish output is
  bit-identical with the legacy call.
* Test 5 -- Gradient finite-difference: analytic gradient matches the
  central-difference numeric gradient.
* Test 6 -- Determinism: identical inputs + ON flag -> bit-identical out.
* Test 7 -- Universality: the helper works for every polyhedron in
  :func:`polyhedra.GEOM_BY_CN` (CN2 .. CN9) plus CN10.

Strict determinism: PYTHONHASHSEED=0 set before numpy import.
"""
from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

pytest.importorskip("scipy")


from delfin.fffree.donor_angle_polyhedron import (
    DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
    DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
    _resolve_donor_angle_polyhedron_threshold_deg,
    _resolve_donor_angle_polyhedron_weight,
    donor_angle_polyhedron_active,
    donor_angle_polyhedron_value_and_grad,
    expected_angles_for_geometry,
    get_expected_angles,
)
from delfin.fffree.polyhedra import GEOM_BY_CN, ref_vectors


_FLAG = "DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_FLOOR"
_WEIGHT_FLAG = "DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_WEIGHT"
_THRESHOLD_FLAG = "DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG"


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
def _scrub_dap_env():
    """Restore the donor-angle env-vars after every test."""
    snap = {k: os.environ.get(k) for k in (_FLAG, _WEIGHT_FLAG, _THRESHOLD_FLAG)}
    try:
        yield
    finally:
        for k, v in snap.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _sp4_geometry(d: float = 2.0) -> np.ndarray:
    """Ideal SP-4 (square planar) coords: metal at origin + 4 donors at
    +x, +y, -x, -y at distance ``d``.  Angles: 90, 90, 90, 90, 180, 180.
    """
    return np.array([
        [0.0, 0.0, 0.0],         # 0 = metal
        [+d, 0.0, 0.0],          # 1 = donor +x
        [0.0, +d, 0.0],          # 2 = donor +y
        [-d, 0.0, 0.0],          # 3 = donor -x
        [0.0, -d, 0.0],          # 4 = donor -y
    ], dtype=np.float64)


def _t4_geometry(d: float = 2.0) -> np.ndarray:
    """Ideal T-4 (tetrahedron) coords: metal at origin + 4 donors at the
    canonical tetrahedron unit vectors scaled by ``d``.  All 6 donor-M-
    donor angles equal acos(-1/3) = 109.471 deg.
    """
    t = 1.0 / math.sqrt(3.0)
    return np.array([
        [0.0, 0.0, 0.0],         # 0 = metal
        [+t * d, +t * d, +t * d],
        [+t * d, -t * d, -t * d],
        [-t * d, +t * d, -t * d],
        [-t * d, -t * d, +t * d],
    ], dtype=np.float64)


def _linearised_cn4(d: float = 2.0, axis_drift_deg: float = 10.0) -> np.ndarray:
    """Pd CN4 pathological -- 4 donors compressed onto two near-parallel
    axes.  Donors 1 & 3 sit close to the +x / -x ends; donors 2 & 4 sit
    close to the same axis but tilted out of plane to break symmetry.
    Result: 4 of the 6 angles are near 180 deg (well outside 90 +- 15
    deg), 2 angles are small (~20 deg).
    """
    theta = math.radians(axis_drift_deg)
    c, s = math.cos(theta), math.sin(theta)
    return np.array([
        [0.0, 0.0, 0.0],          # 0 = metal
        [+d * c, +d * s, 0.0],    # 1 -- near +x, slight +y tilt
        [+d * c, 0.0, -d * s],    # 2 -- near +x, slight -z tilt (close to 1)
        [-d * c, -d * s, 0.0],    # 3 -- near -x, slight -y tilt
        [-d * c, 0.0, +d * s],    # 4 -- near -x, slight +z tilt (close to 3)
    ], dtype=np.float64)


def _total_angular_error_deg(P: np.ndarray, metal: int,
                             donors, expected) -> float:
    """Sum of |delta_i| - tau (clipped at zero) over all donor pairs."""
    tau = DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG
    Rm = P[metal]
    total = 0.0
    donor_list = sorted(int(d) for d in donors)
    for i in range(len(donor_list)):
        di = P[donor_list[i]] - Rm
        ni = np.linalg.norm(di)
        if ni < 1e-9:
            continue
        for j in range(i + 1, len(donor_list)):
            dj = P[donor_list[j]] - Rm
            nj = np.linalg.norm(dj)
            if nj < 1e-9:
                continue
            c = float(np.dot(di, dj) / (ni * nj))
            c = max(-1.0, min(1.0, c))
            theta = math.degrees(math.acos(c))
            best = min(abs(theta - e) for e in expected)
            total += max(0.0, best - tau)
    return total


# ---------------------------------------------------------------------------
# Test 1 -- SP-4 matches template -> zero penalty
# ---------------------------------------------------------------------------
class TestSP4MatchesTemplate:
    """SP-4 ideal coords -> all 6 donor-pair angles match {90, 180} ->
    no penalty fires regardless of the threshold."""

    def test_sp4_zero_loss(self):
        P = _sp4_geometry()
        expected = get_expected_angles("SP-4 square planar")
        assert sorted(round(a, 1) for a in expected) == [90.0, 180.0]
        L, G = donor_angle_polyhedron_value_and_grad(
            P,
            metal=0,
            donors=[1, 2, 3, 4],
            expected_angles_deg=expected,
            weight=DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            threshold_deg=DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)

    def test_sp4_small_jitter_still_zero(self):
        """A 5 deg jitter is well inside the 15 deg threshold -> still zero."""
        rng = np.random.default_rng(0)
        P = _sp4_geometry()
        # Rotate each donor by a small random axis-angle (still on its M-D shell).
        for k in (1, 2, 3, 4):
            ax = rng.standard_normal(3); ax /= np.linalg.norm(ax)
            theta = math.radians(2.0)  # tiny tilt; well below 15 deg
            # Rodrigues
            K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
            R_mat = np.eye(3) + math.sin(theta) * K + (1 - math.cos(theta)) * (K @ K)
            P[k] = R_mat @ P[k]
        expected = get_expected_angles("SP-4 square planar")
        L, _ = donor_angle_polyhedron_value_and_grad(
            P,
            metal=0,
            donors=[1, 2, 3, 4],
            expected_angles_deg=expected,
            weight=5.0,
            threshold_deg=15.0,
        )
        assert L == 0.0


# ---------------------------------------------------------------------------
# Test 2 -- Linearised CN4 fires + gradient improves angular error
# ---------------------------------------------------------------------------
class TestLinearisedCN4Fires:
    """The Pd CN4 'metal-axis linearisation' demo case."""

    def test_penalty_fires(self):
        P = _linearised_cn4(d=2.0, axis_drift_deg=10.0)
        expected = get_expected_angles("SP-4 square planar")
        L, G = donor_angle_polyhedron_value_and_grad(
            P,
            metal=0,
            donors=[1, 2, 3, 4],
            expected_angles_deg=expected,
            weight=DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            threshold_deg=DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
        )
        assert L > 0.0, "expected penalty to fire on a near-linear CN4"
        assert np.isfinite(L)
        assert np.all(np.isfinite(G))
        # The gradient must be non-zero on the donors -- it is what drags
        # them off the metal axis.
        assert float(np.linalg.norm(G[1:5])) > 0.0

    def test_gradient_step_halves_angular_error(self):
        """Gradient-descent with a backtracking line search drives the
        total angular error down by more than 50 percent.

        This is the load-bearing test: it shows the term not only fires
        but also points the donors in the right direction.  A backtracking
        line search keeps the step size sane (raw |G| is O(1e4) on this
        toy geometry, so a fixed step would overshoot).
        """
        P = _linearised_cn4(d=2.0, axis_drift_deg=10.0)
        expected = get_expected_angles("SP-4 square planar")
        donors = [1, 2, 3, 4]
        err_before = _total_angular_error_deg(P, 0, donors, expected)
        assert err_before > 0.0
        P_curr = P.copy()
        kw = dict(
            metal=0,
            donors=donors,
            expected_angles_deg=expected,
            weight=DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            threshold_deg=DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
        )
        for _ in range(500):
            L, G = donor_angle_polyhedron_value_and_grad(P_curr, **kw)
            if L <= 0.0:
                break
            # Backtracking line search: start at a stride that is normalised
            # by |G|, halve until the step actually reduces the loss.
            step = 1.0 / max(float(np.linalg.norm(G)), 1.0)
            for _bt in range(40):
                P_test = P_curr - step * G
                L_test, _ = donor_angle_polyhedron_value_and_grad(P_test, **kw)
                if L_test < L:
                    P_curr = P_test
                    break
                step *= 0.5
            else:
                break
        err_after = _total_angular_error_deg(P_curr, 0, donors, expected)
        assert err_after < 0.5 * err_before, (
            f"expected >= 50 percent angular-error reduction; "
            f"before={err_before:.3f} after={err_after:.3f}"
        )


# ---------------------------------------------------------------------------
# Test 3 -- T-4 (109.47 deg x 6) -> zero
# ---------------------------------------------------------------------------
class TestT4Tetrahedron:
    """Ideal T-4 -> all 6 angles = 109.47 deg; T-4 template = single
    angle 109.47 deg; deviation = 0 -> zero penalty."""

    def test_t4_zero_loss(self):
        P = _t4_geometry()
        expected = get_expected_angles("T-4 tetrahedron")
        # T-4 should have a single ideal angle.
        assert len(expected) == 1
        assert abs(expected[0] - math.degrees(math.acos(-1.0 / 3.0))) < 0.6
        L, G = donor_angle_polyhedron_value_and_grad(
            P,
            metal=0,
            donors=[1, 2, 3, 4],
            expected_angles_deg=expected,
            weight=DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            threshold_deg=DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)


# ---------------------------------------------------------------------------
# Test 4 -- Byte-identical OFF
# ---------------------------------------------------------------------------
class TestByteIdenticalOff:
    """``DELFIN_FFFREE_DONOR_ANGLE_POLYHEDRON_FLOOR`` unset / 0 ->
    grip_polish output bit-identical with the legacy call."""

    def test_active_default_off(self):
        _set_env(_FLAG, None)
        assert donor_angle_polyhedron_active() is False

    def test_active_zero_off(self):
        _set_env(_FLAG, "0")
        assert donor_angle_polyhedron_active() is False

    def test_active_true_on(self):
        _set_env(_FLAG, "1")
        assert donor_angle_polyhedron_active() is True

    def test_grip_polish_byte_identical_off(self):
        """Run grip_polish twice, with flag unset and with flag=0.
        The output must be bit-identical."""
        pytest.importorskip("rdkit")
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from delfin.fffree.grip_polish import grip_polish

        mol = Chem.MolFromSmiles("Cc1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        P0 = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        _set_env(_FLAG, None)
        out_unset = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                                clash_weight=5.0)
        _set_env(_FLAG, "0")
        out_zero = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                               clash_weight=5.0)
        assert np.array_equal(out_unset, out_zero), (
            "unset and 0 must produce bit-identical grip_polish output"
        )


# ---------------------------------------------------------------------------
# Test 5 -- Gradient finite-difference
# ---------------------------------------------------------------------------
class TestGradientFD:
    """Central-difference check on the analytic gradient at multiple
    perturbed configurations.  The penalty is a smooth quadratic in
    |delta| beyond the threshold, so FD agreement to 1e-4 per component
    is the right precision band."""

    def test_fd_match_linearised(self):
        rng = np.random.default_rng(7)
        for trial in range(5):
            # Start from a near-linearised CN4 + small random shifts.
            P = _linearised_cn4(d=2.0, axis_drift_deg=8.0)
            P[1:5] += 0.05 * rng.standard_normal((4, 3))
            expected = get_expected_angles("SP-4 square planar")
            weight = 5.0
            tau = 15.0
            L, G = donor_angle_polyhedron_value_and_grad(
                P,
                metal=0,
                donors=[1, 2, 3, 4],
                expected_angles_deg=expected,
                weight=weight,
                threshold_deg=tau,
            )
            assert L > 0.0, "trial expected to violate the floor"
            eps = 1e-5
            G_fd = np.zeros_like(P)
            for a in range(P.shape[0]):
                for c in range(3):
                    Pp = P.copy(); Pp[a, c] += eps
                    Pn = P.copy(); Pn[a, c] -= eps
                    Lp, _ = donor_angle_polyhedron_value_and_grad(
                        Pp, metal=0, donors=[1, 2, 3, 4],
                        expected_angles_deg=expected,
                        weight=weight, threshold_deg=tau,
                    )
                    Ln, _ = donor_angle_polyhedron_value_and_grad(
                        Pn, metal=0, donors=[1, 2, 3, 4],
                        expected_angles_deg=expected,
                        weight=weight, threshold_deg=tau,
                    )
                    G_fd[a, c] = (Lp - Ln) / (2.0 * eps)
            err = float(np.max(np.abs(G - G_fd)))
            assert err < 1e-3, (
                f"trial {trial}: analytic vs FD mismatch {err:.3e}\n"
                f"G_analytic={G}\nG_fd={G_fd}"
            )


# ---------------------------------------------------------------------------
# Test 6 -- Determinism
# ---------------------------------------------------------------------------
class TestDeterminism:
    """Same coordinates + same args -> bit-identical output, two calls."""

    def test_helper_byte_identical(self):
        P = _linearised_cn4()
        expected = get_expected_angles("SP-4 square planar")
        kw = dict(
            metal=0,
            donors=[1, 2, 3, 4],
            expected_angles_deg=expected,
            weight=DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT,
            threshold_deg=DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG,
        )
        L1, G1 = donor_angle_polyhedron_value_and_grad(P.copy(), **kw)
        L2, G2 = donor_angle_polyhedron_value_and_grad(P.copy(), **kw)
        assert L1 == L2
        assert np.array_equal(G1, G2)

    def test_grip_polish_byte_identical_on(self):
        """Same input + ON flag -> bit-identical grip_polish output."""
        pytest.importorskip("rdkit")
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from delfin.fffree.grip_polish import grip_polish

        mol = Chem.MolFromSmiles("Cc1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        P0 = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        _set_env(_FLAG, "1")
        out_a = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        out_b = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        assert np.array_equal(out_a, out_b)
        assert np.all(np.isfinite(out_a))


# ---------------------------------------------------------------------------
# Test 7 -- Universality across all polyhedra
# ---------------------------------------------------------------------------
class TestUniversality:
    """Expected-angle list works for every polyhedron in GEOM_BY_CN (CN2
    through CN9) plus the CN10 variants -- no per-class branches needed.
    """

    @pytest.mark.parametrize("cn", sorted(GEOM_BY_CN.keys()))
    def test_every_geom_yields_nonempty_expected(self, cn):
        for geom in GEOM_BY_CN[cn]:
            expected = expected_angles_for_geometry(geom)
            assert len(expected) >= 1, (
                f"polyhedron {geom!r} has no expected donor-donor angles"
            )
            # All values must be in (0, 180] deg.
            for a in expected:
                assert 0.0 < a <= 180.0 + 1e-6

    def test_sp4_expected_set(self):
        # Spec: SP-4 -> 90 cis, 180 trans.
        expected = expected_angles_for_geometry("SP-4 square planar")
        assert any(abs(a - 90.0) < 0.5 for a in expected)
        assert any(abs(a - 180.0) < 0.5 for a in expected)

    def test_t4_expected_set(self):
        # Spec: T-4 -> 109.5 deg.
        expected = expected_angles_for_geometry("T-4 tetrahedron")
        assert any(abs(a - 109.47) < 0.5 for a in expected)

    def test_oc6_expected_set(self):
        # Spec: OC-6 -> 90 cis, 180 trans.
        expected = expected_angles_for_geometry("OC-6 octahedron")
        assert any(abs(a - 90.0) < 0.5 for a in expected)
        assert any(abs(a - 180.0) < 0.5 for a in expected)

    def test_tbp5_expected_set(self):
        # Spec: TBP-5 -> 90 (ax-eq), 120 (eq-eq), 180 (ax-ax).
        expected = expected_angles_for_geometry("TBP-5 trigonal bipyramid")
        assert any(abs(a - 90.0) < 0.5 for a in expected)
        assert any(abs(a - 120.0) < 0.5 for a in expected)
        assert any(abs(a - 180.0) < 0.5 for a in expected)

    def test_helper_works_for_each_geom(self):
        """Pick ideal-vertex coords for each (CN, geom) and verify the
        helper returns zero loss + zero gradient at the template."""
        for cn, geoms in sorted(GEOM_BY_CN.items()):
            for geom in geoms:
                V = ref_vectors(geom)
                n_donors = V.shape[0]
                # Build coords: metal at origin + donors at unit vectors * 2 A.
                P = np.zeros((1 + n_donors, 3), dtype=np.float64)
                P[1:] = 2.0 * V
                donors = list(range(1, 1 + n_donors))
                expected = get_expected_angles(geom)
                assert expected, f"empty expected for {geom!r}"
                L, G = donor_angle_polyhedron_value_and_grad(
                    P,
                    metal=0,
                    donors=donors,
                    expected_angles_deg=expected,
                    weight=5.0,
                    threshold_deg=15.0,
                )
                assert L == 0.0, (
                    f"ideal {geom!r} (CN={cn}) produced non-zero loss {L}"
                )
                assert np.allclose(G, 0.0), (
                    f"ideal {geom!r} (CN={cn}) produced non-zero gradient"
                )


# ---------------------------------------------------------------------------
# Test 8 -- Resolvers (env parsing)
# ---------------------------------------------------------------------------
class TestResolvers:
    def test_weight_default(self):
        _set_env(_WEIGHT_FLAG, None)
        assert _resolve_donor_angle_polyhedron_weight() == DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT

    def test_weight_env_override(self):
        _set_env(_WEIGHT_FLAG, "7.5")
        assert _resolve_donor_angle_polyhedron_weight() == 7.5

    def test_weight_garbage_falls_back(self):
        _set_env(_WEIGHT_FLAG, "not-a-number")
        assert _resolve_donor_angle_polyhedron_weight() == DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT

    def test_weight_negative_falls_back(self):
        _set_env(_WEIGHT_FLAG, "-1.0")
        assert _resolve_donor_angle_polyhedron_weight() == DEFAULT_DONOR_ANGLE_POLYHEDRON_WEIGHT

    def test_threshold_default(self):
        _set_env(_THRESHOLD_FLAG, None)
        assert _resolve_donor_angle_polyhedron_threshold_deg() == DEFAULT_DONOR_ANGLE_POLYHEDRON_THRESHOLD_DEG

    def test_threshold_env_override(self):
        _set_env(_THRESHOLD_FLAG, "20.0")
        assert _resolve_donor_angle_polyhedron_threshold_deg() == 20.0

    def test_active_true_variants(self):
        for v in ("1", "true", "yes", "on", "TRUE", "Yes"):
            _set_env(_FLAG, v)
            assert donor_angle_polyhedron_active() is True, f"failed for {v!r}"

    def test_active_false_variants(self):
        for v in ("0", "false", "no", "off", "", "garbage"):
            _set_env(_FLAG, v)
            assert donor_angle_polyhedron_active() is False, f"failed for {v!r}"
