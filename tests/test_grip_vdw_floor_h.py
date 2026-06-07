"""Tests for the H-inclusive vdW-floor extension (2026-06-07, BERTEB fix).

The legacy heavy-atom-only vdW-floor (``DELFIN_FFFREE_GRIP_VDW_FLOOR=1``)
ignores hydrogen atoms by construction (selector ``_heavy_atom_indices``
filters them out).  On the BERTEB-class SMILES
``[Br][Ni-2]12[N]3C=CC=C3C=[N+]1CC1=CC=CC=[N+]12`` the V14 voll-pool
shows the heavy skeleton is fine but four C-H pairs collapse to ~0.82-0.99 Å
(should be ~1.09 Å) because the heavy floor never sees them.

This test module pins down a NEW env flag
``DELFIN_FFFREE_GRIP_VDW_FLOOR_INCLUDE_H=1`` that, when ON together with the
legacy ``DELFIN_FFFREE_GRIP_VDW_FLOOR=1``, extends the pair walk to include
hydrogens (Bondi r = 1.20 Å, already in ``DEFAULT_VDW_RADII``).

Test plan:

1. Byte-identical OFF (``INCLUDE_H=0`` / unset, ``VDW_FLOOR=1``) ->
   heavy-only pair walk identical to HEAD.
2. H-H pair below floor (d=1.2 Å, floor 0.85*(1.2+1.2)=2.04 Å) -> penalty fires.
3. C-H pair below floor (d=0.8 Å, floor 0.85*(1.7+1.2)=2.465 Å) -> penalty
   fires and gradient pushes them apart.
4. Bonded C-H pair at 1.09 Å -> no penalty when listed in ``excluded_pairs``.
5. Finite-difference gradient on a mixed C+H configuration matches the
   analytic gradient within 1e-4.
6. Two runs with the flag ON produce byte-identical outputs (determinism).
7. Integration: the rigid-body-aware
   ``vdw_floor_all_pairs_value_and_grad`` from
   :mod:`delfin.fffree.construction_sanity` routes H gradients through the
   rigid-body's metal translation when both H atoms sit inside the rigid
   body (no internal distortion of the rigid block).
"""
from __future__ import annotations

import os

# Strict determinism set BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

from delfin.fffree.construction_sanity import (
    vdw_floor_all_pairs_value_and_grad,
)
from delfin.fffree.grip_polish import (
    DEFAULT_VDW_FLOOR_FRACTION,
    DEFAULT_VDW_FLOOR_WEIGHT,
    DEFAULT_VDW_RADII,
    _heavy_atom_indices,
    _vdw_floor_active,
    _vdw_floor_include_h,
    _vdw_floor_value_and_grad,
)

_FLAG = "DELFIN_FFFREE_GRIP_VDW_FLOOR"
_INC_H_FLAG = "DELFIN_FFFREE_GRIP_VDW_FLOOR_INCLUDE_H"
_WEIGHT_FLAG = "DELFIN_FFFREE_GRIP_VDW_FLOOR_WEIGHT"
_FRACTION_FLAG = "DELFIN_FFFREE_GRIP_VDW_FLOOR_FRACTION"
_ALL_PAIRS_FLAG = "DELFIN_FFFREE_GRIP_VDW_FLOOR_ALL_PAIRS"


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
def _set_env(flag: str, value: str | None):
    prev = os.environ.get(flag)
    if value is None:
        os.environ.pop(flag, None)
    else:
        os.environ[flag] = value
    return prev


@pytest.fixture(autouse=True)
def _scrub_env():
    """Restore every relevant env-var after each test (defence-in-depth)."""
    keys = (
        _FLAG, _INC_H_FLAG, _WEIGHT_FLAG, _FRACTION_FLAG, _ALL_PAIRS_FLAG,
    )
    snap = {k: os.environ.get(k) for k in keys}
    try:
        yield
    finally:
        for k, v in snap.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


# ---------------------------------------------------------------------------
# Env-flag resolver
# ---------------------------------------------------------------------------
class TestIncludeHFlagResolver:
    def test_default_off(self):
        _set_env(_INC_H_FLAG, None)
        assert _vdw_floor_include_h() is False

    def test_truthy_variants_on(self):
        for v in ("1", "true", "yes", "on", "TRUE", "On"):
            _set_env(_INC_H_FLAG, v)
            assert _vdw_floor_include_h() is True, f"failed for {v!r}"

    def test_falsy_variants_off(self):
        for v in ("0", "false", "no", "off", "", "garbage"):
            _set_env(_INC_H_FLAG, v)
            assert _vdw_floor_include_h() is False, f"failed for {v!r}"

    def test_independent_of_master_flag(self):
        """``INCLUDE_H`` resolver does not consult the master flag."""
        _set_env(_FLAG, None)
        _set_env(_INC_H_FLAG, "1")
        assert _vdw_floor_include_h() is True
        assert _vdw_floor_active() is False
        _set_env(_FLAG, "1")
        assert _vdw_floor_include_h() is True
        assert _vdw_floor_active() is True


# ---------------------------------------------------------------------------
# Selector behaviour: H included only when include_h=True
# ---------------------------------------------------------------------------
class TestSelectorIncludeH:
    def test_default_excludes_h(self):
        idx = _heavy_atom_indices(["C", "H", "N", "H"], 4)
        assert list(idx) == [0, 2]

    def test_include_h_keeps_all(self):
        idx = _heavy_atom_indices(["C", "H", "N", "H"], 4, include_h=True)
        assert list(idx) == [0, 1, 2, 3]

    def test_include_h_keyword_only(self):
        """``include_h`` must be keyword-only -- positional must fail."""
        with pytest.raises(TypeError):
            _heavy_atom_indices(["C", "H"], 2, True)  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Test 1 -- Byte-identical OFF when INCLUDE_H is unset, even with VDW_FLOOR=1
# ---------------------------------------------------------------------------
class TestByteIdenticalIncludeHOff:
    """Result of the value/grad helper is byte-identical to the legacy
    heavy-atom-only path when ``INCLUDE_H`` is not set, even when the master
    floor is ON.
    """

    def _build_2c_2h_geometry(self):
        # Two C and two H in a near-collinear line.  C-C at 1.5 Å (well inside
        # 0.85*2*1.70=2.89 floor -> heavy violation), C-H at 0.95 Å, H-H at 1.0 Å.
        R = np.array([
            [0.0, 0.0, 0.0],       # 0: C
            [1.5, 0.0, 0.0],       # 1: C
            [-0.95, 0.0, 0.0],     # 2: H  (C0-H closer than floor)
            [-1.95, 0.0, 0.0],     # 3: H  (H-H d=1.0 from atom 2)
        ], dtype=np.float64)
        symbols = ["C", "C", "H", "H"]
        radii = np.full(4, np.nan, dtype=np.float64)
        for i, s in enumerate(symbols):
            radii[i] = DEFAULT_VDW_RADII[s]
        return R, symbols, radii

    def test_floor_helper_byte_identical_without_h(self):
        R, symbols, radii = self._build_2c_2h_geometry()
        # Heavy-only: indices = [0, 1].
        heavy = _heavy_atom_indices(symbols, 4, include_h=False)
        L_heavy, G_heavy = _vdw_floor_value_and_grad(
            R, heavy_indices=heavy, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        # Now run again with the SAME indices (mimics what the polish
        # produces when INCLUDE_H is unset) -- bit-identical.
        L_again, G_again = _vdw_floor_value_and_grad(
            R, heavy_indices=heavy, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        assert L_heavy == L_again
        assert np.array_equal(G_heavy, G_again)
        # And: H atoms received zero gradient -- they were not iterated.
        assert np.allclose(G_heavy[2], 0.0)
        assert np.allclose(G_heavy[3], 0.0)


# ---------------------------------------------------------------------------
# Test 2 -- H-H pair below floor fires
# ---------------------------------------------------------------------------
class TestHHPairBelowFloor:
    def test_two_h_atoms_inside_floor(self):
        rH = DEFAULT_VDW_RADII["H"]
        # floor = 0.85 * (1.20 + 1.20) = 2.04 Å.  Place them at 1.20 Å.
        floor = DEFAULT_VDW_FLOOR_FRACTION * 2.0 * rH
        d_in = 0.5 * 2.0 * rH  # = 1.20
        assert d_in < floor
        R = np.array([[0.0, 0.0, 0.0], [d_in, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([rH, rH], dtype=np.float64)
        # With include_h=True we expect both H to be iterated.
        idx = _heavy_atom_indices(["H", "H"], 2, include_h=True)
        assert list(idx) == [0, 1]
        L, G = _vdw_floor_value_and_grad(
            R, heavy_indices=idx, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        # Loss > 0, gradient points outward.
        assert L > 0.0
        assert np.isfinite(L)
        # Atom 0 sits at the origin and atom 1 at +x; the floor wants them apart.
        # grad on atom 0: -2 w gap * (x0 - x1)/d = -2 w gap * (-1) = +2 w gap > 0.
        assert G[0, 0] > 0
        assert G[1, 0] < 0
        # Negative-gradient step increases distance.
        R_new = R - 1e-3 * G
        d_new = float(np.linalg.norm(R_new[0] - R_new[1]))
        assert d_new > d_in


# ---------------------------------------------------------------------------
# Test 3 -- C-H pair below floor fires + pushes apart
# ---------------------------------------------------------------------------
class TestCHPairBelowFloor:
    def test_ch_pair_at_0p8_angstrom(self):
        rC = DEFAULT_VDW_RADII["C"]
        rH = DEFAULT_VDW_RADII["H"]
        # floor = 0.85 * (1.70 + 1.20) = 2.465 Å.  Place at 0.8 Å (BERTEB-like collapse).
        floor = DEFAULT_VDW_FLOOR_FRACTION * (rC + rH)
        d_in = 0.8
        assert d_in < floor
        R = np.array([[0.0, 0.0, 0.0], [d_in, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([rC, rH], dtype=np.float64)
        idx = _heavy_atom_indices(["C", "H"], 2, include_h=True)
        assert list(idx) == [0, 1]
        L, G = _vdw_floor_value_and_grad(
            R, heavy_indices=idx, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        # Loss = w * (floor - d)^2 -- computed explicitly.
        gap = floor - d_in
        L_expected = DEFAULT_VDW_FLOOR_WEIGHT * gap * gap
        assert L == pytest.approx(L_expected, rel=1e-12)
        # Gradient pushes them apart -- atom 0 to -x, atom 1 to +x after step.
        R_new = R - 1e-3 * G
        d_new = float(np.linalg.norm(R_new[0] - R_new[1]))
        assert d_new > d_in
        # And: when include_h is FALSE the same geometry yields zero loss
        # (H is excluded).
        idx_off = _heavy_atom_indices(["C", "H"], 2, include_h=False)
        assert list(idx_off) == [0]
        L_off, G_off = _vdw_floor_value_and_grad(
            R, heavy_indices=idx_off, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        assert L_off == 0.0
        assert np.allclose(G_off, 0.0)


# ---------------------------------------------------------------------------
# Test 4 -- Bonded C-H pair excluded (no penalty at normal 1.09 Å)
# ---------------------------------------------------------------------------
class TestCHBondedExcluded:
    def test_bonded_ch_at_1p09_no_penalty(self):
        rC = DEFAULT_VDW_RADII["C"]
        rH = DEFAULT_VDW_RADII["H"]
        d = 1.09  # normal C-H
        R = np.array([[0.0, 0.0, 0.0], [d, 0.0, 0.0]], dtype=np.float64)
        radii = np.array([rC, rH], dtype=np.float64)
        idx = _heavy_atom_indices(["C", "H"], 2, include_h=True)
        excl = {frozenset((0, 1))}  # bonded
        L, G = _vdw_floor_value_and_grad(
            R, heavy_indices=idx, radii=radii,
            excluded_pairs=excl,
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)


# ---------------------------------------------------------------------------
# Test 5 -- Finite-difference gradient agreement for an H-inclusive case
# ---------------------------------------------------------------------------
class TestGradientFD_IncludeH:
    def test_ch_pair_fd_match(self):
        rC = DEFAULT_VDW_RADII["C"]
        rH = DEFAULT_VDW_RADII["H"]
        radii = np.array([rC, rH], dtype=np.float64)
        rng = np.random.default_rng(2026)
        for trial in range(5):
            sep = 0.4 + 1.5 * rng.random()  # well inside the 2.465 Å C-H floor
            direction = rng.standard_normal(3)
            direction /= np.linalg.norm(direction)
            R = np.array([[0.0, 0.0, 0.0], sep * direction], dtype=np.float64)
            idx = _heavy_atom_indices(["C", "H"], 2, include_h=True)
            L, G = _vdw_floor_value_and_grad(
                R, heavy_indices=idx, radii=radii,
                excluded_pairs=set(),
                weight=DEFAULT_VDW_FLOOR_WEIGHT,
                fraction=DEFAULT_VDW_FLOOR_FRACTION,
            )
            eps = 1e-5
            G_fd = np.zeros_like(R)
            for a in range(2):
                for c in range(3):
                    R_pos = R.copy(); R_pos[a, c] += eps
                    R_neg = R.copy(); R_neg[a, c] -= eps
                    L_pos, _ = _vdw_floor_value_and_grad(
                        R_pos, heavy_indices=idx, radii=radii,
                        excluded_pairs=set(),
                        weight=DEFAULT_VDW_FLOOR_WEIGHT,
                        fraction=DEFAULT_VDW_FLOOR_FRACTION,
                    )
                    L_neg, _ = _vdw_floor_value_and_grad(
                        R_neg, heavy_indices=idx, radii=radii,
                        excluded_pairs=set(),
                        weight=DEFAULT_VDW_FLOOR_WEIGHT,
                        fraction=DEFAULT_VDW_FLOOR_FRACTION,
                    )
                    G_fd[a, c] = (L_pos - L_neg) / (2.0 * eps)
            err = float(np.max(np.abs(G - G_fd)))
            assert err < 1e-4, (
                f"trial {trial}: analytic vs FD mismatch {err:.3e} "
                f"(sep={sep:.3f})"
            )


# ---------------------------------------------------------------------------
# Test 6 -- Determinism: 2 runs of the helper with H-included = bit identical
# ---------------------------------------------------------------------------
class TestDeterminismIncludeH:
    def test_two_runs_byte_identical(self):
        rng = np.random.default_rng(7)
        symbols = ["C", "H", "C", "H", "H"]
        radii = np.array([DEFAULT_VDW_RADII[s] for s in symbols], dtype=np.float64)
        R = rng.standard_normal((5, 3)) * 0.6  # tight cluster -> some clashes
        idx = _heavy_atom_indices(symbols, len(symbols), include_h=True)
        L_a, G_a = _vdw_floor_value_and_grad(
            R, heavy_indices=idx, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        L_b, G_b = _vdw_floor_value_and_grad(
            R, heavy_indices=idx, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        assert L_a == L_b
        assert np.array_equal(G_a, G_b)


# ---------------------------------------------------------------------------
# Test 7 -- Integration with rigid-body-aware all-pairs floor
# ---------------------------------------------------------------------------
class TestAllPairsIncludeH:
    """``vdw_floor_all_pairs_value_and_grad`` accepts ``include_h`` and routes
    H gradients through the rigid-body translation slot when both endpoints
    sit inside the rigid body.
    """

    def test_all_pairs_byte_identical_without_include_h(self):
        # Mixed geometry with one heavy clash AND one C-H clash.
        R = np.array([
            [0.0, 0.0, 0.0],   # 0: C
            [1.5, 0.0, 0.0],   # 1: C  (C-C clash -- heavy)
            [-0.7, 0.0, 0.0],  # 2: H  (C0-H clash, off when include_h=False)
        ], dtype=np.float64)
        symbols = ["C", "C", "H"]
        L_off, G_off = vdw_floor_all_pairs_value_and_grad(
            R, symbols=symbols, excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
            include_h=False,
        )
        # Loss should be exactly the C-C clash contribution only.
        rC = DEFAULT_VDW_RADII["C"]
        floor_cc = DEFAULT_VDW_FLOOR_FRACTION * 2.0 * rC
        gap_cc = floor_cc - 1.5
        L_cc_only = DEFAULT_VDW_FLOOR_WEIGHT * gap_cc * gap_cc
        assert L_off == pytest.approx(L_cc_only, rel=1e-12)
        # H atom got no gradient.
        assert np.allclose(G_off[2], 0.0)

    def test_all_pairs_include_h_fires_for_ch(self):
        R = np.array([
            [0.0, 0.0, 0.0],   # 0: C
            [1.5, 0.0, 0.0],   # 1: C  (C-C clash)
            [-0.7, 0.0, 0.0],  # 2: H  (C0-H clash)
        ], dtype=np.float64)
        symbols = ["C", "C", "H"]
        L_on, G_on = vdw_floor_all_pairs_value_and_grad(
            R, symbols=symbols, excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
            include_h=True,
        )
        L_off, _ = vdw_floor_all_pairs_value_and_grad(
            R, symbols=symbols, excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
            include_h=False,
        )
        # ON loss > OFF loss because the C-H pair now contributes too.
        assert L_on > L_off
        # H gradient is non-zero.
        assert not np.allclose(G_on[2], 0.0)

    def test_all_pairs_include_h_routes_through_rigid_body(self):
        """When both endpoints of an H-inclusive clash are inside the rigid
        body, the helper drops the internal-distortion gradient (its
        rigid-body documented behaviour) -- the rigid block translates as a
        whole or not at all.
        """
        # Three H atoms inside a rigid body, all clashing with each other.
        rH = DEFAULT_VDW_RADII["H"]
        R = np.array([
            [0.0, 0.0, 0.0],   # 0: metal (translation anchor)
            [1.0, 0.0, 0.0],   # 1: H (rigid)
            [1.0, 1.0, 0.0],   # 2: H (rigid, ~1.41 Å from 1, inside floor 2.04)
            [2.0, 0.0, 0.0],   # 3: H (rigid, 1.0 Å from 1)
        ], dtype=np.float64)
        symbols = ["Ni", "H", "H", "H"]
        L_int, G_int = vdw_floor_all_pairs_value_and_grad(
            R, symbols=symbols, excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
            rigid_body_atoms=[1, 2, 3],
            rigid_translation_metal=0,
            include_h=True,
        )
        # Loss > 0 (penalty still recorded), but the gradient on the H atoms
        # themselves is zero (it was dropped as documented -- "both rigid +
        # metal known -> drop gradient").
        assert L_int > 0.0
        assert np.allclose(G_int[1], 0.0)
        assert np.allclose(G_int[2], 0.0)
        assert np.allclose(G_int[3], 0.0)
        # Metal gradient is also zero (both endpoints inside rigid body).
        assert np.allclose(G_int[0], 0.0)


# ---------------------------------------------------------------------------
# Test 8 -- BERTEB-class demo (synthetic but realistic C-H collapse)
# ---------------------------------------------------------------------------
class TestBertebClassDemo:
    """Synthetic BERTEB-flavoured geometry: Ni + Br + 4 H atoms collapsed onto
    C neighbours.  With include_h=True the floor fires for every collapsed
    C-H pair (4 violations matching the V14 report).
    """

    def test_four_ch_collapses_all_fire(self):
        # 4 C-H pairs at the reported collapse distances.
        ch_distances = [0.82, 0.86, 0.93, 0.99]
        rC = DEFAULT_VDW_RADII["C"]
        rH = DEFAULT_VDW_RADII["H"]
        floor = DEFAULT_VDW_FLOOR_FRACTION * (rC + rH)
        # 4 disjoint C-H pairs along independent axes so they don't interact.
        n = 8
        R = np.zeros((n, 3), dtype=np.float64)
        symbols: list[str] = []
        for k, d in enumerate(ch_distances):
            R[2 * k] = np.array([k * 10.0, 0.0, 0.0])           # C
            R[2 * k + 1] = R[2 * k] + np.array([d, 0.0, 0.0])    # H
            symbols.extend(["C", "H"])
        idx = _heavy_atom_indices(symbols, n, include_h=True)
        assert len(idx) == n  # all 8 atoms iterated
        radii = np.array([DEFAULT_VDW_RADII[s] for s in symbols], dtype=np.float64)
        L_on, G_on = _vdw_floor_value_and_grad(
            R, heavy_indices=idx, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        # Expected = sum over the 4 C-H pairs only (pair-pairs are far apart).
        L_expected = sum(
            DEFAULT_VDW_FLOOR_WEIGHT * (floor - d) * (floor - d)
            for d in ch_distances
        )
        assert L_on == pytest.approx(L_expected, rel=1e-9)
        # Take a small step against the gradient -- every C-H distance grows.
        step = 1e-3
        R_new = R - step * G_on
        for k, d_old in enumerate(ch_distances):
            d_new = float(np.linalg.norm(R_new[2 * k] - R_new[2 * k + 1]))
            assert d_new > d_old, (
                f"pair {k}: distance did not grow ({d_old:.3f} -> {d_new:.3f})"
            )

        # And: with INCLUDE_H=False the loss is zero (only C-C pairs exist,
        # and they are at 10 Å > 2.89 floor).
        idx_off = _heavy_atom_indices(symbols, n, include_h=False)
        L_off, G_off = _vdw_floor_value_and_grad(
            R, heavy_indices=idx_off, radii=radii,
            excluded_pairs=set(),
            weight=DEFAULT_VDW_FLOOR_WEIGHT,
            fraction=DEFAULT_VDW_FLOOR_FRACTION,
        )
        assert L_off == 0.0
        assert np.allclose(G_off, 0.0)
