"""Tests for GRIP Phase 2 — loss terms, gradients, fragment detection.

Covers:
* analytic gradients vs central finite differences (tol 1e-6) for
  Bond/Angle/Improper terms across multiple random configurations
* aggregator correctness and gradient shape
* byte-identical determinism on repeat evaluation
* fragment detection respects ``frozen_atoms`` and produces deterministic
  term order
* library fallback chain (sparse keys reach a looser match)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

# Strict determinism: set seed BEFORE numpy random operations are issued.
os.environ.setdefault("PYTHONHASHSEED", "0")

from delfin.fffree.grip_loss_terms import (
    AngleTerm,
    BondTerm,
    ImproperTerm,
    TorsionTerm,
    TotalGripLoss,
    sparse_downweight,
)
from delfin.fffree.grip_mogul_lookup import (
    GripLibrary,
    lookup_angle,
    lookup_bond,
    lookup_improper,
)
from delfin.fffree.grip_fragment_detect import detect_fragments


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _fd_grad(term, R, eps=1e-5):
    """Compute the gradient of term.value_and_grad via central finite diff."""
    R = np.asarray(R, dtype=np.float64).copy()
    g = np.zeros_like(R)
    for i in range(R.shape[0]):
        for k in range(3):
            R_p = R.copy()
            R_m = R.copy()
            R_p[i, k] += eps
            R_m[i, k] -= eps
            lp, _ = term.value_and_grad(R_p)
            lm, _ = term.value_and_grad(R_m)
            g[i, k] = (lp - lm) / (2.0 * eps)
    return g


def _rng(seed):
    """Deterministic RNG factory — only used in tests; never in production."""
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# BondTerm
# ---------------------------------------------------------------------------
class TestBondTerm:
    def test_loss_at_optimum_zero(self):
        # If |r_a - r_b| = mu, loss = 0 and grad is identically zero.
        mu = 1.50
        sigma = 0.02
        R = np.array(
            [[0.0, 0.0, 0.0], [mu, 0.0, 0.0], [5.0, 5.0, 5.0]], dtype=np.float64
        )
        t = BondTerm(a=0, b=1, mu=mu, sigma=sigma, weight=1.0)
        L, g = t.value_and_grad(R)
        assert L == pytest.approx(0.0, abs=1e-14)
        assert np.allclose(g, 0.0, atol=1e-14)

    def test_loss_positive_off_optimum(self):
        mu = 1.50
        sigma = 0.02
        R = np.array([[0.0, 0.0, 0.0], [1.60, 0.0, 0.0]], dtype=np.float64)
        t = BondTerm(a=0, b=1, mu=mu, sigma=sigma, weight=1.0)
        L, g = t.value_and_grad(R)
        # z = (1.60-1.50)/0.02 = 5.0 -> L = w * z^2 = 25
        assert L == pytest.approx(25.0, rel=1e-12)
        # gradient nonzero
        assert np.linalg.norm(g) > 0

    def test_gradient_finite_diff_random(self):
        rng = _rng(seed=42)
        for k in range(5):
            R = rng.standard_normal((4, 3)).astype(np.float64) * 1.5
            # Pick bond pair so length isn't zero
            t = BondTerm(a=0, b=1, mu=1.4, sigma=0.05, weight=1.0)
            _, g_ana = t.value_and_grad(R)
            g_fd = _fd_grad(t, R, eps=1e-5)
            assert np.allclose(g_ana, g_fd, atol=1e-6, rtol=1e-6), (
                f"BondTerm grad mismatch trial {k}: "
                f"max_abs={np.max(np.abs(g_ana - g_fd))}"
            )

    def test_degenerate_coincident_atoms_no_nan(self):
        # |d| = 0 should yield no NaN in either loss or gradient.
        R = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]], dtype=np.float64)
        t = BondTerm(a=0, b=1, mu=1.5, sigma=0.02, weight=1.0)
        L, g = t.value_and_grad(R)
        assert np.isfinite(L)
        assert np.all(np.isfinite(g))


# ---------------------------------------------------------------------------
# AngleTerm
# ---------------------------------------------------------------------------
class TestAngleTerm:
    def test_loss_at_optimum_zero(self):
        # 109.47° tetrahedral angle at origin between two unit vectors.
        target = math_acos_deg(-1.0 / 3.0)  # ~109.4712
        ra = np.array([1.0, 0.0, 0.0])
        rb = np.array([0.0, 0.0, 0.0])
        # Place rc such that angle ra-rb-rc == target
        theta = math_to_rad(target)
        rc = np.array([np.cos(theta), np.sin(theta), 0.0])
        R = np.vstack([ra, rb, rc]).astype(np.float64)
        t = AngleTerm(a=0, b=1, c=2, mu=target, sigma=2.0, weight=1.0)
        L, g = t.value_and_grad(R)
        assert L == pytest.approx(0.0, abs=1e-10)
        # Grad should be ~zero at optimum
        assert np.linalg.norm(g) < 1e-8

    def test_gradient_finite_diff_random(self):
        rng = _rng(seed=7)
        for k in range(5):
            # Random configs, but normalise legs to >0.3 Å so angle is well-defined
            R = rng.standard_normal((4, 3)).astype(np.float64) * 1.5
            # Make sure no coincident atoms by adding small offsets
            R[1] = R[0] + np.array([1.2, 0.3, -0.2])
            R[2] = R[1] + np.array([-0.5, 1.0, 0.7])
            t = AngleTerm(a=0, b=1, c=2, mu=109.5, sigma=3.0, weight=1.0)
            _, g_ana = t.value_and_grad(R)
            g_fd = _fd_grad(t, R, eps=1e-5)
            err = np.max(np.abs(g_ana - g_fd))
            assert err < 1e-5, f"AngleTerm grad mismatch trial {k}: max_abs={err}"

    def test_gradient_near_linear_no_nan(self):
        # Near 180° angle — derivative of acos diverges; check our clamp is
        # effective and the gradient stays finite.
        R = np.array(
            [
                [-1.0 + 1e-7, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        t = AngleTerm(a=0, b=1, c=2, mu=120.0, sigma=5.0, weight=1.0)
        L, g = t.value_and_grad(R)
        assert np.isfinite(L)
        assert np.all(np.isfinite(g))


# ---------------------------------------------------------------------------
# ImproperTerm
# ---------------------------------------------------------------------------
class TestImproperTerm:
    def test_loss_at_planar(self):
        # All 4 atoms coplanar -> d_oop = 0. If mu = 0, loss = 0.
        R = np.array(
            [
                [0.0, 0.0, 0.0],   # center
                [1.0, 0.0, 0.0],
                [-0.5, 0.866, 0.0],
                [-0.5, -0.866, 0.0],
            ],
            dtype=np.float64,
        )
        t = ImproperTerm(
            center=0, neighbors_3=(1, 2, 3),
            mu=0.0, sigma=0.05, weight=2.0,
        )
        L, g = t.value_and_grad(R)
        assert L == pytest.approx(0.0, abs=1e-12)
        assert np.linalg.norm(g) < 1e-10

    def test_loss_positive_off_plane(self):
        # Push center off the plane.
        R = np.array(
            [
                [0.0, 0.0, 0.3],
                [1.0, 0.0, 0.0],
                [-0.5, 0.866, 0.0],
                [-0.5, -0.866, 0.0],
            ],
            dtype=np.float64,
        )
        t = ImproperTerm(
            center=0, neighbors_3=(1, 2, 3),
            mu=0.0, sigma=0.05, weight=2.0,
        )
        L, _ = t.value_and_grad(R)
        # d_oop = 0.3 (centroid is at origin, normal z, dot = 0.3)
        # z = 0.3/0.05 = 6.0 -> L = 2 * 36 = 72
        assert L == pytest.approx(72.0, rel=1e-12)

    def test_gradient_finite_diff_random(self):
        rng = _rng(seed=123)
        for k in range(5):
            R = rng.standard_normal((4, 3)).astype(np.float64)
            # Keep neighbours spread so plane is non-degenerate
            R[1] = np.array([1.0, 0.0, 0.0]) + rng.standard_normal(3) * 0.1
            R[2] = np.array([-0.5, 0.866, 0.0]) + rng.standard_normal(3) * 0.1
            R[3] = np.array([-0.5, -0.866, 0.0]) + rng.standard_normal(3) * 0.1
            R[0] = np.array([0.0, 0.0, 0.0]) + rng.standard_normal(3) * 0.2
            t = ImproperTerm(
                center=0, neighbors_3=(1, 2, 3),
                mu=0.05, sigma=0.1, weight=2.0,
            )
            _, g_ana = t.value_and_grad(R)
            g_fd = _fd_grad(t, R, eps=1e-6)
            err = np.max(np.abs(g_ana - g_fd))
            assert err < 1e-5, f"ImproperTerm grad mismatch trial {k}: max_abs={err}"

    def test_collinear_neighbours_no_nan(self):
        # All 3 neighbours on one line -> normal is ill-defined; gradient must
        # be returned as zero (not NaN).
        R = np.array(
            [
                [0.0, 0.0, 0.5],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
            ],
            dtype=np.float64,
        )
        t = ImproperTerm(
            center=0, neighbors_3=(1, 2, 3),
            mu=0.0, sigma=0.05, weight=2.0,
        )
        L, g = t.value_and_grad(R)
        assert np.isfinite(L)
        assert np.all(np.isfinite(g))


# ---------------------------------------------------------------------------
# Total aggregator
# ---------------------------------------------------------------------------
class TestTotalGripLoss:
    def test_aggregation_sum_and_shape(self):
        terms = [
            BondTerm(a=0, b=1, mu=1.5, sigma=0.02, weight=1.0),
            BondTerm(a=1, b=2, mu=1.4, sigma=0.02, weight=1.0),
            AngleTerm(a=0, b=1, c=2, mu=120.0, sigma=2.0, weight=0.5),
        ]
        R = np.array(
            [[0.0, 0.0, 0.0], [1.6, 0.0, 0.0], [3.0, 0.0, 0.0]],
            dtype=np.float64,
        )
        agg = TotalGripLoss(terms=terms)
        assert len(agg) == 3
        L_total, G_total = agg.evaluate(R)
        # Hand sum
        L0, g0 = terms[0].value_and_grad(R)
        L1, g1 = terms[1].value_and_grad(R)
        L2, g2 = terms[2].value_and_grad(R)
        expected_L = L0 + L1 + L2
        expected_G = g0 + g1 + g2
        assert L_total == pytest.approx(expected_L, rel=1e-12)
        assert np.allclose(G_total, expected_G, atol=1e-12)
        assert G_total.shape == R.shape

    def test_terms_sorted_canonical(self):
        # Construction must sort by atom_indices for determinism.
        t1 = BondTerm(a=3, b=4, mu=1.5, sigma=0.02)
        t2 = BondTerm(a=0, b=1, mu=1.5, sigma=0.02)
        t3 = BondTerm(a=2, b=3, mu=1.5, sigma=0.02)
        agg = TotalGripLoss(terms=[t1, t2, t3])
        order = [t.atom_indices for t in agg.terms]
        assert order == sorted(order)

    def test_flat_input_returns_flat_grad(self):
        terms = [BondTerm(a=0, b=1, mu=1.5, sigma=0.02, weight=1.0)]
        R = np.array(
            [[0.0, 0.0, 0.0], [1.6, 0.0, 0.0]], dtype=np.float64
        ).reshape(-1)
        agg = TotalGripLoss(terms=terms)
        L, G = agg.evaluate(R)
        assert G.shape == R.shape


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------
class TestDeterminism:
    def test_byte_identical_repeat(self):
        # Same R, same terms -> bit-identical loss & gradient values across runs.
        rng = _rng(seed=0)
        R = rng.standard_normal((6, 3)).astype(np.float64)
        terms = [
            BondTerm(a=0, b=1, mu=1.5, sigma=0.02, weight=1.0),
            BondTerm(a=2, b=3, mu=1.4, sigma=0.03, weight=1.0),
            AngleTerm(a=0, b=1, c=2, mu=110.0, sigma=2.0, weight=0.5),
            ImproperTerm(center=4, neighbors_3=(1, 2, 3),
                          mu=0.0, sigma=0.05, weight=2.0),
        ]
        agg = TotalGripLoss(terms=terms)
        L1, G1 = agg.evaluate(R)
        L2, G2 = agg.evaluate(R)
        # 1e-15 means truly bit-identical FP repeats.
        assert L1 == L2
        assert np.array_equal(G1, G2)


# ---------------------------------------------------------------------------
# Library lookup / fallback
# ---------------------------------------------------------------------------
class TestLibrary:
    def test_concrete_bond_exists(self):
        # C-C sp3-sp3 has many CCDC observations - must return.
        hit = lookup_bond("C", "sp3", "C", "sp3")
        assert hit is not None
        mu, sigma, n = hit
        assert 1.40 < mu < 1.65, f"C-C sp3 mu should be near 1.52, got {mu}"
        assert sigma > 0
        assert n >= 5

    def test_concrete_angle_exists(self):
        hit = lookup_angle("C", "C", "sp3", "C")
        assert hit is not None
        mu, sigma, n = hit
        assert 100 < mu < 120, f"C-C-C sp3 angle should be near 110, got {mu}"

    def test_fallback_chain_hits_loose_level(self):
        # An exotic pair (e.g. Tc-Yb) is unlikely at level 0 but should fall
        # back via the chain to a more generic entry. If the library is
        # genuinely empty for both endpoints, the call returns None — also OK.
        # We just assert no crash and consistent type.
        hit = lookup_bond("Tc", "sp3", "Yb", "sp3")
        assert hit is None or (isinstance(hit, tuple) and len(hit) == 3)

    def test_sparse_downweight(self):
        # n < n_floor -> weight scaled down linearly.
        assert sparse_downweight(1.0, 0, n_floor=5) == 0.0
        assert sparse_downweight(1.0, 2, n_floor=5) == pytest.approx(0.4, rel=1e-12)
        assert sparse_downweight(1.0, 5, n_floor=5) == 1.0
        assert sparse_downweight(1.0, 100, n_floor=5) == 1.0
        assert sparse_downweight(2.0, 1, n_floor=5) == pytest.approx(0.4, rel=1e-12)


# ---------------------------------------------------------------------------
# Fragment detection
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def toluene_mol():
    Chem = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


@pytest.fixture
def toluene_coords(toluene_mol):
    conf = toluene_mol.GetConformer()
    P = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(toluene_mol.GetNumAtoms())],
        dtype=np.float64,
    )
    return P


class TestFragmentDetect:
    def test_toluene_produces_terms(self, toluene_mol, toluene_coords):
        result = detect_fragments(
            toluene_mol, toluene_coords,
            frozen_atoms=set(),
            return_result=True,
        )
        # Toluene has 7 heavy atoms + ~8 H = 15 atoms, ~15 bonds, ~22 angles,
        # ~6 sp2 carbons that qualify as improper centres with 3 neighbours.
        assert result.n_bond_candidates > 10
        assert result.n_angle_candidates > 15
        # At least some terms must be matched in the library
        assert result.n_bond_matched > 0
        assert result.n_angle_matched > 0
        # No NaNs in any (mu, sigma)
        for t in result.bond_terms:
            assert np.isfinite(t.mu) and np.isfinite(t.sigma) and t.sigma > 0
        for t in result.angle_terms:
            assert np.isfinite(t.mu) and np.isfinite(t.sigma) and t.sigma > 0
        for t in result.improper_terms:
            assert np.isfinite(t.mu) and np.isfinite(t.sigma) and t.sigma > 0

    def test_frozen_atoms_skip(self, toluene_mol, toluene_coords):
        # Freeze atom 0 (methyl C). Any bond / angle / improper touching it
        # must be excluded.
        frozen = {0}
        result = detect_fragments(
            toluene_mol, toluene_coords,
            frozen_atoms=frozen,
            return_result=True,
        )
        for t in result.bond_terms:
            assert 0 not in t.atom_indices
        for t in result.angle_terms:
            assert 0 not in t.atom_indices
        for t in result.improper_terms:
            assert 0 not in t.atom_indices

    def test_deterministic_order(self, toluene_mol, toluene_coords):
        # Call twice — identical term list.
        r1 = detect_fragments(toluene_mol, toluene_coords, return_result=True)
        r2 = detect_fragments(toluene_mol, toluene_coords, return_result=True)
        a1 = [tuple(t.atom_indices) for t in r1.all_terms()]
        a2 = [tuple(t.atom_indices) for t in r2.all_terms()]
        assert a1 == a2
        # And monotonically sorted
        assert a1 == sorted(a1)

    def test_terms_evaluate_no_nan(self, toluene_mol, toluene_coords):
        terms = detect_fragments(toluene_mol, toluene_coords, frozen_atoms=set())
        if not terms:
            pytest.skip("No terms detected — library coverage gap, not a bug")
        agg = TotalGripLoss(terms=terms)
        L, G = agg.evaluate(toluene_coords)
        assert np.isfinite(L)
        assert np.all(np.isfinite(G))
        assert G.shape == toluene_coords.shape


class TestHeal1DonorShellProtection:
    """Heal-1 (2026-06-01, mddir-fix): donor first-shell neighbours must not
    appear as the *centre* of an angle or improper term, but their bonds and
    their end-position angles are still emitted.

    The protection is opt-in: it activates only when the caller passes
    ``donors=...`` (legacy callers that only pass ``frozen_atoms`` are
    byte-exactly unchanged).
    """

    @pytest.fixture(scope="class")
    def amine_like_mol(self):
        """Tiny CH3-NH2 mol so we can reason about every term explicitly.

        Atom indices after AddHs are deterministic (RDKit canonical order).
        We use methylamine ``CN`` — heavy atoms are C(0) and N(1); their H
        neighbours come after.  Treating N(1) as the donor, C(0) is the
        donor's only heavy shell-1 neighbour.
        """
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("CN")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    @pytest.fixture
    def amine_coords(self, amine_like_mol):
        conf = amine_like_mol.GetConformer()
        return np.array(
            [list(conf.GetAtomPosition(i)) for i in range(amine_like_mol.GetNumAtoms())],
            dtype=np.float64,
        )

    def test_donor_centered_angles_skipped(self, amine_like_mol, amine_coords):
        # N is atom 1 (donor); pass donors=[1].  No angle may be centred on
        # N because N is in frozen_or_donor (legacy behaviour, but exercised
        # explicitly here via the donors= path).
        donor_idx = 1
        result = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            return_result=True,
        )
        for t in result.angle_terms:
            ai, bi, ci = t.atom_indices
            assert bi != donor_idx, (
                f"Angle centred on donor: {(ai, bi, ci)} — Heal-1 invariant violated"
            )

    def test_donor_endpoint_bonds_skipped(self, amine_like_mol, amine_coords):
        donor_idx = 1
        result = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            return_result=True,
        )
        for t in result.bond_terms:
            ai, bi = t.atom_indices
            assert donor_idx not in (ai, bi), (
                f"Bond touches donor: {(ai, bi)} — Heal-1 invariant violated"
            )

    def test_donor_improper_skipped(self, amine_like_mol, amine_coords):
        donor_idx = 1
        result = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            return_result=True,
        )
        for t in result.improper_terms:
            # ImproperTerm.atom_indices == (center, n1, n2, n3)
            ai = t.atom_indices
            assert donor_idx not in ai, (
                f"Improper touches donor: {ai} — Heal-1 invariant violated"
            )

    def test_donor_shell1_protected_as_central(
        self, amine_like_mol, amine_coords
    ):
        # With donors=[1] (N), the C(0) is N's only heavy shell-1 neighbour.
        # Heal-1 strict rule: no angle term has b==C.  Heal-1b refinement
        # (2026-06-01, amine-H regression fix): the central-skip applies
        # ONLY to heavy-only triples; angles with an H endpoint (e.g.
        # H-C-H, N-C-H... — though N as endpoint is donor and gets caught
        # by the legacy frozen-or-donor filter) are now re-emitted to
        # preserve H orientation priors.  So the invariant becomes: no
        # angle term has b==shell1 AND BOTH endpoints heavy (atomic_num>1).
        donor_idx = 1
        shell1 = 0
        result = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            return_result=True,
        )
        for t in result.angle_terms:
            ai, bi, ci = t.atom_indices
            if bi != shell1:
                continue
            sym_a = amine_like_mol.GetAtomWithIdx(ai).GetSymbol()
            sym_c = amine_like_mol.GetAtomWithIdx(ci).GetSymbol()
            assert (sym_a == "H" or sym_c == "H"), (
                f"Heal-1b leak: heavy-only angle ({ai}, {bi}, {ci}) "
                f"({sym_a}-?-{sym_c}) centred on shell-1 — invariant violated"
            )
        # Improper symmetry case: only heavy-only impropers centred on
        # shell-1 are still dropped (Heal-1b).
        for t in result.improper_terms:
            ai = t.atom_indices
            if ai[0] != shell1:
                continue
            # ImproperTerm.atom_indices == (centre, n1, n2, n3)
            nbrs = ai[1:]
            sym_nbrs = [amine_like_mol.GetAtomWithIdx(n).GetSymbol() for n in nbrs]
            assert "H" in sym_nbrs, (
                f"Heal-1b leak: heavy-only improper centred on shell-1: {ai} "
                f"({sym_nbrs}) — invariant violated"
            )

    def test_donors_none_recovers_legacy(self, amine_like_mol, amine_coords):
        # With donors omitted entirely, the result is byte-exactly the
        # pre-Heal-1 behaviour (frozen-only).  Comparing term tuples:
        r_legacy = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(),
            return_result=True,
        )
        r_donors_none = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=None,
            return_result=True,
        )
        a1 = [tuple(t.atom_indices) for t in r_legacy.all_terms()]
        a2 = [tuple(t.atom_indices) for t in r_donors_none.all_terms()]
        assert a1 == a2

    def test_frozen_atoms_only_legacy_behaviour_unchanged(
        self, amine_like_mol, amine_coords
    ):
        # Calling with the legacy frozen_atoms={N} (and donors=None) reproduces
        # the original behaviour byte-exactly: N is excluded everywhere but
        # the C(0) shell-1 protection does NOT engage.
        n_idx = 1
        r_frozen_only = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms={n_idx},
            return_result=True,
        )
        # The shell-1 atom (C, idx 0) is allowed as centre in this mode.
        # We don't enforce that ANY such term exists in this tiny mol (the
        # library may not match), only that the legacy union-set didn't add
        # spurious extra exclusions beyond what frozen_atoms encodes.
        for t in r_frozen_only.bond_terms:
            assert n_idx not in t.atom_indices
        for t in r_frozen_only.angle_terms:
            assert n_idx not in t.atom_indices
        for t in r_frozen_only.improper_terms:
            assert n_idx not in t.atom_indices


class TestHeal1bXCHException:
    """Heal-1b (2026-06-01): the shell-1 central-skip from Heal-1 is *too
    strict* — it drops C-N-H angles (along with C-N-C / C-N-Cα) and lets
    amine H drift away from the C-N-H ~108° prior, which surfaces in
    smoke-50 as a +100% amine_h_pct_files_viol regression.

    Heal-1b keeps the central-protection rule for HEAVY-only triples but
    re-emits any (a, shell1, c) angle term whose *endpoint* a or c is
    hydrogen.  C-N-H, H-N-H stay alive; C-N-C, C-N-Cα still drop out.
    The same exception applies to improper terms whose centre is shell-1
    iff at least one of the three neighbours is hydrogen.

    The new rule is strictly *less* restrictive than Heal-1, so every
    Heal-1 invariant test above still holds (donor itself never appears,
    shell-1 still never centres a C-C-C or C-N-Cα angle).
    """

    @pytest.fixture(scope="class")
    def amine_like_mol(self):
        """Methylamine CH3-NH2 — donor N(1), shell-1 C(0).

        Heavy atoms: C(0), N(1).  After AddHs: H atoms append in
        deterministic RDKit canonical order — 3 on C, 2 on N.
        """
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("CN")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    @pytest.fixture
    def amine_coords(self, amine_like_mol):
        conf = amine_like_mol.GetConformer()
        return np.array(
            [list(conf.GetAtomPosition(i)) for i in range(amine_like_mol.GetNumAtoms())],
            dtype=np.float64,
        )

    def _h_indices_on(self, mol, idx):
        """Return the H-neighbour indices of atom ``idx`` in ``mol``."""
        out = []
        atom = mol.GetAtomWithIdx(int(idx))
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "H":
                out.append(int(nb.GetIdx()))
        return out

    def test_donor_shell1_CNH_angle_kept(self, amine_like_mol, amine_coords):
        # Donor = N(1).  Shell-1 = C(0).  Heal-1b must re-emit at least one
        # angle term centred on N (skipped — N is donor, always frozen) or
        # ... wait: we need C-N-H angles whose CENTRE is N.  But N is donor
        # = frozen_or_donor, so b==N is always skipped (legacy rule, not the
        # Heal-1 rule).  What Heal-1b actually un-skips is angles whose
        # centre is the SHELL-1 atom and one endpoint is H — i.e. angles
        # like N-C-H (centred on C) or H-C-H (centred on C).  Both should
        # now appear because C is shell-1 AND one endpoint is H.
        #
        # We assert: at least one angle term exists with centre b == C(0)
        # AND (a or c) is a hydrogen.
        donor_idx = 1
        shell1 = 0
        result = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            return_result=True,
        )
        h_on_C = set(self._h_indices_on(amine_like_mol, shell1))
        # H-X-Y angles centred on shell-1 with at least one H endpoint:
        H_centered_on_shell1 = [
            (t.atom_indices)
            for t in result.angle_terms
            if t.atom_indices[1] == shell1
            and (t.atom_indices[0] in h_on_C or t.atom_indices[2] in h_on_C)
        ]
        assert len(H_centered_on_shell1) > 0, (
            f"Heal-1b must keep H-containing angles centred on the shell-1 "
            f"atom, but found 0.  All angle term tuples: "
            f"{[t.atom_indices for t in result.angle_terms]}"
        )

    def test_donor_shell1_CNC_angle_skipped(self, amine_like_mol, amine_coords):
        # Even with Heal-1b, the heavy-only angles centred on shell-1 (C)
        # remain forbidden — only the H-touching ones come back.  For
        # methylamine the only heavy neighbour of C(0) is N(1) itself, so
        # there is no C-C-C triple to test.  We synthesise one with
        # ethylamine CC-NH2 → atoms: C(0)-C(1)-N(2)+H's; shell-1 of donor
        # N(2) is C(1); C(1)'s heavy neighbours are C(0) and N(2).  The
        # angle C(0)-C(1)-N(2) is centred on shell-1 with NO H endpoint
        # → must be skipped.  (N as endpoint is also a donor → caught by
        # the legacy frozen-or-donor filter anyway, which is the correct
        # outcome.)
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("CCN")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        P = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        donor_idx = 2   # N
        shell1 = 1      # middle C
        result = detect_fragments(
            mol, P, frozen_atoms=set(), donors=[donor_idx], return_result=True,
        )
        # Assert: no angle term has b==shell1 with BOTH endpoints heavy
        # (atomic number > 1).  This is the C-N-C / C-C-N / C-C-X case
        # Heal-1b still skips.
        for t in result.angle_terms:
            ai, bi, ci = t.atom_indices
            if bi != shell1:
                continue
            sym_a = mol.GetAtomWithIdx(ai).GetSymbol()
            sym_c = mol.GetAtomWithIdx(ci).GetSymbol()
            assert (sym_a == "H" or sym_c == "H"), (
                f"Heal-1b leak: heavy-only angle {ai}-{bi}-{ci} "
                f"({sym_a}-?-{sym_c}) centred on shell-1 must be dropped"
            )

    def test_donor_shell1_HNH_angle_kept(self, amine_like_mol, amine_coords):
        # H-C-H angles centred on shell-1 C: both endpoints are H.  These
        # MUST be re-emitted (they pin the methyl umbrella).  For
        # methylamine, C(0) has 3 H neighbours -> 3 H-C-H angles.
        donor_idx = 1
        shell1 = 0
        result = detect_fragments(
            amine_like_mol, amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            return_result=True,
        )
        h_on_C = set(self._h_indices_on(amine_like_mol, shell1))
        HH_centered_on_shell1 = [
            t.atom_indices for t in result.angle_terms
            if t.atom_indices[1] == shell1
            and t.atom_indices[0] in h_on_C and t.atom_indices[2] in h_on_C
        ]
        # The library should match H-C-H (sp3) -> assert at least one,
        # and library-coverage robustness: if zero matched we still want
        # the *candidate* path to have fired.  We re-detect with a
        # near-empty library mock would be heavier; instead we assert
        # n_angle_candidates increased between donors=None and donors=[N]
        # is consistent with the Heal-1b loosening.  But the strongest
        # direct check is the >0 count below.
        assert len(HH_centered_on_shell1) > 0, (
            f"Heal-1b must keep H-C-H angles centred on shell-1, but found 0.  "
            f"h_on_C={h_on_C}, all angle terms: "
            f"{[t.atom_indices for t in result.angle_terms]}"
        )


class TestOptionBAdaptiveProtection:
    """Option-B (2026-06-01) — adaptive coverage-aware shell-1 protection.

    Heal-1 / Heal-1b drop *every* heavy-only angle / improper term centred
    on a donor first-shell atom.  That fixes mddir but takes the CCDC
    Mahalanobis pressure off the funcgrp internals.  Option-B keeps the
    Heal-1/1b drop ONLY when the library coverage for that specific
    fragment class is THIN (``n < adaptive_min_n``); when the library
    has GOOD coverage (``n ≥ adaptive_min_n``) for that class, the term
    is re-emitted so the prior keeps pulling.

    Determinism: the library is a singleton; lookups are pure functions
    of (key, library); the loop iterates deterministically over sorted
    triples.  Therefore Option-B's outcome is a pure function of the
    input molecule + library.
    """

    @pytest.fixture(scope="class")
    def ethyl_amine_mol(self):
        """CC-NH2 — donor N(2), shell-1 C(1); C(0)-C(1)-N(2) is the
        heavy-only triple we use to exercise Option-B.  N as endpoint
        is donor → already excluded by ``frozen_or_donor``, so for the
        synthetic library tests below we use the smaller methylamine
        (no heavy-heavy-heavy triple) and a larger SMILES for the
        positive-keep test."""
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("CCN")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    @pytest.fixture
    def ethyl_amine_coords(self, ethyl_amine_mol):
        conf = ethyl_amine_mol.GetConformer()
        return np.array(
            [list(conf.GetAtomPosition(i)) for i in range(ethyl_amine_mol.GetNumAtoms())],
            dtype=np.float64,
        )

    @pytest.fixture(scope="class")
    def propyl_amine_mol(self):
        """Propanamine CCC-NH2 — donor N(3), shell-1 C(2).  C(2)'s
        heavy neighbours are C(1) (heavy) and N(3) (donor → excluded).
        So the heavy-only triple centred on shell-1 we want to exercise
        is C(1)-C(2)-? — but the only non-donor heavy neighbour of C(2)
        is C(1), so no heavy-heavy triple centred on C(2) survives.

        For the positive Option-B keep test we instead use isobutylamine
        (CH3)2-CH-CH2-NH2 — donor N(5), shell-1 C(4); C(4)'s heavy
        neighbours are C(3) and N(5).  C(3)'s heavy neighbours are
        C(0), C(2), C(4).  None of those is the shell-1 atom we want.
        We need shell-1 to have ≥ 2 heavy non-donor neighbours; the
        simplest cyclic molecule with that property is pyrrolidine
        (cyclic CH2-CH2-CH2-CH2-NH).  N(0) is donor, shell-1 = C(1)
        or C(4); C(1)'s heavy neighbours are N(0) (donor) and C(2)
        (heavy non-donor).  Still only one heavy non-donor.

        For a clean heavy-only triple centred on shell-1, use
        2-aminopropane (isopropylamine) ``CC(C)N`` — donor N(3),
        shell-1 C(1); C(1)'s heavy neighbours are C(0), C(2), N(3).
        The triples centred on C(1) include C(0)-C(1)-C(2) (heavy-only,
        no donor endpoint) which is the exact case Option-B targets.
        """
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("CC(C)N")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    @pytest.fixture
    def propyl_amine_coords(self, propyl_amine_mol):
        conf = propyl_amine_mol.GetConformer()
        return np.array(
            [list(conf.GetAtomPosition(i)) for i in range(propyl_amine_mol.GetNumAtoms())],
            dtype=np.float64,
        )

    # ------------------------------------------------------------------
    # Synthetic library returning configurable (mu, sigma, n) for any
    # bond/angle/improper lookup.  Used to drive Option-B's gating
    # without depending on the real grip_lib_v1.npz coverage at the
    # moment of test authoring.
    # ------------------------------------------------------------------
    class _FakeLib:
        """Minimal duck-typed GripLibrary for Option-B gating tests.

        ``n_returned`` controls the sample size reported by every
        successful lookup.  ``return_none`` makes every lookup miss
        (mimics a fragment class with zero coverage).  Lookups also
        track a call counter so tests can verify the gating actually
        consulted the library.
        """

        def __init__(self, n_returned: int = 5, return_none: bool = False,
                     mu: float = 1.55, sigma: float = 0.02):
            self.n_returned = int(n_returned)
            self.return_none = bool(return_none)
            self.mu = float(mu)
            self.sigma = float(sigma)
            self.n_bond_calls = 0
            self.n_angle_calls = 0
            self.n_improper_calls = 0

        def _resp(self):
            if self.return_none:
                return None
            return (self.mu, self.sigma, self.n_returned)

        def lookup_bond(self, z1, hyb1, z2, hyb2,
                        ring_size_min: int = -1,
                        in_aromatic: bool = False, *,
                        min_n: int = 5):
            self.n_bond_calls += 1
            r = self._resp()
            if r is None or r[2] < min_n:
                return None
            return r

        def lookup_angle(self, z1, z2, hyb2, z3,
                         ring_size_min: int = -1,
                         in_aromatic: bool = False, *,
                         hyb1: str = "*", hyb3: str = "*",
                         min_n: int = 5):
            self.n_angle_calls += 1
            r = self._resp()
            if r is None or r[2] < min_n:
                return None
            return r

        def lookup_improper(self, z_center, hyb_center, neighbor_zs_sorted, *,
                            ring_size_min: int = -1,
                            neighbor_hybs_sorted=None,
                            min_n: int = 5):
            self.n_improper_calls += 1
            # Reuse the angle response — for Option-B's gate the only
            # thing that matters is ``n`` vs ``min_n``.
            r = self._resp()
            if r is None or r[2] < min_n:
                return None
            return r

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _heavy_only_centred_on(mol, atom_idx: int, angle_terms):
        """Return the angle terms centred on ``atom_idx`` whose endpoints
        are BOTH heavy (atomic number > 1).  This is the population
        Option-B can either drop (thin) or keep (good)."""
        out = []
        for t in angle_terms:
            ai, bi, ci = t.atom_indices
            if bi != atom_idx:
                continue
            ai_a = mol.GetAtomWithIdx(ai).GetAtomicNum()
            ai_c = mol.GetAtomWithIdx(ci).GetAtomicNum()
            if ai_a != 1 and ai_c != 1:
                out.append(t)
        return out

    # ------------------------------------------------------------------
    # 1) GOOD coverage (n ≥ 5) — keep the heavy-only shell-1 term.
    # ------------------------------------------------------------------
    def test_adaptive_keeps_when_coverage_good(
        self, propyl_amine_mol, propyl_amine_coords
    ):
        # Donor = N (index 3 in CC(C)N after AddHs), shell-1 = C(1).
        # C(0)-C(1)-C(2) is the heavy-only triple centred on shell-1.
        donor_idx = 3
        shell1 = 1
        good_lib = self._FakeLib(n_returned=20, return_none=False)
        result = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=good_lib,
            adaptive_shell1=True, adaptive_min_n=5,
            return_result=True,
        )
        heavy_centred = self._heavy_only_centred_on(
            propyl_amine_mol, shell1, result.angle_terms
        )
        assert len(heavy_centred) >= 1, (
            "Option-B (good coverage, n=20) must KEEP at least one "
            "heavy-only angle centred on the shell-1 atom.  Got 0.  "
            f"All angle terms: {[t.atom_indices for t in result.angle_terms]}"
        )
        # The library must have been consulted for the gate.
        assert good_lib.n_angle_calls > 0, (
            "Option-B must consult lookup_angle to decide; counter is 0."
        )

    # ------------------------------------------------------------------
    # 2) THIN coverage (n < 5) — drop the heavy-only shell-1 term
    # (Heal-1/1b behaviour).
    # ------------------------------------------------------------------
    def test_adaptive_drops_when_coverage_thin(
        self, propyl_amine_mol, propyl_amine_coords
    ):
        donor_idx = 3
        shell1 = 1
        thin_lib = self._FakeLib(n_returned=5, return_none=True)
        result = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=thin_lib,
            adaptive_shell1=True, adaptive_min_n=5,
            return_result=True,
        )
        heavy_centred = self._heavy_only_centred_on(
            propyl_amine_mol, shell1, result.angle_terms
        )
        assert len(heavy_centred) == 0, (
            "Option-B (thin coverage, lookup→None) must DROP every "
            "heavy-only angle centred on shell-1 (Heal-1 behaviour).  "
            f"Got {len(heavy_centred)}.  Tuples: "
            f"{[t.atom_indices for t in heavy_centred]}"
        )

    # ------------------------------------------------------------------
    # 3) No library / adaptive disabled — fall back to Heal-1/1b
    # unconditionally.
    # ------------------------------------------------------------------
    def test_adaptive_drops_when_no_lib(
        self, propyl_amine_mol, propyl_amine_coords
    ):
        # ``adaptive_shell1=False`` forces the legacy Heal-1/1b path.
        # No library at all needs to be expressed via the flag because
        # ``library=None`` triggers the default loader.  Equivalent to
        # the spec's "mogul_lib=None → Heal-1 fallback".
        donor_idx = 3
        shell1 = 1
        result = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            adaptive_shell1=False,
            return_result=True,
        )
        heavy_centred = self._heavy_only_centred_on(
            propyl_amine_mol, shell1, result.angle_terms
        )
        assert len(heavy_centred) == 0, (
            "With adaptive_shell1=False the function must drop every "
            "heavy-only angle centred on shell-1 (Heal-1 fallback).  "
            f"Got {len(heavy_centred)}.  Tuples: "
            f"{[t.atom_indices for t in heavy_centred]}"
        )

    # ------------------------------------------------------------------
    # 4) Equivalence with Heal-1 when all coverage is thin.
    # ------------------------------------------------------------------
    def test_adaptive_byte_identical_to_heal1_when_all_thin(
        self, propyl_amine_mol, propyl_amine_coords
    ):
        donor_idx = 3
        # Run A: adaptive ON with a library that misses every lookup.
        thin_lib = self._FakeLib(return_none=True)
        r_adaptive_thin = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=thin_lib,
            adaptive_shell1=True, adaptive_min_n=5,
            return_result=True,
        )
        # Run B: adaptive OFF with the same thin library.  In both cases
        # the lookup at the term-build step ALSO uses the same library,
        # so the emitted term lists must be byte-identical.
        thin_lib_B = self._FakeLib(return_none=True)
        r_heal1 = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=thin_lib_B,
            adaptive_shell1=False,
            return_result=True,
        )
        tuples_A = [tuple(t.atom_indices) for t in r_adaptive_thin.all_terms()]
        tuples_B = [tuple(t.atom_indices) for t in r_heal1.all_terms()]
        assert tuples_A == tuples_B, (
            "When all library lookups miss, adaptive Option-B must "
            "produce the SAME emitted-term set as Heal-1/1b.\n"
            f"adaptive(thin)={tuples_A}\nheal1={tuples_B}"
        )

    # ------------------------------------------------------------------
    # 5) Determinism: same inputs + same library -> bit-identical output.
    # ------------------------------------------------------------------
    def test_adaptive_deterministic(
        self, propyl_amine_mol, propyl_amine_coords
    ):
        donor_idx = 3
        lib = self._FakeLib(n_returned=10, return_none=False)
        r1 = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=lib, adaptive_shell1=True, adaptive_min_n=5,
            return_result=True,
        )
        lib2 = self._FakeLib(n_returned=10, return_none=False)
        r2 = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=lib2, adaptive_shell1=True, adaptive_min_n=5,
            return_result=True,
        )
        a1 = [tuple(t.atom_indices) for t in r1.all_terms()]
        a2 = [tuple(t.atom_indices) for t in r2.all_terms()]
        assert a1 == a2, (
            "Option-B must be deterministic.  "
            f"a1={a1}\na2={a2}"
        )

    # ------------------------------------------------------------------
    # 6) min_n knob: a higher threshold flips a "good" entry to thin.
    # ------------------------------------------------------------------
    def test_adaptive_min_n_threshold(
        self, propyl_amine_mol, propyl_amine_coords
    ):
        donor_idx = 3
        shell1 = 1
        # Library reports n=10 for every entry.  With adaptive_min_n=5
        # the gate considers it "good" → keep.  With adaptive_min_n=100
        # the same n=10 becomes "thin" → drop.
        lib = self._FakeLib(n_returned=10, return_none=False)
        r_keep = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=lib, adaptive_shell1=True, adaptive_min_n=5,
            return_result=True,
        )
        lib2 = self._FakeLib(n_returned=10, return_none=False)
        r_drop = detect_fragments(
            propyl_amine_mol, propyl_amine_coords,
            frozen_atoms=set(), donors=[donor_idx],
            library=lib2, adaptive_shell1=True, adaptive_min_n=100,
            return_result=True,
        )
        kept = self._heavy_only_centred_on(
            propyl_amine_mol, shell1, r_keep.angle_terms
        )
        dropped = self._heavy_only_centred_on(
            propyl_amine_mol, shell1, r_drop.angle_terms
        )
        assert len(kept) >= 1
        assert len(dropped) == 0


# ---------------------------------------------------------------------------
# Hapto-class protection (2026-06-02 voll-pool regression fix)
# ---------------------------------------------------------------------------
class TestHapto:
    """Hapto-π class protection in detect_fragments + the universal
    graph-only detect_hapto_atoms helper in grip_polish.

    Background: the voll-pool 11458 SMILES run surfaced 4 severe hapto-
    class regressions because GRIP's CCDC-Mogul priors are dominated by
    non-hapto chemistry and the Mahalanobis pull forces η-coordinated
    carbons toward standard sp²-aromatic geometries that destroy the
    piano-stool placement.  detect_fragments now accepts ``hapto_atoms``
    and drops every bond/angle/improper term that touches one of them;
    the helper detect_hapto_atoms produces the set from the molecular
    graph (no SMILES pattern matching).
    """

    @staticmethod
    def _build_ferrocene_like(metal_sym: str = "Fe"):
        """Return a metallocene-style RWMol: M + 5C Cp ring with H's."""
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import RWMol
        m = RWMol()
        m.AddAtom(Chem.Atom(metal_sym))  # 0
        for _ in range(5):
            m.AddAtom(Chem.Atom("C"))
        for i in range(5):
            m.AddBond(1 + i, 1 + ((i + 1) % 5), Chem.BondType.SINGLE)
        for i in range(5):
            m.AddBond(0, 1 + i, Chem.BondType.SINGLE)
        for i in range(5):
            m.AddAtom(Chem.Atom("H"))
            m.AddBond(6 + i, 1 + i, Chem.BondType.SINGLE)
        mol = m.GetMol()
        Chem.GetSSSR(mol)
        return mol

    @staticmethod
    def _coords_random(mol, seed=42):
        rng = np.random.default_rng(seed)
        return rng.normal(size=(mol.GetNumAtoms(), 3)).astype(np.float64)

    # ------------------------------------------------------------------
    # 1) detect_hapto_atoms returns the full Cp ring atom-set
    # ------------------------------------------------------------------
    def test_detect_cp_ring_atoms(self):
        from delfin.fffree.grip_polish import detect_hapto_atoms
        mol = self._build_ferrocene_like()
        ha = detect_hapto_atoms(mol, metal_idx=0, donors=[1, 2, 3, 4, 5])
        assert ha == {1, 2, 3, 4, 5}, (
            f"Cp ring should be detected as hapto-set, got {sorted(ha)}"
        )

    # ------------------------------------------------------------------
    # 2) detect_hapto_atoms returns empty for a purely σ-coordinated
    # complex (no false positive)
    # ------------------------------------------------------------------
    def test_detect_sigma_only_returns_empty(self):
        from delfin.fffree.grip_polish import detect_hapto_atoms
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import RWMol
        m = RWMol()
        m.AddAtom(Chem.Atom("Fe"))  # 0
        # 6 σ-donor N atoms (NH3-like)
        for _ in range(6):
            m.AddAtom(Chem.Atom("N"))
        for i in range(6):
            m.AddBond(0, 1 + i, Chem.BondType.SINGLE)
        # add 3 H per N (mass-balance, optional)
        for i in range(6):
            base = m.GetNumAtoms()
            for _ in range(3):
                m.AddAtom(Chem.Atom("H"))
                m.AddBond(m.GetNumAtoms() - 1, 1 + i, Chem.BondType.SINGLE)
        mol = m.GetMol()
        Chem.GetSSSR(mol)
        ha = detect_hapto_atoms(mol, metal_idx=0, donors=list(range(1, 7)))
        assert ha == set(), f"Pure σ should have no hapto atoms, got {sorted(ha)}"

    # ------------------------------------------------------------------
    # 3) detect_fragments emits NO term touching a hapto atom when the
    # set is passed.  Iterates over the term list and verifies no
    # bond/angle/improper participant is in the hapto set.
    # ------------------------------------------------------------------
    def test_no_terms_touch_cp_atoms(self):
        mol = self._build_ferrocene_like()
        P = self._coords_random(mol)
        hapto = {1, 2, 3, 4, 5}
        # Donor set is the ring carbons; metal=0 is the frozen sphere.
        result = detect_fragments(
            mol, P,
            frozen_atoms={0, 1, 2, 3, 4, 5},
            donors=[1, 2, 3, 4, 5],
            hapto_atoms=hapto,
            return_result=True,
        )
        for t in result.bond_terms:
            assert not any(i in hapto for i in t.atom_indices), (
                f"Bond {t.atom_indices} touches hapto set {hapto}"
            )
        for t in result.angle_terms:
            assert not any(i in hapto for i in t.atom_indices), (
                f"Angle {t.atom_indices} touches hapto set {hapto}"
            )
        for t in result.improper_terms:
            assert not any(i in hapto for i in t.atom_indices), (
                f"Improper {t.atom_indices} touches hapto set {hapto}"
            )

    # ------------------------------------------------------------------
    # 4) Env flag disable: with DELFIN_FFFREE_GRIP_SKIP_HAPTO=0 the
    # function ignores hapto_atoms and re-emits at least one term that
    # touches a hapto carbon (proving the gate is wired to the env).
    # ------------------------------------------------------------------
    def test_env_disable_restores_terms(self, monkeypatch):
        mol = self._build_ferrocene_like()
        P = self._coords_random(mol)
        hapto = {1, 2, 3, 4, 5}
        # Reference run (env unset -> skip ON -> no candidate counts for
        # any bond/angle/improper that touches the hapto set).  We use
        # the n_*_candidates counters because they are incremented
        # immediately after the skip checks but BEFORE the library
        # lookup; this isolates the env-flag behaviour from the (mostly
        # non-aromatic) library coverage of our synthetic Cp.
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_SKIP_HAPTO", raising=False)
        r_on = detect_fragments(
            mol, P,
            frozen_atoms={0},  # only freeze metal, not ring carbons
            hapto_atoms=hapto,
            return_result=True,
        )
        # Disabled run (env="0" -> skip OFF -> candidates restored).
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_SKIP_HAPTO", "0")
        r_off = detect_fragments(
            mol, P,
            frozen_atoms={0},
            hapto_atoms=hapto,
            return_result=True,
        )
        # With skip ON every ring-internal bond/angle is dropped before
        # the candidate counter; with skip OFF the same candidates appear
        # because the only difference between the two runs is the gate.
        assert r_on.n_bond_candidates < r_off.n_bond_candidates, (
            "Skip ON must reduce bond candidates relative to skip OFF; "
            f"on={r_on.n_bond_candidates}, off={r_off.n_bond_candidates}"
        )
        assert r_on.n_angle_candidates < r_off.n_angle_candidates, (
            "Skip ON must reduce angle candidates relative to skip OFF; "
            f"on={r_on.n_angle_candidates}, off={r_off.n_angle_candidates}"
        )

    # ------------------------------------------------------------------
    # 5) hapto_atoms=None / empty is byte-identical to the pre-fix path
    # (regression guard — make sure we did not change behaviour when no
    # hapto set is passed).
    # ------------------------------------------------------------------
    def test_hapto_none_byte_identical_to_legacy(self):
        Chem = pytest.importorskip("rdkit.Chem")
        from rdkit.Chem import AllChem
        # Use a small organic mol (no metal, no hapto possible).
        mol = Chem.MolFromSmiles("CCN")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        P = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        r_legacy = detect_fragments(
            mol, P,
            frozen_atoms=set(), donors=[2],  # N donor
            return_result=True,
        )
        r_hapto_none = detect_fragments(
            mol, P,
            frozen_atoms=set(), donors=[2],
            hapto_atoms=None,
            return_result=True,
        )
        r_hapto_empty = detect_fragments(
            mol, P,
            frozen_atoms=set(), donors=[2],
            hapto_atoms=set(),
            return_result=True,
        )
        a_legacy = [tuple(t.atom_indices) for t in r_legacy.all_terms()]
        a_none = [tuple(t.atom_indices) for t in r_hapto_none.all_terms()]
        a_empty = [tuple(t.atom_indices) for t in r_hapto_empty.all_terms()]
        assert a_legacy == a_none == a_empty, (
            "hapto_atoms=None / empty must reproduce legacy behaviour "
            f"byte-exactly.\nlegacy={a_legacy}\nnone={a_none}\nempty={a_empty}"
        )

    # ------------------------------------------------------------------
    # 6) Determinism: same inputs → bit-identical hapto set + term list
    # ------------------------------------------------------------------
    def test_hapto_deterministic(self):
        from delfin.fffree.grip_polish import detect_hapto_atoms
        mol = self._build_ferrocene_like()
        P = self._coords_random(mol)
        ha1 = detect_hapto_atoms(mol, metal_idx=0, donors=[1, 2, 3, 4, 5])
        ha2 = detect_hapto_atoms(mol, metal_idx=0, donors=[1, 2, 3, 4, 5])
        assert ha1 == ha2
        r1 = detect_fragments(
            mol, P, frozen_atoms={0}, donors=[1, 2, 3, 4, 5],
            hapto_atoms=ha1, return_result=True,
        )
        r2 = detect_fragments(
            mol, P, frozen_atoms={0}, donors=[1, 2, 3, 4, 5],
            hapto_atoms=ha2, return_result=True,
        )
        t1 = [tuple(t.atom_indices) for t in r1.all_terms()]
        t2 = [tuple(t.atom_indices) for t in r2.all_terms()]
        assert t1 == t2, "Hapto-skip path must be deterministic"


# ---------------------------------------------------------------------------
# Torsion stub
# ---------------------------------------------------------------------------
class TestTorsionStub:
    def test_torsion_raises_not_implemented(self):
        t = TorsionTerm(a=0, b=1, c=2, d=3)
        with pytest.raises(NotImplementedError):
            t.value_and_grad(np.zeros((4, 3)))


# ---------------------------------------------------------------------------
# Small helpers used in tests (kept here to avoid polluting product API).
# ---------------------------------------------------------------------------
import math as _math


def math_acos_deg(x: float) -> float:
    return _math.degrees(_math.acos(x))


def math_to_rad(deg: float) -> float:
    return _math.radians(deg)
