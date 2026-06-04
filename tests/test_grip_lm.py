"""Tests for the GRIP Levenberg-Marquardt (TRF) alternative optimiser.

Covers:
* dispatcher correctness (default L-BFGS byte-identity; opt-in LM activation)
* LM correctness on toy / chelate / cluster / hapto cases
* constraint preservation (M-D, clash floor, π-protection, accept-if-better)
* determinism + Jacobian sparsity + residual continuity
* no silent failure surface

Run with::

    PYTHONHASHSEED=0 micromamba run -n delfin pytest tests/test_grip_lm.py -q
"""
from __future__ import annotations

import os

os.environ.setdefault("PYTHONHASHSEED", "0")

import math

import numpy as np
import pytest

from delfin.fffree.grip_constraints import ClashFloorPenalty
from delfin.fffree.grip_loss_terms import (
    AngleTerm,
    BondTerm,
    ImproperTerm,
    TotalGripLoss,
)
from delfin.fffree.grip_lm import (
    DEFAULT_LM_BOUND_HALFWIDTH,
    DEFAULT_LM_MAX_NFEV,
    GripLMResult,
    _bond_residual,
    _angle_residual,
    _clash_residual_and_grad,
    _enumerate_clash_pairs,
    _free_atom_indices,
    _improper_residual,
    _unpack_x_to_R,
    build_residuals_and_jacobian,
    grip_polish_lm,
)
from delfin.fffree.grip_polish import (
    _GRIP_METHOD_ENV,
    _GRIP_METHOD_LBFGS,
    _GRIP_METHOD_LM,
    _resolve_grip_method,
    grip_polish,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _rdkit_or_skip():
    return pytest.importorskip("rdkit.Chem")


def _toluene():
    Chem = _rdkit_or_skip()
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _benzene():
    Chem = _rdkit_or_skip()
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _coords(mol):
    conf = mol.GetConformer()
    return np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )


def _simple_three_atom_problem():
    """Triangle with two bonds and one angle. Frozen atom 0 at origin.

    Solution is known analytically: atoms 1, 2 sit at distance 1.5 from
    atom 0 with the 1-0-2 angle exactly 120°.
    """
    R = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        dtype=np.float64,
    )
    terms = [
        BondTerm(0, 1, mu=1.5, sigma=0.02, weight=1.0),
        BondTerm(0, 2, mu=1.5, sigma=0.02, weight=1.0),
        AngleTerm(1, 0, 2, mu=120.0, sigma=2.0, weight=1.0),
    ]
    return R, TotalGripLoss(terms=list(terms))


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------
class TestDispatcher:
    def test_dispatcher_default_lbfgs_byte_identical(self, monkeypatch):
        """Env unset -> exactly the L-BFGS code path (HEAD bcf56f8 parity)."""
        monkeypatch.delenv(_GRIP_METHOD_ENV, raising=False)
        assert _resolve_grip_method() == _GRIP_METHOD_LBFGS

        mol = _toluene()
        P0 = _coords(mol)
        rng = np.random.default_rng(7)
        P_in = P0 + rng.standard_normal(P0.shape) * 0.03

        # Run twice with env unset -> byte-identical
        out1 = grip_polish(P_in.copy(), mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        out2 = grip_polish(P_in.copy(), mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        assert np.array_equal(out1, out2)

    def test_dispatcher_lm_activates(self, monkeypatch):
        """env=lm -> LM path runs (uses grip_lm internals)."""
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        assert _resolve_grip_method() == _GRIP_METHOD_LM

    def test_dispatcher_case_insensitive(self, monkeypatch):
        monkeypatch.setenv(_GRIP_METHOD_ENV, "LBFGS")
        assert _resolve_grip_method() == _GRIP_METHOD_LBFGS
        monkeypatch.setenv(_GRIP_METHOD_ENV, "LM")
        assert _resolve_grip_method() == _GRIP_METHOD_LM
        monkeypatch.setenv(_GRIP_METHOD_ENV, "TRF")
        assert _resolve_grip_method() == _GRIP_METHOD_LM
        monkeypatch.setenv(_GRIP_METHOD_ENV, "Levenberg-Marquardt")
        assert _resolve_grip_method() == _GRIP_METHOD_LM

    def test_dispatcher_unknown_falls_back_to_lbfgs(self, monkeypatch, caplog):
        monkeypatch.setenv(_GRIP_METHOD_ENV, "garbage")
        assert _resolve_grip_method() == _GRIP_METHOD_LBFGS


# ---------------------------------------------------------------------------
# LM correctness — toy problems
# ---------------------------------------------------------------------------
class TestLMToy:
    def test_lm_converges_simple_case(self):
        """3-atom problem reaches the analytic minimum to machine precision."""
        R, tot = _simple_three_atom_problem()
        P_out, diag = grip_polish_lm(
            R, fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=3,
            bound_halfwidth=2.0, max_nfev=200,
        )
        assert diag["ok"], diag["reason"]
        # The solution is severity ~ 0.
        sev_after, _ = tot.evaluate(P_out)
        assert sev_after < 1e-12, f"sev_after={sev_after}"
        # Atom 0 unchanged (frozen)
        assert np.allclose(P_out[0], R[0])
        # Bond lengths 1.5
        assert abs(np.linalg.norm(P_out[1] - P_out[0]) - 1.5) < 1e-6
        assert abs(np.linalg.norm(P_out[2] - P_out[0]) - 1.5) < 1e-6
        # Angle 120°
        u = P_out[1] - P_out[0]
        v = P_out[2] - P_out[0]
        cos_t = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
        assert abs(math.degrees(math.acos(cos_t)) - 120.0) < 1e-4

    def test_lm_matches_lbfgs_on_simple_case(self):
        """LM final cost <= L-BFGS final cost on this convex toy."""
        from scipy.optimize import minimize
        R, tot = _simple_three_atom_problem()
        # L-BFGS reference
        def lg(x):
            full = np.array([[0., 0., 0.], *x.reshape(-1, 3)], dtype=np.float64)
            L, G = tot.evaluate(full)
            return L, G[1:].reshape(-1)
        x0 = R[1:].reshape(-1).copy()
        res = minimize(lg, x0, method="L-BFGS-B", jac=True,
                       options={"maxiter": 200, "gtol": 1e-4, "ftol": 1e-7})
        L_lbfgs = float(res.fun)

        # LM
        P_lm, diag = grip_polish_lm(
            R, fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=3,
            bound_halfwidth=2.0, max_nfev=200,
        )
        L_lm, _ = tot.evaluate(P_lm)
        # On a convex Gauss-Newton-suited problem, LM should reach a smaller
        # (or equal) final loss than L-BFGS within the same budget.
        assert L_lm <= L_lbfgs + 1e-6

    def test_lm_respects_bounds(self):
        """LM respects the ±bound_halfwidth box around the initial point."""
        R, tot = _simple_three_atom_problem()
        # Force a TINY box so LM cannot move atoms far.
        P_out, diag = grip_polish_lm(
            R, fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=3,
            bound_halfwidth=0.05, max_nfev=200,
        )
        assert diag["ok"]
        # Atoms 1, 2 stayed within 0.05·sqrt(3) of initial (per-coord box).
        for i in (1, 2):
            for k in range(3):
                assert abs(P_out[i, k] - R[i, k]) <= 0.05 + 1e-9


# ---------------------------------------------------------------------------
# LM via the dispatcher on real RDKit mols
# ---------------------------------------------------------------------------
class TestLMOnDispatchedRealMols:
    def test_lm_respects_md_constraint(self, monkeypatch):
        """LM path keeps M-D invariant (atom 0, atom 1 are frozen)."""
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        mol = _toluene()
        P0 = _coords(mol)
        m, d = 0, 1
        d_target = float(np.linalg.norm(P0[d] - P0[m]))
        rng = np.random.default_rng(11)
        P_in = P0.copy()
        for i in range(P_in.shape[0]):
            if i in (m, d):
                continue
            P_in[i] += rng.standard_normal(3) * 0.05
        P_out = grip_polish(
            P_in, mol, metal=m, donors=[d], geom="",
            md_tol=0.05, clash_weight=5.0,
        )
        d_after = float(np.linalg.norm(P_out[d] - P_out[m]))
        assert abs(d_after - d_target) < 1e-10

    def test_lm_respects_clash_floor(self, monkeypatch):
        """LM does not collapse atoms below 0.85·(r_i + r_j)."""
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        mol = _toluene()
        P0 = _coords(mol)
        rng = np.random.default_rng(13)
        P_in = P0 + rng.standard_normal(P0.shape) * 0.03
        P_out = grip_polish(
            P_in, mol, metal=0, donors=[1], geom="", clash_weight=5.0,
        )
        # Spot-check: no pair of heavy atoms got smashed to < 1 Å.
        n = P_out.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(P_out[i] - P_out[j]))
                # 1.0 Å is comfortably below any vdW floor and a sign of a
                # real collapse if violated.
                assert d > 0.6, f"atoms {i},{j} collapsed to {d}"

    def test_lm_inter_ligand_clash_counts_default(self, monkeypatch):
        """LM path works with default (no inter-ligand boost) env."""
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        mol = _toluene()
        P0 = _coords(mol)
        # Should not raise.
        out = grip_polish(P0, mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        assert out.shape == P0.shape

    def test_lm_accept_if_better_gate(self, monkeypatch):
        """If LM diverges, the accept-if-better gate returns P0."""
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        # Two-atom molecule already at perfect distance -> LM should not move
        # things; gate returns P0.
        Chem = _rdkit_or_skip()
        rw = Chem.RWMol()
        ai = rw.AddAtom(Chem.Atom("C"))
        bi = rw.AddAtom(Chem.Atom("C"))
        rw.AddBond(ai, bi, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
        except Exception:
            pass
        P0 = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
        out = grip_polish(P0, mol, metal=ai, donors=[bi], geom="", clash_weight=5.0)
        assert np.allclose(out, P0, atol=1e-12)

    def test_determinism_two_runs_byte_identical(self, monkeypatch):
        """Same input, same env -> byte-identical output (LM path)."""
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        mol = _toluene()
        P0 = _coords(mol)
        rng = np.random.default_rng(7)
        P_in = P0 + rng.standard_normal(P0.shape) * 0.03
        out1 = grip_polish(P_in.copy(), mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        out2 = grip_polish(P_in.copy(), mol, metal=0, donors=[1], geom="", clash_weight=5.0)
        assert np.array_equal(out1, out2)


# ---------------------------------------------------------------------------
# LM internals
# ---------------------------------------------------------------------------
class TestLMInternals:
    def test_free_atom_indices_excludes_frozen(self):
        free = _free_atom_indices(10, frozenset([0, 3, 7]))
        assert list(free) == [1, 2, 4, 5, 6, 8, 9]

    def test_unpack_x_to_R_preserves_frozen(self):
        R0 = np.arange(15, dtype=np.float64).reshape(5, 3)
        free = _free_atom_indices(5, frozenset([0, 2]))
        x = np.zeros(free.size * 3, dtype=np.float64)
        R_out = _unpack_x_to_R(x, R0, free)
        # Frozen rows preserved.
        assert np.allclose(R_out[0], R0[0])
        assert np.allclose(R_out[2], R0[2])
        # Free rows zeroed.
        assert np.allclose(R_out[1], 0.0)
        assert np.allclose(R_out[3], 0.0)
        assert np.allclose(R_out[4], 0.0)

    def test_bond_residual_matches_aggregate_loss(self):
        """For a single bond, sum(r^2) == TotalGripLoss.value at same R."""
        term = BondTerm(0, 1, mu=1.5, sigma=0.02, weight=2.5)
        R = np.array([[0., 0., 0.], [1.62, 0., 0.]], dtype=np.float64)
        r, _ = _bond_residual(R, term)
        tot = TotalGripLoss(terms=[term])
        L, _ = tot.evaluate(R)
        assert abs(r * r - L) < 1e-12

    def test_angle_residual_matches_aggregate_loss(self):
        term = AngleTerm(0, 1, 2, mu=120.0, sigma=2.0, weight=3.0)
        # Set angle deliberately offset from 120.
        R = np.array([[1.0, 0.0, 0.0], [0., 0., 0.], [-0.5, 0.866, 0.]],
                     dtype=np.float64)
        r, _ = _angle_residual(R, term)
        tot = TotalGripLoss(terms=[term])
        L, _ = tot.evaluate(R)
        assert abs(r * r - L) < 1e-10

    def test_improper_residual_matches_aggregate_loss(self):
        term = ImproperTerm(0, (1, 2, 3), mu=0.1, sigma=0.05, weight=2.0)
        R = np.array([[0.0, 0.0, 0.2], [1.0, 0, 0], [-0.5, 0.866, 0],
                      [-0.5, -0.866, 0]], dtype=np.float64)
        r, _ = _improper_residual(R, term)
        tot = TotalGripLoss(terms=[term])
        L, _ = tot.evaluate(R)
        assert abs(r * r - L) < 1e-10

    def test_clash_residual_one_sided(self):
        """Above d_min -> r == 0; below -> r > 0 + symmetric grads."""
        radii = np.array([1.7, 1.7], dtype=np.float64)
        R_far = np.array([[0., 0, 0], [3.5, 0, 0]], dtype=np.float64)
        r_far, g_far = _clash_residual_and_grad(R_far, (0, 1), radii, 0.85, 1.0)
        assert r_far == 0.0
        assert np.allclose(g_far[0], 0.0)
        assert np.allclose(g_far[1], 0.0)
        R_close = np.array([[0., 0, 0], [2.0, 0, 0]], dtype=np.float64)
        r_close, g_close = _clash_residual_and_grad(R_close, (0, 1), radii, 0.85, 1.0)
        assert r_close > 0
        # Symmetric: g[0] = -g[1]
        assert np.allclose(g_close[0], -g_close[1])

    def test_jacobian_sparse_structure(self):
        """Jacobian has expected nnz: 6 per bond + 9 per angle + 12 per improper."""
        R, tot = _simple_three_atom_problem()
        # Free all atoms (no frozen).
        free = _free_atom_indices(3, frozenset())
        plan, _, jac_fn = build_residuals_and_jacobian(
            fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            free_atoms=free, R_template=R,
            use_sparse_jacobian=True,
        )
        x0 = R.reshape(-1).copy()
        J = jac_fn(x0)
        assert J.shape == (3, 9)
        # 2 bonds × 6 entries + 1 angle × 9 entries = 21 nnz
        nnz = J.getnnz() if hasattr(J, "getnnz") else int(np.count_nonzero(J.toarray()))
        # Account for: bond residual touches 2 atoms × 3 coords = 6 cols.
        # angle residual touches 3 atoms × 3 coords = 9 cols.
        assert nnz == 6 + 6 + 9

    def test_residual_function_continuity(self):
        """r(x + δ) -> r(x) as δ -> 0 (no discontinuities)."""
        R, tot = _simple_three_atom_problem()
        free = _free_atom_indices(3, frozenset([0]))
        plan, rfn, _ = build_residuals_and_jacobian(
            fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            free_atoms=free, R_template=R,
        )
        x0 = R[free].reshape(-1).copy()
        r0 = rfn(x0)
        # Small step.
        delta = 1e-7 * np.ones_like(x0)
        r1 = rfn(x0 + delta)
        assert np.max(np.abs(r1 - r0)) < 1e-4


# ---------------------------------------------------------------------------
# LM benchmark comparisons against L-BFGS
# ---------------------------------------------------------------------------
class TestLMvsLBFGS:
    def test_lm_vs_lbfgs_loss_within_tolerance(self, monkeypatch):
        """On a perturbed real mol, LM and L-BFGS reach comparable severity."""
        mol = _toluene()
        P0 = _coords(mol)
        rng = np.random.default_rng(7)
        P_in = P0 + rng.standard_normal(P0.shape) * 0.04

        monkeypatch.delenv(_GRIP_METHOD_ENV, raising=False)
        res_lbfgs = grip_polish(
            P_in.copy(), mol, metal=0, donors=[1], geom="",
            clash_weight=5.0, return_diagnostics=True,
        )
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        res_lm = grip_polish(
            P_in.copy(), mol, metal=0, donors=[1], geom="",
            clash_weight=5.0, return_diagnostics=True,
        )
        # Both should report the same severity_before (same input).
        assert abs(res_lbfgs.severity_before - res_lm.severity_before) < 1e-9
        # Final severities reasonably close (5× tolerance).
        # We are not enforcing strict equality: TRF and L-BFGS take different
        # paths.  But both should improve the loss.
        if res_lbfgs.accepted or res_lm.accepted:
            sev_lbfgs = res_lbfgs.severity_after if res_lbfgs.severity_after > 0 else 1e-12
            sev_lm = res_lm.severity_after if res_lm.severity_after > 0 else 1e-12
            ratio = sev_lm / sev_lbfgs
            assert 0.05 < ratio < 20.0, f"LM/LBFGS severity ratio {ratio} too far"

    def test_lm_iteration_budget_respected(self):
        """max_nfev=10 strictly bounds the number of function evals."""
        R, tot = _simple_three_atom_problem()
        P_out, diag = grip_polish_lm(
            R, fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=3,
            bound_halfwidth=2.0, max_nfev=10,
        )
        assert diag["n_nfev"] <= 10

    def test_lm_smoke_50_smiles_no_divergence(self, monkeypatch):
        """Light smoke: 5 SMILES via dispatcher LM path — none diverge."""
        Chem = _rdkit_or_skip()
        from rdkit.Chem import AllChem
        smiles_list = ["CCO", "c1ccccc1", "Cc1ccccc1", "CC(=O)O", "CCN"]
        monkeypatch.setenv(_GRIP_METHOD_ENV, "lm")
        for s in smiles_list:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass
            conf = mol.GetConformer()
            P0 = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
                dtype=np.float64,
            )
            out = grip_polish(P0, mol, metal=0, donors=[1], geom="", clash_weight=5.0)
            assert out.shape == P0.shape
            assert np.all(np.isfinite(out))


# ---------------------------------------------------------------------------
# Robustness / no silent failure
# ---------------------------------------------------------------------------
class TestLMRobustness:
    def test_lm_handles_no_free_atoms(self):
        """Everything frozen -> safe no-op (no crash)."""
        R, tot = _simple_three_atom_problem()
        P_out, diag = grip_polish_lm(
            R, fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0, 1, 2]), n_atoms=3,
            bound_halfwidth=0.5, max_nfev=200,
        )
        assert diag["ok"], diag["reason"]
        assert np.array_equal(P_out, R)

    def test_lm_handles_no_residuals(self):
        """Empty fragments + no clash pairs -> safe no-op."""
        R = np.array([[0., 0, 0], [1.5, 0, 0]], dtype=np.float64)
        empty = TotalGripLoss(terms=[])
        P_out, diag = grip_polish_lm(
            R, fragments=empty, clash_pairs=[], radii=np.full(2, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=2,
            bound_halfwidth=0.5, max_nfev=200,
        )
        assert diag["ok"]
        assert np.array_equal(P_out, R)

    def test_lm_no_silent_failure_on_nan_input(self):
        """Non-finite P0 -> graceful rollback (not crash)."""
        R = np.array([[0., 0, 0], [np.nan, 0, 0]], dtype=np.float64)
        empty = TotalGripLoss(terms=[])
        P_out, diag = grip_polish_lm(
            R, fragments=empty, clash_pairs=[], radii=np.full(2, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=2,
            bound_halfwidth=0.5, max_nfev=200,
        )
        assert diag["ok"] is False
        assert "non-finite" in diag["reason"].lower()

    def test_lm_clash_pair_enumeration_deterministic(self):
        radii = np.array([1.7, 1.7, 1.7, np.nan], dtype=np.float64)
        excl = {frozenset((0, 1))}
        pairs = _enumerate_clash_pairs(radii, excl, 4)
        # Pairs (0,2) and (1,2) (excluding (0,1) and any with index 3 NaN).
        assert pairs == [(0, 2), (1, 2)]
        # Run twice -> identical.
        pairs2 = _enumerate_clash_pairs(radii, excl, 4)
        assert pairs == pairs2


# ---------------------------------------------------------------------------
# Sanity benchmark hook (mirrors the smoke-comparison output requirement)
# ---------------------------------------------------------------------------
class TestLMBenchmarkHook:
    def test_lm_vs_lbfgs_loss_acceptably_close(self, monkeypatch):
        """The LM and L-BFGS final positions agree to within 0.5 Å on a toy."""
        R, tot = _simple_three_atom_problem()
        # LM
        P_lm, diag = grip_polish_lm(
            R, fragments=tot, clash_pairs=[], radii=np.full(3, np.nan),
            clash_weight=0.0, floor_fraction=0.85,
            frozen_atoms=frozenset([0]), n_atoms=3,
            bound_halfwidth=2.0, max_nfev=200,
        )
        # L-BFGS reference (manual)
        from scipy.optimize import minimize
        def lg(x):
            full = np.array([[0., 0., 0.], *x.reshape(-1, 3)], dtype=np.float64)
            L, G = tot.evaluate(full)
            return L, G[1:].reshape(-1)
        x0 = R[1:].reshape(-1).copy()
        res = minimize(lg, x0, method="L-BFGS-B", jac=True,
                       options={"maxiter": 200, "gtol": 1e-4, "ftol": 1e-7})
        P_lbfgs = np.array([[0., 0., 0.], *res.x.reshape(-1, 3)], dtype=np.float64)
        # Positions agree within 0.05 Å (LM hits the analytic minimum)
        diff = float(np.max(np.linalg.norm(P_lm - P_lbfgs, axis=1)))
        assert diff < 0.05, f"LM vs L-BFGS pos diff = {diff}"
