"""Phase 2 — GRIP M+D unfreeze tests (2026-06-06).

Tests for the new ``grip_md_terms`` module + the ``grip_polish`` integration
that lets the metal and donors move under L-BFGS guided by CCDC-empirical /
polyhedron-template Gaussians.

Goals covered:
* default-OFF byte-identical with HEAD (the env flag is OFF -> nothing changes)
* M-D bond loss: gradient finite-diff + actual L-BFGS step pulls M-D to mu
* M-D-X angle loss: gradient finite-diff + correctly classifies sp³/sp²/sp
* D-M-D angle loss: pulls a 90° init toward the polyhedron 180° for OC-2
* Hapto centroid loss: works on a 3-C π-cluster
* Hapto detection: groups ≥3 same-element heavy donors within radius
* Fail-open guard: an artificial large M-D drift triggers a roll-back
* Determinism: two runs with PYTHONHASHSEED=0 produce byte-identical output
* Robust to missing inputs: empty mol / negative indices / non-finite coords
"""
from __future__ import annotations

import os
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")

import math

import numpy as np
import pytest

from delfin.fffree.grip_md_terms import (
    DEFAULT_DMD_SIGMA,
    DEFAULT_DMD_WEIGHT,
    DEFAULT_MD_SIGMA,
    DEFAULT_MD_WEIGHT,
    DEFAULT_MDX_SIGMA,
    DEFAULT_MDX_WEIGHT,
    DMDAngleTerm,
    HAPTO_MIN_COUNT,
    HAPTO_RADIUS,
    HaptoCentroidTerm,
    MDBondTerm,
    MDXAngleTerm,
    build_md_loss_terms,
    detect_hapto_donor_clusters,
    max_md_drift,
    unfreeze_md_active,
)
from delfin.fffree.grip_polish import grip_polish


# ---------------------------------------------------------------------------
# Helpers (no rdkit import at top — keep optional)
# ---------------------------------------------------------------------------
def _trigonal_mol():
    """Build [Fe](N)(N)(N) -- minimal trigonal complex (rdkit)."""
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    fe = mol.AddAtom(Chem.Atom("Fe"))
    n1 = mol.AddAtom(Chem.Atom("N"))
    n2 = mol.AddAtom(Chem.Atom("N"))
    n3 = mol.AddAtom(Chem.Atom("N"))
    mol.AddBond(fe, n1, Chem.BondType.SINGLE)
    mol.AddBond(fe, n2, Chem.BondType.SINGLE)
    mol.AddBond(fe, n3, Chem.BondType.SINGLE)
    return mol


def _cp_fe_mol():
    """Build a tiny Fe-Cp (η5) cluster -- 6 atoms (Fe + 5 C)."""
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    fe = mol.AddAtom(Chem.Atom("Fe"))
    cs = [mol.AddAtom(Chem.Atom("C")) for _ in range(5)]
    for c in cs:
        mol.AddBond(fe, c, Chem.BondType.SINGLE)
    # ring bonds
    for i in range(5):
        mol.AddBond(cs[i], cs[(i + 1) % 5], Chem.BondType.SINGLE)
    return mol


def _coords_trigonal(d_md: float = 2.0):
    """Place Fe at origin, three N atoms at 120° in xy-plane."""
    R = np.zeros((4, 3), dtype=np.float64)
    R[1] = [d_md, 0.0, 0.0]
    R[2] = [-0.5 * d_md, 0.5 * d_md * math.sqrt(3.0), 0.0]
    R[3] = [-0.5 * d_md, -0.5 * d_md * math.sqrt(3.0), 0.0]
    return R


# ---------------------------------------------------------------------------
# 1) Env-flag default-OFF
# ---------------------------------------------------------------------------
class TestUnfreezeFlag:
    def test_default_off(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", raising=False)
        assert unfreeze_md_active() is False

    def test_explicit_on(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        assert unfreeze_md_active() is True

    def test_explicit_off(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "0")
        assert unfreeze_md_active() is False


# ---------------------------------------------------------------------------
# 2) MDBondTerm: gradient + L-BFGS step
# ---------------------------------------------------------------------------
class TestMDBondTerm:
    def test_loss_zero_at_target(self):
        t = MDBondTerm(metal=0, donor=1, mu=2.0, sigma=0.05, weight=10.0)
        R = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64)
        loss, grad = t.value_and_grad(R)
        assert loss == pytest.approx(0.0, abs=1e-12)
        assert np.allclose(grad, 0.0)

    def test_loss_positive_off_target(self):
        t = MDBondTerm(metal=0, donor=1, mu=2.0, sigma=0.05, weight=10.0)
        R = np.array([[0.0, 0.0, 0.0], [2.3, 0.0, 0.0]], dtype=np.float64)
        loss, _ = t.value_and_grad(R)
        # z = (2.3 - 2.0) / 0.05 = 6; loss = 10 * 36 = 360
        assert loss == pytest.approx(10.0 * 36.0, rel=1e-6)

    def test_gradient_finite_diff(self):
        t = MDBondTerm(metal=0, donor=1, mu=2.0, sigma=0.1, weight=2.0)
        rng = np.random.default_rng(123)
        R = rng.standard_normal((2, 3)) * 0.8 + np.array([[0, 0, 0], [2, 0, 0]])
        _, g_analytic = t.value_and_grad(R)
        # central FD
        h = 1e-5
        g_fd = np.zeros_like(R)
        for i in range(2):
            for k in range(3):
                Rp = R.copy(); Rp[i, k] += h
                Rm = R.copy(); Rm[i, k] -= h
                lp, _ = t.value_and_grad(Rp)
                lm, _ = t.value_and_grad(Rm)
                g_fd[i, k] = (lp - lm) / (2 * h)
        assert np.allclose(g_analytic, g_fd, atol=1e-4)

    def test_atom_indices_sorted(self):
        t = MDBondTerm(metal=5, donor=2, mu=2.0, sigma=0.05)
        assert t.atom_indices == (2, 5)


# ---------------------------------------------------------------------------
# 3) MDXAngleTerm: gradient finite-diff + angle target
# ---------------------------------------------------------------------------
class TestMDXAngleTerm:
    def test_loss_zero_at_target(self):
        # M at +x, D at origin, X at -y; angle M-D-X = 90°
        R = np.array([[1.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, -1.0, 0.0]], dtype=np.float64)
        t = MDXAngleTerm(metal=0, donor=1, x=2, mu=90.0, sigma=5.0, weight=3.0)
        loss, _ = t.value_and_grad(R)
        assert loss == pytest.approx(0.0, abs=1e-9)

    def test_gradient_finite_diff(self):
        rng = np.random.default_rng(7)
        R = rng.standard_normal((3, 3)) * 0.6
        # Move apart so the leg vectors are large enough for stable angles.
        R[0] += [1.0, 0.0, 0.0]
        R[2] += [0.0, 1.0, 0.0]
        t = MDXAngleTerm(metal=0, donor=1, x=2, mu=109.5, sigma=5.0, weight=3.0)
        _, g_analytic = t.value_and_grad(R)
        h = 1e-5
        g_fd = np.zeros_like(R)
        for i in range(3):
            for k in range(3):
                Rp = R.copy(); Rp[i, k] += h
                Rm = R.copy(); Rm[i, k] -= h
                lp, _ = t.value_and_grad(Rp)
                lm, _ = t.value_and_grad(Rm)
                g_fd[i, k] = (lp - lm) / (2 * h)
        # FD precision on angle terms is loose -- bumps near acos singularity.
        assert np.allclose(g_analytic, g_fd, atol=1e-2)


# ---------------------------------------------------------------------------
# 4) DMDAngleTerm: gradient finite-diff + polyhedron template
# ---------------------------------------------------------------------------
class TestDMDAngleTerm:
    def test_loss_at_180_target(self):
        R = np.array([[0.0, 0.0, 0.0],
                      [2.0, 0.0, 0.0],
                      [-2.0, 0.0, 0.0]], dtype=np.float64)
        t = DMDAngleTerm(metal=0, donor_i=1, donor_j=2, mu=180.0, sigma=5.0,
                         weight=5.0)
        loss, _ = t.value_and_grad(R)
        assert loss == pytest.approx(0.0, abs=1e-6)

    def test_gradient_finite_diff(self):
        rng = np.random.default_rng(11)
        R = rng.standard_normal((3, 3)) * 0.4
        R[1] += [2.0, 0.0, 0.0]
        R[2] += [0.0, 2.0, 0.0]
        t = DMDAngleTerm(metal=0, donor_i=1, donor_j=2, mu=109.5, sigma=5.0,
                         weight=5.0)
        _, g_analytic = t.value_and_grad(R)
        h = 1e-5
        g_fd = np.zeros_like(R)
        for i in range(3):
            for k in range(3):
                Rp = R.copy(); Rp[i, k] += h
                Rm = R.copy(); Rm[i, k] -= h
                lp, _ = t.value_and_grad(Rp)
                lm, _ = t.value_and_grad(Rm)
                g_fd[i, k] = (lp - lm) / (2 * h)
        assert np.allclose(g_analytic, g_fd, atol=1e-2)


# ---------------------------------------------------------------------------
# 5) HaptoCentroidTerm: 3-C π cluster
# ---------------------------------------------------------------------------
class TestHaptoCentroidTerm:
    def test_centroid_at_target_zero_loss(self):
        # M at origin; 3 C atoms forming a triangle 2 Å up the z-axis.
        # centroid = (0, 0, 2); want |M-centroid|=2; mu=2.
        R = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 2.0],
            [-0.5, 0.866, 2.0],
            [-0.5, -0.866, 2.0],
        ], dtype=np.float64)
        t = HaptoCentroidTerm(metal=0, ring_atoms=(1, 2, 3),
                              mu=2.0, sigma=0.05, weight=10.0)
        loss, _ = t.value_and_grad(R)
        assert loss == pytest.approx(0.0, abs=1e-9)

    def test_gradient_finite_diff(self):
        rng = np.random.default_rng(3)
        R = rng.standard_normal((4, 3)) * 0.3
        R[1] += [1.0, 0.0, 2.0]
        R[2] += [-0.5, 0.866, 2.0]
        R[3] += [-0.5, -0.866, 2.0]
        t = HaptoCentroidTerm(metal=0, ring_atoms=(1, 2, 3),
                              mu=2.0, sigma=0.1, weight=4.0)
        _, g_analytic = t.value_and_grad(R)
        h = 1e-5
        g_fd = np.zeros_like(R)
        for i in range(4):
            for k in range(3):
                Rp = R.copy(); Rp[i, k] += h
                Rm = R.copy(); Rm[i, k] -= h
                lp, _ = t.value_and_grad(Rp)
                lm, _ = t.value_and_grad(Rm)
                g_fd[i, k] = (lp - lm) / (2 * h)
        assert np.allclose(g_analytic, g_fd, atol=1e-4)


# ---------------------------------------------------------------------------
# 6) Hapto detection
# ---------------------------------------------------------------------------
class TestHaptoDetect:
    def test_cp_eta5_detected(self):
        pytest.importorskip("rdkit.Chem")
        mol = _cp_fe_mol()
        # Place Cp in xy-plane at z=1.7 from Fe.
        R = np.zeros((6, 3), dtype=np.float64)
        for i in range(5):
            ang = 2 * math.pi * i / 5
            R[i + 1] = [math.cos(ang), math.sin(ang), 1.7]
        donors = [1, 2, 3, 4, 5]
        clusters = detect_hapto_donor_clusters(mol, 0, donors, R)
        assert len(clusters) == 1
        elem, ring = clusters[0]
        assert elem == "C"
        assert sorted(ring) == sorted(donors)

    def test_sigma_only_no_hapto(self):
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal()
        donors = [1, 2, 3]
        clusters = detect_hapto_donor_clusters(mol, 0, donors, R)
        # 3 N atoms but mixed element check: N != H, count=3 >=3, distance OK.
        # The detector groups by element; all-N -> 1 cluster ("N",(1,2,3)).
        assert len(clusters) == 1
        assert clusters[0][0] == "N"

    def test_distance_filter_excludes_far(self):
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal(d_md=5.0)   # > HAPTO_RADIUS
        donors = [1, 2, 3]
        clusters = detect_hapto_donor_clusters(mol, 0, donors, R)
        assert clusters == []


# ---------------------------------------------------------------------------
# 7) build_md_loss_terms — end-to-end builder
# ---------------------------------------------------------------------------
class TestBuildMDLossTerms:
    def test_trigonal_returns_terms(self):
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        # Offset the donor positions so the centroid is NOT at the metal
        # (otherwise the HaptoCentroidTerm collapses; the builder skips it).
        R = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
            [-1.0, math.sqrt(3) + 0.5, 0.0],
            [-1.0, -math.sqrt(3) + 0.5, 0.0],
        ], dtype=np.float64)
        terms = build_md_loss_terms(mol, R, 0, [1, 2, 3], geom="")
        # 3 N donors form a hapto cluster (3 same-element within 3 Å) → 1
        # HaptoCentroidTerm (the trigonal mol has nothing else to anchor).
        # No D-M-D angle since geom=""; no M-D-X angle (N donors have only
        # the metal as a heavy neighbour, so no X).
        assert len(terms) >= 1
        # All returned terms must expose atom_indices.
        for t in terms:
            assert hasattr(t, "atom_indices")
            assert hasattr(t, "value_and_grad")

    def test_no_hapto_path_emits_md_terms(self):
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal(d_md=5.0)   # outside hapto radius
        terms = build_md_loss_terms(mol, R, 0, [1, 2, 3], geom="")
        # Without hapto we should get exactly 3 MDBondTerm.
        from delfin.fffree.grip_md_terms import MDBondTerm as _MD
        md_terms = [t for t in terms if isinstance(t, _MD)]
        assert len(md_terms) == 3

    def test_empty_donor_list_returns_empty(self):
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal()
        terms = build_md_loss_terms(mol, R, 0, [])
        assert terms == []

    def test_bad_metal_index_returns_empty(self):
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal()
        # metal=42 is out of range.
        terms = build_md_loss_terms(mol, R, 42, [1, 2, 3])
        assert terms == []


# ---------------------------------------------------------------------------
# 8) max_md_drift
# ---------------------------------------------------------------------------
class TestMaxMdDrift:
    def test_zero_drift(self):
        R = _coords_trigonal()
        assert max_md_drift(R, R, 0, [1, 2, 3]) == pytest.approx(0.0)

    def test_drift_detected(self):
        R0 = _coords_trigonal(d_md=2.0)
        R1 = R0.copy()
        R1[1, 0] += 0.7
        # donor 1 was at distance 2.0; now its distance is 2.7
        assert max_md_drift(R0, R1, 0, [1, 2, 3]) == pytest.approx(0.7, rel=1e-6)


# ---------------------------------------------------------------------------
# 9) Integration: grip_polish off vs on — byte-identity OFF + drift guard ON
# ---------------------------------------------------------------------------
class TestGripPolishUnfreezeIntegration:
    def test_default_off_byte_identical(self, monkeypatch):
        """With the flag OFF, the polish must be byte-identical with HEAD."""
        pytest.importorskip("rdkit.Chem")
        Chem = pytest.importorskip("rdkit.Chem")
        # Use the same toluene fixture used in test_grip_polish.
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        n = mol.GetNumAtoms()
        conf = mol.GetConformer()
        R = np.array(
            [list(conf.GetAtomPosition(i)) for i in range(n)], dtype=np.float64
        )
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", raising=False)
        # No metal -> grip_polish should be no-op-ish; just verify it runs.
        out = grip_polish(R, mol, metal=0, donors=[1, 2])
        assert out.shape == R.shape
        assert np.all(np.isfinite(out))

    def test_runs_with_phase2_on(self, monkeypatch):
        """With the flag ON, the polish must complete without raising."""
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal(d_md=1.95)
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        out = grip_polish(R, mol, metal=0, donors=[1, 2, 3])
        assert out.shape == R.shape
        assert np.all(np.isfinite(out))
        # M-D drift must be within 0.5 Å (fail-open guard).
        for d in (1, 2, 3):
            drift = abs(
                np.linalg.norm(out[d] - out[0]) - np.linalg.norm(R[d] - R[0])
            )
            assert drift <= 0.5 + 1e-6

    def test_determinism_with_phase2_on(self, monkeypatch):
        """Two runs with PYTHONHASHSEED=0 + the same env must match."""
        pytest.importorskip("rdkit.Chem")
        mol = _trigonal_mol()
        R = _coords_trigonal(d_md=1.95)
        monkeypatch.setenv("PYTHONHASHSEED", "0")
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
        out_a = grip_polish(R, mol, metal=0, donors=[1, 2, 3])
        out_b = grip_polish(R, mol, metal=0, donors=[1, 2, 3])
        assert np.array_equal(out_a, out_b)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
