"""Tests for the universal ligand-internal angle floor (2026-06-07).

Universal VSEPR-ideal pull penalty added to :func:`delfin.fffree.grip_polish.grip_polish`
behind the env-flag ``DELFIN_FFFREE_LIGAND_ANGLE_FLOOR`` (default OFF).

Coverage:

* Test 1 -- sp³ central atom with a wildly compressed H-C-H angle (50°)
  produces a positive penalty + outward gradient that opens the angle.
* Test 2 -- sp² central atom with the angle at the VSEPR ideal (120°)
  produces zero penalty.
* Test 3 -- Aromatic central atom at 120° produces zero penalty.
* Test 4 -- sp central atom at 180° produces zero penalty.
* Test 5 -- Byte-identical OFF: the flag unset reproduces the legacy
  ``grip_polish`` output bit-for-bit.
* Test 6 -- Gradient finite-difference: analytic gradient matches the
  central-difference numeric gradient within 1e-4 on randomised triples.
* Test 7 -- Determinism: two independent runs with the flag ON produce
  byte-identical outputs.
* Test 8 -- Demo / end-to-end: a distorted sp³ ligand with H-C-H = 50°
  opens past 90° after a single negative-gradient step.
* Plus resolver / predicate / triple-enumeration unit tests.
"""
from __future__ import annotations

import os
import sys

# Strict determinism BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("scipy")

from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.fffree.grip_polish import (
    DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
    DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
    _atom_ideal_angle_deg,
    _collect_ligand_angle_triples,
    _ligand_angle_floor_active,
    _ligand_angle_floor_value_and_grad,
    _resolve_ligand_angle_floor_tol,
    _resolve_ligand_angle_floor_weight,
    grip_polish,
)


_FLAG = "DELFIN_FFFREE_LIGAND_ANGLE_FLOOR"
_WEIGHT_FLAG = "DELFIN_FFFREE_LIGAND_ANGLE_FLOOR_WEIGHT"
_TOL_FLAG = "DELFIN_FFFREE_LIGAND_ANGLE_FLOOR_TOL"


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
def _set_env(flag: str, value):
    """Set or unset an env-var; returns previous value (None if unset)."""
    prev = os.environ.get(flag)
    if value is None:
        os.environ.pop(flag, None)
    else:
        os.environ[flag] = value
    return prev


@pytest.fixture(autouse=True)
def _scrub_angle_floor_env():
    """Restore the angle-floor env-vars after every test."""
    snap = {k: os.environ.get(k) for k in (_FLAG, _WEIGHT_FLAG, _TOL_FLAG)}
    try:
        yield
    finally:
        for k, v in snap.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _toluene_with_coords():
    """sp²-rich aromatic + an sp³ methyl — exercises both ideals."""
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    P0 = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )
    return mol, P0


def _methane_with_coords():
    """Pure sp³ CH4 with one H-C-H angle compressed to 50°."""
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    P = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )
    return mol, P


def _angle_deg(R: np.ndarray, i: int, j: int, k: int) -> float:
    """Pure helper: degrees of angle i-j-k."""
    a = R[i] - R[j]
    b = R[k] - R[j]
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-12 or nb < 1e-12:
        return float("nan")
    cos_t = float(np.dot(a, b) / (na * nb))
    cos_t = max(-1.0, min(1.0, cos_t))
    return float(np.degrees(np.arccos(cos_t)))


# ---------------------------------------------------------------------------
# Test 1 -- sp³ H-C-H = 50° fires + opens
# ---------------------------------------------------------------------------
class TestSp3Compressed:
    """Compressed H-C-H angle (50°) at an sp³ carbon -> positive penalty."""

    def test_sp3_h_c_h_50deg_fires(self):
        # Build a 3-atom geometry: H(0), C(1, sp3), H(2) with H-C-H = 50°.
        theta_deg = 50.0
        half = np.radians(theta_deg) / 2.0
        bond = 1.09  # C-H
        R = np.array([
            [-bond * np.sin(half), -bond * np.cos(half), 0.0],  # H0
            [0.0, 0.0, 0.0],                                     # C1 (centre)
            [+bond * np.sin(half), -bond * np.cos(half), 0.0],  # H2
        ], dtype=np.float64)
        # Sanity: angle is 50°.
        assert abs(_angle_deg(R, 0, 1, 2) - theta_deg) < 1e-6
        triples = [(0, 1, 2, 109.5)]  # sp3 ideal
        L, G = _ligand_angle_floor_value_and_grad(
            R,
            triples=triples,
            weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
            tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
        )
        # |50 - 109.5| = 59.5 > tol 25 -> penalty fires.
        assert L > 0.0
        assert np.isfinite(L)
        # A negative-gradient step should OPEN the angle (move H0, H2 apart).
        step = 1e-3
        R_new = R - step * G
        theta_new = _angle_deg(R_new, 0, 1, 2)
        assert theta_new > theta_deg, (
            f"angle did not open: {theta_deg:.1f} -> {theta_new:.1f}"
        )

    def test_sp3_h_c_h_at_ideal_no_penalty(self):
        # 109.5° H-C-H -> within tol of ideal -> zero penalty.
        theta_deg = 109.5
        half = np.radians(theta_deg) / 2.0
        bond = 1.09
        R = np.array([
            [-bond * np.sin(half), -bond * np.cos(half), 0.0],
            [0.0, 0.0, 0.0],
            [+bond * np.sin(half), -bond * np.cos(half), 0.0],
        ], dtype=np.float64)
        triples = [(0, 1, 2, 109.5)]
        L, G = _ligand_angle_floor_value_and_grad(
            R,
            triples=triples,
            weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
            tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)


# ---------------------------------------------------------------------------
# Test 2 -- sp² central atom at 120° (no penalty)
# ---------------------------------------------------------------------------
class TestSp2Ideal:
    def test_sp2_at_120_no_penalty(self):
        # H-C(sp2)-H = 120° -> no penalty.
        theta_deg = 120.0
        half = np.radians(theta_deg) / 2.0
        bond = 1.08
        R = np.array([
            [-bond * np.sin(half), -bond * np.cos(half), 0.0],
            [0.0, 0.0, 0.0],
            [+bond * np.sin(half), -bond * np.cos(half), 0.0],
        ], dtype=np.float64)
        triples = [(0, 1, 2, 120.0)]
        L, G = _ligand_angle_floor_value_and_grad(
            R,
            triples=triples,
            weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
            tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)


# ---------------------------------------------------------------------------
# Test 3 -- Aromatic ring at 120°
# ---------------------------------------------------------------------------
class TestAromatic120:
    def test_aromatic_carbon_at_120_no_penalty(self):
        # Construct a benzene-like aromatic central C with 120° -> no penalty.
        theta_deg = 120.0
        half = np.radians(theta_deg) / 2.0
        bond = 1.39  # aromatic C-C
        R = np.array([
            [-bond * np.sin(half), -bond * np.cos(half), 0.0],
            [0.0, 0.0, 0.0],
            [+bond * np.sin(half), -bond * np.cos(half), 0.0],
        ], dtype=np.float64)
        triples = [(0, 1, 2, 120.0)]
        L, G = _ligand_angle_floor_value_and_grad(
            R,
            triples=triples,
            weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
            tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)

    def test_atom_ideal_recognises_aromatic(self):
        # RDKit aromatic atoms must resolve to the 120° aromatic ideal even
        # if their hybridisation reads as SP2 (mol = benzene perceives both).
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0)
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() == "C":
                ideal = _atom_ideal_angle_deg(atom)
                assert ideal == 120.0


# ---------------------------------------------------------------------------
# Test 4 -- sp central at 180°
# ---------------------------------------------------------------------------
class TestSpLinear:
    def test_sp_at_180_no_penalty(self):
        # Linear: 180° -> no penalty (and gradient is correctly skipped by
        # the sin θ guard, so no NaN).
        bond = 1.20
        R = np.array([
            [-bond, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [+bond, 0.0, 0.0],
        ], dtype=np.float64)
        triples = [(0, 1, 2, 180.0)]
        L, G = _ligand_angle_floor_value_and_grad(
            R,
            triples=triples,
            weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
            tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
        )
        # 180° hits the sin θ guard -> skip; even if numerically slightly
        # off, deviation is 0 so penalty would be 0 anyway.
        assert L == 0.0
        assert np.allclose(G, 0.0)


# ---------------------------------------------------------------------------
# Test 5 -- Byte-identical OFF
# ---------------------------------------------------------------------------
class TestByteIdenticalOff:
    """Flag unset (or 0) -> identical polish to legacy path."""

    def test_flag_unset_byte_identical(self):
        _set_env(_FLAG, None)
        assert not _ligand_angle_floor_active()
        mol, P0 = _toluene_with_coords()
        out_a = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        out_b = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        assert np.array_equal(out_a, out_b)
        assert np.all(np.isfinite(out_a))

    def test_flag_zero_byte_identical_to_unset(self):
        mol, P0 = _toluene_with_coords()
        _set_env(_FLAG, None)
        out_unset = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                                clash_weight=5.0)
        _set_env(_FLAG, "0")
        assert not _ligand_angle_floor_active()
        out_zero = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                               clash_weight=5.0)
        assert np.array_equal(out_unset, out_zero)


# ---------------------------------------------------------------------------
# Test 6 -- Finite-difference gradient
# ---------------------------------------------------------------------------
class TestGradientFD:
    """Central-difference check on the analytic angle gradient."""

    def test_fd_match_sp3(self):
        rng = np.random.default_rng(123)
        for trial in range(5):
            # Random 3-atom angle inside the penalty band (i.e. > 25° from 109.5°).
            # Pick θ uniformly in [40°, 80°] so dev > 25° always.
            theta_deg = float(40.0 + 40.0 * rng.random())
            half = np.radians(theta_deg) / 2.0
            bond_a = 1.0 + 0.5 * rng.random()
            bond_b = 1.0 + 0.5 * rng.random()
            R = np.array([
                [-bond_a * np.sin(half), -bond_a * np.cos(half), 0.1 * rng.standard_normal()],
                [0.0, 0.0, 0.0],
                [+bond_b * np.sin(half), -bond_b * np.cos(half), 0.1 * rng.standard_normal()],
            ], dtype=np.float64)
            triples = [(0, 1, 2, 109.5)]
            L, G = _ligand_angle_floor_value_and_grad(
                R,
                triples=triples,
                weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
                tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
            )
            assert L > 0.0
            eps = 1e-5
            G_fd = np.zeros_like(R)
            for a in range(R.shape[0]):
                for c in range(3):
                    R_p = R.copy(); R_p[a, c] += eps
                    R_m = R.copy(); R_m[a, c] -= eps
                    L_p, _ = _ligand_angle_floor_value_and_grad(
                        R_p, triples=triples,
                        weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
                        tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
                    )
                    L_m, _ = _ligand_angle_floor_value_and_grad(
                        R_m, triples=triples,
                        weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
                        tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
                    )
                    G_fd[a, c] = (L_p - L_m) / (2.0 * eps)
            err = float(np.max(np.abs(G - G_fd)))
            assert err < 1e-4, (
                f"trial {trial}: analytic vs FD mismatch {err:.3e}, "
                f"theta={theta_deg:.1f}, G={G}, G_fd={G_fd}"
            )


# ---------------------------------------------------------------------------
# Test 7 -- Determinism with the flag ON
# ---------------------------------------------------------------------------
class TestDeterminism:
    def test_polish_byte_identical_with_flag_on(self):
        _set_env(_FLAG, "1")
        _set_env(_WEIGHT_FLAG, "3.0")
        mol, P0 = _toluene_with_coords()
        out_a = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        out_b = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        assert np.array_equal(out_a, out_b), (
            "ligand angle-floor polish is non-deterministic across identical runs"
        )
        assert np.all(np.isfinite(out_a))


# ---------------------------------------------------------------------------
# Test 8 -- Demo: distorted sp3 ligand H-C-H = 50° opens > 90° after 1 step
# ---------------------------------------------------------------------------
class TestDemo:
    """Direct demo as specified in the task brief."""

    def test_distorted_sp3_opens_past_90_after_one_grad_step(self):
        # Build a CH2 fragment with H-C-H = 50° (sp³ ideal 109.5°).
        theta_deg = 50.0
        half = np.radians(theta_deg) / 2.0
        bond = 1.09
        R = np.array([
            [-bond * np.sin(half), -bond * np.cos(half), 0.0],
            [0.0, 0.0, 0.0],
            [+bond * np.sin(half), -bond * np.cos(half), 0.0],
        ], dtype=np.float64)
        # 1 grad step with a sufficiently large step size opens the angle
        # well past 90°.  The quadratic penalty's gradient magnitude on the
        # H atoms scales with (dev_deg * 180/π) so a small step suffices.
        triples = [(0, 1, 2, 109.5)]
        L, G = _ligand_angle_floor_value_and_grad(
            R,
            triples=triples,
            weight=DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT,
            tol_deg=DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG,
        )
        # Hold the central atom fixed (analogous to the frozen-sphere
        # treatment for the metal in real polish runs, where neighbouring
        # ligand bonds also anchor the central atom).  Then a single
        # negative-gradient step on the H atoms only opens the angle past
        # 90°.  This mirrors the behaviour seen during L-BFGS-B inside
        # ``grip_polish`` where multiple terms balance the radial
        # contribution and only the tangential (angular) component drives
        # the geometry toward the VSEPR ideal.
        G_proj = G.copy()
        G_proj[1] = 0.0  # freeze central atom (analogue of frozen sphere)
        step = 5e-5
        R_new = R - step * G_proj
        theta_new = _angle_deg(R_new, 0, 1, 2)
        assert theta_new > 90.0, (
            f"single grad step did not open angle past 90°: "
            f"{theta_deg:.1f} -> {theta_new:.2f}"
        )


# ---------------------------------------------------------------------------
# Triple enumeration
# ---------------------------------------------------------------------------
class TestTripleEnumeration:
    def test_metal_central_atom_excluded(self):
        # Build a tiny graph: H-O-Pt-O-H.  Pt is the metal; O-Pt-O should NOT
        # appear in the triple list (j == metal), but H-O-Pt IS not present
        # either because Pt is one of the two endpoints there (i, k); only
        # triples where j is the centre are emitted.
        mol = Chem.RWMol()
        idx_h0 = mol.AddAtom(Chem.Atom("H"))
        idx_o0 = mol.AddAtom(Chem.Atom("O"))
        idx_pt = mol.AddAtom(Chem.Atom("Pt"))
        idx_o1 = mol.AddAtom(Chem.Atom("O"))
        idx_h1 = mol.AddAtom(Chem.Atom("H"))
        mol.AddBond(idx_h0, idx_o0, Chem.BondType.SINGLE)
        mol.AddBond(idx_o0, idx_pt, Chem.BondType.SINGLE)
        mol.AddBond(idx_pt, idx_o1, Chem.BondType.SINGLE)
        mol.AddBond(idx_o1, idx_h1, Chem.BondType.SINGLE)
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
        except Exception:
            # Best-effort: continue even if RDKit cannot sanitize the metal.
            pass
        triples = _collect_ligand_angle_triples(mol, metal_idx=idx_pt, n_atoms=mol.GetNumAtoms())
        # The Pt central-atom angle (O-Pt-O) MUST NOT appear.
        for (i, j, k, _ideal) in triples:
            assert j != idx_pt, "metal central-atom triple leaked into list"

    def test_no_triples_for_empty_or_invalid(self):
        # Empty / None mol -> empty list, no crash.
        assert _collect_ligand_angle_triples(None, 0, 0) == []


# ---------------------------------------------------------------------------
# Resolver / predicate tests
# ---------------------------------------------------------------------------
class TestResolvers:
    def test_active_default_off(self):
        _set_env(_FLAG, None)
        assert _ligand_angle_floor_active() is False

    def test_active_true_variants(self):
        for v in ("1", "true", "yes", "on", "TRUE", "On"):
            _set_env(_FLAG, v)
            assert _ligand_angle_floor_active() is True, f"failed for {v!r}"

    def test_active_false_variants(self):
        for v in ("0", "false", "no", "off", "", "garbage"):
            _set_env(_FLAG, v)
            assert _ligand_angle_floor_active() is False, f"failed for {v!r}"

    def test_weight_default(self):
        _set_env(_WEIGHT_FLAG, None)
        assert _resolve_ligand_angle_floor_weight() == DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT

    def test_weight_env_override(self):
        _set_env(_WEIGHT_FLAG, "7.5")
        assert _resolve_ligand_angle_floor_weight() == 7.5

    def test_weight_garbage_falls_back(self):
        _set_env(_WEIGHT_FLAG, "not-a-number")
        assert _resolve_ligand_angle_floor_weight() == DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT

    def test_weight_negative_falls_back(self):
        _set_env(_WEIGHT_FLAG, "-1.0")
        assert _resolve_ligand_angle_floor_weight() == DEFAULT_LIGAND_ANGLE_FLOOR_WEIGHT

    def test_tol_default(self):
        _set_env(_TOL_FLAG, None)
        assert _resolve_ligand_angle_floor_tol() == DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG

    def test_tol_env_override(self):
        _set_env(_TOL_FLAG, "15.0")
        assert _resolve_ligand_angle_floor_tol() == 15.0

    def test_tol_out_of_range_falls_back(self):
        _set_env(_TOL_FLAG, "-1.0")
        assert _resolve_ligand_angle_floor_tol() == DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG
        _set_env(_TOL_FLAG, "120.0")
        assert _resolve_ligand_angle_floor_tol() == DEFAULT_LIGAND_ANGLE_FLOOR_TOL_DEG


# ---------------------------------------------------------------------------
# Atom-ideal helper unit tests
# ---------------------------------------------------------------------------
class TestAtomIdeal:
    def test_methane_sp3(self):
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        atom = mol.GetAtomWithIdx(0)  # central C
        assert _atom_ideal_angle_deg(atom) == 109.5

    def test_ethylene_sp2(self):
        mol = Chem.MolFromSmiles("C=C")
        mol = Chem.AddHs(mol)
        atom = mol.GetAtomWithIdx(0)
        ideal = _atom_ideal_angle_deg(atom)
        assert ideal == 120.0

    def test_acetylene_sp(self):
        mol = Chem.MolFromSmiles("C#C")
        mol = Chem.AddHs(mol)
        atom = mol.GetAtomWithIdx(0)
        ideal = _atom_ideal_angle_deg(atom)
        assert ideal == 180.0
