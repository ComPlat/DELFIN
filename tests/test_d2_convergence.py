"""Mission D2 (2026-06-05) — convergence + completeness env-flag tests.

Tests the new env-vars added by Mission D2:

* DELFIN_FFFREE_GRIP_MAX_ITER  -> deeper L-BFGS convergence
* DELFIN_FFFREE_GRIP_GTOL      -> tighter gradient tolerance
* DELFIN_FFFREE_GRIP_RESTARTS  -> multi-restart L-BFGS
* DELFIN_FFFREE_GRIP_RESTART_PERTURB -> restart perturbation amplitude
* DELFIN_FFFREE_GRIP_WEIGHT_{BOND,ANGLE,IMPROPER,TORSION} -> per-class weights

Determinism gate: each restart uses a deterministic perturbation seeded by
``PYTHONHASHSEED`` so 2-run output is byte-identical.

Default-OFF byte-identical: with env unset, every resolver returns the
legacy default value (max_iter=200, gtol=1e-4, restarts=1, weights=(1.0,
0.5, 2.0, 0.2)) and code paths are equivalent to HEAD 8654d8f.
"""
from __future__ import annotations

import os

import numpy as np
import pytest

from delfin.fffree.grip_polish import (
    DEFAULT_GTOL,
    DEFAULT_MAX_ITER,
    DEFAULT_RESTART_PERTURB,
    DEFAULT_RESTARTS,
    _deterministic_perturbation,
    _resolve_gtol,
    _resolve_max_iter,
    _resolve_restart_perturb,
    _resolve_restarts,
)
from delfin.fffree.grip_fragment_detect import (
    WEIGHT_ANGLE_ENV,
    WEIGHT_BOND_ENV,
    WEIGHT_IMPROPER_ENV,
    WEIGHT_TORSION_ENV,
    _apply_weight_env_overrides,
    _resolve_weight_env,
)


# ---------------------------------------------------------------------------
# Fix 1 + 2: max_iter / gtol resolvers
# ---------------------------------------------------------------------------
def test_resolve_max_iter_default(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_MAX_ITER", raising=False)
    assert _resolve_max_iter(None) == DEFAULT_MAX_ITER == 200


def test_resolve_max_iter_env(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_MAX_ITER", "1000")
    assert _resolve_max_iter(None) == 1000


def test_resolve_max_iter_arg_wins(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_MAX_ITER", "1000")
    assert _resolve_max_iter(500) == 500


def test_resolve_max_iter_invalid_env_falls_back(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_MAX_ITER", "not_a_number")
    assert _resolve_max_iter(None) == DEFAULT_MAX_ITER


def test_resolve_max_iter_zero_falls_back(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_MAX_ITER", "0")
    assert _resolve_max_iter(None) == DEFAULT_MAX_ITER


def test_resolve_gtol_default(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_GTOL", raising=False)
    assert _resolve_gtol(None) == pytest.approx(DEFAULT_GTOL)


def test_resolve_gtol_env(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_GTOL", "1e-6")
    assert _resolve_gtol(None) == pytest.approx(1e-6)


def test_resolve_gtol_arg_wins(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_GTOL", "1e-6")
    assert _resolve_gtol(1e-5) == pytest.approx(1e-5)


def test_resolve_gtol_invalid_env_falls_back(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_GTOL", "neg-tol")
    assert _resolve_gtol(None) == pytest.approx(DEFAULT_GTOL)


def test_resolve_gtol_negative_falls_back(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_GTOL", "-1.0")
    assert _resolve_gtol(None) == pytest.approx(DEFAULT_GTOL)


# ---------------------------------------------------------------------------
# Fix 3: restarts + perturbation
# ---------------------------------------------------------------------------
def test_resolve_restarts_default(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_RESTARTS", raising=False)
    assert _resolve_restarts() == DEFAULT_RESTARTS == 1


def test_resolve_restarts_env(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_RESTARTS", "3")
    assert _resolve_restarts() == 3


def test_resolve_restarts_invalid_falls_back(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_RESTARTS", "garbage")
    assert _resolve_restarts() == DEFAULT_RESTARTS


def test_resolve_restart_perturb_default(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_RESTART_PERTURB", raising=False)
    assert _resolve_restart_perturb() == pytest.approx(DEFAULT_RESTART_PERTURB)


def test_resolve_restart_perturb_env(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_RESTART_PERTURB", "0.1")
    assert _resolve_restart_perturb() == pytest.approx(0.1)


def test_deterministic_perturbation_byte_identical(monkeypatch):
    """Same restart_idx + same P_init -> byte-identical output."""
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    rng = np.random.RandomState(42)
    P_init = rng.rand(30, 3).astype(np.float64)
    frozen = np.array([0, 5, 10], dtype=np.int64)
    P_a = _deterministic_perturbation(P_init, 1, 0.05, frozen)
    P_b = _deterministic_perturbation(P_init, 1, 0.05, frozen)
    assert np.array_equal(P_a, P_b)


def test_deterministic_perturbation_zero_restart_identity():
    """restart_idx=0 returns a copy of P_init unchanged."""
    P_init = np.arange(30, dtype=np.float64).reshape(10, 3)
    frozen = np.array([], dtype=np.int64)
    P_a = _deterministic_perturbation(P_init, 0, 0.05, frozen)
    assert np.array_equal(P_a, P_init)
    # is a copy, not the same array
    assert P_a is not P_init


def test_deterministic_perturbation_zero_amplitude():
    """amplitude=0 returns a copy of P_init unchanged."""
    P_init = np.arange(30, dtype=np.float64).reshape(10, 3)
    frozen = np.array([], dtype=np.int64)
    P_a = _deterministic_perturbation(P_init, 5, 0.0, frozen)
    assert np.array_equal(P_a, P_init)


def test_deterministic_perturbation_frozen_unchanged():
    """Frozen atoms (metal + donors) are not perturbed."""
    P_init = np.arange(30, dtype=np.float64).reshape(10, 3)
    frozen = np.array([0, 3, 7], dtype=np.int64)
    P_a = _deterministic_perturbation(P_init, 2, 0.5, frozen)
    for i in frozen:
        assert np.allclose(P_a[i], P_init[i])


def test_deterministic_perturbation_different_restart_idx_differ():
    """Different restart_idx -> different perturbations."""
    P_init = np.arange(30, dtype=np.float64).reshape(10, 3)
    frozen = np.array([], dtype=np.int64)
    P_a = _deterministic_perturbation(P_init, 1, 0.1, frozen)
    P_b = _deterministic_perturbation(P_init, 2, 0.1, frozen)
    assert not np.array_equal(P_a, P_b)


def test_deterministic_perturbation_amplitude_bound():
    """Each coordinate delta is within +/- amplitude."""
    P_init = np.zeros((50, 3), dtype=np.float64)
    frozen = np.array([], dtype=np.int64)
    amp = 0.05
    P_a = _deterministic_perturbation(P_init, 1, amp, frozen)
    assert np.all(np.abs(P_a) <= amp + 1e-12)


# ---------------------------------------------------------------------------
# Fix 6: per-class weight env overrides
# ---------------------------------------------------------------------------
def test_weight_env_default_byte_identical(monkeypatch):
    """No env-flags -> dict unchanged byte-identically."""
    for env in (WEIGHT_BOND_ENV, WEIGHT_ANGLE_ENV,
                WEIGHT_IMPROPER_ENV, WEIGHT_TORSION_ENV):
        monkeypatch.delenv(env, raising=False)
    w = {"bond": 1.0, "angle": 0.5, "improper": 2.0, "torsion": 0.2}
    w_orig = dict(w)
    w_out = _apply_weight_env_overrides(w)
    assert w_out == w_orig


def test_weight_env_bond_override(monkeypatch):
    monkeypatch.setenv(WEIGHT_BOND_ENV, "2.5")
    w = {"bond": 1.0, "angle": 0.5, "improper": 2.0, "torsion": 0.2}
    w_out = _apply_weight_env_overrides(w)
    assert w_out["bond"] == pytest.approx(2.5)
    assert w_out["angle"] == pytest.approx(0.5)  # untouched
    assert w_out["improper"] == pytest.approx(2.0)
    assert w_out["torsion"] == pytest.approx(0.2)


def test_weight_env_angle_override(monkeypatch):
    monkeypatch.setenv(WEIGHT_ANGLE_ENV, "1.5")
    w = {"bond": 1.0, "angle": 0.5, "improper": 2.0, "torsion": 0.2}
    w_out = _apply_weight_env_overrides(w)
    assert w_out["angle"] == pytest.approx(1.5)
    assert w_out["bond"] == pytest.approx(1.0)


def test_weight_env_all_overrides(monkeypatch):
    monkeypatch.setenv(WEIGHT_BOND_ENV, "2.5")
    monkeypatch.setenv(WEIGHT_ANGLE_ENV, "1.5")
    monkeypatch.setenv(WEIGHT_IMPROPER_ENV, "3.0")
    monkeypatch.setenv(WEIGHT_TORSION_ENV, "0.5")
    w = {"bond": 1.0, "angle": 0.5, "improper": 2.0, "torsion": 0.2}
    w_out = _apply_weight_env_overrides(w)
    assert w_out == {"bond": 2.5, "angle": 1.5, "improper": 3.0, "torsion": 0.5}


def test_weight_env_invalid_falls_back(monkeypatch):
    """Invalid env value -> the original / default value is preserved."""
    monkeypatch.setenv(WEIGHT_BOND_ENV, "not_a_float")
    w = {"bond": 1.0}
    assert _resolve_weight_env(WEIGHT_BOND_ENV, w["bond"]) == pytest.approx(1.0)


def test_weight_env_negative_falls_back(monkeypatch):
    """Negative weights are nonsense for a Mahalanobis loss -> fall back."""
    monkeypatch.setenv(WEIGHT_BOND_ENV, "-1.0")
    assert _resolve_weight_env(WEIGHT_BOND_ENV, 1.0) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Fix 4: DELFIN_MAX_ISOMERS env-var lifts max_isomers
# ---------------------------------------------------------------------------
def test_max_isomers_env_lift(monkeypatch):
    """Env-var lifts caller's max_isomers to the env value when env > arg."""
    pytest.importorskip("rdkit")
    monkeypatch.setenv("DELFIN_MAX_ISOMERS", "300")
    # Patch the impl to capture what max_isomers it actually sees
    from delfin import smiles_converter as sc
    seen = {}

    def fake_impl(*args, **kwargs):
        seen["max_isomers"] = kwargs.get("max_isomers")
        # Match real signature: positional max_isomers
        if seen["max_isomers"] is None and len(args) >= 3:
            seen["max_isomers"] = args[2]
        return ([], None)

    monkeypatch.setattr(sc, "_smiles_to_xyz_isomers_impl", fake_impl)
    sc.smiles_to_xyz_isomers("CCO", max_isomers=50)
    assert seen["max_isomers"] == 300, f"env-cap-lift failed: got {seen}"


def test_max_isomers_env_does_not_shrink(monkeypatch):
    """Env-cap NEVER shrinks a larger user-supplied value (only lifts up)."""
    pytest.importorskip("rdkit")
    monkeypatch.setenv("DELFIN_MAX_ISOMERS", "100")
    from delfin import smiles_converter as sc
    seen = {}

    def fake_impl(*args, **kwargs):
        seen["max_isomers"] = kwargs.get("max_isomers")
        if seen["max_isomers"] is None and len(args) >= 3:
            seen["max_isomers"] = args[2]
        return ([], None)

    monkeypatch.setattr(sc, "_smiles_to_xyz_isomers_impl", fake_impl)
    sc.smiles_to_xyz_isomers("CCO", max_isomers=500)
    # 500 > 100, so env cap does NOT shrink it
    assert seen["max_isomers"] == 500, f"env-cap shrunk caller: got {seen}"


def test_max_isomers_env_off_no_change(monkeypatch):
    """Env unset -> caller arg flows through byte-identically."""
    pytest.importorskip("rdkit")
    monkeypatch.delenv("DELFIN_MAX_ISOMERS", raising=False)
    from delfin import smiles_converter as sc
    seen = {}

    def fake_impl(*args, **kwargs):
        seen["max_isomers"] = kwargs.get("max_isomers")
        if seen["max_isomers"] is None and len(args) >= 3:
            seen["max_isomers"] = args[2]
        return ([], None)

    monkeypatch.setattr(sc, "_smiles_to_xyz_isomers_impl", fake_impl)
    sc.smiles_to_xyz_isomers("CCO", max_isomers=50)
    assert seen["max_isomers"] == 50


def test_max_isomers_env_invalid_ignored(monkeypatch):
    """Invalid env value -> caller arg flows through."""
    pytest.importorskip("rdkit")
    monkeypatch.setenv("DELFIN_MAX_ISOMERS", "not_a_number")
    from delfin import smiles_converter as sc
    seen = {}

    def fake_impl(*args, **kwargs):
        seen["max_isomers"] = kwargs.get("max_isomers")
        if seen["max_isomers"] is None and len(args) >= 3:
            seen["max_isomers"] = args[2]
        return ([], None)

    monkeypatch.setattr(sc, "_smiles_to_xyz_isomers_impl", fake_impl)
    sc.smiles_to_xyz_isomers("CCO", max_isomers=50)
    assert seen["max_isomers"] == 50


# ---------------------------------------------------------------------------
# Fix 5: GRACE max_per_isomer env override works (already in place upstream;
# we re-verify the bridge from the env-flag to the resolver).
# ---------------------------------------------------------------------------
def test_grace_max_per_isomer_env_override(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRACE_MAX_PER_ISOMER", "100")
    from delfin.fffree.grace_ensemble import grace_resolve_max_per_isomer
    assert grace_resolve_max_per_isomer(None) == 100


def test_grace_max_per_isomer_arg_wins(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRACE_MAX_PER_ISOMER", "100")
    from delfin.fffree.grace_ensemble import grace_resolve_max_per_isomer
    assert grace_resolve_max_per_isomer(50) == 50


def test_grace_max_isomers_env_override(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_GRACE_MAX_ISOMERS", "300")
    from delfin.fffree.grace_ensemble import grace_resolve_max_isomers
    assert grace_resolve_max_isomers(None) == 300


# ---------------------------------------------------------------------------
# Determinism round-trip: env-set + 2 calls of grip_polish are byte-identical
# ---------------------------------------------------------------------------
def test_grip_polish_with_d2_flags_determinism(monkeypatch):
    """Calling grip_polish twice with the same env + same inputs is byte-id."""
    pytest.importorskip("scipy")
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_MAX_ITER", "1000")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_GTOL", "1e-6")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_RESTARTS", "3")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_WEIGHT_BOND", "2.5")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_WEIGHT_ANGLE", "1.5")

    # Tiny "metal complex" -- octahedral Cu with 6 N donors via NH3
    smi = "[NH3][Cu]([NH3])([NH3])([NH3])([NH3])[NH3]"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        pytest.skip("rdkit could not parse the test SMILES")
    mol = Chem.AddHs(mol)
    res = AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    if res != 0:
        pytest.skip("ETKDG embedding failed -- environment too lean for full polish test")
    conf = mol.GetConformer()
    P = np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )

    # Find Cu + N donors
    metal_idx = next(
        i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "Cu"
    )
    donor_idxs = tuple(
        i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() == "N"
    )

    from delfin.fffree.grip_polish import grip_polish

    P_a = grip_polish(P, mol, metal_idx, donor_idxs, geom="OC-6")
    P_b = grip_polish(P, mol, metal_idx, donor_idxs, geom="OC-6")
    # Byte-identical
    assert np.array_equal(P_a, P_b), "grip_polish is not deterministic under D2 flags"
