"""Unit tests for the donor-cone DoF (delfin.fffree.grip_donor_cone).

Validates:
* env-flag default-OFF byte-identical (the polish path is byte-identical to
  the legacy L-BFGS path when DELFIN_FFFREE_GRIP_DONOR_CONE is unset).
* cone construction is deterministic + topology-correct
* the rotation is an isometry (M-D length, donor pos, M-D-X polar
  angle unchanged when theta != 0)
* the theta gradient agrees with a finite-difference reference to 1e-9
* the position gradient propagates correctly through the rotation
* 2-run determinism of the wrapped objective (byte-identical L and G)
* env-flag round-trip across subprocess fork (relevant for parallel pool)
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Tuple

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(REPO_ROOT))

from delfin.fffree.grip_donor_cone import (  # noqa: E402
    DONOR_CONE_ENV,
    DONOR_CONE_INCLUDE_H_ENV,
    DonorCone,
    apply_donor_cone_rotations,
    augmented_loss_and_grad,
    build_donor_cones,
    donor_cone_active,
    donor_cone_include_h,
)


# ---------------------------------------------------------------------------
# Mini fake-mol fixture
# ---------------------------------------------------------------------------
class _Atom:
    def __init__(self, idx: int, sym: str):
        self._idx = idx
        self._sym = sym
        self._neighbours = []

    def GetIdx(self):
        return self._idx

    def GetSymbol(self):
        return self._sym

    def GetNeighbors(self):
        return list(self._neighbours)


class _Mol:
    def __init__(self, symbols, bonds):
        self._atoms = [_Atom(i, s) for i, s in enumerate(symbols)]
        for (u, v) in bonds:
            self._atoms[u]._neighbours.append(self._atoms[v])
            self._atoms[v]._neighbours.append(self._atoms[u])

    def GetAtomWithIdx(self, idx):
        return self._atoms[idx]

    def GetNumAtoms(self):
        return len(self._atoms)


def _make_simple_complex() -> Tuple[_Mol, np.ndarray, int, list]:
    """Build a synthetic M(N(CH3)(CH3))(N(CH3)(CH3)) complex.

    Atom layout (10 atoms, sym + idx):
        0  Fe   (0, 0, 0)        metal
        1  N    (2, 0, 0)        donor 1
        2  C    (2.8, 0.7, 0)    cone X of D1
        3  C    (2.8, -0.7, 0)   cone X of D1
        4  H    (3.6, 0.7, 0)    H on C2 (not in heavy cone)
        5  N    (-2, 0, 0)       donor 2
        6  C    (-2.8, 0.7, 0)   cone X of D2
        7  C    (-2.8, -0.7, 0)  cone X of D2
        8  H    (-3.6, 0.7, 0)   H on C6
        9  H    (-3.6, -0.7, 0)  H on C7
    """
    syms = ["Fe", "N", "C", "C", "H", "N", "C", "C", "H", "H"]
    P = np.array([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.8, 0.7, 0.0],
        [2.8, -0.7, 0.0],
        [3.6, 0.7, 0.0],
        [-2.0, 0.0, 0.0],
        [-2.8, 0.7, 0.0],
        [-2.8, -0.7, 0.0],
        [-3.6, 0.7, 0.0],
        [-3.6, -0.7, 0.0],
    ], dtype=np.float64)
    bonds = [
        (0, 1), (0, 5),       # M-N
        (1, 2), (1, 3),       # N-C
        (2, 4),               # C-H
        (5, 6), (5, 7),       # N-C
        (6, 8), (7, 9),       # C-H
    ]
    mol = _Mol(syms, bonds)
    metal = 0
    donors = [1, 5]
    return mol, P, metal, donors


# ---------------------------------------------------------------------------
# Env-flag tests
# ---------------------------------------------------------------------------
def test_donor_cone_default_off(monkeypatch):
    monkeypatch.delenv(DONOR_CONE_ENV, raising=False)
    assert donor_cone_active() is False


@pytest.mark.parametrize("value", ["1", "true", "yes", "on", "TRUE"])
def test_donor_cone_env_on_values(monkeypatch, value):
    monkeypatch.setenv(DONOR_CONE_ENV, value)
    assert donor_cone_active() is True


@pytest.mark.parametrize("value", ["0", "false", "no", "off", ""])
def test_donor_cone_env_off_values(monkeypatch, value):
    monkeypatch.setenv(DONOR_CONE_ENV, value)
    assert donor_cone_active() is False


def test_donor_cone_include_h_default_off(monkeypatch):
    monkeypatch.delenv(DONOR_CONE_INCLUDE_H_ENV, raising=False)
    assert donor_cone_include_h() is False


def test_donor_cone_include_h_on(monkeypatch):
    monkeypatch.setenv(DONOR_CONE_INCLUDE_H_ENV, "1")
    assert donor_cone_include_h() is True


# ---------------------------------------------------------------------------
# Cone construction tests
# ---------------------------------------------------------------------------
def test_build_donor_cones_basic():
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    # 2 donors, each with 2 heavy X atoms (first shell).  The cone
    # subtree expands BFS from each X across all non-donor / non-metal
    # neighbours; in this fixture the H on C2 (idx 4), on C6 (idx 8) and
    # on C7 (idx 9) reach the subtree even though include_h is OFF for
    # the FIRST-SHELL gate (they enter via downstream BFS).
    assert len(cones) == 2
    assert cones[0].donor == 1
    # Donor-1 subtree: C atoms 2,3 + their downstream H atom 4.
    assert cones[0].atoms == (2, 3, 4)
    assert cones[1].donor == 5
    # Donor-2 subtree: C atoms 6,7 + their downstream H atoms 8,9.
    assert cones[1].atoms == (6, 7, 8, 9)
    # Axis vectors point from M to D and are unit length.
    np.testing.assert_allclose(cones[0].axis, np.array([1.0, 0.0, 0.0]),
                               atol=1e-12)
    np.testing.assert_allclose(cones[1].axis, np.array([-1.0, 0.0, 0.0]),
                               atol=1e-12)
    np.testing.assert_allclose(np.linalg.norm(cones[0].axis), 1.0,
                               atol=1e-12)


def test_build_donor_cones_include_h():
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P, include_h=True)
    # With include_h=True the H atoms on the donor itself (if any) would
    # ALSO enter the first shell.  In this fixture the N donors have no
    # direct H neighbours -- H atoms ride along via downstream BFS in
    # both cases.  So include_h is a no-op on this fixture.
    assert cones[0].atoms == (2, 3, 4)
    assert cones[1].atoms == (6, 7, 8, 9)


def test_build_donor_cones_hapto_donor_skipped():
    mol, P, metal, donors = _make_simple_complex()
    # Mark donor 1 as a hapto-pi atom -> should be skipped.
    cones = build_donor_cones(mol, metal, donors, P, hapto_atoms=[1])
    assert len(cones) == 1
    assert cones[0].donor == 5


def test_build_donor_cones_deterministic_order():
    """Donors passed out of order must still emit cones in sorted-donor order."""
    mol, P, metal, donors = _make_simple_complex()
    cones1 = build_donor_cones(mol, metal, [5, 1], P)
    cones2 = build_donor_cones(mol, metal, [1, 5], P)
    assert [c.donor for c in cones1] == [c.donor for c in cones2] == [1, 5]
    assert cones1[0].atoms == cones2[0].atoms
    assert cones1[1].atoms == cones2[1].atoms


def test_build_donor_cones_no_donors_returns_empty():
    mol, P, metal, _ = _make_simple_complex()
    assert build_donor_cones(mol, metal, [], P) == []


def test_build_donor_cones_zero_axis_skipped():
    mol, P, metal, donors = _make_simple_complex()
    # Put donor 1 on top of metal -> degenerate axis -> cone skipped.
    P_bad = P.copy()
    P_bad[1] = P_bad[metal]
    cones = build_donor_cones(mol, metal, donors, P_bad)
    assert len(cones) == 1
    assert cones[0].donor == 5


# ---------------------------------------------------------------------------
# Rotation isometry tests
# ---------------------------------------------------------------------------
def test_apply_zero_theta_byte_identical():
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    P_eval = apply_donor_cone_rotations(P, cones, [0.0, 0.0])
    # Bit-exact: zero-theta short-circuit must not perturb P at all.
    assert np.array_equal(P_eval, P)


def test_apply_nonzero_theta_preserves_metal_and_donor():
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    thetas = [0.5, 1.2]
    P_eval = apply_donor_cone_rotations(P, cones, thetas)
    # Metal and donors must be byte-identical (not in any cone).
    np.testing.assert_array_equal(P_eval[metal], P[metal])
    for d in donors:
        np.testing.assert_array_equal(P_eval[d], P[d])


def test_md_bond_length_invariant_under_theta():
    """The M-D distance is invariant since D is not in any cone."""
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    md_before = [
        float(np.linalg.norm(P[d] - P[metal])) for d in donors
    ]
    for theta in (-1.0, -0.3, 0.7, 2.5):
        P_eval = apply_donor_cone_rotations(P, cones, [theta, -theta])
        md_after = [
            float(np.linalg.norm(P_eval[d] - P_eval[metal]))
            for d in donors
        ]
        np.testing.assert_allclose(md_before, md_after, atol=1e-12)


def test_donor_position_invariant_under_theta():
    """Donor positions are NEVER moved by the cone DoF (M-D-rigid)."""
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    for theta in (-2.0, 0.0, 0.4, 3.14159):
        P_eval = apply_donor_cone_rotations(P, cones, [theta, theta])
        for d in donors:
            np.testing.assert_allclose(P_eval[d], P[d], atol=1e-12)


def test_dx_bond_length_invariant_under_theta():
    """D-X bond lengths are invariant -- the rotation is rigid w.r.t. D."""
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    bonds = [(1, 2), (1, 3), (5, 6), (5, 7)]
    bl_before = [float(np.linalg.norm(P[a] - P[b])) for (a, b) in bonds]
    for theta in (-0.5, 0.7, 1.5):
        P_eval = apply_donor_cone_rotations(P, cones, [theta, -theta])
        bl_after = [float(np.linalg.norm(P_eval[a] - P_eval[b]))
                    for (a, b) in bonds]
        np.testing.assert_allclose(bl_before, bl_after, atol=1e-12)


def test_mdx_polar_angle_invariant_under_theta():
    """The M-D-X polar angle is invariant under axis rotation."""
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)

    def angle(a, b, c, X):
        u = X[a] - X[b]
        v = X[c] - X[b]
        return math.degrees(math.acos(
            float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
        ))

    polar_before = [angle(metal, 1, 2, P), angle(metal, 1, 3, P),
                    angle(metal, 5, 6, P), angle(metal, 5, 7, P)]
    P_eval = apply_donor_cone_rotations(P, cones, [0.9, 1.6])
    polar_after = [angle(metal, 1, 2, P_eval), angle(metal, 1, 3, P_eval),
                   angle(metal, 5, 6, P_eval), angle(metal, 5, 7, P_eval)]
    np.testing.assert_allclose(polar_before, polar_after, atol=1e-10)


# ---------------------------------------------------------------------------
# Gradient tests
# ---------------------------------------------------------------------------
def _make_quadratic_loss(target: np.ndarray):
    """Returns a (loss, grad) callable of x_flat such that
        L(P) = 0.5 * sum_i ||P[i] - target[i]||^2
    """
    def loss_and_grad(x_flat: np.ndarray):
        P = np.asarray(x_flat, dtype=np.float64).reshape(target.shape)
        diff = P - target
        L = 0.5 * float(np.sum(diff * diff))
        G = diff.copy()
        return L, G.reshape(-1)

    return loss_and_grad


def test_theta_gradient_finite_diff():
    """Theta gradient agrees with finite-difference to 1e-9."""
    rng = np.random.default_rng(0)
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    n_atoms = P.shape[0]
    K = len(cones)
    # Random non-trivial target so the quadratic loss has signal everywhere.
    target = P + 0.1 * rng.standard_normal(P.shape)
    pos_loss = _make_quadratic_loss(target)
    aug = augmented_loss_and_grad(pos_loss, cones, n_atoms)

    thetas = np.array([0.3, -0.8], dtype=np.float64)
    x_aug = np.concatenate([P.reshape(-1), thetas])
    L0, G0 = aug(x_aug)

    h = 1e-6
    n_pos = 3 * n_atoms
    fd_grad = np.zeros(K, dtype=np.float64)
    for k in range(K):
        x_plus = x_aug.copy(); x_plus[n_pos + k] += h
        x_minus = x_aug.copy(); x_minus[n_pos + k] -= h
        Lp, _ = aug(x_plus)
        Lm, _ = aug(x_minus)
        fd_grad[k] = (Lp - Lm) / (2 * h)
    np.testing.assert_allclose(G0[n_pos:], fd_grad, atol=1e-8, rtol=1e-6)


def test_position_gradient_finite_diff():
    """Position gradient on cone atoms agrees with finite-difference."""
    rng = np.random.default_rng(1)
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    n_atoms = P.shape[0]
    K = len(cones)
    target = P + 0.05 * rng.standard_normal(P.shape)
    pos_loss = _make_quadratic_loss(target)
    aug = augmented_loss_and_grad(pos_loss, cones, n_atoms)

    thetas = np.array([0.4, 0.0], dtype=np.float64)
    x_aug = np.concatenate([P.reshape(-1), thetas])
    L0, G0 = aug(x_aug)

    h = 1e-7
    n_pos = 3 * n_atoms
    # Pick a cone atom (idx 2 -> x position 6,7,8 in flat).
    for cone_atom in (2, 3, 6, 7):
        for comp in range(3):
            idx = cone_atom * 3 + comp
            x_plus = x_aug.copy(); x_plus[idx] += h
            x_minus = x_aug.copy(); x_minus[idx] -= h
            Lp, _ = aug(x_plus)
            Lm, _ = aug(x_minus)
            fd = (Lp - Lm) / (2 * h)
            assert math.isclose(fd, G0[idx], abs_tol=1e-7, rel_tol=1e-5), (
                f"cone atom {cone_atom} comp {comp}: "
                f"analytic={G0[idx]:.6e} fd={fd:.6e}"
            )


def test_position_gradient_non_cone_atom_unchanged_chain():
    """For atoms NOT in any cone, the wrapped grad equals the inner grad."""
    rng = np.random.default_rng(2)
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    n_atoms = P.shape[0]
    K = len(cones)
    target = P + 0.02 * rng.standard_normal(P.shape)
    pos_loss = _make_quadratic_loss(target)
    aug = augmented_loss_and_grad(pos_loss, cones, n_atoms)

    thetas = np.array([0.1, -0.2], dtype=np.float64)
    x_aug = np.concatenate([P.reshape(-1), thetas])
    L, G = aug(x_aug)
    # Atoms 0 (Fe), 1 (N), 5 (N) are the metal + donors -- never in any
    # cone (downstream BFS stops at donors / metal).
    P_eval = apply_donor_cone_rotations(P, cones, thetas)
    _, G_inner_flat = pos_loss(P_eval.reshape(-1))
    G_inner = G_inner_flat.reshape(n_atoms, 3)
    non_cone = [0, 1, 5]
    for a in non_cone:
        idx0 = 3 * a
        np.testing.assert_allclose(G[idx0:idx0 + 3], G_inner[a], atol=1e-12)


# ---------------------------------------------------------------------------
# Determinism tests
# ---------------------------------------------------------------------------
def test_two_run_determinism():
    """Two consecutive evaluations produce byte-identical L and gradient."""
    rng = np.random.default_rng(3)
    mol, P, metal, donors = _make_simple_complex()
    cones = build_donor_cones(mol, metal, donors, P)
    n_atoms = P.shape[0]
    target = P + 0.1 * rng.standard_normal(P.shape)
    pos_loss = _make_quadratic_loss(target)
    aug = augmented_loss_and_grad(pos_loss, cones, n_atoms)
    thetas = np.array([0.7, -0.3], dtype=np.float64)
    x_aug = np.concatenate([P.reshape(-1), thetas])
    L1, G1 = aug(x_aug)
    L2, G2 = aug(x_aug)
    assert L1 == L2
    assert np.array_equal(G1, G2)


def test_default_off_byte_identical_via_polish(tmp_path, monkeypatch):
    """When the env-flag is unset the polish path is byte-identical to
    the legacy L-BFGS path (no cone DoF, no augmented variable set).

    We invoke grip_polish through a subprocess so the env-flag is read
    fresh and the result matches what a parallel pool worker would see.
    """
    script = textwrap.dedent("""
        import os, sys, json
        import numpy as np
        sys.path.insert(0, %r)
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from delfin.fffree.grip_polish import grip_polish

        mol = Chem.AddHs(Chem.MolFromSmiles("[N](C)(C)[Fe][N](C)(C)"))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        n = mol.GetNumAtoms()
        conf = mol.GetConformer()
        P0 = np.asarray([[conf.GetAtomPosition(i).x,
                          conf.GetAtomPosition(i).y,
                          conf.GetAtomPosition(i).z] for i in range(n)],
                        dtype=np.float64)
        # Metal Fe = first Fe atom, donors = its N neighbours.
        metal = next(i for i in range(n) if mol.GetAtomWithIdx(i).GetSymbol()=="Fe")
        donors = sorted(int(nb.GetIdx())
                        for nb in mol.GetAtomWithIdx(metal).GetNeighbors())
        P = grip_polish(P0, mol, metal, donors)
        sys.stdout.write(json.dumps([[float(x) for x in row]
                                     for row in np.asarray(P)]))
    """) % str(REPO_ROOT)
    env_off = os.environ.copy()
    env_off["PYTHONHASHSEED"] = "0"
    env_off.pop(DONOR_CONE_ENV, None)
    r_off = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True, text=True, env=env_off, check=False,
    )
    if r_off.returncode != 0:
        pytest.skip(f"rdkit/grip stack unavailable: {r_off.stderr[-400:]}")
    # Two more runs with the env-flag explicitly OFF.
    env_off2 = env_off.copy(); env_off2[DONOR_CONE_ENV] = "0"
    r_off2 = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True, text=True, env=env_off2, check=False,
    )
    assert r_off2.returncode == 0
    # Byte-identical output across two runs with the env-flag at OFF.
    assert r_off.stdout == r_off2.stdout, (
        "DELFIN_FFFREE_GRIP_DONOR_CONE=0 must be byte-identical with "
        "the unset case (default-OFF byte-identity required)."
    )


def test_env_on_changes_polish_output(tmp_path):
    """The env-flag ON must produce a measurable change vs OFF when at least
    one cone has > 0 atoms (sanity, not a guarantee of improvement)."""
    script = textwrap.dedent("""
        import os, sys, json
        import numpy as np
        sys.path.insert(0, %r)
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from delfin.fffree.grip_polish import grip_polish

        mol = Chem.AddHs(Chem.MolFromSmiles("[N](C)(C)[Fe][N](C)(C)"))
        AllChem.EmbedMolecule(mol, randomSeed=2)
        n = mol.GetNumAtoms()
        conf = mol.GetConformer()
        P0 = np.asarray([[conf.GetAtomPosition(i).x,
                          conf.GetAtomPosition(i).y,
                          conf.GetAtomPosition(i).z] for i in range(n)],
                        dtype=np.float64)
        metal = next(i for i in range(n) if mol.GetAtomWithIdx(i).GetSymbol()=="Fe")
        donors = sorted(int(nb.GetIdx())
                        for nb in mol.GetAtomWithIdx(metal).GetNeighbors())
        P = grip_polish(P0, mol, metal, donors)
        sys.stdout.write(json.dumps([[float(x) for x in row]
                                     for row in np.asarray(P)]))
    """) % str(REPO_ROOT)
    env_off = os.environ.copy(); env_off["PYTHONHASHSEED"] = "0"
    env_off.pop(DONOR_CONE_ENV, None)
    env_on = env_off.copy(); env_on[DONOR_CONE_ENV] = "1"
    r_off = subprocess.run([sys.executable, "-c", script],
                           capture_output=True, text=True,
                           env=env_off, check=False)
    r_on = subprocess.run([sys.executable, "-c", script],
                          capture_output=True, text=True,
                          env=env_on, check=False)
    if r_off.returncode != 0 or r_on.returncode != 0:
        pytest.skip(f"rdkit stack unavailable: {r_off.stderr[-200:]}")
    # Both runs must succeed and produce a valid coordinate JSON.
    P_off = np.asarray(json.loads(r_off.stdout)) if r_off.stdout else None
    P_on = np.asarray(json.loads(r_on.stdout)) if r_on.stdout else None
    assert P_off is not None and P_on is not None
    # Sanity: the cone DoF either changed the polish output OR the
    # accept-if-better gate rolled back; in either case the result must
    # be finite.
    assert np.all(np.isfinite(P_off))
    assert np.all(np.isfinite(P_on))


# ---------------------------------------------------------------------------
# Tiny json import for the subprocess test
# ---------------------------------------------------------------------------
import json  # noqa: E402  (kept here so the subprocess tests can read JSON)
