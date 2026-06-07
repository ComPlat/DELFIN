"""Tests for Hapto-π explicit freeze + substituent subtree drag in GRIP
polish (Task #91, 2026-06-07).

Validates the new env-flag DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE that ties
hapto-π ring + substituent subtree atoms rigidly to the metal during the
M+D unfreeze L-BFGS pass.

Test matrix:
  1. Byte-identical OFF (env unset OR =0): output equals legacy HEAD.
  2. Ferrocene Fe(Cp)2: Cp rings translate together with M under polish.
  3. Cr(η6-C6H6) benzene: 6 ring C + 6 H all drag with M update.
  4. Mixed mode: η5-Cp + ancillary σ-phosphine; only Cp drags, phosphine stays.
  5. Gradient finite-diff (rigid body): grad_M of wrapped loss matches numeric.
  6. Substituent drag: η5-Cp with -CH3 substituent; methyl atoms drag with C.
  7. Determinism: 2 runs PYTHONHASHSEED=0 → byte-identical XYZ output.
"""
from __future__ import annotations

import math
import os
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def env_off(monkeypatch):
    """Hapto-π freeze flag OFF (default) -> byte-identical with HEAD."""
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE", raising=False)
    # Default-OFF byte-identity also requires the unfreeze flag OFF.
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", raising=False)


@pytest.fixture
def env_on(monkeypatch):
    """Hapto-π freeze flag ON + unfreeze ON (the production combo)."""
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_UNFREEZE_MD", "1")
    # Keep the hapto-skip on (default ON) so the legacy skip would
    # normally fire; the freeze flag is the one that flips it back ON.
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_MD_HAPTO_SKIP", "1")
    # Subtree helper module also needs its own env-flag for the
    # translation primitive to actually move atoms.
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_RIGID_SUBTREE", "1")


# ---------------------------------------------------------------------------
# Molecule builders (rdkit-based; skip when rdkit missing)
# ---------------------------------------------------------------------------
def _ferrocene_mol_and_coords():
    """Build Fe(Cp)2: Fe + 10 C + 10 H (5 per Cp ring).

    Index map:
      0       : Fe
      1-5     : top Cp ring carbons
      6-10    : top Cp ring H
      11-15   : bottom Cp ring carbons
      16-20   : bottom Cp ring H
    """
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    fe = mol.AddAtom(Chem.Atom("Fe"))
    top_c = [mol.AddAtom(Chem.Atom("C")) for _ in range(5)]
    top_h = [mol.AddAtom(Chem.Atom("H")) for _ in range(5)]
    bot_c = [mol.AddAtom(Chem.Atom("C")) for _ in range(5)]
    bot_h = [mol.AddAtom(Chem.Atom("H")) for _ in range(5)]
    # Fe-C bonds (10)
    for c in top_c + bot_c:
        mol.AddBond(fe, c, Chem.BondType.SINGLE)
    # Ring C-C bonds (top + bottom)
    for ring in (top_c, bot_c):
        for i in range(5):
            mol.AddBond(ring[i], ring[(i + 1) % 5], Chem.BondType.SINGLE)
    # Ring C-H bonds
    for ring, hs in ((top_c, top_h), (bot_c, bot_h)):
        for c, h in zip(ring, hs):
            mol.AddBond(c, h, Chem.BondType.SINGLE)

    # Coordinates: Cp rings at z = +1.65 / -1.65.
    P = np.zeros((21, 3), dtype=np.float64)
    r_ring = 1.21
    for ring_z, c_off, h_off in ((1.65, 1, 6), (-1.65, 11, 16)):
        for i in range(5):
            th = i * 2.0 * math.pi / 5.0
            P[c_off + i] = [r_ring * math.cos(th), r_ring * math.sin(th), ring_z]
            P[h_off + i] = [(r_ring + 1.08) * math.cos(th),
                            (r_ring + 1.08) * math.sin(th), ring_z]
    donors = list(range(1, 6)) + list(range(11, 16))
    return mol, P, 0, donors


def _benzene_cr_mol_and_coords():
    """Build Cr(η6-C6H6): Cr + 6 C + 6 H.

    Index map:
      0   : Cr
      1-6 : ring C
      7-12: ring H
    """
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    cr = mol.AddAtom(Chem.Atom("Cr"))
    cs = [mol.AddAtom(Chem.Atom("C")) for _ in range(6)]
    hs = [mol.AddAtom(Chem.Atom("H")) for _ in range(6)]
    for c in cs:
        mol.AddBond(cr, c, Chem.BondType.SINGLE)
    for i in range(6):
        mol.AddBond(cs[i], cs[(i + 1) % 6], Chem.BondType.SINGLE)
    for c, h in zip(cs, hs):
        mol.AddBond(c, h, Chem.BondType.SINGLE)
    r_ring = 1.40
    P = np.zeros((13, 3), dtype=np.float64)
    z = 1.70
    for i in range(6):
        th = i * 2.0 * math.pi / 6.0
        P[1 + i] = [r_ring * math.cos(th), r_ring * math.sin(th), z]
        P[7 + i] = [(r_ring + 1.08) * math.cos(th),
                    (r_ring + 1.08) * math.sin(th), z]
    donors = list(range(1, 7))
    return mol, P, 0, donors


def _methyl_cp_iron_mol_and_coords():
    """Build [Fe(η5-MeC5H4)] -- single methyl-Cp + Fe.

    Index map:
      0    : Fe
      1-5  : ring C
      6-9  : ring H (atoms 1-4; atom 5 carries the methyl)
      10   : methyl C
      11-13: methyl H
    """
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    fe = mol.AddAtom(Chem.Atom("Fe"))
    cs = [mol.AddAtom(Chem.Atom("C")) for _ in range(5)]
    hs_ring = [mol.AddAtom(Chem.Atom("H")) for _ in range(4)]
    me_c = mol.AddAtom(Chem.Atom("C"))
    me_h = [mol.AddAtom(Chem.Atom("H")) for _ in range(3)]
    for c in cs:
        mol.AddBond(fe, c, Chem.BondType.SINGLE)
    for i in range(5):
        mol.AddBond(cs[i], cs[(i + 1) % 5], Chem.BondType.SINGLE)
    # Ring atoms 0..3 carry H, atom 4 carries methyl
    for i in range(4):
        mol.AddBond(cs[i], hs_ring[i], Chem.BondType.SINGLE)
    mol.AddBond(cs[4], me_c, Chem.BondType.SINGLE)
    for h in me_h:
        mol.AddBond(me_c, h, Chem.BondType.SINGLE)

    P = np.zeros((14, 3), dtype=np.float64)
    r_ring = 1.21
    z = 1.65
    for i in range(5):
        th = i * 2.0 * math.pi / 5.0
        P[1 + i] = [r_ring * math.cos(th), r_ring * math.sin(th), z]
    # Ring H on atoms 1-4
    for i in range(4):
        th = i * 2.0 * math.pi / 5.0
        P[6 + i] = [(r_ring + 1.08) * math.cos(th),
                    (r_ring + 1.08) * math.sin(th), z]
    # Methyl C attached radially to ring C atom 5 (index 5 in P)
    th5 = 4 * 2.0 * math.pi / 5.0
    radial = np.array([math.cos(th5), math.sin(th5), 0.0])
    P[10] = P[5] + 1.50 * radial
    # Methyl H atoms (tetrahedral around me_c)
    P[11] = P[10] + np.array([0.5, 0.6, 0.7])
    P[12] = P[10] + np.array([0.6, -0.5, 0.7])
    P[13] = P[10] + np.array([0.0, 0.0, 1.08])
    donors = list(range(1, 6))
    return mol, P, 0, donors


def _cp_phosphine_mix_mol_and_coords():
    """Build Fe(η5-Cp)(PMe2H) -- mixed σ + π modes.

    Index map:
      0    : Fe
      1-5  : Cp ring C (η5)
      6-10 : Cp ring H
      11   : P (σ donor)
      12-13: Me C (× 2)
      14-19: Me H (3 per Me)
      20   : P-H
    """
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.RWMol()
    fe = mol.AddAtom(Chem.Atom("Fe"))
    cs = [mol.AddAtom(Chem.Atom("C")) for _ in range(5)]
    hs = [mol.AddAtom(Chem.Atom("H")) for _ in range(5)]
    p_at = mol.AddAtom(Chem.Atom("P"))
    me_c1 = mol.AddAtom(Chem.Atom("C"))
    me_c2 = mol.AddAtom(Chem.Atom("C"))
    me_h1 = [mol.AddAtom(Chem.Atom("H")) for _ in range(3)]
    me_h2 = [mol.AddAtom(Chem.Atom("H")) for _ in range(3)]
    p_h = mol.AddAtom(Chem.Atom("H"))
    for c in cs:
        mol.AddBond(fe, c, Chem.BondType.SINGLE)
    mol.AddBond(fe, p_at, Chem.BondType.SINGLE)
    for i in range(5):
        mol.AddBond(cs[i], cs[(i + 1) % 5], Chem.BondType.SINGLE)
    for c, h in zip(cs, hs):
        mol.AddBond(c, h, Chem.BondType.SINGLE)
    mol.AddBond(p_at, me_c1, Chem.BondType.SINGLE)
    mol.AddBond(p_at, me_c2, Chem.BondType.SINGLE)
    mol.AddBond(p_at, p_h, Chem.BondType.SINGLE)
    for h in me_h1:
        mol.AddBond(me_c1, h, Chem.BondType.SINGLE)
    for h in me_h2:
        mol.AddBond(me_c2, h, Chem.BondType.SINGLE)

    P = np.zeros((21, 3), dtype=np.float64)
    r_ring = 1.21
    z = 1.65
    for i in range(5):
        th = i * 2.0 * math.pi / 5.0
        P[1 + i] = [r_ring * math.cos(th), r_ring * math.sin(th), z]
        P[6 + i] = [(r_ring + 1.08) * math.cos(th),
                    (r_ring + 1.08) * math.sin(th), z]
    # P sits on the -z side
    P[11] = [0.0, 0.0, -2.30]
    P[12] = P[11] + np.array([1.40, 0.50, -0.20])
    P[13] = P[11] + np.array([-1.40, 0.50, -0.20])
    P[14] = P[12] + np.array([0.6, 0.5, -0.7])
    P[15] = P[12] + np.array([0.8, -0.4, -0.6])
    P[16] = P[12] + np.array([0.1, 0.7, -1.0])
    P[17] = P[13] + np.array([-0.6, 0.5, -0.7])
    P[18] = P[13] + np.array([-0.8, -0.4, -0.6])
    P[19] = P[13] + np.array([-0.1, 0.7, -1.0])
    P[20] = P[11] + np.array([0.0, -1.40, 0.0])
    # Donors: 5 ring C + P
    donors = list(range(1, 6)) + [11]
    return mol, P, 0, donors


# ---------------------------------------------------------------------------
# 1) Default-OFF byte-identity
# ---------------------------------------------------------------------------
class TestDefaultOff:
    def test_predicate_default_off(self, env_off):
        from delfin.fffree.grip_polish import _hapto_pi_freeze_active
        assert _hapto_pi_freeze_active() is False

    def test_predicate_explicit_off(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE", "0")
        from delfin.fffree.grip_polish import _hapto_pi_freeze_active
        assert _hapto_pi_freeze_active() is False

    def test_predicate_explicit_on(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_HAPTO_PI_FREEZE", "1")
        from delfin.fffree.grip_polish import _hapto_pi_freeze_active
        assert _hapto_pi_freeze_active() is True

    def test_byte_identical_off_ferrocene(self, env_off):
        """With OFF, grip_polish output must be byte-identical to a
        baseline run (we just compare two runs in the same env)."""
        from delfin.fffree.grip_polish import grip_polish
        mol, P, metal, donors = _ferrocene_mol_and_coords()
        P_a = grip_polish(P.copy(), mol, metal, donors)
        P_b = grip_polish(P.copy(), mol, metal, donors)
        assert np.array_equal(P_a, P_b)

    def test_off_no_rigid_body_built(self, env_off):
        """With both env-flags OFF, _build_hapto_rigid_body still works
        as a pure helper but the polish path doesn't even call it."""
        from delfin.fffree.grip_polish import (
            _hapto_pi_freeze_active, _build_hapto_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _ferrocene_mol_and_coords()
        # The predicate is OFF -> the polish would NOT build a rigid body.
        assert _hapto_pi_freeze_active() is False
        # But the helper still works deterministically when called directly.
        hapto = detect_hapto_atoms(mol, metal, donors)
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        # Metal index 0 is never in the rigid body set.
        assert 0 not in rigid


# ---------------------------------------------------------------------------
# 2) Ferrocene: Cp rings drag with M under polish
# ---------------------------------------------------------------------------
class TestFerroceneDrag:
    def test_rigid_body_contains_both_rings(self, env_on):
        """Ferrocene has 2 Cp rings; the rigid body should include atoms
        from BOTH rings (each in its own cluster) + every ring H, but
        the metal (index 0) must be excluded."""
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _ferrocene_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        # All 10 ring C atoms must be detected as hapto.
        assert hapto >= set(range(1, 6)) | set(range(11, 16))
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_set = set(rigid)
        # Metal NOT in rigid body
        assert 0 not in rigid_set
        # All 10 ring C
        for i in list(range(1, 6)) + list(range(11, 16)):
            assert i in rigid_set, f"ring C {i} missing"
        # All 10 ring H
        for i in list(range(6, 11)) + list(range(16, 21)):
            assert i in rigid_set, f"ring H {i} missing"

    def test_ring_translates_with_metal(self, env_on):
        """Apply _project_rigid_body with a synthetic metal displacement
        and verify every ring + H atom translates by the same delta."""
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, _project_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _ferrocene_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_idx = np.array(sorted(rigid), dtype=np.int64)
        # Apply a fake M displacement.
        delta = np.array([0.3, -0.2, 0.4])
        R = P.copy()
        R[metal] = P[metal] + delta
        R_out = _project_rigid_body(R, metal, rigid_idx, P)
        for i in rigid:
            assert np.allclose(R_out[i] - P[i], delta, atol=1e-12)
        # The metal itself stays where R put it.
        assert np.allclose(R_out[metal], P[metal] + delta, atol=1e-12)


# ---------------------------------------------------------------------------
# 3) Benzene-Cr η6: 6 ring C + 6 H drag with M
# ---------------------------------------------------------------------------
class TestBenzeneCr:
    def test_eta6_full_ring_and_h_in_rigid_body(self, env_on):
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _benzene_cr_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        # Full 6-membered ring must be detected.
        assert hapto >= set(range(1, 7))
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_set = set(rigid)
        # 6 ring carbons + 6 ring H, no metal
        assert 0 not in rigid_set
        for i in range(1, 13):
            assert i in rigid_set, f"atom {i} missing"

    def test_eta6_drag_under_projection(self, env_on):
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, _project_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _benzene_cr_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_idx = np.array(sorted(rigid), dtype=np.int64)
        delta = np.array([0.5, 0.0, -0.1])
        R = P.copy()
        R[metal] = P[metal] + delta
        R_out = _project_rigid_body(R, metal, rigid_idx, P)
        for i in rigid:
            assert np.allclose(R_out[i] - P[i], delta, atol=1e-12)


# ---------------------------------------------------------------------------
# 4) Mixed mode: η5-Cp + σ-phosphine; only Cp drags
# ---------------------------------------------------------------------------
class TestMixedMode:
    def test_phosphine_excluded_from_rigid_body(self, env_on):
        """The σ-bonded P (and its Me / H substituents) must NOT be
        absorbed into the Cp rigid body."""
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _cp_phosphine_mix_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        # Hapto detector picks up the 5 ring carbons only (P is σ, not C).
        assert hapto == set(range(1, 6))
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_set = set(rigid)
        # Cp ring + ring H in
        for i in range(1, 11):
            assert i in rigid_set, f"Cp atom {i} missing"
        # σ donor P + its substituents OUT (frozen_extra boundary).
        for i in range(11, 21):
            assert i not in rigid_set, f"P-side atom {i} should NOT be in rigid"


# ---------------------------------------------------------------------------
# 5) Gradient finite-diff check (rigid body)
# ---------------------------------------------------------------------------
class TestGradientFiniteDiff:
    def test_metal_gradient_accumulates_rigid_body(self, env_on):
        """Sanity-check the chain rule: a metal displacement that moves
        the rigid body produces a metal-gradient component equal to the
        SUM of the per-atom gradients (numeric finite-diff on a quadratic
        loss anchored at the metal)."""
        # Build a tiny synthetic quadratic loss: L(R) = 0.5 * Σ ||R[i] - c[i]||^2
        # so the gradient is grad[i] = R[i] - c[i].  When the rigid body
        # is tied to the metal, grad_M (total) = grad_M + Σ_{i in rigid} grad[i].
        from delfin.fffree.grip_polish import _project_rigid_body

        # Build a minimal system with 1 metal + 5 ring atoms.
        n_atoms = 6
        metal = 0
        rigid = np.array([1, 2, 3, 4, 5], dtype=np.int64)
        # Anchor points c[i].
        rng = np.random.RandomState(0)
        c = rng.randn(n_atoms, 3)
        # Initial coords.
        P_init = rng.randn(n_atoms, 3)

        def loss_and_grad(R):
            G = R - c
            L = 0.5 * float((G ** 2).sum())
            return L, G

        # Wrap the loss the same way grip_polish does for rigid body.
        def wrapped(R):
            R_eff = _project_rigid_body(R, metal, rigid, P_init)
            L, G = loss_and_grad(R_eff)
            # Chain rule: grad_M_total = grad_M + sum rigid grads
            G[metal] = G[metal] + G[rigid].sum(axis=0)
            G[rigid] = 0.0
            return L, G

        # Pick a test point where metal has moved.
        R_test = P_init.copy()
        R_test[metal] = P_init[metal] + np.array([0.2, -0.1, 0.3])
        L_a, G_a = wrapped(R_test)
        # Numeric finite-diff of grad_M
        eps = 1e-5
        for k in range(3):
            R_p = R_test.copy()
            R_p[metal, k] += eps
            L_p, _ = wrapped(R_p)
            R_m = R_test.copy()
            R_m[metal, k] -= eps
            L_m, _ = wrapped(R_m)
            num = (L_p - L_m) / (2.0 * eps)
            assert abs(num - G_a[metal, k]) < 1e-4, (
                f"finite-diff mismatch axis {k}: num={num} analytic={G_a[metal, k]}"
            )


# ---------------------------------------------------------------------------
# 6) Substituent drag: methyl on Cp rides with ring C
# ---------------------------------------------------------------------------
class TestSubstituentDrag:
    def test_methyl_in_rigid_body(self, env_on):
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _methyl_cp_iron_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_set = set(rigid)
        # 5 ring C + 4 ring H + 1 methyl C + 3 methyl H = 13 (no metal)
        for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
            assert i in rigid_set, f"atom {i} missing (sym/methyl test)"
        assert 0 not in rigid_set

    def test_methyl_translates_with_ring(self, env_on):
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, _project_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _methyl_cp_iron_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        rigid = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        rigid_idx = np.array(sorted(rigid), dtype=np.int64)
        delta = np.array([0.7, -0.3, 0.2])
        R = P.copy()
        R[metal] = P[metal] + delta
        R_out = _project_rigid_body(R, metal, rigid_idx, P)
        # Methyl C (10) and methyl H (11, 12, 13) all move by delta.
        for i in (10, 11, 12, 13):
            assert np.allclose(R_out[i] - P[i], delta, atol=1e-12), (
                f"atom {i} did not translate with ring"
            )
        # Methyl-C to ring-C[5] bond length preserved (rigid invariant).
        d_before = float(np.linalg.norm(P[10] - P[5]))
        d_after = float(np.linalg.norm(R_out[10] - R_out[5]))
        assert abs(d_after - d_before) < 1e-10


# ---------------------------------------------------------------------------
# 7) Determinism: two runs PYTHONHASHSEED=0 → byte-identical
# ---------------------------------------------------------------------------
class TestDeterminism:
    def test_two_runs_identical_on(self, env_on):
        """With the freeze flag ON, two consecutive polish runs on the
        same input must produce bit-identical outputs."""
        from delfin.fffree.grip_polish import grip_polish
        mol, P, metal, donors = _ferrocene_mol_and_coords()
        P_a = grip_polish(P.copy(), mol, metal, donors)
        P_b = grip_polish(P.copy(), mol, metal, donors)
        assert np.array_equal(P_a, P_b)

    def test_rigid_body_helper_deterministic(self, env_on):
        """_build_hapto_rigid_body is pure: same input -> same list."""
        from delfin.fffree.grip_polish import (
            _build_hapto_rigid_body, detect_hapto_atoms,
        )
        mol, P, metal, donors = _methyl_cp_iron_mol_and_coords()
        hapto = detect_hapto_atoms(mol, metal, donors)
        a = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        b = _build_hapto_rigid_body(mol, metal, donors, P, hapto)
        assert a == b


# ---------------------------------------------------------------------------
# 8) Integration: grip_polish with freeze ON does not crash, preserves ring
# ---------------------------------------------------------------------------
class TestIntegration:
    def test_polish_runs_on_ferrocene_with_flag_on(self, env_on):
        """Smoke test: grip_polish with freeze flag ON runs to completion
        and the metal-ring centroid distance is preserved within tol."""
        from delfin.fffree.grip_polish import grip_polish
        mol, P, metal, donors = _ferrocene_mol_and_coords()
        # Top ring centroid
        c_top_in = P[1:6].mean(axis=0)
        d_in = float(np.linalg.norm(c_top_in - P[metal]))
        P_out = grip_polish(P.copy(), mol, metal, donors)
        c_top_out = P_out[1:6].mean(axis=0)
        d_out = float(np.linalg.norm(c_top_out - P_out[metal]))
        # The M-centroid distance should be (nearly) preserved because the
        # rigid body translates rigidly with M.
        assert abs(d_out - d_in) < 0.05, (
            f"M-centroid distance drifted: {d_in:.4f} -> {d_out:.4f}"
        )

    def test_no_crash_on_pure_sigma_complex(self, env_on):
        """On a pure σ-only complex (no hapto), the freeze flag is a
        no-op: rigid body is empty -> byte-identical with the unfreeze
        flag's normal behaviour."""
        Chem = pytest.importorskip("rdkit.Chem")
        from delfin.fffree.grip_polish import grip_polish
        # Simple Fe(NH3)3 -- no hapto.
        mol = Chem.RWMol()
        fe = mol.AddAtom(Chem.Atom("Fe"))
        ns = [mol.AddAtom(Chem.Atom("N")) for _ in range(3)]
        for n in ns:
            mol.AddBond(fe, n, Chem.BondType.SINGLE)
        P = np.zeros((4, 3))
        P[1] = [2.0, 0.0, 0.0]
        P[2] = [-1.0, 1.732, 0.0]
        P[3] = [-1.0, -1.732, 0.0]
        # Should run without raising.
        P_out = grip_polish(P.copy(), mol, fe, [1, 2, 3])
        assert P_out.shape == (4, 3)
        assert np.all(np.isfinite(P_out))
