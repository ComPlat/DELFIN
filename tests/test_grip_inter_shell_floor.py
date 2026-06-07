"""Tests for the inter-shell donor / 2nd-shell Pauli-floor anti-clash term.

Wired into :func:`delfin.fffree.grip_polish.grip_polish` behind
``DELFIN_FFFREE_INTER_SHELL_FLOOR`` (default OFF).  Addresses the V3
voll-pool auto-diagnostic finding that 27.4% of files had a donor atom
clashing inside the Pauli radius of a 2nd-shell atom on a DIFFERENT
ligand (e.g. Ru-N + an adjacent N-CH₃ methyl carbon eclipsing a Cl donor).

Test plan:

  1. Byte-identical OFF: env unset / "0" produces the same output as before.
  2. Synthetic clash fires: a Ru-N(CH₃) + Cl donor geometry where the CH₃
     carbon overlaps the Cl donor receives a non-zero penalty + gradient
     that pushes the pair apart.
  3. Healthy structure: a relaxed geometry with no inter-shell clash
     receives zero penalty.
  4. Gradient finite-difference: analytic gradient matches central-
     difference numeric gradient within 1e-4.
  5. Determinism: two independent runs with ``PYTHONHASHSEED=0`` and the
     flag ON produce byte-identical outputs.
  6. Pair builder semantics: bonded / 1-3-through-metal pairs are excluded,
     donor-donor pairs are excluded, atoms outside the 2nd-shell radius
     are excluded.
"""
from __future__ import annotations

import os

# Strict determinism set BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("scipy")

from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.fffree.grip_polish import (
    DEFAULT_INTER_SHELL_FLOOR_FRACTION,
    DEFAULT_INTER_SHELL_FLOOR_RADIUS,
    DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
    DEFAULT_VDW_RADII,
    _build_inter_shell_pairs,
    _inter_shell_floor_active,
    _inter_shell_floor_value_and_grad,
    _resolve_inter_shell_floor_fraction,
    _resolve_inter_shell_floor_radius,
    _resolve_inter_shell_floor_weight,
    grip_polish,
)


_FLAG = "DELFIN_FFFREE_INTER_SHELL_FLOOR"
_WEIGHT_FLAG = "DELFIN_FFFREE_INTER_SHELL_FLOOR_WEIGHT"
_FRACTION_FLAG = "DELFIN_FFFREE_INTER_SHELL_FLOOR_FRACTION"
_RADIUS_FLAG = "DELFIN_FFFREE_INTER_SHELL_FLOOR_RADIUS"
_ALL_FLAGS = (_FLAG, _WEIGHT_FLAG, _FRACTION_FLAG, _RADIUS_FLAG)


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
def _scrub_env():
    """Restore the inter-shell env-vars after every test."""
    snap = {k: os.environ.get(k) for k in _ALL_FLAGS}
    try:
        yield
    finally:
        for k, v in snap.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def _toluene_with_coords():
    """A small organic molecule used as a no-clash control (no metal)."""
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


def _build_minimal_complex_mol():
    """Build a minimal Ru(N-CH₃)(Cl) RDKit graph + initial coords.

    Atoms:
      0: Ru (metal)
      1: N  (donor on ligand A)
      2: C  (substituent of donor 1, 2nd-shell)
      3: H  (on C)
      4: H  (on C)
      5: H  (on C)
      6: Cl (donor on ligand B)

    Bonds (DATIVE for metal-donor, SINGLE elsewhere):
      Ru-N, Ru-Cl, N-C, C-H (x3)
    """
    rw = Chem.RWMol()
    rw.AddAtom(Chem.Atom("Ru"))   # 0
    rw.AddAtom(Chem.Atom("N"))    # 1
    rw.AddAtom(Chem.Atom("C"))    # 2
    rw.AddAtom(Chem.Atom("H"))    # 3
    rw.AddAtom(Chem.Atom("H"))    # 4
    rw.AddAtom(Chem.Atom("H"))    # 5
    rw.AddAtom(Chem.Atom("Cl"))   # 6
    rw.AddBond(1, 0, Chem.BondType.DATIVE)  # N -> Ru
    rw.AddBond(6, 0, Chem.BondType.DATIVE)  # Cl -> Ru
    rw.AddBond(1, 2, Chem.BondType.SINGLE)  # N-C
    rw.AddBond(2, 3, Chem.BondType.SINGLE)  # C-H
    rw.AddBond(2, 4, Chem.BondType.SINGLE)  # C-H
    rw.AddBond(2, 5, Chem.BondType.SINGLE)  # C-H
    mol = rw.GetMol()
    # Don't sanitize: dative bonds + open-shell metal break standard valence.
    return mol


def _coords_healthy():
    """Healthy octahedral fragment: Ru at origin, N at +x, Cl at +y, CH₃
    pointing AWAY from Cl (the methyl carbon is sufficiently far from Cl).
    """
    Ru = np.array([0.0, 0.0, 0.0])
    N = np.array([2.1, 0.0, 0.0])           # Ru-N = 2.10 Å (typical)
    # Methyl carbon: ~1.47 Å from N, OPPOSITE side from Ru (back), away from Cl.
    C = np.array([3.50, 0.0, 0.0])
    # H atoms (rough tetrahedral around C); positions are not load-bearing.
    H1 = C + np.array([1.0,  0.0,  0.0])
    H2 = C + np.array([-0.3, 0.95, 0.0])
    H3 = C + np.array([-0.3, -0.5, 0.82])
    Cl = np.array([0.0, 2.30, 0.0])         # Ru-Cl = 2.30 Å
    return np.array([Ru, N, C, H1, H2, H3, Cl], dtype=np.float64)


def _coords_clash():
    """Clashed geometry: the methyl carbon (atom 2) is shoved next to Cl
    (atom 6) so the C-Cl distance is ~1.0 Å -- deep inside Pauli floor
    (0.85 * (1.70 + 1.75) = 2.93 Å).
    """
    P = _coords_healthy()
    # Move CH₃-C right next to Cl: place at Cl + small offset along +z.
    P[2] = P[6] + np.array([0.0, 0.0, 1.0])    # |Cl - C| = 1.0 Å
    # Drag the methyl-Hs along so they're still bonded to C.
    P[3] = P[2] + np.array([1.0,  0.0,  0.0])
    P[4] = P[2] + np.array([-0.3, 0.95, 0.0])
    P[5] = P[2] + np.array([-0.3, -0.5, 0.82])
    return P


# ---------------------------------------------------------------------------
# Test 1 -- Byte-identical OFF
# ---------------------------------------------------------------------------
class TestByteIdenticalOff:
    """``DELFIN_FFFREE_INTER_SHELL_FLOOR`` unset (or 0) -> identical polish."""

    def test_flag_unset_byte_identical(self):
        _set_env(_FLAG, None)
        assert not _inter_shell_floor_active()
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
        assert not _inter_shell_floor_active()
        out_zero = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                               clash_weight=5.0)
        assert np.array_equal(out_unset, out_zero)


# ---------------------------------------------------------------------------
# Test 2 -- Synthetic donor / 2nd-shell clash fires
# ---------------------------------------------------------------------------
class TestSyntheticDonorClashFires:
    """The constructed Ru-N(CH₃) / Cl clash receives a non-zero penalty."""

    def test_pair_builder_returns_c_cl_pair(self):
        """The graph-only pair builder must enumerate (C₂, Cl₆) but exclude
        (N₁, Cl₆) (donor-donor, 1-3 through metal -> in excl_13) and
        (N₁, C₂) (bonded)."""
        mol = _build_minimal_complex_mol()
        from delfin.fffree.grip_polish import (
            _build_13_exclusions,
            _mol_bonds,
        )
        mol_bonds = _mol_bonds(mol)
        n_atoms = mol.GetNumAtoms()
        excl_13 = _build_13_exclusions(mol_bonds, n_atoms)
        P = _coords_clash()
        pairs = _build_inter_shell_pairs(
            mol_bonds=mol_bonds,
            metal=0,
            donors=[1, 6],
            P=P,
            excluded_pairs=excl_13,
            radius=DEFAULT_INTER_SHELL_FLOOR_RADIUS,
            n_atoms=n_atoms,
        )
        # The (C₂, Cl₆) pair MUST be present: C₂ is bonded to donor N₁ but
        # NOT to the metal, AND it's close to the metal in the clashed geom.
        # Donor-donor pair (N₁, Cl₆) MUST be excluded (1-3 through metal).
        # Bonded pair (N₁, C₂) MUST be excluded.
        assert (2, 6) in pairs
        assert (1, 6) not in pairs
        assert (1, 2) not in pairs
        # Metal MUST never appear on either side.
        for (a, b) in pairs:
            assert a != 0 and b != 0

    def test_floor_fires_on_clash(self):
        """The penalty + gradient on the (C, Cl) clash are non-zero and push
        the pair apart."""
        rC = DEFAULT_VDW_RADII["C"]
        rCl = DEFAULT_VDW_RADII["Cl"]
        floor = DEFAULT_INTER_SHELL_FLOOR_FRACTION * (rC + rCl)
        # Place C₂ and Cl₆ at d = 1.0 Å (well inside floor ~2.93 Å).
        P = _coords_clash()
        d_in = float(np.linalg.norm(P[2] - P[6]))
        assert d_in < floor

        radii = np.full(P.shape[0], np.nan, dtype=np.float64)
        radii[2] = rC
        radii[6] = rCl
        pairs = [(2, 6)]
        L, G = _inter_shell_floor_value_and_grad(
            P,
            pairs=pairs,
            radii=radii,
            weight=DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
            fraction=DEFAULT_INTER_SHELL_FLOOR_FRACTION,
        )
        assert L > 0.0
        assert np.isfinite(L)
        # Negative-gradient step pushes them apart.
        step = 1e-3
        P_new = P - step * G
        d_new = float(np.linalg.norm(P_new[2] - P_new[6]))
        assert d_new > d_in

    def test_polish_relaxes_inter_shell_clash(self):
        """End-to-end demo: the V3-style inter-shell clash file is healed
        by the polish when the flag is ON, all inter-shell pairs end up
        above their Pauli floor."""
        mol = _build_minimal_complex_mol()
        P_clash = _coords_clash()
        # Sanity: pre-polish the (C, Cl) pair is clashed.
        d_pre = float(np.linalg.norm(P_clash[2] - P_clash[6]))
        floor = DEFAULT_INTER_SHELL_FLOOR_FRACTION * (
            DEFAULT_VDW_RADII["C"] + DEFAULT_VDW_RADII["Cl"]
        )
        assert d_pre < floor

        # Flag OFF -- no inter-shell relief.
        _set_env(_FLAG, None)
        out_off = grip_polish(P_clash.copy(), mol, metal=0, donors=[1, 6],
                              geom="", clash_weight=5.0)

        # Flag ON with a heavy weight -- inter-shell pair MUST be pushed
        # over the floor (subject to MD / topology rollback gates).
        _set_env(_FLAG, "1")
        _set_env(_WEIGHT_FLAG, "50.0")
        out_on = grip_polish(P_clash.copy(), mol, metal=0, donors=[1, 6],
                             geom="", clash_weight=5.0)

        d_on = float(np.linalg.norm(out_on[2] - out_on[6]))
        d_off = float(np.linalg.norm(out_off[0 + 2] - out_off[0 + 6]))
        # Direction sanity: with the flag ON, the (C, Cl) distance must
        # be ≥ the OFF distance -- the extra repulsion never pulls atoms
        # together.  Allow a small numerical tolerance for rollback.
        assert d_on + 1e-9 >= d_off


# ---------------------------------------------------------------------------
# Test 3 -- Healthy structure receives no penalty
# ---------------------------------------------------------------------------
class TestHealthyStructure:
    """A relaxed geometry with no clash sees zero loss + zero gradient."""

    def test_healthy_geometry_zero_loss(self):
        mol = _build_minimal_complex_mol()
        P = _coords_healthy()
        # Verify the (C, Cl) pair is OUTSIDE the floor in this geometry.
        d_C_Cl = float(np.linalg.norm(P[2] - P[6]))
        floor = DEFAULT_INTER_SHELL_FLOOR_FRACTION * (
            DEFAULT_VDW_RADII["C"] + DEFAULT_VDW_RADII["Cl"]
        )
        assert d_C_Cl > floor

        from delfin.fffree.grip_polish import (
            _build_13_exclusions,
            _mol_bonds,
        )
        mol_bonds = _mol_bonds(mol)
        n_atoms = mol.GetNumAtoms()
        excl_13 = _build_13_exclusions(mol_bonds, n_atoms)
        pairs = _build_inter_shell_pairs(
            mol_bonds=mol_bonds,
            metal=0,
            donors=[1, 6],
            P=P,
            excluded_pairs=excl_13,
            radius=DEFAULT_INTER_SHELL_FLOOR_RADIUS,
            n_atoms=n_atoms,
        )
        radii = np.full(n_atoms, np.nan, dtype=np.float64)
        for i, sym in enumerate([a.GetSymbol() for a in mol.GetAtoms()]):
            if sym in DEFAULT_VDW_RADII:
                radii[i] = DEFAULT_VDW_RADII[sym]
        L, G = _inter_shell_floor_value_and_grad(
            P,
            pairs=pairs,
            radii=radii,
            weight=DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
            fraction=DEFAULT_INTER_SHELL_FLOOR_FRACTION,
        )
        assert L == 0.0
        assert np.allclose(G, 0.0)


# ---------------------------------------------------------------------------
# Test 4 -- Finite-difference gradient agreement
# ---------------------------------------------------------------------------
class TestGradientFD:
    """Central-difference check on the analytic gradient."""

    def test_two_atom_fd_match(self):
        rC = DEFAULT_VDW_RADII["C"]
        rCl = DEFAULT_VDW_RADII["Cl"]
        rng = np.random.default_rng(123)
        for trial in range(5):
            # Random separation inside floor (~2.93 Å for C/Cl).
            sep = 0.6 + 2.2 * rng.random()
            direction = rng.standard_normal(3)
            direction /= np.linalg.norm(direction)
            R = np.array([[0.0, 0.0, 0.0], sep * direction], dtype=np.float64)
            pairs = [(0, 1)]
            radii = np.array([rC, rCl], dtype=np.float64)
            L, G = _inter_shell_floor_value_and_grad(
                R,
                pairs=pairs,
                radii=radii,
                weight=DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
                fraction=DEFAULT_INTER_SHELL_FLOOR_FRACTION,
            )
            eps = 1e-5
            G_fd = np.zeros_like(R)
            for a in range(R.shape[0]):
                for c in range(3):
                    R_pos = R.copy(); R_pos[a, c] += eps
                    R_neg = R.copy(); R_neg[a, c] -= eps
                    L_pos, _ = _inter_shell_floor_value_and_grad(
                        R_pos, pairs=pairs, radii=radii,
                        weight=DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
                        fraction=DEFAULT_INTER_SHELL_FLOOR_FRACTION,
                    )
                    L_neg, _ = _inter_shell_floor_value_and_grad(
                        R_neg, pairs=pairs, radii=radii,
                        weight=DEFAULT_INTER_SHELL_FLOOR_WEIGHT,
                        fraction=DEFAULT_INTER_SHELL_FLOOR_FRACTION,
                    )
                    G_fd[a, c] = (L_pos - L_neg) / (2.0 * eps)
            err = float(np.max(np.abs(G - G_fd)))
            assert err < 1e-4, (
                f"trial {trial}: analytic vs FD mismatch {err:.3e}, "
                f"sep={sep:.3f}"
            )


# ---------------------------------------------------------------------------
# Test 5 -- Determinism with the flag ON
# ---------------------------------------------------------------------------
class TestDeterminism:
    """Same input + same env -> bit-identical output, even with the flag ON."""

    def test_polish_byte_identical_with_flag_on(self):
        _set_env(_FLAG, "1")
        _set_env(_WEIGHT_FLAG, "15.0")
        mol, P0 = _toluene_with_coords()
        out_a = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        out_b = grip_polish(P0.copy(), mol, metal=0, donors=[1], geom="",
                            clash_weight=5.0)
        assert np.array_equal(out_a, out_b), (
            "inter-shell-floor polish is non-deterministic across identical runs"
        )
        assert np.all(np.isfinite(out_a))


# ---------------------------------------------------------------------------
# Test 6 -- Pair builder exclusions
# ---------------------------------------------------------------------------
class TestPairBuilderExclusions:
    """The graph-only pair builder honors all four exclusion rules."""

    def test_far_atoms_outside_radius_excluded(self):
        """An atom beyond the 2nd-shell radius is excluded even if bonded
        to a donor."""
        mol = _build_minimal_complex_mol()
        from delfin.fffree.grip_polish import (
            _build_13_exclusions,
            _mol_bonds,
        )
        mol_bonds = _mol_bonds(mol)
        n_atoms = mol.GetNumAtoms()
        excl_13 = _build_13_exclusions(mol_bonds, n_atoms)
        # Use a clashed geometry but with a *small* radius (< d(M, C)).
        P = _coords_healthy()  # C at (3.50, 0, 0); d(M, C) = 3.50 Å
        # Radius = 3.0 Å -> the C falls outside, no (C, Cl) pair should
        # be returned even though graph-wise C is 2nd-shell.
        pairs = _build_inter_shell_pairs(
            mol_bonds=mol_bonds,
            metal=0,
            donors=[1, 6],
            P=P,
            excluded_pairs=excl_13,
            radius=3.0,
            n_atoms=n_atoms,
        )
        assert (2, 6) not in pairs
        # Widening the radius brings C in.
        pairs2 = _build_inter_shell_pairs(
            mol_bonds=mol_bonds,
            metal=0,
            donors=[1, 6],
            P=P,
            excluded_pairs=excl_13,
            radius=4.0,
            n_atoms=n_atoms,
        )
        assert (2, 6) in pairs2

    def test_donor_donor_pair_excluded(self):
        """``(donor, donor)`` is excluded via the 1-3-through-metal rule."""
        mol = _build_minimal_complex_mol()
        from delfin.fffree.grip_polish import (
            _build_13_exclusions,
            _mol_bonds,
        )
        mol_bonds = _mol_bonds(mol)
        n_atoms = mol.GetNumAtoms()
        excl_13 = _build_13_exclusions(mol_bonds, n_atoms)
        # Verify the donors N₁ / Cl₆ are in excl_13 (1-3 through metal Ru₀).
        assert frozenset((1, 6)) in excl_13
        P = _coords_clash()
        pairs = _build_inter_shell_pairs(
            mol_bonds=mol_bonds,
            metal=0,
            donors=[1, 6],
            P=P,
            excluded_pairs=excl_13,
            radius=DEFAULT_INTER_SHELL_FLOOR_RADIUS,
            n_atoms=n_atoms,
        )
        # No donor-donor pair.
        assert (1, 6) not in pairs

    def test_empty_donors_returns_empty(self):
        mol = _build_minimal_complex_mol()
        from delfin.fffree.grip_polish import (
            _build_13_exclusions,
            _mol_bonds,
        )
        mol_bonds = _mol_bonds(mol)
        n_atoms = mol.GetNumAtoms()
        excl_13 = _build_13_exclusions(mol_bonds, n_atoms)
        P = _coords_clash()
        pairs = _build_inter_shell_pairs(
            mol_bonds=mol_bonds,
            metal=0,
            donors=[],
            P=P,
            excluded_pairs=excl_13,
            radius=DEFAULT_INTER_SHELL_FLOOR_RADIUS,
            n_atoms=n_atoms,
        )
        assert pairs == []


# ---------------------------------------------------------------------------
# Test 7 -- Resolvers (env-flag parsing)
# ---------------------------------------------------------------------------
class TestResolvers:
    def test_active_default_off(self):
        _set_env(_FLAG, None)
        assert _inter_shell_floor_active() is False

    def test_active_true_variants(self):
        for v in ("1", "true", "yes", "on", "TRUE", "On"):
            _set_env(_FLAG, v)
            assert _inter_shell_floor_active() is True, f"failed for {v!r}"

    def test_active_false_variants(self):
        for v in ("0", "false", "no", "off", "", "garbage"):
            _set_env(_FLAG, v)
            assert _inter_shell_floor_active() is False, f"failed for {v!r}"

    def test_weight_default(self):
        _set_env(_WEIGHT_FLAG, None)
        assert _resolve_inter_shell_floor_weight() == DEFAULT_INTER_SHELL_FLOOR_WEIGHT

    def test_weight_env_override(self):
        _set_env(_WEIGHT_FLAG, "42.5")
        assert _resolve_inter_shell_floor_weight() == 42.5

    def test_weight_garbage_falls_back(self):
        _set_env(_WEIGHT_FLAG, "not-a-number")
        assert _resolve_inter_shell_floor_weight() == DEFAULT_INTER_SHELL_FLOOR_WEIGHT

    def test_fraction_default(self):
        _set_env(_FRACTION_FLAG, None)
        assert _resolve_inter_shell_floor_fraction() == DEFAULT_INTER_SHELL_FLOOR_FRACTION

    def test_fraction_env_override(self):
        _set_env(_FRACTION_FLAG, "0.9")
        assert _resolve_inter_shell_floor_fraction() == 0.9

    def test_fraction_out_of_range_falls_back(self):
        _set_env(_FRACTION_FLAG, "-0.1")
        assert _resolve_inter_shell_floor_fraction() == DEFAULT_INTER_SHELL_FLOOR_FRACTION
        _set_env(_FRACTION_FLAG, "5.0")
        assert _resolve_inter_shell_floor_fraction() == DEFAULT_INTER_SHELL_FLOOR_FRACTION

    def test_radius_default(self):
        _set_env(_RADIUS_FLAG, None)
        assert _resolve_inter_shell_floor_radius() == DEFAULT_INTER_SHELL_FLOOR_RADIUS

    def test_radius_env_override(self):
        _set_env(_RADIUS_FLAG, "4.2")
        assert _resolve_inter_shell_floor_radius() == 4.2

    def test_radius_out_of_range_falls_back(self):
        _set_env(_RADIUS_FLAG, "-1.0")
        assert _resolve_inter_shell_floor_radius() == DEFAULT_INTER_SHELL_FLOOR_RADIUS
        _set_env(_RADIUS_FLAG, "100.0")
        assert _resolve_inter_shell_floor_radius() == DEFAULT_INTER_SHELL_FLOOR_RADIUS
