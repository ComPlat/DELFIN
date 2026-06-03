"""Tests for the class-conditional σ-only GRIP polish mode (2026-06-03).

Validates ``DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE``:

* When the flag is ON, the hapto-protected atom set is expanded to cover
  the WHOLE π-system (aromatic / all-C rings carrying a C donor + H atoms
  on those rings); σ donors (N/O/S/P, isolated C) are NOT expanded.
* When the flag is OFF (default), the σ-only mode helper is NOT called
  and the polish behaves byte-identically with HEAD c03a550.
* The expansion logic is deterministic across calls (no RNG).
* On a structure with no hapto / π carbon-donor, the expansion is a
  no-op (returns the input set verbatim).
"""
from __future__ import annotations

import json
import os

import numpy as np
import pytest


# Defer rdkit import; tests that need RDKit pickle in via fixtures.
@pytest.fixture
def rdkit_chem():
    Chem = pytest.importorskip("rdkit.Chem")
    return Chem


# Patch in the worktree's delfin package via sys.path before importing.
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_WORKTREE_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _WORKTREE_ROOT not in sys.path:
    sys.path.insert(0, _WORKTREE_ROOT)

from delfin.fffree.grip_polish import (  # noqa: E402
    _sigma_only_mode_active,
    detect_hapto_atoms,
    expand_hapto_for_sigma_only,
)


# ---------------------------------------------------------------------------
# 1) Env-flag activation contract
# ---------------------------------------------------------------------------
class TestSigmaOnlyModeFlag:
    def test_default_off_when_unset(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE", raising=False)
        assert _sigma_only_mode_active() is False

    def test_explicit_off(self, monkeypatch):
        for val in ("0", "false", "False", "no", "off", ""):
            monkeypatch.setenv("DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE", val)
            assert _sigma_only_mode_active() is False, f"failed for {val!r}"

    def test_explicit_on(self, monkeypatch):
        for val in ("1", "true", "True", "yes", "on", "YES", "True"):
            monkeypatch.setenv("DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE", val)
            assert _sigma_only_mode_active() is True, f"failed for {val!r}"

    def test_garbage_value_is_off(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE", "maybe")
        assert _sigma_only_mode_active() is False


# ---------------------------------------------------------------------------
# 2) Ferrocene Cp rings are skipped — already by detect_hapto_atoms, but
#    σ-only-mode must NOT shrink that set (always a superset).
# ---------------------------------------------------------------------------
class TestExpandFerroceneCpRings:
    def test_cp_ring_already_in_base_stays_in_expanded(self, rdkit_chem):
        """Build a metal + 2 Cp rings (η5 each).  ``detect_hapto_atoms``
        catches both rings; ``expand_hapto_for_sigma_only`` must include
        them in the output AND also pick up the ring H atoms."""
        Chem = rdkit_chem
        from rdkit.Chem import RWMol
        # Build the molecule by hand: Fe + 2 Cp rings, with H on each
        # ring carbon.
        rw = RWMol()
        # Add Fe.
        fe = rw.AddAtom(Chem.Atom("Fe"))
        # Two Cp rings.
        cp_atoms_top = [rw.AddAtom(Chem.Atom("C")) for _ in range(5)]
        cp_atoms_bot = [rw.AddAtom(Chem.Atom("C")) for _ in range(5)]
        # Each ring C gets one explicit H.
        h_top = [rw.AddAtom(Chem.Atom("H")) for _ in range(5)]
        h_bot = [rw.AddAtom(Chem.Atom("H")) for _ in range(5)]
        # Bond the rings (aromatic).
        for k in range(5):
            rw.AddBond(cp_atoms_top[k], cp_atoms_top[(k + 1) % 5],
                       Chem.BondType.AROMATIC)
            rw.AddBond(cp_atoms_bot[k], cp_atoms_bot[(k + 1) % 5],
                       Chem.BondType.AROMATIC)
            rw.AddBond(cp_atoms_top[k], h_top[k], Chem.BondType.SINGLE)
            rw.AddBond(cp_atoms_bot[k], h_bot[k], Chem.BondType.SINGLE)
        # M-D edges (η5 sand η5).
        for c in cp_atoms_top + cp_atoms_bot:
            rw.AddBond(fe, c, Chem.BondType.SINGLE)
        # Set aromatic + flags so RDKit's SSSR picks up the 5-rings.
        for c in cp_atoms_top + cp_atoms_bot:
            rw.GetAtomWithIdx(c).SetIsAromatic(True)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except Exception:
            # Don't fail the test if sanitisation rejects our toy.
            pass
        donors = cp_atoms_top + cp_atoms_bot
        base = detect_hapto_atoms(mol, fe, donors)
        # Base set must already cover all 10 ring C's.
        assert set(base) >= set(donors)
        # Expanded set is a SUPERSET and includes the H atoms too.
        expanded = expand_hapto_for_sigma_only(mol, fe, donors, base)
        assert set(expanded) >= set(base)
        # Every ring H atom is in the expanded set.
        for h in h_top + h_bot:
            assert h in expanded, f"H atom {h} should be in σ-only expansion"

    def test_expansion_is_idempotent(self, rdkit_chem):
        """Calling expand twice with the second call's input = first call's
        output must return the same set (no further growth)."""
        Chem = rdkit_chem
        from rdkit.Chem import RWMol
        rw = RWMol()
        fe = rw.AddAtom(Chem.Atom("Fe"))
        ring = [rw.AddAtom(Chem.Atom("C")) for _ in range(5)]
        h = [rw.AddAtom(Chem.Atom("H")) for _ in range(5)]
        for k in range(5):
            rw.AddBond(ring[k], ring[(k + 1) % 5], Chem.BondType.AROMATIC)
            rw.AddBond(ring[k], h[k], Chem.BondType.SINGLE)
        for c in ring:
            rw.AddBond(fe, c, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(c).SetIsAromatic(True)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        donors = ring
        base = detect_hapto_atoms(mol, fe, donors)
        e1 = expand_hapto_for_sigma_only(mol, fe, donors, base)
        e2 = expand_hapto_for_sigma_only(mol, fe, donors, e1)
        assert e1 == e2


# ---------------------------------------------------------------------------
# 3) σ-N donor (pyridine N): N stays in the loss, NOT expanded as π.
# ---------------------------------------------------------------------------
class TestPyridineSigmaNDonorIncluded:
    def test_pyridine_N_donor_not_in_expansion(self, rdkit_chem):
        """A pyridine N coordinating to Fe via the lone pair is a σ-donor.
        The σ-only-mode expansion must NOT add the N to the protected set
        (i.e. the N keeps its CCDC bond/angle priors)."""
        Chem = rdkit_chem
        from rdkit.Chem import RWMol
        # Build Fe + pyridine: N at idx 1, C 2..6 around it, H on every C.
        rw = RWMol()
        fe = rw.AddAtom(Chem.Atom("Fe"))
        n_idx = rw.AddAtom(Chem.Atom("N"))
        cs = [rw.AddAtom(Chem.Atom("C")) for _ in range(5)]
        ring = [n_idx] + cs
        hs = [rw.AddAtom(Chem.Atom("H")) for _ in range(5)]
        # Ring bonds (aromatic).
        for k in range(6):
            rw.AddBond(ring[k], ring[(k + 1) % 6], Chem.BondType.AROMATIC)
            rw.GetAtomWithIdx(ring[k]).SetIsAromatic(True)
        for k in range(5):
            rw.AddBond(cs[k], hs[k], Chem.BondType.SINGLE)
        # M-N σ-bond.
        rw.AddBond(fe, n_idx, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        donors = [n_idx]
        # Base hapto set: empty (only 1 C donor and it's the σ-N).
        base = detect_hapto_atoms(mol, fe, donors)
        assert n_idx not in base
        # Expanded: still must not include the N (the donor is σ).
        expanded = expand_hapto_for_sigma_only(mol, fe, donors, base)
        # The pyridine RING (carbons + H) is NOT part of the donor list
        # and so the helper has no C-donor in any ring -> nothing to add.
        assert n_idx not in expanded
        # Ring carbons too should not be added (no donor C in the ring).
        for c in cs:
            assert c not in expanded
        # And the expanded set is identical to the base (empty here).
        assert expanded == base


# ---------------------------------------------------------------------------
# 4) Half-sandwich (η5-Cp + σ-Cl) — Cp ring + ring H added; Cl untouched.
# ---------------------------------------------------------------------------
class TestHalfSandwichExpansion:
    def test_only_pi_ring_and_H_added_not_sigma_Cl(self, rdkit_chem):
        Chem = rdkit_chem
        from rdkit.Chem import RWMol
        rw = RWMol()
        fe = rw.AddAtom(Chem.Atom("Fe"))
        cp = [rw.AddAtom(Chem.Atom("C")) for _ in range(5)]
        h_cp = [rw.AddAtom(Chem.Atom("H")) for _ in range(5)]
        cl_idx = rw.AddAtom(Chem.Atom("Cl"))
        for k in range(5):
            rw.AddBond(cp[k], cp[(k + 1) % 5], Chem.BondType.AROMATIC)
            rw.AddBond(cp[k], h_cp[k], Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(cp[k]).SetIsAromatic(True)
        for c in cp:
            rw.AddBond(fe, c, Chem.BondType.SINGLE)
        rw.AddBond(fe, cl_idx, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        donors = list(cp) + [cl_idx]
        base = detect_hapto_atoms(mol, fe, donors)
        # The Cp ring is detected at η5; Cl is not in the base.
        assert cl_idx not in base
        expanded = expand_hapto_for_sigma_only(mol, fe, donors, base)
        assert cl_idx not in expanded
        # Ring H atoms ARE added.
        for h in h_cp:
            assert h in expanded


# ---------------------------------------------------------------------------
# 5) Determinism: two calls must yield identical sets.
# ---------------------------------------------------------------------------
class TestDeterminism:
    def test_two_runs_identical_set(self, rdkit_chem):
        Chem = rdkit_chem
        from rdkit.Chem import RWMol
        rw = RWMol()
        fe = rw.AddAtom(Chem.Atom("Fe"))
        ring = [rw.AddAtom(Chem.Atom("C")) for _ in range(6)]
        hs = [rw.AddAtom(Chem.Atom("H")) for _ in range(6)]
        for k in range(6):
            rw.AddBond(ring[k], ring[(k + 1) % 6], Chem.BondType.AROMATIC)
            rw.AddBond(ring[k], hs[k], Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(ring[k]).SetIsAromatic(True)
        for c in ring:
            rw.AddBond(fe, c, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        donors = ring
        base = detect_hapto_atoms(mol, fe, donors)
        s1 = sorted(expand_hapto_for_sigma_only(mol, fe, donors, base))
        s2 = sorted(expand_hapto_for_sigma_only(mol, fe, donors, base))
        assert s1 == s2


# ---------------------------------------------------------------------------
# 6) Env OFF byte-identical to HEAD: the grip_polish call with the env var
#    UNSET must produce the same output as calling it without ever loading
#    the σ-only-mode helper.  We assert by checking the helper is a pure
#    superset operation: when ``_sigma_only_mode_active()`` returns False,
#    the polish never calls ``expand_hapto_for_sigma_only``.
# ---------------------------------------------------------------------------
class TestEnvOffByteIdentical:
    def test_polish_does_not_use_expansion_when_env_off(self, monkeypatch, rdkit_chem):
        """When the env flag is OFF, the polish must not consult
        ``expand_hapto_for_sigma_only`` (verified by monkeypatching the
        symbol to raise — if the polish called it, the test would crash).
        """
        Chem = rdkit_chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        P0 = np.array(
            [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE", raising=False)

        called = {"n": 0}

        def trip(*_args, **_kw):  # would crash if called
            called["n"] += 1
            raise RuntimeError("expand_hapto_for_sigma_only was called with env OFF")

        # Patch into the grip_polish module's namespace.
        import delfin.fffree.grip_polish as gp
        monkeypatch.setattr(gp, "expand_hapto_for_sigma_only", trip)

        # Run the polish; must NOT raise.
        from delfin.fffree.grip_polish import grip_polish
        out = grip_polish(P0, mol, metal=0, donors=[1], geom="")
        assert isinstance(out, np.ndarray)
        assert called["n"] == 0, "expand_hapto_for_sigma_only must NOT be called when env OFF"

    def test_polish_does_use_expansion_when_env_on(self, monkeypatch, rdkit_chem):
        """When the env flag is ON, the polish MUST consult the expansion."""
        Chem = rdkit_chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles("Cc1ccccc1")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        P0 = np.array(
            [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
            dtype=np.float64,
        )
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_SIGMA_ONLY_MODE", "1")

        called = {"n": 0}

        # Capture the real function and wrap it to count calls.
        import delfin.fffree.grip_polish as gp
        real_expand = gp.expand_hapto_for_sigma_only

        def counted_expand(*args, **kw):
            called["n"] += 1
            return real_expand(*args, **kw)

        monkeypatch.setattr(gp, "expand_hapto_for_sigma_only", counted_expand)

        from delfin.fffree.grip_polish import grip_polish
        _ = grip_polish(P0, mol, metal=0, donors=[1], geom="")
        assert called["n"] >= 1, "expand_hapto_for_sigma_only should be called when env ON"


# ---------------------------------------------------------------------------
# 7) Defensive: empty donors, no rdkit available, error in expansion -> no
#    crash, returns base set.
# ---------------------------------------------------------------------------
class TestDefensiveFallback:
    def test_empty_donors_returns_input(self, rdkit_chem):
        Chem = rdkit_chem
        from rdkit.Chem import RWMol
        rw = RWMol()
        fe = rw.AddAtom(Chem.Atom("Fe"))
        mol = rw.GetMol()
        base: set = set()
        out = expand_hapto_for_sigma_only(mol, fe, [], base)
        assert out == set()

    def test_none_mol_returns_input(self):
        base = {1, 2, 3}
        out = expand_hapto_for_sigma_only(None, 0, [], base)
        # Mutation safety: the input must not be modified.
        assert out == base
        assert out is not base  # returns a new set
