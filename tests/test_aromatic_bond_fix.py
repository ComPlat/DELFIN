"""Unit tests for the overnight aromatic bond fix.

Covers:
* aromatic_bond_targets.aromatic_ideal / aromatic_mu basic semantics
* aromatic_bond_enforcement byte-identical no-op when flag OFF
* aromatic_bond_enforcement reduces C-C variance on the LUHMOT bug structure
* _bond_decollapse._ideal_bond aromatic-aware path (flag ON) vs legacy (flag OFF)

These tests do NOT require RDKit, Mogul, or CCDC; they run on plain
numpy/python so they can be wired into every CI iteration.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest


def _make_stretched_benzene():
    """Plain benzene with one bond pathologically compressed (1.13 Å)."""
    syms = ["C"] * 6 + ["H"] * 6
    R = 1.40
    th = np.linspace(0.0, 2 * np.pi, 7)[:-1]
    P = np.array(
        [[R * np.cos(t), R * np.sin(t), 0.0] for t in th]
        + [[(R + 1.08) * np.cos(t), (R + 1.08) * np.sin(t), 0.0] for t in th]
    )
    # Smash C0 onto C1: 1.13 Å between them
    move = P[1] - P[0]
    move = move * ((np.linalg.norm(move) - (np.linalg.norm(move) - 1.13)) /
                   np.linalg.norm(move))
    # Move atom 0 toward atom 1 so |P1-P0| = 1.13
    direction = (P[1] - P[0]) / np.linalg.norm(P[1] - P[0])
    new_p0 = P[1] - 1.13 * direction
    P[0] = new_p0
    return syms, P


def test_aromatic_targets_basic():
    from delfin.fffree.aromatic_bond_targets import aromatic_ideal, aromatic_mu

    hit = aromatic_ideal("C", "C")
    assert hit is not None
    mu, sigma, n = hit
    assert 1.38 < mu < 1.41
    assert 0.01 < sigma < 0.10
    assert n > 1_000_000

    # Unordered
    assert aromatic_ideal("N", "C") == aromatic_ideal("C", "N")
    # Missing pair
    assert aromatic_ideal("Cl", "Cl") is None
    # mu helper
    assert abs(aromatic_mu("C", "C") - 1.3992) < 0.001
    assert aromatic_mu("Cl", "Cl") is None
    assert aromatic_mu("Cl", "Cl", fallback=1.7) == 1.7


def test_aromatic_enforcement_byte_identical_when_off(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_AROMATIC_BONDS", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)

    from delfin.fffree.aromatic_bond_enforcement import enforce_aromatic_bonds
    syms, P = _make_stretched_benzene()
    syms_out, P_out = enforce_aromatic_bonds(syms, P, mol=None)
    assert syms_out == list(syms)
    assert np.array_equal(P_out, P)


def test_aromatic_enforcement_reduces_variance_when_on(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_AROMATIC_BONDS", "1")

    from delfin.fffree.aromatic_bond_enforcement import enforce_aromatic_bonds
    syms, P = _make_stretched_benzene()
    # In-ring bonds
    def cc_lens(Q):
        return [float(np.linalg.norm(Q[i] - Q[(i + 1) % 6])) for i in range(6)]
    before = cc_lens(P)
    assert min(before) < 1.20, f"smoke setup: expected min C-C < 1.20, got {min(before):.3f}"
    syms_out, P_out = enforce_aromatic_bonds(syms, P, mol=None)
    after = cc_lens(P_out)
    assert min(after) > min(before), (
        f"Expected min C-C to increase from {min(before):.3f} to >1.13, "
        f"got {min(after):.3f}"
    )
    assert np.std(after) < np.std(before), (
        f"Expected std-C-C drop: before={np.std(before):.4f}  "
        f"after={np.std(after):.4f}"
    )


def test_ideal_bond_aromatic_aware(monkeypatch):
    """_bond_decollapse._ideal_bond returns covalent sum by default and
    CCDC aromatic mean when both gate AND aromatic=True."""
    from delfin._bond_decollapse import _ideal_bond

    # Default (no flag): legacy covalent sum 1.52 regardless of aromatic flag
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_AROMATIC_AWARE", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    assert abs(_ideal_bond("C", "C") - 1.52) < 1e-9
    assert abs(_ideal_bond("C", "C", aromatic=True) - 1.52) < 1e-9, (
        "gate OFF: aromatic=True should still return covalent sum"
    )

    # Flag ON: aromatic=True → 1.40, aromatic=False → 1.52
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_AROMATIC_AWARE", "1")
    assert abs(_ideal_bond("C", "C") - 1.52) < 1e-9
    assert abs(_ideal_bond("C", "C", aromatic=True) - 1.3992) < 0.001
    assert abs(_ideal_bond("C", "N", aromatic=True) - 1.3726) < 0.001


def test_ideal_bond_unknown_pair_fallback(monkeypatch):
    """Unknown aromatic pair → fallback to covalent sum."""
    from delfin._bond_decollapse import _ideal_bond
    monkeypatch.setenv("DELFIN_FFFREE_GRIP_AROMATIC_AWARE", "1")
    # As-As has too few entries (n=18, below MIN_N=50) → fallback expected
    # But Cl-Cl is not in table at all → also fallback.
    res = _ideal_bond("Cl", "Cl", aromatic=True)
    # Covalent: 1.02 + 1.02 = 2.04
    assert abs(res - 2.04) < 1e-9


def test_ideal_bond_pure_track3_activates(monkeypatch):
    """PURE_TRACK3 should activate aromatic awareness via env."""
    from delfin._bond_decollapse import _ideal_bond
    monkeypatch.delenv("DELFIN_FFFREE_GRIP_AROMATIC_AWARE", raising=False)
    monkeypatch.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")
    assert abs(_ideal_bond("C", "C", aromatic=True) - 1.3992) < 0.001


def test_luhmot_bug_fixed(monkeypatch):
    """Smoke regression: load the actual LUHMOT artifact and verify the
    compressed C-C bond gets lifted."""
    monkeypatch.setenv("DELFIN_FFFREE_AROMATIC_BONDS", "1")

    p = Path("/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
             "e9f69af-iso-topo-heal-smoke500/080-LUHMOT.xyz")
    if not p.exists():
        pytest.skip("LUHMOT artifact not present in this environment")

    lines = p.read_text().splitlines()
    n = int(lines[0])
    syms, P = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        syms.append(parts[0])
        P.append([float(x) for x in parts[1:4]])
    P = np.array(P, dtype=float)

    from delfin.fffree.aromatic_bond_enforcement import enforce_aromatic_bonds
    syms_out, P_out = enforce_aromatic_bonds(syms, P, mol=None)

    # Min C-C bond should rise from 1.13 to at least 1.20
    def cc_bonds(Q):
        ds = []
        for i in range(len(syms)):
            if syms[i] != "C":
                continue
            for j in range(i + 1, len(syms)):
                if syms[j] != "C":
                    continue
                d = float(np.linalg.norm(Q[i] - Q[j]))
                if d < 1.65:
                    ds.append(d)
        return ds

    before = cc_bonds(P)
    after = cc_bonds(P_out)
    assert min(before) < 1.15, f"smoke setup: expected min < 1.15, got {min(before):.3f}"
    assert min(after) > min(before), (
        f"LUHMOT min C-C should rise: before={min(before):.3f}, "
        f"after={min(after):.3f}"
    )
    assert np.std(after) < np.std(before), (
        f"LUHMOT C-C std should drop: before={np.std(before):.4f}, "
        f"after={np.std(after):.4f}"
    )


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
