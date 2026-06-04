"""Tests for delfin.fffree.symmetry_equivalence (Iter 6)."""
from __future__ import annotations

import numpy as np
import pytest


def test_benzene_classes():
    """Benzene: 6 C-C ring bonds in one orbit, 6 C-H bonds in another."""
    from delfin.fffree.symmetry_equivalence import find_equivalence_classes
    syms = ["C"] * 6 + ["H"] * 6
    bonds = [(i, (i + 1) % 6) for i in range(6)] + [(i, i + 6) for i in range(6)]
    arom = [True] * 6 + [False] * 6
    classes = find_equivalence_classes(syms, bonds, aromatic_atoms=arom)
    assert len(classes) == 2
    sizes = sorted(len(c) for c in classes)
    assert sizes == [6, 6]


def test_pyridine_classes():
    """Pyridine: aromaticity breaks the 6-fold; we expect:
    - 2× C-N bonds (orbit of size 2)
    - 2× C-C bonds α to N (orbit of size 2)
    - 1× C-C bond γ to N (orbit of size 2 — the para bond)
    - Plus C-H bonds in 3 orbits.

    Total bond classes: 3 (heavy) + 3 (C-H, since N has no H here).
    """
    from delfin.fffree.symmetry_equivalence import find_equivalence_classes
    syms = ["N", "C", "C", "C", "C", "C"] + ["H"] * 5
    # 6-ring + H attached to each C (5 H total, atoms 6..10)
    bonds = [(i, (i + 1) % 6) for i in range(6)]
    for k in range(5):
        bonds.append((k + 1, k + 6))
    arom = [True] * 6 + [False] * 5
    classes = find_equivalence_classes(syms, bonds, aromatic_atoms=arom)
    # We just check the C-N class
    cn_class = None
    for cls in classes:
        if any(syms[a] != syms[b] and {syms[a], syms[b]} == {"C", "N"} for (a, b) in cls):
            cn_class = cls
            break
    assert cn_class is not None
    assert len(cn_class) == 2  # the two C-N edges of pyridine


def test_singleton_class_no_penalty():
    """A single bond with no automorphism partner contributes zero variance."""
    from delfin.fffree.symmetry_equivalence import variance_penalty
    classes = [[(0, 1)]]
    P = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    loss, grad = variance_penalty(P, classes)
    assert loss == 0.0
    assert np.allclose(grad, 0.0)


def test_two_equivalent_bonds_with_variance():
    """Two equivalent bonds, one stretched: variance loss > 0, gradient points
    toward equalisation."""
    from delfin.fffree.symmetry_equivalence import variance_penalty
    classes = [[(0, 1), (2, 3)]]
    P = np.array([
        [0.0, 0.0, 0.0],
        [1.40, 0.0, 0.0],   # bond 0-1: 1.40
        [10.0, 0.0, 0.0],
        [11.30, 0.0, 0.0],  # bond 2-3: 1.30
    ])
    loss, grad = variance_penalty(P, classes, alpha=10.0)
    # Variance = ((1.40 - 1.35)^2 + (1.30 - 1.35)^2) / 2 = 0.0025
    # loss = 10 * 2 * 0.0025 = 0.05
    assert abs(loss - 0.05) < 1e-6
    # Gradient on atom 0 should push it AWAY from atom 1 (lengthen 1-0 bond? no
    # — the 0-1 bond is the LONGER one. To reduce variance the longer bond
    # should shrink, so atom 0 should move toward atom 1.).
    # Actually d - mu = 1.40 - 1.35 = +0.05, coef positive, unit = (P0-P1)/d.
    # P0-P1 = -1.40 in x, so grad[0] in -x direction; grad descent step is -grad,
    # so atom 0 moves in +x — toward atom 1. ✓
    assert grad[0, 0] < 0  # gradient is in -x


def test_classes_lex_sorted_deterministic():
    """Calling twice should give bit-identical class lists."""
    from delfin.fffree.symmetry_equivalence import find_equivalence_classes
    syms = ["C"] * 6 + ["H"] * 6
    bonds = [(i, (i + 1) % 6) for i in range(6)] + [(i, i + 6) for i in range(6)]
    arom = [True] * 6 + [False] * 6
    c1 = find_equivalence_classes(syms, bonds, aromatic_atoms=arom)
    c2 = find_equivalence_classes(syms, bonds, aromatic_atoms=arom)
    assert c1 == c2


def test_is_enabled_default_off(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_SYMMETRY_EQUIVALENCE", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    from delfin.fffree.symmetry_equivalence import is_enabled
    assert is_enabled() is False


def test_is_enabled_via_pure_track3(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_SYMMETRY_EQUIVALENCE", raising=False)
    monkeypatch.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")
    from delfin.fffree.symmetry_equivalence import is_enabled
    assert is_enabled() is True


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
