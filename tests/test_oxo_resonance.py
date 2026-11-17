"""Degenerate terminal groups get their resonance length on the legacy path (DELFIN_FFFREE_OXO_RESONANCE)."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta._oxo_resonance import degenerate_terminal_sets, resonance_xyz  # noqa: E402


def _xyz(syms, P):
    return "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(syms, P)) + "\n"


def test_nitrate_equalised():
    m = Chem.MolFromSmiles("O=[N+]([O-])[O-]", sanitize=False)
    m.UpdatePropertyCache(strict=False)
    sets = degenerate_terminal_sets(m)
    assert len(sets) == 1
    c, term, target = sets[0]
    assert len(term) == 3 and 1.24 < target < 1.31
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    n = c
    P = np.zeros((4, 3))
    dirs = [np.array([1.0, 0, 0]), np.array([-0.5, 0.866, 0]), np.array([-0.5, -0.866, 0])]
    lens = [1.21, 1.45, 1.47]
    for k, o in enumerate(term):
        P[o] = P[n] + dirs[k] * lens[k]
    out = resonance_xyz(_xyz(syms, P), m, sets)
    Q = np.array([[float(v) for v in ln.split()[1:]] for ln in out.strip().splitlines()])
    for o in term:
        assert abs(np.linalg.norm(Q[o] - Q[n]) - target) < 1e-4


def test_carboxylic_acid_untouched():
    m = Chem.MolFromSmiles("CC(=O)O")
    m = Chem.AddHs(m)
    assert degenerate_terminal_sets(m) == []
