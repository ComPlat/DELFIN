"""A ring projected onto its plane keeps its bond lengths (DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH)."""
import math

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta._arom_planarize import correct_xyz  # noqa: E402
from delfin.manta._pi_h_projector import _parse_xyz  # noqa: E402
from delfin.manta.ligand_placement import _snap_aromatic_rings_in_xyz  # noqa: E402

FLAG = "DELFIN_FFFREE_RING_SNAP_KEEP_LENGTH"


def _puckered_benzene_xyz(z=0.25, r=1.39):
    """Benzene, ring atoms alternating +-z (strongly puckered), each H radially at 1.08 A, in the
    atom order of AddHs(c1ccccc1) (6 C, then the 6 H of C0..C5)."""
    C = [np.array([r * math.cos(2 * math.pi * k / 6), r * math.sin(2 * math.pi * k / 6),
                   z * (1 if k % 2 == 0 else -1)]) for k in range(6)]
    H = []
    for c in C:
        u = np.array([c[0], c[1], 0.0]); u /= np.linalg.norm(u)
        H.append(c + 1.08 * u)
    return "\n".join(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for s, p in
                     [("C", c) for c in C] + [("H", h) for h in H]) + "\n"


def _benzene_cc(P):
    return [float(np.linalg.norm(P[k] - P[(k + 1) % 6])) for k in range(6)]


def _benzene_ch(P):
    return [float(np.linalg.norm(P[6 + k] - P[k])) for k in range(6)]


def _coords(xyz):
    return np.array([[float(v) for v in ln.split()[1:4]] for ln in xyz.strip().splitlines()])


def test_snap_projection_shrinks_the_ring_without_the_switch(monkeypatch):
    monkeypatch.delenv(FLAG, raising=False)
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    xyz = _puckered_benzene_xyz()
    P0, P1 = _coords(xyz), _coords(_snap_aromatic_rings_in_xyz(xyz, mol))
    assert np.mean(_benzene_cc(P1)) < np.mean(_benzene_cc(P0)) - 0.03


def test_snap_projection_keeps_ring_and_xh_lengths_with_the_switch(monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    xyz = _puckered_benzene_xyz()
    P0, P1 = _coords(xyz), _coords(_snap_aromatic_rings_in_xyz(xyz, mol))
    assert abs(np.mean(_benzene_cc(P1)) - np.mean(_benzene_cc(P0))) < 1e-6
    assert max(abs(d - 1.08) for d in _benzene_ch(P1)) < 0.03
    z = P1[:6, 2] - P1[:6, 2].mean()
    assert np.abs(z).max() < 1e-6                  # still flat


def test_arom_planarize_keeps_ring_length_with_the_switch(monkeypatch):
    xyz = f"12\nbenzene\n{_puckered_benzene_xyz(z=0.20)}"
    monkeypatch.delenv(FLAG, raising=False)
    _s, p_off, _ = _parse_xyz(correct_xyz(xyz))
    monkeypatch.setenv(FLAG, "1")
    _s, p_on, _ = _parse_xyz(correct_xyz(xyz))
    _s, p0, _ = _parse_xyz(xyz)
    assert np.mean(_benzene_cc(p_off)) < np.mean(_benzene_cc(p0)) - 0.02
    assert abs(np.mean(_benzene_cc(p_on)) - np.mean(_benzene_cc(p0))) < 1e-6
