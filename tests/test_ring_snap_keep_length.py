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


def test_snap_projection_restores_every_bond_not_only_the_mean(monkeypatch):
    """Unequal tilts (one atom lifted far more than the others): the uniform scale alone would
    leave the most-tilted bonds short; every bond comes back to its pre-projection length."""
    monkeypatch.setenv(FLAG, "1")
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    P = _coords(_puckered_benzene_xyz(z=0.04))
    P[0, 2] += 0.45
    P[6, 2] += 0.45                                  # its H rides along
    xyz = "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(["C"] * 6 + ["H"] * 6, P)) + "\n"
    P1 = _coords(_snap_aromatic_rings_in_xyz(xyz, mol))
    assert max(abs(a - b) for a, b in zip(_benzene_cc(P1), _benzene_cc(P))) < 2e-3


def test_sp2_flatten_keeps_the_three_bonds_better_with_the_switch(monkeypatch):
    from rdkit.Chem import AllChem
    from delfin.manta.topology_checks import _flatten_sp2_atoms_xyz
    mol = Chem.AddHs(Chem.MolFromSmiles("C=C(C)C"))
    assert AllChem.EmbedMolecule(mol, randomSeed=3) == 0
    X = mol.GetConformer().GetPositions()
    nb = [n.GetIdx() for n in mol.GetAtomWithIdx(1).GetNeighbors() if n.GetAtomicNum() > 1]
    nrm = np.cross(X[nb[1]] - X[nb[0]], X[nb[2]] - X[nb[0]]); nrm /= np.linalg.norm(nrm)
    X[1] = X[1] + 0.40 * nrm                          # pyramidalise the sp2 centre
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(syms, X)) + "\n"
    L0 = [float(np.linalg.norm(X[1] - X[j])) for j in nb]

    def _err(out):
        Y = _coords(out)
        return max(abs(float(np.linalg.norm(Y[1] - Y[j])) - t) for j, t in zip(nb, L0))

    monkeypatch.delenv(FLAG, raising=False)
    e_off = _err(_flatten_sp2_atoms_xyz(xyz, mol))
    monkeypatch.setenv(FLAG, "1")
    out_on = _flatten_sp2_atoms_xyz(xyz, mol)
    Y = _coords(out_on)
    assert abs(float(np.dot(Y[1] - Y[nb[0]], nrm))) < 1e-6      # planar
    assert _err(out_on) < e_off                       # never worse; a planar centre cannot keep all three
