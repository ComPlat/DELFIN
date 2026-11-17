"""The ligand spin of the legacy topology builder moves the fragment's hydrogens too
(DELFIN_FFFREE_ORIENT_CARRY_H)."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from delfin.manta.ligand_placement import _align_and_orient_ligands  # noqa: E402

FLAG = "DELFIN_FFFREE_ORIENT_CARRY_H"


def _ethoxide_at(donor_pos, chain_dir):
    """An embedded ethoxide whose O sits at donor_pos and whose O->C bond points along chain_dir."""
    m = Chem.AddHs(Chem.MolFromSmiles("CC[O-]"))
    assert AllChem.EmbedMolecule(m, randomSeed=7) == 0
    X = m.GetConformer().GetPositions()
    o = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O")
    c = next(nb.GetIdx() for nb in m.GetAtomWithIdx(o).GetNeighbors())
    X = X - X[o]
    a = X[c] / np.linalg.norm(X[c])
    b = np.asarray(chain_dir, float) / np.linalg.norm(chain_dir)
    v, s, cth = np.cross(a, b), np.linalg.norm(np.cross(a, b)), float(np.dot(a, b))
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.eye(3) + K + K @ K * ((1 - cth) / (s * s)) if s > 1e-9 else np.eye(3)
    X = X @ R.T + np.asarray(donor_pos, float)
    conf = m.GetConformer()
    for i, p in enumerate(X):
        conf.SetAtomPosition(i, p.tolist())
    return m


def _complex():
    """Zn with two ethoxides 90 deg apart whose ethyl chains point into each other (a clash the spin resolves)."""
    l1 = _ethoxide_at((2.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    l2 = _ethoxide_at((0.0, 2.0, 0.0), (1.0, 0.0, 0.0))
    zn = Chem.MolFromSmiles("[Zn+2]")
    zc = Chem.Conformer(1)
    zc.SetAtomPosition(0, (0.0, 0.0, 0.0))
    zn.AddConformer(zc)
    mol = Chem.CombineMols(Chem.CombineMols(l1, l2), zn)
    X = mol.GetConformer().GetPositions()
    metal = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Zn")
    donors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"]
    return mol, X, metal, donors


def _xh(mol, X):
    return {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): float(np.linalg.norm(
        X[b.GetBeginAtomIdx()] - X[b.GetEndAtomIdx()])) for b in mol.GetBonds()
        if 1 in (b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum())}


def _run(mol, X, metal, donors):
    coords = [tuple(p) for p in X]
    _align_and_orient_ligands(coords, mol, metal, donors)
    return np.array(coords, float)


def test_spin_without_the_switch_leaves_the_hydrogens_behind(monkeypatch):
    monkeypatch.delenv(FLAG, raising=False)
    mol, X, metal, donors = _complex()
    Y = _run(mol, X, metal, donors)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    hyd = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    assert np.abs(Y[heavy] - X[heavy]).max() > 0.3          # the spin happened
    assert np.abs(Y[hyd] - X[hyd]).max() == 0.0             # ... and no hydrogen moved
    assert max(_xh(mol, Y).values()) > 1.5                  # torn X-H bonds


def test_spin_with_the_switch_is_a_rigid_motion_and_chooses_the_same_angle(monkeypatch):
    mol, X, metal, donors = _complex()
    monkeypatch.delenv(FLAG, raising=False)
    Y0 = _run(mol, X, metal, donors)
    monkeypatch.setenv(FLAG, "1")
    Y1 = _run(mol, X, metal, donors)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    assert np.abs(Y1[heavy] - Y0[heavy]).max() < 1e-9       # same choice: heavy atoms identical
    before, after = _xh(mol, X), _xh(mol, Y1)
    assert max(abs(after[k] - before[k]) for k in before) < 1e-9    # every X-H bond kept
    D0 = np.linalg.norm(X[:, None] - X[None], axis=2)
    D1 = np.linalg.norm(Y1[:, None] - Y1[None], axis=2)
    for frag in Chem.GetMolFrags(mol):
        if metal in frag:
            continue
        idx = list(frag)
        assert np.abs(D1[np.ix_(idx, idx)] - D0[np.ix_(idx, idx)]).max() < 1e-9   # rigid body
    assert np.abs(Y1[donors] - X[donors]).max() < 1e-9      # donors stay on their vertex
