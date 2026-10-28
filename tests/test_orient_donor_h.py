"""The hydrogens of a one-heavy-neighbour donor turn away from the metal (DELFIN_FFFREE_ORIENT_DONOR_H)."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from delfin.manta.ligand_placement import _align_and_orient_ligands  # noqa: E402


def _zn_methanol():
    m = Chem.AddHs(Chem.MolFromSmiles("CO"))
    assert AllChem.EmbedMolecule(m, randomSeed=3) == 0
    X = m.GetConformer().GetPositions()
    o = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O")
    h = next(nb.GetIdx() for nb in m.GetAtomWithIdx(o).GetNeighbors() if nb.GetAtomicNum() == 1)
    # metal 2.0 A from O along the O-H direction: the H points straight at the metal
    v = (X[h] - X[o]) / np.linalg.norm(X[h] - X[o])
    zn = X[o] + 2.0 * v
    rw = Chem.RWMol(m)
    zi = rw.AddAtom(Chem.Atom(30))
    mol = rw.GetMol()
    return mol, np.vstack([X, zn[None, :]]), zi, [o], h


def _run(mol, X, metal, donors):
    coords = [tuple(p) for p in X]
    _align_and_orient_ligands(coords, mol, metal, donors)
    return np.array(coords, float)


def test_donor_h_off_is_byte_identical(monkeypatch):
    mol, X, zi, donors, h = _zn_methanol()
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_DONOR_H", raising=False)
    ref = _run(mol, X, zi, donors)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_DONOR_H", "0")
    assert np.array_equal(ref, _run(mol, X, zi, donors))


def test_donor_h_turns_away(monkeypatch):
    mol, X, zi, donors, h = _zn_methanol()
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_DONOR_H", raising=False)
    A = _run(mol, X, zi, donors)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_DONOR_H", "1")
    B = _run(mol, X, zi, donors)
    o = donors[0]
    assert np.linalg.norm(B[h] - B[zi]) > np.linalg.norm(A[h] - A[zi])
    assert abs(np.linalg.norm(B[h] - B[o]) - np.linalg.norm(A[h] - A[o])) < 1e-6
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    assert np.allclose(A[heavy], B[heavy])
