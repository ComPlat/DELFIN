"""A chelate spin of the legacy topology builder stays inside the donor lone-pair cone
(DELFIN_FFFREE_ORIENT_CHELATE_LP, 27.10.2026)."""
import math

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from delfin.manta.ligand_placement import _align_and_orient_ligands  # noqa: E402

FLAG = "DELFIN_FFFREE_ORIENT_CHELATE_LP"


def _phen_zn_nh3():
    """Zn(phen)(NH3): phen in the metal plane, an ammonia whose H sit in the phen ring plane
    so that the clash score wants to spin the phen."""
    phen = Chem.AddHs(Chem.MolFromSmiles("c1cnc2c(c1)ccc1cccnc12"))
    assert AllChem.EmbedMolecule(phen, randomSeed=5) == 0
    P = phen.GetConformer().GetPositions()
    ns = [a.GetIdx() for a in phen.GetAtoms() if a.GetSymbol() == "N"]
    lp = np.zeros(3)
    for n in ns:
        for nb in phen.GetAtomWithIdx(n).GetNeighbors():
            v = P[nb.GetIdx()] - P[n]
            lp -= v / np.linalg.norm(v)
    mid = P[ns].mean(axis=0)
    zn = mid + lp / np.linalg.norm(lp) * 1.4
    nh3 = Chem.AddHs(Chem.MolFromSmiles("N"))
    assert AllChem.EmbedMolecule(nh3, randomSeed=5) == 0
    Q = nh3.GetConformer().GetPositions()
    Q = Q - Q[0] + zn + lp / np.linalg.norm(lp) * 2.05
    mol = Chem.RWMol(Chem.CombineMols(phen, nh3))
    zi = mol.AddAtom(Chem.Atom(30))
    m = mol.GetMol()
    X = np.vstack([P, Q, zn[None, :]])
    donors = ns + [phen.GetNumAtoms()]
    return m, X, zi, donors


def _lp_angle(mol, X, metal, d):
    nbs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(d).GetNeighbors() if nb.GetAtomicNum() != 30]
    lp = -sum((X[j] - X[d]) / np.linalg.norm(X[j] - X[d]) for j in nbs)
    dm = X[metal] - X[d]
    return math.degrees(math.acos(np.dot(lp, dm) / np.linalg.norm(lp) / np.linalg.norm(dm)))


def _run(mol, X, metal, donors):
    coords = [tuple(p) for p in X]
    _align_and_orient_ligands(coords, mol, metal, donors)
    return np.array(coords, float)


def test_chelate_lp_off_is_byte_identical(monkeypatch):
    mol, X, metal, donors = _phen_zn_nh3()
    monkeypatch.delenv(FLAG, raising=False)
    ref = _run(mol, X, metal, donors)
    monkeypatch.setenv(FLAG, "0")
    assert np.array_equal(ref, _run(mol, X, metal, donors))


@pytest.mark.parametrize("carry", ["0", "1"])
def test_chelate_lp_keeps_lone_pairs_on_the_metal(monkeypatch, carry):
    mol, X, metal, donors = _phen_zn_nh3()
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", carry)
    monkeypatch.setenv(FLAG, "1")
    Y = _run(mol, X, metal, donors)
    for d in donors[:2]:
        assert np.allclose(Y[d], X[d])
        assert _lp_angle(mol, Y, metal, d) <= max(15.0, _lp_angle(mol, X, metal, d)) + 1e-6
