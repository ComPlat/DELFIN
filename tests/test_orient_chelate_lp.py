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


@pytest.mark.parametrize("carry,mode", [("0", "1"), ("1", "1"), ("0", "2"), ("1", "2")])
def test_chelate_lp_keeps_lone_pairs_on_the_metal(monkeypatch, carry, mode):
    mol, X, metal, donors = _phen_zn_nh3()
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", carry)
    monkeypatch.setenv(FLAG, mode)
    Y = _run(mol, X, metal, donors)
    for d in donors[:2]:
        assert np.allclose(Y[d], X[d])
        assert _lp_angle(mol, Y, metal, d) <= max(15.0, _lp_angle(mol, X, metal, d)) + 1e-6


def test_lp_weight_off_is_byte_identical(monkeypatch):
    mol, X, metal, donors = _phen_zn_nh3()
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_LP_WEIGHT", raising=False)
    ref = _run(mol, X, metal, donors)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_LP_WEIGHT", "0")
    assert np.array_equal(ref, _run(mol, X, metal, donors))


def test_lp_weight_does_not_misalign_the_monodentate(monkeypatch):
    mol, X, metal, donors = _phen_zn_nh3()
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_LP_WEIGHT", raising=False)
    a0 = _lp_angle(mol, _run(mol, X, metal, donors), metal, donors[2])
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_LP_WEIGHT", "5")
    a1 = _lp_angle(mol, _run(mol, X, metal, donors), metal, donors[2])
    assert a1 <= max(a0, 15.0) + 1e-6


def _fold(mol, X, metal, donors):
    d1, d2 = X[donors[0]], X[donors[1]]
    u = (d2 - d1) / np.linalg.norm(d2 - d1)
    mid = 0.5 * (d1 + d2)
    xs = [X[nb.GetIdx()] for d in donors for nb in mol.GetAtomWithIdx(d).GetNeighbors()
          if nb.GetAtomicNum() > 1 and nb.GetAtomicNum() != 30]
    xv, mv = np.mean(xs, axis=0) - mid, X[metal] - mid
    xv, mv = xv - np.dot(xv, u) * u, mv - np.dot(mv, u) * u
    return 180.0 - math.degrees(math.acos(np.dot(xv, mv) / np.linalg.norm(xv) / np.linalg.norm(mv)))


def _spun(mol, X, metal, donors, deg):
    """The phen spun by deg about its N...N axis (folded chelate)."""
    d1, d2 = X[donors[0]], X[donors[1]]
    u = (d2 - d1) / np.linalg.norm(d2 - d1)
    mid = 0.5 * (d1 + d2)
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    Y = X.copy()
    n_phen = donors[2]
    for i in range(n_phen):
        v = X[i] - mid
        Y[i] = mid + v * c + np.cross(u, v) * s + u * np.dot(u, v) * (1 - c)
    return Y


@pytest.mark.parametrize("deg,mode", [(0.0, "1"), (50.0, "1"), (120.0, "1"), (50.0, "2"), (120.0, "2")])
def test_chelate_plane_unfolds_the_chelate(monkeypatch, deg, mode):
    mol, X, metal, donors = _phen_zn_nh3()
    X = _spun(mol, X, metal, donors, deg)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", mode)
    Y = _run(mol, X, metal, donors)
    assert _fold(mol, Y, metal, donors[:2]) <= 20.0 + 1e-6
    for d in donors:
        assert np.allclose(Y[d], X[d])


def test_chelate_plane_off_is_byte_identical(monkeypatch):
    mol, X, metal, donors = _phen_zn_nh3()
    X = _spun(mol, X, metal, donors, 50.0)
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", raising=False)
    ref = _run(mol, X, metal, donors)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", "0")
    assert np.array_equal(ref, _run(mol, X, metal, donors))


@pytest.mark.parametrize("deg", [50.0, 120.0])
def test_chelate_plane_table_window(monkeypatch, tmp_path, deg):
    tab = tmp_path / "fold.tsv"
    tab.write_text("# test table\n5\t1\t12.5\t100\n")
    mol, X, metal, donors = _phen_zn_nh3()
    X = _spun(mol, X, metal, donors, deg)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", "3")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE_TABLE", str(tab))
    Y = _run(mol, X, metal, donors)
    assert _fold(mol, Y, metal, donors[:2]) <= 12.5 + 1e-6


def test_chelate_plane_table_missing_class_is_unconstrained(monkeypatch, tmp_path):
    tab = tmp_path / "fold.tsv"
    tab.write_text("6\t1\t40.0\t100\n")
    mol, X, metal, donors = _phen_zn_nh3()
    X = _spun(mol, X, metal, donors, 50.0)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", raising=False)
    ref = _run(mol, X, metal, donors)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", "3")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE_TABLE", str(tab))
    assert np.array_equal(ref, _run(mol, X, metal, donors))


@pytest.mark.parametrize("deg", [50.0, 120.0])
def test_chelate_plane_lexicographic_window(monkeypatch, tmp_path, deg):
    tab = tmp_path / "fold.tsv"
    tab.write_text("5\t1\t12.5\t100\n")
    mol, X, metal, donors = _phen_zn_nh3()
    X = _spun(mol, X, metal, donors, deg)
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE", "4")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CHELATE_PLANE_TABLE", str(tab))
    Y = _run(mol, X, metal, donors)
    assert _fold(mol, Y, metal, donors[:2]) <= 12.5 + 1e-6


def _md(X, metal, donors):
    return [float(np.linalg.norm(X[d] - X[metal])) for d in donors]


def test_swing_off_is_byte_identical(monkeypatch, tmp_path):
    mol, X, metal, donors = _phen_zn_nh3()
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_SWING", raising=False)
    ref = _run(mol, X, metal, donors)
    tab = tmp_path / "sw.tsv"
    tab.write_text("3\t8.0\t10.0\t100\n")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_SWING", "0")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_SWING_TABLE", str(tab))
    assert np.array_equal(ref, _run(mol, X, metal, donors))


def test_swing_keeps_md_and_bound(monkeypatch, tmp_path):
    mol, X, metal, donors = _tbuo_like()
    tab = tmp_path / "sw.tsv"
    tab.write_text("2\t6.0\t8.0\t100\n")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_CARRY_H", "1")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_SCORE_H", "1")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_SWING", "1")
    monkeypatch.setenv("DELFIN_FFFREE_ORIENT_SWING_TABLE", str(tab))
    Y = _run(mol, X, metal, donors)
    assert np.allclose(_md(Y, metal, donors), _md(X, metal, donors), atol=1e-6)
    for d in donors:
        u, v = X[d] - X[metal], Y[d] - Y[metal]
        ang = math.degrees(math.acos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)))
        assert ang <= 6.0 + 1e-6


def _tbuo_like():
    from test_orient_score_h import _tbuo_pair
    return _tbuo_pair()
