"""A one-substituent donor is seated along its own lone pair, not anti to its substituent
(DELFIN_FFFREE_ORIENT_ONE_SUB_LP).

A synthetic zinc complex, hand-built (no crystal structure): every ligand sits with its donor
on an octahedron vertex and its D->X bond pointing straight away from the metal (M-D-X 180 deg),
which is how the template / fragment seating hands the ligands to _align_and_orient_ligands.

  +x  acetone          O=C(Me)2    sp2 donor, period 2      -> M-O-C 120, metal in the C plane
  +y  thioacetone      S=C(Me)2    sp2 donor, period 3      -> M-S-C 109.47 (cone, six azimuths)
  -x  methanethiolate  [S-]Me      sp3 donor                -> M-S-C 109.47
  -y  trimethylphosphine oxide O=PMe3 (ylidic P=O)          -> stays linear
  +z  acetonitrile     N#CMe       sp substituent            -> stays linear
  -z  methoxide        [O-]Me      sp3 O without H           -> left to the old rule (linear)
"""
import math

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta._lone_pair_space import _perp  # noqa: E402
from delfin.manta.ligand_placement import _align_and_orient_ligands, one_sub_lp_seating  # noqa: E402

FLAG = "DELFIN_FFFREE_ORIENT_ONE_SUB_LP"
M_D = 2.10

LIGANDS = [
    ("acetone", "CC(C)=O", "O", (1, 0, 0)),
    ("thioacetone", "CC(C)=S", "S", (0, 1, 0)),
    ("thiolate", "C[S-]", "S", (-1, 0, 0)),
    ("phosphine_oxide", "CP(C)(C)=O", "O", (0, -1, 0)),
    ("nitrile", "CC#N", "N", (0, 0, 1)),
    ("methoxide", "C[O-]", "O", (0, 0, -1)),
]


def _terminal_donor(m, sym):
    return next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == sym and a.GetDegree() == 1)


def _ligand_coords(m, d, direction):
    """Donor at M_D along ``direction``, its single substituent X further out on the same line,
    X's own substituents spread around that line (trigonal for two, tetrahedral for three,
    straight on for an sp X)."""
    u = np.asarray(direction, float) / np.linalg.norm(direction)
    p = _perp(u)
    q = np.cross(u, p)
    x = m.GetAtomWithIdx(d).GetNeighbors()[0]
    order = m.GetBondBetweenAtoms(d, x.GetIdx()).GetBondTypeAsDouble()
    bl = {("O", 2.0): 1.22, ("S", 2.0): 1.65, ("N", 3.0): 1.15, ("S", 1.0): 1.82,
          ("O", 1.0): 1.41}[(m.GetAtomWithIdx(d).GetSymbol(), order)]
    if x.GetSymbol() == "P":
        bl = 1.50
    X = np.zeros((m.GetNumAtoms(), 3))
    X[d] = M_D * u
    X[x.GetIdx()] = X[d] + bl * u
    rest = [nb.GetIdx() for nb in x.GetNeighbors() if nb.GetIdx() != d]
    if len(rest) == 1 and order == 3.0:
        dirs = [u]
    elif len(rest) == 2:
        a = math.radians(120.0)
        dirs = [-math.cos(a) * u + s * math.sin(a) * p for s in (1.0, -1.0)]
    else:
        a = math.radians(109.47)
        dirs = [-math.cos(a) * u + math.sin(a) * (math.cos(t) * p + math.sin(t) * q)
                for t in np.linspace(0.0, 2 * math.pi, len(rest), endpoint=False)]
    for j, dv in zip(rest, dirs):
        X[j] = X[x.GetIdx()] + 1.50 * np.asarray(dv)
    return X


def _complex():
    mol = Chem.MolFromSmiles("[Zn+2]")
    coords = [np.zeros(3)]
    donors, xs, names = [], [], []
    for name, smi, dsym, direction in LIGANDS:
        lm = Chem.MolFromSmiles(smi)
        d = _terminal_donor(lm, dsym)
        x = lm.GetAtomWithIdx(d).GetNeighbors()[0].GetIdx()
        off = mol.GetNumAtoms()
        coords.extend(list(_ligand_coords(lm, d, direction)))
        mol = Chem.CombineMols(mol, lm)
        donors.append(off + d)
        xs.append(off + x)
        names.append(name)
    rw = Chem.RWMol(mol)
    for d in donors:
        rw.AddBond(0, d, Chem.BondType.SINGLE)     # dative M-D bonds, as in the constructor's graph
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol, np.array(coords, float), donors, xs, names


def _vertex_angle_deg(a, b, c):
    v1, v2 = a - b, c - b
    cs = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cs))))


def _planarity_torsion_deg(p0, p1, p2, p3):
    b1 = (p2 - p1) / np.linalg.norm(p2 - p1)
    v = (p0 - p1) - np.dot(p0 - p1, b1) * b1
    w = (p3 - p2) - np.dot(p3 - p2, b1) * b1
    return math.degrees(math.atan2(float(np.dot(np.cross(b1, v), w)), float(np.dot(v, w))))


def _run(mol, X, donors):
    coords = [tuple(p) for p in X]
    _align_and_orient_ligands(coords, mol, 0, donors)
    return np.array(coords, float)


def _run_bent(mol, X, donors):
    """The bent seating: _align_and_orient_ligands inside one_sub_lp_seating() (the builder's additive
    sibling pass, see topo_isomers)."""
    coords = [tuple(p) for p in X]
    with one_sub_lp_seating():
        _align_and_orient_ligands(coords, mol, 0, donors)
    return np.array(coords, float)


def _mdx(Y, donors, xs, names):
    return {n: _vertex_angle_deg(Y[0], Y[d], Y[x]) for n, d, x in zip(names, donors, xs)}


def test_off_every_one_substituent_donor_stays_linear(monkeypatch):
    """The defect: without the switch the seating keeps M-D-X at 180 deg for every donor."""
    monkeypatch.delenv(FLAG, raising=False)
    mol, X, donors, xs, names = _complex()
    Y = _run(mol, X, donors)
    for n, a in _mdx(Y, donors, xs, names).items():
        assert a > 179.0, (n, a)


def test_the_switch_alone_leaves_a_direct_seating_linear(monkeypatch):
    """The bent seat is an ADDITIVE sibling: with the switch on, the ordinary seating call is unchanged
    (the primary frame is the old one, byte for byte); only the sibling pass bends."""
    mol, X, donors, xs, names = _complex()
    monkeypatch.delenv(FLAG, raising=False)
    Y0 = _run(mol, X, donors)
    monkeypatch.setenv(FLAG, "1")
    Y1 = _run(mol, X, donors)
    assert Y0.tobytes() == Y1.tobytes()


def test_the_seating_context_is_left_clean(monkeypatch):
    mol, X, donors, xs, names = _complex()
    _run_bent(mol, X, donors)
    Y = _run(mol, X, donors)
    for n, a in _mdx(Y, donors, xs, names).items():
        assert a > 179.0, (n, a)


def test_off_is_byte_identical_to_switch_zero(monkeypatch):
    mol, X, donors, xs, names = _complex()
    monkeypatch.delenv(FLAG, raising=False)
    Y0 = _run(mol, X, donors)
    monkeypatch.setenv(FLAG, "0")
    Y1 = _run(mol, X, donors)
    assert Y0.tobytes() == Y1.tobytes()


def test_on_seats_each_donor_by_its_hybridisation(monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    mol, X, donors, xs, names = _complex()
    Y = _run_bent(mol, X, donors)
    a = _mdx(Y, donors, xs, names)
    assert abs(a["acetone"] - 120.0) < 1.0, a
    assert abs(a["thioacetone"] - 109.47) < 1.0, a
    assert abs(a["thiolate"] - 109.47) < 1.0, a
    assert a["phosphine_oxide"] > 179.0, a     # ylidic P=O: the wide seed is kept
    assert a["nitrile"] > 179.0, a             # sp substituent: a linear donor
    assert a["methoxide"] > 179.0, a           # sp3 O without H: outside the scope
    assert np.abs(Y[donors] - X[donors]).max() < 1e-9     # donors stay on their vertices
    assert np.abs(Y[0] - X[0]).max() == 0.0


def test_on_sp2_donor_keeps_the_metal_in_the_plane_of_the_substituents(monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    mol, X, donors, xs, names = _complex()
    Y = _run_bent(mol, X, donors)
    for name in ("acetone",):
        k = names.index(name)
        d, x = donors[k], xs[k]
        r = next(nb.GetIdx() for nb in mol.GetAtomWithIdx(x).GetNeighbors() if nb.GetIdx() != d)
        t = abs(_planarity_torsion_deg(Y[0], Y[d], Y[x], Y[r]))
        assert min(t, 180.0 - t) < 1.0, (name, t)


def test_on_moves_each_ligand_as_a_rigid_body(monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    mol, X, donors, xs, names = _complex()
    Y = _run_bent(mol, X, donors)
    D0 = np.linalg.norm(X[:, None] - X[None], axis=2)
    D1 = np.linalg.norm(Y[:, None] - Y[None], axis=2)
    start = 1
    for _name, smi, _d, _dir in LIGANDS:
        idx = list(range(start, start + Chem.MolFromSmiles(smi).GetNumAtoms()))
        assert np.abs(D1[np.ix_(idx, idx)] - D0[np.ix_(idx, idx)]).max() < 1e-9
        start += len(idx)


def test_on_the_hydrogens_ride_with_the_bent_ligand(monkeypatch):
    """With explicit H the bend is a rigid motion of the whole ligand: no X-H bond is torn
    (left behind, the H read as torn bonds and the topology gate dropped the frame: WACNAQ)."""
    monkeypatch.setenv(FLAG, "1")
    monkeypatch.delenv("DELFIN_FFFREE_ORIENT_CARRY_H", raising=False)
    mol, X, donors, xs, names = _complex()
    k = names.index("acetone")
    d, x = donors[k], xs[k]
    rw = Chem.RWMol(mol)
    Xl = [tuple(p) for p in X]
    for nb in mol.GetAtomWithIdx(x).GetNeighbors():
        if nb.GetIdx() == d:
            continue
        c = nb.GetIdx()
        out = X[c] - X[x]
        out = out / np.linalg.norm(out)
        h = rw.AddAtom(Chem.Atom(1))
        rw.AddBond(c, h, Chem.BondType.SINGLE)
        Xl.append(tuple(X[c] + 1.09 * out))
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    X = np.array(Xl, float)
    Y = _run_bent(mol, X, donors)
    assert abs(_vertex_angle_deg(Y[0], Y[d], Y[x]) - 120.0) < 1.0
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if 1 in (mol.GetAtomWithIdx(i).GetAtomicNum(), mol.GetAtomWithIdx(j).GetAtomicNum()):
            assert abs(np.linalg.norm(Y[i] - Y[j]) - 1.09) < 1e-9


def test_on_is_deterministic(monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    mol, X, donors, xs, names = _complex()
    assert _run_bent(mol, X, donors).tobytes() == _run_bent(mol, X, donors).tobytes()
