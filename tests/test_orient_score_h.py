"""The ligand-spin score of the legacy topology builder counts the hydrogens it carries
(DELFIN_FFFREE_ORIENT_SCORE_H, Ursache 5, 27.10.2026)."""
import numpy as np
import pytest

pytest.importorskip("rdkit.Chem")

from delfin.manta.interpenetration import contact_census  # noqa: E402
from delfin.manta.ml_tables import _METAL_SET  # noqa: E402
from test_orient_carry_h import _complex, _run, _xh  # noqa: E402

CARRY = "DELFIN_FFFREE_ORIENT_CARRY_H"
SCORE = "DELFIN_FFFREE_ORIENT_SCORE_H"


def _tbuo_pair():
    """Zn with two tert-butoxides 90 deg apart: the methyl H of one are only cleared by a spin that sees H."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)(C)[O][Zn][O]C(C)(C)C"))
    assert AllChem.EmbedMolecule(mol, randomSeed=11) == 0
    X = mol.GetConformer().GetPositions()
    zn = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Zn")
    donors = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"]
    X = X - X[zn]
    # put the two O at 90 deg, 1.95 A, and move each fragment rigidly with its O
    targets = [np.array([1.95, 0.0, 0.0]), np.array([0.0, 1.95, 0.0])]
    frag = {}
    for d in donors:
        seen, st = set(), [d]
        while st:
            x = st.pop()
            if x in seen or x == zn:
                continue
            seen.add(x)
            st += [nb.GetIdx() for nb in mol.GetAtomWithIdx(x).GetNeighbors()]
        frag[d] = sorted(seen)
    for d, t in zip(donors, targets):
        X[frag[d]] += t - X[d]
    return mol, X, zn, donors


def _min_inter(mol, X, donors):
    return contact_census(mol, X, _METAL_SET, donors)["inter"][0]


@pytest.mark.parametrize("build", [_complex, _tbuo_pair])
def test_score_h_off_is_byte_identical(monkeypatch, build):
    mol, X, metal, donors = build()
    monkeypatch.setenv(CARRY, "1")
    monkeypatch.delenv(SCORE, raising=False)
    ref = _run(mol, X, metal, donors)
    monkeypatch.setenv(SCORE, "0")
    assert np.array_equal(ref, _run(mol, X, metal, donors))


def test_score_h_needs_carry_h(monkeypatch):
    mol, X, metal, donors = _tbuo_pair()
    monkeypatch.delenv(CARRY, raising=False)
    monkeypatch.delenv(SCORE, raising=False)
    ref = _run(mol, X, metal, donors)
    monkeypatch.setenv(SCORE, "1")
    assert np.array_equal(ref, _run(mol, X, metal, donors))


@pytest.mark.parametrize("build", [_complex, _tbuo_pair])
def test_score_h_never_worse_contact_and_keeps_bonds(monkeypatch, build):
    mol, X, metal, donors = build()
    monkeypatch.setenv(CARRY, "1")
    monkeypatch.delenv(SCORE, raising=False)
    blind = _run(mol, X, metal, donors)
    monkeypatch.setenv(SCORE, "1")
    seeing = _run(mol, X, metal, donors)
    # rigid motion: every X-H bond and every donor position is kept
    for k, v in _xh(mol, X).items():
        assert abs(_xh(mol, seeing)[k] - v) < 1e-6
    for d in donors:
        assert np.allclose(seeing[d], X[d])
    assert _min_inter(mol, seeing, donors) >= _min_inter(mol, blind, donors) - 1e-9


def test_contact_census_sees_inter_fragment_h_overlap():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    one = Chem.AddHs(Chem.MolFromSmiles("C"))
    assert AllChem.EmbedMolecule(one, randomSeed=3) == 0
    P = one.GetConformer().GetPositions()
    mol = Chem.CombineMols(one, one)
    far = np.vstack([P, P + np.array([6.0, 0.0, 0.0])])
    near = np.vstack([P, P + np.array([2.2, 0.0, 0.0])])
    c_far = contact_census(mol, far, _METAL_SET)
    c_near = contact_census(mol, near, _METAL_SET)
    assert c_far["inter"][0] > 1.0
    assert c_near["inter"][0] < 0.6
    assert c_near["bond"][0] > 0.9
