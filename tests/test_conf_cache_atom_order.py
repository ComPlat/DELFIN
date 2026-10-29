"""The ligand conformer memo hands back coordinates in the CALLER's atom order
(DELFIN_FFFREE_CONF_CACHE_ATOM_ORDER, 30.10.2026: RAGLUJ thiocyanate seated through S because the
canonical-SMILES key forgot that the second NCS fragment was cut out as (N, C, S))."""
import numpy as np
import pytest

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

from delfin.manta import assemble_ligand_confs as LC  # noqa: E402


def _frag(smi):
    m = Chem.MolFromSmiles(smi)
    assert m is not None
    return m


@pytest.fixture(autouse=True)
def _clear_cache():
    LC._CONF_CACHE.clear()
    yield
    LC._CONF_CACHE.clear()


def test_off_returns_cached_order_unchanged(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_CONF_CACHE_ATOM_ORDER", raising=False)
    a = LC._ligand_confs_from_mol(_frag("S=C=[N-]"), k=2)
    b = LC._ligand_confs_from_mol(_frag("[N-]=C=S"), k=2)
    assert b is a                                  # historic memo: same object, first order
    assert b[0][0] == "S"                          # atom 0 of the SECOND fragment is N, the memo says S


def test_on_remaps_to_the_callers_order(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_CONF_CACHE_ATOM_ORDER", "1")
    a = LC._ligand_confs_from_mol(_frag("S=C=[N-]"), k=2)
    b = LC._ligand_confs_from_mol(_frag("[N-]=C=S"), k=2)
    assert a[0][:3] == ["S", "C", "N"]
    assert b[0][:3] == ["N", "C", "S"]             # donor index 0 is N again
    assert len(a[1]) == len(b[1])
    for Pa, Pb in zip(a[1], b[1]):                 # the same conformers, permuted (no re-embed)
        assert np.allclose(Pa[[2, 1, 0]], Pb[:3])
    assert b[2].GetAtomWithIdx(0).GetSymbol() == "N"
    assert b[2].GetNumConformers() == a[2].GetNumConformers()
    cid = b[2].GetConformers()[0].GetId()
    assert np.allclose(np.asarray(b[2].GetConformer(cid).GetPositions()), b[1][0])


def test_on_identity_order_returns_the_cached_tuple(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_CONF_CACHE_ATOM_ORDER", "1")
    a = LC._ligand_confs_from_mol(_frag("CS(C)=O"), k=2)
    b = LC._ligand_confs_from_mol(_frag("CS(C)=O"), k=2)
    assert b is a


def test_on_heavier_ligand_keeps_every_bond(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_CONF_CACHE_ATOM_ORDER", "1")
    s1, s2 = "OC(=O)c1ccncc1", "c1cc(C(O)=O)ccn1"          # isonicotinic acid, two atom orders
    a = LC._ligand_confs_from_mol(_frag(s1), k=3)
    b = LC._ligand_confs_from_mol(_frag(s2), k=3)
    mq = Chem.AddHs(_frag(s2))
    assert b[0] == [x.GetSymbol() for x in mq.GetAtoms()]
    for P in b[1]:
        for bd in mq.GetBonds():
            assert np.linalg.norm(P[bd.GetBeginAtomIdx()] - P[bd.GetEndAtomIdx()]) < 1.9
    assert len(b[1]) == len(a[1])


def test_atom_order_map_rejects_a_different_graph():
    m1 = Chem.AddHs(_frag("CCO"))
    m2 = Chem.AddHs(_frag("COC"))
    assert LC._atom_order_map(m1, m2) is None
