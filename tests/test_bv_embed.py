"""Metallacycle embed at bond-valence radii (DELFIN_FFFREE_ML_TARGETS_BV_EMBED)."""
import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta import _ml_targets_bv as BV  # noqa: E402
from delfin.manta.assemble_ligand_embed import _embed_metallacycle  # noqa: E402

# quinolin-8-olate as decompose() hands it over: neutral, H-free O (the quinolinone form)
QUINOLINATE = "O=C1C=CC=C2C=CC=NC12"


def _worst_bond_ratio(mol, syms, X):
    pt = Chem.GetPeriodicTable()
    return min(float(np.linalg.norm(X[b.GetBeginAtomIdx()] - X[b.GetEndAtomIdx()]))
               / (pt.GetRcovalent(syms[b.GetBeginAtomIdx()]) + pt.GetRcovalent(syms[b.GetEndAtomIdx()]))
               for b in mol.GetBonds())


def test_quinolinate_ring_embeds_at_the_target_radii(monkeypatch):
    lig = Chem.MolFromSmiles(QUINOLINATE)
    n = next(a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() == "N")
    dons = [n, 0]
    u = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])]
    tp = [u[0] * 2.14, u[1] * 2.06]
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setattr(BV, "_CTX", {"metal": "Zn", "ox": 2, "cn": 6, "by_class": {}, "by_element": {},
                                     "levels": [], "spins": [], "spin": None})
    monkeypatch.setenv(BV.EMBED_ENV, "1")
    res = _embed_metallacycle(lig, dons, "Zn", donor_target_pos=tp)
    assert res is not None
    syms, confs = res
    for X in confs:
        assert abs(np.linalg.norm(X[n]) - 2.14) < 0.08
        assert abs(np.linalg.norm(X[0]) - 2.06) < 0.08
        assert _worst_bond_ratio(lig, syms, X) > 0.85
    monkeypatch.delenv(BV.EMBED_ENV)
    syms0, confs0 = _embed_metallacycle(lig, dons, "Zn", donor_target_pos=tp)
    assert len(confs0) >= 1          # switch off: the historic embed, whatever it returns


def test_cn5_partner_geometry():
    from delfin.manta.converter_backend import _cn5_opposite_geometry
    assert _cn5_opposite_geometry("TBP-5 trigonal bipyramid") == "SPY-5 square pyramid"
    assert _cn5_opposite_geometry("SPY-5 square pyramid") == "TBP-5 trigonal bipyramid"
    assert _cn5_opposite_geometry("OC-6 octahedron") is None
