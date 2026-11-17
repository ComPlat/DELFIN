"""An additive lone-pair-seated sibling (DELFIN_FFFREE_ORIENT_ONE_SUB_LP) enters only if it is no worse
than its own primary in any geometry class and better in the donor angle.

Measured 05.11.2026 (48 dev systems, fast builds): the siblings carried mpl 15.4 % / bnd 22.5 % / pi 6.8 %
of frames defective against 8.7 / 9.6 / 3.8 % of the old frames -- the stages after the seat relax a bent
ligand into new defects in some systems.  Per class the count of local units beyond the class p99 of
the realism table (frame_realism, the same units as the geometry audit) is compared with the primary's.

Hand-built frames, a synthetic class table; no crystal structure."""
import json

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("networkx")

from delfin.manta import frame_realism as FR  # noqa: E402
from delfin.manta.topo_isomers import _one_sub_lp_sibling_keeps  # noqa: E402

SMI = "Cl[Zn-2](Cl)(Cl)[N+]1=CC=CC=C1"


def test_keeps_a_sibling_that_only_improves_the_donor_angle():
    par = {"pi": 0, "mpl": 1, "dang": 2, "hyb": 0, "con": 0, "bnd": 3}
    sib = {"pi": 0, "mpl": 1, "dang": 0, "hyb": 0, "con": 0, "bnd": 3}
    assert _one_sub_lp_sibling_keeps(sib, par)


@pytest.mark.parametrize("cat", ["pi", "mpl", "hyb", "con", "bnd"])
def test_drops_a_sibling_worse_in_any_other_class(cat):
    par = {"pi": 0, "mpl": 0, "dang": 2, "hyb": 0, "con": 0, "bnd": 0}
    sib = dict(par, dang=0)
    sib[cat] = 1
    assert not _one_sub_lp_sibling_keeps(sib, par)


def test_drops_a_sibling_that_does_not_improve_the_donor_angle():
    par = {"pi": 0, "mpl": 0, "dang": 1, "hyb": 0, "con": 0, "bnd": 0}
    assert not _one_sub_lp_sibling_keeps(dict(par), par)


def test_drops_when_either_frame_cannot_be_measured():
    par = {"pi": 0, "mpl": 0, "dang": 1, "hyb": 0, "con": 0, "bnd": 0}
    assert not _one_sub_lp_sibling_keeps(None, par)
    assert not _one_sub_lp_sibling_keeps(dict(par, dang=0), None)


@pytest.fixture
def table(tmp_path, monkeypatch):
    T = {"fine": {}, "coarse": {},
         "all": {q: {"n": 1000, "med": 0.0, "lo1": 0.0, "hi1": 5.0} for q in FR.CAT}}
    T["all"]["bond_len"] = {"n": 1000, "med": 1.0, "lo1": 0.9, "hi1": 1.1}
    T["all"]["ml_len"] = {"n": 1000, "med": 2.05, "lo1": 1.9, "hi1": 2.2}
    T["all"]["m_contact"] = {"n": 1000, "med": 2.0, "lo1": 1.4, "hi1": 3.0}
    T["all"]["nb_clash"] = {"n": 1000, "med": 0.9, "lo1": 0.7, "hi1": 1.0}
    p = tmp_path / "t.json"
    p.write_text(json.dumps(T))
    monkeypatch.setenv("DELFIN_FRAME_REALISM_TABLE", str(p))
    FR._TABLES.clear()
    FR._CHEM_CACHE.clear()
    return p


def _zn_py(oop_deg):
    from test_frame_realism import _zn_py_frame
    return _zn_py_frame(oop_deg)


def test_category_counts_see_the_metal_out_of_the_donor_plane(table):
    a = FR.category_excess_counts(_zn_py(0), SMI)
    b = FR.category_excess_counts(_zn_py(60), SMI)
    assert a is not None and b is not None
    assert set(a) == {"pi", "mpl", "dang", "hyb", "con", "bnd"}
    assert b["mpl"] > a["mpl"]


def test_category_counts_need_a_smiles_and_a_table(monkeypatch, tmp_path):
    FR.clear_context()
    assert FR.category_excess_counts(_zn_py(0), None) is None
    monkeypatch.setenv("DELFIN_FRAME_REALISM_TABLE", str(tmp_path / "missing.json"))
    FR._TABLES.clear()
    assert FR.category_excess_counts(_zn_py(0), SMI) is None
