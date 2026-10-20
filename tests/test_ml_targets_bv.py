"""Bond-valence M-D targets (DELFIN_FFFREE_ML_TARGETS_BV): off is byte-identical, on reads the table."""
import math

import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta import _ml_targets_bv as BV  # noqa: E402
from delfin.manta import ml_tables as MT  # noqa: E402
from delfin.manta import polyhedra as PH  # noqa: E402
from delfin.manta._oxidation_state import _ox_cut, _ox_cut_full  # noqa: E402

CISPLATIN = "[Cl][Pt-2]([Cl])([NH3+])[NH3+]"
CU_ACAC = "CC1=CC(C)=[O+][Cu-2]2([O]1)[O+]=C(C)C=C(C)[O]2"


@pytest.fixture(autouse=True)
def _bv_env_reset(monkeypatch):
    for k in (BV.FLAG, BV.TABLE_ENV, BV.MIN_N_ENV, BV.LITERATURE_ENV):
        monkeypatch.delenv(k, raising=False)
    BV.clear_context()
    BV._TABLE_CACHE.clear()
    yield
    BV.clear_context()


def _write_bv_table(tmp_path, rows):
    p = tmp_path / "t.tsv"
    lines = ["# test table", "# kind\tmetal\tox\tcn\telement\tclass\tn_struct\tvalue"]
    lines += ["\t".join(str(x) for x in r) for r in rows]
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_cut_full_agrees_with_cut():
    for smi in (CISPLATIN, CU_ACAC, "O=[Re]([Cl])([Cl])[Cl]", "[Cl][Fe-3]([Cl])[Cl]"):
        m = Chem.MolFromSmiles(smi)
        mi = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() in ("Pt", "Cu", "Re", "Fe"))
        full = _ox_cut_full(m, mi)
        assert full is not None
        assert (full[0], full[1]) == _ox_cut(m, mi)


def test_donor_classes_of_known_ligands():
    m = Chem.MolFromSmiles(CISPLATIN)
    metal, ox, donors = BV.complex_donors(m, 1)  # Pt is atom 1
    assert (metal, ox) == ("Pt", 2)
    assert sorted((el, c) for _i, el, c in donors) == [("Cl", "t0x"), ("Cl", "t0x"),
                                                       ("N", "t3l"), ("N", "t3l")]
    m = Chem.MolFromSmiles(CU_ACAC)
    mi = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "Cu")
    metal, ox, donors = BV.complex_donors(m, mi)
    assert (metal, ox) == ("Cu", 2)
    # both acac O are the same class, although the SMILES writes them differently
    assert {c for _i, _el, c in donors} == {"pix"}
    py = Chem.MolFromSmiles("c1ccncc1")
    assert BV.donor_class(py.GetAtomWithIdx(3)) == "pil"
    aqua = Chem.AddHs(Chem.MolFromSmiles("O"))
    assert BV.donor_class(aqua.GetAtomWithIdx(0)) == "t2l"


def test_literature_r0_is_the_okeeffe_brese_formula():
    r0 = BV.r0_literature("Fe", "O")
    ri, ci, rj, cj = 1.16, 1.64, 0.63, 3.15
    want = ri + rj - ri * rj * (math.sqrt(ci) - math.sqrt(cj)) ** 2 / (ci * ri + cj * rj)
    assert r0 == pytest.approx(want)
    assert BV.r0_literature("Fe", "Xx") is None


def test_switch_off_changes_nothing(tmp_path, monkeypatch):
    before = (PH.md_distance("Pt", "N"), PH.md_distance("Pt", "Cl"),
              MT._get_ml_bond_length("Pt", "N"), MT._get_ml_bond_length("Pt", "Cl"))
    monkeypatch.setenv(BV.TABLE_ENV, _write_bv_table(tmp_path, [("D", "Pt", 2, 4, "N", "t3l", 99, 2.051)]))
    BV.set_context(CISPLATIN)
    assert BV.current_context() is None
    after = (PH.md_distance("Pt", "N"), PH.md_distance("Pt", "Cl"),
             MT._get_ml_bond_length("Pt", "N"), MT._get_ml_bond_length("Pt", "Cl"))
    assert before == after


def test_switch_on_reads_the_table_levels(tmp_path, monkeypatch):
    path = _write_bv_table(tmp_path, [
        ("D", "Pt", 2, 4, "N", "t3l", 99, 2.051),
        ("R", "Pt", 2, "*", "Cl", "*", 99, 2.05),       # Cl -> R0 - b ln(2/4)
        ("D", "Pt", 2, 4, "Cl", "t0x", 5, 9.9),         # too few structures: ignored
    ])
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, path)
    BV.set_context(CISPLATIN)
    ctx = BV.current_context()
    assert ctx["metal"] == "Pt" and ctx["ox"] == 2 and ctx["cn"] == 4
    cl = 2.05 - BV.B_BV * math.log(2 / 4)
    assert PH.md_distance("Pt", "N") == pytest.approx(2.051)
    assert PH.md_distance("Pt", "Cl") == pytest.approx(cl)
    assert MT._get_ml_bond_length("Pt", "N") == pytest.approx(2.051)
    assert MT._get_ml_bond_length("Pt", "Cl") == pytest.approx(cl)
    # another metal or element: no answer, historic value
    assert BV.target_for_donor("Pd", "N") is None
    assert BV.target_for_donor("Pt", "P") is None
    BV.clear_context()
    monkeypatch.delenv(BV.FLAG)
    assert MT._get_ml_bond_length("Pt", "N") != pytest.approx(2.051)


def test_class_lookup_from_a_ligand_fragment_atom(tmp_path, monkeypatch):
    path = _write_bv_table(tmp_path, [("D", "Pt", 2, 4, "N", "t3l", 99, 2.051),
                                      ("D", "Pt", 2, 4, "Cl", "t0x", 99, 2.31)])
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, path)
    BV.set_context(CISPLATIN)
    nh3 = Chem.AddHs(Chem.MolFromSmiles("N"))
    assert PH.md_distance("Pt", "N", atom=nh3.GetAtomWithIdx(0), mol=nh3) == pytest.approx(2.051)


def test_unlike_donors_of_one_element_get_no_element_target(tmp_path, monkeypatch):
    # Re(V) oxo + ether O: the per-element question must not average 1.68 and 2.20
    smi = "O=[Re]([Cl])([Cl])([Cl])[O+](C)C"
    m = Chem.MolFromSmiles(smi)
    mi = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "Re")
    _metal, ox, donors = BV.complex_donors(m, mi)
    cls_o = sorted(c for _i, el, c in donors if el == "O")
    rows = [("D", "Re", ox, len(donors), "O", c, 99, v) for c, v in zip(cls_o, (1.68, 2.20))]
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, _write_bv_table(tmp_path, rows))
    BV.set_context(smi)
    ctx = BV.current_context()
    assert ctx is not None and "O" not in ctx["by_element"]
    assert BV.target_for_donor("Re", "O") is None


def test_no_context_for_what_the_cut_declines(monkeypatch):
    monkeypatch.setenv(BV.FLAG, "1")
    BV.set_context("Cl[Cu]Cl.Cl[Cu]Cl")        # two metals
    assert BV.current_context() is None
    BV.set_context("CCO")                      # no metal
    assert BV.current_context() is None
