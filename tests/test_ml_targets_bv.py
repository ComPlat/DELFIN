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
    for k in (BV.FLAG, BV.TABLE_ENV, BV.MIN_N_ENV, BV.LITERATURE_ENV, BV.SCOPE_ENV, BV.SPIN_BOTH_ENV):
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


def test_the_shipped_cod_table_is_the_default_and_reads():
    assert BV.table_path() == BV.SHIPPED_TABLE
    tbl = BV.load_table(BV.SHIPPED_TABLE)
    assert len(tbl) > 100
    # every row is a COD aggregate with its sample size
    assert all(n >= 10 for k, (n, _v) in tbl.items() if k[0] != "P")
    assert 0.05 <= BV.bv_params(tbl)[0] <= 0.45
    t, lev = BV.donor_target(tbl, "Co", 3, 6, "N", "sp3l", 0)
    assert lev is not None and 1.85 < t < 2.05


def test_seat_scope_leaves_every_element_question_alone(tmp_path, monkeypatch):
    path = _write_bv_table(tmp_path, [("D", "Pt", 2, 4, "N", "t3l", 99, 2.051)])
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, path)
    monkeypatch.setenv(BV.SCOPE_ENV, "seat")
    before = MT._get_ml_bond_length("Pt", "N")
    BV.set_context(CISPLATIN)
    nh3 = Chem.AddHs(Chem.MolFromSmiles("N"))
    assert PH.md_distance("Pt", "N", atom=nh3.GetAtomWithIdx(0), mol=nh3) == pytest.approx(2.051)
    assert BV.target_for_donor("Pt", "N") is None
    assert MT._get_ml_bond_length("Pt", "N") == before


FE_PY6 = "[N+]1(=CC=CC=C1)[Fe-4]([N+]1=CC=CC=C1)([N+]1=CC=CC=C1)([N+]1=CC=CC=C1)([N+]1=CC=CC=C1)[N+]1=CC=CC=C1"


def test_spin_states_come_from_the_measured_low_spin_share(tmp_path, monkeypatch):
    rows = [("F", "Fe", 2, "6/6", "*", "*", 99, 0.5),
            ("D", "Fe", 2, 6, "N", "pil", 99, 2.10),
            ("D", "Fe", 2, 6, "N", "pil:ls", 99, 1.97),
            ("D", "Fe", 2, 6, "N", "pil:hs", 99, 2.18)]
    p = tmp_path / "t.tsv"
    p.write_text("\n".join("\t".join(str(x) for x in r) for r in rows) + "\n")
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, str(p))
    BV._TABLE_CACHE.clear()
    BV.set_context(FE_PY6)                                  # spin switch off: spin-agnostic rows
    ctx = BV.current_context()
    assert ctx is not None, "Fe(py)6 must give a context"
    assert ctx["spins"] == [] and BV.target_for_donor("Fe", "N") == pytest.approx(2.10)
    monkeypatch.setenv(BV.SPIN_BOTH_ENV, "1")
    BV.set_context(FE_PY6)
    assert BV.current_context()["spins"] == ["ls", "hs"]
    BV.set_context(FE_PY6, "ls")
    assert BV.target_for_donor("Fe", "N") == pytest.approx(1.97)
    BV.set_context(FE_PY6, "hs")
    assert BV.target_for_donor("Fe", "N") == pytest.approx(2.18)
    BV.clear_context()


def test_the_ligand_field_rule_for_unimodal_cases():
    assert BV.spin_by_rule("Co", 3, [("O", "t2l")] * 6) == "ls"
    assert BV.spin_by_rule("Fe", 2, [("O", "t2l")] * 6) == "hs"
    assert BV.spin_by_rule("Fe", 2, [("N", "pil")] * 6) == "ls"
    assert BV.spin_by_rule("Co", 2, [("N", "pil")] * 6) == "hs"
    assert BV.spin_candidate("Fe", 2, 6) and not BV.spin_candidate("Fe", 2, 4)
    assert not BV.spin_candidate("Ni", 2, 6) and not BV.spin_candidate("Ru", 2, 6)


def test_pins_scope_answers_pins_but_not_checks(tmp_path, monkeypatch):
    path = _write_bv_table(tmp_path, [("D", "Pt", 2, 4, "N", "t3l", 99, 2.051)])
    before = MT._get_ml_bond_length("Pt", "N")
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, path)
    monkeypatch.setenv(BV.SCOPE_ENV, "pins")
    BV.set_context(CISPLATIN)
    assert MT._get_ml_bond_length("Pt", "N") == before             # a check: historic value
    assert MT._pin_ml_bond_length("Pt", "N") == pytest.approx(2.051)  # a pin: the target
    nh3 = Chem.AddHs(Chem.MolFromSmiles("N"))
    assert PH.md_distance("Pt", "N", atom=nh3.GetAtomWithIdx(0), mol=nh3) == pytest.approx(2.051)
    monkeypatch.delenv(BV.FLAG)
    assert MT._pin_ml_bond_length("Pt", "N") == before


def test_fragment_donor_that_lost_its_charge_reads_its_environment(monkeypatch):
    # decompose() hands a quinolinate O over neutral and H-free ('pil'); the complex knows it as
    # 'pix' next to a methanol O ('sp3l', 0.13 A longer, so no per-element target)
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setattr(BV, "_CTX", {"metal": "Zn", "ox": 2, "cn": 6, "spins": [], "spin": None,
                                     "by_class": {("O", "pix"): 2.058, ("O", "sp3l"): 2.192,
                                                  ("N", "pix"): 2.139},
                                     "by_element": {"N": 2.139}, "levels": ["N"]})
    frag = Chem.AddHs(Chem.MolFromSmiles("O=C1C=CC=C2C=CC=NC12"))
    o = frag.GetAtomWithIdx(0)
    assert BV.donor_class(o) == "pil"
    assert BV.target_for_donor("Zn", "O", atom=o) == pytest.approx(2.058)
    meoh = Chem.AddHs(Chem.MolFromSmiles("CO"))
    assert BV.target_for_donor("Zn", "O", atom=meoh.GetAtomWithIdx(1)) == pytest.approx(2.192)
    assert BV.target_for_donor("Zn", "O") is None


def test_lengths_scope_answers_every_length_but_not_checks(tmp_path, monkeypatch):
    path = _write_bv_table(tmp_path, [("D", "Pt", 2, 4, "N", "t3l", 99, 2.051)])
    before = MT._get_ml_bond_length("Pt", "N")
    before_md = PH.md_distance("Pt", "N")
    monkeypatch.setenv(BV.FLAG, "1")
    monkeypatch.setenv(BV.TABLE_ENV, path)
    monkeypatch.setenv(BV.SCOPE_ENV, "lengths")
    BV.set_context(CISPLATIN)
    assert MT._get_ml_bond_length("Pt", "N") == pytest.approx(2.051)   # a builder's length
    assert MT._pin_ml_bond_length("Pt", "N") == pytest.approx(2.051)
    assert MT._check_ml_bond_length("Pt", "N") == before                # a check: historic
    assert PH.md_distance("Pt", "N") == before_md                       # tolerance window: historic
    from delfin.manta import topology_checks as TC
    assert TC._get_ml_bond_length("Pt", "N") == before                  # the checks module reads checks


def test_check_callers_import():
    from delfin.manta import hybrid_fragments, secondary_metal_modules  # noqa: F401
    assert callable(secondary_metal_modules._check_ml_bond_length)
