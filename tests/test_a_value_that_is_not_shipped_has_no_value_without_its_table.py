"""Calibration values that MANTA does not ship are absent without their table.

Four default-off construction paths were calibrated on structure statistics that are
not distributed with DELFIN: per-metal M-C lengths, per-period donor angles with the
measured ring angles and the M-D-X slope, and the element-normalised aromaticity
threshold.  Their numbers are no longer in the code; each path reads them from a file
named by an environment variable (``delfin.manta._unshipped_tables``).  These tests pin
the contract: without the file the path has no value and falls back to its open
reference -- even with its switch on -- and with a file it uses exactly what the file
says.  The numbers written into the temporary tables are arbitrary probe values.
"""
from __future__ import annotations

import math

import pytest

from delfin.manta import _arom_criterion as AC
from delfin.manta import _energy_terms as ET
from delfin.manta import ml_tables as MT
from delfin.manta._unshipped_tables import scalar_table


def _write_probe_table(tmp_path, name, rows):
    path = tmp_path / name
    path.write_text("".join(f"{k}\t{v}\n" for k, v in rows))
    return str(path)


def test_an_unset_variable_gives_an_empty_table(monkeypatch):
    monkeypatch.delenv("DELFIN_TEST_UNSHIPPED", raising=False)
    assert scalar_table("DELFIN_TEST_UNSHIPPED") == {}


def test_a_missing_file_gives_an_empty_table(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_TEST_UNSHIPPED", str(tmp_path / "absent.tsv"))
    assert scalar_table("DELFIN_TEST_UNSHIPPED") == {}


def test_the_table_skips_comments_headers_and_short_rows(monkeypatch, tmp_path):
    path = tmp_path / "t.tsv"
    path.write_text("# comment\nkey\tvalue\n\nonly_a_key\nA\t1.5\textra\nB\t2\n")
    monkeypatch.setenv("DELFIN_TEST_UNSHIPPED", str(path))
    assert scalar_table("DELFIN_TEST_UNSHIPPED") == {"A": 1.5, "B": 2.0}


def test_the_crystal_tables_are_not_in_the_code():
    for module, name in ((MT, "_MC_LEN_CRYSTAL_SOURCE"), (MT, "_MC_LEN_CRYSTAL_MEDIAN"),
                         (ET, "_THETA_BY_PERIOD"), (ET, "_RING_INTERNAL_MEAS"),
                         (ET, "_RING_EXOCYCLIC_MEAS"), (ET, "_PHI_MD_SLOPE"),
                         (AC, "_AROM_RATIO_MAX")):
        assert not hasattr(module, name), f"{module.__name__}.{name} is back in the code"


def test_m_c_length_with_its_switch_on_but_no_table_is_the_switch_off_value(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_MC_LEN_CRYSTAL", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_MC_LEN_TABLE", raising=False)
    off = {m: MT._get_ml_bond_length(m, "C") for m in ("Mn", "Os", "Ti", "La", "Ge", "Fe")}
    monkeypatch.setenv("DELFIN_FFFREE_MC_LEN_CRYSTAL", "1")
    on = {m: MT._get_ml_bond_length(m, "C") for m in off}
    assert on == off


def test_m_c_length_comes_from_the_supplied_table(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_FFFREE_MC_LEN_CRYSTAL", "1")
    monkeypatch.setenv("DELFIN_FFFREE_MC_LEN_TABLE",
                       _write_probe_table(tmp_path, "mc.tsv", [("Os", 1.234)]))
    monkeypatch.delenv("DELFIN_FFFREE_MC_LEN_METALS", raising=False)
    assert MT._get_ml_bond_length("Os", "C") == pytest.approx(1.234)
    monkeypatch.setenv("DELFIN_FFFREE_MC_LEN_METALS", "Re")
    assert MT._get_ml_bond_length("Os", "C") != pytest.approx(1.234)


def test_donor_angles_without_a_table_are_tetrahedral_and_polygonal(monkeypatch):
    monkeypatch.delenv("DELFIN_FFREE_THETA_TABLE", raising=False)
    for period in (2, 3, 4, 5, 6, 7):
        assert ET._period_theta(period) == pytest.approx(109.5)
    assert ET._ring_internal_deg(5) == pytest.approx(108.0)
    assert ET._ring_internal_deg(6) == pytest.approx(120.0)
    assert ET._ring_exocyclic_deg(5, 108.0) == pytest.approx(126.0)
    assert ET._phi_md_slope() == 0.0


def test_donor_angles_come_from_the_supplied_table(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_FFREE_THETA_TABLE", _write_probe_table(tmp_path, "theta.tsv", [
        ("period_3", 101.0), ("ring_internal_6", 117.0), ("ring_exocyclic_6", 121.5),
        ("phi_md_slope", 0.5)]))
    assert ET._period_theta(3) == pytest.approx(101.0)
    assert ET._period_theta(4) == pytest.approx(109.5)
    assert ET._ring_internal_deg(6) == pytest.approx(117.0)
    assert ET._ring_exocyclic_deg(6, 117.0) == pytest.approx(121.5)
    assert ET._phi_md_slope() == pytest.approx(0.5)


def test_signature_angle_of_a_period_3_donor_follows_the_table(monkeypatch, tmp_path):
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.AddHs(Chem.MolFromSmiles("CS(C)[Cu]"))
    s = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "S")
    c1, c2 = [n.GetIdx() for n in mol.GetAtomWithIdx(s).GetNeighbors()
              if n.GetSymbol() == "C"]
    monkeypatch.delenv("DELFIN_FFREE_THETA_TABLE", raising=False)
    assert math.degrees(ET._signature_theta(mol, c1, s, c2)) == pytest.approx(109.5)
    monkeypatch.setenv("DELFIN_FFREE_THETA_TABLE",
                       _write_probe_table(tmp_path, "theta.tsv", [("period_3", 101.0)]))
    assert math.degrees(ET._signature_theta(mol, c1, s, c2)) == pytest.approx(101.0)


_THIOPHENE = (["S", "C", "C", "C", "C"], [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
              [1.71, 1.37, 1.42, 1.37, 1.71])


def test_aromaticity_without_a_threshold_keeps_the_absolute_gate(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_AROM_RATIO_TABLE", raising=False)
    monkeypatch.setenv("DELFIN_FFFREE_AROM_CRITERION_RADII", "0")
    off = AC.ring_rejected_by_length(*_THIOPHENE)
    monkeypatch.setenv("DELFIN_FFFREE_AROM_CRITERION_RADII", "1")
    on = AC.ring_rejected_by_length(*_THIOPHENE)
    assert off is True and on is True


def test_aromaticity_uses_the_supplied_threshold(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_FFFREE_AROM_CRITERION_RADII", "1")
    monkeypatch.setenv("DELFIN_FFFREE_AROM_RATIO_TABLE",
                       _write_probe_table(tmp_path, "arom.tsv", [("ratio_max", 0.95)]))
    assert AC.ring_rejected_by_length(*_THIOPHENE) is False
    monkeypatch.setenv("DELFIN_FFFREE_AROM_RATIO_TABLE",
                       _write_probe_table(tmp_path, "arom2.tsv", [("ratio_max", 0.90)]))
    assert AC.ring_rejected_by_length(*_THIOPHENE) is True
