"""Asked which dye absorbs more in the visible, the operator had the
transition table, chose a threshold and a range by itself, filtered by
hand, and flagged the choice as an assumption in its report -- and named
a field for the brightest visible line as the one change (2026-09-11).
The definitions are stated in the result and the lines derived from
them: first bright, brightest, brightest visible.
"""

from __future__ import annotations

import json

import pytest

from delfin import api
from delfin.ops_server import server as ops

_TDDFT = """\
ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-------------------------------------------------------------------
States    Energy (eV)   Energy (cm-1)   Wavelength (nm)   fosc  D**2  X     Y     Z
-------------------------------------------------------------------
  0-1A  ->  1-3A    2.200000     17744.2     563.6     0.000000    0.0    0.0  0.0  0.0
  0-1A  ->  2-3A    2.650000     21373.7     467.9     0.412300    0.5    0.1  0.2  0.4
  0-1A  ->  3-3A    3.400000     27422.8     364.7     0.021000    1.2    0.0  0.5  1.0
  0-1A  ->  4-3A    4.200000     33875.3     295.2     0.900000    1.2    0.0  0.5  1.0

"""


def _run(tmp_path, text=_TDDFT):
    d = tmp_path / "dye"
    d.mkdir()
    (d / "run.out").write_text(text)
    return str(d)


def test_the_first_bright_line_skips_the_dark_one(tmp_path):
    res = api.extract_excited_states(_run(tmp_path))
    assert res.first_bright["index"] == 1 and res.first_bright["wavelength_nm"] == 467.9
    assert res.first_bright["in_visible"] is True


def test_the_brightest_and_the_brightest_visible_differ_when_the_uv_line_is_stronger(tmp_path):
    res = api.extract_excited_states(_run(tmp_path))
    assert res.brightest["index"] == 3 and res.brightest["fosc"] == 0.9     # 295 nm, UV
    assert res.brightest["in_visible"] is False
    assert res.brightest_visible["index"] == 1                              # 468 nm, 0.41


def test_the_definitions_travel_with_the_answer(tmp_path):
    out = json.loads(ops.tool_extract_excited_states(_run(tmp_path)))
    assert out["bright_threshold_fosc"] == 0.01
    assert list(out["visible_range_nm"]) == [380.0, 780.0]
    assert out["brightest_visible"]["wavelength_nm"] == 467.9
    assert "brightest_visible" in (ops.tool_extract_excited_states.__doc__ or "")


def test_no_bright_line_is_null_not_a_guess(tmp_path):
    dark = _TDDFT.replace("0.412300", "0.000000").replace("0.021000", "0.000000").replace("0.900000", "0.004000")
    res = api.extract_excited_states(_run(tmp_path, dark))
    assert res.n_states == 4
    assert res.first_bright is None and res.brightest is None and res.brightest_visible is None


# --- the overview over many folders --------------------------------------------

_ORB = """
----------------
ORBITAL ENERGIES
----------------

  NO   OCC          E(Eh)            E(eV)
   0   2.0000     -19.246153      -523.7146
   1   2.0000      -0.215000        -5.8505
   2   0.0000      -0.085000        -2.3130
   3   0.0000       0.412345        11.2200

------------------
"""
_FREQ = """
-----------------------
VIBRATIONAL FREQUENCIES
-----------------------
   0:    {m0} cm**-1{imag}
   1:      0.00 cm**-1
   2:    250.50 cm**-1

NORMAL MODES
------------
"""


def _spectra_ws(tmp_path):
    root = tmp_path / "spectra"
    dye = root / "dye"
    dye.mkdir(parents=True)
    (dye / "run.inp").write_text("! B3LYP def2-SVP TDDFT CPCM(water)\n")
    (dye / "run.out").write_text(_ORB + _TDDFT + _FREQ.format(m0="45.20", imag="") + "****ORCA TERMINATED NORMALLY****\n")
    ts = root / "ts_guess"
    ts.mkdir()
    (ts / "run.inp").write_text("! B3LYP def2-SVP Opt Freq\n")
    (ts / "run.out").write_text(_ORB + _FREQ.format(m0="-312.40", imag=" ***imaginary mode***") + "****ORCA TERMINATED NORMALLY****\n")
    return root


def test_the_overview_answers_all_three_questions_in_one_call(tmp_path):
    root = _spectra_ws(tmp_path)
    rows = {r["folder"].rsplit("/", 1)[-1]: r for r in api.extract_spectra_table([str(root / "dye"), str(root / "ts_guess"), str(root / "nowhere")])}
    dye, ts, gone = rows["dye"], rows["ts_guess"], rows["nowhere"]
    assert dye["gap_ev"] == pytest.approx(3.5375, abs=1e-3)
    assert dye["first_bright_nm"] == 467.9 and dye["first_bright_fosc"] == 0.4123
    assert dye["brightest_visible_nm"] == 467.9
    assert dye["n_imag"] == 0 and dye["is_minimum"] is True
    assert dye["method"] == "B3LYP/def2-SVP/water"
    assert ts["n_imag"] == 1 and ts["is_minimum"] is False and ts["most_negative_cm"] == pytest.approx(-312.4)
    assert ts["first_bright_nm"] is None and any("excited states" in n for n in ts["notes"])
    assert gone["gap_ev"] is None and "folder missing" in gone["notes"]


def test_the_overview_tool_is_relative_to_a_root_and_callable(tmp_path):
    root = _spectra_ws(tmp_path)
    from delfin.agent.api_client import _MCP_READONLY_TOOL_BASES
    out = json.loads(ops.tool_extract_spectra_table(f"{root / 'dye'},{root / 'ts_guess'}"))
    assert out["root"] == str(root)
    assert [r["folder"] for r in out["rows"]] == ["dye", "ts_guess"]
    assert "extract_spectra_table" in _MCP_READONLY_TOOL_BASES
    assert any(e["name"] == "extract_spectra_table" for e in api._TOOL_CATALOG)


# --- every value knows its line ---------------------------------------------------------


def _line_of(path, needle):
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        if needle in line:
            return i
    raise AssertionError(needle)


def test_transitions_and_orbitals_know_their_lines(tmp_path):
    root = _spectra_ws(tmp_path)
    out = root / "dye" / "run.out"
    exc = api.extract_excited_states(str(root / "dye"))
    assert exc.first_bright["line"] == _line_of(out, "2-3A    2.650000")
    assert exc.transitions[0].line == _line_of(out, "1-3A    2.200000")
    orb = api.extract_orbital_energies(str(root / "dye"))
    assert orb.homo_line == _line_of(out, "   1   2.0000      -0.215000")
    assert orb.lumo_line == _line_of(out, "   2   0.0000      -0.085000")
    assert orb.block_line == _line_of(out, "ORBITAL ENERGIES")
    imag = api.extract_imaginary_frequencies(str(root / "ts_guess"))
    assert imag.most_negative_line == _line_of(root / "ts_guess" / "run.out", "-312.40 cm**-1")


def test_the_overview_says_whether_the_first_bright_line_is_visible(tmp_path):
    root = _spectra_ws(tmp_path)
    uv = root / "uv_dye"
    uv.mkdir()
    (uv / "run.inp").write_text("! B3LYP def2-SVP TDDFT\n")
    (uv / "run.out").write_text(_ORB + _TDDFT.replace("0.412300", "0.000000").replace("0.021000", "0.000000") + "****ORCA TERMINATED NORMALLY****\n")
    rows = {r["folder"].rsplit("/", 1)[-1]: r for r in api.extract_spectra_table([str(root / "dye"), str(uv)])}
    assert rows["dye"]["first_bright_in_visible"] is True and rows["dye"]["first_bright_line"]
    assert rows["uv_dye"]["first_bright_in_visible"] is False        # the 295 nm line, fosc 0.9
    assert rows["uv_dye"]["brightest_visible_nm"] is None
    assert any("no bright transition in the visible" in n for n in rows["uv_dye"]["notes"])
