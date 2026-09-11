"""Asked which dye absorbs more in the visible, the operator had the
transition table, chose a threshold and a range by itself, filtered by
hand, and flagged the choice as an assumption in its report -- and named
a field for the brightest visible line as the one change (2026-09-11).
The definitions are stated in the result and the lines derived from
them: first bright, brightest, brightest visible.
"""

from __future__ import annotations

import json

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
