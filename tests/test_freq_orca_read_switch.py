"""Call-site switch tests: frequency readers go through delfin.orca_read.

The readers delfin.ir_spectrum.parse_ir_spectrum,
delfin.imag.collect_imaginary_modes and delfin.imag.hessian_frequencies
(the leaves every other frequency reader -- saddle_reason,
imaginary_modes, search_imaginary_mode2, cli_ir_report,
reporting/delfin_collector -- calls) route through delfin.orca_read with
their original bodies as the text fallback. Behaviour is proven
identical by the baseline recorder (scratch/freq_baseline.json, 207
fixture+archive runs, bitwise comparison); ``orca.reader=text`` must
restore the old path exactly.

The ORCA-shaped snippets below are hand-written synthetic test fixtures
(first line carries the DELFIN synthetic marker), not ORCA output.
"""
import pytest

from delfin import imag, ir_spectrum, orca_read

IR_BLOCK = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
-------------------
IR SPECTRUM
-------------------
 Mode freq eps Int T**2 TX TY TZ
-------------------
  6:   89.55  0.000058  0.30  0.000203  (-0.002326 -0.009616 -0.010276)
  7:  113.29  0.000190  0.96  0.000522  ( 0.012132 -0.002599 -0.019194)
  8: 1641.13 12.345600 55.00  1.234500  ( 0.500000  0.600000  0.700000)
-------------------
"""

IMAG_LINES = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
     4:       -831.00 cm**-1   ***imaginary mode***
     5:       -12.50 cm**-1   ***imaginary mode***
     6:        89.55 cm**-1
"""

HESS = """\n$vibrational_frequencies
9
    0                3630.278307
    1                3749.638372
    2                3865.117529
    3                 1595.070406
    4                1595.070406
    5                1595.070406
    6                 -831.000000
    7                  -12.500000
    8                   89.550000
"""


def _run(tmp_path, out_text=IR_BLOCK, hess_text=HESS, base="run"):
    d = tmp_path / "run"
    d.mkdir()
    (d / f"{base}.out").write_text(out_text)
    (d / f"{base}.hess").write_text(hess_text)
    return d


def test_parse_ir_spectrum_routes_through_orca_read(tmp_path):
    d = _run(tmp_path)
    with pytest.MonkeyPatch.context() as mp:
        called = {}

        def fake_read(work_dir, base_name, quantity):
            called["q"] = quantity
            return orca_read.Reading(value=[], unit="cm-1", source="text")

        mp.setattr(orca_read, "read", fake_read)
        modes = ir_spectrum.parse_ir_spectrum(d / "run.out")
    assert called["q"] is orca_read.Quantity.IR_MODES
    assert modes == []


def test_parse_ir_spectrum_falls_back_when_orca_read_has_none(tmp_path):
    d = _run(tmp_path)

    def fake_none(work_dir, base_name, quantity):
        return orca_read.Reading(value=None, unit="", source="none", note="x")

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(orca_read, "read", fake_none)
        modes = ir_spectrum.parse_ir_spectrum(d / "run.out")
    # fallback = the original text parser -> full mode list
    assert [m.mode_number for m in modes] == [6, 7, 8]
    assert modes[0].intensity_km_mol == pytest.approx(0.30)


def test_collect_imaginary_modes_routes_through_orca_read(tmp_path):
    d = _run(tmp_path, out_text=IMAG_LINES)
    with pytest.MonkeyPatch.context() as mp:
        called = {}

        def fake_read(work_dir, base_name, quantity):
            called["q"] = quantity
            return orca_read.Reading(value=[], unit="cm-1", source="text")

        mp.setattr(orca_read, "read", fake_read)
        modes = imag.collect_imaginary_modes(str(d / "run.out"))
    assert called["q"] is orca_read.Quantity.IMAGINARY_MODES
    assert modes == []


def test_collect_imaginary_modes_falls_back(tmp_path):
    d = _run(tmp_path, out_text=IMAG_LINES)

    def fake_none(work_dir, base_name, quantity):
        return orca_read.Reading(value=None, unit="", source="none", note="x")

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(orca_read, "read", fake_none)
        modes = imag.collect_imaginary_modes(str(d / "run.out"))
    assert modes == [(4, -831.0), (5, -12.5)]


def test_hessian_frequencies_routes_through_orca_read(tmp_path):
    d = _run(tmp_path)
    with pytest.MonkeyPatch.context() as mp:
        called = {}

        def fake_read_hess(hess_path, quantity):
            called["q"] = quantity
            return orca_read.Reading(value=[], unit="cm-1", source="text")

        mp.setattr(orca_read, "read_hess", fake_read_hess)
        freqs = imag.hessian_frequencies(d / "run.hess")
    assert called["q"] is orca_read.Quantity.VIBRATIONAL_FREQUENCIES
    assert freqs == []


def test_hessian_frequencies_falls_back(tmp_path):
    d = _run(tmp_path)

    def fake_none(hess_path, quantity):
        return orca_read.Reading(value=None, unit="", source="none", note="x")

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(orca_read, "read_hess", fake_none)
        freqs = imag.hessian_frequencies(d / "run.hess")
    assert freqs[6] == -831.0 and freqs[7] == -12.5 and freqs[8] == 89.55


def test_end_to_end_matches_the_old_bodies(tmp_path):
    """Through the switched readers the values are the old parsers' values."""
    d = _run(tmp_path)
    ir_modes = ir_spectrum.parse_ir_spectrum(d / "run.out")
    assert [(m.frequency_cm1, m.intensity_km_mol) for m in ir_modes] == \
        [(89.55, 0.30), (113.29, 0.96), (1641.13, 55.00)]
    imag_modes = imag.collect_imaginary_modes(str(d / "run.out"))
    assert imag_modes == [(4, -831.0), (5, -12.5)]
    freqs = imag.hessian_frequencies(d / "run.hess")
    assert freqs[6] == -831.0


def test_reader_text_restores_the_old_path_exactly(tmp_path):
    """orca.reader=text: the readers answer from the old parsers alone."""
    d = _run(tmp_path)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(orca_read, "_reader_mode", lambda: "text")
        ir = ir_spectrum.parse_ir_spectrum(d / "run.out")
        imag_modes = imag.collect_imaginary_modes(str(d / "run.out"))
        freqs = imag.hessian_frequencies(d / "run.hess")
    assert [(m.frequency_cm1, m.intensity_km_mol) for m in ir] == \
        [(89.55, 0.30), (113.29, 0.96), (1641.13, 55.00)]
    assert imag_modes == [(4, -831.0), (5, -12.5)]
    assert freqs[6] == -831.0 and freqs[8] == 89.55
    # and the reading still comes from the text reader, not OPI
    reading = orca_read.read(str(d), "run", orca_read.Quantity.IR_MODES)
    assert reading.source == "text"


def test_read_reports_ir_modes_and_imag_modes(tmp_path):
    d = _run(tmp_path)
    r = orca_read.read(str(d), "run", orca_read.Quantity.IR_MODES)
    assert r.source == "text"
    assert [(m.mode_number, m.intensity_km_mol) for m in r.value] == \
        [(6, 0.30), (7, 0.96), (8, 55.00)]
    d2 = _run(tmp_path, out_text=IMAG_LINES, base="imag")
    r2 = orca_read.read(str(d2), "imag", orca_read.Quantity.IMAGINARY_MODES)
    assert r2.value == [(4, -831.0), (5, -12.5)]


def test_read_hess_reports_vibrational_frequencies(tmp_path):
    d = _run(tmp_path)
    r = orca_read.read_hess(d / "run.hess",
                            orca_read.Quantity.VIBRATIONAL_FREQUENCIES)
    assert r.source == "text"
    assert r.value[6] == -831.0


def test_missing_run_never_raises(tmp_path):
    r = orca_read.read(str(tmp_path), "nothing",
                       orca_read.Quantity.IR_MODES)
    assert r.value is None and r.source == "none"
    # missing .out: the switched reader falls back to the old parser
    # (which logs and answers [])
    assert ir_spectrum.parse_ir_spectrum(tmp_path / "nothing.out") == []
    assert imag.collect_imaginary_modes(str(tmp_path / "nothing.out")) == []
