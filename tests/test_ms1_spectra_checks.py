"""MS1: /check reads spectra and flags physically implausible results.

Fixtures are abridged but real-shape ORCA output fragments, following the
precedent of tests/test_result_critic.py and the recorded fixture
tests/fixtures/method_reviewer/water_b3lyp/water.out. Every check must
have a red test (fires) AND a does-not-fire control.

Threshold sources (also in delfin/agent/result_critic.py docstrings):
- Negative IR intensity: Beer-Lambert absorption cross-sections are
  physically non-negative for an absorption spectrum (physical law).
- Mode count 3N-6 / 3N-5: standard vibrational analysis, ORCA manual
  6.1.1 section 4.6 "Vibrational Frequencies".
- <S^2> deviation for excited states: ORCA manual TDDFT output
  description; DELFIN methodology "Spin Contamination and State
  Selection" uses |<S^2> - S(S+1)| as the diagnostic.

Report-only: /check never edits an output, never re-runs, never
loosens a criterion (pinned by tests/test_sc1_thresholds.py).
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import result_critic as rc


_TERM = "ORCA TERMINATED NORMALLY\n"
_ENERGY = "FINAL SINGLE POINT ENERGY     -76.320314120\n"
_ATOMS = "Number of atoms                             ...      3\n"
_FREQ_HDR = "VIBRATIONAL FREQUENCIES\n-----------------------\n"


def _codes(crits):
    return {c.code for c in crits}


def _flag(crits, code):
    return next(c for c in crits if c.code == code)


def _write(tmp_path, name, content) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


# --- recorded-shape fixtures ------------------------------------------------

_WATER_MIN = (
    _TERM + _ENERGY + _ATOMS
    + _FREQ_HDR
    + "   0:         0.00 cm**-1\n"
    + "   1:         0.00 cm**-1\n"
    + "   2:         0.00 cm**-1\n"
    + "   3:         0.00 cm**-1\n"
    + "   4:         0.00 cm**-1\n"
    + "   5:         0.00 cm**-1\n"
    + "   6:      1625.31 cm**-1\n"
    + "   7:      3727.90 cm**-1\n"
    + "   8:      3826.44 cm**-1\n"
    + "\nIR SPECTRUM\n"
    + "-----------\n"
    + " Mode    freq      eps      Int      T**2     TX     TY     TZ\n"
    + "-----------------------------------------------------------------\n"
    + "   6:   1625.31   0.0345    12.30   0.0023   0.12   0.00   0.00\n"
    + "   7:   3727.90   0.5123    45.10   0.0102   0.34   0.01   0.00\n"
    + "   8:   3826.44   0.6011    52.40   0.0110   0.37   0.01   0.00\n"
    + "\nThe first frequency in this output\n"
)


def test_nonlinear_minimum_with_correct_mode_count_does_not_fire(tmp_path):
    """Water, 3 atoms: 3N-6 = 3 real modes; 6 trans/rot zeros excluded.

    Control: the count check must NOT fire on a correct minimum.
    """
    out = _write(tmp_path, "water.out", _WATER_MIN)
    crits = rc.critique_output(out)
    assert "freq-count-mismatch" not in _codes(crits)


def test_nonlinear_minimum_with_wrong_real_mode_count_fires(tmp_path):
    """7 real modes on a 3-atom non-linear minimum (3N-6=6) must flag.

    The count check excludes the 6 zero trans/rot modes before counting.
    """
    bad = _WATER_MIN.replace(
        "   8:      3826.44 cm**-1\n",
        "   8:      3826.44 cm**-1\n   9:      3900.00 cm**-1\n")
    out = _write(tmp_path, "bad.out", bad)
    crits = rc.critique_output(out)
    flag = _flag(crits, "freq-count-mismatch")
    assert flag.level == "error"
    assert "3N-6" in flag.message or "3N" in flag.message


def test_linear_molecule_marker_switches_to_3n_minus_5(tmp_path):
    """A linear 3-atom molecule: 3N-5 = 4 real modes expected.

    CO2-shaped output: ZERO first rotational constant in the
    thermochemistry block (ORCA manual 6.1.1 §4.7 — the real marker;
    s4's F1: no literal "linear molecule" phrase exists in ORCA
    output), 4 real modes, no fire.
    """
    lin = (_TERM + _ENERGY
           + "Number of atoms                             ...      3\n"
           + "Rotational constants in cm-1: 0.000000 0.390721 0.390721\n"
           + _FREQ_HDR
           + "   0:         0.00 cm**-1\n"
           + "   1:         0.00 cm**-1\n"
           + "   2:         0.00 cm**-1\n"
           + "   3:         0.00 cm**-1\n"
           + "   4:      1200.00 cm**-1\n"
           + "   5:      2400.00 cm**-1\n"
           + "   6:      2400.00 cm**-1\n"
           + "   7:      3800.00 cm**-1\n")
    out = _write(tmp_path, "co2.out", lin)
    crits = rc.critique_output(out)
    assert "freq-count-mismatch" not in _codes(crits)


def test_negative_ir_intensity_is_an_error(tmp_path):
    """A negative Int (km/mol) is unphysical for an absorption spectrum
    (Beer-Lambert: cross-section >= 0). Must be flagged as an error.
    """
    bad = _WATER_MIN.replace(
        "   7:   3727.90   0.5123    45.10",
        "   7:   3727.90   0.5123   -45.10")
    out = _write(tmp_path, "negir.out", bad)
    crits = rc.critique_output(out)
    assert "ir-neg-intensity" in _codes(crits)
    assert _flag(crits, "ir-neg-intensity").level == "error"


def test_positive_ir_intensities_do_not_fire(tmp_path):
    """Control: all-positive intensities (water fixture) stay clean."""
    out = _write(tmp_path, "ok.out", _WATER_MIN)
    crits = rc.critique_output(out)
    assert "ir-neg-intensity" not in _codes(crits)


def test_zero_ir_intensity_is_reported_not_flagged_as_error(tmp_path):
    """A zero intensity is physical (symmetry-forbidden mode) — reported
    as an ok-level note, never as an error or warning.
    """
    zero = _WATER_MIN.replace(
        "   6:   1625.31   0.0345    12.30",
        "   6:   1625.31   0.0345     0.00")
    out = _write(tmp_path, "zero.out", zero)
    crits = rc.critique_output(out)
    assert "ir-neg-intensity" not in _codes(crits)


def test_imaginary_mode_listed_in_ir_spectrum_is_flagged(tmp_path):
    """An imaginary mode appearing as a spectrum peak (negative freq row
    in the IR SPECTRUM table) must be flagged: peaks from imaginary
    modes are not physical lines.
    """
    imag_ir = _WATER_MIN.replace(
        "   6:   1625.31   0.0345    12.30   0.0023   0.12   0.00   0.00\n",
        "   6:  -1625.31   0.0345    12.30   0.0023   0.12   0.00   0.00\n")
    out = _write(tmp_path, "imag.out", imag_ir)
    crits = rc.critique_output(out)
    assert "ir-imaginary-as-peak" in _codes(crits)


def test_ir_without_imaginary_rows_does_not_fire(tmp_path):
    """Control: the clean water fixture has no negative-freq IR row."""
    out = _write(tmp_path, "ok.out", _WATER_MIN)
    crits = rc.critique_output(out)
    assert "ir-imaginary-as-peak" not in _codes(crits)

# ---------------------------------------------------------------------------
# Wiring pin (reviewer s4, point 2): the folder parsers read the LARGEST
# .out. /check critiques ONE file. Choice (b): the spectra checks use the
# folder parsers ONLY when the critiqued file is the largest .out in its
# folder; otherwise a visible "spectra-skipped" warning instead of a
# silent divergence between "file flagged" and "file parsed".
# ---------------------------------------------------------------------------


def test_spectra_skipped_when_critiqued_file_is_not_largest(tmp_path):
    big = _write(tmp_path, "big.out", _WATER_MIN + "\n" * 40)
    small = _write(tmp_path, "small.out", _WATER_MIN)
    crits = rc.critique_output(small)
    flag = _flag(crits, "spectra-skipped")
    assert "big.out" in flag.message
    assert not any(c.code.startswith(("ir-", "freq-count", "uv-"))
                   for c in crits)
    # The larger file still gets its own critique when scanned directly.
    assert "spectra-skipped" not in _codes(rc.critique_output(big))
