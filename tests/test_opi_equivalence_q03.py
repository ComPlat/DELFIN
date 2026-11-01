"""DELFIN's ORCA IR reader against OPI's, on synthetic ORCA 6.1.1 outputs.

DELFIN's reader of record for IR is :func:`delfin.ir_spectrum.parse_ir_spectrum`,
which greps the ``IR SPECTRUM`` block out of the ``.out`` text.  ORCA's own
structured interface OPI (``orca-pi`` 2.0.0) reads the *same* printed block
through ``Output.get_ir()``.  The equivalence test compares both on a synthetic
ORCA output written into ``tmp_path`` -- real ORCA output is never committed
(# SYNTHETIC rule); the values were recorded from a real H2O opt+freq run and
hand-copied into the block below.

Measured (see scratch/report_q03.md):
  * on the hand-recorded IR block ``parse_ir_spectrum`` and OPI ``get_ir``
    agree value-for-value on mode count, wavenumber (cm^-1), eps
    (L/(mol*cm)) and intensity (km/mol); both read the same printed lines, so
    the tolerance is the printed precision.
  * on a synthetic NOT-converged run there is no IR block: DELFIN answers
    ``[]`` and OPI ``None``, and OPI's ``scf_converged()`` is False.  A value
    from an unconverged run is meaningless, and neither reader hands one back.

OPI has NO reader for the .hess frequency list, the VIBRATIONAL FREQUENCIES
section or the ***imaginary mode*** lines, so delfin.imag's readers have no
OPI counterpart (they keep the text as reader of record; see the report).

CI does not install the ``[orca]`` extra, so the real-OPI comparison is gated
with ``pytest.importorskip("opi")``, following tests/test_opi_is_a_witness_
when_it_is_there.py.  The deterministic reader test uses no OPI.
"""

from __future__ import annotations

import pytest

from delfin.ir_spectrum import parse_ir_spectrum

# Tolerances in the units ORCA prints (wavenumber cm^-1, eps L/(mol*cm),
# intensity km/mol).  ORCA prints 2, 6 and 2 significant decimals
# respectively; a reader that only agrees to the printed precision is still
# the same value, so these are tighter than any physical uncertainty.
WAVENUMBER_TOL = 0.01        # cm^-1
EPS_TOL = 1e-6               # L/(mol*cm)
INTENSITY_TOL = 0.01         # km/mol

# The IR block as ORCA 6.1.1 prints it for the 3 IR-active vibrations of
# water (recorded from a real opt+freq run; the banner above it has been
# shortened, the mode lines are verbatim).
IR_BLOCK = (
    "-----------\n"
    "IR SPECTRUM\n"
    "-----------\n"
    "\n"
    " Mode   freq       eps      Int      T**2         TX        TY        TZ\n"
    "       cm**-1   L/(mol*cm) km/mol    a.u.\n"
    "----------------------------------------------------------------------------\n"
    "   6:   1637.64   0.010942   55.30  0.002085  ( 0.000000  0.000000  0.045663)\n"
    "   7:   3787.07   0.000924    4.67  0.000076  (-0.000000  0.000000  0.008728)\n"
    "   8:   3882.05   0.005245   26.51  0.000422  ( 0.000000 -0.020533  0.000000)\n"
    "\n"
)

# Minimal header the way ORCA 6.1.1 prints it: the readers that OPI and
# DELFIN share only need the banner, the SCF SUCCESS marker and the IR block.
SUCCESSFUL_FREQ_RUN = (
    "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
    "\n"
    "                                 *****************\n"
    "                                 * O   R   C   A *\n"
    "                                 *****************\n"
    "\n"
    "                         Program Version 6.1.1  -  RELEASE   -\n"
    "\n"
    "-----------------------\n"
    "ORCA SCF CALCULATIONS\n"
    "-----------------------\n"
    "SUCCESS\n"
    "\n"
    + IR_BLOCK
    + "\n"
    "                   **** ORCA TERMINATED NORMALLY ****\n"
)

NOT_CONVERGED_RUN = (
    "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
    "\n"
    "                         Program Version 6.1.1  -  RELEASE   -\n"
    "\n"
    "-----------------------\n"
    "ORCA SCF CALCULATIONS\n"
    "-----------------------\n"
    "The SCF was not converged in 200 iterations\n"
    "NOT CONVERGED\n"
    "\n"
)


def test_ir_readers_give_the_same_values_from_the_same_lines(tmp_path):
    """Both readers parse the printed 'IR SPECTRUM' block value-for-value.

    Runs DELFIN's reader against the hand-recorded block (no OPI needed), so
    the tolerance rationale is checked even where OPI is not installed.
    """
    out = tmp_path / "h2o_freq.out"
    out.write_text(SUCCESSFUL_FREQ_RUN)

    modes = parse_ir_spectrum(out)
    assert len(modes) == 3
    expected = {
        6: (1637.64, 0.010942, 55.30),
        7: (3787.07, 0.000924, 4.67),
        8: (3882.05, 0.005245, 26.51),
    }
    by_mode = {m.mode_number: m for m in modes}
    for mode_no, (freq, eps, intens) in expected.items():
        m = by_mode[mode_no]
        assert abs(m.frequency_cm1 - freq) < WAVENUMBER_TOL
        assert abs(m.epsilon - eps) < EPS_TOL
        assert abs(m.intensity_km_mol - intens) < INTENSITY_TOL


def test_delfin_and_opi_agree_on_the_synthetic_run(tmp_path):
    """DELFIN's text reader and OPI's get_ir return the same modes/values.

    The synthetic output is written into tmp_path (never tests/fixtures/),
    because OPI constructs an Output with a working_dir and must not write
    anywhere near the repository.
    """
    pytest.importorskip("opi")

    from opi.output.core import Output

    out = tmp_path / "h2o_freq.out"
    out.write_text(SUCCESSFUL_FREQ_RUN)

    delfin_modes = {m.mode_number: m for m in parse_ir_spectrum(out)}
    opi = Output(basename=out.stem, working_dir=tmp_path)
    opi_modes = opi.get_ir()

    assert opi_modes, "OPI must read the synthetic IR block"
    assert len(delfin_modes) == len(opi_modes)
    for mode_no, opi_m in opi_modes.items():
        d_m = delfin_modes[mode_no]
        assert abs(d_m.frequency_cm1 - opi_m.wavenumber) < WAVENUMBER_TOL, (
            f"mode {mode_no} wavenumber: DELFIN={d_m.frequency_cm1} "
            f"OPI={opi_m.wavenumber}"
        )
        assert abs(d_m.epsilon - opi_m.eps) < EPS_TOL, (
            f"mode {mode_no} eps: DELFIN={d_m.epsilon} OPI={opi_m.eps}"
        )
        assert abs(d_m.intensity_km_mol - opi_m.intensity) < INTENSITY_TOL, (
            f"mode {mode_no} intensity: DELFIN={d_m.intensity_km_mol} "
            f"OPI={opi_m.intensity}"
        )


def test_not_converged_run_reports_no_ir_to_either_reader(tmp_path):
    """A run whose SCF did not converge has no IR block: DELFIN says [] and
    OPI says None, and OPI's convergence flag says the run is bad -- so no
    reader can hand back a bogus mode list."""
    out = tmp_path / "scf_failed.out"
    out.write_text(NOT_CONVERGED_RUN)

    assert parse_ir_spectrum(out) == [], (
        "DELFIN must not invent IR modes for an unconverged run"
    )

    try:
        from opi.output.core import Output
    except Exception:
        # OPI absent (CI without the [orca] extra): the DELFIN half above
        # still holds.
        return

    opi = Output(basename=out.stem, working_dir=tmp_path)
    assert opi.scf_converged() is False
    assert opi.get_ir() is None
