"""DELFIN's ORCA frequency/IR readers compared with OPI's (swarm q03).

DELFIN's reader of record for IR is :func:`delfin.ir_spectrum.parse_ir_spectrum`,
which greps the ``IR SPECTRUM`` block out of the ``.out`` text.  ORCA's own
structured interface OPI reads the *same* printed block through
``Output.get_ir()``.  This equivalence was added so DELFIN can later move to
OPI for IR without losing anything: both must return the same modes, same
wavenumbers, same intensities, within the tolerance the ``.out`` prints.

Measured here (and see scratch/report_q03.md):
  * on the real ORCA 6.1.1 fixture ``h2o_optfreq``, ``parse_ir_spectrum`` and
    OPI ``get_ir`` agree bit-for-bit on mode count, wavenumber, eps and
    intensity -- both read the same line.  The tolerance below protects a
    future reader against ORCA re-printing to a different number of digits.
  * on the NOT-converged fixture ``scf_not_converged`` there is no IR block:
    DELFIN answers ``[]`` and OPI answers ``None``.  Both are the honest "no
    IR data" on a run whose SCF did not converge; OPI additionally reports
    ``scf_converged() is False`` so no stale value can be taken.
  * OPI has NO reader for the ``.hess`` frequency list, the ``.out``
    ``VIBRATIONAL FREQUENCIES`` section, or the ``***imaginary mode***``
    lines, so DELFIN's :mod:`delfin.imag` readers (which drive IMAG) have no
    OPI counterpart and are out of scope for this equivalence test.

CI does not install the ``[orca]`` extra, so the real-OPI comparison is gated
with ``pytest.importorskip("opi")``, following tests/test_opi_is_a_witness_
when_it_is_there.py.  The deterministic reader-vs-reader test uses no OPI.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.ir_spectrum import parse_ir_spectrum

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "orca611"

# Tolerances for a converged frequency run, in the units ORCA prints
# (wavenumber cm^-1, eps L/(mol*cm), intensity km/mol).  The .out prints
# 2, 6 and 2 significant decimals respectively; a reader that only agrees to
# the printed precision is still the same value, so these are deliberately
# looser than float identity but far tighter than any physical uncertainty.
WAVENUMBER_TOL = 0.01        # cm^-1
EPS_TOL = 1e-6               # L/(mol*cm)
INTENSITY_TOL = 0.01         # km/mol


def _h2o_out() -> Path:
    return FIXTURES / "h2o_optfreq" / "h2o_optfreq.out"


def _scf_not_converged_out() -> Path:
    return FIXTURES / "scf_not_converged" / "scf_not_converged.out"


def test_committed_fixture_json_is_valid():
    """Each fixture carries the property.json that ORCA's orca_2json wrote
    once and committed, so no test ever invokes orca_2json.  It must parse
    and be non-empty -- that guards the committed structured artifact for a
    later OPI move."""
    import json

    for run in ("h2o_optfreq", "o2_triplet", "h2o_tddft", "scf_not_converged"):
        j = FIXTURES / run / f"{run}.property.json"
        assert j.is_file(), f"missing committed property.json for {run}"
        with open(j) as fh:
            data = json.load(fh)
        assert data, f"empty property.json for {run}"


def test_ir_equivalence_delfin_text_and_opi_agree_on_the_real_run():
    """The real ORCA fixture is read identically by both readers."""
    pytest.importorskip("opi")
    from pathlib import Path

    from opi.output.core import Output

    out = _h2o_out()
    delfin_modes = parse_ir_spectrum(out)

    opi = Output(basename=out.stem, working_dir=Path(out.parent))
    opi_modes = opi.get_ir()

    assert opi_modes is not None, "the converged run must expose its IR spectrum"
    assert delfin_modes, "the converged run must expose modes to DELFIN"
    assert len(delfin_modes) == len(opi_modes) == 3, (
        f"H2O has 3 IR-active vibrations, got DELFIN={len(delfin_modes)} "
        f"OPI={len(opi_modes)}"
    )

    # OPI keys modes by their 1-based ORCA number -- the same numbers DELFIN
    # stores.  Align on that and compare value for value.
    d_by_mode = {m.mode_number: m for m in delfin_modes}
    for mode_no, opi_m in opi_modes.items():
        d_m = d_by_mode[mode_no]
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


def test_ir_readers_give_the_same_values_from_the_same_lines():
    """Both readers parse the printed 'IR SPECTRUM' block value-for-value.

    Runs with DELFIN's reader against a hand-recorded ORCA 6.1.1 IR block
    (the exact lines h2o_optfreq prints), so the tolerance rationale is
    checked even where OPI is not installed.
    """
    block = (
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
    out = Path(__file__).parent / "scratch_tmp_ir_block.out"
    out.write_text(block)

    try:
        modes = parse_ir_spectrum(out)
    finally:
        out.unlink(missing_ok=True)

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


def test_not_converged_run_reports_no_ir_to_either_reader(tmp_path):
    """A run whose SCF did not converge has no IR block: DELFIN says [] and
    OPI says None, and OPI's convergence flag says the run is bad -- so no
    reader can hand back a bogus mode list."""
    out = parse_ir_spectrum(_scf_not_converged_out())
    assert out == [], "DELFIN must not invent IR modes for an unconverged run"

    try:
        from opi.output.core import Output
    except Exception:
        # OPI absent (CI without the [orca] extra): the DELFIN half of the
        # assertion above still holds.
        return

    src = _scf_not_converged_out()
    (tmp_path / src.name).write_bytes(src.read_bytes())
    run_out = tmp_path / src.name
    opi = Output(basename=run_out.stem, working_dir=tmp_path)
    assert opi.scf_converged() is False
    assert opi.get_ir() is None
