"""An ESD rate constant that is not a result says so, and the time grid is ORCA's.

Measured on ORCA 6.1.1, formaldehyde, the job DELFIN wrote (MAXTIME 12000
a.u., 290 fs) against ORCA's own window (2934 fs for LINEW 50 cm-1):

    ISC S1>T1 (Ms=0)   5.698e7 s-1   ->  2.460e7 s-1  (4x and 8x the old window: 2.4606e7, 2.4598e7)
    IC  S1>S0         -8.064e-3 s-1  ->  3.795e-3 s-1 (ORCA: "negative rates are unphysical")

At 290 fs the correlation function still had 6.5 % of its amplitude
(exp(-LINEW t) with LINEW as a Hartree frequency).  More points on the same
window did not help (-7.82e-3); a longer window did.  Every archived CONTROL
file carries the old pair NPOINTS 131072 / MAXTIME 12000, and every archived
ISC of a TADF emitter (Emitter8) was computed on that window.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.common.esd_numerics import TRUNCATION_LIMIT, rate_problem, window_remainders
from delfin.esd_input_generator import _esd_grid_lines
from delfin.esd_results import ICResult, ISCResult, _parse_ic_output, _parse_isc_output, esd_rate_problem


def _esd_output(maximum_time_fs: float, rate: str, *, homogeneous: float = 50.0, inhomogeneous=None,
                negative_warning: bool = False, kk_error: bool = False, kind: str = "internal conversion") -> str:
    lines = ["------------------------------------------------------------------------------",
             "                           ORCA EXCITED STATE DYNAMICS",
             "------------------------------------------------------------------------------"]
    if kk_error:
        lines.append("Error (ORCA_ESD): The K*K value is too large (>30)!")
    lines.append(f"Homogeneous linewidth is:\t\t\t{homogeneous:.2f} cm-1")
    if inhomogeneous is not None:
        lines.append(f"Inhomogeneous linewidth is:\t\t\t{inhomogeneous:.2f} cm-1")
    lines += ["Number of points:\t\t\t\t131072",
              f"Maximum time:\t\t\t\t\t{maximum_time_fs:.2f} fs",
              "Temperature used:\t\t\t\t298.15 K",
              "0-0 energy difference:\t\t\t\t30296.80 cm-1",
              f"The calculated {kind} rate constant is\t{rate} s-1"]
    if negative_warning:
        lines.append("WARNING: negative rates are unphysical! It means something went wrong with the CorrFunc integration.")
    lines.append("                             ****ORCA TERMINATED NORMALLY****")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- the check

def test_orcas_own_window_is_converged_and_the_old_one_is_not():
    # ORCA prints "Cutoff for the correlation function: 1.00e-12"; its window reaches it
    assert window_remainders(_esd_output(2933.77, "3.79e-03"))[0] == pytest.approx(1.0e-12, rel=0.05)
    assert window_remainders(_esd_output(290.27, "5.70e+07"))[0] == pytest.approx(0.065, rel=0.02)
    assert window_remainders(_esd_output(1160.3, "2.46e+07"))[0] < TRUNCATION_LIMIT

    assert rate_problem(_esd_output(2933.77, "3.79e-03")) is None
    assert "6.5%" in rate_problem(_esd_output(290.27, "5.70e+07", kind="ISC"))


def test_a_negative_rate_is_not_a_result_even_on_a_long_window():
    problem = rate_problem(_esd_output(2933.77, "-8.06e-03", negative_warning=True))

    assert problem.startswith("unphysical")


def test_a_voigt_profile_decays_by_both_widths():
    # the Gaussian part alone takes the remainder far below the limit
    assert window_remainders(_esd_output(290.27, "1.0e+05", inhomogeneous=250.0))[0] < 1e-30


def test_the_parsed_result_carries_the_problem(tmp_path):
    ic = tmp_path / "S1_S0_IC.out"
    ic.write_text(_esd_output(290.27, "-8.063829e-03", negative_warning=True))
    isc = tmp_path / "S1_T1_ISC_ms0.out"
    isc.write_text(_esd_output(290.27, "5.698207e+07", kind="ISC", kk_error=True))

    assert _parse_ic_output(ic).problem.startswith("unphysical")
    parsed = _parse_isc_output(isc)
    assert "not converged" in parsed.problem and "K*K" in parsed.problem
    assert esd_rate_problem(_esd_output(2933.77, "3.79e-03")) is None


def test_the_reports_mark_what_is_not_a_result(tmp_path):
    from delfin.esd_results import ESDSummary
    from delfin.reporting.esd_report import generate_esd_report

    ok = ICResult(3.79e-3, 298.15, 30296.8, tmp_path / "a.out")
    bad = ICResult(-8.06e-3, 298.15, 30296.8, tmp_path / "b.out", problem="unphysical: negative")
    summary = ESDSummary(states={}, isc={}, ic={"S1>S0": ok, "S2>S0": bad})
    out = tmp_path / "ESD.txt"
    generate_esd_report(summary, out)
    text = out.read_text()

    assert "S1>S0: 3.790000e-03 s^-1" in text and "S1>S0: 3.790000e-03 s^-1 (T=298.15 K, Δ0-0=30296.80 cm^-1)\n" in text
    assert "S2>S0: -8.060000e-03 s^-1 (T=298.15 K, Δ0-0=30296.80 cm^-1)  [NOT A RESULT: unphysical: negative]" in text


def test_older_results_keep_their_fields():
    # positional construction as the existing callers do it still works
    assert ISCResult(1.0, 298.15, 1.0, None, None, None, Path("x")).problem is None


# ----------------------------------------------------------------- the grid

@pytest.mark.parametrize("config, lines", [
    ({}, []),
    ({"ESD_NPOINTS": "auto", "ESD_MAXTIME": "auto"}, []),
    ({"ESD_NPOINTS": 131072, "ESD_MAXTIME": 12000}, []),          # every archived CONTROL file
    ({"ESD_NPOINTS": 262144, "ESD_MAXTIME": 60000}, ["  NPOINTS         262144", "  MAXTIME         60000"]),
    ({"ESD_NPOINTS": "auto", "ESD_MAXTIME": 80000}, ["  MAXTIME         80000"]),
])
def test_the_time_grid_is_left_to_orca_unless_chosen(config, lines):
    assert _esd_grid_lines(config) == lines


def test_isc_and_ic_inputs_carry_no_grid_from_the_template(tmp_path):
    from delfin.config import _load_template_defaults
    from delfin.esd_input_generator import create_ic_input, create_isc_input

    config = {**_load_template_defaults(), "functional": "PBE0", "main_basisset": "def2-SVP", "PAL": 4,
              "maxcore": 1000, "solvent": "water", "ESD_modul": "yes"}
    for state in ("S0", "S1", "T1"):
        (tmp_path / f"{state}.xyz").write_text("2\n\nH 0 0 0\nH 0 0 0.74\n")
    isc = Path(create_isc_input("S1>T1", tmp_path, 0, "water", [], "def2-SVP", "def2-TZVP", config, trootssl=0))
    ic = Path(create_ic_input("S1>S0", tmp_path, 0, "water", [], "def2-SVP", "def2-TZVP", config))

    for text in (isc.read_text(), ic.read_text()):
        assert "%ESD" in text
        assert "NPOINTS" not in text and "MAXTIME" not in text


# -------------------------------------------------- ISC from the sublevels

def test_the_observed_isc_rate_sums_from_a_singlet_and_averages_from_a_triplet():
    """ORCA manual 6.1.1, 5.5.4: singlet-to-triplet the sublevel rates add up;
    triplet-to-singlet 'the observed rate constant is the average, not the
    sum' -- each sublevel holds a third of the triplet.  The reports summed
    both, so every RISC rate (the TADF figure of merit) came out 3x too fast.
    Formaldehyde, dashboard run: T1>S1 sublevels 3.74e6, 3.74e6, 7.40e7."""
    from delfin.esd_results import observed_isc_rate

    sub = {1: 3.744717e6, -1: 3.744717e6, 0: 7.404665e7}
    assert observed_isc_rate("S1", sub) == pytest.approx(8.153608e7)
    assert observed_isc_rate("T1", sub) == pytest.approx(8.153608e7 / 3)
    assert observed_isc_rate("T1", {0: 3.0e7}) == pytest.approx(3.0e7)        # one computed: all equal
    assert observed_isc_rate("S1", {0: 1.0e7, 1: 2.0e6}) == pytest.approx(1.4e7)  # Ms=-1 mirrors Ms=+1
    assert observed_isc_rate("T1", {}) is None


def test_every_report_gives_the_same_observed_rate(tmp_path):
    from delfin.esd_results import ESDSummary
    from delfin.reporting.delfin_collector import _calculate_total_isc_rate, parse_esd_summary
    from delfin.reporting.esd_report import generate_esd_report

    def isc(rate):
        return ISCResult(rate, 298.15, 324.24, None, None, None, tmp_path / "x.out")

    summary = ESDSummary(isc={"T1>S1(Ms=+1)": isc(3.0e6), "T1>S1(Ms=-1)": isc(3.0e6), "T1>S1(Ms=0)": isc(6.0e7)})
    generate_esd_report(summary, tmp_path / "ESD.txt")
    parsed = parse_esd_summary(tmp_path)["isc"]["T1_S1"]

    assert parsed["total_rate_s1"] == pytest.approx(2.2e7, rel=1e-6)
    assert _calculate_total_isc_rate(parsed["ms_components"], "T1") == pytest.approx(2.2e7)
    assert _calculate_total_isc_rate(parsed["ms_components"], "S1") == pytest.approx(6.6e7)


def test_a_report_rebuilt_for_an_older_run_corrects_and_flags_it(tmp_path):
    """An ESD.txt written before 2026-09 summed T1>S1 and carried no flags;
    the JSON (and DELFIN.docx) takes the rate from the sublevels and the
    flags from the ORCA outputs themselves."""
    from delfin.reporting.delfin_collector import collect_esd_data

    (tmp_path / "CONTROL.txt").write_text("ESD_modul=yes\nstates=S1,T1\nISCs=T1>S1\n")
    (tmp_path / "ESD.txt").write_text(
        "ISC rate constants (s^-1):\n"
        "  T1>S1(Ms=+1): 3.000000e+06 s^-1 (T=298.15 K, Δ0-0=324.24 cm^-1)\n"
        "  T1>S1(Ms=-1): 3.000000e+06 s^-1 (T=298.15 K, Δ0-0=324.24 cm^-1)\n"
        "  T1>S1(Ms=0): 6.000000e+07 s^-1 (T=298.15 K, Δ0-0=324.24 cm^-1)\n"
        "  T1>S1 (total): 6.600000e+07 s^-1\n")
    esd = tmp_path / "ESD"
    esd.mkdir()
    for ms, rate in (("ms0", "6.0e+07"), ("msp1", "3.0e+06"), ("msm1", "3.0e+06")):
        (esd / f"T1_S1_ISC_{ms}.out").write_text(_esd_output(290.27, rate, kind="ISC"))

    entry = collect_esd_data(tmp_path)["intersystem_crossing"]["T1_S1"]

    assert entry["total_rate_s1"] == pytest.approx(2.2e7)
    assert "not converged" in entry["problem"]
