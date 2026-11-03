"""Result checker for chemistry answers.

``delfin.tools.result_checks.check_run`` validates a finished ORCA run
against the failure modes known from autonomous chemistry agents:
values guessed instead of computed (a reported number with no source is
flagged as not computed), overconfident analysis (job-type consistency:
a minimum has 0 imaginary modes, a TS exactly 1 and its mode is named;
SCF and geometry converged; <S**2> deviation within a stated bound for
UKS; a T1 diagnostic for DLPNO-CCSD(T) reported and flagged above the
manual's threshold), and plausibility against reference values when the
task gives them.

This file runs under the CI interpreter (no OPI). The checker reads
from real fixture runs and from constructed stubs; values are taken
through :mod:`delfin.orca_read` and small text readers.
"""
from pathlib import Path

import pytest

from delfin.tools.result_checks import T1_THRESHOLD, check_run

# The real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2): real output is never committed. Every
# test that reads one SKIPS without the store; constructed stubs cover the
# bad cases without it, so this file still answers in plain CI.
from test_lc1_fixtures_helper import private_fixtures_or_skip  # noqa: E402


def _copy_fixture(name: str, tmp_path: Path) -> Path:
    """Copy a real-ORCA fixture run from the private store into tmp."""
    src = private_fixtures_or_skip(f"tests/fixtures/orca611/{name}")
    dst = tmp_path / name
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        (dst / f.name).write_bytes(f.read_bytes())
    return dst


def _stub(tmp_path: Path, name: str, out_body: str,
          inp_text: str = "#\n") -> Path:
    """Write a minimal ORCA-like run and return the work directory."""
    d = tmp_path / name
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{name}.out").write_text(out_body)
    (d / f"{name}.inp").write_text(inp_text)
    return d


# ---------------------------------------------------------------------------
# job-type consistency
# ---------------------------------------------------------------------------

def test_minimum_with_zero_imaginary_modes_passes(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(d, "h2o_optfreq", job_type="minimum",
                       reported={"e_tot": {"value": -76.321269021507,
                                           "source": str(d / "out")}})
    assert report["ok"] is True
    job = report["job_type"]
    assert job["imag_freq_count"] == 0
    assert job["consistent"] is True


def test_missing_scf_marker_is_reported(tmp_path):
    d = _stub(tmp_path, "bad",
              "ORCA TERMINATED NORMALLY\n"
              "FINAL SINGLE POINT ENERGY       -76.321269021507\n",
              inp_text="! B3LYP\n* xyz 0 1\nO 0 0 0\n*\n")
    report = check_run(d, "bad", job_type="sp")
    assert report["ok"] is False
    conv = report["convergence"]
    assert conv["scf_converged"] is None
    assert any(c["check"] == "scf_converged" and not c["passed"]
               for c in report["checks"])


def test_geometry_converged_enforced_for_opt_runs(tmp_path):
    d = _stub(tmp_path, "geom_bad",
              "SCF CONVERGED AFTER  8 ITERATIONS\n"
              "ORCA TERMINATED NORMALLY\n",
              inp_text="#\n")
    report = check_run(d, "geom_bad", job_type="minimum")
    conv = report["convergence"]
    assert conv["geometry_converged"] is not True
    assert not report["ok"]


def test_sp_run_geometry_untested_is_not_a_failure(tmp_path):
    d = _stub(tmp_path, "sp_ok",
              "SCF CONVERGED AFTER  8 ITERATIONS\n"
              "ORCA TERMINATED NORMALLY\n",
              inp_text="#\n")
    report = check_run(d, "sp_ok", job_type="sp")
    assert report["convergence"]["geometry_converged"] is None
    assert report["ok"] is True


def test_ts_with_exactly_one_named_mode_passes(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "\n"
        "                           --------------------\n"
        "                           VIBRATIONAL FREQUENCIES\n"
        "                           --------------------\n"
        "\n"
        "  4:   -123.45 cm**-1  ***imaginary mode***\n"
        "\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "ts_ok", body, inp_text="#\n")
    report = check_run(d, "ts_ok", job_type="ts")
    job = report["job_type"]
    assert job["consistent"] is True
    assert job["imag_freq_count"] == 1
    assert job["imag_modes"] == [(4, -123.45)]
    assert report["ok"] is True


def test_ts_with_zero_imaginary_modes_is_rejected(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "\n"
        "  1:    3800.00 cm**-1\n"
        "\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "ts_zero", body, inp_text="#\n")
    report = check_run(d, "ts_zero", job_type="ts")
    job = report["job_type"]
    assert job["consistent"] is False
    assert job["imag_freq_count"] == 0
    assert not report["ok"]
    assert any(c["check"] == "job_type" and not c["passed"]
               for c in report["checks"])


def test_ts_with_two_imaginary_modes_is_rejected(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "\n"
        "  4:   -120.00 cm**-1  ***imaginary mode***\n"
        "  5:    -90.00 cm**-1  ***imaginary mode***\n"
        "\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "ts_two", body, inp_text="#\n")
    report = check_run(d, "ts_two", job_type="ts")
    job = report["job_type"]
    assert job["consistent"] is False
    assert job["imag_freq_count"] == 2
    assert not report["ok"]


def test_minimum_with_an_imaginary_mode_is_rejected(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "\n"
        "  4:   -123.45 cm**-1  ***imaginary mode***\n"
        "\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "min_imag", body, inp_text="#\n")
    report = check_run(d, "min_imag", job_type="minimum")
    job = report["job_type"]
    assert job["consistent"] is False
    assert not report["ok"]


# ---------------------------------------------------------------------------
# spin: <S**2> deviation within a stated bound for UKS
# ---------------------------------------------------------------------------

def test_uks_s2_deviation_within_bound_passes(tmp_path):
    d = _copy_fixture("o2_triplet", tmp_path)
    report = check_run(d, "o2_triplet", job_type="sp",
                       s2_bound=0.05,
                       reported={"e_tot": {"value": -150.051820514846,
                                           "source": str(d / "out")}})
    spin = report["spin"]
    assert spin["multiplicity"] == 3
    assert spin["s2_deviation"] == pytest.approx(0.007206, abs=1e-5)
    assert spin["within_bound"] is True
    assert report["ok"] is True


def test_uks_s2_deviation_over_bound_fails(tmp_path):
    d = _copy_fixture("o2_triplet", tmp_path)
    report = check_run(d, "o2_triplet", job_type="sp", s2_bound=0.001)
    spin = report["spin"]
    assert spin["s2_deviation"] == pytest.approx(0.007206, abs=1e-5)
    assert spin["within_bound"] is False
    assert not report["ok"]


def test_uks_without_a_bound_is_reported_untested(tmp_path):
    d = _copy_fixture("o2_triplet", tmp_path)
    report = check_run(d, "o2_triplet", job_type="sp")
    spin = report["spin"]
    assert spin["s2_deviation"] == pytest.approx(0.007206, abs=1e-5)
    assert spin["s2_bound"] is None
    assert spin["within_bound"] is None
    assert any(c["check"] == "s2_bound" and c["passed"] is False
               and c["untested"] is True
               for c in report["checks"])


# ---------------------------------------------------------------------------
# T1 diagnostic for DLPNO-CCSD(T)
# ---------------------------------------------------------------------------

def test_t1_missing_on_dlpno_run_is_reported(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "dlpno_no_t1", body,
              inp_text="#\n")
    report = check_run(d, "dlpno_no_t1", job_type="sp", method="DLPNO-CCSD(T)")
    t1 = report["t1_diagnostic"]
    assert t1["present"] is False
    assert not report["ok"]
    assert any(c["check"] == "t1_diagnostic" and not c["passed"]
               for c in report["checks"])


def test_t1_below_threshold_is_reported_not_flagged(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "T1 diagnostic                              ...      0.006287935\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "dlpno_ok", body, inp_text="#\n")
    report = check_run(d, "dlpno_ok", job_type="sp", method="DLPNO-CCSD(T)")
    t1 = report["t1_diagnostic"]
    assert t1["present"] is True
    assert t1["value"] == pytest.approx(0.006287935, abs=1e-9)
    assert t1["threshold"] == T1_THRESHOLD
    assert t1["flagged"] is False
    assert report["ok"] is True


def test_t1_dlpno_variant_line_is_read(tmp_path):
    # ORCA also prints the diagnostic as "T1(DLPNO) : 0.0104"; the
    # reader must recognise both spellings, not only "T1 diagnostic ...".
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "T1(DLPNO) : 0.0104\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "dlpno_variant", body, inp_text="#\n")
    report = check_run(d, "dlpno_variant", job_type="sp", method="DLPNO-CCSD(T)")
    t1 = report["t1_diagnostic"]
    assert t1["present"] is True
    assert t1["value"] == pytest.approx(0.0104, abs=1e-9)
    assert t1["flagged"] is False
    assert report["ok"] is True


def test_committed_dlpno_synth_fixture_t1_read_and_scf_flagged(tmp_path):
    # The committed synthetic DLPNO-CCSD(T1) fixture (orca611, "CONSTRUCTED
    # fixture" header): T1(DLPNO) : 0.0104, no "SCF CONVERGED AFTER" line
    # but "ORCA SCF ... Total Energy" converged output, terminated normally.
    # The checker reads the T1 value and does NOT wave the unstated SCF
    # through: a normally-terminated run must state its SCF outcome, so
    # the run is reported not-ok with the scf_converged check failing.
    src = (Path(__file__).parent / "fixtures" / "orca611"
           / "butane_dlpno_ccsdt_pvtz")
    d = tmp_path / "butane_dlpno_ccsdt_pvtz"
    d.mkdir(parents=True)
    for f in src.iterdir():
        (d / f.name).write_bytes(f.read_bytes())
    report = check_run(d, "butane_dlpno_ccsdt_pvtz", job_type="sp",
                       method="DLPNO-CCSD(T)")
    t1 = report["t1_diagnostic"]
    assert t1["present"] is True
    assert t1["value"] == pytest.approx(0.0104, abs=1e-9)
    assert t1["flagged"] is False
    assert report["ok"] is False
    assert any(c["check"] == "scf_converged" and not c["passed"]
               for c in report["checks"])


def test_t1_above_threshold_is_flagged(tmp_path):
    body = (
        "SCF CONVERGED AFTER  8 ITERATIONS\n"
        "T1 diagnostic                              ...      0.045100000\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    d = _stub(tmp_path, "dlpno_hi", body, inp_text="#\n")
    report = check_run(d, "dlpno_hi", job_type="sp", method="DLPNO-CCSD(T)")
    t1 = report["t1_diagnostic"]
    assert t1["present"] is True
    assert t1["flagged"] is True
    # A run with a high T1 is reported and flagged, not silently passed.
    assert any(c["check"] == "t1_diagnostic" and c["flagged"] is True
               for c in report["checks"])


def test_t1_not_required_for_non_dlpno(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(d, "h2o_optfreq", job_type="minimum", method="B3LYP")
    assert report["t1_diagnostic"]["required"] is False
    assert "present" not in report["t1_diagnostic"] or report["ok"]


# ---------------------------------------------------------------------------
# traceability: a value with no source is flagged as not computed
# ---------------------------------------------------------------------------

def test_value_without_source_is_not_computed(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(
        d, "h2o_optfreq", job_type="minimum",
        reported={"e_tot": {"value": -76.321269021507, "source": None}},
    )
    trace = report["traceability"]
    assert trace["e_tot"]["not_computed"] is True
    assert not report["ok"]
    assert any(c["check"] == "traceability" and not c["passed"]
               for c in report["checks"])


def test_value_with_source_is_computed(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(
        d, "h2o_optfreq", job_type="minimum",
        reported={"e_tot": {"value": -76.321269021507,
                            "source": str(d / "out")}},
    )
    assert report["traceability"]["e_tot"]["not_computed"] is False
    assert report["ok"] is True


def test_reported_missing_totally_is_not_computed(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(d, "h2o_optfreq", job_type="minimum",
                       reported={"gibbs": {"value": -76.32}})
    assert report["traceability"]["gibbs"]["not_computed"] is True


# ---------------------------------------------------------------------------
# plausibility against reference values
# ---------------------------------------------------------------------------

def test_value_within_reference_tolerance_passes(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(
        d, "h2o_optfreq", job_type="minimum",
        reported={"e_tot": {"value": -76.321269021507, "source": str(d / "out")}},
        references={"e_tot": {"expected": -76.321269021507, "tolerance": 1e-6}},
    )
    plaus = report["plausibility"]["e_tot"]
    assert plaus["within_tolerance"] is True
    assert report["ok"] is True


def test_value_outside_reference_tolerance_fails(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(
        d, "h2o_optfreq", job_type="minimum",
        reported={"e_tot": {"value": -76.0, "source": str(d / "out")}},
        references={"e_tot": {"expected": -76.321269021507, "tolerance": 1e-6}},
    )
    plaus = report["plausibility"]["e_tot"]
    assert plaus["within_tolerance"] is False
    assert not report["ok"]


def test_reference_with_no_reported_value_fails(tmp_path):
    d = _copy_fixture("h2o_optfreq", tmp_path)
    report = check_run(
        d, "h2o_optfreq", job_type="minimum",
        reported={},
        references={"e_tot": {"expected": -76.321269021507, "tolerance": 1e-6}},
    )
    plaus = report["plausibility"]["e_tot"]
    assert plaus["within_tolerance"] is False
    assert not report["ok"]
