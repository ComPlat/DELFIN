"""The committed ESD grading fixture: synthetic text shaped like job 8181368.

SYNTHETIC — not produced by ORCA.  Every ORCA-format string below is
hand-written and carries the marker; no ORCA-produced file is
committed (licence hold, operator 2026-10-31).  The shapes come from
the REAL practice run (SLURM job 8181368, formaldehyde, PBE0-D4/
def2-SVP gas, outputs under .gate/practice/): state order clean,
ISC rate problem-free and window-decayed, ISC rate 25x BELOW the
in-tree reference (FAIL — carried as FAIL with its source, never a
tuned PASS, per QS instruction), adiabatic Delta E_ST 0.889 eV
outside the experimental band (FAIL — the F6 regression case),
fluorescence and phosphorescence parsed with no problem.

The skip-if-absent test at the end reads the REAL run from
.gate/practice when present, so the real records are checked without
sharing them.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.common.esd_numerics import (
    est_reference_comparison,
    isc_reference_comparison,
    rate_problem,
    window_remainders,
)

# --- synthetic ORCA-format pieces (marked) ---------------------------------

ISC_OUT_SYNTH = """SYNTHETIC — not produced by ORCA. Practice-run shape, job 8181368.
Homogeneous linewidth is: 50 cm-1
Temperature used: 298.15 K
0-0 energy difference: 35301.0 cm-1
The calculated ISC rate constant is 9.634758E+05 s^-1
"""

STATE_LINES = {
    # FSPE values from the real run (Job 8181368); synthetic marker text
    # around them, no ORCA file committed.
    "S0": -114.284737350784,
    "S1": -114.145811100034,
    "T1": -114.178402699168,
}


def _state_out(state: str) -> str:
    return (
        f"SYNTHETIC — not produced by ORCA. {state} final single point.\n"
        f"Final Single Point Energy   {STATE_LINES[state]:.12f} Eh\n"
    )


def _delta_est_ev() -> float:
    # Convention of est_reference_comparison: Delta E_ST = E(S1) - E(T1),
    # TADF-favourable (T1 below S1) is POSITIVE.
    eh2ev = 27.211386245988
    return (STATE_LINES["S1"] - STATE_LINES["T1"]) * eh2ev


# --- the grading chain on synthetic records --------------------------------


def test_isc_rate_matches_the_run_shape_and_has_no_problem():
    assert rate_problem(ISC_OUT_SYNTH) is None
    # Window decayed: 50 cm-1 Lorentzian, no MAXTIME line -> ORCA's own
    # window, no remainder computable, i.e. no truncation evidence.
    assert window_remainders(ISC_OUT_SYNTH) == []


def test_isc_vs_reference_fails_25x_below():
    record = isc_reference_comparison(9.634758e5, None, True, "")
    assert record["verdict"] == "FAIL"
    # The value carries the ratio and the reference with its source.
    assert "0.04" in record["value"]
    assert "2.46e+07" in record["value"]


def test_delta_est_fails_the_band_the_loose_check_would_have_passed():
    record = est_reference_comparison(_delta_est_ev())
    assert record["verdict"] == "FAIL"
    assert abs(_delta_est_ev() - 0.889) < 0.01


def test_full_grade_records_carry_source_and_verdict():
    # One assembled record set as the grader emits it (shape contract).
    records = [
        isc_reference_comparison(9.634758e5, None, True, ""),
        est_reference_comparison(_delta_est_ev()),
    ]
    assert all(r["criterion"] and r["verdict"] in ("PASS", "FAIL", "NOT GRADED")
               for r in records)
    assert all(r["source"] for r in records)
    # No tuned PASS: both real-run verdicts are FAIL with sources.
    assert json.dumps(records).count('"FAIL"') == 2


# --- skip-if-absent: the REAL run from .gate/ -------------------------------

PRACTICE = Path(__file__).resolve().parents[1] / ".gate" / "practice"


@pytest.mark.skipif(
    not (PRACTICE / "ESD" / "S1_T1_ISC_ms0.out").exists(),
    reason="real practice run (job 8181368) not present in .gate/practice",
)
def test_real_run_isc_record_is_problem_free():
    from delfin.esd_results import collect_esd_results

    summary = collect_esd_results(
        PRACTICE / "ESD", ["S1", "T1"], ["S1>T1"], [], {}
    )
    (key, record), = summary.isc.items()
    assert record.rate is not None and record.problem is None
    comparison = isc_reference_comparison(
        record.rate, record.problem, True, ""
    )
    # The real run's ISC grades FAIL 25x below the reference — kept,
    # never tuned.
    assert comparison["verdict"] == "FAIL"
