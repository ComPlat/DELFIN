"""A run's own files are graded with DELFIN's parsers, level-gated.

SYNTHETIC — not produced by ORCA. Small hand-written text in ORCA's
ESD output format drives the REAL product path — collect_esd_results
(parsing), isc_level_matches (level gate), isc_reference_comparison
(reference verdict) and window_remainders (truncation) — the exact
chain the practice-run grader (.gate/practice/grade_es1.py,
uncommitted per the operator's ruling) calls for job 8181368.

Old decision being replaced (rule 6, reviewer finding F5): the test
first pinned the .gate path of the grader script (9e419fa51), then an
un-landed accept module; it now exercises only committed product
code.
"""

from __future__ import annotations

from pathlib import Path

from delfin.common.esd_numerics import (
    isc_level_matches,
    isc_reference_comparison,
    rate_problem,
    window_remainders,
)
from delfin.esd_results import collect_esd_results

# Hand-written ORCA-format ESD output text, marked synthetic.
ISC_OUT = """SYNTHETIC — not produced by ORCA. Fixture for the grading script.
ISC calculation
Homogeneous linewidth is: 50 cm-1
Temperature used: 298.15 K
0-0 energy difference: 35301.0 cm-1
The calculated ISC rate constant is 2.5E+07 s^-1
"""

ISC_TRUNCATED_OUT = """SYNTHETIC — not produced by ORCA. Fixture for the grading script.
ISC calculation
Homogeneous linewidth is: 50 cm-1
Temperature used: 298.15 K
0-0 energy difference: 35301.0 cm-1
The calculated intersystem crossing rate constant is 2.5E+07 s^-1
"""

ISC_MISMATCH_OUT = """SYNTHETIC — not produced by ORCA. Fixture for the grading script.
ISC calculation
Homogeneous linewidth is: 120 cm-1
Temperature used: 298.15 K
0-0 energy difference: 35301.0 cm-1
The calculated intersystem crossing rate constant is 2.5E+07 s^-1
"""

MATCHED_CONTROL = {
    "functional": "PBE0",
    "disp_corr": "D4",
    "main_basisset": "def2-SVP",
    "solvent": "",
}


def _run_dir(tmp_path: Path, isc_text: str) -> Path:
    (tmp_path / "CONTROL.txt").write_text(
        "method=classic\nfunctional=PBE0\ndisp_corr=D4\n"
        "main_basisset=def2-SVP\ncharge=0\nmultiplicity=1\nsolvent=\n"
        "input_file=input.txt\nPAL=2\nmaxcore=2000\nESD_modul=yes\n"
        "ESD_modus=TDDFT\nESD_states=[S0,S1,T1]\nESD_ISC=S1>T1\nTROOTSSL=0\n"
    )
    esd = tmp_path / "ESD"
    esd.mkdir(exist_ok=True)
    for state, energy in (("S0", -114.15000), ("S1", -113.99000), ("T1", -114.04000)):
        (esd / f"{state}.out").write_text(
            "SYNTHETIC — not produced by ORCA.\n"
            f"TOTAL ENERGY      {energy:.8f} Eh\n"
            "The current calculation has used the line:\n"
        )
    (esd / "S1_T1_ISC_ms0.out").write_text(isc_text)
    return tmp_path


def _grade(isc_text: str, tmp_path: Path) -> dict:
    """The grader's chain on a synthetic run folder."""
    run = _run_dir(tmp_path, isc_text)
    summary = collect_esd_results(
        run / "ESD", ["S1", "T1"], ["S1>T1"], [], MATCHED_CONTROL
    )
    text = (run / "ESD" / "S1_T1_ISC_ms0.out").read_text()
    rec = next(iter(summary.isc.values()))
    ok, reason = isc_level_matches(MATCHED_CONTROL, text)
    return isc_reference_comparison(rec.rate, rate_problem(text), ok, reason)


def test_level_matched_synthetic_compares(tmp_path):
    record = _grade(ISC_OUT, tmp_path)
    assert record["verdict"] == "PASS", record
    assert "in-tree measurement" in record["source"]


def test_functional_mismatch_refuses(tmp_path):
    ok, reason = isc_level_matches(
        {**MATCHED_CONTROL, "functional": "B3LYP"}, ISC_OUT
    )
    assert ok is False
    assert "functional=B3LYP" in reason
    record = isc_reference_comparison(2.5e7, None, ok, reason)
    assert record["verdict"] == "NOT GRADED"
    assert "level differs" in record["value"]


def test_linewidth_mismatch_refuses(tmp_path):
    ok, reason = isc_level_matches(MATCHED_CONTROL, ISC_MISMATCH_OUT)
    assert ok is False
    assert "LINEW=120" in reason


def test_truncated_window_is_flagged_not_graded(tmp_path):
    # A window that ends while the correlation function still has
    # amplitude: window_remainders computes the remainder from ORCA's
    # "Maximum time:" line and the damping, rate_problem turns a
    # remainder above TRUNCATION_LIMIT into a NOT-A-RESULT, and the
    # comparison is then NOT GRADED — a truncated window never reads
    # as a result.
    text = ISC_OUT.replace(
        "The calculated ISC rate constant is 2.5E+07 s^-1",
        "Maximum time: 100 fs\n"
        "The calculated ISC rate constant is 2.5E+07 s^-1",
    )
    rems = window_remainders(text)
    assert rems and max(rems) > 1e-4, rems
    problem = rate_problem(text)
    assert problem is not None and "not converged" in problem
    record = isc_reference_comparison(2.5e7, problem, True, "")
    assert record["verdict"] == "NOT GRADED"
    assert "not converged" in record["value"]
