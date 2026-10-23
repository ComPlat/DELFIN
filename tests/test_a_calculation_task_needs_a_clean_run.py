"""Red control: a calculation task can complete with no evidence at all.

Wave-10 / s1, phase 2. check_completion_claim knows a test-task class
(_TEST_TASK_RE -> a green run in the window is required) but no
calculation-task class: a subject like "Berechne die Geometrie von
Wasser" or "Optimize the ligand and report the energy" completes as
plain 'unchecked' -- no run, no result-critic check, nothing. A
calculation task should need the same kind of proof a test task needs.

Kept as the permanent acceptance tests for the calc-task class; the
check only fires when a calc ledger exists, so code tasks ("Optimiere
den Import") stay untouched. Backward compatibility (no calcs argument
-> unchanged verdicts) is pinned by test_task_evidence_characterization.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent.task_evidence import check_completion_claim  # noqa: E402


_GREEN_CALCS = [
    {"ts": 100.0, "folder": "calc/water", "outcome": "succeeded",
     "worst": "ok"},
]
_FAILED_CALCS = [
    {"ts": 100.0, "folder": "calc/water", "outcome": "failed (exit code 1)",
     "worst": "error"},
]


def test_calc_task_with_a_clean_run_verifies():
    r = check_completion_claim(
        "Berechne die Geometrie von Wasser",
        changes=[], observed=[], window_start=0.0,
        calcs=_GREEN_CALCS)
    assert r["verdict"] == "verified", r
    assert r["kind"] == "calc", r


def test_calc_task_with_only_failed_runs_is_unmet():
    r = check_completion_claim(
        "Optimize the ligand and report the energy",
        changes=[], observed=[], window_start=0.0,
        calcs=_FAILED_CALCS)
    assert r["verdict"] == "unmet", r


def test_calc_task_without_a_calc_ledger_stays_unchecked():
    r = check_completion_claim(
        "Berechne die Geometrie von Wasser",
        changes=[], observed=[], window_start=0.0)
    assert r["verdict"] == "unchecked", r


def test_non_calc_subject_is_untouched():
    r = check_completion_claim(
        "Schreibe den Bericht als Markdown",
        changes=[{"path": "report.md", "ts": 50.0, "created": True}],
        observed=[], window_start=0.0)
    assert r["verdict"] == "verified", r
