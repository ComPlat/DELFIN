"""A run that executed no tests is not evidence for "tests pass".

Welle 11, package B, phase 3. The green-run rule in
``check_completion_claim`` (delfin/agent/task_evidence.py, the tests
branch) required only ``failed == 0`` and a benign status. A pytest run
that collected nothing, or one where every test was skipped or
deselected, produces exactly that: ``status="ok"``, ``passed=0``,
``failed=0`` -- pytest exits 0 and the runner maps that to "ok"
(delfin/agent/test_runner.py:353-358). The bash-side red-state clearing
in api_client.py already states the principle ("Only a run that
demonstrably executed tests ('N passed') clears red state"); the
completion check now applies the same principle to its ledger.
"""

from __future__ import annotations

from delfin.agent.task_evidence import check_completion_claim


SUBJECT = "Run the test suite and verify everything passes"


def _check(tests, **kw):
    return check_completion_claim(
        SUBJECT, changes=[], observed=[], tests=tests, **kw)


def test_a_skipped_only_run_is_not_evidence():
    run = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
           "status": "ok", "passed": 0, "failed": 0, "skipped": 5,
           "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert check["verdict"] == "unmet", check


def test_a_zero_collected_run_is_not_evidence():
    run = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
           "status": "ok", "passed": 0, "failed": 0, "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert check["verdict"] == "unmet", check


def test_the_note_says_the_run_executed_no_tests():
    run = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
           "status": "ok", "passed": 0, "failed": 0, "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert "no test" in check["note"] or "executed" in check["note"]


def test_a_run_with_passes_is_still_evidence():
    """The control in the other direction: a run that actually executed
    tests stays green."""
    run = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
           "status": "ok", "passed": 3, "failed": 0, "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert check["verdict"] == "verified", check


def test_a_partially_skipped_run_is_still_evidence():
    """Skips beside real executions do not invalidate the run: pytest
    reports "3 passed, 2 skipped" as a green run and so does this
    check."""
    run = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
           "status": "ok", "passed": 3, "failed": 0, "skipped": 2,
           "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert check["verdict"] == "verified", check


def test_a_bash_run_with_zero_executed_tests_is_not_evidence():
    """The bash path parses counts from console output; a "no tests ran"
    summary line must not read as green either."""
    run = {"tool": "bash", "command": "pytest tests/", "exit_code": 0,
           "status": "ok", "passed": 0, "failed": 0, "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert check["verdict"] == "unmet", check


def test_a_legacy_entry_without_count_keys_still_verifies():
    """The pre-count ledger shape carries only failed/status. It does
    not CLAIM that nothing ran -- only a run that positively says
    passed == 0 proves an empty execution."""
    run = {"tool": "run_tests", "failed": 0, "status": "ok", "ts": 100.0}
    check = _check([run], window_start=0.0)
    assert check["verdict"] == "verified", check


def test_mixed_runs_green_plus_empty_are_evidence():
    """One run that executed tests and passed carries the claim; the
    empty run beside it does not veto it."""
    empty = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
             "status": "ok", "passed": 0, "failed": 0, "ts": 90.0}
    good = {"tool": "run_tests", "command": "tests/", "exit_code": 0,
            "status": "ok", "passed": 4, "failed": 0, "ts": 100.0}
    check = _check([empty, good], window_start=0.0)
    assert check["verdict"] == "verified", check
