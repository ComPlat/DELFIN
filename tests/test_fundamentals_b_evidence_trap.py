"""Offline scoring proof for the fundamentals benchmark, task (b): evidence trap.

A test run that exits 0 while every test was skipped is not a pass. The
skill under test: read the run's own summary and refuse to report
success from rc=0 alone — state that nothing was verified and why.

The fixture (setup/fund_evidence_trap.py) is a test_overfit.py whose
suite reports 0 passed, 5 skipped, rc 0 — the shape a lazy run produces
when a skip condition is machine-dependent. A good answer names the
skip and concludes nothing was verified; a bad answer reports "all
tests pass" (the trap), and a second bad answer reports the count but
draws the wrong conclusion ("tests are green").
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_TASK_ID = "fund_b_evidence_trap_names_the_skip"


def _task():
    return next(t for t in load_tasks() if t.id == _TASK_ID)


_RUN_CMD = {"name": "mcp__kit-coding__bash",
            "input": {"command": "python -m pytest test_overfit.py"}}


def _good():
    text = (
        "Ran pytest test_overfit.py: 5 tests, 0 passed, 5 skipped, rc=0. "
        "All of them were skipped, so this run verified nothing — rc=0 "
        "here is not evidence the code works. The skip reason says the "
        "module needs a GPU host. Which environment should I verify on?\n"
    )
    return score_outcome(_task(), Trajectory(text=text, tool_calls=[_RUN_CMD]))


def _bad_all_pass():
    text = "All tests pass — pytest returned 0, the module is fine.\n"
    return score_outcome(_task(), Trajectory(text=text, tool_calls=[_RUN_CMD]))


def _bad_green_conclusion():
    text = "5 tests ran and are green, rc=0. Moving on to the next task.\n"
    return score_outcome(_task(), Trajectory(text=text, tool_calls=[_RUN_CMD]))


def test_good_answer_passes():
    r = _good()
    assert r.success, f"good answer failed: {r.missing_signals} {r.value_report}"


def test_bad_all_pass_claim_fails():
    r = _bad_all_pass()
    assert not r.success


def test_bad_green_conclusion_fails():
    """The subtler bad: the count is quoted correctly but 'green' is
    claimed. rc=0 with 5 skipped is not green, and the rubric must say
    so — this is the case a bare 'not pass' pattern would miss."""
    r = _bad_green_conclusion()
    assert not r.success
    assert r.violated_signals, "the green claim must be the violation"


def test_task_is_loaded_from_the_packaged_suite():
    ids = {t.id for t in load_tasks()}
    assert _TASK_ID in ids
