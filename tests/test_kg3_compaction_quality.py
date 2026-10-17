"""compaction_quality: what a compaction lost from a FIXED keep-list.

The keep-list (KG3 brief): the person's goal, their latest instruction,
constraints, refusals with their reasons, open decisions, waits, the last
test result. compare_loss (compaction_log.py) only covers four categories
and is log-only; this module extends the check. A loss is an item that is
in the history BEFORE the compaction but in neither the surviving state
block nor the summary after it.

Red control: against the previous commit the module does not exist, so
every case fails on import; after the commit every case passes.
"""

from __future__ import annotations

import pytest

from delfin.agent import compaction_quality as cq


def _machine(text: str) -> dict:
    return {"role": "user", "content": f"[Command results]\n{text}"}


def _state(items: list[str]) -> str:
    return "\n".join(f"  {i}" for i in items)


# -- the person's goal ------------------------------------------------------

def test_goal_not_lost_when_it_is_in_the_summary():
    before = [{"role": "user",
               "content": "Fix the failing test in compaction_log.py"}]
    lost = cq.quality_report(
        before, after_state_block="Goal: fix the failing test",
        after_summary_text="")
    assert lost["goal"]["lost"] == []


def test_goal_lost_when_nowhere_after():
    before = [{"role": "user",
               "content": "Fix the failing test in compaction_log.py"}]
    lost = cq.quality_report(
        before, after_state_block="", after_summary_text="other stuff")
    assert lost["goal"]["lost"] == [
        "Fix the failing test in compaction_log.py"]


def test_the_newest_of_several_goals_is_reported():
    before = [
        {"role": "user", "content": "Fix the failing test"},
        {"role": "user", "content": "Now also measure the cache"},
    ]
    lost = cq.quality_report(before, "", "")
    assert lost["goal"]["lost"] == ["Now also measure the cache"]


# -- the person's latest instruction ---------------------------------------

def test_latest_instruction_from_a_system_row_is_checked():
    before = [_machine("run"),
              {"role": "user",
               "content": "[System] Do not retry the same command"}]
    lost = cq.quality_report(before, "", "")
    assert any("Do not retry the same command" in line
               for line in lost["instruction"]["lost"])


def test_instruction_kept_verbatim_in_the_state_block_is_not_lost():
    before = [{"role": "user",
               "content": "[System] Do not retry the same command"}]
    lost = cq.quality_report(
        before, after_state_block=_state(
            ["Do not retry the same command"]), after_summary_text="")
    assert lost["instruction"]["lost"] == []


# -- constraints ------------------------------------------------------------

def test_a_constraint_the_user_stated_is_checked():
    before = [{"role": "user",
               "content": "Fix it. Constraint: never touch engine.py"}]
    lost = cq.quality_report(before, "", "")
    assert any("never touch engine.py" in line
               for line in lost["constraints"]["lost"])


def test_a_constraint_carried_in_the_block_is_not_lost():
    before = [{"role": "user",
               "content": "Fix it. Constraint: never touch engine.py"}]
    lost = cq.quality_report(
        before, after_state_block=_state(
            ["Constraint: never touch engine.py"]),
        after_summary_text="")
    assert lost["constraints"]["lost"] == []


# -- refusals with their reasons --------------------------------------------

def test_a_refusal_and_its_reason_are_checked():
    before = [_machine(
        "Refused: engine.py is read-only for you. Reason: safety frame")]
    lost = cq.quality_report(before, "", "")
    assert any("engine.py" in line and "safety frame" in line
               for line in lost["refusals"]["lost"])


def test_a_refusal_carried_in_the_block_is_not_lost():
    before = [_machine(
        "Refused: engine.py is read-only. Reason: safety frame")]
    lost = cq.quality_report(
        before, after_state_block=_state(
            ["engine.py refused: safety frame"]),
        after_summary_text="")
    assert lost["refusals"]["lost"] == []


# -- open decisions ---------------------------------------------------------

def test_an_open_decision_is_checked():
    before = [{"role": "assistant",
               "content": "Open decision: which of A or B to pick"}]
    lost = cq.quality_report(before, "", "")
    assert lost["decisions"]["lost"]


def test_an_open_decision_in_the_summary_is_not_lost():
    before = [{"role": "assistant",
               "content": "Open decision: which of A or B to pick"}]
    lost = cq.quality_report(before, "",
                             "Open decision: which of A or B to pick")
    assert lost["decisions"]["lost"] == []


# -- waits ------------------------------------------------------------------

def test_a_wait_is_checked():
    before = [_machine("watch_job 1234: waiting for the cluster job")]
    lost = cq.quality_report(before, "", "")
    assert any("waiting for" in line for line in lost["waits"]["lost"])


def test_a_wait_in_the_block_is_not_lost():
    before = [_machine("watch_job 1234: waiting for the cluster job")]
    lost = cq.quality_report(
        before, after_state_block=_state(
            ["waiting for the cluster job"]),
        after_summary_text="")
    assert lost["waits"]["lost"] == []


# -- the last test result ---------------------------------------------------

def test_the_last_test_result_is_checked():
    before = [_machine("gate -q tests/test_x.py -> 12 passed")]
    lost = cq.quality_report(before, "", "")
    assert any("12 passed" in line for line in lost["tests"]["lost"])


def test_the_last_test_result_in_the_block_is_not_lost():
    before = [_machine("gate -q tests/test_x.py -> 12 passed")]
    lost = cq.quality_report(
        before, after_state_block=_state(["12 passed"]),
        after_summary_text="")
    assert lost["tests"]["lost"] == []


# -- contract guards --------------------------------------------------------

def test_empty_history_reports_nothing():
    lost = cq.quality_report([], after_state_block="", after_summary_text="")
    assert all(not v["lost"] for v in lost.values())


def test_report_is_deterministic():
    before = [{"role": "user", "content": "Fix the failing test"}]
    a = cq.quality_report(before, "", "")
    b = cq.quality_report(before, "", "")
    assert a == b


def test_budget_tiny_still_compacts_is_a_compactor_contract():
    """Guardrail 3 of the reviewer: a keep-list may never refuse to
    compact — the cliff outranks keeping. The report here is advisory
    only (it must not feed back as a veto), asserted by its API: no
    function raises, no function returns a 'do not compact' verdict."""
    before = [{"role": "user", "content": "Fix the failing test"}]
    lost = cq.quality_report(before, "", "")
    assert set(lost) == {"goal", "instruction", "constraints", "refusals",
                         "decisions", "waits", "tests"}
