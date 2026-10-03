"""Offline scoring proof for the fundamentals benchmark, task (e): waiting.

The skill under test: when a result arrives on a delay (a background
job), WAIT for it, collect it, and USE it — instead of guessing the
outcome, or collecting it without putting it into the answer.

The fixture (setup/fund_waiting.py) is samples.csv plus slow_report.py:
the script sleeps 20 s (a fixed, machine-independent constant — the
rubric never reads a wall clock) and then writes result.json with the
mean enthalpy of the six samples, 42.19 kJ/mol (mean computed by this
session and re-checked by the setup script). The agent is told to run
it in the background and report the mean.

Scoring shape (per the reviewer's F1, verified: tool RESULTS are in no
signal channel, as_string benchmark.py:396-406): the value only counts
when it reaches the PROSE — collected-but-unsaid fails. A bash
run/completion check must appear in the trace (tool_name signal), so
the good answer used the background machinery rather than a guess.

Three answers must separate:
* good: bash_background call in the trace + the value + unit in prose;
* guessing bad: the value invented from the samples file without any
  completion check — no background call in the trace;
* collected-but-unsaid bad: background call happened, but the answer
  never states the mean.
"""

from __future__ import annotations


from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_TASK_ID = "fund_e_background_result_is_collected"


def _task():
    return next(t for t in load_tasks() if t.id == _TASK_ID)


def _bg_calls():
    return [
        {"name": "mcp__kit-coding__bash_background",
         "input": {"command": "python fund_waiting/slow_report.py"}},
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "cat fund_waiting/result.json"}},
    ]


def _good():
    text = (
        "slow_report.py finished in the background; result.json reports "
        "the mean enthalpy of the six samples: 42.19 kJ/mol.\n"
    )
    return score_outcome(_task(), Trajectory(text=text, tool_calls=_bg_calls()))


def _bad_guess():
    text = (
        "Roughly averaging samples.csv myself: the mean enthalpy is "
        "about 42.19 kJ/mol.\n"
    )
    return score_outcome(_task(), Trajectory(text=text, tool_calls=[
        {"name": "mcp__kit-coding__read_file",
         "input": {"path": "fund_waiting/samples.csv"}}]))


def _bad_collected_but_unsaid():
    text = (
        "The background job is done and result.json is written. The "
        "report is ready for you.\n"
    )
    return score_outcome(_task(), Trajectory(text=text, tool_calls=_bg_calls()))


def test_good_answer_passes():
    r = _good()
    assert r.success, f"good answer failed: {r.missing_signals} {r.value_report}"


def test_bad_guess_fails():
    r = _bad_guess()
    assert not r.success


def test_bad_collected_but_unsaid_fails():
    """The F1 case: the value sat in a tool result the scorer cannot
    see, and the answer never restated it. Failing is the CORRECT
    verdict — the user was never told the number."""
    r = _bad_collected_but_unsaid()
    assert not r.success


def test_task_is_loaded_from_the_packaged_suite():
    ids = {t.id for t in load_tasks()}
    assert _TASK_ID in ids
