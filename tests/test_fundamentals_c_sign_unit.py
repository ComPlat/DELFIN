"""Offline scoring proof for the fundamentals benchmark, task (c): sign/unit.

The skill under test: carry a sign through a unit conversion and say
what it means. The fixture gives two total energies in hartree, E(A) =
-154.772300 Eh and E(B) = -154.793100 Eh, and the conversion factor
(1 Eh = 627.5094740631 kcal/mol). The defined difference is

    dE = E(B) - E(A) = -0.0208 Eh = -13.052 kcal/mol,

negative because B is LOWER — B is the more stable isomer by
13.052 kcal/mol. Every number here was computed by this session
(python: EB-EA = -0.0208; -0.0208 * 627.5094740631 = -13.052) and is
re-checked by the setup script, so the rubric grades measured values.

Three answers must separate:
* good: the signed difference in kcal/mol AND "B is more stable";
* sign-flip bad: +13.05 kcal/mol, "A is more stable" — the magnitude
  matches but the interpretation signal is missing;
* unit bad: "-13.05 kJ/mol" (the kJ figure would be -54.61, so this
  wrong-unit answer states the right NUMBER with the wrong unit) —
  caught by the scoped forbidden pattern, because the figure and the
  interpretation both look right.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome

_TASK_ID = "fund_c_sign_and_unit_of_a_difference"


def _task():
    return next(t for t in load_tasks() if t.id == _TASK_ID)


def _good():
    text = (
        "dE = E(B) - E(A) = -154.793100 - (-154.772300) Eh = -0.0208 Eh. "
        "With 1 Eh = 627.509 kcal/mol that is dE = -13.052 kcal/mol "
        "(=-54.61 kJ/mol). The difference is negative, so B has the "
        "lower energy: B is the more stable isomer, by 13.052 kcal/mol.\n"
    )
    calls = [{"name": "mcp__kit-coding__read_file",
              "input": {"path": "fund_sign_unit/energies.txt"}}]
    return score_outcome(_task(), Trajectory(text=text, tool_calls=calls))


def _bad_sign_flip():
    text = (
        "dE = E(B) - E(A) = +13.052 kcal/mol, so A is the more stable "
        "isomer by 13.052 kcal/mol.\n"
    )
    return score_outcome(
        _task(), Trajectory(text=text, tool_calls=[
            {"name": "mcp__kit-coding__read_file",
             "input": {"path": "fund_sign_unit/energies.txt"}}]))


def _bad_wrong_unit():
    text = (
        "The difference is -13.05 kJ/mol; B is the more stable isomer.\n"
    )
    return score_outcome(
        _task(), Trajectory(text=text, tool_calls=[
            {"name": "mcp__kit-coding__read_file",
             "input": {"path": "fund_sign_unit/energies.txt"}}]))


def test_good_answer_passes():
    r = _good()
    assert r.success, f"good answer failed: {r.missing_signals} {r.value_report}"


def test_bad_sign_flip_fails():
    r = _bad_sign_flip()
    assert not r.success
    # The magnitude is right; what must fail is the direction — either
    # the interpretation signal (B more stable) or the value verdict.
    assert r.missing_signals or any(
        v.endswith(":wrong") or v.endswith(":absent")
        for v in r.value_report.values()), r.value_report


def test_bad_wrong_unit_fails():
    r = _bad_wrong_unit()
    assert not r.success
    assert r.violated_signals, "the kJ-as-answer claim must be the violation"


def test_task_is_loaded_from_the_packaged_suite():
    ids = {t.id for t in load_tasks()}
    assert _TASK_ID in ids
