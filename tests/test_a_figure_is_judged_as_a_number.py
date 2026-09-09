"""A rubric that reads digits cannot read a number.

Every figure this suite checks was a regex over characters: `6\\.?070[,.]55`
for a total, `23[.,]?716` for the one it must not be. Three things that
cannot do:

* it has to be written for the separators the answer happens to use.
  1.234,50 and 1,234.50 and 1234.5 are one amount, and the pattern sees
  three unrelated strings — so a correct answer fails for writing its
  figure the way its document writes it, which is the error this suite
  has actually made, repeatedly;
* it cannot express a tolerance, so a rubric either demands a rounding
  the model was never told to use, or is loosened until it stops
  discriminating;
* it cannot say WHY it missed. A wrong total and an answer that never
  arrived produce the same red line, and they need opposite repairs.

``expected_values`` asks the question directly: this figure, within this
tolerance, in whatever convention it was written. The verdict is
matched / wrong / absent, and the run summary counts the three apart —
because an absent figure is usually not arithmetic, it is an empty turn
or a truncation, and counting it as a model error is how a harness fault
gets filed against the model.

The permissiveness is deliberate and bounded: an ambiguous token like
"1.265" is read BOTH ways, since the only cost is failing to catch a
wrong answer, while the stricter reading fails correct ones.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark import (
    ExpectedValue, Task, Trajectory, numbers_in, readings_of, score_outcome,
    summarise_run,
)


# ---------------------------------------------------------------------------
# Reading a written figure
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("written,expected", [
    ("6.070,55", 6070.55),      # German
    ("6,070.55", 6070.55),      # English
    ("6070.55", 6070.55),       # plain
    ("23.716,25", 23716.25),
    ("1.234.567,89", 1234567.89),
    ("1265", 1265.0),
    ("0,00", 0.0),
    ("-12,00", -12.0),
])
def test_one_amount_in_every_convention_reads_as_that_amount(written, expected):
    assert expected in readings_of(written), (written, readings_of(written))


def test_an_ambiguous_token_is_read_both_ways():
    """Nothing in '1.265' says whether it is a thousand or one and a bit."""
    assert set(readings_of("1.265")) == {1.265, 1265.0}


def test_a_separator_that_cannot_be_grouping_is_not_read_as_grouping():
    """'1,5' is one and a half. Fifteen would need '1,500'."""
    assert readings_of("1,5") == (1.5,)
    assert 15.0 not in readings_of("1,5")


def test_the_figures_are_found_inside_a_sentence():
    found = numbers_in("Im März 2026 wurden 6.070,55 € auf 5 von 25 Zeilen "
                       "gebucht.")
    for value in (2026.0, 6070.55, 5.0, 25.0):
        assert value in found, (value, found)


def test_a_record_reference_is_not_a_figure():
    """Administrative answers are full of them, and without this an
    answer that merely NAMES beleg R-014 satisfies an expected 14."""
    assert numbers_in("Beleg R-014 und Sheet2") == ()
    assert numbers_in("Vorgang A-2026-003") == ()
    # A date is a date, not three figures.
    assert numbers_in("am 2026-03-19") == (2026.0,)


def test_a_sign_still_works_where_a_sign_can_stand():
    for text in ("-12,00 EUR", "Saldo: -12,00", "(-3,5)"):
        assert any(v < 0 for v in numbers_in(text)), text


# ---------------------------------------------------------------------------
# The verdict, and why it is three-valued
# ---------------------------------------------------------------------------

def _judge(text, value=6070.55, **kw):
    return ExpectedValue(value=value, **kw).judge(text)


def test_the_right_figure_matches_however_it_is_written():
    for text in ("Im März: 6.070,55 €.", "March total 6,070.55 EUR",
                 "6070.55", "Summe: 6.070,55"):
        assert _judge(text) == "matched", text


def test_a_different_figure_is_wrong_not_absent():
    assert _judge("Im März: 23.716,25 €.") == "wrong"


def test_no_figure_at_all_is_absent_not_wrong():
    """The distinction the summary is built on."""
    assert _judge("Das konnte ich nicht ermitteln.") == "absent"
    assert _judge("[empty turn] The backend ended this turn without any "
                  "answer text.") in ("absent", "wrong")


def test_the_tolerance_is_relative_and_has_an_absolute_floor():
    assert _judge("6070.60", tolerance=0.001) == "matched"
    assert _judge("6100.00", tolerance=0.001) == "wrong"
    # A figure of zero needs the absolute floor; a relative margin on 0 is 0.
    assert ExpectedValue(value=0.0, absolute=0.01).judge("0,004") == "matched"
    assert ExpectedValue(value=0.0, absolute=0.01).judge("0,5") == "wrong"


# ---------------------------------------------------------------------------
# Scoring and the run summary
# ---------------------------------------------------------------------------

def _task(**kw):
    return Task(id="probe", task_class="office", mode="office", prompt="p",
                expected_values=(ExpectedValue(value=6070.55, label="maerz",
                                               **kw),))


def test_a_required_figure_decides_the_task():
    assert score_outcome(_task(), Trajectory(text="6.070,55 €")).success
    assert not score_outcome(_task(), Trajectory(text="23.716,25 €")).success


def test_an_optional_figure_does_not():
    assert score_outcome(_task(optional=True),
                         Trajectory(text="keine Zahl hier")).success


def test_the_result_records_why_it_missed():
    wrong = score_outcome(_task(), Trajectory(text="23.716,25 €"))
    absent = score_outcome(_task(), Trajectory(text="weiß ich nicht"))
    assert list(wrong.value_report.values()) == ["wrong"]
    assert list(absent.value_report.values()) == ["absent"]
    # And the label travels, so a run with several figures is readable.
    assert all("maerz" in k for k in wrong.value_report)


def test_the_summary_counts_the_two_failures_apart():
    rows = [
        score_outcome(_task(), Trajectory(text="6.070,55 €")).__dict__,
        score_outcome(_task(), Trajectory(text="23.716,25 €")).__dict__,
        score_outcome(_task(), Trajectory(text="weiß ich nicht")).__dict__,
    ]
    s = summarise_run(rows)
    assert s["values_matched"] == 1
    assert s["values_wrong"] == 1
    assert s["values_absent"] == 1


def test_a_run_with_no_figures_reports_zeroes_not_a_crash():
    plain = Task(id="p", task_class="office", mode="office", prompt="p")
    s = summarise_run([score_outcome(plain, Trajectory(text="ok")).__dict__])
    assert s["values_matched"] == s["values_wrong"] == s["values_absent"] == 0


# ---------------------------------------------------------------------------
# The task that now uses it
# ---------------------------------------------------------------------------

def test_the_month_total_is_asked_for_as_a_value():
    from delfin.agent.benchmark import load_tasks

    task = next(t for t in load_tasks()
                if t.id == "office_one_month_is_totalled_by_the_tool")
    assert task.expected_values, "the figure went back to being a pattern"
    figure = task.expected_values[0]
    assert figure.value == pytest.approx(6070.55)
    # The answer in either convention passes; the year total does not.
    for text in ("Im März 2026: 6.070,55 €.", "March 2026: 6,070.55 EUR"):
        assert figure.judge(text) == "matched", text
    assert figure.judge("Im März 2026: 23.716,25 €.") == "wrong"


def test_the_figure_comes_from_the_tool_not_from_hand_arithmetic():
    """Same discipline as before: the rubric's number is what the tool
    returns, re-derived here rather than trusted."""
    from pathlib import Path

    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
    from delfin.agent.benchmark import load_tasks

    workspace = Path("tests/fixtures/office_workspace")
    if not (workspace / "Buchungen_2026.xlsx").exists():
        pytest.skip("the workbook fixture is not built in this checkout")
    perms = KitToolPermissions(workspace=str(workspace))
    perms.mode = "acceptEdits"
    perms.task_session_id = "figure-judge"
    out = _DocToolExecutor()._execute_sum_column(
        {"path": "Buchungen_2026.xlsx", "column": "Betrag",
         "date_column": "Datum", "period": "2026-03"}, perms)

    task = next(t for t in load_tasks()
                if t.id == "office_one_month_is_totalled_by_the_tool")
    assert task.expected_values[0].judge(out) == "matched", out
