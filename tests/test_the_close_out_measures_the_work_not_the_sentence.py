"""Closing out a task list is what the task is named for.

`office_task_list_is_closed_out` asks for three steps — sum, reconcile,
check the form — planned as a task list and worked off. Its forbidden
signal banned the sentence "alle drei Aufgaben erledigt".

Measured 2026-09-08 on kit.glm-5.3, sixteen tool calls including
sum_column and compare_tables: it summed with the n/a row named,
reconciled with three concrete discrepancies, checked the form — and
said so. q=41, failed, for the sentence.

The gap the ban was reaching for was real: the reconciliation, the second
of the three steps, had no expected signal, so a run that did the other
two could close the list out having skipped it. That is now measured
directly. R-002 is booked at 289,90 and invoiced at 298,90 and R-009
exists only in the invoices — neither is visible without doing the
comparison.

The ban also carried `optional: true`, which the scorer reads only for
EXPECTED signals. On a forbidden one it is a no-op, so the soft signal
its author wrote was scored as a hard one.
"""

from __future__ import annotations

import re

import yaml

from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_TASK = "office_task_list_is_closed_out"


def _task() -> dict:
    for name in ("tasks_office.yaml", "tasks.yaml"):
        data = yaml.safe_load(
            (_ROOT / "delfin" / "agent" / "pack" / "benchmark" / name)
            .read_text(encoding="utf-8")) or {}
        for t in data.get("tasks", []):
            if t["id"] == _TASK:
                return t
    raise AssertionError(f"{_TASK} not found")


def _text_patterns():
    return [re.compile(s["pattern"]) for s in _task()["expected_signals"]
            if s.get("against") == "text"]


def test_the_reconciliation_has_a_signal_of_its_own():
    """Its substance is only discoverable by doing it."""
    joined = " ".join(s["pattern"] for s in _task()["expected_signals"])
    assert "R-002" in joined or "R-009" in joined


def test_an_answer_that_did_the_work_passes_every_text_signal():
    answer = (
        "Monatsabschluss Juni — alle drei Aufgaben erledigt.\n"
        "1. Summe 1.998,40 € über 5 von 6 Zeilen: R-004 hat im Betrag n/a.\n"
        "2. Abgleich: Differenz R-002 (289,90 gebucht vs 298,90 berechnet); "
        "R-009 nur in den Rechnungen.\n"
        "3. kostenstellen_roh.csv: Zeile 4 hat eine leere Spalte Leitung.\n")
    for rx in _text_patterns():
        assert rx.search(answer), rx.pattern


def test_an_answer_that_skipped_the_reconciliation_does_not():
    answer = (
        "Alle drei Aufgaben erledigt.\n"
        "1. Summe 1.998,40 €, R-004 ohne Betrag.\n"
        "3. kostenstellen_roh.csv: eine Zeile hat eine leere Spalte.\n")
    missed = [rx.pattern for rx in _text_patterns() if not rx.search(answer)]
    assert missed, "a report that skipped step 2 passed every signal"


def test_saying_the_work_is_done_is_not_itself_forbidden():
    assert _task().get("forbidden_signals") == []


def test_optional_is_a_no_op_on_a_forbidden_signal():
    """Stated so nobody writes another one expecting it to soften a ban.

    score_outcome consults `optional` only while walking the EXPECTED
    signals; a forbidden match flips success whatever the flag says.
    """
    from delfin.agent.benchmark import Signal, Task, Trajectory, score_outcome

    task = Task(id="t", task_class="office", mode="office", prompt="p",
                expected_signals=(),
                forbidden_signals=(Signal(pattern="(?i)verboten",
                                          against="text", optional=True),),
                max_duration_s=60.0, max_cost_usd=0.1, max_tool_calls=5)
    traj = Trajectory(text="das ist verboten", tool_calls=[], duration_s=1.0,
                      cost_usd=0.0, input_tokens=1, output_tokens=1)
    res = score_outcome(task, traj, model="m")
    assert res.violated_signals == ["t.forbidden[0]"]
    assert res.success is False
