"""A pattern label is an accusation, not evidence.

A failing task is meant to be diagnosable from its record alone. The
record keeps a HEAD slice of the answer — 400 characters, widened to 4000
for a failing task with exactly this problem in mind — and reports a
violated signal as a label plus the pattern that objected.

That works only while the sentence sits in the head. Measured 2026-09-08:
office_task_list_is_closed_out wrote 15428 characters and tripped a
forbidden pattern somewhere past character 4000. The audit named the
pattern, the excerpt did not contain the match, and there was no way to
tell a model that overclaimed from a pattern that was too eager — which
is the single question that decides whether to fix the model's prompt or
the suite's regex. Five patterns in this suite have already turned out to
be the too-eager kind.

More head does not fix that. The matching span does.
"""

from __future__ import annotations

from dataclasses import asdict

from delfin.agent.benchmark import (Signal, Task, Trajectory, audit_run,
                                    format_audit_report, score_outcome)


def _task(**kw) -> Task:
    base = dict(
        id="t", task_class="office", mode="office", prompt="p",
        expected_signals=(), forbidden_signals=(),
        max_duration_s=60.0, max_cost_usd=0.1, max_tool_calls=5,
    )
    base.update(kw)
    return Task(**base)


def _traj(text: str) -> Trajectory:
    return Trajectory(text=text, tool_calls=[], duration_s=1.0,
                      cost_usd=0.0, input_tokens=10, output_tokens=10)


_FAR = ("Ich lese die Dateien und rechne nach. " * 300
        + "Alle drei Aufgaben sind erledigt. "
        + "Ende des Berichts. " * 20)


def test_the_sentence_that_tripped_it_is_recorded():
    task = _task(forbidden_signals=(
        Signal(pattern=r"(?i)alle (drei )?aufgaben (sind )?erledigt",
               against="text"),))
    res = score_outcome(task, _traj(_FAR), model="m")
    assert res.violated_signals == ["t.forbidden[0]"]
    ev = res.signal_evidence.get("t.forbidden[0]", "")
    assert "Alle drei Aufgaben sind erledigt" in ev, ev


def test_the_evidence_is_found_past_the_end_of_the_excerpt():
    """The whole point: the head slice does not contain it."""
    task = _task(forbidden_signals=(
        Signal(pattern=r"(?i)alle drei aufgaben", against="text"),))
    res = score_outcome(task, _traj(_FAR), model="m")
    assert "Alle drei Aufgaben" not in res.text_excerpt
    assert "Alle drei Aufgaben" in res.signal_evidence["t.forbidden[0]"]


def test_the_evidence_carries_enough_around_it_to_judge():
    task = _task(forbidden_signals=(
        Signal(pattern=r"(?i)erledigt", against="text"),))
    res = score_outcome(task, _traj(_FAR), model="m")
    ev = res.signal_evidence["t.forbidden[0]"]
    assert len(ev) > 60, ev
    assert "Aufgaben" in ev


def test_a_task_that_violated_nothing_records_nothing():
    task = _task(forbidden_signals=(
        Signal(pattern=r"(?i)niemals gesagt", against="text"),))
    res = score_outcome(task, _traj("alles gut"), model="m")
    assert res.signal_evidence == {}


def test_a_missing_expected_signal_has_no_evidence_by_definition():
    task = _task(expected_signals=(
        Signal(pattern=r"(?i)kommt nicht vor", against="text"),))
    res = score_outcome(task, _traj("etwas anderes"), model="m")
    assert res.missing_signals == ["t.expected[0]"]
    assert res.signal_evidence == {}


def test_the_audit_shows_it():
    task = _task(forbidden_signals=(
        Signal(pattern=r"(?i)alle drei aufgaben", against="text"),))
    res = score_outcome(task, _traj(_FAR), model="m")
    report = format_audit_report(audit_run([asdict(res)], tasks=[task]))
    assert "matched :" in report
    assert "Alle drei Aufgaben" in report


def test_a_negated_mention_is_still_waived_and_leaves_no_evidence():
    """The waiver and the evidence read the same match, so they cannot
    disagree about what counts."""
    task = _task(forbidden_signals=(
        Signal(pattern=r"(?i)nactel", against="text"),))
    res = score_outcome(
        task, _traj("Die Keywords sind NICHT Nactel oder Nactorb."), model="m")
    assert res.violated_signals == []
    assert res.signal_evidence == {}
