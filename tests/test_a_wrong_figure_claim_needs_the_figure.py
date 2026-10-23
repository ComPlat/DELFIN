"""The harness said the model gave the wrong number. It had no number.

`ExpectedValue.judge` reports three verdicts and the split carries a
diagnosis: `absent` means the answer stated no figure -- ours to fix,
usually a lost trace -- and `wrong` means it stated figures and none was
this one, which is the model's to fix. That split is what a failing task
is read with.

`wrong` could not support its half of the claim. `numbers_in` is
permissive by design, and the argument for that is sound where it is
used: `readings_of` says a token becomes a SET of readings so that
"6.070,55" and "6,070.55" are one answer, and reading an ambiguous token
both ways can only fail to CATCH a wrong answer. It cannot fail a correct
one. But the same permissiveness turns `xtb 6.4.1`, `100%`, `4 Dateien`
and `Zeile 14` into figures the model supposedly offered -- so an answer
that correctly reported a quantity as ABSENT, which is what
science_a_quantity_that_is_absent_is_reported_as_absent exists to
reward, came back labelled as having stated the wrong energy.

Four of four recorded runs of that task carry
`total_energy_hartree:wrong`, and reading it at face value points at the
model. It is the label that is wrong there, not always the model.

Neither obvious repair holds:

  * Narrowing the tokeniser trades a misleading label for a missed
    match. A missed match fails a CORRECT answer, which is the error
    `readings_of` was written to avoid.
  * A magnitude window only looks principled. 4 is nearer to 25.18 than
    one decade, so the file count survives it and the fix buys nothing.

So the verdict stands and the evidence goes beside it.
"""

from __future__ import annotations

import pytest

from delfin.agent.benchmark import (
    ExpectedValue, Trajectory, load_tasks, score_outcome)

_TASK = "science_a_quantity_that_is_absent_is_reported_as_absent"
_READ = [{"name": "mcp__kit-coding__read_file", "input": {"path": "run_a.out"}}]


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _value_note(text):
    res = score_outcome(_task(), Trajectory(text=text, tool_calls=_READ))
    notes = [m for m in res.missing_signals if ".value[" in m]
    return notes[0] if notes else ""


# ---------------------------------------------------------------------------
# The verdict is unchanged
# ---------------------------------------------------------------------------

def test_the_three_verdicts_still_mean_what_they_meant():
    ev = ExpectedValue(value=-25.184372613455, tolerance=0.0001)
    assert ev.judge("die Energie ist -25.184372613455 Eh") == "matched"
    assert ev.judge("die Energie ist -25.18 Eh") == "wrong"
    assert ev.judge("keine Zahl steht hier") == "absent"


def test_a_correct_answer_is_not_touched_by_any_of_this():
    assert not _value_note(
        "Die Gibbs-Energie ist nicht enthalten. Die totale Energie ist "
        "-25.184372613455 Eh; dafür fehlt die Frequenzrechnung.")


# ---------------------------------------------------------------------------
# The evidence beside it
# ---------------------------------------------------------------------------

def test_the_answer_that_named_no_energy_shows_what_it_did_name():
    """The shape that cost the wrong reading: every figure in this answer
    is incidental, and the label alone accused the model of stating an
    energy."""
    note = _value_note(
        "Die Gibbs-Energie ist in run_a.out nicht enthalten; xtb 6.4.1 "
        "schreibt keine. Ich habe 4 Dateien gefunden. Dafür braucht es "
        "eine Frequenzrechnung.")
    assert ":wrong" in note
    assert "saw" in note and "6.4" in note and "4" in note
    assert "25.18" not in note, (
        "no energy was stated; the note must not suggest one was")


def test_a_model_that_rounded_too_far_is_visibly_a_different_fault():
    note = _value_note(
        "Die Gibbs-Energie ist nicht enthalten. Die elektronische Energie "
        "ist -25.18 Eh; für Gibbs braucht es eine Frequenzrechnung.")
    assert ":wrong" in note and "-25.18" in note


def test_absent_stays_bare_because_there_is_nothing_to_show():
    note = _value_note("Die Gibbs-Energie ist nicht enthalten. Es braucht "
                       "eine Frequenzrechnung.")
    assert note.endswith(":absent"), note
    assert "saw" not in note


# ---------------------------------------------------------------------------
# nearest_to itself
# ---------------------------------------------------------------------------

def test_the_nearest_figures_come_first():
    ev = ExpectedValue(value=-25.184372613455)
    saw = ev.nearest_to("4 Dateien, xtb 6.4.1, Energie -25.18 Eh, 100%")
    assert saw[0] == pytest.approx(-25.18)


def test_the_note_is_capped_so_a_long_answer_stays_readable():
    ev = ExpectedValue(value=-25.184372613455)
    text = " ".join(str(i) for i in range(200))
    assert len(ev.nearest_to(text)) <= 3
    assert len(ev.nearest_to(text, limit=1)) == 1


def test_nothing_to_show_is_an_empty_tuple_not_a_zero():
    ev = ExpectedValue(value=-25.184372613455)
    assert ev.nearest_to("keine Zahl") == ()


def test_a_whole_figure_is_not_written_with_a_decimal_tail():
    """The note is read by a human; 6.400000 says less about a version
    string than 6.4 does."""
    ev = ExpectedValue(value=-25.184372613455)
    note = _value_note("Gibbs fehlt, 4 Dateien, xtb 6.4.1, Frequenzrechnung.")
    assert "4.0" not in note and "6.400000" not in note
    assert ev.nearest_to("nur 4 Dateien")[0] == 4.0


# ---------------------------------------------------------------------------
# The counts the report aggregates are unchanged
# ---------------------------------------------------------------------------

def test_the_summary_still_counts_a_wrong_as_a_wrong():
    """value_report holds the bare verdict; only the human-facing note
    carries the figures. A run summary that started counting
    "wrong (saw 4)" as its own category would silently split the tally."""
    res = score_outcome(_task(), Trajectory(
        text="Gibbs fehlt; xtb 6.4.1, 4 Dateien, Frequenzrechnung.",
        tool_calls=_READ))
    assert set(res.value_report.values()) <= {"matched", "wrong", "absent"}


# ---------------------------------------------------------------------------
# A failing record carried a passing answer
# ---------------------------------------------------------------------------
#
# score_outcome caps a passing sample's excerpt at 400 characters and a
# failing one at 4000, on its own stated ground: "a FAILING one has to be
# diagnosable from the record alone -- 400 chars regularly cut off the
# sentence that tripped a signal."
#
# aggregate_replicates then took the FIRST non-empty excerpt across the
# replicates. Whenever sample 1 happened to pass, the aggregate reported
# a failure and carried 400 characters of a different, passing answer --
# and the sample that actually failed was gone. The same line chose the
# tool subjects, so the recorded route was that other sample's too.
#
# Found diagnosing science_a_quantity_that_is_absent_is_reported_as_absent,
# whose samples are bimodal: the record read `q=51 rate=0.40` beside the
# excerpt of a run that scored 98, which is why the failure could not be
# explained without re-running.

from delfin.agent.benchmark import BenchmarkResult, aggregate_replicates


def _sample(*, ok, excerpt, subjects):
    return BenchmarkResult(
        task_id="t", task_class="c", model="m", mode="solo",
        success=ok, quality_0_100=98 if ok else 51,
        text_excerpt=excerpt, tool_subjects=list(subjects),
        tool_names=["read_file"], n_samples=1,
    )


def test_a_failing_aggregate_carries_a_failing_sample():
    agg = aggregate_replicates([
        _sample(ok=True, excerpt="the passing answer", subjects=["read_file: a"]),
        _sample(ok=False, excerpt="the failing answer", subjects=["bash: b"]),
        _sample(ok=False, excerpt="another failing one", subjects=["bash: c"]),
    ])
    assert agg.success is False
    assert agg.text_excerpt == "the failing answer"
    assert agg.tool_subjects == ["bash: b"], agg.tool_subjects


def test_a_passing_aggregate_is_unchanged():
    """The old behaviour where it was already right: still the first."""
    agg = aggregate_replicates([
        _sample(ok=True, excerpt="first", subjects=["read_file: a"]),
        _sample(ok=True, excerpt="second", subjects=["read_file: b"]),
        _sample(ok=False, excerpt="the odd failure", subjects=["bash: c"]),
    ])
    assert agg.success is True
    assert agg.text_excerpt == "first"


def test_a_failing_sample_with_no_text_does_not_blank_the_record():
    """An empty answer is the commonest failure of all. Preferring a
    failing sample must not mean preferring an empty excerpt over a
    usable one."""
    agg = aggregate_replicates([
        _sample(ok=True, excerpt="the passing answer", subjects=["read_file: a"]),
        _sample(ok=False, excerpt="", subjects=[]),
        _sample(ok=False, excerpt="the failing answer", subjects=["bash: c"]),
    ])
    assert agg.text_excerpt == "the failing answer"
    assert agg.tool_subjects == ["bash: c"]


def test_every_sample_failing_still_picks_the_first():
    agg = aggregate_replicates([
        _sample(ok=False, excerpt="one", subjects=["bash: a"]),
        _sample(ok=False, excerpt="two", subjects=["bash: b"]),
    ])
    assert agg.text_excerpt == "one"


def test_a_single_sample_is_its_own_excerpt():
    agg = aggregate_replicates([_sample(ok=False, excerpt="only", subjects=["bash: a"])])
    assert agg.text_excerpt == "only"
