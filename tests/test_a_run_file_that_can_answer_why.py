"""The recorded run kept tool names and dropped what they touched.

Reading finished runs is the cheapest way to find where a model
struggles: the models have already told us, in the traces, and nobody has
to spend a token to hear it again. On 2026-09-09 six framework faults
came out of existing run files that way.

It only works if the record carries the subject. It did not. One task
took fifty tool calls to total a column and the file said `edit_file x12,
bash x11` — which files, which commands, it could not say. The cause was
a prompt sentence sending column arithmetic to a shell, and the only
reason it was found is that the console log still existed. Once that is
gone the run file cannot answer the question its own numbers raise.

`tool_subjects` is `name: subject` per call, in order, using the same
notion of "the argument that identifies a call" that `_lead_with_subject`
already applies for matching. Never the content of a write: that is what
made keeping the inputs unthinkable in the first place, and it is why
this costs about a hundred bytes a task rather than a megabyte.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.benchmark import (
    Task, Trajectory, aggregate_replicates, score_outcome, tool_subject,
)


def _task():
    return Task(id="probe", task_class="science_analysis", mode="solo",
                prompt="x")


# ---------------------------------------------------------------------------
# What a subject is
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("call,expected", [
    ({"name": "mcp__kit-coding__bash", "input": {"command": "ls -la"}},
     "bash: ls -la"),
    ({"name": "write_file", "input": {"path": "summary.py"}},
     "write_file: summary.py"),
    ({"name": "sum_column", "input": {"path": "b.csv", "column": "Betrag"}},
     "sum_column: b.csv"),
    ({"name": "grep_file", "input": {"pattern": "HOMO", "path": "a.out"}},
     "grep_file: a.out"),
])
def test_the_subject_is_the_argument_that_identifies_the_call(call, expected):
    assert tool_subject(call) == expected


def test_a_write_never_carries_its_content():
    """The reason the arguments were dropped in the first place."""
    subject = tool_subject(
        {"name": "write_file",
         "input": {"path": "s.py", "content": "x" * 100_000}})
    assert subject == "write_file: s.py"
    assert "xxxx" not in subject


def test_a_call_with_no_identifying_argument_is_still_distinguishable():
    subject = tool_subject({"name": "task_update",
                            "input": {"status": "done", "note": "ok"}})
    assert subject.startswith("task_update: ")
    assert "status" in subject


def test_a_long_subject_is_truncated():
    subject = tool_subject(
        {"name": "bash", "input": {"command": "echo " + "a" * 500}})
    assert len(subject) < 120


def test_whitespace_is_flattened_so_one_call_is_one_line():
    subject = tool_subject(
        {"name": "bash", "input": {"command": "line one\n  line two"}})
    assert "\n" not in subject
    assert "line one line two" in subject


# ---------------------------------------------------------------------------
# What reaches the file
# ---------------------------------------------------------------------------

def test_the_route_is_recorded_in_order():
    traj = Trajectory(text="fertig", tool_calls=[
        {"name": "mcp__kit-coding__list_files", "input": {"path": "."}},
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "awk -F';' '{s+=$5}END{print s}' b.csv"}},
        {"name": "write_file", "input": {"path": "summary.py",
                                         "content": "x" * 9000}},
    ])
    res = score_outcome(_task(), traj)
    assert res.tool_subjects == [
        "list_files: .",
        "bash: awk -F';' '{s+=$5}END{print s}' b.csv",
        "write_file: summary.py",
    ]


def test_the_detour_is_visible_without_the_console_log():
    """The case this exists for: a shell one-liner doing arithmetic a
    tool does, which the names alone cannot show."""
    res = score_outcome(_task(), Trajectory(text="1998.40", tool_calls=[
        {"name": "bash",
         "input": {"command": "awk -F';' '$2 ~ /06.2026/ {s+=$5} END{print s}' b.csv"}},
    ]))
    joined = " ".join(res.tool_subjects)
    assert "awk" in joined and "06.2026" in joined


def test_it_costs_about_a_hundred_bytes():
    """Cheap enough that no run has to weigh whether to keep it."""
    res = score_outcome(_task(), Trajectory(text="x", tool_calls=[
        {"name": "read_file", "input": {"path": f"run_{i}.out"}}
        for i in range(10)]))
    assert len(json.dumps(res.tool_subjects)) < 500


def test_a_run_with_no_tool_calls_records_nothing():
    res = score_outcome(_task(), Trajectory(text="nur Text"))
    assert res.tool_subjects == []


# ---------------------------------------------------------------------------
# Repeats
# ---------------------------------------------------------------------------

def test_repeats_keep_one_route_rather_than_interleaving_three():
    """A union of three replicates would read as a route nobody took."""
    def _run(marker):
        return score_outcome(_task(), Trajectory(text="x", tool_calls=[
            {"name": "bash", "input": {"command": f"echo {marker}"}}]))

    agg = aggregate_replicates([_run("a"), _run("b"), _run("c")])
    assert agg.tool_subjects == ["bash: echo a"]


def test_the_first_sample_that_called_anything_is_the_one_kept():
    empty = score_outcome(_task(), Trajectory(text="x"))
    called = score_outcome(_task(), Trajectory(text="x", tool_calls=[
        {"name": "read_file", "input": {"path": "a.out"}}]))
    agg = aggregate_replicates([empty, called])
    assert agg.tool_subjects == ["read_file: a.out"]
