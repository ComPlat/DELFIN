"""The office benchmark suite: loadable, well-formed, and actually run.

A benchmark suite that has to be registered somewhere is a suite that
gets written and then silently never runs. These tests pin that the
office tasks are discovered by existing, that their rubrics compile, and
that their fixture workspace is restored between attempts — without
that, a later attempt scores against rows an earlier one changed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks

_FIXTURES = (Path(__file__).resolve().parent / "fixtures" / "office_workspace")


def _office_tasks():
    return [t for t in load_tasks() if getattr(t, "task_class", "") == "office"]


def test_the_office_suite_is_discovered_without_being_registered():
    assert len(_office_tasks()) >= 6


def test_hand_written_suites_are_picked_up_too():
    """The glob used to cover only auto-generated siblings."""
    import inspect

    from delfin.agent import benchmark

    source = inspect.getsource(benchmark.load_tasks)
    assert 'glob("tasks_*.yaml")' in source


def test_every_office_task_runs_in_office_mode():
    for task in _office_tasks():
        assert task.mode == "office", task.id


def test_task_ids_stay_unique_across_all_suites():
    ids = [t.id for t in load_tasks()]
    assert len(ids) == len(set(ids))


def test_every_rubric_pattern_compiles():
    for task in _office_tasks():
        for signal in list(task.expected_signals) + list(task.forbidden_signals):
            re.compile(signal.pattern)


def test_each_task_carries_a_budget():
    for task in _office_tasks():
        assert task.max_cost_usd > 0, task.id
        assert task.max_duration_s > 0, task.id
        assert task.max_tool_calls > 0, task.id


def test_the_fixture_files_exist():
    for name in ("buchungen.csv", "rechnungen.csv", "inventar.csv",
                 "kostenstellen_roh.csv"):
        assert (_FIXTURES / name).is_file(), name


def test_the_fixture_workspace_is_restored_between_attempts():
    """Office tasks read and edit the same tables. Without a restore a
    later attempt scores against what an earlier one changed."""
    from delfin.agent import benchmark_runner

    rels = [str(p) for p in benchmark_runner._BEHAVIOR_WS_RELS]
    assert any("office_workspace" in r for r in rels)


def test_the_fixtures_carry_the_traps_the_rubrics_check():
    """A rubric that no longer matches its data passes for the wrong
    reason, so the two are pinned against each other."""
    bookings = (_FIXTURES / "buchungen.csv").read_bytes().decode("cp1252")
    assert "n/a" in bookings                     # a row with no amount
    assert bookings.count("R-005") == 2          # a duplicate key
    assert "1.234,50" in bookings                # decimal comma

    invoices = (_FIXTURES / "rechnungen.csv").read_text(encoding="utf-8")
    assert "Belegnummer" in invoices             # a different vocabulary
    assert "298.90" in invoices                  # the transposition
    assert "R-009" in invoices                   # only on this side

    inventory = (_FIXTURES / "inventar.csv").read_text(encoding="utf-8")
    assert "8.986" in inventory                  # ambiguous throughout

    raw = (_FIXTURES / "kostenstellen_roh.csv").read_bytes().decode("cp1252")
    assert raw.splitlines()[3].count(";") == 3   # one field too many


def test_the_ambiguous_column_really_cannot_be_decided():
    """The task forbids a total. That is only fair if nothing in the
    column settles the reading."""
    from delfin.agent import office

    values = [line.split(",")[2] for line in
              (_FIXTURES / "inventar.csv").read_text(
                  encoding="utf-8").splitlines()[1:] if line]
    assert office.detect_number_convention(values)[0] == office.AMBIGUOUS


def test_no_office_task_expects_a_tool_that_does_not_exist():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI

    known = {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}
    for task in _office_tasks():
        for signal in task.expected_signals:
            if getattr(signal, "against", "") == "tool_name":
                assert signal.pattern in known, (
                    f"{task.id} expects a tool named {signal.pattern!r}")


def test_office_tasks_start_inside_the_office_folder():
    """Office mode is defined by working inside one folder, so a task in
    that mode has to start there. Pointed at the repository root, the
    fixtures sat three levels down and whether the model found them was
    luck — one run cited the full path, the next answered "the file is
    not in the working directory", and that scored as a failure of a
    rubric it never got to exercise."""
    import inspect

    from delfin.agent import benchmark_runner

    source = inspect.getsource(benchmark_runner._build_engine) if hasattr(
        benchmark_runner, "_build_engine") else inspect.getsource(
        benchmark_runner)
    assert 'if (mode or "") == "office":' in source
    assert '"office_workspace"' in source


def test_the_office_prompts_name_files_the_way_a_user_would():
    """Bare file names, not fixture paths — the workspace is what makes
    them resolvable, and a prompt carrying a repository path would be
    testing something no user ever types."""
    for task in _office_tasks():
        assert "tests/fixtures" not in task.prompt, task.id
