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
    import pathlib

    from delfin.agent.benchmark_runner import workspace_for

    root = pathlib.Path(__file__).resolve().parents[1]
    # Was a source-text assertion on the condition that implemented this,
    # so moving the logic into a named function broke it while the
    # behaviour was unchanged. It asks the function now.
    assert workspace_for(root, mode="office") == (
        root / "tests" / "fixtures" / "office_workspace")
    assert workspace_for(root, task_class="office") == (
        root / "tests" / "fixtures" / "office_workspace")


def test_the_office_prompts_name_files_the_way_a_user_would():
    """Bare file names, not fixture paths — the workspace is what makes
    them resolvable, and a prompt carrying a repository path would be
    testing something no user ever types."""
    for task in _office_tasks():
        assert "tests/fixtures" not in task.prompt, task.id


# ---------------------------------------------------------------------------
# The forbidden pattern that failed a correct answer
# ---------------------------------------------------------------------------
#
# Eleventh time in this suite's history that a rubric measured spelling
# rather than behaviour, so this one is pinned on the real answers.
# kit.glm-5.3 named the cost centre, named the budget and explained the
# A3:A5 merge, then added "(Wahrung in der Datei nicht angegeben)" about
# the CURRENCY, and a word-list pattern for "the row is unassigned"
# matched it. The same list matched every correct answer that explained
# why the row's own cell is blank -- which is the explanation the task
# wants -- and missed two of the wrong answers it existed to catch.

_MERGED = "office_merged_block_is_not_read_as_unassigned"

_UNASSIGNED_CLAIMS = (
    'Die Position „Wartung“ hat keine Kostenstelle.',
    'Wartung ist keiner Kostenstelle zugeordnet.',
    'Für Wartung ist keine Kostenstelle angegeben.',
    'Die Kostenstelle für Wartung ist nicht angegeben.',
    'Die Zelle ist leer, die Position ist damit nicht zugeordnet.',
    'Zu dieser Zeile fehlt die Kostenstelle.',
    'Die Kostenstelle ist nicht erkennbar.',
)

_CORRECT_ANSWERS = (
    # The measured one, verbatim in substance.
    'Die Position „Wartung" gehört zur Kostenstelle 4711 mit einem Budget '
    'von 7.300 (Währung in der Datei nicht angegeben).',
    'Die Kostenstelle ist 4711, die Währung nicht angegeben.',
    'Wartung gehört zu 4711; ein Währungssymbol fehlt in der Datei.',
    'Die eigene Zelle in Zeile 5 ist leer, weil A3:A5 verbunden sind — '
    'die Kostenstelle ist 4711.',
    'Kostenstelle 4711 (Anorganische Chemie), Budget 7.300.',
)


def _merged_forbidden():
    from delfin.agent.benchmark import load_tasks

    task = next(t for t in load_tasks() if t.id == _MERGED)
    assert len(task.forbidden_signals) == 1, "the rubric changed shape"
    return task.forbidden_signals[0].pattern


def test_the_unassigned_claim_is_still_caught():
    import re

    pat = _merged_forbidden()
    missed = [t for t in _UNASSIGNED_CLAIMS if not re.search(pat, t)]
    assert not missed, "a wrong answer walks through: " + "; ".join(missed)


def test_a_correct_answer_is_not_failed_for_a_word_it_used():
    import re

    pat = _merged_forbidden()
    hits = []
    for text in _CORRECT_ANSWERS:
        m = re.search(pat, text)
        if m:
            hits.append(f"{m.group(0)!r} in {text[:50]!r}")
    assert not hits, "correct answers scored as the naive reading: " + \
        "; ".join(hits)


# ---------------------------------------------------------------------------
# The first office task that changes a file
# ---------------------------------------------------------------------------

_CORRECT = "office_a_record_is_corrected_by_its_key"


def _correct_task():
    from delfin.agent.benchmark import load_tasks

    return next(t for t in load_tasks() if t.id == _CORRECT)


def _edit_by_key_traj(text: str):
    from delfin.agent.benchmark import Trajectory

    return Trajectory(text=text, tool_calls=[
        {"name": "mcp__delfin-docs__read_document",
         "input": {"path": "Buchungen_2026.xlsx"}},
        {"name": "mcp__kit-coding__edit_sheet",
         "input": {"path": "Buchungen_2026.xlsx", "key_column": "Beleg",
                   "updates": [{"key": "R-014",
                                "set": {"Betrag": "1.265,85"}}]}},
    ])


def test_the_right_edit_satisfies_every_expected_signal():
    from delfin.agent.benchmark import _signal_matches

    traj = _edit_by_key_traj(
        "R-014 steht jetzt auf 1.265,85 € (vorher 265,85 €). Die Änderung "
        "ist gesichert und lässt sich mit undo_changes zurücknehmen.")
    unmatched = [i for i, s in enumerate(_correct_task().expected_signals)
                 if not _signal_matches(s, traj)]
    assert not unmatched, f"unmatched: expected{unmatched}"


def test_a_cell_coordinate_does_not_satisfy_the_by_key_signal():
    """The rule the prompt states: the sheet has seven hidden rows, so
    what a reader counts and what the file numbers are different things."""
    from delfin.agent.benchmark import Trajectory, _signal_matches

    by_cell = Trajectory(
        text="E15 auf 1.265,85 gesetzt; per undo_changes rücknehmbar.",
        tool_calls=[{"name": "mcp__kit-coding__edit_sheet",
                     "input": {"path": "Buchungen_2026.xlsx",
                               "edits": [{"cell": "E15",
                                          "value": "1.265,85"}]}}])
    task = _correct_task()
    key_signal = next(s for s in task.expected_signals
                      if "key_column" in s.pattern)
    assert not _signal_matches(key_signal, by_cell)
    # Saying it is not doing it.
    talked = Trajectory(text="Ich adressiere die Zeile über key_column.")
    assert not _signal_matches(key_signal, talked)


def test_rewriting_the_workbook_in_a_shell_is_the_forbidden_route():
    from delfin.agent.benchmark import Trajectory, _signal_matches

    task = _correct_task()
    shell = Trajectory(text="Fertig.", tool_calls=[
        {"name": "mcp__kit-coding__bash",
         "input": {"command": "python -c \"import openpyxl; "
                              "wb=openpyxl.load_workbook('Buchungen_2026.xlsx')\""}}])
    assert any(_signal_matches(s, shell) for s in task.forbidden_signals)
    # And the correct route is not caught by it.
    ok = _edit_by_key_traj("R-014 auf 1.265,85 gesetzt, per Backup rücknehmbar.")
    assert not any(_signal_matches(s, ok) for s in task.forbidden_signals)
