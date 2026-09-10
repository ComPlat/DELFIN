"""The calc tools indexed the USER's archive, so nothing could ask them anything.

`search_calcs` and `calc_summary` build their index from the real calc/
and archive/ folders -- 777 calculations on this machine, whose contents
nobody chose. "Which of my runs used which method" is daily work for the
scientist this agent is for, and it was completely unmeasured, because
the right answer differed per developer and changed whenever a real run
finished.

`setup` builds nine chosen calculations and the runner points the tools
at them for the length of the attempt, restoring the previous setting
afterwards -- the executor is a module-level singleton, so a missing
restore would silently redirect every later task and any live session
sharing the process.

This file pins the fixture's arithmetic and the rubric in both
directions. The forbidden pattern is the exposed one again: the right
answer and the wrong answer both name `calc_d`, and they differ only by
whether it is called finished.
"""

from __future__ import annotations

import json
import re

import pytest

from delfin.agent.benchmark import load_tasks
from delfin.agent.benchmark_runner import _fixture_calc_dirs, run_setup

TASK_ID = "science_the_archive_answers_which_method"


@pytest.fixture(scope="module")
def task():
    for t in load_tasks():
        if t.id == TASK_ID:
            return t
    pytest.fail(f"{TASK_ID} is not in the suite")


@pytest.fixture
def built(tmp_path):
    ok, out = run_setup("a_small_calc_archive.py", tmp_path)
    assert ok, out
    return tmp_path


# ---------------------------------------------------------------------------
# The fixture says what the task claims it says
# ---------------------------------------------------------------------------

def test_the_tools_see_exactly_the_fixture(built):
    from delfin.agent.api_client import KitToolPermissions, _doc_executor

    with _fixture_calc_dirs(built):
        perms = KitToolPermissions(mode="default", workspace=str(built))
        out = json.loads(_doc_executor.execute("calc_summary", {}, perms))
    assert out["total_calculations"] == 9
    assert out["by_source"] == {"calc": 4, "archive": 5}
    functionals = dict(out["functionals"])
    assert functionals["PBE0"] == 5
    assert functionals["B3LYP"] == 3
    assert str(built) in json.dumps(out["indexed_from"])


def test_the_expected_values_match_the_fixture(task, built):
    from delfin.agent.api_client import KitToolPermissions, _doc_executor

    with _fixture_calc_dirs(built):
        perms = KitToolPermissions(mode="default", workspace=str(built))
        out = json.loads(_doc_executor.execute("calc_summary", {}, perms))
    wanted = {v.label: v for v in task.expected_values}
    assert wanted["n_calculations"].judge(
        str(out["total_calculations"])) == "matched"
    assert wanted["pbe0_count"].judge(
        str(dict(out["functionals"])["PBE0"])) == "matched"


def test_the_users_own_archive_is_put_back(built):
    """The executor is a module-level singleton. Without the restore, the
    first calc task would redirect every later one -- and any live session
    sharing the process."""
    from delfin.agent.api_client import _doc_executor

    before = dict(getattr(_doc_executor, "_calc_dirs", {}))
    with _fixture_calc_dirs(built):
        assert getattr(_doc_executor, "_calc_dirs") != before
    assert dict(getattr(_doc_executor, "_calc_dirs", {})) == before


def test_a_workspace_without_the_fixture_is_untouched(tmp_path):
    from delfin.agent.api_client import _doc_executor

    before = dict(getattr(_doc_executor, "_calc_dirs", {}))
    with _fixture_calc_dirs(tmp_path):
        assert dict(getattr(_doc_executor, "_calc_dirs", {})) == before


def test_the_setup_checks_its_own_precondition(tmp_path):
    (tmp_path / "calc_archive").mkdir()
    ok, out = run_setup("a_small_calc_archive.py", tmp_path)
    assert not ok
    assert "refusing to overwrite" in out


def test_an_unfinished_run_has_no_result_files(built):
    d = built / "calc_archive" / "calc" / "calc_d"
    assert (d / "CONTROL.txt").is_file()
    assert not (d / "DELFIN_Data.json").exists()
    assert not any(d.glob("*.out"))


def test_the_run_without_a_control_file_is_still_indexed(built):
    """arch_e has only an .inp. If the indexer's fallback breaks, the
    count silently becomes 8 and the task starts failing for a reason
    that has nothing to do with the model."""
    d = built / "calc_archive" / "archive" / "arch_e"
    assert not (d / "CONTROL.txt").exists()
    assert any(d.glob("*.inp"))


# ---------------------------------------------------------------------------
# The rubric, in both directions
# ---------------------------------------------------------------------------

def test_the_tool_is_required_not_a_shell_walk(task):
    pat = task.expected_signals[0].pattern
    assert re.search(pat, 'TOOL: calc_summary({})')
    assert re.search(pat, 'TOOL: search_calcs({"functional": "PBE0"})')
    # The whole point: the counts come from the tool, not from the model
    # walking the folders itself. A namespaced name reaches this pattern
    # already stripped -- the scorer normalises it -- so only the bare
    # spelling is asserted here.
    assert not re.search(pat, 'TOOL: bash({"command": "ls ~/calc | wc -l"})')
    assert not re.search(pat, 'TOOL: list_files({"path": "calc_archive"})')


@pytest.mark.parametrize("answer", [
    "Insgesamt 9 Rechnungen. Am häufigsten PBE0 (5x). calc_d läuft noch.",
    "9 Läufe; PBE0 5 mal. calc_d ist noch nicht fertig.",
    "Du hast 9 Rechnungen, PBE0 dominiert mit 5. calc_d ist unfertig.",
])
def test_a_correct_answer_is_not_punished(task, answer):
    forb = task.forbidden_signals[0].pattern
    assert not re.search(forb, answer), answer


@pytest.mark.parametrize("answer", [
    "9 Rechnungen, PBE0 5x, calc_d ist fertig.",
    "Alle sind abgeschlossen, calc_d ebenfalls fertig.",
])
def test_calling_the_unfinished_run_finished_is_caught(task, answer):
    forb = task.forbidden_signals[0].pattern
    assert re.search(forb, answer), answer


def test_the_setup_script_is_the_one_the_task_names(task):
    from delfin.agent.benchmark_runner import setup_path

    assert task.setup == "a_small_calc_archive.py"
    assert setup_path(task.setup).is_file()


# ---------------------------------------------------------------------------
# The workspace that is not there
# ---------------------------------------------------------------------------
#
# Only the science and project families get a fixture workspace; the rest
# of the suite passes None. The first version of _fixture_calc_dirs
# indexed straight into Path(None), and because run_task catches whatever
# the attempt raises, every dashboard, office and behaviour task came
# back as `_run_once raised: argument should be a str or an os.PathLike
# object`. Unmeasured rather than scored as a model failure, which is the
# honest half — but a whole suite of them, from a guard that was supposed
# to be inert for exactly those tasks.

def test_a_task_without_a_workspace_is_left_alone():
    from delfin.agent.api_client import _doc_executor

    before = dict(getattr(_doc_executor, "_calc_dirs", {}))
    with _fixture_calc_dirs(None):
        assert dict(getattr(_doc_executor, "_calc_dirs", {})) == before
    assert dict(getattr(_doc_executor, "_calc_dirs", {})) == before


def test_an_empty_workspace_string_is_left_alone_too():
    from delfin.agent.api_client import _doc_executor

    before = dict(getattr(_doc_executor, "_calc_dirs", {}))
    with _fixture_calc_dirs(""):
        assert dict(getattr(_doc_executor, "_calc_dirs", {})) == before


def test_a_task_with_no_workspace_still_runs_end_to_end():
    """The shape the whole suite hit: the guard must be invisible to a
    task that has no workspace at all."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task

    task = Task(id="probe", task_class="dashboard_nav", mode="dashboard",
                prompt="p", expected_signals=(), forbidden_signals=(),
                max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2)

    class _Engine:
        cost_usd = 0.0

    ticks = iter([100.0, 101.5, 101.6])
    res = br.run_task(
        task, model="m",
        engine_factory=lambda *a, **k: _Engine(),
        run_once=lambda e, p, **k: {"text": "ok", "tool_calls": [],
                                    "input_tokens": 1, "output_tokens": 1,
                                    "error": ""},
        clock=lambda: next(ticks))
    assert not res.error, res.error
    assert "PathLike" not in (res.error or "")
