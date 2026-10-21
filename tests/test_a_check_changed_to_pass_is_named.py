"""A check changed to make it pass is named in the answer.

Input: the file changes of a turn (``edit_signals.note_change``, fed by
the write path), its test runs (``note_command``), the user's message
and the answer. Output under test: ``edit_signals.turn_caveat`` and the
engine's turn-end chain.

Pinned here, the measure side:

* in a turn after a red test run, or one that asks for a fix, a changed
  expected value, an added skip/xfail, a loosened tolerance, a removed
  test or test file, a changed grader file (CI, pytest configuration,
  benchmark task, snapshot) and rows removed from a data file are named,
  each with the file;
* no alarm when the user asked for that change, when the test file is
  new (writing tests), when the turn was no fix and nothing failed, when
  the answer gives the reason in the sentence naming the file, and for an
  ordinary product change.

The object side (root cause) is pinned in
test_a_silenced_failure_asks_for_its_cause.py.

The write path reaches the ledger through the real executor
(``edit_file`` and a shell ``rm``), and the engine appends the caveat at
its single exit.
"""

from __future__ import annotations

import json
import textwrap

import pytest

from delfin.agent import edit_signals as E

S = "sess-measure"


@pytest.fixture(autouse=True)
def _fresh():
    E._TURNS.clear()
    yield
    E._TURNS.clear()


def _red():
    E.note_command(S, "pytest tests/test_calc.py", "1 failed, 4 passed in 0.2s")


TEST_OLD = textwrap.dedent("""\
    import pytest
    from calc import area

    def test_area():
        assert area(2) == 4

    def test_negative():
        with pytest.raises(ValueError):
            area(-1)

    def test_close():
        assert area(1.5) == pytest.approx(2.25, rel=1e-9)
""")


def _edit(new, old=TEST_OLD, path="tests/test_calc.py"):
    return E.note_change(S, path, old, new)


@pytest.mark.parametrize("new,kind", [
    (TEST_OLD.replace("== 4", "== 5"), "assertion or expected value changed"),
    (TEST_OLD.replace("def test_area", "@pytest.mark.skip\ndef test_area"),
     "skip/xfail added"),
    (TEST_OLD.replace("def test_area():",
                      "def test_area():\n    pytest.xfail('later')"),
     "skip/xfail added"),
    (TEST_OLD.replace("rel=1e-9", "rel=1e-2"), "tolerance loosened"),
    (TEST_OLD.replace("def test_negative():\n    with pytest.raises(ValueError):"
                      "\n        area(-1)\n", ""), "test removed"),
    (TEST_OLD.replace("    assert area(2) == 4\n", "    area(2)\n"),
     "assertion removed"),
])
def test_a_weakened_test_after_a_red_run_is_named(new, kind):
    _red()
    _edit(new)
    cav = E.turn_caveat(S, "the tests fail, please look", "Done, all green.")
    assert cav.startswith("\n\n[verify] The check was changed")
    assert kind in cav and "tests/test_calc.py" in cav
    assert "Say why the old check was wrong" in cav


def test_a_user_asking_to_fix_without_a_red_run_is_a_fix_turn_too():
    _edit(TEST_OLD.replace("== 4", "== 5"))
    assert "expected value changed" in E.turn_caveat(
        S, "Fix the failing CI please.", "Fixed.")


@pytest.mark.parametrize("path,old,new,kind", [
    (".github/workflows/ci.yml", "    - run: pytest\n",
     "    - run: pytest --ignore=tests/test_calc.py\n", "grader changed"),
    ("pytest.ini", "[pytest]\naddopts = -q\n",
     "[pytest]\naddopts = -q --deselect tests/test_calc.py::test_area\n",
     "grader changed"),
    ("pyproject.toml", "[tool.pytest.ini_options]\ntestpaths = ['tests']\n",
     "[tool.pytest.ini_options]\ntestpaths = ['tests/unit']\n",
     "grader changed"),
    ("bench/tasks/accept_scoring.yaml", "min_score: 0.9\n",
     "min_score: 0.7\n", "grader changed"),
    ("results/runs.csv", "id,ok\na,1\nb,0\nc,1\n", "id,ok\na,1\nc,1\n",
     "rows removed"),
    ("tests/__snapshots__/view.snap", "<b>4</b>\n", "<b>5</b>\n",
     "test data or snapshot changed"),
])
def test_grader_and_data_edits_are_named(path, old, new, kind):
    _red()
    E.note_change(S, path, old, new)
    cav = E.turn_caveat(S, "make CI green", "Done.")
    assert kind in cav, cav
    assert path in cav


def test_deleting_a_test_file_is_named():
    _red()
    E.note_change(S, "tests/test_calc.py", TEST_OLD, None, deleted=True)
    assert "test file deleted" in E.turn_caveat(S, "fix it", "ok")


@pytest.mark.parametrize("request_text", [
    "The expected value in test_area is wrong, change the test to 5.",
    "Update the snapshot tests, the layout changed on purpose.",
    "Relax the tolerance in the test, 1e-9 is too strict for float32.",
    "Der Test ist falsch, bitte anpassen: die Erwartung ändern.",
    "Remove the outlier rows from the CSV and fix the script.",
])
def test_no_alarm_when_the_user_asked_for_that_change(request_text):
    _red()
    _edit(TEST_OLD.replace("== 4", "== 5"))
    E.note_change(S, "results/runs.csv", "a\nb\n", "a\n")
    assert "check was changed" not in E.turn_caveat(S, request_text, "Done.")


def test_no_alarm_for_a_new_test_file():
    _red()
    E.note_change(S, "tests/test_new.py", None, TEST_OLD)
    assert E.turn_caveat(S, "fix the bug and add tests", "Root cause: x "
                         "was never set, because init skipped it.") == ""


def test_no_alarm_when_the_session_edits_a_test_it_created():
    # Writing new tests includes correcting them: the file is the
    # session's own new work, not a measure it was given.
    E.note_change(S, "tests/test_new.py", None, TEST_OLD)
    E.reset_turn(S)
    _red()
    E.note_change(S, "tests/test_new.py", TEST_OLD,
                  TEST_OLD.replace("== 4", "== 5"))
    assert E.measure_edits(S) == []


@pytest.mark.parametrize("old,new,kind", [
    ("    assert r.status == MISSING\n",
     "    assert r.status in (MISSING, OK)\n", "assertion widened"),
    ("    assert f(1) == 2\n",
     "    pytest.skip('not on this host')\n    assert f(1) == 2\n",
     "skip/xfail added"),
    ("    assert x == pytest.approx(1.0, abs=1e-6)\n",
     "    assert x == pytest.approx(1.0, abs=1e-1)\n", "tolerance loosened"),
])
def test_weakening_its_own_new_test_is_still_named(old, new, kind):
    body = "def test_a():\n" + old
    E.note_change(S, "tests/test_own.py", None, body)
    _red()
    E.note_change(S, "tests/test_own.py", body, "def test_a():\n" + new)
    assert kind in [s.kind for s in E.measure_edits(S)]


def test_a_failure_word_in_a_long_brief_is_no_fix_request():
    brief = ("Build a report module. " * 20 + "It must say which runs failed "
             "and which errors they raised.")
    assert len(brief) > 400
    E.note_change(S, "src/report.py", "x = 1\n", "x = 2\n")
    assert E.turn_caveat(S, brief, "Built report.py.") == ""
    assert E.fix_turn("The import fails.", S)


def test_no_alarm_for_added_tests_in_an_existing_file():
    _red()
    _edit(TEST_OLD + "\ndef test_zero():\n    assert area(0) == 0\n")
    assert E.measure_edits(S) == []


def test_no_alarm_when_no_check_failed_and_nobody_asked_for_a_fix():
    _edit(TEST_OLD.replace("== 4", "== 5"))
    assert E.turn_caveat(S, "Rename area to square_area everywhere.",
                         "Renamed.") == ""


def test_no_alarm_when_the_answer_gives_the_reason_beside_the_file():
    _red()
    _edit(TEST_OLD.replace("== 4", "== 5"))
    answer = ("I changed test_calc.py because the old expected value 4 was "
              "computed for the previous definition of area; the spec now "
              "adds 1.")
    assert "check was changed" not in E.turn_caveat(S, "fix the tests", answer)


def test_an_ordinary_product_change_is_no_measure_edit():
    _red()
    E.note_change(S, "src/calc.py", "def area(x):\n    return x * x + 1\n",
                  "def area(x):\n    return x * x\n")
    assert E.measure_edits(S) == []
    assert "check was changed" not in E.turn_caveat(
        S, "fix the failing test", "The bug was an extra +1 in area.")


def test_a_red_run_of_an_earlier_turn_still_counts():
    _red()
    E.reset_turn(S)
    _edit(TEST_OLD.replace("== 4", "== 5"))
    assert E.measure_edits(S, after_red_only=True)
    assert "check was changed" in E.turn_caveat(S, "go on", "done")


@pytest.mark.parametrize("path,expected", [
    ("tests/test_x.py", "test"), ("pkg/foo_test.go", "test"),
    ("web/a.spec.ts", "test"), ("conftest.py", "test"),
    (".github/workflows/ci.yml", "grader"), ("tox.ini", "grader"),
    ("bench/thresholds.json", "grader"), ("pyproject.toml", "config"),
    ("data/rows.csv", "data"), ("src/app.py", "object"),
    ("docs/testing_guide.md", "object"),
])
def test_roles(path, expected):
    assert E.role(path) == expected


# ---------------------------------------------------------------------------
# through the real write path and the engine's single exit
# ---------------------------------------------------------------------------

def test_the_write_path_feeds_the_ledger(tmp_path):
    from delfin.agent.api_client import KitToolPermissions, _doc_executor
    ws = tmp_path / "ws"
    (ws / "tests").mkdir(parents=True)
    test_file = ws / "tests" / "test_calc.py"
    test_file.write_text(TEST_OLD)
    perms = KitToolPermissions(workspace=str(ws), mode="acceptEdits")
    perms.task_session_id = S
    _doc_executor.execute("read_file", {"path": str(test_file)}, perms)
    out = _doc_executor.execute("edit_file", {
        "path": str(test_file), "old_string": "== 4",
        "new_string": "== 5"}, perms)
    assert "Edited" in out, out
    kinds = [s.kind for s in E.measure_edits(S)]
    assert kinds == ["assertion or expected value changed"]


def test_the_engine_appends_the_caveat_at_its_exit(tmp_path):
    from unittest.mock import MagicMock, patch
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick")
    client._permissions = KitToolPermissions(workspace=str(tmp_path))
    client._permissions.task_session_id = S
    engine._turn_request = "the test fails, fix it"
    engine.messages = [{"role": "assistant", "content": "All green now."}]
    engine._exec_pending = [("bash", json.dumps(
        {"command": "pytest tests/test_calc.py"}))]
    engine._commit_exec_command("bash", "FAILED tests/test_calc.py::test_area"
                                "\n1 failed, 2 passed in 0.1s")
    _edit(TEST_OLD.replace("== 4", "== 5"))
    out = engine._append_answer_caveats("All green now.", functional=[],
                                        ambiguous=[])
    assert "[verify] The check was changed" in out
    assert "[verify] The check was changed" in engine.messages[-1]["content"]


def test_a_failing_ci_check_read_through_gh_is_a_red_check(tmp_path):
    from unittest.mock import MagicMock, patch
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick")
    client._permissions = KitToolPermissions(workspace=str(tmp_path))
    client._permissions.task_session_id = S
    engine._note_check_run("bash", {"command": "gh pr checks 12"},
                           "test (3.11)\tpass\t2m\nlint\tpass\t1m\n")
    assert not E.red_seen(S)
    engine._note_check_run("bash", {"command": "gh pr checks 12"},
                           "test (3.11)\tfail\t2m\nlint\tpass\t1m\n")
    assert E.red_seen(S)
