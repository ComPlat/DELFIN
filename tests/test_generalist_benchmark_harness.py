"""Harness support for the ``generalist`` benchmark category.

Four mechanisms, each pinned here against the real loader, scorer, runner
and CLI parser (no model call):

* ``prompt_en``: a task carries an English prompt beside the German one;
  ``english_variant`` turns it into a task with its own id (``<id>.en``)
  so the two languages never share a history row or a baseline.
* ``select_tasks``: ``--task`` accepts shell-style globs
  (``generalist_*``); a pattern that matched nothing is returned so the
  CLI can say so.
* ``QUALITY: <0-100>`` printed by an acceptance script is the graded
  part of the verdict. For a task with no expected signals it replaces
  the flat 15/20 routing share, and it is kept on the result
  (``verify_quality``) and through replicate aggregation (median).
* ``workspace_for(task_class="generalist")`` is the generalist fixture,
  and that fixture is restored after every attempt.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import benchmark as bm
from delfin.agent import benchmark_runner as br

_ROOT = Path(__file__).resolve().parents[1]


def _task(**kw) -> bm.Task:
    base = dict(id="generalist_x", task_class="generalist", mode="solo",
                prompt="Mach X.", verify="x.py")
    base.update(kw)
    return bm.Task(**base)


# --- prompt_en ------------------------------------------------------------

def test_a_task_file_carries_an_english_prompt(tmp_path):
    p = tmp_path / "t.yaml"
    p.write_text(
        "tasks:\n"
        "  - id: g1\n"
        "    task_class: generalist\n"
        "    prompt: Mach X.\n"
        "    prompt_en: Do X.\n", encoding="utf-8")
    (t,) = bm.load_tasks(p)
    assert t.prompt == "Mach X."
    assert t.prompt_en == "Do X."


def test_the_english_variant_has_its_own_id_and_prompt():
    t = _task(prompt_en="Do X.")
    en = bm.english_variant(t)
    assert en is not None
    assert en.id == "generalist_x.en"
    assert en.prompt == "Do X."
    # Everything else is the same task: same fixture, same verdict.
    assert (en.setup, en.verify, en.task_class, en.max_tool_calls) == (
        t.setup, t.verify, t.task_class, t.max_tool_calls)


def test_a_task_without_an_english_prompt_has_no_variant():
    assert bm.english_variant(_task()) is None


# --- select_tasks ---------------------------------------------------------

def test_a_glob_selects_a_family_and_an_exact_id_still_works():
    tasks = [_task(id="generalist_a"), _task(id="generalist_b"),
             _task(id="gen_other")]
    got, unmatched = bm.select_tasks(tasks, ["generalist_*", "gen_other"])
    assert [t.id for t in got] == ["generalist_a", "generalist_b",
                                   "gen_other"]
    assert unmatched == []


def test_a_pattern_that_matched_nothing_is_reported():
    got, unmatched = bm.select_tasks([_task(id="generalist_a")],
                                     ["generalist_a", "nosuch_*", "typo"])
    assert [t.id for t in got] == ["generalist_a"]
    assert unmatched == ["nosuch_*", "typo"]


def test_no_pattern_selects_everything():
    tasks = [_task(id="a"), _task(id="b")]
    assert bm.select_tasks(tasks, []) == (tasks, [])


def test_comma_separated_patterns_are_split():
    tasks = [_task(id="generalist_a"), _task(id="x")]
    got, _ = bm.select_tasks(tasks, ["generalist_*,x"])
    assert [t.id for t in got] == ["generalist_a", "x"]


def test_the_cli_accepts_a_language_flag():
    from delfin.agent import cli
    parser = cli.build_parser()
    ns = parser.parse_args(["bench", "run", "--model", "m",
                            "--task", "generalist_*", "--lang", "en"])
    assert ns.lang == "en"
    ns = parser.parse_args(["bench", "run", "--model", "m"])
    assert ns.lang == "de"


# --- QUALITY from the acceptance script -----------------------------------

def _traj(ok, out):
    tr = bm.Trajectory(text="done", duration_s=1.0)
    tr.verify_ok = ok
    tr.verify_output = out
    return tr


def test_the_acceptance_quality_is_kept_and_grades_the_routing_share():
    r = bm.score_outcome(_task(), _traj(True, "PASS a\nQUALITY: 40\nOK"))
    assert r.verify_quality == 40
    assert r.components["routing_pts"] == 8          # 20 * 40 / 100
    assert r.success is True


def test_a_failed_acceptance_keeps_its_partial_quality():
    r = bm.score_outcome(_task(), _traj(False, "FAIL a\nQUALITY: 75\nFAIL"))
    assert r.success is False
    assert r.verify_quality == 75
    assert r.components["routing_pts"] == 15


def test_without_a_quality_line_nothing_changes():
    r = bm.score_outcome(_task(), _traj(True, "OK"))
    assert r.verify_quality is None
    assert r.components["routing_pts"] == 15


def test_a_task_with_signals_keeps_the_signal_share():
    t = _task(expected_signals=(bm.Signal(pattern="done", against="text"),))
    r = bm.score_outcome(t, _traj(True, "QUALITY: 0\nOK"))
    assert r.verify_quality == 0
    assert r.components["routing_pts"] == 15


@pytest.mark.parametrize("line", ["QUALITY: 101", "QUALITY: -3",
                                  "QUALITY: abc", "quality 50"])
def test_a_malformed_quality_line_is_ignored(line):
    r = bm.score_outcome(_task(), _traj(True, line))
    assert r.verify_quality is None


def test_replicates_keep_the_median_acceptance_quality():
    rs = [bm.score_outcome(_task(), _traj(ok, f"QUALITY: {q}"))
          for ok, q in ((True, 100), (False, 20), (True, 80))]
    agg = bm.aggregate_replicates(rs)
    assert agg.verify_quality == 80
    assert agg.per_run_verify_quality == [100, 20, 80]


# --- the fixture workspace ------------------------------------------------

def test_the_generalist_class_writes_into_its_own_fixture():
    ws = br.workspace_for(_ROOT, task_class="generalist")
    assert ws == _ROOT / "tests" / "fixtures" / "generalist_workspace"


def test_the_generalist_fixture_is_restored_after_every_attempt():
    rel = Path("tests") / "fixtures" / "generalist_workspace"
    assert rel in br._BEHAVIOR_WS_RELS


# --- process signals --------------------------------------------------------

def _call(name, **inp):
    return {"name": name, "input": inp}


def _ptraj(calls, results):
    tr = bm.Trajectory(text="done", duration_s=1.0, tool_calls=calls,
                       tool_results=results)
    tr.verify_ok = True
    tr.verify_output = "QUALITY: 100\nOK"
    return tr


def test_a_run_that_planned_edited_and_tested_last_says_so():
    calls = [_call("task_create", subject="plan"),
             _call("read_file", file_path="textkit/src/textkit/slug.py"),
             _call("bash", command="cd textkit && python -m pytest -q"),
             _call("write_file", file_path="textkit/src/textkit/words.py"),
             _call("bash", command="cd textkit && python -m pytest -q")]
    results = ["ok", "text", '{"exit_code": 1, "stdout": "1 failed"}', "ok",
               '{"exit_code": 0, "stdout": "3 passed"}']
    p = bm.process_signals(_task(max_tool_calls=10), _ptraj(calls, results))
    assert p["planned"] == 1.0 and p["scouted"] == 1.0
    assert p["tested"] == 1.0 and p["tested_after_last_edit"] == 1.0
    assert (p["test_runs"], p["edits"], p["tool_calls"]) == (2.0, 1.0, 5.0)
    assert p["nonzero_exits"] == 1.0 and p["tool_errors"] == 0.0
    assert p["tool_budget_used"] == 0.5


def test_an_edit_after_the_last_test_run_is_not_a_tested_end_state():
    calls = [_call("bash", command="cargo test"),
             _call("edit_file", file_path="unitconv/src/lib.rs"),
             _call("bash", command="cat unitconv/src/lib.rs")]
    results = ['{"exit_code": 0}', '{"error": "old_string not found"}', "x"]
    p = bm.process_signals(_task(), _ptraj(calls, results))
    assert p["tested"] == 1.0
    assert p["tested_after_last_edit"] == 0.0
    assert p["planned"] == 0.0
    assert p["tool_errors"] == 1.0


@pytest.mark.parametrize("cmd", [
    "npm test", "npm run test", "node --no-warnings --test test/a.test.ts",
    "go test ./...", "cd stats && cargo test --offline",
    "python3 -m pytest tests/", "pytest -q"])
def test_every_toolchain_test_command_counts(cmd):
    p = bm.process_signals(_task(), _ptraj([_call("bash", command=cmd)],
                                           ['{"exit_code": 0}']))
    assert p["test_runs"] == 1.0, cmd


def test_other_classes_get_no_process_signals():
    t = _task(task_class="generic_project")
    assert bm.process_signals(t, _ptraj([], [])) == {}


def test_process_signals_reach_the_result_and_average_over_replicates():
    t = _task()
    a = bm.score_outcome(t, _ptraj([_call("bash", command="pytest")],
                                   ['{"exit_code": 0}']))
    b = bm.score_outcome(t, _ptraj([], []))
    assert a.process["tested"] == 1.0 and b.process["tested"] == 0.0
    agg = bm.aggregate_replicates([a, b])
    assert agg.process["tested"] == 0.5


# --- toolchains -----------------------------------------------------------

def _fake_acceptance(monkeypatch, tmp_path, body: str) -> None:
    script = tmp_path / "acc.py"
    script.write_text(body, encoding="utf-8")
    monkeypatch.setattr(br, "acceptance_path", lambda name, root=None: script)


def test_a_missing_toolchain_is_not_run_not_failed(monkeypatch, tmp_path):
    _fake_acceptance(monkeypatch, tmp_path,
                     "import sys\nprint('toolchain unavailable: cargo',"
                     " file=sys.stderr)\nsys.exit(3)\n")
    ok, out = br.run_acceptance("x.py", tmp_path, root=tmp_path)
    assert ok is None
    assert "toolchain unavailable" in out


def test_exit_three_without_the_phrase_is_an_ordinary_failure(monkeypatch,
                                                               tmp_path):
    _fake_acceptance(monkeypatch, tmp_path, "import sys\nsys.exit(3)\n")
    ok, _ = br.run_acceptance("x.py", tmp_path, root=tmp_path)
    assert ok is False


def test_the_toolchain_is_first_on_path_for_a_generalist_attempt(
        monkeypatch, tmp_path):
    tc = tmp_path / "tc-bin"
    tc.mkdir()
    monkeypatch.setenv("DELFIN_TOOLCHAIN_BIN", str(tc))
    monkeypatch.setenv("PATH", "/usr/bin")
    with br._toolchain_on_path("generalist"):
        assert os.environ["PATH"].split(os.pathsep)[0] == str(tc)
    assert os.environ["PATH"] == "/usr/bin"
    with br._toolchain_on_path("generic_project"):
        assert os.environ["PATH"] == "/usr/bin"


def test_no_toolchain_directory_changes_nothing(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_TOOLCHAIN_BIN", str(tmp_path / "absent"))
    monkeypatch.setenv("PATH", "/usr/bin")
    assert br.toolchain_bin() is None
    with br._toolchain_on_path("generalist"):
        assert os.environ["PATH"] == "/usr/bin"


@pytest.mark.parametrize("setup", [
    "generalist_rust_feature.py", "generalist_rust_red_to_green.py",
    "generalist_go_feature.py", "generalist_go_red_to_green.py",
    "generalist_ts_feature.py"])
def test_a_setup_without_its_toolchain_says_so_and_builds_nothing(
        setup, tmp_path):
    """Run on every host: PATH and the toolchain directory are emptied,
    so the setup must take the not-run route whatever is installed."""
    from conftest import child_env
    env = child_env(tmp_path)
    empty = tmp_path / "empty-bin"
    empty.mkdir()
    env["PATH"] = str(empty)
    env["DELFIN_TOOLCHAIN_BIN"] = str(empty)
    ws = tmp_path / "ws"
    ws.mkdir()
    proc = subprocess.run(
        [sys.executable, str(br.setup_path(setup, _ROOT)), str(ws)],
        capture_output=True, text=True, timeout=60, env=env, cwd=str(ws))
    assert proc.returncode == 3, proc.stderr
    assert "toolchain unavailable" in proc.stderr
    assert not any(ws.iterdir())
