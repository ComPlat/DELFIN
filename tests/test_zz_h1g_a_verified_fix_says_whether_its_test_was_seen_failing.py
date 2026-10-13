"""A verified test task says whether its own test was ever seen failing.

A test that passes with the fix and would also pass without it proves
nothing about the fix. The completion check already says "verified" when a
fresh green run exists; this adds one sentence to that verdict when the
session's durable tool trace never shows the changed test failing first.

Everything here goes through the real path with the real writers: the file
changes are recorded by ``change_journal.record_change``, the test runs by
``tool_trace.record`` (which stores the tool input as TEXT, not as a dict,
and caps the output), and the verdict is read from the document executor's
``_check_completion`` -- the function that runs when the model completes a
task. No hand-made ledger shapes.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent import agent_tasks, change_journal, tool_trace
from delfin.agent.api_client import KitToolPermissions, _doc_executor

SID = "sess-h1g-end-to-end"
GATE = "/repo/.delfin/gate-tools/gate tests/test_thing.py -q"
RED = ("tests/test_thing.py F.\nFAILED tests/test_thing.py::test_it - "
       "assert 1 == 2\n1 failed, 1 passed in 0.10s")
GREEN = "tests/test_thing.py ..\n2 passed in 0.10s"


@pytest.fixture
def ws(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    w = tmp_path / "ws"
    (w / "tests").mkdir(parents=True)
    (w / "pkg").mkdir()
    agent_tasks._STORES.clear()
    return w


def _write(ws: Path, rel: str) -> None:
    p = ws / rel
    p.write_text("x = 1\n", encoding="utf-8")
    assert change_journal.record_change(
        SID, tool="write_file", path=str(p), old_text=None,
        new_text="x = 1\n") is not None


def _ran(command: str, output: str, ok: bool,
         tool: str = "mcp__kit-coding__bash") -> None:
    tool_trace.record(SID, tool=tool,
                      tool_input={"command": command, "description": "run"},
                      output=output, ok=ok)
    time.sleep(0.02)


def _verdict(ws: Path) -> dict:
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits",
                               task_session_id=SID)
    out = json.loads(_doc_executor.execute(
        "task_create", {"subject": "Run the tests for the fix"}, perms))
    tid = out["task"]["id"]
    _doc_executor.execute(
        "task_update", {"task_id": tid, "status": "in_progress"}, perms)
    live = [{"tool": "bash", "command": GATE, "exit_code": 0, "status": "ok",
             "passed": 2, "failed": 0, "ts": time.time()}]
    perms.evidence_provider = lambda: {"observed": set(), "tests": live}
    res = _doc_executor._check_completion(tid, perms)
    assert isinstance(res, dict), res
    return res


def _fix(ws: Path) -> None:
    _write(ws, "pkg/thing.py")
    _write(ws, "tests/test_thing.py")


def test_the_trace_writer_stores_the_input_as_text(ws):
    """The shape this whole file stands on: if the writer ever stores a
    dict again, the tests below would silently test another shape."""
    _ran(GATE, RED, False)
    entry = tool_trace.read(SID)[-1]
    assert isinstance(entry["input"], str) and "gate-tools/gate" in entry["input"]


def test_red_then_green_through_the_gate_tool_adds_nothing(ws):
    _fix(ws)
    _ran(GATE, RED, False)
    _ran(GATE, GREEN, True)
    res = _verdict(ws)
    assert res["verdict"] == "verified", res
    assert "no control run" not in (res.get("note") or ""), res


def test_green_only_says_the_test_was_never_seen_failing(ws):
    _fix(ws)
    _ran(GATE, GREEN, True)
    res = _verdict(ws)
    assert res["verdict"] == "verified", res
    assert "never seen failing" in (res.get("note") or ""), res
    assert "tests/test_thing.py" in res["note"]


def test_a_changed_test_that_never_ran_is_named(ws):
    _fix(ws)
    res = _verdict(ws)
    assert "was changed but never run" in (res.get("note") or ""), res


def test_green_first_and_red_afterwards_is_not_a_control(ws):
    _fix(ws)
    node = "/repo/.delfin/gate-tools/gate tests/test_thing.py::test_it -q"
    _ran(node, "1 passed in 0.1s", True)
    _ran(node, "FAILED tests/test_thing.py::test_it\n1 failed in 0.1s", False)
    note = _verdict(ws).get("note") or ""
    assert "no control run" in note, note


def test_a_red_run_on_another_file_is_not_this_fixs_control(ws):
    _fix(ws)
    _ran("/repo/.delfin/gate-tools/gate tests/test_other.py -q",
         "FAILED tests/test_other.py::test_b\n1 failed in 0.1s", False)
    _ran(GATE, GREEN, True)
    assert "never seen failing" in (_verdict(ws).get("note") or "")


def test_a_red_run_whose_summary_was_cut_off_still_counts_as_red(ws):
    """The trace keeps 2000 characters of output; a long failing run loses
    its closing '1 failed' line and is red by its exit status alone."""
    _fix(ws)
    _ran(GATE, "tests/test_thing.py F.\n" + "E   detail\n" * 400, False)
    _ran(GATE, GREEN, True)
    assert "no control run" not in (_verdict(ws).get("note") or "")


def test_the_sentence_never_changes_the_verdict(ws):
    """It is information for the reader, not a second gate: a task with a
    fresh green run stays verified whatever the sentence says."""
    _fix(ws)
    _ran(GATE, GREEN, True)
    assert _verdict(ws)["verdict"] == "verified"


def test_a_documentation_change_gets_no_sentence(ws):
    _write(ws, "README.md")
    _ran(GATE, GREEN, True)
    assert "no control run" not in (_verdict(ws).get("note") or "")


# -- the sentence on the model's side, and run shapes the trace really holds --

def _complete(ws: Path) -> str:
    """The tool result the model receives for task_update(completed)."""
    perms = KitToolPermissions(workspace=ws, mode="acceptEdits",
                               task_session_id=SID)
    tid = json.loads(_doc_executor.execute(
        "task_create", {"subject": "Run the tests for the fix"},
        perms))["task"]["id"]
    _doc_executor.execute(
        "task_update", {"task_id": tid, "status": "in_progress"}, perms)
    live = [{"tool": "bash", "command": GATE, "exit_code": 0, "status": "ok",
             "passed": 2, "failed": 0, "ts": time.time()}]
    perms.evidence_provider = lambda: {"observed": set(), "tests": live}
    return _doc_executor.execute(
        "task_update", {"task_id": tid, "status": "completed"}, perms)


def test_the_sentence_reaches_the_model_in_the_task_update_result(ws):
    _fix(ws)
    _ran(GATE, GREEN, True)
    out = json.loads(_complete(ws))
    assert out["task"]["verified"] == "verified", out
    assert "never seen failing" in out.get("note", ""), out
    assert "never seen failing" in out["task"]["verify_note"], out


def test_a_controlled_fix_adds_no_note_to_the_task_update_result(ws):
    _fix(ws)
    _ran(GATE, RED, False)
    _ran(GATE, GREEN, True)
    out = json.loads(_complete(ws))
    assert "note" not in out, out


def _run_tests(target: str, *, passed: int, failed: int, errors: int = 0,
               nodes=()) -> None:
    status = "failed" if (failed or errors) else "ok"
    tool_trace.record(
        SID, tool="mcp__kit-coding__run_tests",
        tool_input={"target": target},
        output=json.dumps({
            "status": status, "framework": "pytest",
            "exit_code": 1 if status == "failed" else 0,
            "summary": {"passed": passed, "failed": failed,
                        "errors": errors},
            "failures": [{"node_id": n, "message": "assert 1 == 2"}
                         for n in nodes]}),
        ok=True)
    time.sleep(0.02)


def test_the_run_tests_tool_in_its_recorded_shape_counts(ws):
    """run_tests is recorded with its own arguments ({"target": ...}) and a
    JSON result; it has no "command" and no pytest console line."""
    _fix(ws)
    _run_tests("tests/test_thing.py", passed=1, failed=1,
               nodes=["tests/test_thing.py::test_it"])
    _run_tests("tests/test_thing.py", passed=2, failed=0)
    assert "no control run" not in (_verdict(ws).get("note") or "")


def test_the_run_tests_tool_green_only_is_never_seen_failing(ws):
    _fix(ws)
    _run_tests("tests/test_thing.py", passed=2, failed=0)
    assert "never seen failing" in (_verdict(ws).get("note") or "")


def test_a_command_that_only_prints_a_summary_is_not_a_test_run(ws):
    """No runner token: an echo or a cat of a file that contains pytest
    words is not a run, whatever it prints."""
    _fix(ws)
    _ran("echo FAILED tests/test_thing.py::test_it 1 failed in"
         " tests/test_thing.py",
         "FAILED tests/test_thing.py::test_it\n1 failed", False)
    _ran("cat tests/test_thing.py", "# 2 passed\n", True)
    _ran(GATE, GREEN, True)
    assert "never seen failing" in (_verdict(ws).get("note") or "")


def test_a_collection_error_is_not_the_test_failing(ws):
    """An ImportError at collection means the test never ran; the red is
    not the test's assertion failing without the fix."""
    _fix(ws)
    _ran(GATE, "ERROR tests/test_thing.py - ImportError: cannot import "
               "name 'helper'\n!!! Interrupted: 1 error during collection"
               " !!!\n1 error in 0.20s", False)
    _ran(GATE, GREEN, True)
    assert "never seen failing" in (_verdict(ws).get("note") or "")


def test_a_green_run_of_another_selection_does_not_control_a_red_node(ws):
    """-k selects; a green -k test_other run says nothing about test_it."""
    _fix(ws)
    _ran(GATE, RED, False)
    _ran(GATE + " -k test_other", "1 passed, 1 deselected in 0.1s", True)
    assert "no control run" in (_verdict(ws).get("note") or "")


def test_the_same_selection_red_then_green_is_a_control(ws):
    _fix(ws)
    _ran(GATE + " -k test_it", RED, False)
    _ran(GATE + " -k test_it", "1 passed, 1 deselected in 0.1s", True)
    assert "no control run" not in (_verdict(ws).get("note") or "")


def test_a_flag_before_the_target_does_not_hide_the_target(ws):
    """-q takes no value; the token after it is the test file."""
    _fix(ws)
    _ran("python -m pytest -q tests/test_thing.py", GREEN, True)
    note = _verdict(ws).get("note") or ""
    assert "never seen failing" in note, note
    assert "never run" not in note, note


def test_a_test_edited_after_its_red_run_is_not_a_control(ws):
    """The test that went green is not the one seen red when it changed in
    between. The journal stamps whole seconds, hence the waits."""
    _fix(ws)
    time.sleep(1.1)
    _ran(GATE, RED, False)
    time.sleep(1.1)
    _write(ws, "tests/test_thing.py")
    _ran(GATE, GREEN, True)
    note = _verdict(ws).get("note") or ""
    assert "changed after it was seen failing" in note, note
