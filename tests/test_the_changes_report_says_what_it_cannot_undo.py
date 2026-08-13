""""What did you change?" and "undo that" have to answer the same list.

The report was built from the audit log alone, and undoability was
implied by the tool's NAME. That was wrong in both directions:

* A shell command, a ``fill_series`` run and every MCP file mutator
  (audited under the BARE native tool name, indistinguishable from a
  native write) appeared under "Files written" exactly like a journalled
  write. The user read the list, asked for an undo, and got a subset —
  with no statement of what was never covered.
* The record's ``path`` was read from the ``path`` argument only, while
  the write tools also accept ``file_path``/``filename``/``file``/
  ``target`` — the aliases that exist for the open-weights models this
  project targets. A ``write_file(file_path=…)`` logged an EMPTY path,
  and the report answered "No recorded changes" for a file that had been
  written and journalled.

Two more the same records got wrong: a command that ran and FAILED was
logged ``ok`` because the decision described the gate's verdict on the
attempt, and the report deliberately ignored rotated history, so on the
first day of a new ISO week it reported "No recorded changes" without
saying that rotation had happened.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin.agent import audit_log as al
from delfin.agent import change_journal as cj
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    return tmp_path / "home"


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


def _log(home: Path) -> Path:
    p = home / ".delfin" / "audit.log"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def _perms(ws: Path, sid: str) -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode="default")
    p.task_session_id = sid
    return p


# ---------------------------------------------------------------------------
# The path the executor used, whatever the model called the argument
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("key", ["path", "file_path", "filename", "file", "target"])
def test_a_written_file_is_reported_whatever_the_argument_was_called(
        home, ws, key):
    ex = _DocToolExecutor()
    perms = _perms(ws, "s-alias")
    args = {key: "note.txt", "content": "hello\n"}
    result = ex._execute_write_file(args, perms)
    assert "File created" in result, result
    ex._audit_call("write_file", args, perms, result)

    report = al.build_changes_report("s-alias", log_path=_log(home))
    paths = [f["path"] for f in report["files_written"]]
    assert paths == [str(ws / "note.txt")], report
    assert "No recorded changes" not in al.format_changes_report(report)


def test_the_reported_path_lines_up_with_the_journal(home, ws):
    ex = _DocToolExecutor()
    perms = _perms(ws, "s-match")
    args = {"file_path": "note.txt", "content": "hello\n"}
    ex._audit_call("write_file", args, perms,
                   ex._execute_write_file(args, perms))

    report = al.build_changes_report("s-match", log_path=_log(home))
    assert report["files_written"][0]["undoable"] is True
    assert report["unjournalled"] == []


# ---------------------------------------------------------------------------
# Undoable comes from the journal, not from the tool's name
# ---------------------------------------------------------------------------

def test_a_write_with_no_pre_image_is_marked_not_undoable(home, ws):
    """An MCP write is audited under the bare native tool name."""
    al.append(al.make_record(
        tool="write_file", decision="ok", session_id="s-mcp",
        path=str(ws / "remote.txt")), log_path=_log(home))

    report = al.build_changes_report("s-mcp", log_path=_log(home))
    assert report["files_written"][0]["undoable"] is False
    assert report["unjournalled"] == [str(ws / "remote.txt")]
    text = al.format_changes_report(report)
    assert "NOT undoable" in text
    assert "Not undoable (1)" in text


def test_the_journalled_and_the_unjournalled_are_told_apart(home, ws):
    journalled = ws / "native.txt"
    journalled.write_text("agent\n", encoding="utf-8")
    cj.record_change("s-mix", tool="write_file", path=journalled,
                     old_text="user\n", new_text="agent\n")
    for p in (journalled, ws / "by_shell.txt"):
        al.append(al.make_record(
            tool="write_file", decision="ok", session_id="s-mix",
            path=str(p)), log_path=_log(home))

    report = al.build_changes_report("s-mix", log_path=_log(home))
    marks = {f["path"]: f["undoable"] for f in report["files_written"]}
    assert marks[str(journalled)] is True
    assert marks[str(ws / "by_shell.txt")] is False


def test_an_entry_without_a_restorable_pre_image_is_not_promised(home, ws):
    p = ws / "big.txt"
    p.write_text("small\n", encoding="utf-8")
    cj.record_change("s-trunc", tool="write_file", path=p,
                     old_text="x" * (cj.MAX_PRE_IMAGE_BYTES + 10),
                     new_text="small\n")
    al.append(al.make_record(
        tool="write_file", decision="ok", session_id="s-trunc",
        path=str(p)), log_path=_log(home))

    report = al.build_changes_report("s-trunc", log_path=_log(home))
    assert report["files_written"][0]["undoable"] is False, (
        "a truncated pre-image is a record, not a restorable file")


def test_without_a_session_undoability_is_unknown_not_promised(home, ws):
    al.append(al.make_record(
        tool="write_file", decision="ok", session_id="",
        path=str(ws / "x.txt")), log_path=_log(home))
    report = al.build_changes_report(log_path=_log(home))
    assert report["files_written"][0]["undoable"] is None
    assert "NOT undoable" not in al.format_changes_report(report)


def test_every_write_tool_either_journals_or_is_declared_unjournalled():
    """The gap has to be a declaration somebody maintains.

    A write route that reaches the user's files without a pre-image is
    allowed to exist — ``bash`` cannot be journalled — but it may not be
    presented as undoable, and adding a new one must not be able to slip
    through silently.
    """
    src = (Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    captured = set(re.findall(
        r'_capture(?:_binary|_raw)?_change\(\s*"([a-z_]+)"', src))
    # undo_changes journals through change_journal.revert itself.
    captured.add("undo_changes")

    tools = set(al._WRITE_TOOLS) | set(al._COMMAND_TOOLS)
    undeclared = sorted(
        t for t in tools
        if t not in captured and t not in al._UNJOURNALLED_WRITE_TOOLS)
    assert not undeclared, (
        "these write tools neither journal a pre-image nor are declared "
        f"unjournalled: {undeclared}")

    # And the declaration must not claim a tool that does journal.
    wrongly_declared = sorted(al._UNJOURNALLED_WRITE_TOOLS & captured)
    assert not wrongly_declared, wrongly_declared


def test_the_native_writers_really_do_journal(home, ws):
    """The declaration above is only worth what the executor does."""
    ex = _DocToolExecutor()
    perms = _perms(ws, "s-real")
    (ws / "e.txt").write_text("one\ntwo\n", encoding="utf-8")
    perms.read_tracker[str(ws / "e.txt")] = (ws / "e.txt").stat().st_mtime

    ex._execute_write_file({"path": "w.txt", "content": "x\n"}, perms)
    ex._execute_edit_file(
        {"path": "e.txt", "old_string": "two", "new_string": "TWO"}, perms)
    perms.read_tracker[str(ws / "e.txt")] = (ws / "e.txt").stat().st_mtime
    ex._execute_multi_edit(
        {"path": "e.txt", "edits": [{"old_string": "one", "new_string": "ONE"}]},
        perms)

    tools = [r["tool"] for r in cj.list_changes("s-real")]
    assert tools == ["write_file", "edit_file", "multi_edit"], tools


# ---------------------------------------------------------------------------
# Attempts versus outcomes
# ---------------------------------------------------------------------------

def test_a_command_that_ran_and_failed_is_not_logged_ok(home, ws):
    ex = _DocToolExecutor()
    perms = _perms(ws, "s-exit")
    result = ex._execute_bash({"command": "exit 7"}, perms)
    ex._audit_call("bash", {"command": "exit 7"}, perms, result)

    rec = al.read_last_n(1, log_path=_log(home))[0]
    assert rec["decision"] == "error"
    assert rec["exit_code"] == 7
    text = al.format_changes_report(
        al.build_changes_report("s-exit", log_path=_log(home)))
    assert "[exit 7]" in text
    assert "exit 7" in text


def test_a_command_that_succeeded_is_still_logged_ok(home, ws):
    ex = _DocToolExecutor()
    perms = _perms(ws, "s-ok")
    result = ex._execute_bash({"command": "true"}, perms)
    ex._audit_call("bash", {"command": "true"}, perms, result)

    rec = al.read_last_n(1, log_path=_log(home))[0]
    assert rec["decision"] == "ok"
    assert rec["exit_code"] == 0
    assert "[exit" not in al.format_changes_report(
        al.build_changes_report("s-ok", log_path=_log(home)))


def test_a_background_job_gets_a_completion_record(home, ws):
    """The launch says a job started; nothing said how it ended."""
    ex = _DocToolExecutor()
    perms = _perms(ws, "s-bg")
    started = json.loads(
        ex._execute_bash_background({"command": "exit 3"}, perms))
    ex._audit_call("bash_background", {"command": "exit 3"}, perms,
                   json.dumps(started))
    job_id = started["job_id"]

    import time
    for _ in range(200):
        payload = json.loads(ex._execute_bash_output({"job_id": job_id}))
        ex._audit_call("bash_output", {"job_id": job_id}, perms,
                       json.dumps(payload))
        if payload.get("running") is False:
            break
        time.sleep(0.02)

    records = al.read_last_n(50, log_path=_log(home))
    done = [r for r in records if r.get("event") == "completed"]
    assert done, "a background job still finishes without a record"
    assert done[-1]["exit_code"] == 3
    assert done[-1]["decision"] == "error"
    assert len(done) == 1, "polling must not append a record per poll"


# ---------------------------------------------------------------------------
# Rotated history
# ---------------------------------------------------------------------------

def test_the_first_day_of_a_new_week_still_finds_last_weeks_changes(home, ws):
    log = _log(home)
    rotated = log.with_name("audit-2026-W31.log")
    rotated.write_text(json.dumps({
        "ts": "2026-08-01T09:00:00Z", "session_id": "s-old",
        "tool": "write_file", "decision": "ok", "path": str(ws / "old.txt"),
    }) + "\n", encoding="utf-8")

    report = al.build_changes_report("s-old", log_path=log)
    assert report["files_written"], "a whole week of history was invisible"
    assert report["window"]["rotated"] == "audit-2026-W31.log"
    assert "incl. rotated audit-2026-W31.log" in al.format_changes_report(report)


def test_an_empty_answer_says_whether_history_was_consulted(home, ws):
    text = al.format_changes_report(
        al.build_changes_report("s-none", log_path=_log(home)))
    assert "No recorded changes" in text
    assert "rotated" in text.lower()


def test_the_active_log_alone_is_enough_when_it_answers(home, ws):
    log = _log(home)
    log.with_name("audit-2026-W31.log").write_text(json.dumps({
        "ts": "2026-08-01T09:00:00Z", "session_id": "s-cur",
        "tool": "write_file", "decision": "ok", "path": str(ws / "old.txt"),
    }) + "\n", encoding="utf-8")
    al.append(al.make_record(
        tool="write_file", decision="ok", session_id="s-cur",
        path=str(ws / "new.txt")), log_path=log)

    report = al.build_changes_report("s-cur", log_path=log)
    assert [f["path"] for f in report["files_written"]] == [str(ws / "new.txt")]
    assert report["window"]["rotated"] == ""
