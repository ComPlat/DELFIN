"""list_changes_made answered with another agent's writes.

The tool exists so "what did you change?" is answered from the record
rather than from memory -- its own description says so. Driven on
2026-09-09 with no session id attached, it answered with the newest 200
audit records from EVERY session: another agent's writes in the same
shared checkout, two benchmark worktrees, four unrelated temp
directories, and a notebook edited by a probe half an hour earlier.

``build_changes_report(None)`` means "no filter", and
``_execute_list_changes`` passed None whenever ``task_session_id`` was
empty. An agent reading that output would tell the user, citing its audit
log, that it had edited files it never opened. A fabricated claim with a
tool call behind it is worse than one without: the grounding rules are
satisfied, and the answer is still false.

Two changes, because one is not enough:

* the workspace is passed as a filter, which removes the absolute strays;
* when there is no session id the report SAYS so, in its first line.

The second is not belt-and-braces. The workspace filter cannot close the
gap alone: file tools record workspace-RELATIVE paths, so a relative path
from a different workspace is indistinguishable from one of ours, and
``_under_workspace`` lets it through by design. A caveat is the honest
answer where a filter cannot be one -- and it goes first, because under a
long list it is a line nobody reaches.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def ws():
    with tempfile.TemporaryDirectory(prefix="chg-") as tmp:
        yield Path(tmp)


def _perms(ws, sid=""):
    p = A.KitToolPermissions(mode="bypassPermissions", workspace=str(ws))
    p.task_session_id = sid
    return p


def _report(perms) -> str:
    return str(A._doc_executor.execute("list_changes_made", {}, perms))


def test_a_session_scoped_report_holds_only_this_sessions_work(ws):
    import uuid
    sid = uuid.uuid4().hex
    perms = _perms(ws, sid)
    A._doc_executor.execute(
        "write_file", {"path": "mine.py", "content": "x = 1\n"}, perms)

    other = _perms(ws, uuid.uuid4().hex)
    A._doc_executor.execute(
        "write_file", {"path": "theirs.py", "content": "y = 2\n"}, other)

    out = _report(perms)
    assert "mine.py" in out
    assert "theirs.py" not in out
    assert "NOT SESSION-SCOPED" not in out


def test_without_a_session_the_report_says_it_is_not_scoped(ws):
    perms = _perms(ws)
    A._doc_executor.execute(
        "write_file", {"path": "a.py", "content": "x = 1\n"}, perms)
    out = _report(perms)
    assert out.startswith("NOT SESSION-SCOPED")
    assert "not necessarily your own work" in out
    assert "unless you can point at the call that did" in out


def test_the_caveat_comes_before_the_list_not_after(ws):
    perms = _perms(ws)
    A._doc_executor.execute(
        "write_file", {"path": "a.py", "content": "x = 1\n"}, perms)
    out = _report(perms)
    if "### Changes made" in out:
        assert out.index("NOT SESSION-SCOPED") < out.index("### Changes made")


def test_an_absolute_path_from_another_workspace_is_filtered(ws, tmp_path):
    """The half a filter CAN do."""
    from delfin.agent import audit_log as _al
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    report = _al.build_changes_report(
        None, workspace=str(ws),
        log_path=_write_log(tmp_path, [
            {"ts": "2026-09-09T10:00:00Z", "session_id": "s1",
             "tool": "write_file", "decision": "ok",
             "path": str(elsewhere / "theirs.py")},
            {"ts": "2026-09-09T10:00:01Z", "session_id": "s1",
             "tool": "write_file", "decision": "ok",
             "path": str(ws / "mine.py")},
        ]))
    paths = [f["path"] for f in report["files_written"]]
    assert str(ws / "mine.py") in paths
    assert str(elsewhere / "theirs.py") not in paths


def _write_log(tmp_path: Path, records: list[dict]) -> Path:
    import json
    p = tmp_path / "audit.log"
    p.write_text("\n".join(json.dumps(r) for r in records) + "\n",
                 encoding="utf-8")
    return p


def test_the_tool_still_answers_rather_than_refusing(ws):
    """Refusing outright was the other option and is the wrong one: the
    engine sets a session id in normal use, and where it does not, an
    unscoped list with an honest label still helps more than nothing."""
    perms = _perms(ws)
    A._doc_executor.execute(
        "write_file", {"path": "a.py", "content": "x = 1\n"}, perms)
    out = _report(perms)
    assert "a.py" in out
    assert "error" not in out.lower()[:40]
