"""One look says which sessions are open and what each is waiting on.

``delfin-agent sessions`` listed what ran BEFORE, which is misleading
when five are running now: the supervisor reads a table of history and
concludes nothing is live.

Assembling the live picture on 2026-09-20 took `tmux capture-pane` per
pane, a `git log` per worktree and a grep through the audit log -- and
the pane capture only says a session is waiting, never on what, because
the dialog is cut to 24 lines and to the pane's width.

Both halves exist now: a terminal session announces its presence, and a
question waiting at one is published whole. This joins them.

  who is open        from presence, so a dead session is not listed
  what it waits on   from the published question, by session key
  where it works     the branch, because five worktrees look alike
"""

from __future__ import annotations

import pytest

from delfin.agent import cli_resume as CR
from delfin.agent import session_presence as pres
from delfin.agent import terminal_confirm as tc


@pytest.fixture
def world(tmp_path, monkeypatch):
    monkeypatch.setattr(pres, "_DIR", tmp_path / "presence")
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(pres, "_last_written", {})
    return tmp_path


def _open(key, workspace, title=""):
    pres.announce(key, session_id=f"id-{key}", title=title,
                  workspace=str(workspace))


def _waiting(key, tool="bash", preview="x"):
    broker = tc.TerminalConfirmBroker(session_id=f"id-{key}", session_key=key)
    return broker, broker._enqueue(
        tc.ConfirmRequest(kind=tc.CONFIRM, tool=tool,
                          args={"command": "c"}, preview=preview))


class TestTheOverview:
    def test_an_open_session_is_listed(self, world):
        _open("runde2-s1", world)
        rows = CR.open_now()
        assert [r["key"] for r in rows] == ["runde2-s1"]

    def test_nothing_open_is_an_empty_list(self, world):
        assert CR.open_now() == []

    def test_what_it_waits_on_is_joined_in(self, world):
        _open("runde2-s1", world)
        _waiting("runde2-s1", tool="edit_file")
        row = CR.open_now()[0]
        assert row["waiting"]["tool"] == "edit_file"

    def test_a_session_waiting_on_nothing_says_so(self, world):
        _open("runde2-s2", world)
        assert CR.open_now()[0]["waiting"] is None

    def test_the_question_goes_to_the_right_session(self, world):
        _open("runde2-s1", world)
        _open("runde2-s2", world)
        _waiting("runde2-s2", tool="write_file")
        by_key = {r["key"]: r["waiting"] for r in CR.open_now()}
        assert by_key["runde2-s1"] is None
        assert by_key["runde2-s2"]["tool"] == "write_file"

    def test_a_protected_question_is_marked(self, world):
        _open("runde2-s1", world)
        _waiting("runde2-s1", preview="[SELF-MODIFICATION GUARD]\nd")
        assert CR.open_now()[0]["waiting"]["protected"] is True

    def test_an_answered_question_leaves_the_row_clean(self, world):
        _open("runde2-s1", world)
        broker, req = _waiting("runde2-s1")
        broker.resolve(req, True)
        assert CR.open_now()[0]["waiting"] is None


class TestTheRendering:
    def test_it_names_the_key_and_the_wait(self, world):
        _open("runde2-s1", world)
        _waiting("runde2-s1", tool="edit_file",
                 preview="[SELF-MODIFICATION GUARD]\nd")
        text = CR.render_open(CR.open_now())
        assert "runde2-s1" in text
        assert "edit_file" in text
        assert "PROTECTED" in text

    def test_a_quiet_session_renders_without_a_wait(self, world):
        _open("runde2-s2", world)
        text = CR.render_open(CR.open_now())
        assert "runde2-s2" in text
        assert "PROTECTED" not in text

    def test_nothing_open_renders_as_nothing(self, world):
        assert CR.render_open([]) == ""

    def test_the_overview_never_raises(self, monkeypatch):
        monkeypatch.setattr(pres, "open_sessions",
                            lambda **k: (_ for _ in ()).throw(OSError("gone")))
        assert CR.open_now() == []
