"""A terminal session can be reached without typing into its terminal.

``session_messages.send`` wrote a file that nobody read: ``take`` is
called in the dashboard and nowhere else. A ``delfin-agent chat`` session
announced no presence either, so it did not appear in ``open_sessions``
and the ``session_message`` tool answered "only inside an open dashboard
session".

The only way in was the terminal itself, and that way has a trap: an
approval dialog reads single keys, so a message pasted while one is up
does not queue -- its first character ANSWERS the dialog. Supervising
five sessions on 2026-09-20 meant checking every pane for an open dialog
before typing, and the operator's recall channel, which the written
instructions prescribe, reached nothing at all.

The delivery point is the one the job wake-up already uses: the idle
prompt, polled between reads. Never during a turn, never while a dialog
holds the reader -- so the trap cannot arise -- and never while something
is typed, because a draft belongs to whoever typed it.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_messages as msgs


class _Theme:
    def dim(self, t): return t
    def bold(self, t): return t
    def red(self, t): return t


class _Transcript:
    def __init__(self): self.lines = []
    theme = _Theme()
    def chrome(self, line): self.lines.append(line)


class _Opts:
    def __init__(self, name=""): self.session_name = name


def _agent(name="runde2-s1", session_id="abcdef0123456789"):
    from delfin.agent.repl import TerminalAgent
    a = object.__new__(TerminalAgent)
    a.opts = _Opts(name)
    a.engine = type("E", (), {"session_id": session_id})()
    a.transcript = _Transcript()
    return a


@pytest.fixture
def inbox(tmp_path, monkeypatch):
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "session_inbox")
    return tmp_path


class TestTheKey:
    def test_the_name_is_the_address(self):
        assert _agent(name="runde2-s1")._presence_key() == "runde2-s1"

    def test_without_a_name_the_session_id_serves(self):
        assert _agent(name="")._presence_key() == "abcdef01"


class TestDelivery:
    def test_a_waiting_message_reaches_the_prompt(self, inbox):
        msgs.send("runde2-s1", "stop searching in the home directory")
        got = _agent()._operator_messages("")
        assert "stop searching in the home directory" in got

    def test_an_empty_inbox_is_silent(self, inbox):
        assert _agent()._operator_messages("") == ""

    def test_a_draft_holds_the_message_back(self, inbox):
        msgs.send("runde2-s1", "later")
        a = _agent()
        assert a._operator_messages("half a thought") == ""
        assert a._operator_messages("") != "", "it waits, it is not dropped"

    def test_a_message_is_taken_only_once(self, inbox):
        msgs.send("runde2-s1", "once")
        a = _agent()
        assert a._operator_messages("") != ""
        assert a._operator_messages("") == ""

    def test_the_sender_and_the_text_are_shown_on_the_transcript(self, inbox):
        # The message becomes the next prompt, and a prompt is not
        # echoed: without this line nobody watching the pane can see
        # what was delivered.
        msgs.send("runde2-s1", "stop searching in the home directory",
                  from_title="the operator")
        a = _agent()
        a._operator_messages("")
        shown = " ".join(a.transcript.lines)
        assert "operator" in shown
        assert "stop searching in the home directory" in shown

    def test_a_session_without_a_key_takes_nothing(self, inbox):
        msgs.send("", "nowhere")
        a = _agent(name="", session_id="")
        assert a._operator_messages("") == ""

    def test_a_broken_inbox_never_takes_the_prompt_with_it(self, inbox, monkeypatch):
        def boom(_key): raise OSError("inbox is on a dead mount")
        monkeypatch.setattr(msgs, "take", boom)
        assert _agent()._operator_messages("") == ""


class TestPresence:
    def test_the_session_announces_itself(self, tmp_path, monkeypatch):
        from delfin.agent import session_presence as pres
        monkeypatch.setattr(pres, "_DIR", tmp_path / "presence")
        a = _agent()
        a.opts.cwd = tmp_path
        a._announce_presence()
        assert any(r.get("key") == "runde2-s1" for r in pres.open_sessions())

    def test_it_withdraws_when_it_leaves(self, tmp_path, monkeypatch):
        from delfin.agent import session_presence as pres
        monkeypatch.setattr(pres, "_DIR", tmp_path / "presence")
        a = _agent()
        a.opts.cwd = tmp_path
        a._announce_presence()
        a._withdraw_presence()
        assert not [r for r in pres.open_sessions() if r.get("key") == "runde2-s1"]

    def test_announcing_never_raises(self, monkeypatch):
        from delfin.agent import session_presence as pres
        monkeypatch.setattr(pres, "announce",
                            lambda *a, **k: (_ for _ in ()).throw(OSError("no")))
        a = _agent()
        a.opts.cwd = "/nonexistent"
        a._announce_presence()          # must not raise
