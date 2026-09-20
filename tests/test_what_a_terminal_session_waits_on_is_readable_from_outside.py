"""What a terminal session is waiting on can be read from outside it.

An approval prompt is rendered into a pane: the body is capped at 24
lines and every line is cut to the pane's width. On 2026-09-20 that put
an operator in front of

    $ grep -rn "mcp_isolation" ~/.delfin_…/proc /usr/bin/true' && echo BWRAP_OK

with the middle of a compound command missing, and in front of a change
to the self-modification guard whose diff said "… 30 more lines". The
module's own comment says why credentials are masked rather than cut:
"an approval you cannot read is one you cannot give." The width does the
same damage as the masking it avoids.

Widening the window afterwards does not help -- the block is already
drawn. And from outside, the only way to see that a session was waiting
at all was to capture its pane and grep for the frame character.

So the request is also published, whole, where a supervisor can read it:
the full preview, untruncated, with the session it belongs to.

  it is INFORMATION, not a second door   the answer is given at the
                                         terminal, and the record says so
  it disappears when the question does   a stale question read as open is
                                         worse than none at all
  publishing may never break the prompt  a supervisor's convenience does
                                         not get to fail an approval
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import terminal_confirm as tc


@pytest.fixture
def room(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "terminal_confirmations")
    return tmp_path


LONG = "\n".join(f"line {i:03d} " + "x" * 200 for i in range(60))


def _broker(**kw):
    return tc.TerminalConfirmBroker(session_id="s-123", session_key="runde2-s1",
                                    **kw)


def _ask(broker, preview=LONG, tool="bash"):
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool=tool,
                            args={"command": "echo hi"}, preview=preview)
    return broker._enqueue(req)


class TestPublishing:
    def test_a_waiting_request_shows_up(self, room):
        _ask(_broker())
        assert len(tc.pending_at_terminals()) == 1

    def test_the_preview_is_not_cut(self, room):
        _ask(_broker())
        rec = tc.pending_at_terminals()[0]
        assert rec["preview"] == LONG
        assert len(rec["preview"].splitlines()) == 60

    def test_it_says_which_session(self, room):
        _ask(_broker())
        rec = tc.pending_at_terminals()[0]
        assert rec["session_key"] == "runde2-s1"
        assert rec["session_id"] == "s-123"

    def test_it_says_where_it_can_be_answered(self, room):
        # The first version said "terminal", which was right while there
        # was no checked way in from outside. There is one now --
        # file_confirm's, with all of its checks -- so the record says
        # both ends, and the terminal still wins a race.
        _ask(_broker())
        assert tc.pending_at_terminals()[0]["answer_at"] == "terminal or supervisor"

    def test_an_answer_takes_the_question_away(self, room):
        b = _broker()
        req = _ask(b)
        b.resolve(req, True)
        assert tc.pending_at_terminals() == []

    def test_an_abort_takes_every_question_away(self, room):
        b = _broker()
        _ask(b)
        _ask(b)
        assert len(tc.pending_at_terminals()) == 2
        b.abort_all()
        assert tc.pending_at_terminals() == []

    def test_two_sessions_do_not_overwrite_each_other(self, room):
        _ask(_broker())
        _ask(tc.TerminalConfirmBroker(session_id="s-999",
                                      session_key="runde2-s2"))
        keys = sorted(r["session_key"] for r in tc.pending_at_terminals())
        assert keys == ["runde2-s1", "runde2-s2"]

    def test_the_record_is_owner_only(self, room):
        _ask(_broker())
        import os
        f = next((room / "terminal_confirmations").rglob("*.request.json"))
        assert oct(os.stat(f).st_mode & 0o777) == "0o600"


class TestItNeverBreaksTheApproval:
    def test_a_failing_publish_is_silent(self, room, monkeypatch):
        monkeypatch.setattr(tc, "_publish_pending",
                            lambda *a, **k: (_ for _ in ()).throw(OSError("full")))
        b = _broker()
        req = _ask(b)                      # must not raise
        assert req.resolved is False

    def test_a_failing_withdraw_is_silent(self, room, monkeypatch):
        b = _broker()
        req = _ask(b)
        monkeypatch.setattr(tc, "_withdraw_pending",
                            lambda *a, **k: (_ for _ in ()).throw(OSError("gone")))
        assert b.resolve(req, True) is True

    def test_a_broker_without_a_session_still_works(self, room):
        b = tc.TerminalConfirmBroker()
        req = _ask(b)
        assert b.resolve(req, True) is True

    def test_an_unreadable_room_reads_as_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "nope" / "deeper")
        assert tc.pending_at_terminals() == []
