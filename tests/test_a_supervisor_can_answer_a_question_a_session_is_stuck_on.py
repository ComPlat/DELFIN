"""A supervisor can answer a question a terminal session is stuck on.

Reading it was the first half and is done: the question is published
whole, because the pane cuts it to 24 lines and to its own width.
Answering was still a keystroke into that pane -- and an approval dialog
reads single keys, so a supervisor who typed while one was up answered
it blind.

The decision this reverses is the one the first version of the published
record stated in the record itself: "answer_at: terminal". That was
right while there was no checked way in. There is one now, and it is not
a new one: file_confirm already answers headless sessions, and every
careful part of it applies unchanged -- the answer may not be a symlink,
must belong to this user, must not be group- or other-writable, must
name this very question, and must not predate it.

  the terminal still wins a race   resolve() is the single point that
                                   decides, and the first answer takes
                                   it; a later one is a no-op
  it is recorded, loudly           an approval given from outside the
                                   pane is a security event, so the
                                   panel shows it like any other
  reading beats guessing           the supervisor decides on the WHOLE
                                   preview, which is more than the pane
                                   was ever able to show
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import file_confirm as FC
from delfin.agent import security_events
from delfin.agent import terminal_confirm as tc


@pytest.fixture
def room(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    security_events.clear()
    return tmp_path / "pending"


def _stuck(preview="[SELF-MODIFICATION GUARD]\nthe whole diff"):
    broker = tc.TerminalConfirmBroker(session_id="s-1", session_key="runde2-s1",
                                      timeout_s=6.0, poll_s=0.02)
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="edit_file",
                            args={"command": "c"}, preview=preview)
    out: dict = {}

    def _run():
        broker._enqueue(req)
        out["decision"] = broker._wait(req)

    t = threading.Thread(target=_run, daemon=True)
    t.start()
    for _ in range(200):
        if tc.pending_at_terminals():
            break
        time.sleep(0.01)
    return broker, req, t, out


def _the_question():
    rows = tc.pending_at_terminals()
    assert rows, "nothing was published"
    return rows[0]["id"]


class TestAnsweringFromOutside:
    def test_an_approval_reaches_the_waiting_session(self, room):
        broker, req, t, out = _stuck()
        assert tc.answer_waiting(_the_question(), True, by="operator")
        t.join(6)
        assert out["decision"] is True

    def test_a_refusal_reaches_it_too(self, room):
        broker, req, t, out = _stuck()
        assert tc.answer_waiting(_the_question(), False, by="operator")
        t.join(6)
        assert out["decision"] is False

    def test_the_question_leaves_the_listing(self, room):
        broker, req, t, out = _stuck()
        tc.answer_waiting(_the_question(), True)
        t.join(6)
        assert tc.pending_at_terminals() == []

    def test_it_is_recorded_as_a_security_event(self, room):
        broker, req, t, out = _stuck()
        tc.answer_waiting(_the_question(), True, by="operator")
        t.join(6)
        kinds = [e.kind for e in security_events.recent(20)]
        assert "approval_from_outside" in kinds

    def test_the_kind_has_a_label(self):
        assert "approval_from_outside" in security_events.known_kinds()

    def test_an_answer_to_a_question_nobody_asked_is_refused(self, room):
        assert tc.answer_waiting("1789-deadbeef", True) is False


class TestTheTerminalStillWins:
    def test_the_first_answer_takes_it(self, room):
        broker, req, t, out = _stuck()
        assert broker.resolve(req, "terminal-said-this") is True
        t.join(6)
        # The outside answer arrives after the question is gone.
        assert tc.answer_waiting(_the_question() if tc.pending_at_terminals()
                                 else "1789-gone", True) is False
        assert out["decision"] == "terminal-said-this"


class TestTheCarefulPartsStillApply:
    def test_a_symlinked_answer_is_ignored(self, room, tmp_path):
        broker, req, t, out = _stuck()
        rid = _the_question()
        target = tmp_path / "elsewhere.json"
        target.write_text('{"id": "%s", "decision": "approve"}' % rid,
                          encoding="utf-8")
        (room / f"{rid}.answer.json").symlink_to(target)
        t.join(8)
        assert out["decision"] is False, "an expired question refuses"

    def test_an_answer_naming_another_question_is_ignored(self, room):
        broker, req, t, out = _stuck()
        rid = _the_question()
        FC.answer(rid, True, room=room, by="x")
        (room / f"{rid}.answer.json").write_text(
            '{"id": "somebody-else", "decision": "approve"}', encoding="utf-8")
        t.join(8)
        assert out["decision"] is False
