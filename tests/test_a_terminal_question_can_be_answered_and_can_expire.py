"""A question a terminal session asks can be answered from outside, and can expire.

Found running wave 14 (2026-10-09): every choice question a terminal
session asked (ask_user_question) was unanswerable from outside --
`delfin-agent approvals answer` failed with "no request file", because the
answer was written to the headless room while terminal sessions publish
their questions in their own. And a terminal session had no timeout at
all: a session nobody watched waited for hours on one question. Now the
answer goes beside the question, and `chat --confirm-timeout` lets an
unanswered question expire as "not now" -- never as a refusal.
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import approval_answers as aa
from delfin.agent import file_confirm as FC
from delfin.agent import security_events
from delfin.agent import terminal_confirm as tc


@pytest.fixture
def rooms(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "terminal")
    monkeypatch.setattr(FC, "requests_dir", lambda: tmp_path / "headless")
    security_events.clear()
    return tmp_path


def _ask(timeout_s=6.0):
    broker = tc.TerminalConfirmBroker(session_id="s-1", session_key="wave-s14",
                                      timeout_s=timeout_s, poll_s=0.02)
    req = tc.ConfirmRequest(
        kind=tc.ASK, tool="ask_user_question", args={}, preview="When?",
        payload={"question": "When should it check?",
                 "options": [{"label": "Before asking"},
                             {"label": "After approval"},
                             {"label": "Both"}]})
    out: dict = {}

    def _run():
        broker._enqueue(req)
        out["decision"] = broker._wait(req)
        # Per thread, like the gate reads it: off __self__ in the asking thread.
        out["timed_out"] = broker.last_timed_out

    t = threading.Thread(target=_run, daemon=True)
    t.start()
    for _ in range(300):
        if tc.pending_at_terminals():
            break
        time.sleep(0.01)
    return broker, req, t, out


def test_an_answer_from_outside_reaches_a_terminal_question(rooms):
    broker, req, t, out = _ask()
    rows = tc.pending_at_terminals()
    assert rows, "the question was not published"
    picks = aa.answer(rows[0]["id"], "3", by="operator")
    t.join(timeout=5)
    assert picks == ["Both"]
    assert not t.is_alive(), "the session is still waiting"
    assert out["decision"] == {"answers": ["Both"]}
    # It went beside the question, not into the headless room.
    assert not (rooms / "headless").exists() or not any(
        (rooms / "headless").glob("*.answer.json"))


def test_an_unanswered_question_expires_as_not_now(rooms):
    broker, req, t, out = _ask(timeout_s=0.3)
    t.join(timeout=5)
    assert not t.is_alive()
    assert req.expired is True
    assert out["timed_out"] is True


def test_the_chat_command_takes_a_confirm_timeout():
    from delfin.agent import cli
    parser = cli.build_parser()
    args = parser.parse_args(["chat", "--confirm-timeout", "900", "hi"])
    assert args.confirm_timeout == 900.0
    args = parser.parse_args(["chat", "hi"])
    assert args.confirm_timeout == 0.0
