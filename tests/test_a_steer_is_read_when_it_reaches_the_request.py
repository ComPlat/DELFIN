"""A session message is read when it reaches the model's request, and a
blocked task is stale only against the last USER message.

Two defects of the first wiring, each pinned here against the real objects:

* The steer path marked the message read when ``engine.steer`` QUEUED it.
  A turn that stops before its next round hands the queued text back as
  "Not delivered" -- the receipt said ``read`` for text the model never
  received. The read mark now rides on the steer (``on_delivered``) and
  runs only in a drain that puts the text into a request.
* The stale reference was the session file's ``updated_at``. Every turn
  end saves the session, so a block set during a turn was older than the
  save the moment the turn ended: every fresh block read "(stale)". The
  reference is now ``engine.last_user_message_at``, set when a user turn
  starts.
"""

from __future__ import annotations

import threading
import time
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import agent_tasks, session_messages as msgs
from delfin.agent import task_ticker as TT
from delfin.agent.api_client import OpenAIClient


def _client():
    c = OpenAIClient.__new__(OpenAIClient)
    c._steer_lock = threading.Lock()
    c._steer_msgs = []
    return c


@pytest.fixture
def inbox(tmp_path, monkeypatch):
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "session_inbox")
    return tmp_path


class _Theme:
    def dim(self, t):
        return t


class _Transcript:
    theme = _Theme()

    def __init__(self):
        self.lines = []

    def chrome(self, line):
        self.lines.append(line)


class _Opts:
    session_name = "target-s1"

    def __init__(self, cwd=None):
        self.cwd = cwd


def _terminal(engine, cwd=None):
    from delfin.agent.repl import TerminalAgent
    a = object.__new__(TerminalAgent)
    a.opts = _Opts(cwd)
    a.engine = engine
    a.transcript = _Transcript()
    a._mail_steer_last = 0.0
    return a


def _engine_with_real_queue(tmp_path):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                          model="kit.qwen3.5-397b-A17b", mode="solo")
    eng.client = _client()
    return eng


# -- the client queue ------------------------------------------------------

def test_the_callback_runs_on_a_delivering_drain():
    c = _client()
    seen = []
    c.push_steer("look at this", on_delivered=lambda: seen.append(1))
    assert seen == []
    assert c._drain_steer() == ["look at this"]
    assert seen == [1]


def test_the_callback_is_dropped_by_a_stop_drain():
    c = _client()
    seen = []
    c.push_steer("look at this", on_delivered=lambda: seen.append(1))
    assert c._drain_steer(delivered=False) == ["look at this"]
    assert c._drain_steer() == []
    assert seen == []


def test_a_failing_callback_does_not_lose_the_text():
    c = _client()

    def boom():
        raise OSError("receipts unwritable")
    c.push_steer("still delivered", on_delivered=boom)
    assert c._drain_steer() == ["still delivered"]


# -- the terminal's running-turn mail --------------------------------------

def test_a_steered_message_is_not_read_while_it_is_only_queued(inbox):
    sent = msgs.send("target-s1", "steer into this turn", from_key="ops")
    eng = _engine_with_real_queue(inbox)
    _terminal(eng)._steer_operator_mail()
    assert eng.client._steer_msgs, "the mail must be queued for the turn"
    assert msgs.status(sent["id"], "ops") == "delivered"


def test_a_steered_message_is_read_once_the_request_takes_it(inbox):
    sent = msgs.send("target-s1", "steer into this turn", from_key="ops")
    eng = _engine_with_real_queue(inbox)
    _terminal(eng)._steer_operator_mail()
    eng.client._drain_steer()
    assert msgs.status(sent["id"], "ops") == "read"


def test_a_steered_message_handed_back_by_a_stop_is_not_read(inbox):
    sent = msgs.send("target-s1", "steer into this turn", from_key="ops")
    eng = _engine_with_real_queue(inbox)
    _terminal(eng)._steer_operator_mail()
    note = eng._hand_back_undelivered_steer()
    assert "Not delivered" in note
    assert msgs.status(sent["id"], "ops") != "read"


# -- the dashboard's delivery ----------------------------------------------

def test_the_dashboard_steer_carries_the_read_mark_instead_of_setting_it():
    """Driving _deliver_session_messages needs the whole dashboard; the
    shape is asserted on its source: the steer branch hands the read mark
    to push_steer and does not mark read itself."""
    from pathlib import Path
    import delfin.dashboard.tab_agent as tab
    src = Path(tab.__file__).read_text(encoding="utf-8")
    i = src.index("def _deliver_session_messages")
    body = src[i:src.index("\n    def ", i + 10)]
    assert "push_steer(text, on_delivered=_read)" in body
    steer_line = body.index("push_steer(text, on_delivered=_read)")
    after = body[steer_line:]
    assert after.index("elif _send_on_its_own(") < after.index("_read()")


# -- the stale reference ---------------------------------------------------

@pytest.fixture
def store(tmp_path):
    agent_tasks._STORES.clear()
    return agent_tasks.get_store(tmp_path)


def test_a_fresh_block_is_not_stale_after_the_turn_end_save(
        store, tmp_path, monkeypatch):
    """The turn that blocked the task ends and the session is saved; nobody
    has said anything since. The terminal's task report must not call the
    block stale."""
    from delfin.agent import session_store
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sessions")
    (tmp_path / "sessions").mkdir()
    eng = MagicMock()
    eng.session_id = "s1"
    eng.last_user_message_at = time.time() - 30.0   # the message that began it
    t = store.create("run the opt", session_id="s1")
    store.update(t["id"], status="blocked",
                 blocked_reason="the user's choice of basis")
    time.sleep(1.1)                                  # the turn goes on
    session_store.save_session("s1", chat_messages=[
        {"role": "user", "content": "optimise it"}], workspace=str(tmp_path))
    a = _terminal(eng, cwd=tmp_path)
    a._show_tasks = True
    a._report_tasks()
    shown = "\n".join(a.transcript.lines)
    assert "waiting on the user's choice of basis" in shown
    assert "(stale)" not in shown


def test_a_block_before_the_last_user_message_is_stale(store, tmp_path):
    eng = MagicMock()
    eng.session_id = "s1"
    t = store.create("run the opt", session_id="s1")
    store.update(t["id"], status="blocked", blocked_reason="the basis")
    eng.last_user_message_at = time.time() + 5.0
    a = _terminal(eng, cwd=tmp_path)
    a._show_tasks = True
    a._report_tasks()
    assert "(stale)" in "\n".join(a.transcript.lines)


def test_a_block_in_the_same_second_as_the_message_is_not_stale(store, tmp_path):
    """The store keeps whole seconds: a block written at 12:00:00.9 reads as
    12:00:00 and must not be called older than a message at 12:00:00.5."""
    t = store.create("x", session_id="s")
    store.update(t["id"], status="blocked", blocked_reason="r")
    floor = TT._block_epoch(store.get(t["id"]))
    row = TT.rows(tmp_path, session_id="s", last_message_at=floor + 0.5)
    assert "(stale)" not in row[0]["label"]
    row = TT.rows(tmp_path, session_id="s", last_message_at=floor + 1.0)
    assert "(stale)" in row[0]["label"]


def test_a_user_turn_records_when_it_began(tmp_path):
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                          model="kit.qwen3.5-397b-A17b", mode="solo")

    def answer(*a, **kw):
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", output_tokens=1, cost_usd=0.0)
    eng.client.stream_message = MagicMock(side_effect=answer)
    before = time.time()
    eng.stream_response("hallo")
    after = time.time()
    assert before <= float(eng.last_user_message_at) <= after
