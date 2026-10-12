"""An approval request carries the id of the session that raised it.

Input: two sessions in one dashboard, one ``KitConfirmBroker`` (the tab
builds exactly one), an approval window that expires in session C.
Output: the inbox entry names C; session A's turn-opening reminder
(``AgentEngine._remind_of_waiting_approvals``) says nothing about it;
session C's does.

Defect pinned: the dashboard never set ``KitConfirmBroker.session_id``,
so every dashboard entry was unattributed, and the reminder counted an
unattributed entry younger than one hour as its own. A session whose
own edit_file calls had all passed was told an edit_file request of
another session was still waiting for it, and acted on it.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest


class _Client:
    def __init__(self):
        self.notes: list[str] = []

    def push_run_note(self, text):
        self.notes.append(text)


def _engine(sid: str):
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng.client = _Client()
    eng.session_id = sid
    return eng


@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", lambda: h)
    from delfin.agent import terminal_confirm as tc
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "room")
    return h


def _dashboard_broker(state: dict):
    """The broker as ``_ensure_kit_broker`` binds it: following whatever
    engine the tab currently holds."""
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=0.05)
    broker.follow_session(
        lambda: getattr(state.get("engine"), "session_id", ""))
    return broker


def _expire(broker, path="a.py"):
    broker._away_since = None          # each request gets its own window
    assert broker.callback("edit_file", {"path": path}, "") is False
    assert broker.last_timed_out is True


def _pending():
    from delfin.agent import attention
    return attention.list_pending("confirm_pending")


def test_another_sessions_expired_request_is_not_mentioned(home):
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=0.05)   # as built before
    _expire(broker)
    assert len(_pending()) == 1
    a = _engine("session-a")
    a._remind_of_waiting_approvals()
    assert a.client.notes == [], a.client.notes


def test_the_session_that_asked_is_reminded_and_no_other(home):
    a, c = _engine("session-a"), _engine("session-c")
    state = {"engine": c}
    broker = _dashboard_broker(state)
    _expire(broker, "c.py")
    assert [e["session_id"] for e in _pending()] == ["session-c"]

    state["engine"] = a                    # the user switches session
    a._remind_of_waiting_approvals()
    assert a.client.notes == []
    c._remind_of_waiting_approvals()
    assert len(c.client.notes) == 1 and "edit_file" in c.client.notes[0]


def test_a_switch_and_a_new_id_are_followed(home):
    a, c = _engine("session-a"), _engine("session-c")
    state = {"engine": c}
    broker = _dashboard_broker(state)
    _expire(broker, "1.py")
    state["engine"] = a
    _expire(broker, "2.py")
    a.session_id = "session-a2"           # the engine minted a new id
    _expire(broker, "3.py")
    state["engine"] = None                 # no engine: honest about it
    _expire(broker, "4.py")
    owners = sorted(str(e.get("session_id") or "") for e in _pending())
    assert owners == ["", "session-a", "session-a2", "session-c"]


def test_the_room_record_names_the_session_too(home):
    state = {"engine": _engine("session-c")}
    broker = _dashboard_broker(state)
    seen = {}
    from delfin.agent import terminal_confirm as tc
    real = tc.publish_waiting

    def _spy(**kw):
        seen.update(kw)
        return real(**kw)
    import unittest.mock as um
    with um.patch.object(tc, "publish_waiting", _spy):
        _expire(broker)
    assert seen.get("session_id") == "session-c"


def test_a_set_id_still_works_without_a_source(home):
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=0.05)
    broker.session_id = "fixed"
    _expire(broker)
    assert [e["session_id"] for e in _pending()] == ["fixed"]


def test_a_failing_source_does_not_break_the_dialog(home):
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=0.05)

    def _boom():
        raise RuntimeError("engine gone")
    broker.follow_session(_boom)
    _expire(broker)
    assert [e["session_id"] for e in _pending()] == [""]


def test_the_dashboard_binds_the_broker_to_the_current_engine():
    from delfin.dashboard import tab_agent as T
    src = inspect.getsource(T)
    i = src.index("def _ensure_kit_broker()")
    body = src[i:src.index("return state.get(\"_kit_confirm_broker\")", i)]
    assert "broker.follow_session(" in body
    assert 'state.get("engine")' in body[body.index("broker.follow_session("):]
