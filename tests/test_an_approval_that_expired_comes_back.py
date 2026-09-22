"""An expired approval is not a refusal — and nothing said when to ask again.

The window is 300 s. When it closes unanswered the model is told the
right thing — absence, not a denial, carry on and do not retry now —
and the request is parked in the USER's inbox (`/attention`). Nobody
tells the AGENT when the user comes back, so a task whose only remaining
step needs that approval simply stops.

Measured on 2026-09-18: one session lost the last quarter of its work to
exactly this, twice. The user was away with four sessions running; the
second request expired without even being shown, because the broker
latches "away" for the length of the window.

A turn now opens by saying what is still waiting — once per request, so
a user who is still away is not nagged every turn.
"""

from __future__ import annotations

import pytest


class _Client:
    def __init__(self):
        self.notes: list[str] = []

    def push_run_note(self, text):
        self.notes.append(text)


SESSION = "this-session"


@pytest.fixture()
def engine(monkeypatch):
    from delfin.agent.engine import AgentEngine as DelfinAgent
    eng = DelfinAgent.__new__(DelfinAgent)
    eng.client = _Client()
    eng.session_id = SESSION
    return eng


def _mine(eid, title):
    return {"id": eid, "title": title, "session_id": SESSION}


def _pending(monkeypatch, items):
    from delfin.agent import attention as att
    monkeypatch.setattr(att, "list_pending", lambda kind=None: list(items))


def test_a_waiting_approval_opens_the_turn(engine, monkeypatch):
    _pending(monkeypatch, [_mine("e1", "edit delfin/dashboard/tab_agent.py")])
    engine._remind_of_waiting_approvals()
    assert engine.client.notes, "the turn was never told"
    note = engine.client.notes[0]
    assert "still waiting" in note
    assert "tab_agent.py" in note, "it names what is waiting"
    assert "absence, not a refusal" in note
    assert "ask again now" in note


def test_it_says_nothing_when_nothing_waits(engine, monkeypatch):
    _pending(monkeypatch, [])
    engine._remind_of_waiting_approvals()
    assert engine.client.notes == []


def test_the_same_request_is_mentioned_once(engine, monkeypatch):
    """A user who is still away must not be nagged every turn."""
    _pending(monkeypatch, [_mine("e1", "edit x.py")])
    engine._remind_of_waiting_approvals()
    engine._remind_of_waiting_approvals()
    assert len(engine.client.notes) == 1


def test_a_new_request_is_mentioned(engine, monkeypatch):
    _pending(monkeypatch, [_mine("e1", "edit x.py")])
    engine._remind_of_waiting_approvals()
    _pending(monkeypatch, [_mine("e1", "edit x.py"), _mine("e2", "edit y.py")])
    engine._remind_of_waiting_approvals()
    assert len(engine.client.notes) == 2
    assert "y.py" in engine.client.notes[1]
    assert "x.py" not in engine.client.notes[1], "only what is new"


def test_it_never_raises(engine, monkeypatch):
    from delfin.agent import attention as att

    def _boom(kind=None):
        raise RuntimeError("no inbox")
    monkeypatch.setattr(att, "list_pending", _boom)
    engine._remind_of_waiting_approvals()      # must not raise
    assert engine.client.notes == []


def test_a_turn_asks_before_it_starts():
    """Wired into stream_response, not left as a helper nobody calls."""
    import ast
    import inspect
    from delfin.agent import engine as E

    tree = ast.parse(inspect.getsource(E))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "stream_response")
    assert "_remind_of_waiting_approvals" in ast.unparse(fn)


# -- whose request is it, anyway --------------------------------------------

def test_another_sessions_request_is_not_mine(engine, monkeypatch):
    """Reported from the field: a fresh session that had said nothing but
    "Hallo" was told seven requests were waiting for it, all of them from
    sessions that had ended hours before. It spent a turn working that
    out."""
    _pending(monkeypatch, [{"id": "e1", "title": "edit x.py",
                            "session_id": "some-other-session"}])
    engine._remind_of_waiting_approvals()
    assert engine.client.notes == []


def test_an_old_unattributed_request_is_not_mine(engine, monkeypatch):
    """Entries from before requests carried a session say nothing about
    whose they are, so they count only while they are recent."""
    import time
    _pending(monkeypatch, [{"id": "e1", "title": "edit x.py",
                            "created_at": time.time() - 4 * 3600}])
    engine._remind_of_waiting_approvals()
    assert engine.client.notes == []


def test_a_recent_unattributed_request_still_counts(engine, monkeypatch):
    import time
    _pending(monkeypatch, [{"id": "e1", "title": "edit x.py",
                            "created_at": time.time() - 60}])
    engine._remind_of_waiting_approvals()
    assert engine.client.notes, "a request from minutes ago is plausibly mine"
