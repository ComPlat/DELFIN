"""The terminal runs at most three "announced but not done" follow-ups in a row.

Wave-12 sessions ended turns with "Let me …" and sat idle; the engine now
records one follow-up note per such turn and the terminal runs it. A model
that ONLY ever announces must not loop: after _ANNOUNCE_FOLLOWUP_CAP
follow-ups in a row the session waits for real input, and a paused session
gets none at all.
"""
from __future__ import annotations

from delfin.agent import repl as R
from delfin.agent import session_pause
from delfin.agent.engine import AgentEngine


class _AlwaysAnnouncing:
    """An engine whose every turn announced work again."""

    def __init__(self):
        self.cleared = 0
        self.pending_turn_continuation = "continue the announced step"

    def clear_turn_continuation(self):
        self.cleared += 1
        # the next turn announces again
        self.pending_turn_continuation = "continue the announced step"


def _agent(engine, monkeypatch):
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    a.engine = engine
    monkeypatch.setattr(a, "_presence_key", lambda: "sess-test", raising=False)
    monkeypatch.setattr(session_pause, "wake_blocked", lambda key: False)
    return a


def test_an_always_announcing_model_gets_at_most_the_cap(monkeypatch):
    eng = _AlwaysAnnouncing()
    a = _agent(eng, monkeypatch)
    a._announce_followups = 0
    got = [a._continue_after_announcement() for _ in range(6)]
    cap = R.TerminalAgent._ANNOUNCE_FOLLOWUP_CAP
    assert sum(1 for g in got if g) == cap, got
    assert all(not g for g in got[cap:])
    assert eng.cleared == 6, "the note must be consumed every time"


def test_real_input_rearms_the_cap(monkeypatch):
    a = _agent(_AlwaysAnnouncing(), monkeypatch)
    a._announce_followups = R.TerminalAgent._ANNOUNCE_FOLLOWUP_CAP
    assert a._continue_after_announcement() == ""
    a._announce_followups = 0          # what run() does on a real-input turn
    assert a._continue_after_announcement()


def test_a_paused_session_gets_no_follow_up(monkeypatch):
    a = _agent(_AlwaysAnnouncing(), monkeypatch)
    monkeypatch.setattr(session_pause, "wake_blocked", lambda key: True)
    a._announce_followups = 0
    assert a._continue_after_announcement() == ""


def test_no_note_means_no_follow_up(monkeypatch):
    eng = object.__new__(AgentEngine)
    a = _agent(eng, monkeypatch)
    a._announce_followups = 0
    assert a._continue_after_announcement() == ""


def test_the_engine_latch_records_one_note_for_an_announcing_answer(monkeypatch):
    eng = object.__new__(AgentEngine)
    monkeypatch.setattr(AgentEngine, "_open_task_count", lambda self: 2)
    eng._note_turn_continuation("Let me write the next test.")
    assert eng.pending_turn_continuation
    first = eng.pending_turn_continuation
    eng._note_turn_continuation("I will now run the suite.")
    assert eng.pending_turn_continuation == first, "latch must hold"
    eng.clear_turn_continuation()
    assert eng.pending_turn_continuation == ""


def test_no_open_tasks_means_no_note(monkeypatch):
    eng = object.__new__(AgentEngine)
    monkeypatch.setattr(AgentEngine, "_open_task_count", lambda self: 0)
    eng._note_turn_continuation("Let me write the next test.")
    assert eng.pending_turn_continuation == ""
