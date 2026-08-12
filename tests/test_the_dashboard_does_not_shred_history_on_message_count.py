"""A second compactor ran after every dashboard turn and kept six messages.

``_check_auto_compact`` fired unconditionally at the end of every turn.
Over 30 messages it kept the last 6 and, in place of everything else, a
200-character excerpt of only the newest 8 of the ones it dropped. On a
31-message session that is everything before index 17, gone. It computed
a token total and never used it, so it ignored the budget entirely and
fired at roughly 15% of the window -- re-creating, in the dashboard, the
exact failure the engine's own compaction comment records as fixed
("compacted at 15% full -> agent confused").

It also ignored pinned messages, wrote no archive, and left no elided
record, so nothing it cut had any retrieval path at all. The ``/compact``
command had the same shape.

The engine already compacts: inside the turn, under token pressure, with
pins exempt, the pre-compaction transcript archived and a
``history_get('elided:<ref>')`` marker left behind for every body it
shortens. There is no argument for a second one, so both are gone and
``/compact`` routes to the engine's.
"""

from __future__ import annotations

import inspect
import pathlib

from delfin.agent.engine import AgentEngine

_TAB = pathlib.Path(
    inspect.getfile(__import__("delfin.dashboard.tab_agent",
                               fromlist=["x"]))).read_text(encoding="utf-8")


def _bare_engine() -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.context_window_tokens = 100_000
    eng.auto_compact_pct = 0.95
    eng.last_compaction_info = None
    eng.session_id = ""
    eng.backend = "api"
    eng.client = None
    eng.current_role_index = 0
    eng.route = ["solo_agent"]
    return eng


# ---------------------------------------------------------------------------
# Both compactors are gone
# ---------------------------------------------------------------------------

def test_the_post_turn_compactor_is_gone():
    assert "_check_auto_compact" not in _TAB, (
        "the message-count compactor is back: it shreds a 31-message "
        "session at 15% of the window, with no archive and nothing to "
        "retrieve any of it from")


def test_no_dashboard_code_rebuilds_the_history_from_an_excerpt():
    """The signature of both: replace engine.messages with a
    hand-assembled summary + a fixed tail."""
    assert "[Auto-compacted" not in _TAB
    assert "[Context summary of" not in _TAB


def test_the_compact_command_routes_to_the_engine():
    assert "engine._compact_history(force=True)" in _TAB


# ---------------------------------------------------------------------------
# The engine's compactor, which the command now uses
# ---------------------------------------------------------------------------

def test_a_short_session_at_low_pressure_is_left_alone():
    """31 messages at 15% of the window: the old compactor's trigger."""
    eng = _bare_engine()
    eng.messages = [{"role": "user" if i % 2 == 0 else "assistant",
                     "content": f"message {i}"} for i in range(31)]
    before = list(eng.messages)
    eng._compact_history()
    assert eng.messages == before


def test_a_forced_compaction_still_protects_pinned_messages():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "PINNED GOAL", "_pinned": True}]
    eng.messages += [{"role": "assistant", "content": f"answer {i}"}
                     for i in range(20)]
    eng._compact_history(force=True)
    assert any(m.get("_pinned") for m in eng.messages), (
        "the pinned message was compacted away")
    assert any(m.get("content") == "PINNED GOAL" for m in eng.messages)


def test_a_forced_compaction_archives_what_it_removes(monkeypatch):
    archived = {}

    def _archive(session_id, msgs, info=None):
        archived["msgs"] = list(msgs)

    monkeypatch.setattr(
        "delfin.agent.session_store.archive_pre_compaction_transcript",
        _archive)
    eng = _bare_engine()
    eng.session_id = "s1"
    eng.messages = [{"role": "assistant", "content": f"answer {i}"}
                    for i in range(20)]
    eng._compact_history(force=True)
    assert len(archived.get("msgs") or []) >= 15, (
        "compaction removed messages without archiving them")


def test_a_forced_compaction_on_a_two_message_session_does_nothing():
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "hi"},
                    {"role": "assistant", "content": "hello"}]
    eng._compact_history(force=True)
    assert len(eng.messages) == 2
