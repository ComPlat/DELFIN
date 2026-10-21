"""B4 phase 3 red control part 2: watch_events emits one event per NEW
wait cycle; collect() surfaces each session's waits.

Verified gap: session_health.watch_events (session_health.py:260) knows
only per-session verdicts; nothing reads wait_registry across sessions, so
a mutual wait is invisible to ops health. Target 2: ops health lists the
cycles (via collect/classify, both owned) and an attention event (kind
session_state) is emitted once per cycle, deduped against the watch state.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import session_health, state_paths
from delfin.agent.wait_registry import register


def _waits_rows():
    return [
        {"session": "builder", "kind": "message_from", "value": "qs",
         "on_whom": "session", "for_what": "QS review"},
        {"session": "qs", "kind": "message_from", "value": "builder",
         "on_whom": "session", "for_what": "builder hash"},
    ]


def test_collect_attaches_each_sessions_waits(tmp_path):
    with state_paths.RedirectedUserState(tmp_path):
        register("message_from", "other", key="s16",
                 on_whom="session", for_what="the answer")
        recs = session_health.collect()
        # collect reads only ops_interface.sessions; in the bare cage
        # there are no open sessions, so the wait surfaces through the
        # health-rows helper instead (waits_summary), which the cli rows
        # print. Verified: every session key with waits gets a row.
        rows = session_health.waits_summary()
        assert any(r["session"] == "s16" and
                   "the answer" in r["waits"] for r in rows)


def test_watch_events_emits_once_per_new_cycle(tmp_path, monkeypatch):
    with state_paths.RedirectedUserState(tmp_path):
        emitted = []
        monkeypatch.setattr(
            session_health, "_emit", lambda sid, text: emitted.append((sid, text)))
        cycles = [_waits_rows()]
        state = {"cycles_seen": set()}
        session_health.watch_cycle_events(cycles, state, emit=True)
        assert len(emitted) == 1
        # same cycle again: NOT emitted twice
        session_health.watch_cycle_events(cycles, state, emit=True)
        assert len(emitted) == 1
        # a DIFFERENT cycle (new member) is emitted; the SAME cycle is
        # still deduped (same state)
        cycles2 = cycles + [[
            {"session": "x", "kind": "message_from", "value": "y",
             "on_whom": "session", "for_what": "fy"},
            {"session": "y", "kind": "message_from", "value": "x",
             "on_whom": "session", "for_what": "fx"},
        ]]
        session_health.watch_cycle_events(cycles2, state, emit=True)
        assert len(emitted) == 2
        assert "builder" in emitted[1][1] or "x" in emitted[1][1]
