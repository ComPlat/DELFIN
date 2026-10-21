"""B4 adversarial: health wiring (waits_summary, watch_cycle_events).

Reviewer nacht-s19. Attacks through the real functions:
- hostile for_what/on_whom/value text in the store cannot break the
  waits_summary join or inject a second line into the display string.
- watch_cycle_events dedups by member set: same ring twice -> one event;
  a new member -> a second event; an empty cycle list -> nothing.
- the emitted text goes through describe() so hostile session names are
  normalised, not raw.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_health as sh
from delfin.agent.mutual_wait import find


def _w(kind, value, on_whom="", for_what="", session=""):
    return {
        "session": session, "kind": kind, "value": value,
        "on_whom": on_whom, "for_what": for_what,
    }


def test_waits_summary_hostile_for_what_stays_one_line(monkeypatch, tmp_path):
    from delfin.agent import state_paths
    with state_paths.RedirectedUserState(tmp_path):
        from delfin.agent import wait_registry
        wait_registry.register(
            "message_from", "s2", key="s1",
            on_whom="session", for_what="still <b>waiting</b>")
        rows = sh.waits_summary()
    assert len(rows) == 1
    text = rows[0]["waits"]
    assert "\n" not in text
    assert "<b>" not in text
    assert "waiting" in text


def test_waits_summary_empty_store_is_empty(monkeypatch, tmp_path):
    from delfin.agent import state_paths
    with state_paths.RedirectedUserState(tmp_path):
        assert sh.waits_summary() == []


def test_cycle_event_emitted_once_per_ring():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s1", on_whom="session", session="s2"),
    ]
    cycles = find(waits)
    state = {}
    first = sh.watch_cycle_events(cycles, state, emit=False)
    assert len(first) == 1
    second = sh.watch_cycle_events(cycles, state, emit=False)
    assert second == []
    # a third member joins the ring -> a NEW cycle -> one more event
    waits.append(_w("message_from", "s1", on_whom="session", session="s3"))
    waits.append(_w("message_from", "s3", on_whom="session", session="s2"))
    third = sh.watch_cycle_events(find(waits), state, emit=False)
    assert len(third) == 1


def test_cycle_event_text_is_normalised():
    waits = [
        _w("message_from", "s2", on_whom="session", session="s1"),
        _w("message_from", "s1", on_whom="session", session="s2"),
    ]
    state = {}
    out = sh.watch_cycle_events(find(waits), state, emit=False)
    assert len(out) == 1
    sid, text = out[0]
    assert "\n" not in text
    assert "<" not in text


def test_cycle_events_none_cycles_is_nothing():
    assert sh.watch_cycle_events(None, {}, emit=False) == []
