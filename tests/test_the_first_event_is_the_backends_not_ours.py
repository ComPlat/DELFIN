"""``first_event_ms`` is stamped by the backend's first event, not by ours.

The client emits a "waiting" notice every _WAIT_TICK_MAX_S while nothing
has arrived. It is an event on the same stream, so the stamp that asks
"did the transport deliver anything at all?" was answered by our own
clock: first_event_ms sat at ~10040 on 19 of 22 turns of one session and
11 of 18 of another (2026-10-08), while the real first token took 9-72 s
(ttft_ms). The instrument reported its cadence.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import turn_metrics as tm
from delfin.agent.api_client import StreamEvent
from tests.test_a_crashed_turn_was_logged_as_a_silent_backend import (  # noqa: F401
    _engine, _last, agent_tree, home)


def test_a_waiting_tick_does_not_stamp_the_first_event(agent_tree, home):
    def slow(system, messages, **kw):
        yield StreamEvent(type="waiting", text="waiting 10 s")
        time.sleep(0.25)
        yield StreamEvent(type="message_start", input_tokens=900)
        yield StreamEvent(type="text_delta", text="hi")

    eng = _engine(agent_tree, slow)
    eng.stream_response("hi")
    entry = _last(eng)
    assert entry["first_event_ms"] is not None
    assert entry["first_event_ms"] >= 200, (
        "the first event was stamped at our own waiting tick, before the "
        f"backend sent anything: {entry['first_event_ms']} ms")


def test_a_backend_event_still_stamps_it(agent_tree, home):
    def quick(system, messages, **kw):
        yield StreamEvent(type="message_start", input_tokens=900)

    eng = _engine(agent_tree, quick)
    eng.stream_response("hi")
    assert _last(eng)["first_event_ms"] is not None


def test_only_waiting_ticks_leave_it_unstamped(agent_tree, home):
    """A turn that saw nothing but our own ticks heard nothing."""
    def nothing(system, messages, **kw):
        yield StreamEvent(type="waiting", text="waiting 10 s")

    eng = _engine(agent_tree, nothing)
    eng.stream_response("hi")
    entry = _last(eng)
    assert entry["first_event_ms"] is None
    assert tm.silence_kind(entry) != "model"
