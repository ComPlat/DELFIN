"""A mid-session nudge to keep memory — driven by work, never by count.

Package 7, phase 3. Hermes nudges the agent at intervals to persist
what it learned, instead of only at session end. The one hard rule
carried over from DELFIN's own history: the trigger is WORK (tool
calls, tokens produced since the last nudge), NEVER the number of
messages — a message-counter trigger once fired at 15% context fill.

These tests pin ``delfin.agent.memory_nudge``:

- enough work since the last nudge (tool calls AND estimated tokens
  both past their thresholds) -> the nudge text, which asks "is
  anything here worth keeping?" and names /memorize;
- not enough work -> "" ;
- many messages with no tool work -> "" (the anti message-count rule,
  stated as a test);
- at most ONE nudge per turn: the second call within the same turn
  returns "" even after more work;
- ``agent.memory_nudge: false`` in user settings switches it off;
- the decision function is pure w.r.t. its state object — it updates
  the state in place (last-nudge counters) and returns text or "".

Nothing here writes to disk or to the real ~/.delfin.
"""

from __future__ import annotations

import pytest


def _state(**kw):
    from delfin.agent import memory_nudge
    return memory_nudge.NudgeState(**kw)


def _maybe(state, *, tool_calls=0, chars=0, turn=None, settings=None):
    from delfin.agent import memory_nudge
    return memory_nudge.maybe_nudge(
        state, tool_calls=tool_calls, chars=chars, turn=turn,
        settings=settings)


def test_enough_work_triggers_the_nudge():
    from delfin.agent import memory_nudge as mn
    out = _maybe(_state(), tool_calls=mn.NUDGE_TOOL_CALLS + 1,
                 chars=mn.NUDGE_CHARS * 5)
    assert out and "worth keeping" in out
    assert "/memorize" in out


def test_too_little_work_no_nudge():
    from delfin.agent import memory_nudge as mn
    assert _maybe(_state(),
                  tool_calls=mn.NUDGE_TOOL_CALLS - 1,
                  chars=mn.NUDGE_CHARS * 5) == ""
    assert _maybe(_state(),
                  tool_calls=mn.NUDGE_TOOL_CALLS + 1,
                  chars=mn.NUDGE_CHARS - 1) == ""


def test_messages_without_tool_work_never_nudge():
    """The anti rule: a message-count trigger fired at 15% fill once.

    A thousand 'messages' (whatever they are) with no tool calls and
    no produced text is zero WORK — the nudge must stay silent.
    """
    from delfin.agent import memory_nudge as mn
    assert _maybe(_state(),
                  tool_calls=0, chars=0) == ""


def test_at_most_one_nudge_per_turn():
    from delfin.agent import memory_nudge as mn
    state = _state()
    first = _maybe(state, tool_calls=mn.NUDGE_TOOL_CALLS + 5,
                   chars=mn.NUDGE_CHARS * 10, turn=7)
    assert first
    # More work arrives in the SAME turn: still silent, and the
    # counters keep accumulating for the next turn's decision.
    second = _maybe(state, tool_calls=mn.NUDGE_TOOL_CALLS + 5,
                    chars=mn.NUDGE_CHARS * 10, turn=7)
    assert second == ""
    # A NEW turn with fresh work nudges again.
    third = _maybe(state, tool_calls=mn.NUDGE_TOOL_CALLS + 5,
                   chars=mn.NUDGE_CHARS * 10, turn=8)
    assert third


def test_setting_switches_the_nudge_off():
    from delfin.agent import memory_nudge as mn
    assert _maybe(_state(),
                  tool_calls=mn.NUDGE_TOOL_CALLS + 5,
                  chars=mn.NUDGE_CHARS * 10,
                  settings={"agent": {"memory_nudge": False}}) == ""


def test_state_tracks_work_since_last_nudge():
    from delfin.agent import memory_nudge as mn
    state = _state()
    _maybe(state, tool_calls=3, chars=200)
    assert state.tool_calls_since_nudge == 3
    assert state.chars_since_nudge == 200
    out = _maybe(state, tool_calls=mn.NUDGE_TOOL_CALLS,
                 chars=mn.NUDGE_CHARS * 5)
    assert out
    assert state.tool_calls_since_nudge == 0
    assert state.chars_since_nudge == 0
