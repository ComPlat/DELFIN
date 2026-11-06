"""LK1 — the DSML re-request policy as a pure, bounded state object.

The 06.11 defect: a wave builder wrote its ``git commit`` as text twice in
one user turn. The first leak was re-requested and ran; the second leak in
the SAME user turn hit the one-shot latch (``_dsml_leak_spent`` in
``AgentEngine``) and the turn ended — the session sat idle for ~4 hours with
uncommitted work. One re-request per user turn is too strict for a long
autonomous turn with dozens of tool calls; the latch exists to stop a loop,
not to stop recovery after real progress.

This package moves the re-request *policy* into a small, pure, tested state
object (``delfin.agent.leak_recovery.DSMLReplayPolicy``):

* the FIRST re-request of a user turn is always allowed when the leaked
  call dominates the answer — exactly today's behaviour;
* a SECOND (or later) re-request in the same turn is allowed ONLY when at
  least one REAL tool call has completed since the last re-request
  (progress) AND fewer than ``agent.leak_recovery.max_per_turn``
  re-requests happened this turn (default 3);
* two leaks back-to-back without a real tool call between them are
  annotated and stopped, exactly as today (no loop).

All bounds live here and are tested here; the protected engine.py / 
api_client.py hooks only *call* these decisions.

Red control: the module does not exist yet, so importing it fails.
"""
from __future__ import annotations

from delfin.agent.leak_recovery import DSMLReplayPolicy

DEFAULT_MAX = 3


def test_first_rerequest_is_allowed_even_without_prior_progress():
    """The very first re-request of a turn must succeed on a dominant leak —
    this is today's single-leak behaviour preserved."""
    p = DSMLReplayPolicy()
    assert p.allowed_for(dominates=True) is True
    assert p.used == 0
    p.spend()
    assert p.used == 1


def test_first_rerequest_requires_dominance():
    """A non-dominant answer (prose, instructive markup) is never a
    re-request, matching leaked_tool_call_dominates."""
    p = DSMLReplayPolicy()
    assert p.allowed_for(dominates=False) is False
    assert p.used == 0  # a refused re-request is not spent


def test_second_rerequest_needs_a_real_tool_call_in_between():
    """The core of the fix: after the first re-request, a second one is
    allowed only when a real tool call completed since."""
    p = DSMLReplayPolicy()
    p.spend()                       # leak A re-requested
    # no real tool call since -> a second leak is NOT re-requested
    assert p.allowed_for(dominates=True) is False
    # a real tool call now provides the progress
    p.note_real_tool_call()
    assert p.allowed_for(dominates=True) is True


def test_two_leaks_in_a_row_annotate_and_stop():
    """leak A -> re-request -> leak A again (no real call between): stop,
    exactly as today. This is the never-loops core."""
    p = DSMLReplayPolicy()
    assert p.allowed_for(True) is True
    p.spend()
    assert p.allowed_for(True) is False  # no progress -> no second re-request


def test_real_tool_call_progress_consumed_by_the_next_rerequest():
    """Once a second re-request is granted, that progress is spent: another
    leak right after it (without a new real call) is again refused."""
    p = DSMLReplayPolicy()
    p.spend()                       # 1
    p.note_real_tool_call()
    assert p.allowed_for(True) is True
    p.spend()                       # 2 consumes the progress
    assert p.allowed_for(True) is False  # nothing ran since -> refused


def test_the_turn_is_bounded_by_max_per_turn():
    """A model that keeps making a real call between leaks is still bounded:
    after max_per_turn re-requests the turn stops no matter the progress."""
    p = DSMLReplayPolicy()          # max 3
    for _ in range(DEFAULT_MAX):
        p.spend()
        p.note_real_tool_call()
    assert p.used == DEFAULT_MAX
    # even with fresh progress, the per-turn bound refuses
    p.note_real_tool_call()
    assert p.allowed_for(True) is False


def test_max_per_turn_is_configurable():
    p = DSMLReplayPolicy(max_per_turn=1)
    assert p.allowed_for(True) is True
    p.spend()
    p.note_real_tool_call()
    assert p.allowed_for(True) is False  # bound of 1 reached


def test_exhausted_reports_when_no_more_rerequests_are_left():
    p = DSMLReplayPolicy()
    for _ in range(DEFAULT_MAX):
        p.spend()
    assert p.exhausted is True
    assert p.used == DEFAULT_MAX


def test_note_real_tool_call_alone_grants_nothing_on_the_first():
    """Progress alone never grants a re-request on a non-dominant answer."""
    p = DSMLReplayPolicy()
    p.note_real_tool_call()
    assert p.allowed_for(dominates=False) is False
