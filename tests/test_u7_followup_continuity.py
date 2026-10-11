"""U7 / phase 3 — the continuation fires even on a blocked-only task list.

s18 regression guard: a turn that ended "Let me …" used to get NO follow-up
turn when the task list held only blocked tasks, leaving the session idle
until poked. The root cause was that the detector first ASKED whether
non-blocked open work remained and only then checked the announcement; a
blocked-only list read as "nothing to do" and the follow-up was suppressed.

The fix (turn_continuation, committed fdf55623): ``should_continue`` decides
purely on the announcement — the open-task count is carried into the note but
no longer gates it — so an announced action fires a continuation for ANY open
list, blocked-only included. This file PINS that behaviour so the s18 hole
cannot be re-introduced by a later count-gate. It is a regression guard, not
a red-to-green fix: it is green on this base.
"""

from __future__ import annotations

from delfin.agent.turn_continuation import should_continue

_ANNOUNCED = "Let me write the failing test first."
_HANDBACK = "I'll wait for your reply."

#: A count that could only be made up of blocked tasks has no actionable
#: (pending / in_progress) step — the exact shape that suppressed the
#: follow-up before the fix.
_BLOCKED_ONLY_SCALE = 3


def test_an_announced_action_fires_on_a_blocked_only_list():
    ok, note = should_continue(_ANNOUNCED, _BLOCKED_ONLY_SCALE)
    assert ok is True
    assert "Continue the announced work now" in note


def test_it_also_fires_when_the_caller_reports_no_open_tasks():
    ok, _ = should_continue(_ANNOUNCED, 0)
    assert ok is True


def test_the_note_tells_a_blocked_agent_what_to_do_instead_of_nagging():
    _, note = should_continue(_ANNOUNCED, _BLOCKED_ONLY_SCALE)
    # A blocked-only list is exactly the case where the follow-up needs the
    # "if it cannot proceed, mark the tasks blocked and say what is missing"
    # instruction — otherwise the continued turn only re-announces.
    assert "mark the tasks blocked" in note


def test_a_handback_is_not_forced_to_continue_even_on_a_blocked_list():
    ok, note = should_continue(_HANDBACK, _BLOCKED_ONLY_SCALE)
    assert ok is False
    assert note == ""


def test_the_return_shape_is_stable_for_every_open_count():
    for n in (0, 1, 3, 10):
        out = should_continue(_ANNOUNCED, n)
        assert isinstance(out, tuple) and len(out) == 2
