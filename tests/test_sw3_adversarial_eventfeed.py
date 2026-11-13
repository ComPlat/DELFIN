"""Adversarial tests for SW3 phase 1 (reviewer s21).

Target: delfin/agent/swarm/lead.py LeadEventFeed deduplication.

The task (package SW3, target 2) requires that "member finished a
phase" and "approval waiting" events are DELIVERED to the lead so it
reacts to them. The feed's dedup key is (kind, package, member) -- it
omits the event's own payload (phase, item, ...). A member that moves
through several phases therefore has every phase change AFTER the first
silently dropped, because they all share one dedup key. A second,
distinct event for the same member is not a retry of the first -- it
is new information the lead must see.
"""

from delfin.agent.swarm import lead
from delfin.agent.swarm import status as swarm_status_mod


def _feed_and_watch(**kw):
    feed = lead.LeadEventFeed(**kw)
    seen = []
    feed.on_deliver = seen.append
    off = feed.attach()
    return feed, seen, off


def test_two_distinct_phase_changes_for_one_member_are_both_delivered():
    """A member moving building -> in_review -> qs_ok raises a phase
    change per transition. Each is new information; the feed must
    deliver all of them, not collapse them onto the first key."""
    feed, seen, off = _feed_and_watch()
    try:
        swarm_status_mod.emit_event(
            {"kind": "phase_change", "package": "core",
             "member": "core-builder", "phase": "building"})
        swarm_status_mod.emit_event(
            {"kind": "phase_change", "package": "core",
             "member": "core-builder", "phase": "in_review"})
        assert len(seen) == 2, seen
        # The second delivery carries the NEW phase, not a stale echo.
        assert "in_review" in seen[1], seen
    finally:
        off()


def test_two_separate_approval_waits_for_one_member_are_both_delivered():
    """A distinct approval wait (a second queued item for the same
    member) is not a retry of the first; the lead must see both."""
    feed, seen, off = _feed_and_watch()
    try:
        swarm_status_mod.emit_event(
            {"kind": "approval_waiting", "package": "core",
             "member": "core-builder", "item": "perm"})
        swarm_status_mod.emit_event(
            {"kind": "approval_waiting", "package": "core",
             "member": "core-builder", "item": "push"})
        assert len(seen) == 2, seen
    finally:
        off()
