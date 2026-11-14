"""Adversarial tests, round 3: LeadEventFeed identity and detach.

Reviewer nacht-s21, package SW3. Attacks on delfin/agent/swarm/lead.py:

1. Unicode look-alike payloads: _event_identity folds
   json.dumps(rest, sort_keys=True) into the key. Two events whose
   payloads differ only by look-alike characters (Cyrillic 'a' vs Latin
   'a') must be DELIVERED SEPARATELY -- over-broad key folding (the F1
   defect the builder fixed once already) would swallow one of them.
   format_event renders only kind/package/member/phase/pr, so the
   assertions read member= from the rendered text.
   Pinned green on the current code.
2. detach() must stop delivery through the REAL bus path: feed.attach()
   subscribes to status.emit_event (lead.py:242 -> status.py:107);
   emit_event delivers to every subscriber (status.py:92); after
   detach() a new emit_event must not reach on_deliver, and a second
   detach() is a no-op, not an error.
"""
from delfin.agent.swarm import lead as lead_mod
from delfin.agent.swarm import status as status_mod


def _event(kind, member, note):
    return {"kind": kind, "member": member, "note": note}


def test_unicode_lookalike_payloads_are_distinct_events():
    """Cyrillic-a and Latin-a member names must both be delivered."""
    delivered = []
    feed = lead_mod.LeadEventFeed(on_deliver=delivered.append)
    feed._on_event(_event("phase_change", "m\u04301", "cyrillic a"))
    feed._on_event(_event("phase_change", "ma1", "latin a"))
    assert len(delivered) == 2, delivered
    texts = "\n".join(delivered)
    # format_event renders member= but drops the 'note' key: assert on
    # the rendered member values, not the note.
    assert "member=m\u04301" in texts
    assert "member=ma1" in texts


def test_detach_stops_delivery_and_is_idempotent():
    """attach -> emit delivered; detach -> silent; detach again -> no error."""
    delivered = []
    feed = lead_mod.LeadEventFeed(on_deliver=delivered.append)
    off = feed.attach()
    try:
        status_mod.emit_event(_event("session_stalled", "m1", "first"))
        assert len(delivered) == 1, delivered
    finally:
        off()
        feed.detach()  # second unsubscribe must be a no-op, not an error
        feed.detach()
    # After a real detach, an emit_event on the module-level bus must not
    # reach this feed. A distinct event, so _seen does not dedup it.
    status_mod.emit_event(_event("session_stalled", "m1", "after detach"))
    assert len(delivered) == 1, delivered
