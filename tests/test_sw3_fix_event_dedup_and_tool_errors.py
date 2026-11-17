"""SW3 phase-1 follow-up: two reviewer-confirmed defects, red first.

Adversarial controls that pin the two fixes on delfin/agent/swarm/lead.py:

F1 -- event loss in LeadEventFeed. The dedup key pinned
     (kind, package, member) collapses DISTINCT events of the same kind
     for the same member: a member's second ``phase_change`` (a different
     phase) and a second ``approval_waiting`` (for a different approval)
     are silently dropped, so the lead never learns about them. The dedup
     identity must include the payload's discriminator fields (phase, pr)
     so a retry of the SAME event dedups while a DISTINCT fact delivers.

F2 -- error contract. A lead tool must return a structured dict, never
     raise. Two escapes:
     - ``request_triage_answer`` lets SW1's ``TriageError`` (an
       escalation) propagate out of the tool instead of returning a
       structured refusal;
     - ``request_swarm_start`` calls ``plan_launch`` on its unconfirmed
       branch, so a bad spec raises ``SpecError`` instead of returning the
       safe dry-run/ask verdict. The unconfirmed branch is supposed to
       never crash -- it is the "present the plan and ask" path.
"""

from __future__ import annotations

import pytest

from delfin.agent.swarm import lead
from delfin.agent.swarm import status as swarm_status_mod


def event(kind: str, package: str = "core", member: str = "core-builder",
          **extra) -> dict:
    return {"kind": kind, "package": package, "member": member, **extra}


def write_events(feed: lead.LeadEventFeed, *events) -> list[str]:
    seen: list[str] = []
    feed.on_deliver = seen.append
    for e in events:
        swarm_status_mod.emit_event(e)
    return seen


class TestEventFeedDeliversDistinctEvents:
    """A retry of the same event dedups; a DISTINCT fact must not."""

    def test_two_phase_changes_of_one_member_both_deliver(self):
        # A member first enters 'building', then 'in_review'. Both are
        # separate facts the lead must see; the old (kind, package,
        # member) key dropped the second.
        clock = [0.0]
        feed = lead.LeadEventFeed(max_per_window=20, window_s=60.0,
                                  now=lambda: clock[0])
        off = feed.attach()
        try:
            seen = write_events(
                feed,
                event("phase_change", phase="building"),
                event("phase_change", phase="in_review"),
            )
            assert len(seen) == 2, seen
            assert any("phase=building" in s for s in seen), seen
            assert any("phase=in_review" in s for s in seen), seen
        finally:
            off()

    def test_two_approval_waiting_of_one_member_both_deliver(self):
        # Two distinct approvals queue for the same member (different
        # item ids are the only discriminator); both must deliver.
        clock = [0.0]
        feed = lead.LeadEventFeed(max_per_window=20, window_s=60.0,
                                  now=lambda: clock[0])
        off = feed.attach()
        try:
            seen = write_events(
                feed,
                event("approval_waiting", item="approval-1"),
                event("approval_waiting", item="approval-2"),
            )
            assert len(seen) == 2, seen
        finally:
            off()

    def test_ci_red_for_two_prs_both_deliver(self):
        clock = [0.0]
        feed = lead.LeadEventFeed(max_per_window=20, window_s=60.0,
                                  now=lambda: clock[0])
        off = feed.attach()
        try:
            seen = write_events(
                feed,
                event("ci_red", pr=611),
                event("ci_red", pr=612),
            )
            assert len(seen) == 2, seen
        finally:
            off()


class TestLeadToolsNeverRaise:
    """A lead tool returns a structured dict, never a traceback."""

    def test_triage_answer_refuses_an_escalation_as_data(self):
        # An escalation item must come back as a structured refusal, not
        # as SW1's TriageError escaping the tool.
        escalation = {
            "id": "a1",
            "kind": "session_message",
            "from": "outsider",          # not a roster member -> escalate
            "member": "core-builder",
            "roster": ["core-builder"],  # the swarm's member identities
        }
        result = lead.request_triage_answer(escalation, answer="nope",
                                            confirmed=True)
        assert isinstance(result, dict), result
        assert result.get("answered") is False
        assert "action" in result
        assert "refused" in result
        assert result["refused"]  # the reason names it goes to the user

    def test_triage_answer_unconfirmed_is_a_refusal(self):
        result = lead.request_triage_answer({"id": "x"}, answer="y")
        assert result["action"] == "ask_user"
        assert result.get("answered") is False

    def test_swarm_start_with_bad_spec_unconfirmed_returns_refusal(self):
        # The unconfirmed path must never crash on a bad spec: it is the
        # "present the plan and ask the user" path.
        bad = {"goal": "no packages declared"}
        result = lead.request_swarm_start(bad)
        assert isinstance(result, dict), result
        # Either an ask_user (good spec refused) or a structured refusal
        # (bad spec) -- never an exception, never a started member.
        assert not result.get("started")
        assert result.get("action") in ("ask_user", "refused")

    def test_triage_answer_confirmed_auto_passes_the_flag_through(self):
        # Sanity: a genuinely auto-answerable item still reaches SW1's
        # triage_answer (the real path), not a local re-implementation.
        # classify reads 'from', requires the member in the roster, and
        # auto-answers a test_output read by its own member.
        auto = {
            "id": "a2",
            "kind": "test_output",
            "from": "core-builder",
            "member": "core-builder",
            "roster": ["core-builder"],
        }
        result = lead.request_triage_answer(auto, answer="ok",
                                            confirmed=True)
        assert result.get("answered") is True, result
