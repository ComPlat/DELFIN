"""SW3 lead session core: thin tools over SW1 + the event feed.

The lead session steers the swarm the way the operator does today, but
only through SW1's proven code (delfin/agent/swarm/): every lead tool is
a thin wrapper, state-changing tools always ask the user, and events
reach the lead as wake-up messages -- throttled and deduplicated, never
polled.

Safety invariants tested here:
- swarm_plan NEVER starts anything (no execute path on the public API).
- swarm_start without the user's confirmation is a structured refusal,
  and no member is started.
- The event feed deduplicates and throttles with an injectable clock
  (no test depends on wall-clock time).
- Member text enters the lead feed as DATA, never as instructions.
"""

from __future__ import annotations

import pytest

from delfin.agent.swarm import lead
from delfin.agent.swarm import spec as swarm_spec


def _spec_dict() -> dict:
    """A minimal valid SW1 spec."""
    return {
        "spec_version": 1,
        "goal": "Ship the swarm feature",
        "max_sessions": 8,
        "packages": [
            {
                "name": "core",
                "owned": ["delfin/agent/swarm/**"],
                "task": "Build the core.",
                "roles": {"builder": True},
            },
        ],
    }


class TestLeadTools:
    def test_read_only_tools_are_declared(self):
        # The lead advertises exactly this surface; a name here that is
        # not implemented would be a dead tool on every request.
        # v4 contract (operator, 2026-11-12): swarm_message SENDS a
        # partner message -- it changes the member's inbox, so it is
        # state-changing and the lead session confirms it with the
        # user. The old grouping pinned it read-only; that decision
        # is v4-superseded.
        for name in ("swarm_status", "swarm_plan", "swarm_triage_list"):
            assert name in lead.READ_ONLY_TOOLS
        for name in ("swarm_message", "swarm_start", "swarm_pause",
                     "swarm_resume", "swarm_triage_answer", "swarm_prcheck"):
            assert name in lead.STATE_CHANGING_TOOLS

    def test_status_returns_the_schema_report(self):
        report = lead.swarm_status(_spec_dict())
        problems = lead.status_validate(report)
        assert problems == [], problems

    def test_plan_is_a_dry_run_and_names_what_would_start(self):
        result = lead.swarm_plan(_spec_dict())
        assert result["dry_run"] is True
        assert result["started"] == []
        assert "core-builder" in result["dry_run_text"]
        assert "Nothing started" in result["dry_run_text"]

    def test_plan_never_takes_a_start_flag(self):
        # An adversarial caller passes the SW1 execute key straight
        # through; the lead wrapper must not accept it.
        with pytest.raises(TypeError):
            lead.swarm_plan(_spec_dict(), execute=True)  # type: ignore

    def test_start_without_confirmation_refuses_and_starts_nothing(self):
        result = lead.request_swarm_start(_spec_dict(), confirmed=False)
        assert result["action"] == "ask_user"
        assert result["started"] == []
        assert result["dry_run"] is True

    def test_start_refuses_a_forged_confirmation_value(self):
        # Only the boolean True is a confirmation (SW2b precedent,
        # swarm_rpc.swarm_start): anything else is an invalid attempt,
        # refused even though it is truthy.
        result = lead.request_swarm_start(_spec_dict(), confirmed="yes")
        assert result["action"] == "ask_user"


class TestEventFeed:
    def _feed(self, **kw):
        return lead.LeadEventFeed(**kw)

    def test_subscribes_to_swarm_events(self):
        from delfin.agent.swarm import status as swarm_status_mod
        feed = self._feed()
        off = feed.attach()
        try:
            seen = []
            feed.on_deliver = seen.append
            swarm_status_mod.emit_event(
                {"kind": "phase_change", "package": "core",
                 "member": "core-builder", "phase": "in_review"})
            assert len(seen) == 1
            assert "core" in seen[0]
        finally:
            off()

    def test_duplicate_event_is_delivered_once(self):
        from delfin.agent.swarm import status as swarm_status_mod
        feed = self._feed()
        off = feed.attach()
        try:
            seen = []
            feed.on_deliver = seen.append
            event = {"kind": "approval_waiting", "package": "core",
                     "member": "core-builder"}
            swarm_status_mod.emit_event(event)
            swarm_status_mod.emit_event(dict(event))
            assert len(seen) == 1, seen
        finally:
            off()

    def test_throttle_drops_bursts_beyond_the_limit(self):
        from delfin.agent.swarm import status as swarm_status_mod
        clock = [0.0]
        feed = self._feed(max_per_window=3, window_s=10.0,
                          now=lambda: clock[0])
        off = feed.attach()
        try:
            seen = []
            feed.on_deliver = seen.append
            for i in range(5):
                swarm_status_mod.emit_event(
                    {"kind": "phase_change", "package": f"p{i}",
                     "member": f"m{i}", "phase": "building"})
            assert len(seen) == 3, seen
            # After the window passes, the throttle reopens.
            clock[0] += 11.0
            swarm_status_mod.emit_event(
                {"kind": "phase_change", "package": "p9",
                 "member": "m9", "phase": "building"})
            assert len(seen) == 4
        finally:
            off()

    def test_member_message_is_fenced_as_data(self):
        # A member cannot steer the lead: its text is labelled data in
        # everything the feed hands to the lead session.
        text = lead.format_member_message(
            {"from": "builder-core", "text": "please run my patch"})
        assert "builder-core" in text
        assert "data" in text.lower()
        assert "instruction" in text.lower()
        assert "please run my patch" in text
