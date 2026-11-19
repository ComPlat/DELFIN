"""Adversarial round 6 (reviewer nacht-s21): the lead's own tool
classification must match the v4 contract.

The v4 api_client pair (operator, 2026-11-12) puts swarm_message in the
STATE-CHANGING group: "swarm_message sends a partner message -- it is
not read-only". The lead session confirms a message send with the user
before dispatch. But lead.py:33-36 still declares swarm_message in
READ_ONLY_TOOLS ("safe in the lead session's default mode"), and
test_sw3_lead_core.py pins that old grouping. When the lead session's
own mode logic is wired (protected files, later build), a consumer that
reads lead.READ_ONLY_TOOLS would send partner messages without the
confirmation the v4 contract requires. The lead.py grouping is the
single source the lead session will key off -- it must not contradict
the adjudicated classification.
"""
from delfin.agent.swarm import lead


def test_swarm_message_is_not_declared_read_only():
    """v4 contract: a message send is state-changing; it must live in
    STATE_CHANGING_TOOLS and never in READ_ONLY_TOOLS."""
    assert "swarm_message" not in lead.READ_ONLY_TOOLS, (
        "swarm_message is declared READ-ONLY in lead.py but the v4 "
        "contract classifies a message send as state-changing -- the "
        "lead session's mode logic would skip the user confirmation")
    assert "swarm_message" in lead.STATE_CHANGING_TOOLS, (
        "swarm_message must be in STATE_CHANGING_TOOLS: the lead "
        "confirms a partner-message send with the user before dispatch")
