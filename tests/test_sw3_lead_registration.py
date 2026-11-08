"""SW3 phase 3: lead tool names are registered so the lead session's
preset can list them (red control — see delfin/agent/swarm/lead.py for
the tool-name truth, api_client.py for where they become known).

The phase-2 preset (delfin/agent/pack/agents/swarm-lead_subagent.md)
lists only generic plan-mode tools; it cannot yet list the swarm_* lead
tools because those names are not part of _known_preset_tool_names().
This registration makes them known. The change lives in the PROTECTED
file delfin/agent/api_client.py (`_DELFIN_ONLY_TOOL_NAMES`), delivered
as .gate/px/orig_api_client.py + .gate/px/new_api_client.py and built
by the operator. This file is the red control: it must FAIL before the
registration exists and PASS after it is applied.
"""

import pytest


def _lead_tool_names() -> set[str]:
    """The lead tools the phase-2 preset must be able to list once the
    registration lands. Source of truth: lead.py READ_ONLY_TOOLS and
    STATE_CHANGING_TOOLS (lines 33-41)."""
    return {
        # read-only lead tools (safe in default mode)
        "swarm_status", "swarm_plan", "swarm_message", "swarm_triage_list",
        # state-changing lead tools (always ASK the user)
        "swarm_start", "swarm_pause", "swarm_resume", "swarm_triage_answer",
        "swarm_prcheck",
    }


def _known_preset_tool_names():
    from delfin.agent.subagents import _known_preset_tool_names as _k
    return _k()


@pytest.mark.parametrize(
    "tool",
    ["swarm_status", "swarm_plan", "swarm_message", "swarm_triage_list",
     "swarm_start", "swarm_pause", "swarm_resume", "swarm_triage_answer",
     "swarm_prcheck"],
)
def test_lead_tool_is_known_after_registration(tool):
    """Each lead tool name must be part of the known preset-tool surface
    once the operator applies the api_client registration."""
    assert tool in _known_preset_tool_names(), (
        f"{tool} is not known to _known_preset_tool_names(); the SW3 "
        f"phase-3 registration in api_client.py _DELFIN_ONLY_TOOL_NAMES "
        f"has not been applied yet."
    )


def test_known_surface_covers_the_whole_lead_toolset():
    """The union of the registered names must account for every lead
    tool that lead.py declares, so the preset can expose the lead's
    full documented surface."""
    known = _known_preset_tool_names()
    missing = sorted(_lead_tool_names() - set(known))
    assert not missing, f"lead tools not yet registered: {missing}"
