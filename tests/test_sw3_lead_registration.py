"""SW3 phase 3 / v4: the lead tool names are KNOWN to the preset-tool
surface through the LEAD-ONLY set, not _DELFIN_ONLY_TOOL_NAMES.

v4 contract (operator, 2026-11-12): the nine swarm_* lead tools are not
advertised tools of the normal agent -- putting them into
_DELFIN_ONLY_TOOL_NAMES broke the schema budget, the catalogue
measurement, the DELFIN-workspace filter and the office denylist. They
now live in _SWARM_LEAD_TOOL_NAMES (delfin/agent/api_client.py), which
feeds only the typo check (unknown_tool_names) and the T4 grounding
_known set. _known_preset_tool_names() unions those known tools, so the
SA0 lead preset can list them.

This file is the red control: it must FAIL before the operator build
applies the v4 api_client copy and PASS after it.
"""

import pytest


def _lead_tool_names() -> set[str]:
    """The lead tools the SA0 preset must be able to list once the v4
    registration lands. Source of truth: lead.py READ_ONLY_TOOLS and
    STATE_CHANGING_TOOLS. The SA0 preset is a PLANNING role (phase-2
    safety design, tests/test_sw3_lead_prompt.py): it lists only the
    read-only lead tools; the state-changing ones are the lead
    session's own surface, confirmed with the user there."""
    return {
        # read-only lead tools (safe in default mode)
        "swarm_status", "swarm_plan", "swarm_message", "swarm_triage_list",
        # state-changing lead tools (always ASK the user)
        "swarm_start", "swarm_pause", "swarm_resume", "swarm_triage_answer",
        "swarm_prcheck",
    }


def _preset_lead_tool_names() -> set[str]:
    """What the SA0 preset lists: the READ-ONLY lead tools only."""
    return {"swarm_status", "swarm_plan", "swarm_triage_list"}


def _known_preset_tool_names():
    from delfin.agent.subagents import _known_preset_tool_names as _k
    return _k()


def _delfin_only_tool_names():
    from delfin.agent.api_client import _DELFIN_ONLY_TOOL_NAMES
    return set(_DELFIN_ONLY_TOOL_NAMES)


@pytest.mark.parametrize(
    "tool",
    ["swarm_status", "swarm_plan", "swarm_triage_list"],
)
def test_lead_tool_is_known_after_registration(tool):
    """Each READ-ONLY lead tool must be part of the known preset-tool
    surface once the operator applies the v4 api_client copy. The SA0
    preset is a planning role: it lists the read-only lead tools only
    (phase-2 safety design); the state-changing ones are the lead
    session's own confirmed surface."""
    assert tool in _known_preset_tool_names(), (
        f"{tool} is not known to _known_preset_tool_names(); the SW3 "
        f"v4 lead registration (_SWARM_LEAD_TOOL_NAMES feeding the "
        f"unknown/grounding known sets) has not been applied yet."
    )


def test_known_surface_covers_the_whole_lead_toolset():
    """The union of the registered names must account for every lead
    tool that lead.py declares, so the preset can expose the lead's
    full documented surface."""
    known = _known_preset_tool_names()
    missing = sorted(_preset_lead_tool_names() - set(known))
    assert not missing, f"lead tools not yet registered: {missing}"


def test_lead_tools_stay_out_of_the_delfin_only_set():
    """v4 contract: the lead tools are NOT in _DELFIN_ONLY_TOOL_NAMES --
    they are not advertised tools of the normal agent, and the schema
    budget / catalogue measurement / workspace filter / office denylist
    all key off that set."""
    leaked = sorted(_lead_tool_names() & _delfin_only_tool_names())
    assert not leaked, (
        f"lead tools leaked into _DELFIN_ONLY_TOOL_NAMES: {leaked} -- "
        "they must live only in _SWARM_LEAD_TOOL_NAMES (v4 contract)")
