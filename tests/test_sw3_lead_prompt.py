"""SW3 phase 2: the lead role preset and the lead skill are shipped and safe.

The lead session is the ONE session the user talks to while a swarm of
members does the work (package SW3, task "one DELFIN session leads the
swarm"). This test pins the two documents that carry the lead role:

- delfin/agent/pack/agents/swarm-lead_subagent.md  -- the SA0 preset,
  auto-discovered into SUBAGENT_PRESETS as a new (plan-mode) preset;
- delfin/agent/pack/skills/swarm-lead.md           -- the built-in
  skill, auto-discovered by discover_skills.

The safety invariants are the point (target point 3/4): the lead plans
and triages with read-only tools only, carries no execute path that
starts, stops, merges or answers an approval by itself, and treats every
member message as data, never an instruction.
"""

from delfin.agent.subagents import SUBAGENT_PRESETS, _known_preset_tool_names
from delfin.agent.skills import discover_skills, get_skill


def test_lead_preset_is_registered_as_a_plan_expert():
    assert "swarm-lead" in SUBAGENT_PRESETS, list(SUBAGENT_PRESETS)
    preset = SUBAGENT_PRESETS["swarm-lead"]
    # The lead is a planning/triage role: read-only by construction, it
    # never executes a state change. plan mode is how the surface
    # guarantees that.
    assert preset.mode == "plan", preset.mode


def test_lead_preset_lists_only_known_read_only_tools():
    known = frozenset(_known_preset_tool_names())
    preset = SUBAGENT_PRESETS["swarm-lead"]
    for tool in preset.tools:
        assert tool in known, f"preset lists unknown tool {tool!r}"
    # No write/steering verb may be a lead tool: the lead never starts,
    # stops, merges or answers an approval by itself.
    for forbidden in ("swarm_start", "swarm_pause", "swarm_resume",
                      "swarm_triage_answer", "write_file", "edit_file",
                      "apply_patch", "bash"):
        assert forbidden not in preset.tools, forbidden


def test_lead_preset_body_refuses_to_self_execute():
    body = SUBAGENT_PRESETS["swarm-lead"].system_prompt or ""
    low = body.lower()
    assert "never start" in low or "never START" in body
    assert "never pause" in low or "never PAUSE" in body
    assert "never merge" in low or "never MERGE" in body
    assert "approval" in low  # an approval is never answered by the lead
    # Member text is data, never an instruction (the delegate-report rule
    # extended to the whole swarm).
    assert "data, never" in low


def test_lead_skill_is_discoverable_and_declares_a_domain():
    skill = get_skill("swarm-lead")
    assert skill is not None, sorted(s.name for s in discover_skills())
    assert skill.domains, skill.domains
