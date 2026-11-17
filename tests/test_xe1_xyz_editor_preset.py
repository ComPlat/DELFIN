"""XE1: xyz-editor preset validates, and the ops-server geometry tools become
known to validate_preset (Path C, operator).

subagents._known_preset_tool_names must DERIVE the ops-server's structure /
geometry tool names (bare + mcp__delfin-ops__namespaced) from the server
module, and must keep refusing an unknown name. The acceptance test below is
RED on this branch until the operator applies the .gate/px subtitle patch;
the bogus-name test is a guard that must stay green (it would catch a
catalogue that started accepting arbitrary names).
"""
from __future__ import annotations

from pathlib import Path

from delfin.agent.subagents import validate_preset

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "agents"


def _preset_problems() -> list[str]:
    return validate_preset(str(_PACK / "xyz-editor_subagent.md"))


def test_xyz_editor_preset_passes_validate_preset_once_geometry_is_known():
    """The real xyz-editor preset (which lists the ops-server geometry tools)
    validates once the catalogue knows them (operator Path C patch). RED on
    this branch until that patch lands; GREEN after."""
    problems = _preset_problems()
    assert problems == [], problems


def test_xyz_editor_refuses_a_made_up_geometry_tool():
    """Guard (operator): a made-up mcp__delfin-ops__geometry_bogus is ALWAYS
    refused — the catalogue may never accept an arbitrary name. Green now and
    after the Path C patch."""
    md = """---
name: xyz-editor
description: edits xyz files
mode: default
tools: [mcp__delfin-ops__geometry_bogus]
write_globs: [structure_edit_output/**]
may_delegate_to: []
department: structure
enabled_setting: agent.departments.structure
---
body
"""
    problems = validate_preset(md)
    assert any("geometry_bogus" in p for p in problems), problems


def test_xyz_editor_preset_defaults_to_off():
    """The structure department ships behind enabled_setting, default OFF:
    nothing new is offered to the model until the person turns it on."""
    from delfin.agent.skills import _parse_frontmatter

    meta, _ = _parse_frontmatter(
        Path(_PACK / "xyz-editor_subagent.md").read_text(encoding="utf-8")
    )
    assert meta.get("enabled_setting") == "agent.departments.structure"
    # default OFF: the meta carries only the setting name, and the registry
    # does not load a department whose setting is off.
    assert "enabled" not in meta or meta.get("enabled", "").lower() != "true"
