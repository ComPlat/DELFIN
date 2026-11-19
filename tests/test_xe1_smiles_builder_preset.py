"""XE1: smiles-builder preset validates, refuses a made-up tool, and ships
default OFF behind the structure-department setting.

The smiles-builder expert reviews and builds metal-complex SMILES in
DELFIN's standard (MANTA) convention: N directly bonded to a metal is an
amido ligand, O[M] is hydroxo, an ammine is written [NH3] or N->[M] as a
dative charge-separated bond. The preset never rewrites the user's SMILES;
it FLAGS nonstandard forms and only builds a 3D structure from a corrected
SMILES when asked.
"""
from __future__ import annotations

from pathlib import Path

from delfin.agent.subagents import validate_preset

_PACK = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "agents"


def _preset_problems() -> list[str]:
    return validate_preset(str(_PACK / "smiles-builder_subagent.md"))


def test_smiles_builder_preset_passes_validate_preset():
    """The real smiles-builder preset (which lists known read + geometry
    build tools) validates with no problems."""
    problems = _preset_problems()
    assert problems == [], problems


def test_smiles_builder_refuses_a_made_up_tool():
    """Guard: a made-up tool name is always refused — the preset may only
    list known read/geometry/report tools."""
    md = """---
name: smells-builder
description: builds SMILES
mode: default
tools: [mcp__delfin-ops__geometry_smell_bad]
write_globs: []
may_delegate_to: []
department: structure
enabled_setting: agent.departments.structure
---
body
"""
    problems = validate_preset(md)
    assert any("geometry_smell_bad" in p for p in problems), problems


def test_smiles_builder_defaults_to_off():
    """The structure department ships behind enabled_setting, default OFF:
    nothing new is offered to the model until the person turns it on."""
    from delfin.agent.skills import _parse_frontmatter

    meta, _ = _parse_frontmatter(
        Path(_PACK / "smiles-builder_subagent.md").read_text(encoding="utf-8")
    )
    assert meta.get("enabled_setting") == "agent.departments.structure"
    assert "enabled" not in meta or meta.get("enabled", "").lower() != "true"
