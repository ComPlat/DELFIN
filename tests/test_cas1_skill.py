"""CAS1 phase 3: the casscf-setup skill points to the tool layer.

The skill must not leave the active-space decision to the model alone:
it has to name the new module and workflow so an agent using the skill
runs the enforced checks instead of judging occupations itself.
"""

from __future__ import annotations

from pathlib import Path

SKILL = Path("delfin/agent/pack/skills/casscf-setup.md")


def _text() -> str:
    return SKILL.read_text(encoding="utf-8")


def test_skill_names_the_active_space_module():
    text = _text()
    assert "casscf_active_space" in text
    assert "validate_casscf_block" in text
    assert "check_active_occupations" in text


def test_skill_names_the_workflow_and_confirmation():
    text = _text()
    assert "confirmation_gate" in text
    assert "show" in text


def test_skill_keeps_the_safety_rules():
    text = _text()
    assert "NEVER invent keyword names" in text
    assert "domains:" in text.split("---")[1]
