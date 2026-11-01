"""wave14/w04 RULE 1: no cost to other abilities.

The stable system-prompt head must be identical for three non-chemistry
tasks (Python bugfix, Rust feature, docs edit) and must not carry ORCA
knowledge, the new tool descriptions, or a chemistry reminder.
"""
from __future__ import annotations

import pytest

TASKS = (
    "fix the failing test in foo.py",
    "add a SIMD kernel to the Rust encoder crate",
    "clarify the installation section of README.md",
)


def _prompt(task_text: str, role: str, mode: str) -> str:
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    return loader.build_system_prompt(
        role_id=role, mode_id=mode, route=[role], task_text=task_text,
        session_key=f"identity-w04-{role}-{len(task_text)}")


@pytest.mark.parametrize("role,mode", [("solo_agent", "solo"),
                                       ("dashboard_agent", "dashboard")])
def test_prompts_identical_across_non_chemistry_tasks(role, mode):
    prompts = [_prompt(t, role, mode) for t in TASKS]
    assert prompts[0] == prompts[1] == prompts[2], (
        f"{role}: a non-chemistry task changed the prompt")


@pytest.mark.parametrize("role,mode", [("solo_agent", "solo"),
                                       ("dashboard_agent", "dashboard")])
def test_no_orca_knowledge_in_the_head(role, mode):
    text = _prompt(TASKS[0], role, mode)
    for marker in ("orca_explain", "orca_check_input",
                   "workflow_candidates", "3conv", "NoSplitRIJ",
                   "SLOWCONV", "TIGHTSCF"):
        assert marker not in text, (role, marker)
