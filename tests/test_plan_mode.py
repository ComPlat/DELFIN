"""Tests for plan-mode wiring: addendum injection + dropdown surfacing.

Plan-mode is solo_agent role + permission_profile="plan" (read-only)
+ a markdown addendum that tells the model to use ExitPlanMode for
approval. The addendum file lives at
``delfin/agent/pack/shared/plan_mode_addendum.md`` and the prompt
loader picks it up only when ``mode_id == "plan"``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader


def test_plan_mode_addendum_file_exists():
    p = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack" / "shared" / "plan_mode_addendum.md"
    assert p.is_file(), f"missing addendum: {p}"
    body = p.read_text(encoding="utf-8")
    assert "Plan Mode" in body
    assert "ExitPlanMode" in body


def test_solo_prompt_includes_plan_addendum_when_mode_is_plan():
    loader = PromptLoader()
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="plan",
        mode_description="plan mode test",
        task_text="figure out how to add feature X",
    )
    assert "Plan Mode" in prompt
    assert "ExitPlanMode" in prompt


def test_solo_prompt_skips_plan_addendum_in_other_modes():
    loader = PromptLoader()
    for mode in ("solo", "dashboard", "quick"):
        prompt = loader.build_system_prompt(
            role_id="solo_agent",
            mode_id=mode,
            mode_description=f"{mode} mode test",
            task_text="anything",
        )
        # The addendum's distinctive heading must not leak into non-plan modes
        assert "# Plan Mode\n" not in prompt, mode


def test_plan_mode_addendum_documents_exit_plan_mode_handoff():
    p = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack" / "shared" / "plan_mode_addendum.md"
    body = p.read_text(encoding="utf-8")
    # Key contract elements must be present
    assert "ExitPlanMode" in body
    assert "approve" in body.lower()
    assert "acceptEdits" in body
