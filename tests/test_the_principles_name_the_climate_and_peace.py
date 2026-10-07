"""The principles' climate-and-peace commitment reaches the model.

Added 2026-10-07 at the maintainer's request. The principles already
committed to "the great problems of our time" and "protecting the planet
and its life"; this makes the orientation explicit: the climate crisis is
named among those problems, and the work serves these ends peacefully,
through sound science, for the common good and the living natural world.

Pinned by MEANING in the BUILT prompt, for every role — the same bar the
other principle commitments are held to. The agent also enforces the
file's digest at startup (principles_guard), so this guards the wording's
PRESENCE in what the model receives, while the guard protects the file
from silent change; the two are different questions.
"""

from __future__ import annotations

import pytest

from delfin.agent.prompt_loader import PromptLoader

_TERMS = ("climate crisis", "peacefully", "sound science",
          "common good", "living natural world")


@pytest.mark.parametrize("role, mode", [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
    ("office_agent", "office"),
])
def test_the_commitment_reaches_every_role(role, mode):
    prompt = PromptLoader().build_system_prompt(
        role_id=role, mode_id=mode, task_text="tidy the workspace")
    missing = [t for t in _TERMS if t not in prompt]
    assert not missing, (
        f"{role} no longer carries the climate/peace commitment: {missing}")


def test_it_is_an_orientation_not_a_ranking():
    """Deliberately NOT "foremost": a system prompt orients the work, it
    does not rank the planet's crises. The word was dropped in review on
    2026-10-07; this keeps it out so the choice is not quietly undone."""
    from pathlib import Path
    src = (Path(__file__).resolve().parents[1]
           / "delfin" / "agent" / "pack" / "shared"
           / "principles_addendum.md").read_text(encoding="utf-8")
    assert "climate crisis" in src
    assert "foremost" not in src.lower(), (
        "the principles rank the climate crisis above every other problem; "
        "that was considered and dropped — the prompt orients the work "
        "rather than ranking crises")
