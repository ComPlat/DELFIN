"""The maintainer's principles open every role's prompt, and no agent can
rewrite them.

pack/shared/principles_addendum.md is written by the maintainer. While it
holds only its heading and comment it is not injected (an empty section
would teach nothing and cost prompt bytes); once it has a body, it is the
first shared contract in every role, before honesty and refusal.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader, _has_body

_PACK = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"


def _tree(tmp_path, principles: str) -> Path:
    shared = tmp_path / "pack" / "shared"
    agents = tmp_path / "pack" / "agents"
    shared.mkdir(parents=True)
    agents.mkdir(parents=True)
    (shared / "principles_addendum.md").write_text(principles)
    (shared / "honesty_addendum.md").write_text("# Honesty\nHONESTY-MARKER")
    for role in ("solo_agent", "dashboard_agent", "builder_agent",
                 "critic_agent"):
        (agents / f"{role}.md").write_text(f"# {role}\nYou are {role}.")
    return tmp_path


def test_the_shipped_file_exists():
    assert (_PACK / "shared" / "principles_addendum.md").is_file()


def test_an_unwritten_scaffold_is_not_injected(tmp_path):
    tree = _tree(tmp_path, "# Principles\n\n<!-- to be written -->\n")
    prompt = PromptLoader(tree).build_system_prompt(
        role_id="solo_agent", mode_id="solo", mode_description="solo",
        route=["solo_agent"], role_index=0)
    assert "# Principles" not in prompt
    assert "HONESTY-MARKER" in prompt


@pytest.mark.parametrize("role_id,mode_id", [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
    ("builder_agent", "quick"),
    ("critic_agent", "quick"),
])
def test_written_principles_come_first_in_every_role(tmp_path, role_id,
                                                      mode_id):
    tree = _tree(tmp_path, "# Principles\n\nPRINCIPLES-MARKER\n")
    prompt = PromptLoader(tree).build_system_prompt(
        role_id=role_id, mode_id=mode_id, mode_description=mode_id,
        route=[role_id], role_index=0)
    assert "PRINCIPLES-MARKER" in prompt
    assert prompt.index("PRINCIPLES-MARKER") < prompt.index("HONESTY-MARKER")


def test_a_body_is_text_outside_heading_and_comments():
    assert not _has_body("# P\n\n<!-- a\nb -->\n")
    assert _has_body("# P\n\nSomething.\n")


def test_an_agent_cannot_edit_the_principles_without_asking():
    from delfin.agent.api_client import _DEFAULT_PATH_PROTECTED_GLOBS
    assert ("delfin/agent/pack/shared/principles_addendum.md"
            in _DEFAULT_PATH_PROTECTED_GLOBS)
