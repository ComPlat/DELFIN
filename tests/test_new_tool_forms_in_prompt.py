"""Controls that the new tool forms are named in the built prompt.

Red on the previous commit: the reading-tools section of
solo_agent.md still steered the agent to offset/limit reads and said
nothing about append — so the patterns Welle 4 counted (sed -n
"$(grep …)" and cat >> file <<EOF) had no in-prompt alternative.
After the edit the three new forms (around=, mode="append") must be
present BOTH in the source markdown and in the prompt the loader
actually builds for a solo session, and the token budget of the
source file must still hold (that one is also pinned in
tests/test_prompt_token_budget.py; this file re-checks it because the
section edit sits exactly on that budget's edge).
"""

from __future__ import annotations

import re
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_SOLO = _REPO / "delfin" / "agent" / "pack" / "agents" / "solo_agent.md"


def _markers_stripped(text: str) -> str:
    return re.sub(r"^<!--\s*module:[a-zA-Z0-9_-]+\s*-->\s*$\n?", "",
                  text, flags=re.M)


def test_solo_prompt_names_the_around_read():
    text = _SOLO.read_text()
    assert "around=<pattern>" in text


def test_solo_prompt_names_the_append_mode():
    text = _SOLO.read_text()
    assert 'mode="append"' in text
    assert "cat >> file <<EOF" in text  # named as the thing NOT to do


def test_the_old_offset_recipe_is_gone():
    # The section must not still tell the agent to use offset/limit as
    # the way to read around a grep hit -- that advice is what drove
    # the $( ) idiom. offset/limit itself stays valid elsewhere.
    text = _SOLO.read_text()
    assert "`read_file` with offset/limit" not in text


def test_built_prompt_carries_the_new_forms(monkeypatch):
    from delfin import user_settings
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(
        user_settings, "load_settings",
        lambda *a, **k: {"agent": {"slim_prompt": True}})
    built = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="fix the failing test in foo.py",
        session_key="new-forms-1")
    assert "around=<pattern>" in built
    assert 'mode="append"' in built


def test_solo_budget_still_holds():
    # 10735 is solo_agent.md's pinned budget in
    # tests/test_prompt_token_budget.py; the section rewrite was sized
    # to be budget-neutral, and this guard keeps a future re-edit of
    # this exact section honest.
    tokens = (len(_markers_stripped(_SOLO.read_text())) + 3) // 4
    assert tokens <= 10735, f"{tokens} tokens (>10735). Trim first."
