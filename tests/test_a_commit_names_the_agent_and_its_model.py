"""A commit names the agent and the model that wrote it.

Asked for on 2026-09-15: the way Claude Code ends its commits with a
co-author line naming itself and its model, the DELFIN agent ends its own
with ``Co-Authored-By: DELFIN-Agent (<model>)``. On by default, off with
``agent.commit_coauthor = false``, and the user's own instructions or
remembered preferences win over it.
"""

from __future__ import annotations

import pathlib

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]


def _solo_prompt(model: str) -> str:
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = _ROOT
    return loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="commit the fix", session_key="attribution-test",
        model=model)


def test_the_prompt_names_the_model_in_the_trailer():
    prompt = _solo_prompt("kit.glm-5.3")
    assert "Co-Authored-By: DELFIN-Agent (kit.glm-5.3)" in prompt
    assert "remembered preferences win" in prompt


def test_the_setting_turns_it_off(monkeypatch):
    import delfin.user_settings as us
    monkeypatch.setattr(
        us, "load_settings", lambda *a, **k: {"agent": {"commit_coauthor": False}})
    assert "Co-Authored-By: DELFIN-Agent" not in _solo_prompt("kit.glm-5.3")
