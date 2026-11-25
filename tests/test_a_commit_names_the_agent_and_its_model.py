"""A commit names the agent and the model that wrote it.

Asked for on 2026-09-15: a commit an agent makes ends with a co-author
line naming the agent and the model behind it,
``Co-Authored-By: DELFIN-Agent (<model>)``. On by default, off with
``agent.commit_coauthor = false``, and the user's own instructions or
remembered preferences win over it.

What this file checks is the PROMPT: that the model is told to write the
line. The mechanical half -- the rewrite that appends it whatever the
model does, and the call that arms it -- is in
tests/test_a_commit_carries_the_model_that_wrote_it.py, because for
months the prompt was the only half there was: nothing called
``set_commit_trailer_model``, so the rewrite never fired and the line
appeared only when a model remembered to type it.
"""

from __future__ import annotations

import pathlib

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def _an_ordinary_session(monkeypatch):
    """Every case here documents the ORDINARY session's prompt.

    The block a session is given is decided (since the w48 landing rules) by
    ``commit_attribution.commit_coauthor_for``, which reads
    ``DELFIN_SWARM_RUN``: run this file from inside a wave without pinning it
    and the five prompt cases below fail -- not because the attribution broke
    but because a member is deliberately NOT given the line. The pin says
    which session each case is about; the member case is its own test at the
    bottom. Before the member decision existed the code path did not read the
    environment, so the omission of this fixture measured nothing.
    """
    monkeypatch.delenv("DELFIN_SWARM_RUN", raising=False)


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


def test_the_address_is_where_the_credit_is_decided(monkeypatch):
    """GitHub attaches a co-author to an account by its email. The
    default belongs to none on purpose — a deployment that wants the
    avatar and the contributor credit names its own bot address."""
    import delfin.user_settings as us
    monkeypatch.setattr(
        us, "load_settings",
        lambda *a, **k: {"agent": {
            "commit_coauthor_email":
                "12345678+delfin-agent@users.noreply.github.com"}})
    prompt = _solo_prompt("kit.glm-5.3")
    assert ("Co-Authored-By: DELFIN-Agent (kit.glm-5.3) "
            "<12345678+delfin-agent@users.noreply.github.com>") in prompt


def test_without_a_setting_it_belongs_to_nobody(monkeypatch):
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {"agent": {}})
    assert "<noreply@delfin-agent.invalid>" in _solo_prompt("kit.glm-5.3")


def test_the_prompt_does_not_say_both_things():
    """A learned rule said "No Co-Authored-By in commits" while the
    attribution block asked for exactly that trailer — the same prompt,
    two instructions. It meant no THIRD-PARTY attribution."""
    prompt = _solo_prompt("kit.glm-5.3")
    assert "No Co-Authored-By in commits" not in prompt
    assert "Co-Authored-By: DELFIN-Agent (kit.glm-5.3)" in prompt


def test_the_shipped_address_is_the_agents_own_account():
    """Its noreply form: the trailer links to the machine account and no
    personal address enters the history."""
    prompt = _solo_prompt("kit.glm-5.3")
    assert "@users.noreply.github.com>" in prompt
    assert "@gmail.com" not in prompt


def test_a_wave_member_is_not_asked_to_sign_its_commits(monkeypatch):
    """A member's commits are merged as the wave's work, so the line is public
    text on a history nobody can reword: the prompt gives it nothing. Red on
    the base, where the block read only `agent.commit_coauthor` and every
    member was told to sign."""
    import delfin.user_settings as us
    monkeypatch.setenv("DELFIN_SWARM_RUN", "w48")
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {"agent": {}})
    assert "Co-Authored-By: DELFIN-Agent" not in _solo_prompt("kit.glm-5.3")


def test_a_deployment_can_have_its_members_sign_anyway(monkeypatch):
    """Off, not gone: one setting restores the line for a wave."""
    import delfin.user_settings as us
    monkeypatch.setenv("DELFIN_SWARM_RUN", "w48")
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {
        "agent": {"swarm_commit_coauthor": True}})
    assert "Co-Authored-By: DELFIN-Agent (kit.glm-5.3)" in _solo_prompt("kit.glm-5.3")
