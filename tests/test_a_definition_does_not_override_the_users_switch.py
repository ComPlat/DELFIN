"""A subagent definition's `model: cheap` is a default; the user's switch
agent.subagents.cheap_tier=false still wins. Only a call may ask past it.

And every correction prompt the harness sends opens with the one phrase the
model has learned to read as "not the user": "Automatic check, not a
message from the user". The auto-verify prompt said "[Harness check -- not
a message from the user]" (review 2026-09-16).
"""
import inspect

from delfin.agent import api_client as A
from delfin.agent import subagents as SA


class _Client:
    model = "big-model"
    _provider = "kit"
    _base_url = ""


def test_a_preset_cheap_tier_respects_the_users_off_switch(monkeypatch):
    monkeypatch.setattr(SA, "_cheap_tier_enabled", lambda: False)
    monkeypatch.setitem(SA.SUBAGENT_PRESETS, "deep-review", SA.SubagentPreset(
        name="deep-review", description="d", system_prompt="s", mode="plan",
        model="cheap"))
    assert SA._resolve_subagent_model(_Client(), "deep-review") == ("big-model", "parent")


def test_a_call_may_still_ask_for_the_cheap_tier(monkeypatch):
    monkeypatch.setattr(SA, "_cheap_tier_enabled", lambda: False)
    monkeypatch.setattr("delfin.agent.model_routing.tier_model", lambda p, t, s=None: "small-model")
    monkeypatch.setattr("delfin.agent.model_routing.is_known_broken", lambda m: False)
    import types
    caps = types.SimpleNamespace(supports_tools=True)
    monkeypatch.setattr("delfin.agent.model_capabilities.resolve", lambda *a, **k: caps)
    assert SA._resolve_subagent_model(_Client(), "explore", "cheap") == ("small-model", "cheap")


def test_the_auto_verify_prompt_uses_the_agreed_phrase():
    src = inspect.getsource(A.OpenAIClient.stream_message)
    assert "[Verify] Automatic check, not a message from the " in src
    assert "Harness check" not in src
