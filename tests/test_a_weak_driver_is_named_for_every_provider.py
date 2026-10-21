"""A model too small to drive an agent is named, whatever serves it.

Session 78c36f8c (2026-10-19) ran an hour of coding on ``azure.gpt-5-nano``
served through the KIT endpoint. The only check was KIT's, keyed on the
``kit.`` prefix, so nothing was said: no tool-format warning, no hint at a
better model -- and the terminal never ran a check at all.
"""

from __future__ import annotations

import pytest

from delfin.agent import model_capabilities as mc


def _caps(model, window=200_000, tools=True):
    return mc.ModelCapabilities(model=model, provider="kit",
                                context_window=window, supports_tools=tools)


@pytest.mark.parametrize("model", [
    "azure.gpt-5-nano", "gpt-5-nano", "gemini-2.5-flash-lite",
    "some-vendor/tiny-chat", "qwen2.5-coder:7b",
])
def test_a_small_tier_is_named(model):
    warn = mc.agent_fitness("openai", model, _caps(model))
    assert model in warn and "Better" in warn


@pytest.mark.parametrize("model", [
    "gpt-5-mini", "kit.deepseek-v4-flash", "claude-opus-5-5",
    "kit.glm-5.3", "qwen3-coder-480b", "llama-3.1-nanotron-70b",
])
def test_a_capable_model_is_not(model):
    assert mc.agent_fitness("openai", model, _caps(model)) == ""


def test_no_tools_or_a_small_window_is_named():
    assert "no tool support" in mc.agent_fitness(
        "ollama", "x", _caps("x", tools=False))
    assert "small context window" in mc.agent_fitness(
        "ollama", "x", _caps("x", window=8_000))


def test_the_kit_endpoint_points_at_its_best_model():
    assert mc.KIT_BEST_MODEL in mc.agent_fitness(
        "kit", "azure.gpt-5-nano", _caps("azure.gpt-5-nano"))


def test_preflight_warns_for_a_model_kits_own_check_passed(monkeypatch):
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: _caps("azure.gpt-5-nano"))
    ok, msg = mc.preflight("kit", "azure.gpt-5-nano")
    assert ok and "smallest tier" in msg


def test_the_terminal_says_it_on_switch(monkeypatch):
    from types import SimpleNamespace
    from delfin.agent import repl_commands as RC
    client = SimpleNamespace(_provider="kit", _base_url="")
    monkeypatch.setattr(mc, "resolve", lambda p, m, *a, **k: _caps(m))
    assert "smallest tier" in RC._fitness_note(client, "azure.gpt-5-nano")
    assert RC._fitness_note(client, "kit.glm-5.3") == ""


def _banner(tmp_path, model):
    from delfin.agent import cli as c

    class _Report:
        git = type("G", (), {"is_repo": False, "branch": "", "dirty": ()})()

    class _Engine:
        client = type("C", (), {"model": model})()
        provider = "kit"
        mode = "solo"
        session_id = ""
        kit_permissions = type("P", (), {"mode": "default"})()

    return c._startup_banner(_Engine(), _Report(), tmp_path)


def test_the_terminal_says_it_at_start_without_asking_the_endpoint(
        monkeypatch, tmp_path):
    def _no_network(*a, **k):
        raise AssertionError("the banner must not resolve over the network")

    monkeypatch.setattr(mc, "resolve", _no_network)
    assert "smallest tier" in _banner(tmp_path, "azure.gpt-5-nano")
    assert "weak driver" not in _banner(tmp_path, "kit.deepseek-v4-flash")
