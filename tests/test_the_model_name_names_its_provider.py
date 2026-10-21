"""A model id that names its provider starts that provider when none is set.

Input: ``cli._build_engine`` with ``--model`` and no ``--provider``, in a
HOME without settings. Output: the engine is constructed with the provider
``pricing.infer_provider`` reads from the id (``kit.`` / ``azure.`` -> kit,
``name:tag`` -> ollama, ``gpt-`` / ``o1`` / ``o3`` / ``o4`` -> openai,
``claude`` -> claude), and the start notes name the choice. An explicit
``--provider`` and a provider from the settings file are never replaced;
an id that names no provider leaves the provider empty, as before.

Before: the provider stayed empty, ``create_client`` fell through to the
Anthropic backend, and ``--model kit.deepseek-v4-flash`` failed there with
an error about a key it was never meant to use.
"""
from __future__ import annotations

import argparse

import pytest

from delfin.agent import cli as agent_cli


class _Recorder:
    last: dict = {}

    def __init__(self, **kwargs):
        type(self).last = dict(kwargs)


@pytest.fixture
def recorded(monkeypatch, tmp_path):
    import delfin.agent.engine as engine_mod
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    monkeypatch.setattr(engine_mod, "AgentEngine", _Recorder)
    _Recorder.last = {}
    return _Recorder


def _args(**over) -> argparse.Namespace:
    base = dict(backend="", provider="", model="", mode="", cwd="",
                effort="", permission_mode="", settings_defaults=True)
    base.update(over)
    return argparse.Namespace(**base)


@pytest.mark.parametrize("model,provider", [
    ("kit.deepseek-v4-flash", "kit"),
    ("azure.gpt-5.4", "kit"),
    ("qwen3-coder:32b", "ollama"),
    ("gpt-4o", "openai"),
    ("claude-sonnet-4", "claude"),
])
def test_the_provider_follows_the_model_name(recorded, tmp_path, model,
                                             provider):
    args = _args(model=model, cwd=str(tmp_path))
    agent_cli._build_engine(args)
    assert recorded.last.get("provider") == provider
    notes = agent_cli._provider_notes(args)
    assert len(notes) == 1
    assert provider in notes[0] and model in notes[0]


def test_an_explicit_provider_always_wins(recorded, tmp_path):
    args = _args(model="kit.deepseek-v4-flash", provider="openai",
                 cwd=str(tmp_path))
    agent_cli._build_engine(args)
    assert recorded.last.get("provider") == "openai"
    assert agent_cli._provider_notes(args) == []


def test_the_settings_provider_is_not_replaced(recorded, monkeypatch,
                                               tmp_path):
    monkeypatch.setattr(agent_cli, "_engine_defaults",
                        lambda: {"provider": "ollama"})
    args = _args(model="kit.deepseek-v4-flash", cwd=str(tmp_path))
    agent_cli._build_engine(args)
    assert recorded.last.get("provider") == "ollama"
    assert agent_cli._provider_notes(args) == []


def test_a_model_that_names_no_provider_changes_nothing(recorded, tmp_path):
    args = _args(model="some-local-model", cwd=str(tmp_path))
    agent_cli._build_engine(args)
    assert recorded.last.get("provider") == ""
    assert agent_cli._provider_notes(args) == []


def test_run_infers_too(recorded, tmp_path):
    """`run` keeps out of the settings file, but a model id it was given
    still names its provider."""
    args = _args(model="kit.deepseek-v4-flash", settings_defaults=False,
                 cwd=str(tmp_path))
    agent_cli._build_engine(args)
    assert recorded.last.get("provider") == "kit"


def test_the_interactive_start_carries_the_note(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()

    class _Engine:
        def set_kit_permission_mode(self, *a, **k):
            pass

    def _build(ns):
        agent_cli._resolve_provider(ns)
        return _Engine()

    monkeypatch.setattr(agent_cli, "_build_engine", _build)
    monkeypatch.setattr(agent_cli, "_open_session", lambda *a, **k: True)
    monkeypatch.setattr(agent_cli, "_apply_launch_grants", lambda *a: None)
    monkeypatch.setattr(agent_cli, "_apply_isolation_flags", lambda a: "")
    import delfin.agent.isolation_posture as ip
    import delfin.agent.run_limits as rl
    monkeypatch.setattr(ip, "apply_at_start", lambda *a, **k: [])
    monkeypatch.setattr(rl, "apply_unattended", lambda *a, **k: [])
    args = _args(model="kit.deepseek-v4-flash", cwd=str(ws))
    import os
    before = os.getcwd()
    try:
        start = agent_cli._start_interactive_session(args, None, ws,
                                                     can_ask=False)
    finally:
        os.chdir(before)
    assert any("kit" in n and "kit.deepseek-v4-flash" in n
               for n in start.start_notes), start.start_notes
