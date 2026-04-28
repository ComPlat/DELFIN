"""Tests for AgentEngine._compact_history CLI-restart behaviour.

S5 introduced the user-settings opt-out ``agent.compact_resets_cli``.
Default True (kill on compact) keeps the historical behaviour;
False keeps the subprocess alive and only mutates local messages.
"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.engine import AgentEngine


@pytest.fixture
def agent_tree(tmp_path):
    """Minimal pack tree so AgentEngine can boot.

    Mirrors the layout used in tests/test_agent_prompt_loader.py
    (pack/ for roles, pack_lite/ for modes manifest).
    """
    import textwrap

    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()

    (shared / "delfin_context.md").write_text("# DELFIN Context\nTest.")
    (shared / "work_cycle_rules.md").write_text("# Work Cycle Rules\nRule 1.")
    (shared / "universal_input_template.md").write_text("# Input Template")
    (shared / "minimal_final_verdict.md").write_text("# Verdict")
    (agents / "solo_agent.md").write_text("# Solo Agent\nYou are solo.")

    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# Mode: quick\nQuick mode.")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        recommended_default_mode: quick
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - solo_agent
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _build_engine(agent_tree, *, n_messages=20):
    """Boot an engine with a fake CLI client and a long history."""
    fake_client = MagicMock()
    fake_client.session_id = "sess-1"
    with patch("delfin.agent.engine.create_client", return_value=fake_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )
    # current_role is derived from the mode's route — solo_agent mode
    # gives us solo_agent as current_role automatically.
    # Long alternating history (must exceed _COMPACTION_THRESHOLD = 12)
    for i in range(n_messages):
        role = "user" if i % 2 == 0 else "assistant"
        engine.messages.append({"role": role, "content": f"msg {i}"})
    return engine


def _settings_dir_with(setting_value, tmp_path, monkeypatch):
    """Point load_settings at a temp file with the given compact_resets_cli value."""
    cfg_dir = tmp_path / "settings"
    cfg_dir.mkdir(exist_ok=True)
    cfg_path = cfg_dir / ".delfin_settings.json"
    cfg_path.write_text(json.dumps({"agent": {"compact_resets_cli": setting_value}}))
    monkeypatch.setenv("HOME", str(cfg_dir))
    return cfg_path


def test_compact_kills_cli_when_setting_true_default(agent_tree):
    """Default behaviour: kill the CLI subprocess after compaction."""
    engine = _build_engine(agent_tree, n_messages=20)
    engine._compact_history()
    engine.client.kill.assert_called_once()


def test_compact_does_not_kill_cli_when_setting_false(
    agent_tree, tmp_path, monkeypatch,
):
    """compact_resets_cli=False keeps the subprocess alive."""
    _settings_dir_with(False, tmp_path, monkeypatch)
    engine = _build_engine(agent_tree, n_messages=20)
    engine._compact_history()
    engine.client.kill.assert_not_called()


def test_compact_local_messages_always_trimmed(agent_tree, tmp_path, monkeypatch):
    """Whether or not the CLI was killed, local messages are always compacted."""
    _settings_dir_with(False, tmp_path, monkeypatch)
    engine = _build_engine(agent_tree, n_messages=20)
    engine._compact_history()
    # Compaction → 1 summary user + 1 ack assistant + last 4 originals = 6
    assert len(engine.messages) == 6
    assert "Conversation summary" in engine.messages[0]["content"]


def test_compact_below_threshold_is_noop(agent_tree):
    """Short conversations don't compact at all."""
    engine = _build_engine(agent_tree, n_messages=4)
    engine._compact_history()
    engine.client.kill.assert_not_called()
    assert len(engine.messages) == 4
