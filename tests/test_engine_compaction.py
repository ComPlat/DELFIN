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
    # Compaction is token-driven (the legacy message-count trigger was
    # removed — it fired at ~15% window full). Shrink the window so the
    # long history below crosses the auto-compact threshold; the 12-message
    # floor still protects short conversations from compacting at all.
    engine.context_window_tokens = 10
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


# ---------------------------------------------------------------------------
# Re-compaction fidelity: a prior summary must not be re-truncated to 400 chars
# (that silently dropped any fact sitting deep in it). Exercises the compaction
# methods directly via __new__ so no agent-pack is needed.
# ---------------------------------------------------------------------------


def _bare_engine():
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 10      # force compaction (tiny window)
    eng.auto_compact_pct = 0.95
    eng.backend = "cli"                 # extractive path (no LLM call)
    eng.client = None
    eng.last_compaction_info = {}
    eng.session_id = ""
    return eng


def test_prior_summary_survives_second_compaction():
    """A load-bearing fact sitting DEEP in the first summary must survive when
    a second compaction folds that summary back in. Regression: the prior
    summary was treated as a 400-char user goal, dropping everything past it."""
    eng = _bare_engine()
    eng.messages.append({"role": "user", "content": "ORIGINAL_GOAL: build the EOS solver"})
    for i in range(5):
        eng.messages.append({"role": "user", "content": f"step {i}: " + "work " * 18})
        eng.messages.append({"role": "assistant", "content": f"did {i}: " + "done " * 18})
    eng.messages.append({"role": "user",
                         "content": "DEEP_CONSTRAINT: bug reports group-readable to qmchem_shared only"})
    for i in range(5, 8):
        eng.messages.append({"role": "user", "content": f"step {i}: " + "work " * 18})
        eng.messages.append({"role": "assistant", "content": f"did {i}: " + "done " * 18})

    eng._compact_history()
    block1 = eng.messages[0]["content"]
    assert block1.startswith("[Conversation summary")
    # precondition: the constraint really is deep (past the old 400-char cut)
    assert "DEEP_CONSTRAINT" in block1 and block1.find("DEEP_CONSTRAINT") > 400

    for i in range(12):
        eng.messages.append({"role": "user" if i % 2 == 0 else "assistant",
                             "content": f"later {i}: " + "more " * 10})
    eng._compact_history()
    block2 = eng.messages[0]["content"]
    assert "ORIGINAL_GOAL" in block2      # original intent preserved
    assert "DEEP_CONSTRAINT" in block2    # deep fact no longer truncated away
