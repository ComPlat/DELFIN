"""A failed mid-turn checkpoint write is said once, not swallowed.

The engine writes a cheap mid-turn crash checkpoint (throttled) during
long tool loops so a SIGKILL costs the last rounds, not the whole turn.
The write was wrapped in ``except Exception: pass`` — best-effort stays
right (a failed checkpoint must never break the turn), but silent was
wrong: a run whose disk was full, or whose sessions dir was unwritable,
lost the crash insurance for the WHOLE turn while looking perfectly
guarded. Nothing anywhere said so.

Contract pinned here: a failing checkpoint write produces one stream
notice per turn (not one per throttled write, so a full disk cannot
spam), a successful write stays silent, and a failed write never breaks
the turn.
"""

from __future__ import annotations

import pytest
from unittest.mock import MagicMock, patch

from delfin.agent import session_store as ss
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def sessions_dir(monkeypatch, tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    return d


@pytest.fixture
def agent_tree(tmp_path):
    """Minimal agent pack tree (same shape as test_agent_engine.py)."""
    import textwrap
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()
    (shared / "delfin_context.md").write_text("# Context")
    (shared / "work_cycle_rules.md").write_text("# Rules")
    (shared / "goal_decomposition_rules.md").write_text("# Goal Decomposition")
    (shared / "universal_input_template.md").write_text("")
    (shared / "minimal_final_verdict.md").write_text("")
    (agents / "session_manager.md").write_text("# Session Manager")
    (agents / "builder_agent.md").write_text("# Builder Agent")
    (agents / "test_agent.md").write_text("# Test Agent")
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
              - builder_agent
              - test_agent
    """))
    return tmp_path


def _run_turn(notices: list[str], sessions_dir, agent_tree, failing: bool):
    """Stream one 11-round turn, recording notices. Returns the engine's
    returned answer text (stream_response returns the full text)."""
    from delfin.agent.engine import AgentEngine

    def fake_stream(system, messages, max_tokens=4096, session_id="",
                    thinking_budget=0):
        yield StreamEvent(type="session_init", text="ckpt-sid-1")
        yield StreamEvent(type="message_start", input_tokens=100)
        yield StreamEvent(type="text_delta", text="working... ")
        for i in range(11):
            yield StreamEvent(type="tool_use", tool_name="Bash",
                              tool_input='{"command": "ls"}')
            yield StreamEvent(type="tool_result", tool_name="Bash",
                              tool_output=f"round {i}")
        yield StreamEvent(type="text_delta", text="done")
        yield StreamEvent(type="message_delta", output_tokens=10,
                          cost_usd=0.01)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=fake_stream)

    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)

    def on_notice(text):
        notices.append(text)

    real_save = ss.save_turn_checkpoint

    def maybe_failing_save(sid, payload):
        if failing:
            raise OSError("simulated checkpoint write failure")
        return real_save(sid, payload)

    with patch.object(ss, "save_turn_checkpoint", maybe_failing_save):
        answer = engine.stream_response("do a long refactor",
                                        on_notice=on_notice)
    return engine, notices, answer


def test_failing_checkpoint_write_says_once_per_turn(sessions_dir,
                                                     agent_tree):
    """11 tool rounds → one throttled checkpoint write that fails: the
    notice is emitted exactly once, and the turn still completes."""
    notices: list[str] = []
    _, _, answer = _run_turn(notices, sessions_dir, agent_tree, failing=True)
    assert len([n for n in notices if "checkpoint" in n]) == 1
    assert "done" in answer


def test_healthy_checkpoint_write_stays_silent(sessions_dir, agent_tree):
    """A successful write gains no machinery speech."""
    notices: list[str] = []
    _run_turn(notices, sessions_dir, agent_tree, failing=False)
    assert not [n for n in notices if "checkpoint" in n]


def test_failing_write_does_not_break_the_turn(sessions_dir, agent_tree):
    """The failure must not raise or truncate the answer."""
    notices: list[str] = []
    _, _, answer = _run_turn(notices, sessions_dir, agent_tree, failing=True)
    # The turn produced its full answer despite the failing write.
    assert "done" in answer
