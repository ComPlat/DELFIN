"""Tests for delfin.agent.engine (with mocked backends)."""

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    """Create a minimal agent directory tree."""
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()

    (shared / "delfin_context.md").write_text("# Context")
    (shared / "work_cycle_rules.md").write_text("# Rules")
    (shared / "universal_input_template.md").write_text("")
    (shared / "minimal_final_verdict.md").write_text("")
    (agents / "session_manager.md").write_text("# Session Manager")
    (agents / "builder_agent.md").write_text("# Builder Agent")
    (agents / "test_agent.md").write_text("# Test Agent")

    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# quick mode")

    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - session_manager
              - builder_agent
              - test_agent
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)

    return tmp_path


@pytest.fixture
def mock_client():
    """Create a mock client that yields StreamEvents."""
    from delfin.agent.api_client import StreamEvent

    def fake_stream(system, messages, max_tokens=4096, session_id="", thinking_budget=0):
        yield StreamEvent(type="session_init", text="test-session-123")
        yield StreamEvent(type="message_start", input_tokens=100)
        yield StreamEvent(type="text_delta", text="Hello ")
        yield StreamEvent(type="text_delta", text="from ")
        yield StreamEvent(type="text_delta", text="Claude!")
        yield StreamEvent(type="message_delta", output_tokens=50, cost_usd=0.01)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=fake_stream)
    return client


def test_engine_init(agent_tree, mock_client):
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree,
            backend="cli",
            mode="quick",
            pack_dir=agent_tree,
        )
    assert engine.mode == "quick"
    assert engine.route == ["session_manager", "builder_agent", "test_agent"]
    assert engine.current_role == "session_manager"
    assert engine.current_role_index == 0
    assert not engine.is_cycle_complete
    assert engine.backend == "cli"


def test_engine_stream_response(agent_tree, mock_client):
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree
        )

    chunks = []
    response = engine.stream_response(
        "Hello agent",
        on_token=lambda t: chunks.append(t),
    )

    assert response == "Hello from Claude!"
    assert chunks == ["Hello ", "from ", "Claude!"]
    assert len(engine.messages) == 2
    assert engine.messages[0]["role"] == "user"
    assert engine.messages[1]["role"] == "assistant"
    assert engine.token_usage["input"] == 100
    assert engine.token_usage["output"] == 50
    assert engine.cost_usd == 0.01
    assert engine.session_id == "test-session-123"


def test_engine_advance_role(agent_tree, mock_client):
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree
        )

    engine.stream_response("Plan the work")

    assert engine.current_role == "session_manager"
    has_next = engine.advance_role()
    assert has_next is True
    assert engine.current_role == "builder_agent"
    assert "session_manager" in engine.role_outputs

    engine.stream_response("Build it")
    has_next = engine.advance_role()
    assert has_next is True
    assert engine.current_role == "test_agent"

    engine.stream_response("Test it")
    has_next = engine.advance_role()
    assert has_next is False
    assert engine.is_cycle_complete


def test_engine_reset_cycle(agent_tree, mock_client):
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree
        )

    engine.stream_response("Hello")
    engine.advance_role()

    engine.reset_cycle()
    assert engine.current_role_index == 0
    assert engine.messages == []
    assert engine.role_outputs == {}
    assert engine.token_usage == {"input": 0, "output": 0}
    assert engine.cost_usd == 0.0


def test_engine_get_status(agent_tree, mock_client):
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree
        )

    status = engine.get_status()
    assert status["mode"] == "quick"
    assert status["backend"] == "cli"
    assert status["role"] == "session_manager"
    assert status["role_index"] == 0
    assert status["role_total"] == 3
    assert status["cycle_complete"] is False


def test_engine_request_stop(agent_tree, mock_client):
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    def stream_that_stops_mid(system, messages, max_tokens=4096, session_id="", thinking_budget=0):
        for i in range(10):
            if i == 2:
                engine._stop_requested = True
            yield StreamEvent(type="text_delta", text=f"chunk{i} ")
        yield StreamEvent(type="message_delta", output_tokens=10)

    mock_client.stream_message = MagicMock(side_effect=stream_that_stops_mid)

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree
        )

    response = engine.stream_response("Hello")
    assert "chunk0" in response
    assert "chunk1" in response
    assert "chunk9" not in response


def test_engine_available_modes(agent_tree, mock_client):
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree
        )

    modes = engine.available_modes()
    assert "quick" in modes


def test_cli_client_init():
    """Test CLIClient can be created when claude is in PATH."""
    from delfin.agent.api_client import CLIClient

    with patch("shutil.which", return_value="/usr/bin/claude"):
        client = CLIClient(model="sonnet")
    assert client.model == "sonnet"


def test_cli_client_not_found():
    """Test CLIClient raises when claude is not in PATH."""
    from delfin.agent.api_client import CLIClient

    with patch("shutil.which", return_value=None):
        with pytest.raises(FileNotFoundError, match="Claude Code CLI not found"):
            CLIClient()


def test_create_client_cli():
    """Test factory creates CLIClient for backend='cli'."""
    from delfin.agent.api_client import CLIClient, create_client

    with patch("shutil.which", return_value="/usr/bin/claude"):
        client = create_client(backend="cli", model="sonnet")
    assert isinstance(client, CLIClient)


def test_engine_tool_use_callback(agent_tree):
    """Test that tool_use events are forwarded to the on_tool_use callback."""
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    def stream_with_tools(system, messages, max_tokens=4096, session_id="", thinking_budget=0):
        yield StreamEvent(type="tool_use", tool_name="Read", tool_input='{"file_path": "test.py"}')
        yield StreamEvent(type="text_delta", text="File contents here.")
        yield StreamEvent(type="message_delta", output_tokens=20, cost_usd=0.005)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=stream_with_tools)

    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree)

    tools_used = []
    response = engine.stream_response(
        "Read test.py",
        on_tool_use=lambda name, inp: tools_used.append((name, inp)),
    )
    assert response == "File contents here."
    assert len(tools_used) == 1
    assert tools_used[0][0] == "Read"


def test_engine_session_persistence(agent_tree):
    """Test that session_id is preserved across calls."""
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    call_count = {"n": 0}

    def stream_fn(system, messages, max_tokens=4096, session_id="", thinking_budget=0):
        call_count["n"] += 1
        if call_count["n"] == 1:
            yield StreamEvent(type="session_init", text="session-abc-123")
        yield StreamEvent(type="text_delta", text=f"Response {call_count['n']}")
        yield StreamEvent(type="message_delta", output_tokens=10)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=stream_fn)

    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree)

    engine.stream_response("First message")
    assert engine.session_id == "session-abc-123"

    # Second call should pass the session_id
    engine.stream_response("Second message")
    second_call = client.stream_message.call_args_list[1]
    assert second_call.kwargs.get("session_id") == "session-abc-123" or \
        (len(second_call.args) > 3 and second_call.args[3] == "session-abc-123")

    # Reset should clear session_id
    engine.reset_cycle()
    assert engine.session_id == ""


def test_engine_export_state(agent_tree, mock_client):
    """Test that export_state captures all relevant engine state."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree)

    engine.stream_response("Hello")
    engine.advance_role()

    exported = engine.export_state()
    assert exported["mode"] == "quick"
    assert exported["role_index"] == 1
    assert exported["route"] == ["session_manager", "builder_agent", "test_agent"]
    assert "session_manager" in exported["role_outputs"]
    assert len(exported["engine_messages"]) == 2
    assert exported["session_id"] == "test-session-123"
    assert exported["cost_usd"] == 0.01


def test_engine_restore_state(agent_tree, mock_client):
    """Test that restore_state rebuilds engine from saved data."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree)

    saved = {
        "mode": "quick",
        "role_index": 2,
        "role_outputs": {"session_manager": "plan", "builder_agent": "code"},
        "engine_messages": [
            {"role": "user", "content": "Hi"},
            {"role": "assistant", "content": "Hello"},
        ],
        "token_usage": {"input": 5000, "output": 1200},
        "cost_usd": 0.15,
        "session_id": "restored-session-456",
    }
    engine.restore_state(saved)

    assert engine.current_role == "test_agent"
    assert engine.current_role_index == 2
    assert engine.session_id == "restored-session-456"
    assert engine.cost_usd == 0.15
    assert engine.token_usage["input"] == 5000
    assert len(engine.messages) == 2
    assert "session_manager" in engine.role_outputs

    # Verify streaming still works after restore
    response = engine.stream_response("Continue")
    assert response == "Hello from Claude!"


def test_retry_from_builder(agent_tree, mock_client):
    """Test that retry_from_builder rewinds to the builder step."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )

    # Walk through: session_manager → builder → test
    engine.stream_response("Plan")
    engine.advance_role()  # → builder
    engine.stream_response("Build")
    engine.advance_role()  # → test
    engine.stream_response("Test FAILED: 2 errors")
    assert engine.current_role == "test_agent"

    # Retry should rewind to builder
    assert engine.retry_from_builder() is True
    assert engine.current_role == "builder_agent"
    assert "test_agent" in engine.role_outputs


def test_suggest_mode_detects_cluster():
    """Test that cluster files trigger cluster mode suggestion."""
    from delfin.agent.engine import AgentEngine

    result = AgentEngine.suggest_mode(
        "Fix the bug in delfin/dashboard/backend_slurm.py", "quick"
    )
    assert result == "cluster"


def test_suggest_mode_detects_reviewed():
    """Test that cli.py triggers reviewed mode suggestion."""
    from delfin.agent.engine import AgentEngine

    result = AgentEngine.suggest_mode(
        "Refactor delfin/cli.py argument parsing", "quick"
    )
    assert result == "reviewed"


def test_suggest_mode_no_escalation_needed():
    """Test that unrelated files don't trigger escalation."""
    from delfin.agent.engine import AgentEngine

    result = AgentEngine.suggest_mode(
        "Fix a typo in the README", "quick"
    )
    assert result is None


def test_suggest_mode_already_high_enough():
    """Test that no suggestion if current mode is already sufficient."""
    from delfin.agent.engine import AgentEngine

    result = AgentEngine.suggest_mode(
        "Fix delfin/cli.py", "cluster"
    )
    assert result is None  # cluster > reviewed, no escalation


def test_suggest_mode_cluster_over_reviewed():
    """Test that cluster wins when both cluster and reviewed files mentioned."""
    from delfin.agent.engine import AgentEngine

    result = AgentEngine.suggest_mode(
        "Change delfin/cli.py and backend_slurm.py together", "quick"
    )
    assert result == "cluster"


def test_thinking_budget_for_role():
    """Test that different roles get different thinking budgets."""
    from delfin.agent.engine import AgentEngine

    builder = AgentEngine.thinking_budget_for_role("builder_agent")
    session = AgentEngine.thinking_budget_for_role("session_manager")
    critic = AgentEngine.thinking_budget_for_role("critic_agent")

    assert builder > session  # builder needs more thinking
    assert builder > critic  # builder needs the most
    assert session > 0
    assert critic > 0


def test_build_handoff_message(agent_tree, mock_client):
    """Test that handoff messages include prior outputs and task."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )

    # Session Manager phase
    engine.stream_response("Fix the bug in cli.py")
    engine.advance_role()  # → builder_agent

    msg = engine.build_handoff_message("Fix the bug in cli.py")
    assert "Fix the bug in cli.py" in msg
    assert "session_manager" in msg  # prior output referenced
    assert "Implement" in msg  # builder-specific instruction


def test_build_handoff_message_test_agent(agent_tree, mock_client):
    """Test that test agent handoff includes pytest instruction."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )

    engine.stream_response("Plan")
    engine.advance_role()
    engine.stream_response("Build")
    engine.advance_role()  # → test_agent

    msg = engine.build_handoff_message("Fix the bug")
    assert "pytest" in msg
    assert "criterion" in msg


def test_extract_acceptance_criteria():
    """Test parsing acceptance criteria from Session Manager output."""
    from delfin.agent.engine import AgentEngine

    sm_output = """## PLAN

**Task:** Fix the bug
**Class:** bugfix

### Acceptance criteria
1. The function returns correct values for edge cases
2. All existing tests pass
3. No regressions in CLI behavior

### Execution plan
1. Read the file
"""
    criteria = AgentEngine.extract_acceptance_criteria(sm_output)
    assert len(criteria) == 3
    assert "edge cases" in criteria[0]
    assert "existing tests pass" in criteria[1]
    assert "regressions" in criteria[2]


def test_extract_acceptance_criteria_empty():
    """Test extraction from output without criteria section."""
    from delfin.agent.engine import AgentEngine

    criteria = AgentEngine.extract_acceptance_criteria("Just some text")
    assert criteria == []


def test_pipeline_status(agent_tree, mock_client):
    """Test pipeline_status returns correct step states."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )

    status = engine.pipeline_status()
    assert len(status) == 3
    assert status[0]["status"] == "active"
    assert status[1]["status"] == "pending"
    assert status[2]["status"] == "pending"

    engine.stream_response("Plan")
    engine.advance_role()

    status = engine.pipeline_status()
    assert status[0]["status"] == "done"
    assert status[1]["status"] == "active"
    assert status[2]["status"] == "pending"


def test_compact_for_next_role(agent_tree, mock_client):
    """Test that compaction clears messages but preserves role outputs."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )

    engine.stream_response("Plan")
    engine.advance_role()
    assert len(engine.messages) == 2
    assert "session_manager" in engine.role_outputs

    engine.compact_for_next_role()
    assert engine.messages == []
    assert "session_manager" in engine.role_outputs  # preserved


def test_create_client_api():
    """Test factory creates APIClient for backend='api'."""
    pytest.importorskip("anthropic")
    from delfin.agent.api_client import APIClient, create_client

    client = create_client(backend="api", api_key="test-key", model="claude-sonnet-4-20250514")
    assert isinstance(client, APIClient)
