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
    (shared / "goal_decomposition_rules.md").write_text("# Goal Decomposition")
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


def test_extract_stage_gates():
    """Test that stage gates are parsed from the Session Manager plan."""
    from delfin.agent.engine import AgentEngine

    plan = textwrap.dedent("""\
        ## PLAN
        **Task:** Improve agent handoffs
        **Class:** agent architecture
        **Risk:** medium
        **Mode:** reviewed

        ### Acceptance criteria
        1. Builder gets the locked goal

        ### Stage gates
        1. Reproduce the orchestration weakness — Exit evidence: failing handoff identified
        2. Patch engine handoffs — Exit evidence: builder sees stage gates
        3. Add tests — Exit evidence: focused suite passes
    """)

    gates = AgentEngine.extract_stage_gates(plan)
    assert len(gates) == 3
    assert "Reproduce the orchestration weakness" in gates[0]
    assert "Patch engine handoffs" in gates[1]


def test_validate_session_manager_output_requires_goal_lock_and_stage_gates():
    """Test that incomplete Session Manager plans are rejected."""
    from delfin.agent.engine import AgentEngine

    incomplete = textwrap.dedent("""\
        ## PLAN
        **Task:** Improve agent reliability
        **Class:** agent architecture
        **Risk:** medium
        **Mode:** quick

        ### Affected files
        - `delfin/agent/engine.py` — update handoff logic

        ### Scope
        - Improve handoffs

        ### Acceptance criteria
        1. Builder sees the plan

        ### Execution plan
        1. Patch engine

        ### Known risks
        - Prompt drift

        **confidence:** medium
        **reason:** Missing fields on purpose
    """)

    errors = AgentEngine.validate_role_output("session_manager", incomplete)
    assert "missing ### Goal lock" in errors
    assert "missing ### Stage gates" in errors
    assert "goal lock must include primary goal, success metric/oracle, and wrong proxy" in errors


def test_validate_test_agent_output_requires_stage_gate_verification():
    """Test that incomplete Test Agent reports are rejected."""
    from delfin.agent.engine import AgentEngine

    incomplete = textwrap.dedent("""\
        ## TEST REPORT

        **Test command:** `python -m pytest tests/ -q`
        **Result:** 10 passed, 0 failed

        **Acceptance criteria verification:**
        1. Builder sees plan — PASS — verified manually

        **Regression check:**
        - None found

        **confidence:** high
        **reason:** Looked good
        **status:** approve
        **key findings:** none
        **open risks:** none
        **recommended next step:** done
    """)

    errors = AgentEngine.validate_role_output("test_agent", incomplete)
    assert "missing **Stage gate verification:**" in errors
    assert "missing **New tests written:**" in errors


def test_validate_reviewer_output_accepts_skip_and_question():
    """Test that reviewer special control outputs are allowed."""
    from delfin.agent.engine import AgentEngine

    assert AgentEngine.validate_role_output("reviewer_agent", "SKIP — trivial change.") == []
    assert AgentEngine.validate_role_output("reviewer_agent", "QUESTION: should I block on this?") == []


def test_evaluate_role_gate_pauses_on_incomplete_session_manager_plan():
    """Test that incomplete Session Manager plans trigger a schema gate."""
    from delfin.agent.engine import AgentEngine

    incomplete_plan = textwrap.dedent("""\
        ## PLAN
        **Task:** Fix something
        **Class:** bugfix
        **Risk:** low
        **Mode:** quick

        ### Affected files
        - `foo.py`

        ### Scope
        - Fix the thing

        ### Acceptance criteria
        1. It works

        ### Execution plan
        1. Do it

        ### Known risks
        - None

        **confidence:** high
        **reason:** simple
    """)

    action, gate_type, message = AgentEngine.evaluate_role_gate("session_manager", incomplete_plan)
    assert action == "pause"
    assert gate_type == "schema"
    assert "incomplete" in message.lower()
    assert "Goal lock" in message or "Stage gates" in message


def test_evaluate_role_gate_skips_conversational_session_manager():
    """Test that conversational SM responses don't trigger a schema gate."""
    from delfin.agent.engine import AgentEngine

    greeting = "Hallo! Welcome to DELFIN. What would you like to work on?"
    action, gate_type, message = AgentEngine.evaluate_role_gate("session_manager", greeting)
    assert action == "continue"
    assert gate_type == ""


def test_is_conversational_detects_greetings():
    """Test that greetings are detected as conversational."""
    from delfin.agent.engine import AgentEngine

    assert AgentEngine.is_conversational("session_manager", "Hello! How can I help?") is True
    assert AgentEngine.is_conversational("session_manager", "## PLAN\n**Task:** Fix bug") is False
    assert AgentEngine.is_conversational("session_manager", "") is False


def test_validate_role_output_accepts_conversational():
    """Test that conversational SM output passes validation."""
    from delfin.agent.engine import AgentEngine

    errors = AgentEngine.validate_role_output("session_manager", "Hi there! What's the task?")
    assert errors == []


def test_evaluate_role_gate_continues_on_complete_session_manager_plan():
    """Test that a complete Session Manager plan passes the gate."""
    from delfin.agent.engine import AgentEngine

    complete_plan = textwrap.dedent("""\
        ## PLAN
        **Task:** Fix something
        **Class:** bugfix
        **Risk:** low
        **Mode:** quick

        ### Affected files
        - `foo.py` — fix the thing

        ### Goal lock
        - Primary goal: fix the actual bug
        - Success metric / oracle: test passes with correct output
        - Wrong proxy to avoid: silencing the error without fixing root cause

        ### Scope
        - Fix the thing

        ### Out of scope
        - Refactoring

        ### Acceptance criteria
        1. It works correctly

        ### Stage gates
        1. Reproduce the bug — Exit evidence: failing test
        2. Fix the root cause — Exit evidence: test passes

        ### Execution plan
        1. Read the file
        2. Fix the bug

        ### Known risks
        - None

        **confidence:** high
        **reason:** simple fix
    """)

    action, gate_type, message = AgentEngine.evaluate_role_gate("session_manager", complete_plan)
    assert action == "continue"
    assert gate_type == ""


def test_evaluate_role_gate_continues_on_approve_with_risks():
    """approve_with_risks should auto-continue — pausing wastes tokens."""
    from delfin.agent.engine import AgentEngine

    report = textwrap.dedent("""\
        ## RESEARCH REPORT
        **confidence:** high
        **reason:** enough sources
        **status:** approve_with_risks
        **key findings:** proxy is weak
        **recommended next step:** tighten the oracle
    """)

    action, gate_type, message = AgentEngine.evaluate_role_gate("research_agent", report)
    assert action == "continue"

    # But reject should still pause
    reject_report = report.replace("approve_with_risks", "reject")
    action2, gate_type2, _ = AgentEngine.evaluate_role_gate("research_agent", reject_report)
    assert action2 == "pause"
    assert gate_type2 == "risk"


def test_evaluate_role_gate_pauses_on_builder_partial_or_blocked():
    """Test that Builder cannot silently pass partial or blocked work onward."""
    from delfin.agent.engine import AgentEngine

    report = textwrap.dedent("""\
        ## BUILD REPORT

        **Changes made:**
        1. `delfin/agent/engine.py` — added something

        **Stage gate status:**
        1. Lock goal — DONE — extracted from plan
        2. Carry gates into handoff — PARTIAL — builder sees them but reviewer does not

        **Critic/Reviewer/Runtime findings addressed:**
        - none

        **Tests run:**
        - `pytest -q` — passed

        **Acceptance criteria status:**
        1. Builder sees goal lock — DONE
        2. Reviewer sees goal lock — BLOCKED

        **Remaining work:**
        - fix reviewer handoff

        **confidence:** medium
        **reason:** still incomplete
        **status:** approve_with_risks
        **open risks:** reviewer still blind
        **recommended next step:** verify reviewer handoff
    """)

    action, gate_type, message = AgentEngine.evaluate_role_gate("builder_agent", report)
    assert action == "pause"
    assert gate_type == "partial"
    assert ("blocked items" in message.lower()) or ("partial items" in message.lower())


def test_evaluate_role_gate_pauses_on_reviewer_goal_lock_issues():
    """Test that reviewer goal-lock failures are treated as a gate."""
    from delfin.agent.engine import AgentEngine

    report = textwrap.dedent("""\
        ## CODE REVIEW

        **Files reviewed:**
        1. `delfin/agent/engine.py` — handoff changes

        **Findings:**
        1. [MINOR] `delfin/agent/engine.py:1` — none — no fix

        **Goal-lock check:**
        - ISSUES — builder solved a weaker proxy

        **Verdict:** PASS
        **Summary:** Logic compiles but the wrong goal was optimized.

        **confidence:** high
        **reason:** direct diff review
    """)

    action, gate_type, message = AgentEngine.evaluate_role_gate("reviewer_agent", report)
    assert action == "pause"
    assert gate_type == "goal-lock"
    assert "goal-lock" in message.lower()


def test_builder_handoff_includes_locked_contract(agent_tree, mock_client):
    """Test that builder handoff includes the locked plan contract and stage gates."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree)

    engine.role_outputs["session_manager"] = textwrap.dedent("""\
        ## PLAN
        **Task:** Fix wrong success metric in agent cycle
        **Class:** agent architecture
        **Risk:** medium
        **Mode:** quick

        ### Goal lock
        - Primary goal: preserve the user's actual target
        - Success metric / oracle: downstream handoffs carry goal, scope, and gates
        - Wrong proxy to avoid: "pipeline completed" without solving the right problem

        ### Scope
        - Update engine handoffs

        ### Out of scope
        - Rewriting the dashboard UI

        ### Acceptance criteria
        1. Builder sees locked goal details

        ### Stage gates
        1. Parse goal lock from Session Manager plan — Exit evidence: handoff contains goal lock
        2. Pass stage gates to Builder — Exit evidence: handoff contains stage gates
    """)
    engine.current_role_index = 1  # builder_agent

    handoff = engine.build_handoff_message("Improve the agent architecture")
    assert "Locked plan contract:" in handoff
    assert "Primary goal: preserve the user's actual target" in handoff
    assert "Stage gates:" in handoff
    assert "Do not silently substitute an easier proxy metric" in handoff


def test_reviewer_handoff_is_role_specific(agent_tree, mock_client):
    """Test that reviewer_agent receives a dedicated review handoff."""
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree)

    engine.route = ["session_manager", "builder_agent", "reviewer_agent", "test_agent"]
    engine.current_role_index = 2
    engine.role_outputs["session_manager"] = textwrap.dedent("""\
        ## PLAN
        **Task:** Review real changes against locked goal
        **Class:** agent architecture
        **Risk:** medium
        **Mode:** reviewed

        ### Goal lock
        - Primary goal: catch goal drift in review

        ### Acceptance criteria
        1. Reviewer checks for goal drift

        ### Stage gates
        1. Add reviewer handoff — Exit evidence: reviewer sees locked contract
    """)
    engine.role_outputs["builder_agent"] = "## BUILD REPORT\nImplemented the new handoff."

    handoff = engine.build_handoff_message("Improve reviewed mode")
    assert "Review the implemented changes for this task" in handoff
    assert "Check correctness, regressions, and whether the Builder stayed aligned" in handoff
    assert "Primary goal: catch goal drift in review" in handoff


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


def test_recommend_task_route_prefers_dashboard_for_dashboard_ops():
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "Please set CONTROL key METHOD to r2scan-3c and submit the job in the dashboard.",
        "dashboard",
    )
    assert decision["mode"] == "dashboard"
    assert decision["task_class"] == "dashboard"
    assert decision["intent"] == "operate"


def test_recommend_task_route_prefers_research_for_chemistry_questions():
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "Which DFT functional and basis set should I use for a Ni metal complex redox question?",
        "dashboard",
    )
    assert decision["mode"] == "research"
    assert decision["task_class"] == "chemistry"
    assert decision["intent"] == "research"


def test_recommend_task_route_prefers_solo_for_code_questions():
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "How does delfin/pipeline.py handle retry logic and where is it wired into the codebase?",
        "dashboard",
    )
    assert decision["mode"] == "solo"
    assert decision["task_class"] == "coding"
    assert decision["intent"] == "question"


def test_recommend_task_route_escalates_cluster_for_runtime_changes():
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "Fix restart handling in delfin/dashboard/backend_slurm.py and scratch recovery logic.",
        "quick",
    )
    assert decision["mode"] == "cluster"
    assert decision["risk_flags"]["cluster"] is True


def test_recommend_task_route_escalates_reviewed_for_api_semantics():
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "Change CONTROL validation and public API semantics for result parsing.",
        "quick",
    )
    assert decision["mode"] == "reviewed"
    assert decision["intent"] == "change"


def test_recommend_task_route_chemistry_code_change_goes_reviewed():
    """Chemistry + code change verbs should route to reviewed (not quick)."""
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "Fix the CREST conformer search implementation for metal complexes.",
        "quick",
    )
    assert decision["mode"] == "reviewed"
    assert decision["task_class"] == "chemistry"
    assert decision["intent"] == "change"


def test_recommend_task_route_crest_keyword_recognized():
    """Verify that new chemistry keywords like CREST are detected."""
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "What is the best CREST conformer search protocol for transition metals?",
        "dashboard",
    )
    assert decision["task_class"] == "chemistry"
    assert decision["mode"] == "research"


def test_recommend_task_route_occupier_code_change():
    """Changes to OCCUPIER workflow files should escalate to reviewed."""
    from delfin.agent.engine import AgentEngine

    decision = AgentEngine.recommend_task_route(
        "Refactor the OCCUPIER auto tree logic in delfin/occupier_auto.py.",
        "quick",
    )
    assert decision["mode"] == "reviewed"


def test_suggest_mode_escalates_for_chemistry_files():
    """Chemistry workflow files should trigger reviewed mode suggestion."""
    from delfin.agent.engine import AgentEngine

    assert AgentEngine.suggest_mode("Fix delfin/esd_module.py", "quick") == "reviewed"
    assert AgentEngine.suggest_mode("Fix delfin/xtb_crest.py", "quick") == "reviewed"
    assert AgentEngine.suggest_mode("Fix delfin/calculators.py", "quick") == "reviewed"


def test_research_agent_uses_sonnet():
    """Research agent should use sonnet for chemistry method synthesis."""
    from delfin.agent.engine import AgentEngine

    assert AgentEngine.model_for_role("research_agent") == "sonnet"


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
