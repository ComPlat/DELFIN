"""The per-turn steering blocks have to reach the turn they steer.

Every block in the interactive prompt -- the open-task reminder, the run
budget wind-down, a late answer to a parked question -- is built once,
while the system prompt for the turn is assembled, and then frozen until
the turn ends. That is the wrong moment for all of them:

  * The reminder exists so a finished step gets checked off. At turn
    start the agent has not created the task yet, so the reminder is
    empty; by the time the task exists, the prompt that would mention it
    is already sent.
  * The budget wind-down fires at 80% spend. A turn that crosses 80% in
    round 30 reads about it only if the user sends another message.
  * A benchmark task is a single turn, so on this project's own
    workspaces these blocks measured empty in EVERY case.

The fix is the channel the background-job events already use: ask the
engine between tool rounds and inject what changed. These tests hold the
mechanism to what makes it worth its tokens -- it must carry a change,
it must not repeat itself, and it must be bounded.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E


@pytest.fixture
def office_engine(tmp_path):
    """An office-mode engine on a scratch workspace, with no live client."""
    ws = tmp_path / "buero"
    ws.mkdir()
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(
            repo_dir=ws, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="office")
    # Turn start: this is what freezes the blocks for the rest of the turn.
    eng._build_current_system_prompt("", task_text="Wie viel im Juni?")
    return eng


def _open_task(engine, subject: str) -> None:
    from delfin.agent.agent_tasks import get_store
    perms = engine.kit_permissions
    get_store(perms.workspace).create(
        subject, session_id=getattr(perms, "task_session_id", "") or "")


# ---------------------------------------------------------------------------
# It carries what changed
# ---------------------------------------------------------------------------

def test_a_task_created_mid_turn_reaches_the_model(office_engine):
    """The case the reminder exists for, and the one it never covered."""
    assert office_engine._drain_turn_steering() == []

    _open_task(office_engine, "Rechnungen abgleichen")

    blocks = office_engine._drain_turn_steering()
    assert blocks, "a task created during the turn never reached the model"
    assert "Rechnungen abgleichen" in "\n".join(blocks)


def test_the_same_block_is_not_sent_twice(office_engine):
    _open_task(office_engine, "Rechnungen abgleichen")
    first = office_engine._drain_turn_steering()
    assert first
    assert office_engine._drain_turn_steering() == [], (
        "an unchanged block was re-sent, paying tokens to repeat itself")


def test_a_block_already_in_the_system_prompt_is_not_repeated(tmp_path):
    """The turn-start prompt and the refresh must not both carry it."""
    ws = tmp_path / "buero"
    ws.mkdir()
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(
            repo_dir=ws, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="office")
    _open_task(eng, "Schon offen vor dem Turn")
    prompt = eng._build_current_system_prompt("", task_text="Weiter")
    assert "Schon offen vor dem Turn" in prompt
    assert eng._drain_turn_steering() == []


def test_the_refresh_is_capped(office_engine):
    """A flapping list must not be able to spend the turn on itself."""
    for i in range(E._MAX_STEERING_REFRESHES + 4):
        _open_task(office_engine, f"Aufgabe {i}")
        office_engine._drain_turn_steering()
    _open_task(office_engine, "Eine mehr")
    assert office_engine._drain_turn_steering() == []
    assert office_engine._steering_refreshes <= E._MAX_STEERING_REFRESHES


def test_a_new_turn_starts_the_budget_over(office_engine):
    for i in range(E._MAX_STEERING_REFRESHES + 2):
        _open_task(office_engine, f"Aufgabe {i}")
        office_engine._drain_turn_steering()
    office_engine._build_current_system_prompt("", task_text="Nächste Frage")
    assert office_engine._steering_refreshes == 0
    _open_task(office_engine, "Nach dem Turnwechsel")
    assert office_engine._drain_turn_steering()


# ---------------------------------------------------------------------------
# It stays out of where it does not belong
# ---------------------------------------------------------------------------

def test_pipeline_roles_get_no_steering(tmp_path):
    """Scripted roles run one step and keep no list of their own."""
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="quick")
    assert eng.current_role not in E._STEERING_ROLES
    _open_task(eng, "Nicht für diese Rolle")
    assert eng._drain_turn_steering() == []


def test_the_office_role_is_included():
    """Office works on real records for hours; it needs the list most."""
    assert "office_agent" in E._STEERING_ROLES


def test_the_static_blocks_are_not_refreshed():
    """Re-sending them would cost tokens to repeat the system prompt."""
    assert "project_dir" not in E._MID_TURN_STEERING_KEYS
    assert "backend_parity" not in E._MID_TURN_STEERING_KEYS


def test_the_context_snapshot_speaks_only_when_it_warns(office_engine, monkeypatch):
    """A percentage that ticks up every round is noise, not steering."""
    monkeypatch.setattr(E.AgentEngine, "_build_context_status_block",
                        lambda self: "- Context: 41% of window used")
    monkeypatch.setattr(office_engine, "_steering_blocks",
                        lambda role: [("context_status",
                                       office_engine._build_context_status_block())])
    assert office_engine._drain_turn_steering() == []
    monkeypatch.setattr(E.AgentEngine, "_build_context_status_block",
                        lambda self: "- Context: 84% — WARNING: nearing auto-compaction")
    assert office_engine._drain_turn_steering()


def test_a_broken_block_builder_does_not_end_the_turn(office_engine, monkeypatch):
    def boom(role):
        raise RuntimeError("store unreadable")
    monkeypatch.setattr(office_engine, "_steering_blocks", boom)
    assert office_engine._drain_turn_steering() == []


# ---------------------------------------------------------------------------
# The loop side
# ---------------------------------------------------------------------------

def test_the_loop_asks_only_when_the_engine_offered():
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    assert client._drain_turn_steering() == []


def test_the_loop_passes_the_blocks_through():
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client.steering_provider = lambda: ["[live update]\n- one", "", None]
    assert client._drain_turn_steering() == ["[live update]\n- one"]


def test_a_failing_callback_does_not_end_the_turn():
    client = A.OpenAIClient.__new__(A.OpenAIClient)

    def boom():
        raise RuntimeError("engine gone")

    client.steering_provider = boom
    assert client._drain_turn_steering() == []


def test_the_engine_installs_the_callback(office_engine):
    """Without the wiring the loop has nothing to ask."""
    from delfin.agent.api_client import StreamEvent

    def fake_stream(**kwargs):
        yield StreamEvent(type="text_delta", text="ok")

    office_engine.client.stream_message = MagicMock(side_effect=fake_stream)
    office_engine.stream_response("Was steht an?")
    assert callable(getattr(office_engine.client, "steering_provider", None))


def test_the_refresh_can_be_switched_off(office_engine, monkeypatch):
    """So a benchmark arm can measure this mechanism against its absence."""
    _open_task(office_engine, "Rechnungen abgleichen")
    monkeypatch.setenv("DELFIN_TURN_STEERING", "0")
    assert office_engine._drain_turn_steering() == []
    monkeypatch.delenv("DELFIN_TURN_STEERING")
    assert office_engine._drain_turn_steering()
