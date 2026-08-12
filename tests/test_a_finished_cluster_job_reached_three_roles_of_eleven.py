"""A finished background job was announced to three roles out of eleven.

The per-turn blocks were all gated on ``_STEERING_ROLES`` -- solo,
dashboard, office. That gate earns its place for the blocks recomputed
from state every turn: a scripted pipeline role keeps no task list of its
own, has no run budget to wind down and is pinned to its folder by
construction, so those blocks would cost it tokens to say nothing.

It buys nothing for the drain-backed blocks. Those are empty unless
something actually finished, they are consumed exactly once by their own
store, and a builder or test role that submitted a cluster job needs the
completion exactly as much as a solo one does -- more, since it has
fewer rounds to spend polling for it.

So the gate is dropped for the drain-backed keys only: finished jobs,
finished background delegates, and a late answer to a question that had
timed out. Everything else stays gated.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E


def _real_perms(engine, workspace) -> None:
    """Replace the mock client's mock permissions with real ones.

    ``engine.kit_permissions`` reads ``client._permissions``; on a mock
    client that is a mock, and a task store built from a mock workspace
    resolves relative to the process CWD -- inside the checkout.
    """
    engine.client._permissions = A.KitToolPermissions(workspace=workspace)


@pytest.fixture
def pipeline_engine(tmp_path):
    """An engine running a scripted role -- not one of _STEERING_ROLES."""
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")
    _real_perms(eng, tmp_path)
    eng.route = ["builder_agent"]
    eng.current_role_index = 0
    assert eng.current_role not in E._STEERING_ROLES
    return eng


def _finished_job(monkeypatch, jid="4711"):
    """One drained completion event, exactly once, as the store behaves."""
    seen = {"n": 0}

    def _drain(ws):
        seen["n"] += 1
        return [{"job_id": jid, "command": "sbatch opt.slurm", "exit_code": 0,
                 "runtime_s": 3600.0, "stdout_tail": "ORCA TERMINATED NORMALLY",
                 "stderr_tail": ""}] if seen["n"] == 1 else []

    monkeypatch.setattr("delfin.agent.bash_jobs.drain_all_finished_events",
                        _drain)
    monkeypatch.setattr("delfin.agent.job_monitor.check_agent_jobs",
                        lambda ws: [])
    return seen


# ---------------------------------------------------------------------------
# The drain-backed blocks reach every role
# ---------------------------------------------------------------------------

def test_a_scripted_role_is_told_its_job_finished(pipeline_engine, monkeypatch):
    _finished_job(monkeypatch)
    prompt = pipeline_engine._build_current_system_prompt("", task_text="baue")
    assert "4711" in prompt, (
        "a builder role was never told the cluster job it submitted had ended")


def test_a_scripted_role_gets_it_mid_turn_too(pipeline_engine, monkeypatch):
    _finished_job(monkeypatch, jid="8123")
    blocks = pipeline_engine._drain_turn_steering()
    assert any("8123" in b for b in blocks)


def test_the_event_is_still_announced_only_once(pipeline_engine, monkeypatch):
    _finished_job(monkeypatch, jid="9001")
    first = pipeline_engine._drain_turn_steering()
    assert any("9001" in b for b in first)
    assert pipeline_engine._drain_turn_steering() == []


def test_nothing_finished_costs_a_scripted_role_nothing(
        pipeline_engine, monkeypatch):
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_all_finished_events",
                        lambda ws: [])
    monkeypatch.setattr("delfin.agent.job_monitor.check_agent_jobs",
                        lambda ws: [])
    assert pipeline_engine._drain_turn_steering() == []


# ---------------------------------------------------------------------------
# The gate stays where it earns its place
# ---------------------------------------------------------------------------

def test_a_scripted_role_still_gets_no_task_reminder(pipeline_engine):
    """It keeps no list of its own -- the reason the gate exists."""
    from delfin.agent.agent_tasks import get_store
    perms = pipeline_engine.kit_permissions
    get_store(perms.workspace).create(
        "not for this role",
        session_id=getattr(perms, "task_session_id", "") or "")
    assert pipeline_engine._drain_turn_steering() == []
    prompt = pipeline_engine._build_current_system_prompt("", task_text="baue")
    assert "not for this role" not in prompt


def test_only_the_drain_backed_keys_lost_the_gate():
    assert set(E._DRAINED_STEERING_KEYS) == {
        "finished_jobs", "finished_subagents", "answered"}
    for key in ("open_tasks", "budget", "context_status", "project_dir",
                "unmet_delegation", "unmet_tasklist"):
        assert key not in E._DRAINED_STEERING_KEYS


def test_the_interactive_roles_keep_everything(tmp_path, monkeypatch):
    """Dropping the gate for some keys must not drop it for the rest."""
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")
    _real_perms(eng, tmp_path)
    from delfin.agent.agent_tasks import get_store
    perms = eng.kit_permissions
    get_store(perms.workspace).create(
        "Wire the parser", session_id=getattr(perms, "task_session_id", "") or "")
    assert "Wire the parser" in eng._build_current_system_prompt("", task_text="go")
