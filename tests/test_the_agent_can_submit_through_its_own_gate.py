"""DELFIN's own job tools are reachable, and every call is the user's decision.

Input: the agent's built-in registration of the ``delfin-ops`` MCP server
and the agent's MCP gate. Field report 2026-10-22 (KIT cluster): the user
authorised the cluster jobs in the prompt; ``submit_calculation``,
``run_orca_input`` and ``run_application`` all answered
``mutation_blocked`` and the DFT refinement never ran.

The ops server refuses a mutating tool unless the process that started it
set ``DELFIN_OPS_ALLOW_MUTATE`` -- consent from the host, never from a tool
argument (tests/test_a_mutating_ops_tool_cannot_grant_itself_permission.py).
Nothing in DELFIN ever set it, so no route existed at all. The agent is
that host, and its gate decides every non-read-only MCP call before it
reaches the server: the approval dialog in the default modes, a refusal
when nobody can be asked, a refusal in plan mode, and in bypassPermissions
the user's standing grant with the classifier's second opinion.

After: the agent's built-in ``delfin-ops`` is started with the grant, and
this file pins that the gate stands in front of every mutating ops tool.
"""
from __future__ import annotations

import pytest

from delfin.agent import mcp_client

MUTATING = [
    "submit_calculation", "cancel_calculation", "rename_calc_folder",
    "create_calc_folder", "move_calc_folder", "move_to_archive",
    "delete_calc_folder", "kill_all_user_jobs", "prepare_recalc",
    "run_calc_option", "cleanup", "stop", "pipeline_prepare",
    "pipeline_run", "run_orca_input", "co2", "tadf_xtb", "hyperpol",
]


def test_the_agents_ops_server_is_started_with_the_grant():
    env = mcp_client._BUILTIN_SERVERS["delfin-ops"].get("env") or {}
    assert env.get("DELFIN_OPS_ALLOW_MUTATE") == "1"
    # Only the ops server: nothing else is handed the grant.
    for name, cfg in mcp_client._BUILTIN_SERVERS.items():
        if name != "delfin-ops":
            assert "DELFIN_OPS_ALLOW_MUTATE" not in (cfg.get("env") or {})


def test_the_grant_reaches_the_server_process():
    from delfin.agent.api_client import _configured_child_env
    env = _configured_child_env(mcp_client._BUILTIN_SERVERS["delfin-ops"]
                                .get("env"))
    assert env.get("DELFIN_OPS_ALLOW_MUTATE") == "1"


@pytest.fixture
def gate(tmp_path, monkeypatch):
    from delfin.agent import api_client as A
    executor = A._doc_executor
    seen = []
    monkeypatch.setattr(executor, "_classifier_second_opinion",
                        lambda perms, name, args, mode: seen.append(name))

    def run(name, mode, confirm=None):
        perms = A.KitToolPermissions(workspace=tmp_path, mode=mode,
                                     agent_role="solo_agent")
        perms.confirm_callback = confirm
        return executor._gate_mcp_tool(f"mcp__delfin-ops__{name}",
                                       {"folder": str(tmp_path)}, perms)
    run.seen = seen
    return run


@pytest.mark.parametrize("name", MUTATING)
def test_with_nobody_to_ask_a_mutating_ops_call_is_refused(gate, name):
    assert gate(name, "default") is not None


@pytest.mark.parametrize("name", MUTATING)
def test_in_plan_mode_a_mutating_ops_call_is_refused(gate, name):
    assert gate(name, "plan") is not None


@pytest.mark.parametrize("name", MUTATING)
def test_the_user_is_asked_and_their_answer_decides(gate, name):
    asked = []
    yes = gate(name, "default", confirm=lambda *a, **k: asked.append(1) or True)
    no = gate(name, "default", confirm=lambda *a, **k: asked.append(1) or False)
    assert yes is None and no is not None and len(asked) == 2


def test_in_bypass_the_second_opinion_is_consulted(gate):
    assert gate("submit_calculation", "bypassPermissions") is None
    assert gate.seen == ["mcp__delfin-ops__submit_calculation"]
