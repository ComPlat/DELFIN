"""``check_only`` on an MCP call is a word in the arguments, not a dry run.

Input: the agent's MCP gate (``_gate_mcp_tool`` and the side-effect gate
behind it). Both returned "allowed" for any MCP tool whose arguments
carried ``check_only: true`` -- in plan mode and when nobody could be
asked -- on the reasoning that a dry run changes nothing. No tool of the
built-in servers (delfin-ops, delfin-tools) has such a parameter, and
an MCP server drops an argument its schema does not name (FastMCP's
argument model ignores extra keys), so the call ran for real. Measured
2026-10-22 on this branch: ``mcp__delfin-ops__submit_calculation`` with
``check_only: true`` was refused without the flag and allowed with it,
in plan mode and in default mode with no dialog. With the host grant on
the built-in ops server, that is a submission nobody approved.

After: an MCP call is judged on its name, whatever its arguments say.
The native ``apply_patch`` keeps its exemption, because its executor
runs ``git apply --check`` and writes nothing when the flag is set.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A

CALLS = [
    "mcp__delfin-ops__submit_calculation",
    "mcp__delfin-ops__run_orca_input",
    "mcp__delfin-tools__run_application",
    "mcp__github__create_pull_request",
]


@pytest.fixture
def gate(tmp_path):
    def run(name, mode, args):
        perms = A.KitToolPermissions(workspace=tmp_path, mode=mode,
                                     agent_role="solo_agent")
        perms.confirm_callback = None
        return A._doc_executor._gate_mcp_tool(name, args, perms)
    return run


@pytest.mark.parametrize("name", CALLS)
@pytest.mark.parametrize("mode", ["plan", "default"])
def test_a_check_only_flag_does_not_open_the_gate(gate, name, mode):
    args = {"folder": "/tmp/x", "check_only": True}
    assert gate(name, mode, args) is not None


@pytest.mark.parametrize("name", CALLS)
def test_plan_mode_on_the_dispatch_path_ignores_the_flag(tmp_path, name):
    perms = A.KitToolPermissions(workspace=tmp_path, mode="plan",
                                 agent_role="solo_agent")
    out = A._doc_executor.execute(name, {"check_only": True}, perms)
    assert "plan mode" in out


def test_the_native_patch_check_is_still_a_dry_run(tmp_path):
    perms = A.KitToolPermissions(workspace=tmp_path, mode="plan",
                                 agent_role="solo_agent")
    out = A._doc_executor.execute(
        "apply_patch", {"diff": "--- a/x\n+++ b/x\n", "check_only": True},
        perms)
    assert "plan mode (read-only)" not in out
