"""A program may start a calculation, and the gate judges every one.

Part B. Code Mode lets the model write ONE program instead of one turn
per tool call, and the case it pays off on is the multi-step one -- the
field's own rule of thumb. A four-step workflow was seven model rounds
through individual calls; as a program it is one.

The two runners are side-effecting and allowed here for ONE reason:
the gate judges them now. A call from a program reaches
``KitToolExecutor.execute`` -- the same single public entry a direct
call takes, with the deny rules, plan mode, sandbox, role and shape
checks ahead of it -- and from there
``api_client._mcp_calculation_verdict`` ->
``calc_policy.judge_application_run``, which applies the SAME check a
``bash orca job.inp`` gets.

THE ORDER MATTERS AND IS THE POINT. Before that check existed this
would have been a hole: the verdict returned None for run_application
by design, which in an unattended session meant "refused for want of a
dialog", so a program could not have run one -- and the day it could,
it would have run unjudged. The gate first, then the door.

WHAT MAKES A LOOP SAFE, measured rather than assumed: the gate judges
EVERY call. Three good calls produced three verdicts and three bad ones
three refusals; nothing is cached, so a program iterating over ten
applications is judged ten times.
"""
from __future__ import annotations

import tempfile

import pytest

from delfin.agent import api_client as ac
from delfin.agent import calc_policy as cp
from delfin.agent import code_mode as cm

PIPE = "mcp__delfin-tools__"
APP = "conformer_dG"


@pytest.fixture()
def gate():
    ex = ac._DocToolExecutor.__new__(ac._DocToolExecutor)
    ex._classifier_second_opinion = lambda *a, **k: None
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    perms.confirm_callback = None

    def ask(args):
        return ex._gate_side_effecting_mcp_tool(
            PIPE + "run_application", "run_application", args, perms)

    return ask


# --- the door ---------------------------------------------------------

@pytest.mark.parametrize("tool", ["run_application", "submit_application"])
def test_a_runner_is_allowed_from_a_program(tool):
    assert cm._check_allowed(PIPE + tool) is None
    assert PIPE + tool in cm.RUNNER_TOOLS


def test_the_runners_are_the_only_side_effecting_tools_allowed():
    """Named as a set so an addition has to be deliberate."""
    side_effecting = {t for t in cm.READ_ONLY_TOOLS
                      if t.startswith(PIPE)} - cm.READ_ONLY_PIPELINE_TOOLS
    assert side_effecting == cm.RUNNER_TOOLS


@pytest.mark.parametrize("tool", ["register_module", "save_application",
                                  "cancel_run"])
def test_a_side_effecting_tool_that_is_not_a_calculation_stays_refused(tool):
    """judge_application_run does not judge these, so nothing would
    check them."""
    assert cm._check_allowed(PIPE + tool) is not None
    assert ac._mcp_calculation_verdict(tool, {"name": APP}) is None


@pytest.mark.parametrize("tool", ["write_file", "bash", "apply_patch"])
def test_a_native_write_tool_stays_refused(tool):
    assert cm._check_allowed(tool) is not None


# --- the gate still decides ------------------------------------------

def test_a_checked_run_is_admitted(gate):
    assert gate({"name": APP, "cores": 1}) is None


def test_a_run_over_the_limit_is_refused_with_its_reason(gate):
    err = gate({"name": APP, "cores": 9999})
    assert err is not None and "9999 cores" in err


def test_an_unreadable_application_is_refused(gate):
    assert gate({"name": "no-such-application", "cores": 1}) is not None


# --- what makes a loop safe ------------------------------------------

def test_the_gate_judges_every_call_and_caches_nothing(gate, monkeypatch):
    """A program that iterates must be judged per iteration. An
    admission remembered across calls would mean one check and ten
    runs."""
    seen: list = []
    real = cp.judge_application_run

    def counting(*a, **k):
        seen.append(a[:2])
        return real(*a, **k)

    monkeypatch.setattr(cp, "judge_application_run", counting)
    for cores in (1, 1, 1):
        gate({"name": APP, "cores": cores})
    assert len(seen) == 3, f"{len(seen)} verdicts for 3 calls"


def test_a_refusal_is_not_remembered_as_an_admission(gate):
    """Three bad calls, three refusals -- not one refusal and two
    passes."""
    errs = [gate({"name": APP, "cores": 9999}) for _ in range(3)]
    assert all(e is not None for e in errs)


# --- both enforcement points hold ------------------------------------

def test_the_server_side_check_refuses_a_fabricated_payload():
    """The RpcServer checks the allow-list BEFORE the executor, so a
    payload naming a write tool is refused whatever the client stub
    claims."""
    class _Executor:
        def __init__(self):
            self.calls = []

        def execute(self, name, arguments, permissions=None):
            self.calls.append(name)
            return "{}"

    ex = _Executor()
    server = cm.RpcServer(ex)
    out = server.handle_payload('{"tool": "write_file", "arguments": {}}')
    assert "write_file" in out
    assert ex.calls == [], "the executor was touched by a refused call"


def test_an_allowed_runner_reaches_the_one_mcp_dispatcher():
    """A permitted MCP name goes to ``api_client.execute_mcp_call`` --
    the one dispatcher the agent loop uses -- and NOT to
    ``executor.execute``, which answers "Unknown tool" for every MCP
    name. That was the whole defect: the allow-list let the name past
    and the executor refused it as unknown.
    """
    import inspect
    src = inspect.getsource(cm.Connection.handle)
    assert "execute_mcp_call(" in src
    assert 'name.startswith("mcp__")' in src, (
        "an MCP name must be routed by its namespace, not by a list "
        "that would have to be kept in step with the allow-list")


# --- the refusal text tells the truth now ----------------------------

def test_the_refusal_no_longer_claims_code_mode_is_read_only():
    """It was "Code Mode is read-only (phase 1)", which stopped being
    true with this change. A refusal that misstates the rule sends the
    reader looking for a workaround."""
    refusal = cm._check_allowed(PIPE + "register_module")
    assert "read-only (phase 1)" not in refusal
    assert "run_application" in refusal
    assert "judged by the gate" in refusal
