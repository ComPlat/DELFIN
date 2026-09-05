"""The stall watchdog blamed the provider for the agent's own tool.

The dashboard ends a turn after 120 s without a stream event and tells
the user the endpoint went quiet. But a dispatched tool call produces no
stream event either, and every blocking tool has a budget far longer than
120 s: ``bash_status`` blocks server-side for its wait cap, a foreground
``bash`` runs to the permission object's maximum, a sub-agent runs to its
wall budget -- and ANY tool can be parked at a confirmation dialog for as
long as the broker waits for the user.

So a tool that was working correctly was reported as a dead provider, the
turn was killed underneath it, and the approval the user gave a moment
later was consumed by a broker whose turn no longer existed.

The silence is now attributed: while a dispatched tool is inside its own
deadline the kill waits for THAT deadline, and when a tool does blow its
budget the message names the tool instead of the provider.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

from delfin.dashboard import tab_agent as T

_SRC = pathlib.Path(inspect.getfile(T)).read_text(encoding="utf-8")


def _fn(name: str) -> ast.AST:
    for node in ast.walk(ast.parse(_SRC)):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found")


def _calls(node: ast.AST) -> set[str]:
    return {getattr(n.func, "id", "") or getattr(n.func, "attr", "")
            for n in ast.walk(node) if isinstance(n, ast.Call)}


# ---------------------------------------------------------------------------
# the budgets come from the tools, not from a guess
# ---------------------------------------------------------------------------

def test_a_confirmation_prompt_outlives_the_stall_budget():
    """Any tool can sit at the dialog for the broker's whole timeout; the
    120 s stall budget killed the turn a third of the way in."""
    from delfin.agent.kit_confirm import KitConfirmBroker

    wait = inspect.signature(
        KitConfirmBroker.__init__).parameters["default_timeout_s"].default
    assert T._tool_stall_budget_s("write_file") >= wait
    assert T._tool_stall_budget_s("write_file") > 120.0


def test_a_bash_status_wait_gets_its_server_side_cap():
    from delfin.agent.api_client import _BASH_STATUS_WAIT_CAP_S

    assert T._tool_stall_budget_s("bash_status") > _BASH_STATUS_WAIT_CAP_S


def test_a_foreground_bash_gets_the_permission_objects_maximum():
    from delfin.agent.api_client import KitToolPermissions

    assert (T._tool_stall_budget_s("bash")
            > KitToolPermissions.bash_max_timeout_s)


def test_a_subagent_gets_its_wall_budget():
    from delfin.agent.subagents import _MAX_WALL_S

    assert T._tool_stall_budget_s("subagent") > _MAX_WALL_S


def test_the_namespaced_name_is_the_same_tool():
    """The KIT backend emits mcp__kit-coding__bash_status; a budget table
    that only knows the bare name would give it the floor."""
    assert (T._tool_stall_budget_s("mcp__kit-coding__bash_status")
            == T._tool_stall_budget_s("bash_status"))


# ---------------------------------------------------------------------------
# what is in flight
# ---------------------------------------------------------------------------

def test_nothing_in_flight_excuses_nothing():
    assert T._longest_pending_tool({}, 100.0) is None
    assert T._longest_pending_tool(None, 100.0) is None


def test_the_tool_with_the_most_budget_left_decides():
    """Parallel calls: the fan-out's sub-agent must not be judged by the
    read_file that started beside it."""
    now = 10_000.0
    name, _waited, _budget = T._longest_pending_tool(
        {"read_file": now - 10.0, "subagent": now - 10.0}, now)
    assert name == "subagent"


def test_a_tool_over_its_own_budget_stops_excusing_the_silence():
    now = 10_000.0
    _name, waited, budget = T._longest_pending_tool(
        {"read_file": now - 5_000.0}, now)
    assert waited > budget


def test_a_junk_timestamp_is_ignored_rather_than_raising():
    assert T._longest_pending_tool({"read_file": None}, 10_000.0) is None


# ---------------------------------------------------------------------------
# the decision
# ---------------------------------------------------------------------------

def test_a_running_tool_holds_the_kill_off():
    """120 s of stall budget, 30 s spent, and a confirmation dialog open:
    the turn used to die here."""
    pending = ("write_file", 130.0, T._tool_stall_budget_s("write_file"))
    assert T._stall_remaining_s(elapsed=130.0, budget=120.0,
                                pending=pending) > 0


def test_the_wait_is_the_tools_deadline_not_the_providers():
    pending = ("bash_status", 10.0, 600.0)
    assert T._stall_remaining_s(
        elapsed=10.0, budget=120.0, pending=pending) == 590.0


def test_silence_with_no_tool_running_is_still_a_stall():
    assert T._stall_remaining_s(elapsed=130.0, budget=120.0,
                                pending=None) < 0


def test_a_tool_past_its_deadline_no_longer_protects_the_turn():
    pending = ("read_file", 5_000.0, 330.0)
    assert T._stall_remaining_s(elapsed=5_000.0, budget=120.0,
                                pending=pending) < 0


# ---------------------------------------------------------------------------
# ...and it is wired into both halves of the watchdog
# ---------------------------------------------------------------------------

def test_both_watchdog_halves_ask_what_is_running():
    for name in ("_check_kill", "_check_stale"):
        assert "_longest_pending_tool" in _calls(_fn(name)), name


def test_the_kill_waits_for_the_tool_deadline():
    assert "_stall_remaining_s" in _calls(_fn("_check_kill"))


def test_the_dispatch_and_the_result_bracket_the_tool():
    """Without both ends the marker is either never set or never cleared,
    and 'never cleared' hands the rest of the turn the tool's budget."""
    assert '_tool_inflight' in ast.dump(_fn("_on_tool_use"))
    assert '_tool_inflight' in ast.dump(_fn("_on_tool_result"))


def test_the_marker_is_cleared_before_an_empty_result_returns():
    """_on_tool_result returns early when the result body is empty; the
    clear has to happen first or the tool stays 'running' forever."""
    body = _fn("_on_tool_result").body
    clear_at = min(i for i, stmt in enumerate(body)
                   if "_tool_inflight" in ast.dump(stmt))
    return_at = min(i for i, stmt in enumerate(body)
                    if isinstance(stmt, ast.If)
                    and any(isinstance(s, ast.Return) for s in stmt.body))
    assert clear_at < return_at


def test_a_streaming_provider_clears_the_marker():
    """Tokens mean the tools are answered -- a result event with an empty
    body reaches no callback at all."""
    for name in ("_on_token", "_on_thinking"):
        assert "_tool_inflight" in ast.dump(_fn(name)), name


def test_the_kill_message_names_the_tool_rather_than_the_provider():
    assert "the tool" in _SRC and "This is the tool, not the" in _SRC
