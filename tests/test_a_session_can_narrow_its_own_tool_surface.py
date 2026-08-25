"""One session's own tool allow list and deny list, and where they bite.

``--allowed-tools`` used to travel to the subprocess CLI backend alone and
was dropped for every provider that builds an ``OpenAIClient`` — which is
the only client with a tool loop of its own. The flag was announced and
changed nothing on the backend the project runs on.

The fix has one property worth more than the flag: the lists live on the
permissions object the EXECUTOR is handed, so an excluded tool is REFUSED
when called. Filtering only the advertised list would be
advertising-as-security — the model would not see the tool and anything
else that reaches the executor (an MCP backend, a replayed call, a
sub-agent) would still run it. Several tests below call the executor
directly with an excluded name for exactly that reason.

The rest of the file pins the three semantics that are silent when they are
backwards: an empty allow list is no restriction, deny wins over allow, and
neither list can widen what the role already forbids.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli as agent_cli
from delfin.agent.api_client import (
    KitToolPermissions,
    ToolSurfaceContext,
    _DocToolExecutor,
    _session_tool_refusal,
    _tool_denied_for_session,
    advertisable_tools,
    tool_unavailable_reason,
    unknown_tool_names,
)


@pytest.fixture
def workspace(tmp_path) -> Path:
    (tmp_path / "a.txt").write_text("hello\n", encoding="utf-8")
    return tmp_path


def _perms(workspace, **over) -> KitToolPermissions:
    base = dict(workspace=workspace, mode="default")
    base.update(over)
    return KitToolPermissions(**base)


def _error_of(result: str) -> str:
    try:
        return json.loads(result).get("error", "")
    except Exception:
        return ""


def _args(**over) -> argparse.Namespace:
    base = dict(backend="", provider="", model="", mode="", cwd="",
                effort="", permission_mode="", settings_defaults=False,
                allowed_tools="", disallowed_tools="")
    base.update(over)
    return argparse.Namespace(**base)


# ---------------------------------------------------------------------------
# The executor refuses. Not the advertised list — the executor.
# ---------------------------------------------------------------------------

def test_a_denied_tool_is_refused_when_called_not_merely_hidden(workspace):
    """The defect this whole design exists to avoid.

    The executor is called directly, the way an MCP backend, a replayed
    tool call or a sub-agent would reach it — with no advertised list
    anywhere in the path. A deny list that only filtered the schema sent to
    the model would let this call through, and the tool would run.
    """
    ex = _DocToolExecutor()
    open_perms = _perms(workspace)
    assert "hello" in ex.execute("read_file", {"path": "a.txt"}, open_perms)

    denied = _perms(workspace, session_denied_tools=frozenset({"read_file"}))
    result = ex.execute("read_file", {"path": "a.txt"}, denied)
    assert "hello" not in result
    assert "--disallowed-tools" in _error_of(result)


def test_a_tool_left_off_the_allow_list_is_refused_by_the_executor(workspace):
    ex = _DocToolExecutor()
    perms = _perms(workspace, session_allowed_tools=frozenset({"read_file"}))
    assert "hello" in ex.execute("read_file", {"path": "a.txt"}, perms)

    result = ex.execute("write_file", {"path": "b.txt", "content": "x"}, perms)
    assert "--allowed-tools" in _error_of(result)
    assert not (workspace / "b.txt").exists(), (
        "the write ran despite being off the session's allow list")


def test_an_mcp_backend_is_not_a_way_around_the_session_lists(workspace):
    """MCP names are dispatched without passing through ``execute``.

    So ``--disallowed-tools bash`` would stop the native tool and leave
    ``mcp__<server>__bash`` running the identical command. The gate that
    guards that path carries the same check, judged on the bare name.
    """
    ex = _DocToolExecutor()
    perms = _perms(workspace, session_denied_tools=frozenset({"bash"}))
    blocked = ex._gate_mcp_tool("mcp__kit-coding__bash", {"command": "ls"},
                                perms)
    assert blocked and "not a way around" in blocked


def test_a_sub_agent_inherits_the_lists_it_was_started_under(workspace):
    """``dataclasses.replace`` copies field values, which is the mechanism.

    A child that could be spawned without them would be a way around the
    parent's restriction rather than a delegation of it.
    """
    parent = _perms(workspace, session_denied_tools=frozenset({"bash"}),
                    session_allowed_tools=frozenset({"read_file", "bash"}))
    child = dataclasses.replace(parent, subagent_depth=1)
    assert child.session_denied_tools == frozenset({"bash"})
    assert _tool_denied_for_session(child, "bash")


# ---------------------------------------------------------------------------
# The advertised surface is DERIVED from that refusal, never a second policy
# ---------------------------------------------------------------------------

def test_the_advertised_surface_drops_what_the_executor_would_refuse(workspace):
    ctx = ToolSurfaceContext(session_denied=frozenset({"read_file"}))
    names = {t["function"]["name"]
             for t in advertisable_tools(A._DOC_TOOLS_OPENAI, ctx)}
    assert "read_file" not in names
    assert "write_file" in names


def test_an_allow_list_leaves_only_what_it_names_on_the_surface():
    ctx = ToolSurfaceContext(
        session_allowed=frozenset({"read_file", "grep_file"}))
    names = {t["function"]["name"]
             for t in advertisable_tools(A._DOC_TOOLS_OPENAI, ctx)}
    assert names <= {"read_file", "grep_file"}
    assert "read_file" in names


def test_nothing_refused_by_the_executor_is_still_advertised(workspace):
    """The invariant the derivation exists for, over the whole catalogue.

    Advertising must stay a strict SUBSET of what may execute. The one
    permitted asymmetry is the other way round — the meta calls the harness
    needs to finish a turn stay executable under an allow list without
    being advertised, so a shortfall here can only ever be safe.
    """
    allowed = frozenset({"read_file", "grep_file"})
    denied = frozenset({"grep_file"})
    perms = _perms(workspace, session_allowed_tools=allowed,
                   session_denied_tools=denied)
    ctx = ToolSurfaceContext(session_allowed=allowed, session_denied=denied)
    for tool in A._DOC_TOOLS_OPENAI:
        name = tool["function"]["name"]
        refused = _tool_denied_for_session(perms, name) is not None
        hidden = tool_unavailable_reason(name, ctx) is not None
        if refused:
            assert hidden, (
                f"{name} would be refused by the executor and is still "
                "advertised — the model is offered a tool it cannot use")


# ---------------------------------------------------------------------------
# The three semantics that are silent when they are backwards
# ---------------------------------------------------------------------------

def test_an_empty_allow_list_means_no_restriction_not_no_tools(workspace):
    """Read the other way round, this disables the agent without a word."""
    perms = _perms(workspace)
    assert perms.session_allowed_tools == frozenset()
    for name in ("read_file", "write_file", "bash"):
        assert _tool_denied_for_session(perms, name) is None
    names = {t["function"]["name"]
             for t in advertisable_tools(A._DOC_TOOLS_OPENAI,
                                         ToolSurfaceContext())}
    assert "bash" in names and "read_file" in names


def test_deny_wins_when_a_name_is_on_both_lists(workspace):
    perms = _perms(workspace,
                   session_allowed_tools=frozenset({"read_file", "bash"}),
                   session_denied_tools=frozenset({"bash"}))
    assert _tool_denied_for_session(perms, "read_file") is None
    assert "--disallowed-tools" in (_tool_denied_for_session(perms, "bash") or "")


def test_the_lists_narrow_and_never_widen_a_roles_surface(workspace):
    """A flag that could widen a role is a privilege escalation.

    ``dashboard_agent`` may not run ``bash`` at all. Naming it on the
    session allow list must not grant it — the role gate is evaluated
    separately and still refuses.
    """
    perms = _perms(workspace, agent_role="dashboard_agent",
                   session_allowed_tools=frozenset({"bash", "read_file"}))
    result = _DocToolExecutor().execute("bash", {"command": "echo hi"}, perms)
    assert "not available to the 'dashboard_agent' role" in _error_of(result)

    ctx = ToolSurfaceContext(role="dashboard_agent",
                             session_allowed=frozenset({"bash"}))
    assert tool_unavailable_reason("bash", ctx) is not None
    names = {t["function"]["name"]
             for t in advertisable_tools(A._DOC_TOOLS_OPENAI, ctx)}
    assert "bash" not in names


def test_a_deny_list_narrows_a_role_further_than_the_role_does(workspace):
    """The other direction of the same rule: narrowing is always allowed."""
    perms = _perms(workspace, agent_role="dashboard_agent",
                   session_denied_tools=frozenset({"search_docs"}))
    assert _tool_denied_for_session(perms, "search_docs") is not None


def test_an_allow_list_does_not_break_the_calls_the_harness_needs(workspace):
    """``_ALWAYS_ALLOWED_TOOLS`` is exempt from an allow list, as for a role.

    A session that names three working tools did not mean to make it
    impossible to submit a plan or return a sub-agent result. An explicit
    deny still holds over them — that is a decision someone typed out.
    """
    perms = _perms(workspace, session_allowed_tools=frozenset({"read_file"}))
    for meta in sorted(A._ALWAYS_ALLOWED_TOOLS):
        assert _tool_denied_for_session(perms, meta) is None, meta

    denied = _perms(workspace,
                    session_denied_tools=frozenset({"exit_plan_mode"}))
    assert _tool_denied_for_session(denied, "exit_plan_mode") is not None


def test_a_namespaced_call_is_judged_by_the_tool_it_names():
    assert _session_tool_refusal(
        "mcp__kit-coding__bash", denied=frozenset({"bash"})) is not None
    assert _session_tool_refusal(
        "mcp__kit-coding__read_file",
        allowed=frozenset({"read_file"})) is None


# ---------------------------------------------------------------------------
# A name that matches no tool is a typo the user hears about
# ---------------------------------------------------------------------------

def test_a_name_that_matches_no_tool_is_reported():
    stray = unknown_tool_names(["read_file", "reed_file", "bahs"])
    assert stray == ["bahs", "reed_file"]


def test_an_mcp_name_is_never_called_a_typo():
    """MCP servers are discovered inside the turn; a guess would be a false
    alarm about a name that is very likely correct."""
    assert unknown_tool_names(["mcp__github__create_issue"]) == []


def test_the_startup_notice_names_the_stray_tool():
    engine = type("E", (), {
        "client": type("C", (), {"model": "m"})(),
        "provider": "kit", "backend": "api"})()
    notes = agent_cli._bounding_notices(
        _args(allowed_tools="read_file,reed_file"), engine)
    assert any("reed_file" in n and "typo" in n for n in notes)


# ---------------------------------------------------------------------------
# The flags, and the client that now takes them
# ---------------------------------------------------------------------------

def test_both_flags_are_offered_on_the_command_line():
    args = agent_cli.build_parser().parse_args(
        ["chat", "--allowed-tools", "read_file,bash",
         "--disallowed-tools", "web_fetch"])
    assert args.allowed_tools == "read_file,bash"
    assert args.disallowed_tools == "web_fetch"


def test_the_flag_values_reach_the_object_the_executor_reads(tmp_path):
    """End to end onto the permissions, which is where enforcement lives."""
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client._permissions = _perms(tmp_path)
    engine = type("E", (), {"client": client, "provider": "kit",
                            "backend": "api"})()
    agent_cli._apply_tool_surface(
        engine, _args(allowed_tools="read_file, bash ",
                      disallowed_tools="bash"))
    assert client._permissions.session_allowed_tools == {"read_file", "bash"}
    assert client._permissions.session_denied_tools == {"bash"}
    assert _tool_denied_for_session(client._permissions, "bash") is not None


def test_a_client_without_a_tool_loop_stores_nothing_and_says_so(tmp_path):
    """The plain OpenAI provider builds no permissions object, so there is
    nothing to narrow — and the run has to be told rather than left
    believing the surface was cut."""
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client._permissions = None
    engine = type("E", (), {"client": client, "provider": "openai",
                            "backend": "api"})()
    agent_cli._apply_tool_surface(engine, _args(allowed_tools="read_file"))
    assert client.enforced_tool_surface() == ((), ())
    notes = agent_cli._bounding_notices(
        _args(allowed_tools="read_file"), engine)
    assert any("--allowed-tools REQUESTED but not applied" in n for n in notes)


def test_the_cannot_enforce_notice_no_longer_fires_for_the_tool_loop(tmp_path):
    """The backend this project runs on CAN enforce it now, so the notice
    that said otherwise must stop firing — a stale warning trains the user
    to ignore the ones that are true."""
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client._permissions = _perms(tmp_path)
    engine = type("E", (), {"client": client, "provider": "kit",
                            "backend": "api"})()
    args = _args(allowed_tools="read_file", disallowed_tools="bash")
    agent_cli._apply_tool_surface(engine, args)
    assert agent_cli._bounding_notices(args, engine) == []


def test_the_subprocess_backend_keeps_its_allow_list_and_admits_no_deny_list():
    """Its tool loop runs in a child process, so the only lever is the
    command line — an allow list it has, a deny list it has not."""
    client = A.CLIClient.__new__(A.CLIClient)
    client.allowed_tools = ["read_file", "bash"]
    engine = type("E", (), {"client": client, "provider": "",
                            "backend": "cli"})()
    notes = agent_cli._bounding_notices(
        _args(allowed_tools="read_file,bash"), engine)
    assert notes == []
    notes = agent_cli._bounding_notices(
        _args(allowed_tools="read_file,bash", disallowed_tools="bash"), engine)
    assert any("--disallowed-tools REQUESTED but not applied" in n
               for n in notes)


def test_a_run_that_asks_for_nothing_narrows_nothing(tmp_path):
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client._permissions = _perms(tmp_path)
    engine = type("E", (), {"client": client, "provider": "kit",
                            "backend": "api"})()
    agent_cli._apply_tool_surface(engine, _args())
    assert client.enforced_tool_surface() == ((), ())
    assert agent_cli._bounding_notices(_args(), engine) == []


def test_whitespace_only_lists_are_no_restriction(tmp_path):
    """`--allowed-tools " , "` is a user who typed nothing, not a user who
    asked for a session with no tools at all."""
    perms = _perms(tmp_path, session_allowed_tools=" , ")
    assert perms.session_allowed_tools == frozenset()
    assert _tool_denied_for_session(perms, "bash") is None
