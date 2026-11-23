"""The allow-list let the name past and the executor said "Unknown tool".

RED FIRST, and this file is why. Two earlier changes added DELFIN's
pipeline tools to Code Mode's allow-list, their own tests passed, and
the capability did not exist: ``_DocToolExecutor.execute`` dispatches
NATIVE tools only. Measured before the refactor:

    execute("mcp__delfin-tools__get_guide")      -> "Unknown tool"
    execute("mcp__delfin-tools__resolve_spec")   -> "Unknown tool"
    execute("mcp__delfin-tools__run_application")-> "Unknown tool"
    execute("search_calcs")                      -> works

MCP tools are dispatched in the AGENT LOOP, not in the executor, and
that branch says so itself: "MCP tools skip _DocToolExecutor.execute(),
so run the SAME PreToolUse hooks + content-gate + audit/PostToolUse
here". Four steps the RPC server never had.

So the tests here assert a CALL, never a membership. A test that checks
an allow-list proves that a name is in a set; this one proves a program
gets an answer. That is the difference between the two attempts, and
the reason this file is the first commit of the change rather than the
last.
"""
from __future__ import annotations

import json
import tempfile

import pytest

from delfin.agent import api_client as ac
from delfin.agent import code_mode as cm

PIPE = "mcp__delfin-tools__"


@pytest.fixture()
def server():
    """An RpcServer over the REAL executor and real permissions.

    Not a recording double: a double would answer whatever it is asked
    and prove nothing about dispatch, which is exactly how the no-op
    passed its tests.
    """
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    perms.confirm_callback = None
    return cm.RpcServer(ac._DocToolExecutor(), permissions=perms)


def _call(server, tool, **arguments):
    return json.loads(server.handle_payload(
        json.dumps({"tool": tool, "arguments": arguments})))


# --- the call, not the membership -------------------------------------

def test_a_program_gets_an_answer_from_get_guide(server):
    """The whole point: an answer, not "Unknown tool"."""
    reply = _call(server, PIPE + "get_guide")
    assert "Unknown tool" not in json.dumps(reply), reply
    assert reply, "an empty reply is not an answer"


def test_a_program_gets_an_answer_from_catalog(server):
    reply = _call(server, PIPE + "catalog")
    assert "Unknown tool" not in json.dumps(reply), reply


def test_a_native_tool_still_answers(server):
    """The route that worked must keep working."""
    reply = _call(server, "search_calcs", query="nothing-matches-this")
    assert "Unknown tool" not in json.dumps(reply), reply


def test_an_unknown_tool_still_says_so(server):
    """The error that proved the no-op must still be reachable for a
    name that really does not exist."""
    reply = _call(server, PIPE + "no_such_tool_at_all")
    assert "error" in reply


# --- the gate travels with the dispatch -------------------------------

def test_a_run_over_the_limit_is_refused_through_the_program_route(server):
    """The same refusal a direct call gets, with its reason. If the
    program route reached the registry without the gate, this would
    start a 9999-core run."""
    reply = _call(server, PIPE + "run_application",
                  name="conformer_dG", cores=9999)
    text = json.dumps(reply)
    assert "9999" in text or "cores" in text, text


def test_a_refused_write_tool_never_reaches_the_dispatcher(server):
    reply = _call(server, "write_file", path="x", content="y")
    assert "error" in reply
    assert "directly" in reply["error"].lower()


# --- and over the real socket, in a real child process ---------------

def _short_workspace():
    import os
    import tempfile
    bases = [b for b in ("/tmp", tempfile.gettempdir()) if os.path.isdir(b)]
    return tempfile.mkdtemp(prefix="cmpipe_", dir=min(bases, key=len))


def test_a_real_program_reaches_the_pipeline_over_the_socket(tmp_path):
    """The documented route, end to end: the hook starts the server, the
    child imports the stub and calls get_guide(), the answer comes back.

    This is the test the two earlier attempts needed and did not have.
    They asserted that a name was in a frozenset; a frozenset cannot
    tell you that a program gets an answer.
    """
    import json as _json
    import os
    import subprocess
    import sys

    workspace = _short_workspace()
    with open(os.path.join(workspace, "tools_api.py"), "w",
              encoding="utf-8") as fh:
        fh.write(cm.generate_stub())
    program = os.path.join(workspace, "prog.py")
    with open(program, "w", encoding="utf-8") as fh:
        fh.write("import json\n"
                 "import tools_api as api\n"
                 "out = api.get_guide()\n"
                 "print('UNKNOWN' if 'Unknown tool' in out else 'ANSWERED')\n")

    perms = ac.KitToolPermissions(mode="default", workspace=workspace)
    perms.confirm_callback = None
    # conftest.child_env, not os.environ: a child a test starts resolves
    # its own ~/.delfin sinks at import time against whatever HOME says,
    # and the suite's in-memory redirects stop at the process boundary.
    # Measured 2026-09-22: nine confirmation requests from suite children
    # landed in the operator's real ~/.delfin/terminal_confirmations. The
    # repo's guard test caught this one before it shipped.
    from conftest import child_env
    env = dict(child_env(tmp_path))
    env.update(ac._code_mode_child_env(
        "python3 -m delfin.agent.code_mode run prog.py",
        workspace, perms, ac._DocToolExecutor()) or {})
    if "DELFIN_CODE_MODE_RPC" not in env:
        pytest.fail("the hook exported no socket; the route is not wired")
    env["PYTHONPATH"] = os.pathsep.join(
        [workspace, env.get("PYTHONPATH", "")]).strip(os.pathsep)

    proc = subprocess.run(
        [sys.executable, "-m", "delfin.agent.code_mode", "run", program],
        cwd=workspace, env=env, capture_output=True, text=True, timeout=180)
    assert proc.returncode == 0, proc.stderr[-800:]
    assert "ANSWERED" in proc.stdout, (proc.stdout, proc.stderr[-400:])


# --- one dispatcher means one behaviour ------------------------------

def test_the_loop_and_the_rpc_call_the_same_dispatcher():
    """Extracted, not copied. Two dispatch paths would mean the gate
    holds on one of them, and which one would depend on where a call
    came from."""
    import inspect
    loop = inspect.getsource(ac.KitAgentClient) \
        if hasattr(ac, "KitAgentClient") else inspect.getsource(ac)
    assert "execute_mcp_call(" in loop
    rpc = inspect.getsource(cm.Connection.handle)
    assert "execute_mcp_call(" in rpc
    # and the four steps live in one place
    disp = inspect.getsource(ac.execute_mcp_call)
    for step in ("_mcp_call_refusal", "get_registry", "_wrap_untrusted",
                 "_run_post_tool_hooks"):
        assert step in disp, step


def test_a_blocking_hook_ends_a_call_from_a_program(monkeypatch):
    """The PreToolUse hooks answer FIRST, from a program as from the
    loop. If they did not, a hook written to stop `bash` would stop it
    in a turn and not in a program."""
    calls: list = []

    def _refuse(name, args, permissions):
        calls.append(name)
        return json.dumps({"error": "a hook blocked this"})

    monkeypatch.setattr(ac, "_mcp_call_refusal", _refuse)
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    srv = cm.RpcServer(ac._DocToolExecutor(), permissions=perms)
    reply = json.loads(srv.handle_payload(json.dumps(
        {"tool": PIPE + "get_guide", "arguments": {}})))
    assert "a hook blocked this" in json.dumps(reply)
    assert calls == [PIPE + "get_guide"]


def test_the_post_hooks_see_the_un_namespaced_name(monkeypatch):
    """An audit rule on `bash` must cover `mcp__kit-coding__bash`, which
    is why the base name is what the hooks get -- on both routes."""
    seen: list = []

    class _Exec(ac._DocToolExecutor):
        def _run_post_tool_hooks(self, base, args, perms, result):
            seen.append(base)

    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    perms.confirm_callback = None
    srv = cm.RpcServer(_Exec(), permissions=perms)
    srv.handle_payload(json.dumps(
        {"tool": PIPE + "get_guide", "arguments": {}}))
    assert seen == ["get_guide"], seen


def test_a_result_from_our_own_server_is_still_wrapped_as_data():
    """The distinction this change keeps: a tool DESCRIPTION is this
    repository's text and is advertised plainly (first-party), and a
    RESULT is data about the world whoever served it."""
    import inspect
    assert "_wrap_untrusted" in inspect.getsource(ac.execute_mcp_call)
