"""A server that could not be started reported as healthy.

``MCPServer.list_tools`` fails CLOSED by contract: a missing binary, a
server that exits during ``initialize``, a refused HTTP request, a timeout
— each records the cause in ``last_error`` and returns ``[]``. It never
raises. Both consumers decided by catching an exception:

    unreachable_servers()   for name, server in ...:
                                try: server.list_tools(); continue
                                except Exception: ...append

    doctor._check_mcp()     try: srv.list_tools()
                            except Exception: unreachable.append(name)

so neither branch ever ran. Measured on a server whose command does not
exist:

    list_tools          -> []
    last_error          -> "command not found: …"
    unreachable_servers -> []

and the doctor printed "configured + reachable" while the dashboard printed
"0 tools, no problems". The same held for the real doc server started
against a missing index: it exited, ``last_error`` said "EOF from server",
and nothing reported it.

The verdict now comes from ``last_error`` after the call.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import mcp_client as MC


# ---------------------------------------------------------------------------
# The real failure mode: a spawn that fails, with no exception anywhere
# ---------------------------------------------------------------------------

def test_a_missing_binary_is_reported_although_nothing_raised():
    reg = MC.MCPRegistry()
    reg.servers["dead"] = MC.MCPServer(
        name="dead", command="delfin-no-such-binary-xyz", args=[])

    # The premise: it fails without raising, and says why in last_error.
    assert reg.servers["dead"].list_tools() == []
    assert "command not found" in reg.servers["dead"].last_error

    problems = MC.unreachable_servers(reg)
    assert problems, "a server that could not be started was not reported"
    assert "dead" in problems[0]
    assert "command not found" in problems[0]


def test_a_server_that_exits_on_startup_is_reported():
    """The doc server against a missing index used to be exactly this."""
    import sys
    reg = MC.MCPRegistry()
    reg.servers["exits"] = MC.MCPServer(
        name="exits", command=sys.executable,
        args=["-c", "import sys; sys.exit(1)"])
    assert reg.servers["exits"].list_tools() == []
    problems = MC.unreachable_servers(reg)
    assert problems and "exits" in problems[0]


def test_an_http_server_that_cannot_be_reached_is_reported(monkeypatch):
    def _boom(req, timeout=None):
        raise OSError("connection refused")
    monkeypatch.setattr(MC.urllib.request, "urlopen", _boom)
    reg = MC.MCPRegistry()
    reg.servers["remote"] = MC.MCPServer(
        name="remote", command="", transport="http",
        url="https://down.example/mcp")
    assert reg.servers["remote"].list_tools() == []   # still fails closed
    problems = MC.unreachable_servers(reg)
    assert problems and "remote" in problems[0]


# ---------------------------------------------------------------------------
# The legitimate cases still pass
# ---------------------------------------------------------------------------

class _Healthy:
    name = "ok"
    last_error = ""

    def list_tools(self):
        return [MC.MCPTool(server="ok", name="a", description="", schema={})]


class _Silent:
    """Reachable, offers nothing, has nothing to complain about."""
    name = "empty"
    last_error = ""

    def list_tools(self):
        return []


class _Recovered:
    """Answered this time; an error from an earlier call is history."""
    name = "flaky"
    last_error = "timed out on a previous call"

    def list_tools(self):
        return [MC.MCPTool(server="flaky", name="a", description="",
                           schema={})]


@pytest.mark.parametrize("server", [_Healthy(), _Silent(), _Recovered()])
def test_a_working_server_is_not_reported_as_a_problem(server):
    reg = MC.MCPRegistry(servers={server.name: server}, workspace=None,
                         loaded=True)
    assert MC.unreachable_servers(reg) == []


def test_reporting_never_raises():
    class _Wild:
        name = "wild"

        @property
        def last_error(self):
            raise RuntimeError("even this")

        def list_tools(self):
            raise RuntimeError("and this")

    reg = MC.MCPRegistry(servers={"wild": _Wild()}, workspace=None,
                         loaded=True)
    assert isinstance(MC.unreachable_servers(reg), list)


# ---------------------------------------------------------------------------
# The doctor decides the same way, from the same function
# ---------------------------------------------------------------------------

def test_the_doctor_reports_an_unreachable_server(tmp_path, monkeypatch):
    from delfin.agent import doctor

    cfg = {"servers": {"dead": {"command": "delfin-no-such-binary-xyz"}}}
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "mcp_servers.json").write_text(json.dumps(cfg))
    monkeypatch.setattr(MC, "_user_config_path",
                        lambda: tmp_path / "nonexistent.json")
    # Only the configured server, so the built-ins do not spawn here.
    monkeypatch.setattr(MC, "_BUILTIN_SERVERS", {})
    # A workspace's servers are withheld until the user trusts the
    # directory. This test is about REPORTING a server that cannot be
    # started, so it grants trust and then asks what the doctor says.
    from delfin.agent import workspace_trust as _wt
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    _wt.trust_workspace(tmp_path, (_wt.KIND_MCP_SERVERS,), actor="user")

    rows = doctor._check_mcp({"workspace": str(tmp_path), "fast": False})
    assert rows, "the doctor said nothing about MCP"
    assert any(r["status"] != doctor.PASS for r in rows), (
        "a server that cannot be started was reported as reachable")
    assert any("unreachable" in r["detail"] for r in rows)


def test_the_doctor_still_passes_a_reachable_server(tmp_path, monkeypatch):
    import sys
    from delfin.agent import doctor

    script = tmp_path / "srv.py"
    script.write_text(
        "import json,sys\n"
        "for line in sys.stdin:\n"
        "    msg = json.loads(line)\n"
        "    if 'id' not in msg:\n"
        "        continue\n"
        "    if msg['method'] == 'initialize':\n"
        "        r = {'protocolVersion': '2025-06-18', 'capabilities': {}}\n"
        "    elif msg['method'] == 'tools/list':\n"
        "        r = {'tools': [{'name': 'echo', 'inputSchema': "
        "{'type': 'object'}}]}\n"
        "    else:\n"
        "        r = {}\n"
        "    sys.stdout.write(json.dumps("
        "{'jsonrpc': '2.0', 'id': msg['id'], 'result': r}) + '\\n')\n"
        "    sys.stdout.flush()\n",
        encoding="utf-8")
    cfg = {"servers": {"live": {"command": sys.executable,
                                "args": [str(script)]}}}
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "mcp_servers.json").write_text(json.dumps(cfg))
    monkeypatch.setattr(MC, "_user_config_path",
                        lambda: tmp_path / "nonexistent.json")
    monkeypatch.setattr(MC, "_BUILTIN_SERVERS", {})
    from delfin.agent import workspace_trust as _wt
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    _wt.trust_workspace(tmp_path, (_wt.KIND_MCP_SERVERS,), actor="user")

    rows = doctor._check_mcp({"workspace": str(tmp_path), "fast": False})
    assert rows and all(r["status"] == doctor.PASS for r in rows), rows


# ---------------------------------------------------------------------------
# Reported as dead is not enough: it has to say WHY
# ---------------------------------------------------------------------------
#
# The reason a broken built-in went unnoticed for as long as it did. The
# child was spawned with ``stderr=subprocess.DEVNULL``, so when mcp 2.0.0
# removed ``mcp.server.fastmcp`` and all three built-ins exited 1 on their
# import line, the traceback naming the moved module was discarded and the
# caller was told "server not running" — which reads as a server nobody
# configured rather than as an install that cannot work.

def _dies_with(script: str, name: str = "dying") -> "MC.MCPServer":
    import sys
    return MC.MCPServer(name=name, command=sys.executable, args=["-c", script])


def test_a_server_that_dies_on_an_import_says_which_import():
    """The live defect, in the shape it actually had."""
    srv = _dies_with("import mcp.server.fastmcp_gone_xyz")
    try:
        assert srv.list_tools() == []
        # The FIRST call, not a second one: this is the pass that runs at
        # startup and whose verdict reaches the doctor and the /mcp listing.
        assert "ModuleNotFoundError" in srv.last_error, srv.last_error
        assert "fastmcp_gone_xyz" in srv.last_error, srv.last_error
        assert "exited with code 1" in srv.last_error, srv.last_error
    finally:
        srv.stop()


def test_the_reason_survives_into_the_report_the_user_reads():
    srv = _dies_with("import sys; print('boom: the index is missing',"
                     " file=sys.stderr); sys.exit(2)")
    reg = MC.MCPRegistry()
    reg.servers["dying"] = srv
    try:
        problems = MC.unreachable_servers(reg)
        assert problems, "a server that exited was not reported at all"
        assert "boom: the index is missing" in problems[0], problems
    finally:
        srv.stop()


def test_a_healthy_server_that_writes_to_stderr_is_not_a_failure():
    """Draining stderr must not turn logging into an error.

    An SDK server logs its startup banner there as a matter of course, so a
    rule of "anything on stderr is a problem" would report every working
    server as broken.
    """
    import sys
    srv = MC.MCPServer(name="chatty", command=sys.executable, args=["-c", (
        "import json, sys\n"
        "print('starting up', file=sys.stderr, flush=True)\n"
        "for line in sys.stdin:\n"
        "    msg = json.loads(line)\n"
        "    print('served a request', file=sys.stderr, flush=True)\n"
        "    if msg.get('id') is None:\n"
        "        continue\n"
        "    r = ({'tools': [{'name': 'echo', 'inputSchema': "
        "{'type': 'object'}}]} if msg['method'] == 'tools/list' else {})\n"
        "    sys.stdout.write(json.dumps({'jsonrpc': '2.0', "
        "'id': msg['id'], 'result': r}) + '\\n')\n"
        "    sys.stdout.flush()\n")])
    reg = MC.MCPRegistry()
    reg.servers["chatty"] = srv
    try:
        assert [t.name for t in srv.list_tools()] == ["echo"]
        assert srv.last_error == "", srv.last_error
        assert MC.unreachable_servers(reg) == []
    finally:
        srv.stop()


def test_a_server_that_floods_stderr_is_held_and_reported_bounded():
    """A server is an arbitrary program and may write without end.

    Measured: 200,000 lines of 500 characters — about 100 MB — offered to a
    client that must neither buffer it nor block the child by refusing to
    read it.
    """
    srv = _dies_with(
        "import sys\n"
        "for _ in range(200000): print('x' * 500, file=sys.stderr)\n"
        "sys.exit(3)", name="flood")
    try:
        assert srv.list_tools() == []
        assert "exited with code 3" in srv.last_error
        assert len(srv.last_error) < 4000, len(srv.last_error)
        assert len(srv._stderr_tail_lines) <= MC._STDERR_KEEP_LINES
    finally:
        srv.stop()


def test_a_restart_does_not_report_the_previous_process_dying_words():
    """One server object outlives the processes it spawns."""
    import sys
    srv = _dies_with("import sys; print('the first failure',"
                     " file=sys.stderr); sys.exit(1)")
    try:
        srv.list_tools()
        assert "the first failure" in srv.last_error
        srv.stop()

        srv.args = ["-c", "import sys; print('a different failure',"
                          " file=sys.stderr); sys.exit(1)"]
        srv.command = sys.executable
        srv.list_tools()
        assert "a different failure" in srv.last_error
        assert "the first failure" not in srv.last_error, srv.last_error
    finally:
        srv.stop()


def test_a_stopped_server_is_not_blamed_on_the_output_of_a_healthy_one():
    """The failure mode the stderr capture could have introduced.

    A stopped server has no exit status to attach output to — the handle
    was dropped deliberately, not lost to a death — and its ring still
    holds whatever it logged while it was working. Reporting the startup
    banner of a server the user shut down, as the reason a later call
    failed, would be a worse answer than the bare one it replaced.
    """
    import sys
    srv = MC.MCPServer(name="banner", command=sys.executable, args=["-c", (
        "import sys, time\n"
        "print('INFO  server listening on stdio', file=sys.stderr,"
        " flush=True)\n"
        "time.sleep(30)\n")])
    srv.start()
    try:
        assert srv.proc is not None
        srv.stop()
        resp = srv._send("tools/list")
        assert resp["error"]["message"] == "server not running", resp
    finally:
        srv.stop()
