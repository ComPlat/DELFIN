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
