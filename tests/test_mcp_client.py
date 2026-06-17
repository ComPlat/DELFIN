"""Tests for delfin.agent.mcp_client (MCP stdio client)."""

from __future__ import annotations

import json
import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

from delfin.agent import mcp_client as M


# A tiny in-process MCP-shaped server: reads JSON-RPC lines from stdin,
# answers initialize / tools/list / tools/call.
_FAKE_SERVER = textwrap.dedent('''
    import json, sys
    while True:
        line = sys.stdin.readline()
        if not line:
            break
        line = line.strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except Exception:
            continue
        method = req.get("method")
        rid = req.get("id")
        if method == "initialize":
            sys.stdout.write(json.dumps({
                "jsonrpc": "2.0",
                "id": rid,
                "result": {
                    "protocolVersion": "2025-06-18",
                    "capabilities": {},
                    "serverInfo": {"name": "fake", "version": "0.1"},
                },
            }) + "\\n"); sys.stdout.flush()
        elif method == "notifications/initialized":
            pass    # no response
        elif method == "tools/list":
            sys.stdout.write(json.dumps({
                "jsonrpc": "2.0",
                "id": rid,
                "result": {"tools": [
                    {"name": "echo",
                     "description": "Echo input",
                     "inputSchema": {
                         "type": "object",
                         "properties": {"text": {"type": "string"}},
                     }},
                    {"name": "add",
                     "description": "Add a+b",
                     "inputSchema": {"type": "object"}},
                ]},
            }) + "\\n"); sys.stdout.flush()
        elif method == "tools/call":
            tname = req["params"]["name"]
            args = req["params"].get("arguments", {})
            if tname == "echo":
                payload = args.get("text", "")
            elif tname == "add":
                payload = str(int(args.get("a", 0)) + int(args.get("b", 0)))
            else:
                payload = "?"
            sys.stdout.write(json.dumps({
                "jsonrpc": "2.0",
                "id": rid,
                "result": {"content": [{"type": "text", "text": payload}]},
            }) + "\\n"); sys.stdout.flush()
''')


@pytest.fixture
def fake_server_script():
    with tempfile.NamedTemporaryFile(
        "w", suffix=".py", delete=False,
    ) as f:
        f.write(_FAKE_SERVER)
        f.flush()
        yield Path(f.name)


@pytest.fixture(autouse=True)
def _reset_registry():
    M.reset_registry()
    yield
    M.reset_registry()


def test_load_configs_merges_user_and_project(monkeypatch):
    with tempfile.TemporaryDirectory() as user_d, \
         tempfile.TemporaryDirectory() as proj_d:
        user_path = Path(user_d) / "mcp_servers.json"
        user_path.write_text(json.dumps({
            "servers": {"user_srv": {"command": "echo"}},
        }))
        proj_path = Path(proj_d) / ".delfin" / "mcp_servers.json"
        proj_path.parent.mkdir(parents=True)
        proj_path.write_text(json.dumps({
            "servers": {"proj_srv": {"command": "echo"}},
        }))
        monkeypatch.setattr(M, "_user_config_path", lambda: user_path)
        configs = M._load_configs(Path(proj_d))
        assert "user_srv" in configs
        assert "proj_srv" in configs


def test_load_configs_skips_disabled(monkeypatch):
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "mcp_servers.json"
        path.write_text(json.dumps({
            "servers": {
                "off": {"command": "echo", "enabled": False},
                "on": {"command": "echo"},
            },
        }))
        monkeypatch.setattr(M, "_user_config_path", lambda: path)
        configs = M._load_configs(None)
        assert "on" in configs and "off" not in configs
        # built-in DELFIN servers are present by default alongside user config
        assert "delfin-tools" in configs


def test_builtin_delfin_tools_default_and_disable(monkeypatch):
    """delfin-tools ships as a built-in default (so Pipeline mode works without
    manual per-user setup) but can be disabled via the user config."""
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "mcp_servers.json"
        monkeypatch.setattr(M, "_user_config_path", lambda: path)
        assert "delfin-tools" in M._load_configs(None)          # default on
        path.write_text(json.dumps(
            {"servers": {"delfin-tools": {"enabled": False}}}))
        assert "delfin-tools" not in M._load_configs(None)      # disable-able


def test_server_initialize_and_list_tools(fake_server_script):
    srv = M.MCPServer(
        name="test", command=sys.executable,
        args=[str(fake_server_script)],
    )
    try:
        assert srv.initialize()
        tools = srv.list_tools()
        names = sorted(t.name for t in tools)
        assert names == ["add", "echo"]
        assert {t.namespaced_name for t in tools} == {
            "mcp__test__echo", "mcp__test__add",
        }
    finally:
        srv.stop()


def test_server_call_tool_returns_text(fake_server_script):
    srv = M.MCPServer(
        name="t", command=sys.executable,
        args=[str(fake_server_script)],
    )
    try:
        assert srv.initialize()
        result = srv.call_tool("echo", {"text": "hello"})
        assert "hello" in result
        result2 = srv.call_tool("add", {"a": 2, "b": 3})
        assert result2.strip() == "5"
    finally:
        srv.stop()


def test_server_handles_missing_command():
    srv = M.MCPServer(name="x", command="/bin/no-such-command")
    assert not srv.initialize()
    assert srv.proc is None or (srv.last_error and "command" in srv.last_error.lower())


def test_registry_load_and_call(monkeypatch, fake_server_script):
    with tempfile.TemporaryDirectory() as d:
        cfg = Path(d) / "mcp_servers.json"
        cfg.write_text(json.dumps({
            "servers": {
                "fake": {
                    "command": sys.executable,
                    "args": [str(fake_server_script)],
                },
            },
        }))
        monkeypatch.setattr(M, "_user_config_path", lambda: cfg)
        reg = M.get_registry(None)
        try:
            tools = reg.discover_all()
            assert len(tools) >= 2
            names = {t.namespaced_name for t in tools}
            assert "mcp__fake__echo" in names

            out = reg.call("mcp__fake__echo", {"text": "hi"})
            assert "hi" in out
        finally:
            M.reset_registry()


def test_registry_unknown_tool():
    reg = M.MCPRegistry()
    out = reg.call("mcp__unknown__nope", {})
    payload = json.loads(out)
    assert "error" in payload


def test_malformed_namespaced_name():
    reg = M.MCPRegistry()
    out = reg.call("not_mcp", {})
    payload = json.loads(out)
    assert "error" in payload

    out2 = reg.call("mcp__noundsep", {})
    payload = json.loads(out2)
    assert "error" in payload
