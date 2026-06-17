"""Tests for the MCP Streamable-HTTP transport.

No live HTTP server: ``urllib.request.urlopen`` is monkeypatched with a fake
MCP endpoint that speaks JSON-RPC and can answer either as ``application/json``
or ``text/event-stream`` (SSE). Covers transport selection, the initialize +
Mcp-Session-Id lifecycle, tools/list + tools/call over both response shapes,
and registry wiring.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M


# ---------------------------------------------------------------------------
# Fake HTTP MCP endpoint
# ---------------------------------------------------------------------------


class _FakeResp:
    def __init__(self, body: bytes, headers: dict):
        self._body = body
        self.headers = headers

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _result_for(method: str) -> dict:
    if method == "initialize":
        return {"protocolVersion": "2025-06-18",
                "serverInfo": {"name": "fake", "version": "1"},
                "capabilities": {"tools": {}}}
    if method == "tools/list":
        return {"tools": [{"name": "echo", "description": "Echo it",
                           "inputSchema": {"type": "object"}}]}
    if method == "tools/call":
        return {"content": [{"type": "text", "text": "hello-from-http"}]}
    if method == "resources/list":
        return {"resources": [
            {"uri": "file:///notes.md", "name": "Notes",
             "description": "Project notes", "mimeType": "text/markdown"},
            {"uri": "mem://state", "name": "State"},
        ]}
    if method == "resources/read":
        return {"contents": [
            {"uri": "file:///notes.md", "mimeType": "text/markdown",
             "text": "# Notes\nhello"},
        ]}
    if method == "prompts/list":
        return {"prompts": [
            {"name": "summarize", "description": "Summarize text",
             "arguments": [{"name": "text", "required": True}]},
        ]}
    if method == "prompts/get":
        return {"description": "Summarize",
                "messages": [
                    {"role": "user",
                     "content": {"type": "text", "text": "Summarize this: X"}},
                ]}
    return {}


def make_fake_endpoint(*, sse: bool = False, session: str = "sess-123"):
    """Return (urlopen_fn, calls) where calls records (method, session_hdr)."""
    calls: list[tuple] = []

    def _urlopen(req, timeout=None):
        body = json.loads(req.data.decode("utf-8"))
        method = body.get("method", "")
        # urllib title-cases header keys → "Mcp-session-id"
        sess_hdr = req.headers.get("Mcp-session-id")
        calls.append((method, sess_hdr))
        if "id" not in body:                       # notification → 202, no body
            return _FakeResp(b"", {})
        payload = {"jsonrpc": "2.0", "id": body["id"],
                   "result": _result_for(method)}
        headers = {"Mcp-Session-Id": session}
        if sse:
            headers["Content-Type"] = "text/event-stream"
            data = f"event: message\ndata: {json.dumps(payload)}\n\n".encode()
            return _FakeResp(data, headers)
        headers["Content-Type"] = "application/json"
        return _FakeResp(json.dumps(payload).encode("utf-8"), headers)

    return _urlopen, calls


def _http_server(monkeypatch, **kw):
    fake, calls = make_fake_endpoint(**kw)
    monkeypatch.setattr(M.urllib.request, "urlopen", fake)
    srv = M.MCPServer(name="remote", command="", transport="http",
                      url="https://mcp.example.com/mcp")
    return srv, calls


# ---------------------------------------------------------------------------
# Transport selection from config
# ---------------------------------------------------------------------------


def test_server_from_config_http_by_type():
    srv = M._server_from_config("r", {"type": "http",
                                      "url": "https://x/mcp",
                                      "headers": {"Authorization": "Bearer t"}})
    assert srv.transport == "http"
    assert srv.url == "https://x/mcp"
    assert srv.headers["Authorization"] == "Bearer t"


def test_server_from_config_http_by_url_only():
    srv = M._server_from_config("r", {"url": "https://x/mcp"})
    assert srv.transport == "http"


def test_server_from_config_stdio_default():
    srv = M._server_from_config("s", {"command": "echo", "args": ["hi"]})
    assert srv.transport == "stdio"
    assert srv.command == "echo"


# ---------------------------------------------------------------------------
# initialize + tools over HTTP (JSON responses)
# ---------------------------------------------------------------------------


def test_http_initialize_and_list_tools(monkeypatch):
    srv, calls = _http_server(monkeypatch)
    tools = srv.list_tools()
    assert [t.name for t in tools] == ["echo"]
    assert tools[0].namespaced_name == "mcp__remote__echo"
    # The handshake happened: initialize, then the initialized notification.
    methods = [m for m, _ in calls]
    assert methods[0] == "initialize"
    assert "notifications/initialized" in methods


def test_http_session_id_captured_and_echoed(monkeypatch):
    srv, calls = _http_server(monkeypatch, session="abc-789")
    srv.list_tools()
    assert srv.session_id == "abc-789"
    # The initialize call had no session yet; later calls must echo it.
    init_call = next(c for c in calls if c[0] == "initialize")
    later_calls = [c for c in calls if c[0] == "tools/list"]
    assert init_call[1] is None
    assert later_calls and later_calls[0][1] == "abc-789"


def test_http_call_tool(monkeypatch):
    srv, _ = _http_server(monkeypatch)
    out = srv.call_tool("echo", {"text": "x"})
    assert out == "hello-from-http"


# ---------------------------------------------------------------------------
# SSE responses
# ---------------------------------------------------------------------------


def test_http_list_tools_via_sse(monkeypatch):
    srv, _ = _http_server(monkeypatch, sse=True)
    tools = srv.list_tools()
    assert [t.name for t in tools] == ["echo"]


def test_http_call_tool_via_sse(monkeypatch):
    srv, _ = _http_server(monkeypatch, sse=True)
    assert srv.call_tool("echo", {}) == "hello-from-http"


def test_extract_jsonrpc_from_sse_matches_id():
    raw = (
        "event: message\n"
        'data: {"jsonrpc":"2.0","id":7,"result":{"ok":true}}\n\n'
    )
    msg = M._extract_jsonrpc_from_sse(raw, 7)
    assert msg["result"]["ok"] is True


def test_extract_jsonrpc_from_sse_no_payload():
    msg = M._extract_jsonrpc_from_sse("event: ping\n\n", 1)
    assert "error" in msg


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def test_http_request_failure_is_graceful(monkeypatch):
    def _boom(req, timeout=None):
        raise OSError("connection refused")
    monkeypatch.setattr(M.urllib.request, "urlopen", _boom)
    srv = M.MCPServer(name="r", command="", transport="http",
                      url="https://down.example/mcp")
    assert srv.list_tools() == []          # fails closed, no exception
    assert "http request failed" in srv.last_error


# ---------------------------------------------------------------------------
# Registry wiring
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Resources + prompts
# ---------------------------------------------------------------------------


def test_http_list_and_read_resources(monkeypatch):
    srv, _ = _http_server(monkeypatch)
    resources = srv.list_resources()
    assert [r.uri for r in resources] == ["file:///notes.md", "mem://state"]
    assert resources[0].name == "Notes"
    assert resources[0].mime_type == "text/markdown"
    body = srv.read_resource("file:///notes.md")
    assert "hello" in body


def test_http_list_and_get_prompts(monkeypatch):
    srv, _ = _http_server(monkeypatch)
    prompts = srv.list_prompts()
    assert [p.name for p in prompts] == ["summarize"]
    assert prompts[0].namespaced_name == "mcp__remote__summarize"
    assert prompts[0].arguments[0]["name"] == "text"
    rendered = srv.get_prompt("summarize", {"text": "X"})
    assert "Summarize this: X" in rendered


def test_registry_resources_and_prompts(monkeypatch, tmp_path):
    fake, _ = make_fake_endpoint()
    monkeypatch.setattr(M.urllib.request, "urlopen", fake)
    monkeypatch.setattr(M, "_user_config_path",
                        lambda: tmp_path / "none.json")
    cfg = {"servers": {"remote": {"url": "https://mcp.example.com/mcp"}}}
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "mcp_servers.json").write_text(json.dumps(cfg))
    reg = M.MCPRegistry()
    reg.load(workspace=tmp_path)
    assert [r.uri for r in reg.discover_resources()] == \
        ["file:///notes.md", "mem://state"]
    assert [p.namespaced_name for p in reg.discover_prompts()] == \
        ["mcp__remote__summarize"]
    assert "hello" in reg.read_resource("remote", "file:///notes.md")
    assert "Summarize this: X" in reg.get_prompt(
        "mcp__remote__summarize", {"text": "X"})


def test_registry_loads_http_and_discovers(monkeypatch, tmp_path):
    fake, _ = make_fake_endpoint()
    monkeypatch.setattr(M.urllib.request, "urlopen", fake)
    cfg = {"servers": {"remote": {"type": "http",
                                  "url": "https://mcp.example.com/mcp"},
                       # disable the built-in default so this stays focused on
                       # the HTTP server under test (no real subprocess spawn).
                       "delfin-tools": {"enabled": False}}}
    (tmp_path / ".delfin").mkdir()
    (tmp_path / ".delfin" / "mcp_servers.json").write_text(json.dumps(cfg))
    monkeypatch.setattr(M, "_user_config_path",
                        lambda: tmp_path / "nonexistent.json")
    reg = M.MCPRegistry()
    reg.load(workspace=tmp_path)
    assert "remote" in reg.servers
    assert reg.servers["remote"].transport == "http"
    tools = reg.discover_all()
    assert [t.namespaced_name for t in tools] == ["mcp__remote__echo"]
    # And dispatch routes back through the HTTP transport.
    assert reg.call("mcp__remote__echo", {}) == "hello-from-http"
