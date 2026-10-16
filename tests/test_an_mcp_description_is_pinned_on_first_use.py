"""An MCP tool's description is pinned on first use; a change waits for
the user.

Input: an MCP server (a fake registry) whose tools carry harmless
descriptions, discovered by the real client loop (the tool surface of
each request), across sessions that share one trust store under
``tmp_path``. Between sessions a description, a schema or the tool list
changes. Output: the ``tools`` of each request, the questions the
session asks, the security events, the MCP gate's answer to a call, and
``delfin-agent mcp-pins``.

Semantics:

* the first time a server is seen its tools are pinned: a SHA-256 of
  each tool's description and input schema, kept in the trust store
  (owner-only, outside the workspace), recorded as ``mcp_pinned``;
* a pinned tool whose hash is unchanged is advertised as before;
* a changed description or schema, or a new tool on a known server, is
  withheld from the surface and its calls are refused, until the user
  accepts it: an attended session asks once per change and process
  (``mcp_description_accepted`` or ``mcp_description_changed``); an
  unattended one withholds it with a blocked ``mcp_description_changed``
  event, and ``delfin-agent mcp-pins accept`` accepts it from the
  person's shell;
* DELFIN's own built-in servers are not pinned (their text ships with
  DELFIN).
"""
from __future__ import annotations

import json
import stat
from pathlib import Path
from unittest.mock import patch

import pytest

import delfin.agent.api_client as A
from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc
from delfin.agent import security_events as SE

MARKER = "DELFIN_AGENT_SESSION"


class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None
        self.reasoning_content = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.index = 0
        self.delta = delta
        self.finish_reason = finish


class _Chunk:
    def __init__(self, choices):
        self.choices = choices
        self.usage = None


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _Registry:
    def __init__(self, tools=(), sources=None):
        self._tools = list(tools)
        self._sources = dict(sources or {})

    def discover_all(self):
        return list(self._tools)

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []

    def source_of(self, name):
        return self._sources.get(name, "/home/u/.delfin/mcp_servers.json")


def _tool(desc="Look up a compound by name.", name="lookup", server="chem",
          schema=None):
    return M.MCPTool(server=server, name=name, description=desc,
                     schema=schema or {"type": "object", "properties": {
                         "name": {"type": "string"}}})


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.delenv(MARKER, raising=False)
    from delfin.agent import trust_store
    monkeypatch.setattr(trust_store, "store_dir",
                        lambda: tmp_path / "trust_store")
    trust_store.forget_cached_keys()
    _forget()
    SE.clear()
    w = tmp_path / "ws"
    w.mkdir()
    yield w
    _forget()
    SE.clear()


def _session(ws, tools, *, mode="bypassPermissions", answer=None,
             sources=None, requests=2):
    """One session: ``requests`` requests, each assembling the surface.
    Returns ({tool name: description} per request, questions asked)."""
    asked: list = []
    callback = None
    if answer is not None:
        def callback(name, args, preview):
            asked.append((name, args, preview))
            return answer
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    registry = _Registry(tools, sources)
    patches = [patch.object(mc, "resolve", lambda *a, **k: caps),
               patch.object(M, "get_registry", lambda *a, **k: registry),
               patch.object(A.time, "sleep", lambda *a, **k: None)]
    for p in patches:
        p.start()
    seen: list = []
    try:
        client = A.create_client(backend="api", provider="ollama",
                                 model="qwen2.5-coder:7b", cwd=str(ws))
        perms = A.KitToolPermissions(workspace=ws, mode=mode,
                                     confirm_callback=callback)
        client.set_permissions(perms)

        def _create(**kwargs):
            seen.append(kwargs)
            return _Stream([_Chunk([_Choice(_Delta(content="done"),
                                            finish="stop")])])
        client.client.chat.completions.create = _create
        for _ in range(requests):
            list(client.stream_message(
                "sys", [{"role": "user", "content": "go"}], max_tokens=64))
    finally:
        for p in patches:
            p.stop()
    out = []
    for kw in seen:
        out.append({t["function"]["name"]: t["function"]["description"]
                    for t in (kw.get("tools") or [])
                    if t["function"]["name"].startswith("mcp__")})
    return out, asked, perms


def _forget():
    """A new process: the in-memory answers of the last one are gone."""
    getattr(M, "forget_pin_decisions", lambda: None)()


def _kinds():
    return [e.kind for e in SE.recent(200)]


def test_first_use_pins_and_advertises(ws, tmp_path):
    surfaces, asked, _ = _session(ws, [_tool()])
    assert all("mcp__chem__lookup" in s for s in surfaces)
    assert asked == []
    assert "mcp_pinned" in _kinds()
    store = tmp_path / "trust_store" / "mcp_pins.json"
    assert store.exists()
    assert stat.S_IMODE(store.stat().st_mode) == 0o600
    assert not str(store).startswith(str(ws))
    data = json.loads(store.read_text())
    assert "lookup" in json.dumps(data)


def test_an_unchanged_description_is_advertised_without_a_question(ws):
    _session(ws, [_tool()])
    _forget()
    surfaces, asked, _ = _session(ws, [_tool()], mode="default", answer=True)
    assert all("mcp__chem__lookup" in s for s in surfaces)
    assert asked == []


def test_a_changed_description_is_withheld_unattended(ws):
    _session(ws, [_tool()])
    _forget()
    SE.clear()
    changed = _tool("Look up a compound by its CAS number.")
    surfaces, asked, perms = _session(ws, [changed])
    assert all("mcp__chem__lookup" not in s for s in surfaces), surfaces
    events = [e for e in SE.recent(200)
              if e.kind == "mcp_description_changed"]
    assert len(events) == 1 and events[0].blocked
    refusal = A._doc_executor._gate_mcp_tool(
        "mcp__chem__lookup", {"name": "x"}, perms)
    assert refusal and "accept" in refusal, refusal


@pytest.mark.parametrize("answer", [True, False])
def test_an_attended_session_asks_once_per_change(ws, answer):
    _session(ws, [_tool()])
    _forget()
    SE.clear()
    changed = _tool("Look up a compound by its CAS number.")
    surfaces, asked, _ = _session(ws, [changed], mode="default",
                                  answer=answer, requests=3)
    assert len(asked) == 1, asked
    assert "CAS number" in asked[0][2]
    if answer:
        assert all("CAS number" in s["mcp__chem__lookup"] for s in surfaces)
        assert "mcp_description_accepted" in _kinds()
        _forget()
        again, asked2, _ = _session(ws, [changed], mode="default",
                                    answer=True)
        assert asked2 == [] and "mcp__chem__lookup" in again[0]
    else:
        assert all("mcp__chem__lookup" not in s for s in surfaces)
        assert "mcp_description_changed" in _kinds()


def test_a_new_tool_on_a_known_server_waits(ws):
    _session(ws, [_tool()])
    _forget()
    added = _tool("List the stored spectra.", name="spectra")
    surfaces, _, _ = _session(ws, [_tool(), added])
    assert all("mcp__chem__lookup" in s for s in surfaces)
    assert all("mcp__chem__spectra" not in s for s in surfaces)


def test_a_schema_change_counts_as_a_change(ws):
    _session(ws, [_tool()])
    _forget()
    wider = _tool(schema={"type": "object", "properties": {
        "name": {"type": "string", "description": "the compound name"}}})
    surfaces, _, _ = _session(ws, [wider])
    assert all("mcp__chem__lookup" not in s for s in surfaces)


def test_delfins_own_servers_are_not_pinned(ws):
    """A server whose entry is DELFIN's own (source "built-in"); the
    name is neutral because the default tool set hides delfin-tools."""
    own = _tool(server="own", name="calc")
    _session(ws, [own], sources={"own": "built-in"})
    _forget()
    changed = _tool("A newer text of the shipped tool.", server="own",
                    name="calc")
    surfaces, asked, _ = _session(ws, [changed], mode="default", answer=False,
                                  sources={"own": "built-in"})
    assert asked == []
    assert all("mcp__own__calc" in s for s in surfaces), surfaces


def test_the_person_accepts_a_change_from_the_shell(ws, capsys, monkeypatch):
    from delfin.agent import cli as agent_cli
    _session(ws, [_tool()])
    _forget()
    changed = _tool("Look up a compound by its CAS number.")
    _session(ws, [changed])
    assert agent_cli.main(["mcp-pins"]) == 0
    listing = capsys.readouterr().out
    assert "chem" in listing and "lookup" in listing and "changed" in listing
    monkeypatch.setenv(MARKER, "99")
    assert agent_cli.main(["mcp-pins", "accept", "chem", "lookup"]) == 2
    monkeypatch.delenv(MARKER)
    assert agent_cli.main(["mcp-pins", "accept", "chem", "lookup"]) == 0
    _forget()
    surfaces, _, _ = _session(ws, [changed])
    assert all("CAS number" in s["mcp__chem__lookup"] for s in surfaces)


def test_the_event_kinds_have_labels():
    for kind in ("mcp_pinned", "mcp_description_changed",
                 "mcp_description_accepted"):
        assert kind in SE.known_kinds()
