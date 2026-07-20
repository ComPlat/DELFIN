"""MCP tools honour PreToolUse hooks + the audit log.

MCP tools dispatch straight to the registry, skipping _DocToolExecutor.execute()
where hooks and the audit log live. So a user's PreToolUse hook on `bash` never
fired for `mcp__kit-coding__bash`, and code-modifying MCP calls went unaudited.
Now the MCP path runs the same hooks + audit, keyed on the base tool name.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as ac
from delfin.agent import mcp_client as M
from delfin.agent.api_client import KitToolPermissions, create_client, _doc_executor


def _write_hook(ws, matcher, command):
    d = ws / ".delfin"
    d.mkdir(parents=True, exist_ok=True)
    (d / "settings.json").write_text(json.dumps({"hooks": {
        "PreToolUse": [{"matcher": matcher,
                        "hooks": [{"type": "command", "command": command}]}]
    }}), encoding="utf-8")


# --- helper units -----------------------------------------------------------

def test_pre_hook_helper_returns_block_reason(tmp_path):
    _write_hook(tmp_path, "bash", "echo nope >&2; exit 1")
    perms = KitToolPermissions(workspace=tmp_path, mode="default")
    reason = _doc_executor._run_pre_tool_hooks("bash", {"command": "ls"}, perms)
    assert reason is not None


def test_pre_hook_helper_runs_programmatic_callback(tmp_path):
    seen = []
    perms = KitToolPermissions(workspace=tmp_path, mode="default")
    perms.pre_tool_hook = lambda n, a: seen.append((n, a))
    _doc_executor._run_pre_tool_hooks("write_file", {"path": "x"}, perms)
    assert seen and seen[0][0] == "write_file"


def test_post_hook_helper_invokes_audit(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(_doc_executor, "_audit_call",
                        lambda n, a, p, r: calls.append(n))
    perms = KitToolPermissions(workspace=tmp_path, mode="default")
    _doc_executor._run_post_tool_hooks(
        "write_file", {"path": "x"}, perms, '{"ok":1}')
    assert calls == ["write_file"]   # audited under the base name


# --- MCP dispatch integration (fake stream) ---------------------------------

class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def test_mcp_call_blocked_by_pretooluse_hook(tmp_path, monkeypatch):
    _write_hook(tmp_path, "bash", "echo nope >&2; exit 1")
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))

    streams = [
        _Stream([
            _Chunk([_Choice(_Delta(tool_calls=[_TCDelta(
                0, id="c1", name="mcp__kit-coding__bash",
                arguments='{"command": "ls"}')]))]),
            _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
        ]),
        _Stream([_Chunk([_Choice(_Delta(content="done"), finish="stop")],
                        usage=_Usage())]),
    ]
    calls = {"n": 0}

    def _create(**kwargs):
        s = streams[min(calls["n"], len(streams) - 1)]
        calls["n"] += 1
        return s

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    results = [e for e in events if getattr(e, "type", "") == "tool_result"]
    assert results, "no tool_result emitted"
    # The MCP bash call was blocked by the PreToolUse hook, keyed on `bash`.
    assert any("blocked_by_hook" in (getattr(e, "tool_output", "") or "")
               for e in results)
