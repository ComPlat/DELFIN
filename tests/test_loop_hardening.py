"""Tool-loop robustness: a single bad tool call must not crash the whole turn.

Drives the real OpenAIClient ReAct loop with a faked OpenAI stream (no network)
and asserts the loop degrades gracefully instead of raising out of the
generator (which would discard every round's progress).
"""

from __future__ import annotations

import pytest

from delfin.agent import mcp_client as M


# --- Fake OpenAI streaming primitives (mirrors test_mcp_meta_tools) ---------

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


def _tool_round(tool_name, args_json, tc_id="c1"):
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TCDelta(0, id=tc_id, name=tool_name, arguments=args_json)]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final():
    return _Stream([
        _Chunk([_Choice(_Delta(content="done"), finish="stop")],
               usage=_Usage()),
    ])


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture
def client(monkeypatch, tmp_path):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import create_client
    return create_client(backend="api", provider="ollama",
                         model="qwen2.5-coder:7b", cwd=str(tmp_path))


def _drive(client, streams):
    captured = []
    calls = {"n": 0}

    def _create(**kwargs):
        captured.append(kwargs)
        s = streams[min(calls["n"], len(streams) - 1)]
        calls["n"] += 1
        return s

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    return events, captured


def _tool_results(events):
    return [e for e in events if getattr(e, "type", "") == "tool_result"]


# --- A5: valid-JSON-but-not-an-object args must not crash -------------------

@pytest.mark.parametrize("bad_args", ['"true"', "[1,2,3]", "42", '"a string"'])
def test_non_object_tool_args_do_not_crash_turn(client, bad_args):
    # Before the fix: json.loads → non-dict → handler .get() → AttributeError
    # escapes the generator and the whole turn dies.
    events, _ = _drive(client, [_tool_round("read_file", bad_args), _final()])
    results = _tool_results(events)
    assert results, "tool_result was never produced (turn crashed?)"
    # The turn completed to the final text.
    assert any(getattr(e, "type", "") == "text_delta"
               and "done" in (getattr(e, "text", "") or "") for e in events)


# --- A4: a handler exception must degrade to a recoverable error ------------

def test_tool_handler_exception_does_not_crash_turn(client, monkeypatch):
    from delfin.agent import api_client as ac

    def _boom(name, args, permissions=None):
        raise RuntimeError("handler blew up")

    monkeypatch.setattr(ac._doc_executor, "execute", _boom)
    events, _ = _drive(client, [_tool_round("read_file", "{}"), _final()])
    results = _tool_results(events)
    assert results, "no tool_result — the exception escaped the turn"
    assert any("failed" in (getattr(e, "tool_output", "") or "").lower()
               for e in results)


# --- A6: an empty streamed tool_call id is backfilled + kept consistent -----

def test_empty_tool_call_id_is_backfilled(client):
    _, captured = _drive(
        client, [_tool_round("read_file", "{}", tc_id=None), _final()])
    assert len(captured) >= 2, "second round never happened"
    round2_msgs = captured[1].get("messages") or []
    asst = [m for m in round2_msgs
            if m.get("role") == "assistant" and m.get("tool_calls")]
    tool = [m for m in round2_msgs if m.get("role") == "tool"]
    assert asst and tool, "assistant/tool messages missing from round 2"
    tc_id = asst[0]["tool_calls"][0]["id"]
    assert tc_id, "tool_call id was left empty"
    # the tool result must reference the SAME (backfilled) id
    assert any(m.get("tool_call_id") == tc_id for m in tool)
