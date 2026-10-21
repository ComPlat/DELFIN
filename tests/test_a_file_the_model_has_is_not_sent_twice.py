"""A file the model already has, unchanged and in view, is not sent twice.

Session 78c36f8c (2026-10-19): one 441-line HTML file read whole again
and again, the requests growing past 400k input tokens. The short answer
is given only while the earlier result row is still in the request, byte
for byte, and the file's size and mtime are unchanged.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import mcp_client as M
from delfin.agent.reread import ReadLedger


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "game.html").write_text("<script>\nconst BX=4;\n</script>\n")
    return tmp_path


def _row(cid, content):
    return {"role": "tool", "tool_call_id": cid, "content": content}


def test_an_unchanged_read_still_in_view_gets_the_note(ws):
    led = ReadLedger()
    led.record({"path": "game.html"}, ws, "c1", "BODY")
    note = led.repeat_note({"path": "game.html"}, ws, [_row("c1", "BODY")])
    assert note and "unchanged" in note and "game.html" in note


def test_an_absolute_and_a_relative_path_are_the_same_file(ws):
    led = ReadLedger()
    led.record({"path": "game.html"}, ws, "c1", "BODY")
    assert led.repeat_note({"path": str(ws / "game.html")}, ws,
                           [_row("c1", "BODY")])


def test_a_changed_file_is_read_again(ws):
    led = ReadLedger()
    led.record({"path": "game.html"}, ws, "c1", "BODY")
    (ws / "game.html").write_text("<script>\nconst BX=10;\n</script>\n")
    assert led.repeat_note({"path": "game.html"}, ws, [_row("c1", "BODY")]) \
        is None


def test_a_same_size_edit_in_the_same_second_is_read_again(ws):
    led = ReadLedger()
    led.record({"path": "game.html"}, ws, "c1", "BODY")
    st = os.stat(ws / "game.html")
    (ws / "game.html").write_text("<script>\nconst BX=5;\n</script>\n")
    os.utime(ws / "game.html", ns=(st.st_atime_ns, st.st_mtime_ns + 1))
    assert led.repeat_note({"path": "game.html"}, ws, [_row("c1", "BODY")]) \
        is None


@pytest.mark.parametrize("messages", [
    [],                                   # dropped by a compaction
    [_row("c1", "BODY[... trimmed]")],    # rewritten by a context edit
])
def test_a_result_no_longer_in_view_is_read_again(ws, messages):
    led = ReadLedger()
    led.record({"path": "game.html"}, ws, "c1", "BODY")
    assert led.repeat_note({"path": "game.html"}, ws, messages) is None


def test_other_lines_are_read(ws):
    led = ReadLedger()
    led.record({"path": "game.html"}, ws, "c1", "BODY")
    assert led.repeat_note({"path": "game.html", "offset": 2}, ws,
                           [_row("c1", "BODY")]) is None


# -- through the tool loop --------------------------------------------------

class _Fn:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TC:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _Fn(name, arguments)


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


def _read_call(cid):
    return _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TC(
        0, id=cid, name="read_file",
        arguments=json.dumps({"path": "game.html"}))]), finish="tool_calls")],
        usage=_Usage())])


def _say(text):
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


@pytest.fixture
def client(monkeypatch, ws):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import KitToolPermissions, create_client
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(ws))
    cl.set_permissions(KitToolPermissions(
        workspace=ws, mode="default", task_session_id="sess-rr"))
    return cl


def _tool_rows(sent):
    return [str(m.get("content")) for m in sent[-1] if m.get("role") == "tool"]


def test_the_second_read_in_a_turn_gets_the_note(client):
    sent = []
    streams = [_read_call("c1"), _read_call("c2"), _say("BX ist 4.")]

    def _create(**kw):
        sent.append(list(kw.get("messages") or []))
        return streams[min(len(sent) - 1, len(streams) - 1)]

    client.client.chat.completions.create = _create
    list(client.stream_message(
        "sys", [{"role": "user", "content": "was ist BX?"}], max_tokens=100))
    rows = _tool_rows(sent)
    assert len(rows) == 2
    assert "const BX=4" in rows[0]
    assert "unchanged" in rows[1] and "const BX=4" not in rows[1]
