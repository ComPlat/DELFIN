"""A message typed into a running turn stays in the history.

Input: a turn whose tool round is followed by a typed steer and, in the
same drain, a harness row (here a run note). The loop appends both as
``user`` rows, so the turn's rows hold two consecutive user rows.

Semantics pinned here: the next turn's alternation repair
(``AgentEngine._sanitize_messages``) merges consecutive user rows into
one row, the parts joined with a blank line, in order, each intact. It
used to keep the newer row and drop the older one, so the typed steer
was gone from the next request, from ``engine.messages`` and from every
saved session after it.

Driven through the real OpenAI-compatible tool loop and the real engine
turn; only ``chat.completions.create`` is a stub.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import mcp_client as M


STEER = "nimm bitte nur mod1, nicht mod2"
NOTE = "[system] The permission mode changed to plan while you worked."


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


def _tool_round(path):
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TCDelta(0, id="c1", name="read_file",
                     arguments=json.dumps({"path": path}))]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final(text):
    return _Stream([
        _Chunk([_Choice(_Delta(content=text), finish="stop")],
               usage=_Usage()),
    ])


@pytest.fixture
def engine(tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "mod1.py").write_text("x = 1\n")
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_finished_events",
                        lambda ws: [])
    from delfin.agent.api_client import KitToolPermissions, create_client
    from delfin.agent.engine import AgentEngine
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(ws))
    cl.set_permissions(KitToolPermissions(
        workspace=ws, mode="default", task_session_id="sess-1"))
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = AgentEngine(repo_dir=ws, backend="api", provider="ollama",
                          model="qwen2.5-coder:7b", mode="solo")
    eng.client = cl
    return eng, ws


def _drive(eng, ws):
    """Turn 1: a tool round, during which a steer is typed and a run note
    queued; then an answer. Turn 2: a plain question. Returns every
    request's messages."""
    requests: list[list[dict]] = []
    script = [_tool_round(str(ws / "mod1.py")), _final("mod1 gelesen."),
              _final("ok")]

    def _create(**kwargs):
        requests.append([dict(m) for m in kwargs.get("messages") or []])
        if len(requests) == 1:
            assert eng.steer(STEER)
            eng.client.push_run_note(NOTE)
        return script[len(requests) - 1]

    eng.client.client.chat.completions.create = _create
    eng.stream_response("lies die Module")
    after_turn = [dict(m) for m in eng.messages]
    eng.stream_response("und jetzt?")
    return requests, after_turn


def _texts(rows):
    return "\n".join(str(m.get("content") or "") for m in rows
                     if m.get("role") == "user")


def test_the_typed_steer_reaches_the_next_turns_request(engine):
    eng, ws = engine
    requests, _after = _drive(eng, ws)
    assert len(requests) == 3
    assert STEER in _texts(requests[1]) and NOTE in _texts(requests[1])
    assert STEER in _texts(requests[2]), (
        "the typed steer was dropped from the next turn's request")
    assert NOTE in _texts(requests[2])


def test_the_typed_steer_stays_in_the_history(engine):
    eng, ws = engine
    _requests, after_turn = _drive(eng, ws)
    assert STEER in _texts(after_turn)
    assert STEER in _texts(eng.messages)
    assert NOTE in _texts(eng.messages)


def test_the_merged_row_keeps_both_parts_in_order(engine):
    eng, ws = engine
    _drive(eng, ws)
    row = next(m for m in eng.messages
               if m.get("role") == "user" and STEER in str(m.get("content")))
    text = row["content"]
    assert NOTE in text and text.index(STEER) < text.index(NOTE)
    roles = [m.get("role") for m in eng.messages]
    assert all(not (a == b == "user") for a, b in zip(roles, roles[1:]))


def test_the_tool_round_stays_paired(engine):
    eng, ws = engine
    _drive(eng, ws)
    msgs = eng.messages
    i = next(i for i, m in enumerate(msgs) if m.get("tool_calls"))
    assert msgs[i + 1]["role"] == "tool"
    assert msgs[i + 1]["tool_call_id"] == msgs[i]["tool_calls"][0]["id"]


def test_the_typed_steer_survives_save_and_restore(engine, tmp_path):
    eng, ws = engine
    _drive(eng, ws)
    state = json.loads(json.dumps(eng.export_state(), default=str))
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        other = AgentEngine(repo_dir=ws, backend="api", provider="ollama",
                            model="qwen2.5-coder:7b", mode="solo")
    other.restore_state(state)
    assert STEER in _texts(other.messages)
    other._sanitize_messages()
    assert STEER in _texts(other.messages)


# -- the repair itself -----------------------------------------------------

def _bare():
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    return eng


def test_an_alternating_history_is_left_as_the_same_objects():
    eng = _bare()
    hist = [{"role": "user", "content": "q"},
            {"role": "assistant", "content": "",
             "tool_calls": [{"id": "c1", "type": "function",
                             "function": {"name": "read_file",
                                          "arguments": "{}"}}]},
            {"role": "tool", "tool_call_id": "c1", "content": "x"},
            {"role": "user", "content": "steer"},
            {"role": "assistant", "content": "a"}]
    eng.messages = list(hist)
    before = json.dumps(hist)
    eng._sanitize_messages()
    assert all(a is b for a, b in zip(eng.messages, hist))
    assert json.dumps(eng.messages) == before


def test_list_content_is_merged_as_parts():
    eng = _bare()
    img = {"type": "image_url", "image_url": {"url": "data:x"}}
    eng.messages = [{"role": "assistant", "content": "a"},
                    {"role": "user", "content": "typed"},
                    {"role": "user", "content": [
                        {"type": "text", "text": "[Image you opened: p.png]"},
                        img]}]
    eng._sanitize_messages()
    assert len(eng.messages) == 2
    parts = eng.messages[1]["content"]
    assert parts[0] == {"type": "text", "text": "typed"}
    assert img in parts


def test_a_harness_mark_survives_only_when_every_part_carries_it():
    eng = _bare()
    eng.messages = [{"role": "assistant", "content": "a"},
                    {"role": "user", "content": "typed"},
                    {"role": "user", "content": "note", "_harness": True}]
    eng._sanitize_messages()
    assert "_harness" not in eng.messages[1]
    eng.messages = [{"role": "assistant", "content": "a"},
                    {"role": "user", "content": "n1", "_harness": True},
                    {"role": "user", "content": "n2", "_harness": True}]
    eng._sanitize_messages()
    assert eng.messages[1].get("_harness") is True


def test_a_pinned_row_is_kept_verbatim_on_either_side():
    eng = _bare()
    eng.messages = [{"role": "assistant", "content": "a"},
                    {"role": "user", "content": "typed"},
                    {"role": "user", "content": "pinned", "_pinned": True}]
    eng._sanitize_messages()
    users = [m for m in eng.messages if m.get("role") == "user"]
    assert [u["content"] for u in users][-1] == "pinned"
    assert any(u["content"] == "typed" for u in users)
    roles = [m.get("role") for m in eng.messages]
    assert all(not (a == b) for a, b in zip(roles, roles[1:]))
