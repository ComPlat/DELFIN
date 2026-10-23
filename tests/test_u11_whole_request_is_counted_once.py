"""What the whole-request guard in ``api_client._compact_between_rounds``
counts, and what its compaction must keep sendable.

Input: the round rows the client sends (``api_messages``), the planning
window (tokens), ``system_chars`` and ``schema_chars``. The guard sums the
request in chars (4 chars/token) and compacts older rounds above 95 % of the
window. Pinned here:

* The system prompt is row 0 of ``api_messages`` (role system/developer);
  ``system_chars`` is booked only when no such row exists, so the prompt is
  counted once. Measured on the client loop: a 50 000-char prompt arrived as
  row 0 AND as ``system_chars=50000``.
* Tool-call arguments are request text: a turn that writes files carries
  its payload there, with ``content=None``.
* Advertised tool schemas ride on every request outside the rows: 41 244
  chars for the 69 tools a client without permissions sends (~10 k tokens,
  7.9 % of a 131 072 window).
* A compaction never splits an assistant message with parallel tool calls
  from any of its results: an orphaned ``tool`` row is rejected by the API.
"""

from __future__ import annotations

import json
import threading
import types

import pytest

from delfin.agent import api_client
from delfin.agent import in_turn_compaction

_WINDOW = 131_072
_BUDGET = int(_WINDOW * 4 * 0.95)          # 498 073 chars


@pytest.fixture(autouse=True)
def _no_log(monkeypatch):
    monkeypatch.setattr("delfin.agent.compaction_log.record_compaction",
                        lambda *a, **k: None)


def _round(i, *, tool_chars=100, args_chars=2, content="r"):
    return [
        {"role": "assistant", "content": content,
         "tool_calls": [{"id": f"c{i}", "type": "function",
                         "function": {"name": "write_file",
                                      "arguments": "x" * args_chars}}]},
        {"role": "tool", "tool_call_id": f"c{i}", "content": "x" * tool_chars},
    ]


def _paired(msgs):
    calls = {tc["id"] for m in msgs for tc in (m.get("tool_calls") or [])}
    answers = {m["tool_call_id"] for m in msgs if m.get("role") == "tool"}
    return calls == answers


def test_the_system_row_is_not_counted_a_second_time():
    # 260 000 once is under the 498 073-char line; twice it is over.
    msgs = [{"role": "system", "content": "s" * 260_000},
            {"role": "user", "content": "goal"}]
    for i in range(4):
        msgs += _round(i)
    before = [dict(m) for m in msgs]
    note = api_client._compact_between_rounds(
        msgs, 10**9, 0, "s", planning_window=_WINDOW, system_chars=260_000)
    assert note == ""
    assert msgs == before


def test_system_chars_count_when_the_rows_carry_no_system_row():
    msgs = [{"role": "user", "content": "goal"}]
    for i in range(4):
        msgs += _round(i, tool_chars=2_000)
    note = api_client._compact_between_rounds(
        msgs, 10**9, 0, "s", planning_window=_WINDOW, system_chars=510_000)
    assert "summarised" in note


def test_tool_call_arguments_are_counted():
    # 4 x 150 000 chars of write_file payload, content None, tiny results.
    msgs = [{"role": "user", "content": "goal"}]
    for i in range(4):
        msgs += _round(i, args_chars=150_000, content=None)
    note = api_client._compact_between_rounds(
        msgs, 10**9, 0, "s", planning_window=_WINDOW)
    assert "summarised" in note
    assert _paired(msgs)


def test_tool_schemas_are_counted():
    msgs = [{"role": "user", "content": "goal"}]
    for i in range(4):
        msgs += _round(i, tool_chars=115_000)    # 460 000 chars of results
    note = api_client._compact_between_rounds(
        msgs, 10**9, 0, "s", planning_window=_WINDOW, schema_chars=41_244)
    assert "summarised" in note


def test_parallel_calls_keep_every_result_paired():
    msgs = [{"role": "user", "content": "goal"}]
    for i in range(3):
        msgs += _round(i, content="x" * 170_000)
    msgs.append({"role": "assistant", "content": None, "tool_calls": [
        {"id": "p1", "type": "function",
         "function": {"name": "read_file", "arguments": "{}"}},
        {"id": "p2", "type": "function",
         "function": {"name": "read_file", "arguments": "{}"}}]})
    msgs.append({"role": "tool", "tool_call_id": "p1", "content": "a" * 10})
    msgs.append({"role": "tool", "tool_call_id": "p2", "content": "b" * 10})
    note = api_client._compact_between_rounds(
        msgs, 10**9, 0, "s", planning_window=_WINDOW)
    assert "summarised" in note
    assert _paired(msgs)
    assert [m.get("tool_call_id") for m in msgs[-2:]] == ["p1", "p2"]


def test_a_parallel_group_is_compacted_whole_or_kept_whole():
    msgs = [{"role": "user", "content": "goal"},
            {"role": "assistant", "content": None, "tool_calls": [
                {"id": "a", "type": "function",
                 "function": {"name": "bash", "arguments": "{}"}},
                {"id": "b", "type": "function",
                 "function": {"name": "bash", "arguments": "{}"}}]},
            {"role": "tool", "tool_call_id": "a", "content": "x" * 5_000},
            {"role": "tool", "tool_call_id": "b", "content": "x" * 5_000}]
    for i in range(3):
        msgs += _round(i, tool_chars=100)
    new, rec = in_turn_compaction.compact(msgs, 6_000)
    assert rec is not None
    assert _paired(new)


def _client(system):
    tc = types.SimpleNamespace(index=0, id="c1", function=types.SimpleNamespace(
        name="read_file", arguments='{"path": "x"}'))
    calls = types.SimpleNamespace(usage=None, choices=[types.SimpleNamespace(
        delta=types.SimpleNamespace(content=None, tool_calls=[tc],
                                    reasoning_content=None),
        finish_reason="tool_calls")])

    class _Stream:
        def __iter__(self):
            return iter((calls,))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    _Stub.sent.append(kw)
                    return _Stream()
    _Stub.sent = []
    c = api_client.OpenAIClient.__new__(api_client.OpenAIClient)
    c.client = _Stub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://example.invalid/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    c.metrics_session = "s"
    c._system_prompt_chars = len(system)
    return c, _Stub


def test_the_client_loop_books_schemas_and_the_prompt_once(monkeypatch):
    seen = []
    real = api_client._compact_between_rounds

    def spy(msgs, budget, base, sid, **kw):
        seen.append((msgs[0]["role"], kw))
        return real(msgs, budget, base, sid, **kw)

    monkeypatch.setattr(api_client, "_compact_between_rounds", spy)
    system = "s" * 50_000
    client, stub = _client(system)
    try:
        for _ev in client.stream_message(
                messages=[{"role": "user", "content": "hi"}],
                system=system, max_tokens=64):
            if len(seen) >= 2:
                break
    except Exception:
        pass
    assert seen, "the round loop never reached the guard"
    role, kw = seen[0]
    assert role == "system"
    assert int(kw.get("system_chars", 0) or 0) == 0
    tools = stub.sent[0].get("tools") if stub.sent else None
    assert tools
    assert kw.get("schema_chars") == len(json.dumps(tools))
