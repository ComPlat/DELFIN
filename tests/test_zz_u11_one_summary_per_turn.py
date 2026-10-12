"""One in-turn summary row per request, however long the turn runs.

Input: ``OpenAIClient.stream_message`` against a stub endpoint that calls
``read_file`` on a new 80 kB file every round, for 32 rounds, through both
compaction routes (the tool-char budget, and the whole-request guard of
#158). Output measured per request: the number of rows that start with
the in-turn summary header, and the chars of those rows.

Measured on main 6d6295af: the summaries accumulate, one per compaction.
The last request of the tool route held 30 summary rows (8510 chars, +283
chars per round), the whole route 15 (6965 chars, +465 per compaction), so
the compaction grows the request it exists to shrink and, over hundreds of
rounds, the summaries alone outgrow the window.
"""

from __future__ import annotations

import copy
import json
import types

import pytest

from delfin.agent import api_client
from delfin.agent import in_turn_compaction as itc

HEADER = itc._SUMMARY_HEADER


def _drive(tmp_path, monkeypatch, route, rounds=32):
    """Every request the loop sent, as a deep copy of its rows."""
    ws = tmp_path / "ws"
    ws.mkdir()
    for k in range(rounds):
        (ws / f"f{k}.txt").write_text((f"line {k} of text " * 8 + "\n") * 600)
    seen = []

    class _Stream:
        def __init__(self, i):
            self.i = i

        def __iter__(self):
            if self.i >= rounds:
                delta = types.SimpleNamespace(content="done", tool_calls=None,
                                              reasoning_content=None)
                fin = "stop"
            else:
                delta = types.SimpleNamespace(
                    content="y" * (60_000 if route == "whole" else 10),
                    reasoning_content=None,
                    tool_calls=[types.SimpleNamespace(
                        index=0, id=f"c{self.i}",
                        function=types.SimpleNamespace(
                            name="read_file", arguments=json.dumps(
                                {"path": str(ws / f"f{self.i}.txt")})))])
                fin = "tool_calls"
            return iter((types.SimpleNamespace(usage=None, choices=[
                types.SimpleNamespace(delta=delta, finish_reason=fin)]),))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    seen.append(copy.deepcopy(kw.get("messages")))
                    return _Stream(len(seen) - 1)

    if route == "tool":
        monkeypatch.setattr(api_client, "_tool_context_char_budget",
                            lambda caps, model="": 12_000)
    else:
        monkeypatch.setattr(api_client, "_tool_context_char_budget",
                            lambda caps, model="": 10_000_000)
        monkeypatch.setattr(api_client, "_planning_window_for_turn",
                            lambda self, caps: 60_000)
    client = api_client.OpenAIClient(
        api_key="x", model="kit.glm-5.3", provider="kit",
        base_url="https://ki-toolbox.scc.kit.edu/api/v1",
        permissions=api_client.KitToolPermissions(
            workspace=ws, mode="bypassPermissions", task_session_id="S"))
    client.client = _Stub()
    for _ev in client.stream_message(
            messages=[{"role": "user", "content": "go"}], system="t",
            max_tokens=64):
        pass
    return seen


def _summaries(rows):
    return [str(m.get("content") or "") for m in rows
            if m.get("role") == "user"
            and str(m.get("content") or "").startswith(HEADER)]


def _orphans(rows):
    """Tool rows whose call is not in an earlier assistant row, and calls
    with no tool row after them."""
    calls, answered, bad = set(), set(), 0
    for m in rows:
        for tc in m.get("tool_calls") or []:
            calls.add(str(tc.get("id")))
        if m.get("role") == "tool":
            cid = str(m.get("tool_call_id"))
            if cid not in calls:
                bad += 1
            answered.add(cid)
    return bad + len(calls - answered)


@pytest.mark.parametrize("route", ["tool", "whole"])
def test_a_long_turn_carries_one_summary_row(tmp_path, monkeypatch, route):
    seen = _drive(tmp_path, monkeypatch, route)
    assert len(seen) == 33
    counts = [len(_summaries(rows)) for rows in seen]
    assert max(counts) >= 1, "the turn never compacted"
    assert max(counts) == 1, f"summary rows per request: {counts}"
    # The fold is bounded: the one summary never outgrows the cap.
    sizes = [sum(len(s) for s in _summaries(rows)) for rows in seen]
    assert max(sizes) <= itc._SUMMARY_MAX_CHARS, sizes
    # The newest compacted round is named in the one summary: the fold
    # keeps the newest information when it has to drop.
    kept = {str(m.get("tool_call_id")) for m in seen[-1]
            if m.get("role") == "tool"}
    newest = max((i for i in range(32) if f"c{i}" not in kept), default=None)
    assert newest is not None
    assert f"round result c{newest}:" in _summaries(seen[-1])[0]
    for rows in seen:
        assert _orphans(rows) == 0
