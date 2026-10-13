"""A TUI session reports its turn as events a window can draw as it happens.

Driven through the real worker process with the scripted engine: the order
of a turn's events, tool results landing on their own call, the answer's
deltas adding up to the text that is kept, and a stdout that carries the
protocol and nothing else.
"""

from __future__ import annotations

import json

from delfin.agent.service.protocol import DURABLE_EVENTS
from tui_worker_client import WorkerClient


def test_the_worker_introduces_the_session_before_anything_else(tmp_path):
    with WorkerClient(tmp_path) as w:
        info = w.init_result["session"]
        assert w.init_result["protocolVersion"] == 1
        assert info["provider"] == "fake" and info["workspace"] == str(w.workspace)
        assert w.types()[:3] == ["worker.ready", "session.info", "session.state"]
        assert w.events[2]["state"] == "idle"


def test_a_turn_arrives_in_order_and_ends_idle(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("hello there")
        types = w.types(start, end + 1)
        # Shown the moment it arrives; busy once the loop takes it.
        assert types[:3] == ["user.message", "session.state", "turn.started"]
        assert w.events[start + 1]["state"] == "running"
        assert types[-1] == "turn.ended"
        order = ["reasoning.started", "reasoning.ended", "tool.started",
                 "tool.ended", "text.started", "text.ended"]
        assert [t for t in types if t in order] == order
        w.wait_for(lambda e: e["type"] == "session.state" and e["state"] == "idle",
                   after=end)
        turn_end = w.events[end]
        assert turn_end["toolCalls"] == 1
        assert turn_end["inputTokens"] > 0 and turn_end.get("error") is None


def test_a_tool_result_lands_on_the_call_it_belongs_to(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("hello")
        started = [e for e in w.events[start:end] if e["type"] == "tool.started"]
        ended = [e for e in w.events[start:end] if e["type"] == "tool.ended"]
        assert [e["callId"] for e in started] == [e["callId"] for e in ended]
        assert started[0]["input"] == {"path": "README.md"}
        assert ended[0]["meta"]["ok"] is True


def test_the_kept_text_is_exactly_what_was_streamed(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("hello")
        part = next(e for e in w.events[start:end] if e["type"] == "text.ended")
        deltas = [e["delta"] for e in w.events[start:end]
                  if e["type"] == "text.delta" and e["partId"] == part["partId"]]
        assert deltas and "".join(deltas) == part["text"]
        assert "**one** README" in part["text"]


def test_durable_events_are_numbered_in_the_order_they_were_sent(tmp_path):
    with WorkerClient(tmp_path) as w:
        w.turn("hello")
        w.turn("hello again")
        durable = [e for e in w.events if e["type"] in DURABLE_EVENTS]
        seqs = [e["seq"] for e in durable]
        assert seqs == list(range(1, len(seqs) + 1))
        # Live-only events are never replayed, so they carry no number.
        assert not any("seq" in e for e in w.events
                       if e["type"] not in DURABLE_EVENTS)
        session = w.init_result["session"]["sessionId"]
        assert {e["session"] for e in w.events} == {session}


def test_a_failing_turn_says_why_and_the_session_carries_on(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("fail now")
        assert "scripted failure" in w.events[end]["error"]
        assert any(e["type"] == "error" for e in w.events[start:end])
        _, end2 = w.turn("hello")
        assert w.events[end2].get("error") is None


def test_stdout_carries_the_protocol_and_nothing_else(tmp_path):
    with WorkerClient(tmp_path) as w:
        w.turn("be noisy")
        w.turn("hello")
    for line in w.raw_lines:
        json.loads(line)              # every byte on stdout is protocol
    err = w.stderr()
    assert "stray print from inside a turn" in err
    assert "stray output from a child process" in err


def test_shutdown_ends_the_session_with_a_last_event(tmp_path):
    w = WorkerClient(tmp_path)
    w.turn("hello")
    assert w.close() == 0
    assert w.events[-1]["type"] == "session.ended"
    assert w.events[-1]["exitCode"] == 0


def test_a_commands_output_keeps_its_lines_in_one_notice(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "/help"})
        notice, at = w.wait_for("notice", after=start)
        assert "\n" in notice["text"]
        # Nothing of /help is split into a second notice.
        w.wait_for(lambda e: e["type"] == "session.state"
                   and e["state"] == "idle", after=at)
        assert [e["type"] for e in w.events[start:]].count("notice") == 1
        assert not any(e["type"] == "turn.started" for e in w.events[start:])
