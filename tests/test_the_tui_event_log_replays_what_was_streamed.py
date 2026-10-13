"""The session's event log gives a window that attaches late what it missed.

Durable events are written in the order they were sent, with the numbers
they were sent with; a replay of the last turns starts at a user message
and carries the session's header; the file is the user's alone.
"""

from __future__ import annotations

import json
import stat

from delfin.agent.service import eventlog
from delfin.agent.service.protocol import DURABLE_EVENTS
from tui_worker_client import WorkerClient


def test_numbers_continue_across_a_restart_and_only_durable_events_are_kept(tmp_path):
    path = tmp_path / "s.jsonl"
    log = eventlog.EventLog(path)
    for kind in ("session.info", "text.delta", "user.message", "usage", "turn.ended"):
        log.stamp({"type": kind})
    kept = eventlog.read_events(path)
    assert [e["type"] for e in kept] == ["session.info", "user.message", "turn.ended"]
    assert [e["seq"] for e in kept] == [1, 3, 5]
    again = eventlog.EventLog(path)
    assert again.stamp({"type": "notice"})["seq"] == 6


def test_a_line_still_being_written_is_not_read(tmp_path):
    path = tmp_path / "s.jsonl"
    path.write_text('{"type": "notice", "seq": 1}\n{"type": "notice", "se')
    assert eventlog.read_events(path) == [{"type": "notice", "seq": 1}]


def test_a_replay_starts_at_a_user_message_with_the_header_first():
    events = [{"type": "session.info", "seq": 1}]
    for n in range(1, 4):
        events += [{"type": "user.message", "seq": 10 * n},
                   {"type": "turn.ended", "seq": 10 * n + 1}]
    picked = eventlog.tail_turns(events, 2)
    assert [e["seq"] for e in picked] == [1, 20, 21, 30, 31]
    assert eventlog.tail_turns(events, 0) == events
    assert eventlog.tail_turns(events, 99) == events


def test_the_replay_is_what_a_window_watching_from_the_start_kept(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        w.turn("hello")
        start = len(w.events)
        w.request("session.prompt", {"text": "edit it"})
        asked, _ = w.wait_for("permission.asked", after=start)
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        w.wait_for("turn.ended", after=start)
        w.turn("hello again")
        replay = w.request("session.replay", {"turns": 2})["events"]
    live_durable = [e for e in w.events if e["type"] in DURABLE_EVENTS]
    assert replay[0]["type"] == "session.info"
    body = replay[1:]
    assert body[0]["type"] == "user.message" and body[0]["text"] == "edit it"
    by_seq = {e["seq"]: e for e in live_durable}
    for event in body:
        assert by_seq[event["seq"]] == event
    assert [e["seq"] for e in body] == \
        [e["seq"] for e in live_durable if e["seq"] >= body[0]["seq"]
         and e["type"] != "session.ended"][:len(body)]


def test_the_log_is_readable_by_its_owner_only(tmp_path):
    with WorkerClient(tmp_path) as w:
        w.turn("hello")
        session = w.init_result["session"]["sessionId"]
    path = w.home / ".delfin" / "agent_service" / "sessions" / f"{session}.jsonl"
    assert path.exists(), list((w.home / ".delfin").rglob("*.jsonl"))
    assert stat.S_IMODE(path.stat().st_mode) & 0o077 == 0
    assert stat.S_IMODE(path.parent.stat().st_mode) & 0o077 == 0
    for line in path.read_text().splitlines():
        assert json.loads(line)["type"] in DURABLE_EVENTS
