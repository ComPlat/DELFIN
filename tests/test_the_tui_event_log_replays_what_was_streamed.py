"""The session's event log gives a window that attaches late what it missed.

Durable events are written in the order they were sent, with the numbers
they were sent with, and no number is ever given out twice -- not across a
restart either; a replay of the last turns starts where a turn started and
carries the session's header as it is now; the file is the user's alone.
"""

from __future__ import annotations

import json
import stat

from delfin.agent.service import eventlog
from delfin.agent.service.protocol import DURABLE_EVENTS
from tui_worker_client import WorkerClient


def test_only_durable_events_are_numbered_and_kept(tmp_path):
    path = tmp_path / "s.jsonl"
    log = eventlog.EventLog(path)
    stamped = [log.stamp({"type": kind}) for kind in
               ("session.info", "text.delta", "user.message", "usage",
                "turn.ended")]
    assert [e.get("seq") for e in stamped] == [1, None, 2, None, 3]
    kept = eventlog.read_events(path)
    assert [e["type"] for e in kept] == ["session.info", "user.message", "turn.ended"]
    assert [e["seq"] for e in kept] == [1, 2, 3]


def test_a_number_is_never_given_out_twice_across_a_restart(tmp_path):
    path = tmp_path / "s.jsonl"
    log = eventlog.EventLog(path)
    for kind in ("notice", "text.delta", "text.delta", "notice"):
        log.stamp({"type": kind})
    log.close()
    again = eventlog.EventLog(path)
    assert again.stamp({"type": "notice"})["seq"] == 3


def test_a_line_still_being_written_is_not_read(tmp_path):
    path = tmp_path / "s.jsonl"
    path.write_text('{"type": "notice", "seq": 1}\n{"type": "notice", "se')
    assert eventlog.read_events(path) == [{"type": "notice", "seq": 1}]


def test_a_torn_last_line_does_not_swallow_the_next_runs_first_event(tmp_path):
    path = tmp_path / "s.jsonl"
    path.write_text('{"type": "notice", "seq": 1}\n{"type": "notice", "se')
    log = eventlog.EventLog(path)
    log.stamp({"type": "session.info"})
    log.close()
    kept = eventlog.read_events(path)
    assert [(e["type"], e["seq"]) for e in kept] == \
        [("notice", 1), ("session.info", 2)]


def test_text_that_is_not_utf8_is_kept_not_fatal(tmp_path):
    import os
    path = tmp_path / "s.jsonl"
    log = eventlog.EventLog(path)
    odd = os.fsdecode(b"caf\xe9.xyz")
    log.stamp({"type": "notice", "text": odd})
    assert log.writable
    assert eventlog.read_events(path)[0]["text"] == odd


def _turn(n: int, *before: dict) -> list[dict]:
    return [*before, {"type": "user.message", "seq": 10 * n, "source": "user"},
            {"type": "turn.started", "seq": 10 * n + 1},
            {"type": "turn.ended", "seq": 10 * n + 5}]


def test_a_replay_starts_where_a_turn_started_with_the_header_first():
    events = [{"type": "session.info", "seq": 1}, *_turn(1), *_turn(2), *_turn(3)]
    picked = eventlog.tail_turns(events, 2)
    assert [e["seq"] for e in picked] == [1, 20, 21, 25, 30, 31, 35]
    assert eventlog.tail_turns(events, 0) == events
    assert eventlog.tail_turns(events, 99) == events


def test_a_message_queued_during_a_turn_is_not_a_turn_of_its_own():
    events = [{"type": "session.info", "seq": 1}, *_turn(1)]
    # Turn 2 runs; a message is queued during it, and becomes turn 3.
    events += [{"type": "user.message", "seq": 20, "source": "user"},
               {"type": "turn.started", "seq": 21},
               {"type": "user.message", "seq": 22, "source": "queued"},
               {"type": "turn.ended", "seq": 25},
               {"type": "turn.started", "seq": 31},
               {"type": "turn.ended", "seq": 35}]
    assert [e["seq"] for e in eventlog.tail_turns(events, 1)] == [1, 31, 35]
    assert [e["seq"] for e in eventlog.tail_turns(events, 2)] == \
        [1, 20, 21, 22, 25, 31, 35]


def test_the_header_given_replaces_the_logged_ones():
    events = [{"type": "session.info", "seq": 1, "approval": "plan"}, *_turn(1),
              {"type": "session.info", "seq": 15, "approval": "plan"}, *_turn(2)]
    now = {"type": "session.info", "approval": "acceptEdits"}
    picked = eventlog.tail_turns(events, 1, header=now)
    assert picked[0] is now
    assert [e.get("seq") for e in picked] == [None, 20, 21, 25]


def test_the_replay_is_what_a_window_watching_from_the_start_kept(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        w.turn("hello")
        start = len(w.events)
        w.request("session.prompt", {"text": "edit it"})
        asked, _ = w.wait_for("permission.asked", after=start)
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        w.wait_for("turn.ended", after=start)
        w.turn("hello again")
        w.request("session.setMode", {"mode": "acceptEdits"})
        replay = w.request("session.replay", {"turns": 2})["events"]
    live_durable = [e for e in w.events if e["type"] in DURABLE_EVENTS]
    header = replay[0]
    assert header["type"] == "session.info"
    assert header["approval"] == "acceptEdits"        # as it is now
    body = replay[1:]
    assert body[0]["type"] == "user.message" and body[0]["text"] == "edit it"
    by_seq = {e["seq"]: e for e in live_durable}
    for event in body:
        assert by_seq[event["seq"]] == event
    assert [e["seq"] for e in body] == \
        [e["seq"] for e in live_durable if e["seq"] >= body[0]["seq"]
         and e["type"] not in ("session.ended", "session.info")][:len(body)]


def test_the_log_is_readable_by_its_owner_only(tmp_path):
    with WorkerClient(tmp_path) as w:
        w.turn("hello")
    path = w.log_path()
    assert path.exists(), list(tmp_path.rglob("*.jsonl"))
    assert not list(w.home.rglob("*.jsonl"))       # nothing in the real home
    assert stat.S_IMODE(path.stat().st_mode) & 0o077 == 0
    for parent in path.parents[:3]:
        assert stat.S_IMODE(parent.stat().st_mode) & 0o077 == 0, parent
    for line in path.read_text().splitlines():
        assert json.loads(line)["type"] in DURABLE_EVENTS


def test_a_resumed_session_numbers_its_turns_on(tmp_path, monkeypatch):
    import io

    from delfin.agent.service import protocol
    from delfin.agent.service.fake_engine import FakeEngine
    from delfin.agent.service.worker import Worker

    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    engine = FakeEngine(str(tmp_path))
    log = eventlog.EventLog(eventlog.log_path(engine.session_id))
    for n in (1, 2):
        log.stamp({"type": "turn.started", "turnId": f"t{n}"})
    log.close()
    worker = Worker(engine, writer=protocol.LineWriter(io.BytesIO()),
                    reader=io.BytesIO(), workspace=tmp_path, banner=[],
                    notices=[], fake=True, save=False)
    assert worker.agent._turn_seq == 2      # the next turn is t3
