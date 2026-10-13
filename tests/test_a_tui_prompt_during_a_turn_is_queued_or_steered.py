"""What a window sends while a turn runs: queued, steered, or a cancel.

The REPL's rules, kept: Enter during a turn queues the message for when the
turn ends; steering sends it into the turn the model is running; Esc ends
the turn. The approval posture changes as shift+tab changes it -- never
onto bypass.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent.service import protocol
from tui_worker_client import WorkerClient, WorkerError


def test_a_message_sent_during_a_turn_waits_for_the_next(tmp_path):
    # Sent the moment the turn is announced, again and again: the window
    # between "turn.started" and the engine's own turn flag is where a
    # message used to become a turn of its own.
    with WorkerClient(tmp_path, delay=0.0) as w:
        for n in range(6):
            start = len(w.events)
            w.request("session.prompt", {"text": f"brief please ({n})"})
            w.wait_for("turn.started", after=start)
            got = w.request("session.prompt", {"text": f"and then this ({n})"})
            assert got == {"accepted": True, "queued": True}
            queued, _ = w.wait_for(lambda e: e["type"] == "user.message"
                                   and e["source"] == "queued", after=start)
            assert queued["text"] == f"and then this ({n})"
            _, first_end = w.wait_for("turn.ended", after=start)
            _, second_start = w.wait_for("turn.started", after=first_end)
            _, second_end = w.wait_for("turn.ended", after=second_start)
            # Shown once, when it was queued -- not again when it runs.
            shown = [e for e in w.events[start:second_end]
                     if e["type"] == "user.message"
                     and e["text"] == f"and then this ({n})"]
            assert len(shown) == 1
            w.wait_for(lambda e: e["type"] == "session.state"
                       and e["state"] == "idle", after=second_end)


def test_a_steered_message_reaches_the_running_turn(tmp_path):
    with WorkerClient(tmp_path) as w:
        for n in range(5):
            start = len(w.events)
            w.request("session.prompt", {"text": "steer test"})
            w.wait_for("turn.started", after=start)
            got = w.request("session.prompt", {"text": f"use the SMD model {n}",
                                               "how": "steer"})
            assert got["steered"] is True
            _, end = w.wait_for("turn.ended", after=start)
            assert f"Noted mid-turn: use the SMD model {n}" in w.text_of(start, end)
            assert any(e["type"] == "user.message" and e["source"] == "steer"
                       for e in w.events[start:end])


def test_steering_with_no_turn_running_is_an_ordinary_message(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        got = w.request("session.prompt", {"text": "hello", "how": "steer"})
        assert got == {"accepted": True}
        w.wait_for("turn.ended", after=start)


def test_a_cancel_right_after_the_prompt_still_ends_that_turn(tmp_path):
    # The engine clears its stop flag when its turn begins; a cancel that
    # lands before that must hold anyway.
    full = len("Working through a long calculation step by step. " * 4)
    with WorkerClient(tmp_path, delay=0.05) as w:
        for n in range(5):
            start = len(w.events)
            w.send("session.prompt", {"text": f"slow please ({n})"})
            got = w.request("session.cancel")
            assert got["stopping"] is True
            _, end = w.wait_for("turn.ended", after=start)
            said = "".join(e["delta"] for e in w.events[start:end]
                           if e["type"] == "text.delta")
            assert len(said) < full
            w.wait_for(lambda e: e["type"] == "session.state"
                       and e["state"] == "idle", after=start)
        _, end = w.turn("hello")          # the next turn is not cancelled
        assert w.events[end].get("error") is None


def test_a_cancel_with_nothing_running_stops_nothing(tmp_path):
    with WorkerClient(tmp_path) as w:
        assert w.request("session.cancel") == {"stopping": False, "refused": 0}
        _, end = w.turn("hello")
        assert w.events[end].get("error") is None


def test_cancel_ends_the_turn_early(tmp_path):
    with WorkerClient(tmp_path, delay=0.05) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("text.started", after=start)
        w.request("session.cancel")
        turn_end, _ = w.wait_for("turn.ended", after=start)
        full = len("Working through a long calculation step by step. " * 4)
        text = "".join(e["delta"] for e in w.events[start:] if e["type"] == "text.delta")
        assert len(text) < full
        assert turn_end["turnId"]


def test_the_posture_moves_along_the_ladder_but_never_onto_bypass(tmp_path):
    with WorkerClient(tmp_path) as w:
        assert w.init_result["session"]["approval"] == "plan"
        start = len(w.events)
        assert w.request("session.setMode", {"mode": "acceptEdits"}) == \
            {"approval": "acceptEdits"}
        updated, _ = w.wait_for("session.updated", after=start)
        assert updated["approval"] == "acceptEdits"
        with pytest.raises(WorkerError) as refused:
            w.request("session.setMode", {"mode": "bypassPermissions"})
        assert refused.value.code == protocol.INVALID_PARAMS
        assert w.request("session.status")["approval"] == "acceptEdits"


def test_nothing_is_accepted_before_initialize(tmp_path):
    w = WorkerClient(tmp_path, initialize=False)
    try:
        with pytest.raises(WorkerError) as early:
            w.request("session.prompt", {"text": "hi"})
        assert early.value.code == protocol.NOT_ALLOWED
        with pytest.raises(WorkerError) as old:
            w.request("initialize", {"protocolVersion": 0})
        assert old.value.code == protocol.VERSION_MISMATCH
        w.request("initialize", {"protocolVersion": 1})
        w.turn("hello")
    finally:
        w.close()


def test_a_prompt_that_looks_like_the_start_of_a_block_is_still_one_message(tmp_path):
    with WorkerClient(tmp_path) as w:
        for text in ("a path ends in a backslash C:\\", '"""',
                     "first line\nsecond line\\\nthird"):
            start, end = w.turn(text)
            sent = next(e for e in w.events[start:end] if e["type"] == "user.message")
            assert sent["text"] == text.strip()
            assert w.events[end].get("error") is None


def test_the_window_going_away_ends_the_session_cleanly(tmp_path):
    w = WorkerClient(tmp_path)
    w.turn("hello")
    w.proc.stdin.close()            # the window side hangs up: EOF
    assert w.ended() == 129         # a hangup, as for a closed terminal
    assert w.events[-1]["type"] == "session.ended"
    assert w.events[-1]["reason"] == "eof"


def test_the_window_going_away_mid_turn_ends_it_within_the_deadline(tmp_path):
    w = WorkerClient(tmp_path, delay=0.05)
    start = len(w.events)
    w.request("session.prompt", {"text": "slow please"})
    w.wait_for("text.started", after=start)
    began = time.monotonic()
    w.proc.stdin.close()
    assert w.ended() == 129
    assert time.monotonic() - began < 15
    assert w.events[-1]["type"] == "session.ended"


def test_shutdown_drops_what_was_queued_and_runs_nothing_more(tmp_path):
    marker = tmp_path / "workspace" / "ran.txt"
    w = WorkerClient(tmp_path, delay=0.05)
    start = len(w.events)
    w.request("session.prompt", {"text": "slow please"})
    w.wait_for("turn.started", after=start)
    w.request("session.prompt", {"text": "hello after"})
    w.request("session.prompt", {"text": f"!touch {marker}"})
    assert w.close() == 0
    assert [e["type"] for e in w.events[start:]].count("turn.started") == 1
    said = "".join(e["delta"] for e in w.events[start:]
                   if e["type"] == "text.delta")
    assert len(said) < len("Working through a long calculation step by step. " * 4)
    assert not marker.exists()
    assert w.events[-1]["type"] == "session.ended"


def test_a_prompt_after_shutdown_is_refused(tmp_path):
    w = WorkerClient(tmp_path, delay=0.05)
    start = len(w.events)
    w.request("session.prompt", {"text": "slow please"})
    w.wait_for("turn.started", after=start)
    w.request("session.shutdown")
    with pytest.raises(WorkerError) as late:
        w.request("session.prompt", {"text": "one more"})
    assert late.value.code == protocol.NOT_ALLOWED
    assert w.ended() == 0


def _outage(w: WorkerClient) -> int:
    """A turn that died on an outage; returns the index of its notice."""
    start = len(w.events)
    w.request("session.prompt", {"text": "outage please"})
    _, end = w.wait_for("turn.ended", after=start)
    assert "temporarily unavailable" in w.events[end]["error"]
    _, at = w.wait_for(lambda e: e["type"] == "notice"
                       and "retrying in" in e["text"], after=end)
    return at


def test_a_message_ends_the_outage_wait_and_goes_out_at_once(tmp_path):
    with WorkerClient(tmp_path) as w:
        at = _outage(w)
        # Waiting out the outage is the session working on its own.
        assert w.states(at)[-1:] in ([], ["running"])
        assert w.request("session.status")["state"] == "running"
        began = time.monotonic()
        got = w.request("session.prompt", {"text": "hello instead"})
        assert got == {"accepted": True, "queued": True}
        started, _ = w.wait_for("turn.started", after=at)
        assert time.monotonic() - began < 20        # not the 60 s pause
        _, end = w.wait_for("turn.ended", after=at)
        assert w.events[end].get("error") is None
        # The continuation the outage would have sent is gone.
        assert not any(e["type"] == "user.message" and e["source"] == "harness"
                       for e in w.events[at:end])


def test_a_cancel_ends_the_outage_wait(tmp_path):
    with WorkerClient(tmp_path) as w:
        at = _outage(w)
        assert w.request("session.cancel")["stopping"] is True
        w.wait_for(lambda e: e["type"] == "session.state"
                   and e["state"] == "idle", after=at, timeout=20)
        assert not any(e["type"] == "turn.started" for e in w.events[at:])


def test_shutdown_during_the_outage_wait_is_prompt(tmp_path):
    w = WorkerClient(tmp_path)
    _outage(w)
    began = time.monotonic()
    assert w.close() == 0
    assert time.monotonic() - began < 20             # not the 60 s pause


def _in_process_worker(tmp_path, monkeypatch):
    import io

    from delfin.agent.service.fake_engine import FakeEngine
    from delfin.agent.service.worker import Worker

    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    out = io.BytesIO()
    worker = Worker(FakeEngine(str(tmp_path)), writer=protocol.LineWriter(out),
                    reader=io.BytesIO(), workspace=tmp_path, banner=[],
                    notices=[], fake=True, save=False)
    return worker, out


def test_a_line_queued_on_the_way_to_the_prompt_is_read_there(tmp_path, monkeypatch):
    # The loop looks for queued lines, then goes to the prompt; a line
    # queued in between must be read at the prompt, and not shown twice.
    import json

    worker, out = _in_process_worker(tmp_path, monkeypatch)
    agent = worker.agent
    agent.idle = False                    # between a turn and the prompt
    agent.queued.append("late line")
    assert agent._line_ready() is True
    assert agent._read_next() == "late line"
    assert agent.queued == []
    agent._show_user_input("late line")   # shown already, when it was queued
    agent._show_user_input("next line")   # an ordinary one is shown
    shown = [json.loads(line)["params"] for line in out.getvalue().splitlines()]
    assert [e["text"] for e in shown if e["type"] == "user.message"] == ["next line"]


def test_a_message_sent_as_a_turn_ends_is_not_left_behind(tmp_path):
    # Each prompt goes out the moment the last turn says it ended: while
    # the loop is between that turn and the prompt, where a queued line
    # could miss both the look for queued lines and the idle read.
    with WorkerClient(tmp_path, delay=0.0) as w:
        first = len(w.events)
        for n in range(20):
            w.turn(f"hello {n}")
        shown = [e["text"] for e in w.events[first:] if e["type"] == "user.message"]
        assert shown == [f"hello {n}" for n in range(20)]


def test_exit_typed_as_a_prompt_ends_the_session(tmp_path):
    w = WorkerClient(tmp_path)
    w.turn("hello")
    w.request("session.prompt", {"text": "/exit"})
    assert w.ended() == 0
    assert w.events[-1]["type"] == "session.ended"
    assert w.events[-1]["reason"] == "exit"
