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
    with WorkerClient(tmp_path, delay=0.02) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("turn.started", after=start)
        got = w.request("session.prompt", {"text": "and then this"})
        assert got == {"accepted": True, "queued": True}
        queued, _ = w.wait_for(lambda e: e["type"] == "user.message"
                               and e["source"] == "queued", after=start)
        assert queued["text"] == "and then this"
        _, first_end = w.wait_for("turn.ended", after=start)
        _, second_start = w.wait_for("turn.started", after=first_end)
        _, second_end = w.wait_for("turn.ended", after=second_start)
        # The queued message is shown once, when it was queued -- not again.
        shown = [e for e in w.events[start:second_end]
                 if e["type"] == "user.message" and e["text"] == "and then this"]
        assert len(shown) == 1


def test_a_steered_message_reaches_the_running_turn(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "steer test"})
        w.wait_for("turn.started", after=start)
        time.sleep(0.2)
        got = w.request("session.prompt", {"text": "use the SMD model",
                                           "how": "steer"})
        assert got["steered"] is True
        _, end = w.wait_for("turn.ended", after=start)
        text = "".join(e["text"] for e in w.events[start:end] if e["type"] == "text.ended")
        assert "Noted mid-turn: use the SMD model" in text
        assert any(e["type"] == "user.message" and e["source"] == "steer"
                   for e in w.events[start:end])


def test_cancel_ends_the_turn_early(tmp_path):
    with WorkerClient(tmp_path, delay=0.05) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("text.started", after=start)
        began = time.monotonic()
        w.request("session.cancel")
        turn_end, _ = w.wait_for("turn.ended", after=start)
        assert time.monotonic() - began < 3.0
        full = len("Working through a long calculation step by step. " * 12)
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
    assert w.proc.wait(timeout=20) == 0
    w.wait_for("session.ended")
