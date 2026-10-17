"""A session that only announces work, turn after turn, stops after three.

Input: a model whose every answer announces work it does not do ("Let me
resolve the next task.") and makes no tool call, while a task stays open
and no question to the person is pending. Driven through the real
``AgentEngine.stream_response``.

Before: each such turn left a follow-up note and the session was nudged
again and again (wave 16, 2026-10-16: 35 nudges in one session, 20 in
another), because an announcing answer counted as progress and reset the
stall guard's streak every turn.

Pinned here:

* three announcing turns in a row stop the follow-ups;
* two do not;
* a line the person types lifts the stop and starts a fresh count.

``tests/test_zz_h1e_*`` pins the other half: a turn whose answer is a
leaked tool call advances the guard once, under its own kind, so the two
kinds never reset each other.
"""
from __future__ import annotations

from unittest.mock import patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


class _Client:
    """A model that answers with announced-but-not-done intent and makes no
    tool call — so the turn ends with open work and no question pending."""

    calls = 0

    def stream_message(self, *_a, **_k):
        type(self).calls += 1
        yield StreamEvent(type="text_delta", text="Let me resolve the next task.")
        yield StreamEvent(type="message_delta", stop_reason="end_turn",
                          input_tokens=1, output_tokens=1)


@pytest.fixture
def make(tmp_path):
    def _make():
        client = _Client()
        with patch("delfin.agent.engine.create_client", return_value=client):
            eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api",
                                provider="kit", model="kit.qwen3.5-397b-A17b",
                                mode="solo")
        # The turn ends with work still open, and no question is pending.
        eng._open_task_count = lambda: 1
        return eng, client
    return _make


def _run_three(make):
    eng, _ = make()
    for i in range(3):
        eng.stream_response(f"turn {i}")
        # The continuation latch records one note per turn; drop it so each of
        # the three turns is judged afresh (as the harness does between
        # follow-ups).
        eng.clear_turn_continuation()
    return eng


def test_three_open_work_turns_in_a_row_stop(make):
    eng = _run_three(make)
    assert eng._stall_guard_stopped, (
        "three turns in a row that end with open work and no pending "
        "question must stop, like three identical dead turns")


def test_two_open_work_turns_do_not_stop(make):
    eng, _ = make()
    eng.stream_response("go")
    eng.clear_turn_continuation()
    eng.stream_response("harness follow-up 1")
    assert eng._stall_guard_stopped == ""


def test_a_question_is_not_a_dead_ending(make):
    """A turn that ends by ASKING the person is not the dead kind W2 guards:
    there is a pending question, so it must not count toward the stop."""
    eng = _run_three(make)                     # three open-work-dead turns
    assert eng._stall_guard_stopped            # ... stop
    eng.note_typed("I will do it")             # a person answers -> fresh count
    assert eng._stall_guard_stopped == ""
    eng._open_task_count = lambda: 1
    for i in range(2):
        eng.stream_response(f"again {i}")
        eng.clear_turn_continuation()
    assert eng._stall_guard_stopped == ""      # two, not three, after the answer
