"""What a broken turn must not destroy, and what Stop must actually stop.

Three defects at the seams between the engine, the loop and the caller.

THE QUESTION DISAPPEARED. On an exception the engine pops the user
message only when nothing streamed; when text HAD streamed it keeps the
user message and re-raises before the assistant message is appended. The
history then ends on a user message with no answer. The next send appends
another user message, and the alternation sanitiser resolves two
consecutive user messages by keeping the NEWEST -- so "continue" silently
replaced "refactor the parser", and the model answers against a history
that never mentioned the task. The chat still looked intact, because the
dashboard prints the partial text it received.

STOP DID NOT STOP. The engine checks its flag between STREAM EVENTS. The
tool loop that runs between those events never looked at it -- the name
does not occur in the client module at all. During a ten-minute command
or a stalled stream, Stop set a flag nobody read.

TWO TURNS COULD SHARE ONE ENGINE. Nothing refused a second turn while the
first was still running: both appended to the same message list, both
mutated the same per-turn cost state, and both ran tools in the same
workspace. The guard belonged on the engine, which owns the conversation,
not on a widget flag that a Stop had already cleared.
"""

from __future__ import annotations

import threading
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")


def _roles(engine):
    return [m.get("role") for m in engine.messages]


# ---------------------------------------------------------------------------
# The question survives
# ---------------------------------------------------------------------------

def test_a_turn_that_dies_after_streaming_keeps_both_sides(eng):
    def half_then_die(**kw):
        yield StreamEvent(type="text_delta", text="I read the parser and ")
        raise RuntimeError("connection reset")

    eng.client.stream_message = MagicMock(side_effect=half_then_die)
    with pytest.raises(RuntimeError):
        eng.stream_response("refactor the ORCA parser")

    assert _roles(eng) == ["user", "assistant"], (
        "the turn ended on a user message with no answer; the next send "
        "will replace the question")
    assert "parser" in eng.messages[0]["content"]
    assert "I read the parser and" in eng.messages[1]["content"]


def test_the_partial_answer_is_marked_as_incomplete(eng):
    def half_then_die(**kw):
        yield StreamEvent(type="text_delta", text="partial")
        raise RuntimeError("connection reset")

    eng.client.stream_message = MagicMock(side_effect=half_then_die)
    with pytest.raises(RuntimeError):
        eng.stream_response("do the thing")
    assert "partial" in eng.messages[-1]["content"]
    # The model must not read its own truncated text as a finished answer.
    assert "[" in eng.messages[-1]["content"]


def test_the_next_message_no_longer_replaces_the_question(eng):
    def half_then_die(**kw):
        yield StreamEvent(type="text_delta", text="started")
        raise RuntimeError("connection reset")

    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=half_then_die)
    with pytest.raises(RuntimeError):
        eng.stream_response("refactor the ORCA parser")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("continue")

    joined = " ".join(str(m.get("content")) for m in eng.messages)
    assert "refactor the ORCA parser" in joined, (
        "the original task was dropped by the alternation sanitiser")


def test_a_turn_that_dies_before_streaming_leaves_no_orphan(eng):
    def die(**kw):
        raise RuntimeError("refused")
        yield  # pragma: no cover

    eng.client.stream_message = MagicMock(side_effect=die)
    with pytest.raises(RuntimeError):
        eng.stream_response("hello")
    assert eng.messages == [], "an unanswered user message was left behind"


# ---------------------------------------------------------------------------
# Stop reaches the loop
# ---------------------------------------------------------------------------

def test_the_engine_hands_the_loop_a_way_to_see_a_stop(eng):
    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("hi")
    probe = getattr(eng.client, "should_stop", None)
    assert callable(probe), "the loop has no way to observe a stop"
    assert probe() is False
    eng.request_stop()
    assert probe() is True
    eng.clear_stop()
    assert probe() is False


def test_the_loop_checks_it_between_rounds():
    """A flag nobody reads is not a stop."""
    import pathlib
    src = pathlib.Path(
        __import__("delfin.agent.api_client", fromlist=["x"]).__file__
    ).read_text(encoding="utf-8")
    code = "\n".join(l for l in src.splitlines()
                     if not l.lstrip().startswith("#"))
    assert "should_stop" in code


def test_a_client_without_the_probe_still_runs():
    from delfin.agent import api_client as A
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    assert client._stop_was_requested() is False


def test_a_broken_probe_does_not_end_the_turn():
    from delfin.agent import api_client as A
    client = A.OpenAIClient.__new__(A.OpenAIClient)

    def boom():
        raise RuntimeError("engine gone")

    client.should_stop = boom
    assert client._stop_was_requested() is False


# ---------------------------------------------------------------------------
# One conversation, one turn
# ---------------------------------------------------------------------------

def test_a_second_turn_cannot_start_while_the_first_runs(eng):
    """Two workers on one engine interleave the history and the cost state."""
    started = threading.Event()
    release = threading.Event()

    def slow(**kw):
        started.set()
        release.wait(timeout=5)
        yield StreamEvent(type="text_delta", text="done")

    eng.client.stream_message = MagicMock(side_effect=slow)
    out: list[str] = []
    t = threading.Thread(target=lambda: out.append(eng.stream_response("first")))
    t.start()
    assert started.wait(timeout=5)

    second = eng.stream_response("second")
    assert "already" in second.lower() or "running" in second.lower(), second
    release.set()
    t.join(timeout=10)
    assert out and out[0] == "done"
    # The refused turn must not have left its message in the history.
    assert [m.get("content") for m in eng.messages if m.get("role") == "user"] \
        == ["first"]


def test_the_guard_clears_after_a_failed_turn(eng):
    def die(**kw):
        raise RuntimeError("boom")
        yield  # pragma: no cover

    eng.client.stream_message = MagicMock(side_effect=die)
    with pytest.raises(RuntimeError):
        eng.stream_response("first")

    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fine)
    assert eng.stream_response("second") == "ok", (
        "a failed turn left the engine permanently refusing new ones")
