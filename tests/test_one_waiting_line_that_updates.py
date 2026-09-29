"""The wait indicator is one line that updates, not a line per tick.

Input: a round in which the endpoint sends no byte for a while, so
``stream_message`` yields ``StreamEvent(type="waiting")`` repeatedly.
Output: each surface shows the elapsed time ONCE, updated in place.

Semantics. The producer's cadence is a tenth of the request deadline,
capped at ten seconds (api_client._WAIT_TICK_MIN_S/_MAX_S). With the
default 600 s deadline that is a waiting event every ten seconds, so a
round that stalls to the deadline yields sixty of them. Every surface
appended one message per event, and the chat filled with sixty copies of
the same sentence differing only in the number.

The event is progress, not history: it REPLACES the one before it. The
engine therefore routes it to ``on_wait`` when the caller offers one, and
only falls back to ``on_notice`` for callers that do not -- a fallback,
not a default, because a caller with no in-place surface still has to see
that something is happening.

Alternative considered and rejected: dropping every nth event in the
producer. It cuts the count without fixing the shape (the survivors still
append), and it makes the deadline check and the display share one
cadence, so shortening the display interval shortens the timeout.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

import delfin.agent.engine as E
from delfin.agent.api_client import StreamEvent
from delfin.agent.cli_stream import StreamRenderer


@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")


def _stream(*events):
    def _gen(**kwargs):
        for event in events:
            yield event
        yield StreamEvent(type="message_delta", input_tokens=5,
                          output_tokens=7, cost_usd=0.0)
    return _gen


_TICKS = tuple(
    StreamEvent(type="waiting", text=f"waiting for the model … {n}s")
    for n in (10, 20, 30))


def test_the_waiting_event_goes_to_on_wait_when_there_is_one(eng):
    """A caller with an in-place surface gets the ticks there, and its
    notice channel stays free for the harness speech that IS history."""
    eng.client.stream_message = _stream(
        *_TICKS, StreamEvent(type="notice", text="retrying 2/3 in 4s"),
        StreamEvent(type="text_delta", text="the answer"))
    waits: list[str] = []
    notices: list[str] = []
    answer: list[str] = []
    eng.stream_response("go", on_token=answer.append,
                        on_notice=notices.append, on_wait=waits.append)

    assert waits == [e.text for e in _TICKS]
    assert notices == ["retrying 2/3 in 4s"], (
        "a retry banner is history and stays on the notice channel; only "
        "the wait tick moves")
    assert "".join(answer).strip() == "the answer"


def test_without_on_wait_the_ticks_still_reach_the_notice_channel(eng):
    """The fallback. A caller that has no in-place surface must not go
    blind: it keeps exactly the behaviour it had before on_wait existed."""
    eng.client.stream_message = _stream(
        *_TICKS, StreamEvent(type="text_delta", text="the answer"))
    notices: list[str] = []
    answer: list[str] = []
    eng.stream_response("go", on_token=answer.append,
                        on_notice=notices.append)

    assert notices == [e.text for e in _TICKS]
    assert "".join(answer).strip() == "the answer"


def test_a_wait_event_never_reaches_the_answer(eng):
    """Unchanged guarantee, restated here because this file adds a
    second route to the same text: a caller passing neither callback
    must not find the wait line in what it recorded as the answer."""
    eng.client.stream_message = _stream(
        *_TICKS, StreamEvent(type="text_delta", text="the answer"))
    answer: list[str] = []
    eng.stream_response("go", on_token=answer.append)
    assert "".join(answer).strip() == "the answer"


def test_a_terminal_rewrites_the_line_instead_of_adding_one():
    """On a tty the renderer returns the tick as a carriage return plus
    an erase, so the cursor stays on one row. Three ticks, three writes,
    one line on the screen."""
    events = [{"type": "wait", "text": t.text} for t in _TICKS]
    out = list(StreamRenderer(events, is_tty=True))

    assert len(out) == 3
    for line, tick in zip(out, _TICKS):
        assert line.startswith("\r"), "the tick must return to column one"
        assert "\x1b[K" in line, "the previous, longer tick must be erased"
        assert tick.text in line
        assert not line.endswith("\n"), (
            "a newline would leave the old tick standing above the new one")


def test_a_pipe_gets_one_wait_line_and_not_a_flood():
    """A pipe has no cursor to move. The tick is progress and not a
    record, so the first one says the endpoint is slow and the rest are
    dropped rather than written sixty times into a log."""
    events = [{"type": "wait", "text": t.text} for t in _TICKS]
    out = list(StreamRenderer(events, is_tty=False))

    assert out == [_TICKS[0].text]


def test_a_notice_is_still_appended_on_both(eng):
    """The collapse applies to the wait tick alone. Two retry banners
    are two events and stay two lines, on a tty and in a pipe."""
    events = [{"type": "notice", "text": "retrying 1/3 in 2s"},
              {"type": "notice", "text": "retrying 2/3 in 4s"}]
    assert len(list(StreamRenderer(events, is_tty=False))) == 2
    assert len(list(StreamRenderer(events, is_tty=True))) == 2
