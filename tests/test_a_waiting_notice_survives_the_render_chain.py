"""The waiting indicator survives the whole render chain.

M3 made stream_message emit StreamEvent(type="waiting", ...) while the
endpoint has not answered. That event is useless unless every layer
between the client and the screen carries it through without folding it
into the answer text:

* engine.stream_response must route "waiting" to on_notice (not
  on_token — a notice is not the model's answer, and a turn whose
  recorded output was three retry banners was once scored as a model
  that answered badly).
* cli_stream.StreamRenderer must print a "notice" event, not drop it.
* cli._run_once must pass on_notice through to the engine.
"""

from unittest.mock import MagicMock, patch

import delfin.agent.engine as E
from delfin.agent.cli_stream import StreamRenderer
from delfin.agent.api_client import StreamEvent

import pytest


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


def test_the_waiting_event_goes_to_on_notice_not_the_answer(eng):
    """'waiting for the model … 10s' is harness speech: shown to the
    user on the notice channel, never mixed into the answer text."""
    eng.client.stream_message = _stream(
        StreamEvent(type="waiting", text="waiting for the model … 10s"),
        StreamEvent(type="text_delta", text="the answer"),
    )
    answer: list[str] = []
    notices: list[str] = []
    eng.stream_response(
        "go", on_token=answer.append, on_notice=notices.append)

    assert "".join(answer).strip() == "the answer"
    assert notices == ["waiting for the model … 10s"]


def test_without_on_notice_the_waiting_line_stays_out_of_the_answer(eng):
    """The non-breaking guarantee: a caller that passes no on_notice
    (old dashboard, benchmark, scheduler) must not suddenly find the
    waiting line in the answer text — the fallback drops nothing into
    on_token that was not already there."""
    eng.client.stream_message = _stream(
        StreamEvent(type="waiting", text="waiting for the model … 10s"),
        StreamEvent(type="text_delta", text="the answer"),
    )
    answer: list[str] = []
    eng.stream_response("go", on_token=answer.append)
    assert "".join(answer).strip() == "the answer"


def test_stream_renderer_prints_a_notice_event():
    lines = list(StreamRenderer(
        [{"type": "notice", "text": "waiting for the model … 10s"}],
        is_tty=False))
    assert lines == ["waiting for the model … 10s"], (
        "the waiting notice must be shown, not dropped")


def test_stream_renderer_drops_an_empty_notice():
    lines = list(StreamRenderer(
        [{"type": "notice", "text": "   "}], is_tty=False))
    assert lines == []
