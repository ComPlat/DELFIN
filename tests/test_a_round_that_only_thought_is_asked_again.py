"""A round that produced only reasoning is a round with no answer.

The streamed request to an OpenAI-compatible backend sometimes ends
having produced nothing usable. `api_client` already retries such a round
once, non-streaming, because the identical request answers when asked the
other way — measured at the protocol level on kit.glm-5.3, twelve times
out of twelve.

The retry refused to fire when the reasoning channel had carried
something, on the ground that a model thinking out loud is different from
a model producing nothing. It is not different in the way that matters:
reasoning is not the answer, it is not handed back to the caller, and
what the engine does with such a round is declare the whole turn empty,
take the user's message out of the history and tell them to send it
again.

Captured 2026-09-09 from a suite run: 36577 characters of system prompt,
one 44-character user message, 16 tools advertised — 34 characters of
reasoning, no text, no tool call, and `[empty turn]` reported for a
question the model could answer. Replaying that exact request answered
correctly six times out of six, streamed and non-streamed alike. The
request was fine; the response was empty once. That is why a retry is the
cure.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _should_retry_empty_round


def test_a_round_that_produced_only_reasoning_is_re_asked():
    """The case this exists for. Reasoning is not an argument against."""
    assert _should_retry_empty_round(True, [], {})


@pytest.mark.parametrize("text_chunks, tool_calls", [
    (["Gespeichert."], {}),
    ([], {"0": {"name": "remember"}}),
    (["part"], {"0": {"name": "bash"}}),
])
def test_a_round_that_produced_something_is_left_alone(text_chunks, tool_calls):
    """Content and a tool call are the two things a round can contribute;
    with either, there is something to hand back and a retry would throw
    it away."""
    assert not _should_retry_empty_round(True, text_chunks, tool_calls)


def test_a_non_streamed_round_is_not_retried():
    """The retry works by asking the other way. A non-streamed round has
    no other way — it would be the same request twice."""
    assert not _should_retry_empty_round(False, [], {})


def test_the_retry_is_wired_to_the_predicate():
    """The branch lives inside a long generator; if the call site stops
    using the predicate, the tests above stop describing anything.

    What the retry DOES is proven with a stub in
    tests/test_a_stream_that_said_nothing_is_asked_again.py — this only
    keeps the two from drifting apart.
    """
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client)
    assert "_should_retry_empty_round(" in src
    # ...and the fallback it guards is still the non-streaming one.
    idx = src.index("_should_retry_empty_round(\n")
    window = src[idx:idx + 700]
    assert '_nk["stream"] = False' in window, window[:200]
