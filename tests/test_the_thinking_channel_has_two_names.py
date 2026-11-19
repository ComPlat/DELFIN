"""A reasoning channel is read by its content, not by one field name.

OpenAI-compatible endpoints stream a model's thinking in a field beside
``content``, and they do not agree on what to call it. Ollama (qwen3,
gpt-oss) and vLLM (r1) send ``reasoning_content``; the Helmholtz
Blablador gateway sends ``reasoning``.

Measured against Blablador on 2026-11-19, model alias-fast
(GPT-OSS-120b), one prompt with max_tokens=64: 28 stream chunks, of
which 24 carried ``reasoning`` and exactly one carried ``content``. A
reader that knows only the first spelling therefore sees nothing at all
in 24 of 28 chunks and drops the model's thinking.

Two consequences, and the second is the one that shows:

* the time-to-first-token marker never trips on thinking, so a model
  that thinks for seconds before answering is recorded as having taken
  that long to start;
* the thinking is not forwarded, so a round that spent its whole output
  budget reasoning reaches the caller as silence.

That second case is reachable without any unusual setting. Also
measured, same model and prompt: with max_tokens at 8, 16 or 32 the
response is HTTP 200 with ``content: null`` and ``finish_reason:
"length"`` -- the reasoning field alone consumed the budget. Only at 64
and above does ``content`` arrive. So on a short budget the answer is
empty by construction, and whether the person sees *why* depends
entirely on this field being read.
"""

from __future__ import annotations

import inspect

import pytest

from delfin.agent import api_client as ac


class _Delta:
    """A stream delta carrying exactly the attributes it was given.

    Not a dict: the reader uses getattr, because what the SDK hands back
    is an object whose absent fields raise rather than return None.
    """

    def __init__(self, **fields):
        for key, value in fields.items():
            setattr(self, key, value)


def test_the_ollama_and_vllm_spelling_is_read():
    assert ac._reasoning_delta(_Delta(reasoning_content="thinking")) == "thinking"


def test_the_blablador_spelling_is_read():
    assert ac._reasoning_delta(_Delta(reasoning="thinking")) == "thinking"


def test_a_delta_with_neither_is_empty_not_an_error():
    assert ac._reasoning_delta(_Delta(content="the answer")) == ""


def test_a_missing_delta_is_empty_rather_than_raising():
    """A chunk without a delta is normal at the end of a stream, and the
    reader must never be the thing that breaks one."""
    assert ac._reasoning_delta(None) == ""


@pytest.mark.parametrize("value", [None, "", 0, [], {}])
def test_an_empty_value_does_not_count_as_thinking(value):
    assert ac._reasoning_delta(_Delta(reasoning=value,
                                      reasoning_content=value)) == ""


def test_both_spellings_present_takes_the_first_non_empty():
    """No backend sends both. If one ever does, the behaviour is defined
    rather than dependent on attribute order."""
    assert ac._reasoning_delta(
        _Delta(reasoning_content="", reasoning="second")) == "second"
    assert ac._reasoning_delta(
        _Delta(reasoning_content="first", reasoning="second")) == "first"


def test_a_value_that_is_not_a_string_is_still_forwarded_as_one():
    """Forwarded into a StreamEvent whose text field is a str."""
    assert ac._reasoning_delta(_Delta(reasoning=42)) == "42"


def test_the_stream_reads_the_channel_through_the_helper():
    """The wiring, not the helper: a second hand-written getattr in the
    stream loop would keep exactly the defect this closes, and it would
    pass every test above.
    """
    src = inspect.getsource(ac)
    start = src.index("# Reasoning channel:")
    block = src[start:start + 700]

    assert "_reasoning_delta(delta)" in block, (
        "the stream loop must read the channel through the helper")
    assert 'getattr(delta, "reasoning_content"' not in block, (
        "a hand-written read of one spelling is what this closes")


def test_the_first_token_marker_counts_thinking_as_a_token():
    """Thinking is the first thing a reasoning model emits. A marker that
    ignores it measures the wrong thing on every such model."""
    src = inspect.getsource(ac)
    start = src.index("_round_ttft is None")
    block = src[start:start + 200]

    assert "_reasoning_delta(delta)" in block, (
        "the prefill ends at the first token of any kind, thinking "
        "included")


def test_both_spellings_are_declared_in_one_place():
    """Two readers with two lists is the shape that drifts."""
    assert ac._REASONING_FIELDS == ("reasoning_content", "reasoning")
