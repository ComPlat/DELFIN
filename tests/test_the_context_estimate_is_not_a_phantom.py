"""The estimate that drives compaction measured the wrong thing.

The tool-loop backends emit one `message_start` per ROUND, and every
round after the first carries that turn's accumulated tool results. Those
results live only in the client's local request list -- they never enter
`self.messages` and are gone by the time the next request is built. The
engine kept the LAST one as its ground-truth floor, so the floor was the
intra-turn peak.

Measured on a real loop before the fix: true next-request size 1,155
tokens, floor kept 17,634. Fifteen times. The next turn's compaction then
saw pressure that did not exist and trimmed every message older than the
last four -- on a conversation occupying two percent of the window, every
turn, for the rest of the session. Above the summary threshold it also
bought an LLM summarisation call to destroy nine messages worth 1.5% of
the window. And the same estimate feeds the self-monitoring block, which
told the agent it was at 89% and should wind down and delegate while it
held one token of history.

The docstring already said the floor was "the last real count for THIS
context". Round N's prompt is this context plus ephemeral tool results,
so the code contradicted its own contract.

The direct Anthropic backend has no tool loop and was never affected;
every OpenAI-compatible provider and the CLI backend were.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


def _engine(tmp_path, events, mode="solo"):
    def stream(**kw):
        for ev in events:
            yield ev

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode=mode)
    return eng


# ---------------------------------------------------------------------------
# The floor is the request that will be repeated, not the one that grew
# ---------------------------------------------------------------------------

def test_the_floor_comes_from_the_first_round(tmp_path):
    eng = _engine(tmp_path, [
        StreamEvent(type="message_start", input_tokens=1_155),
        StreamEvent(type="message_start", input_tokens=9_898),
        StreamEvent(type="message_start", input_tokens=17_634),
        StreamEvent(type="text_delta", text="done"),
        StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0),
    ])
    eng.stream_response("read sixteen files")
    assert eng._last_input_tokens == 1_155


def test_billing_still_counts_every_round(tmp_path):
    """Each round IS a separately billed request; only the floor changed."""
    eng = _engine(tmp_path, [
        StreamEvent(type="message_start", input_tokens=1_155),
        StreamEvent(type="message_start", input_tokens=9_898),
        StreamEvent(type="message_start", input_tokens=17_634),
        StreamEvent(type="text_delta", text="done"),
        StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0),
    ])
    eng.stream_response("read sixteen files")
    assert eng.token_usage["input"] == 1_155 + 9_898 + 17_634


def test_the_next_turn_takes_its_own_floor(tmp_path):
    """The latch is per turn, not per session."""
    eng = _engine(tmp_path, [
        StreamEvent(type="message_start", input_tokens=1_000),
        StreamEvent(type="message_start", input_tokens=50_000),
        StreamEvent(type="text_delta", text="a"),
        StreamEvent(type="message_delta", output_tokens=1, cost_usd=0.0),
    ])
    eng.stream_response("first")
    assert eng._last_input_tokens == 1_000

    def stream2(**kw):
        yield StreamEvent(type="message_start", input_tokens=2_000)
        yield StreamEvent(type="text_delta", text="b")
        yield StreamEvent(type="message_delta", output_tokens=1, cost_usd=0.0)

    eng.client.stream_message = MagicMock(side_effect=stream2)
    eng.stream_response("second")
    assert eng._last_input_tokens == 2_000


def test_a_long_tool_turn_does_not_shred_a_short_conversation(tmp_path):
    """The behaviour the number was destroying."""
    events = [StreamEvent(type="message_start", input_tokens=n)
              for n in (900, 40_000, 90_000, 150_000)]
    events += [StreamEvent(type="text_delta", text="ok"),
               StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)]
    eng = _engine(tmp_path, events)
    eng.context_window_tokens = 200_000
    eng.stream_response("a turn with many tool rounds")
    assert eng._estimate_context_tokens() < 20_000, (
        "the estimate still carries the intra-turn peak")


# ---------------------------------------------------------------------------
# The distiller does not run where its only path is the lossy one
# ---------------------------------------------------------------------------

def test_it_is_off_on_a_provider_it_cannot_call(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode="solo")
    assert eng._distiller is not None
    assert eng._distiller.enabled is False


def test_it_stays_on_where_it_works(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api",
                            provider="claude", model="haiku", mode="solo")
    assert eng._distiller.enabled is True


def test_it_never_builds_a_client_for_another_provider():
    """Not merely broken: with an Anthropic key exported for unrelated
    reasons it would ship this engine's project memory, external memory
    and repo map to a provider the user never selected."""
    from delfin.agent.context_distiller import ContextDistiller

    d = ContextDistiller(enabled=True, provider="kit")
    with patch("anthropic.Anthropic") as spy:
        assert d._api_distill("some context") == ""
        assert spy.call_count == 0


def test_text_it_cannot_classify_is_left_alone():
    from delfin.agent.context_distiller import _split_layer0

    prompt = "role rules with no section marker anywhere " * 50
    layer0, rest = _split_layer0(prompt)
    assert layer0 == prompt and rest == ""
