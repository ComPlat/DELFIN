"""The system prompt changed every turn, so the prefix cache missed every turn.

The endpoint caches the longest common prefix of system prompt and messages,
and the system prompt comes first. One of its blocks -- the context status,
with a usage counter that moves every turn -- was rebuilt into it each turn,
so the first request of every turn was cold for the whole history behind it.
Measured on a live engine: two plain turns diverged at 60,503 of 61,425
prompt characters, inside that block. Three cluster sessions read 3.5M, 2.8M
and 10.2M input tokens for 10, 8 and 5 turns.

The volatile blocks ride at the end of the turn's user message now, kept on
the stored message as a private key so later requests repeat the bytes the
endpoint already holds, and dropped from older messages only when a trim
rewrites the history anyway. Driven through the real engine with a client
that records what it is sent.
"""

from __future__ import annotations

import copy
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def eng(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        e = E.AgentEngine(
            repo_dir=ws, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")
    # Real permissions on the scratch workspace: a mock's workspace puts
    # the task store under the process CWD, inside the checkout.
    e.client._permissions = A.KitToolPermissions(workspace=ws)
    return e


def _record(eng) -> list[dict]:
    """A client that keeps every request's system prompt and messages."""
    seen: list[dict] = []

    def _gen(**kw):
        seen.append({"system": kw.get("system"),
                     "messages": copy.deepcopy(kw.get("messages"))})
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", input_tokens=5,
                          output_tokens=7, cost_usd=0.0)

    eng.client.stream_message = _gen
    return seen


def _last_user(messages: list[dict]) -> dict:
    return next(m for m in reversed(messages) if m.get("role") == "user")


# ---------------------------------------------------------------------------
# The prefix
# ---------------------------------------------------------------------------

def test_two_plain_turns_share_the_system_prompt_byte_for_byte(eng):
    seen = _record(eng)
    eng.stream_response("Hallo")
    eng.stream_response("Weiter")
    assert len(seen) == 2
    assert seen[0]["system"] == seen[1]["system"], (
        "the system prompt still changes between two plain turns, so the "
        "endpoint's cache of it (and of every message after it) is lost")
    assert "Current usage" not in seen[1]["system"]


def test_the_blocks_still_reach_the_model_at_the_end_of_the_message(eng):
    seen = _record(eng)
    eng.stream_response("Hallo")
    eng.stream_response("Weiter")
    tail = _last_user(seen[1]["messages"])["content"]
    assert tail.startswith("Weiter\n\n"), tail[:60]
    assert "# Context status" in tail and "Current usage: 3 msgs" in tail
    # Each turn's own numbers: this is what used to move the prompt.
    first = _last_user(seen[0]["messages"])["content"]
    assert "Current usage: 1 msgs" in first


def test_an_earlier_message_is_sent_as_it_was_cached(eng):
    """The second request repeats the first turn's message, tail and all:
    that is the prefix the endpoint has, and it must not move."""
    seen = _record(eng)
    eng.stream_response("Hallo")
    eng.stream_response("Weiter")
    assert seen[1]["messages"][0] == seen[0]["messages"][0]


def test_the_stored_text_is_what_the_user_wrote(eng):
    """The tail is a private key, not part of the content: the triggers,
    the language matchers and the transcript read the content."""
    _record(eng)
    eng.stream_response("Hallo")
    stored = _last_user(eng.messages)
    assert stored["content"] == "Hallo"
    assert "# Context status" in stored["_steer"]
    # And the private key never reaches a provider as a field.
    assert all("_steer" not in m for m in eng._wire_messages())


def test_the_dashboard_state_rides_in_the_tail_too(eng):
    seen = _record(eng)
    eng.set_live_state("calc_dir: /scratch/run_042")
    eng.stream_response("Hallo")
    assert "calc_dir: /scratch/run_042" not in seen[0]["system"]
    assert "calc_dir: /scratch/run_042" in _last_user(
        seen[0]["messages"])["content"]


def test_a_rebuild_within_the_turn_keeps_the_tail_it_has(eng):
    """A continuation rebuilds the prompt mid-turn; the message's bytes
    must not change under the endpoint's cache."""
    _record(eng)
    eng.stream_response("Hallo")
    before = _last_user(eng.messages)["_steer"]
    eng.messages.append({"role": "user", "content": "x"})
    eng.messages[-1]["_steer"] = "frozen"
    eng._build_current_system_prompt("", task_text="x")
    assert eng.messages[-1]["_steer"] == "frozen"
    assert _last_user(eng.messages[:-1])["_steer"] == before


# ---------------------------------------------------------------------------
# Outside a turn, and the shapes
# ---------------------------------------------------------------------------

def test_a_prompt_built_with_no_message_still_carries_the_blocks(eng):
    prompt = eng._build_current_system_prompt("", task_text="x")
    assert "# Context status" in prompt


def test_list_content_gets_the_tail_as_a_text_block(eng):
    eng.messages = [{"role": "user", "_steer": "TAIL", "content": [
        {"type": "text", "text": "look"},
        {"type": "image_url", "image_url": {"url": "data:,"}}]}]
    wire = eng._wire_messages()
    assert wire[0]["content"][-1] == {"type": "text", "text": "TAIL"}
    assert len(wire[0]["content"]) == 3
    assert "_steer" not in wire[0]


def test_the_estimate_counts_the_tail():
    eng = E.AgentEngine.__new__(E.AgentEngine)
    eng.messages = [{"role": "user", "content": "hi"}]
    base = eng._estimate_context_tokens()
    eng.messages[0]["_steer"] = "T" * 4000
    assert eng._estimate_context_tokens() == base + 1000


# ---------------------------------------------------------------------------
# The cost side: the tails do not pile up past a trim
# ---------------------------------------------------------------------------

def test_a_trim_drops_the_stale_tails_and_keeps_the_newest(eng):
    _record(eng)
    for text in ("eins", "zwei", "drei"):
        eng.stream_response(text)
    users = [m for m in eng.messages if m.get("role") == "user"]
    assert all("_steer" in m for m in users)
    # Pressure between the sliding line and the compaction cliff.
    est = eng._estimate_context_tokens()
    eng.context_window_tokens = int(est / 0.80)
    eng.auto_compact_pct = 0.95
    assert eng._should_slide() and not eng._should_auto_compact()
    eng._compact_history()
    assert "_steer" not in users[0] and "_steer" not in users[1]
    assert "_steer" in users[2], "the newest turn lost its own tail"


def test_dropping_tails_credits_the_input_floor():
    eng = E.AgentEngine.__new__(E.AgentEngine)
    eng.messages = [
        {"role": "user", "content": "a", "_steer": "x" * 300},
        {"role": "assistant", "content": "b"},
        {"role": "user", "content": "c", "_steer": "y" * 100},
    ]
    eng._trimmed_chars_since_floor = 0
    assert eng._drop_stale_steering_tails() == 300
    assert eng._trimmed_chars_since_floor == 300
    assert eng.messages[2]["_steer"] == "y" * 100
    assert eng._drop_stale_steering_tails() == 0
