"""Context-truth accounting: the compaction ladder must run on real
provider token counts on ALL backends, and in-place trims must be
visible to its stop conditions.

Two verified defects covered here:
1. The OpenAI-compatible client never emitted message_start, so on
   KIT/Ollama the engine's input accounting and compaction floor stayed 0
   and the 70%/95% triggers ran on a schema-blind chars//4 guess.
2. The floor (_last_input_tokens) is a PRE-trim snapshot and was never
   credited for chars removed by the sliding-window trim / hard clear, so
   the deterministic stages could not terminate and could never avert the
   paid LLM summary.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# quick mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                          mode="quick", pack_dir=agent_tree)


def _client_from(events):
    client = MagicMock()
    client.stream_message = MagicMock(
        side_effect=lambda *a, **k: iter(events))
    return client


# ---------------------------------------------------------------------------
# Engine: message_start / message_delta accounting
# ---------------------------------------------------------------------------

def test_message_start_feeds_floor_and_input(agent_tree):
    from delfin.agent.api_client import StreamEvent
    ev = [
        StreamEvent(type="message_start", input_tokens=1234),
        StreamEvent(type="text_delta", text="ok"),
        StreamEvent(type="message_delta", input_tokens=99999,
                    output_tokens=5, cost_usd=0.0),
    ]
    engine = _engine(agent_tree, _client_from(ev))
    engine.stream_response("hi")
    assert engine.token_usage["input"] == 1234       # delta input NOT added
    assert engine._last_input_tokens == 1234


def test_message_delta_fallback_counts_input_without_message_start(agent_tree):
    """Backends that report usage only in the final delta still get input
    accounting — but the compaction floor must NOT be set from the
    cumulative delta value."""
    from delfin.agent.api_client import StreamEvent
    ev = [
        StreamEvent(type="text_delta", text="ok"),
        StreamEvent(type="message_delta", input_tokens=777,
                    output_tokens=5, cost_usd=0.0),
    ]
    engine = _engine(agent_tree, _client_from(ev))
    engine.stream_response("hi")
    assert engine.token_usage["input"] == 777
    assert engine._last_input_tokens == 0            # floor untouched


def test_per_round_message_starts_sum_like_billing(agent_tree):
    """A multi-round tool turn emits one message_start per request. The
    engine SUMS them, because each round is a separately billed request --
    and takes its compaction floor from the FIRST one only.

    The floor used to track the latest round, which made it the intra-turn
    peak: every later round carries that turn's accumulated tool results,
    and those never enter self.messages. Measured on a real loop, the true
    next-request size was 1,155 tokens while the floor kept 17,634. The
    next turn's compaction then trimmed every older message on a
    conversation occupying 2% of the window, and the self-monitoring block
    told the agent it was at 89% and should wind down while it held one
    token of history."""
    from delfin.agent.api_client import StreamEvent
    ev = [
        StreamEvent(type="message_start", input_tokens=1000),
        StreamEvent(type="message_start", input_tokens=1600),
        StreamEvent(type="text_delta", text="ok"),
        StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0),
    ]
    engine = _engine(agent_tree, _client_from(ev))
    engine.stream_response("hi")
    assert engine.token_usage["input"] == 2600
    assert engine._last_input_tokens == 1000


# ---------------------------------------------------------------------------
# Engine: floor credit for in-place trims
# ---------------------------------------------------------------------------

def _quiet_engine(agent_tree):
    ev = []
    return _engine(agent_tree, _client_from(ev))


def test_slide_window_trim_credits_the_floor(agent_tree):
    engine = _quiet_engine(agent_tree)
    engine.context_window_tokens = 10_000
    big = "x" * 40_000
    engine.messages = (
        [{"role": "user", "content": "goal"},
         {"role": "assistant", "content": big}]
        + [{"role": "user", "content": f"u{i}"} for i in range(engine._KEEP_RECENT)]
    )
    engine._last_input_tokens = 20_000     # pre-trim provider snapshot
    before = engine._estimate_context_tokens()
    assert before == 20_000                # floor is binding
    trimmed = engine._slide_window_trim()
    assert trimmed >= 1
    assert engine._trimmed_chars_since_floor > 0
    after = engine._estimate_context_tokens()
    assert after < before                  # trim is now visible to the estimator


def test_hard_clear_credits_the_floor(agent_tree):
    engine = _quiet_engine(agent_tree)
    engine.context_window_tokens = 10_000
    big = "y" * 30_000
    msgs = [{"role": "assistant", "content": big}]
    engine.messages = msgs + [{"role": "user", "content": "goal"}]
    engine._last_input_tokens = 20_000
    cleared = engine._hard_clear_old_tool_results(msgs)
    assert cleared == 1
    assert engine._estimate_context_tokens() < 20_000


def test_fresh_message_start_resets_trim_credit(agent_tree):
    from delfin.agent.api_client import StreamEvent
    ev = [
        StreamEvent(type="message_start", input_tokens=5000),
        StreamEvent(type="text_delta", text="ok"),
        StreamEvent(type="message_delta", output_tokens=1, cost_usd=0.0),
    ]
    engine = _engine(agent_tree, _client_from(ev))
    engine._trimmed_chars_since_floor = 12_345
    engine.stream_response("hi")
    # New provider count already reflects earlier trims -> credit resets.
    assert engine._trimmed_chars_since_floor == 0
    assert engine._last_input_tokens == 5000


# ---------------------------------------------------------------------------
# OpenAI-compatible client: per-round message_start emission
# ---------------------------------------------------------------------------

class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    def __init__(self, prompt=0, completion=0):
        self.prompt_tokens = prompt
        self.completion_tokens = completion


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


@pytest.fixture
def openai_client(monkeypatch, tmp_path):
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _FakeRegistry())
    import delfin.agent.api_client as A
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)
    return A.create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))


def test_openai_client_emits_message_start_with_round_prompt_tokens(openai_client):
    def _create(**kwargs):
        return _Stream([_Chunk(
            [_Choice(_Delta(content="done"), finish="stop")],
            usage=_Usage(prompt=4321, completion=7))])

    openai_client.client.chat.completions.create = _create
    events = list(openai_client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=64))
    starts = [e for e in events if e.type == "message_start"]
    assert len(starts) == 1
    assert starts[0].input_tokens == 4321
    # cached stays on message_delta only (no double count via message_start)
    assert all((getattr(e, "cached_tokens", 0) or 0) == 0 for e in starts)
