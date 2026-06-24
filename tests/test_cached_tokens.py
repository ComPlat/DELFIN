"""Capture prompt-cache hits (cached_tokens) on the OpenAI/KIT + Anthropic
paths so DELFIN can SEE how much input is served free from the prefix cache —
the foundation for the caching/efficiency work."""

from __future__ import annotations

from delfin.agent.api_client import _cached_tokens_of, StreamEvent


class _Det:
    def __init__(self, c):
        self.cached_tokens = c


class _Usage:
    def __init__(self, det=None):
        if det is not None:
            self.prompt_tokens_details = det


def test_extracts_cached_tokens():
    assert _cached_tokens_of(_Usage(_Det(1234))) == 1234


def test_zero_when_field_absent():
    assert _cached_tokens_of(_Usage()) == 0          # no prompt_tokens_details
    assert _cached_tokens_of(_Usage(_Det(None))) == 0
    assert _cached_tokens_of(None) == 0


def test_stream_event_has_cached_field():
    ev = StreamEvent(type="message_delta", input_tokens=100, cached_tokens=80)
    assert ev.cached_tokens == 80
    assert StreamEvent(type="message_delta").cached_tokens == 0   # default


def test_capture_is_wired_in_both_paths():
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    # OpenAI/KIT streaming + non-streaming capture
    assert src.count("_total_cached += _cached_tokens_of(") >= 2
    # Anthropic path surfaces cache_read
    assert "cached_tokens=cache_read" in src
    # OpenAI path surfaces the accumulator on the final delta
    assert "cached_tokens=_total_cached" in src


def test_engine_accumulates_and_exposes_cached():
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "engine.py").read_text(encoding="utf-8")
    assert 'self.token_usage = {"input": 0, "output": 0, "cached": 0}' in src
    assert 'self.token_usage.get(\n                            "cached"' in src \
        or 'self.token_usage.get("cached"' in src
    assert '"cached_tokens": self.token_usage.get("cached", 0)' in src
