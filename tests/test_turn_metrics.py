"""Per-turn timing metrics: record/read, stall detection, summary."""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import turn_metrics as tm


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(tm, "_DIR", tmp_path / ".delfin" / "turn_metrics")
    return tmp_path


def test_record_and_read(home):
    tm.record("s1", provider="kit", model="azure.gpt-5.4", role="dashboard_agent",
              total_ms=4200, ttft_ms=3300, output_chars=6, tool_calls=0)
    tm.record("s1", provider="ollama", model="qwen3-vl:4b", role="solo_agent",
              total_ms=12000, ttft_ms=800, output_chars=50, tool_calls=3)
    e = tm.read("s1")
    assert len(e) == 2
    assert e[0]["model"] == "azure.gpt-5.4" and e[0]["ttft_ms"] == 3300
    assert e[1]["tool_calls"] == 3


def test_ttft_none_is_preserved(home):
    # A turn that emitted nothing (pure stall / error) records ttft=None.
    tm.record("s2", model="m", total_ms=90000, ttft_ms=None, output_chars=0)
    assert tm.read("s2")[0]["ttft_ms"] is None


def test_is_stall_flags_backend_wait(home):
    stall = {"ttft_ms": 89000, "tool_calls": 0, "output_chars": 6}
    fast = {"ttft_ms": 800, "tool_calls": 0, "output_chars": 6}
    busy = {"ttft_ms": 89000, "tool_calls": 4, "output_chars": 6}   # tools ran
    chatty = {"ttft_ms": 89000, "tool_calls": 0, "output_chars": 5000}  # generated a lot
    assert tm.is_stall(stall) is True
    assert tm.is_stall(fast) is False
    assert tm.is_stall(busy) is False
    assert tm.is_stall(chatty) is False
    assert tm.is_stall({"ttft_ms": None, "tool_calls": 0, "output_chars": 0}) is False


def test_format_summary_marks_stall(home):
    tm.record("s3", model="azure.gpt-5.4", total_ms=92700, ttft_ms=92000,
              output_chars=6, tool_calls=0)            # the 92.7s "Hallo" shape
    tm.record("s3", model="azure.gpt-5.4", total_ms=4000, ttft_ms=3300,
              output_chars=20, tool_calls=0)
    out = tm.format_summary(tm.read("s3"))
    assert "backend-stall" in out
    assert out.count("backend-stall") == 1         # only the slow one flagged


def test_read_empty_session(home):
    assert tm.read("nope") == []
