"""Per-turn timing metrics: record/read, stall detection, summary."""

from __future__ import annotations

import json
import time
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


# ---------------------------------------------------------------------------
# aggregate_turn_stats — windowed cross-session roll-up
# ---------------------------------------------------------------------------

def test_aggregate_counts_stalls_stops_and_windows(home):
    now = time.time()
    # session a: one normal turn, one backend stall
    tm.record("a", model="m", total_ms=4000, ttft_ms=1000,
              output_chars=50, tool_calls=1)
    tm.record("a", model="m", total_ms=90000, ttft_ms=89000,
              output_chars=6, tool_calls=0)             # stall shape
    # session b: a stopped turn + one 30-day-old entry (outside window)
    tm.record("b", model="m", total_ms=5000, ttft_ms=2000,
              output_chars=80, tool_calls=0, stopped=True)
    old = {"ts": now - 30 * 86400, "model": "m", "total_ms": 1,
           "ttft_ms": 1, "output_chars": 1, "tool_calls": 0,
           "stopped": True}
    with tm.metrics_path("b").open("a", encoding="utf-8") as f:
        f.write(json.dumps(old) + "\n")
    s = tm.aggregate_turn_stats(7, now=now)
    assert s["turns"] == 3                              # old entry excluded
    assert s["stalls"] == 1
    assert s["stopped_count"] == 1
    assert s["avg_ttft_ms"] == int(round((1000 + 89000 + 2000) / 3))
    assert s["p90_ttft_ms"] == 89000


def test_aggregate_skips_ttft_none_but_counts_turn(home):
    tm.record("s", model="m", total_ms=90000, ttft_ms=None, output_chars=0)
    tm.record("s", model="m", total_ms=2000, ttft_ms=500, output_chars=20)
    s = tm.aggregate_turn_stats()
    assert s["turns"] == 2
    assert s["avg_ttft_ms"] == 500                      # None excluded from avg


def test_aggregate_tolerates_corrupt_lines(home):
    tm.record("s", model="m", total_ms=100, ttft_ms=50, output_chars=5)
    with tm.metrics_path("s").open("a", encoding="utf-8") as f:
        f.write("garbage\n[1, 2]\n")
    s = tm.aggregate_turn_stats()
    assert s["turns"] == 1


def test_aggregate_empty_and_missing_dir_never_raise(home):
    # "crashes" joined the roll-up when turns that RAISED stopped being
    # counted as backend stalls: without its own counter, a run of
    # crashing turns would have made this report quieter than before.
    #
    # "never_started" and "ttft_sample" joined it for the same reason,
    # one layer along. A turn that waited out a silent backend was only
    # reachable through the public predicate, so the roll-up -- which is
    # what a person actually reads -- could not say it had happened. And
    # a mean over three ttft samples out of ninety turns reads exactly
    # like a mean over ninety; the sample size is the difference between
    # a number and a number that means something.
    zeros = {"turns": 0, "avg_ttft_ms": 0, "p90_ttft_ms": 0,
             "stalls": 0, "crashes": 0, "stopped_count": 0,
             "never_started": 0, "ttft_sample": 0}
    assert tm.aggregate_turn_stats() == zeros           # dir doesn't exist yet
    tm.record("s", model="m", total_ms=1)
    tm.metrics_path("s").write_text("", encoding="utf-8")
    assert tm.aggregate_turn_stats() == zeros           # empty file
