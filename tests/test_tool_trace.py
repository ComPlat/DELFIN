"""Tool-call trace + bug-report integration."""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent import tool_trace as tt


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(tt, "_DIR", tmp_path / ".delfin" / "tool_traces")
    return tmp_path


def test_record_and_read(home):
    tt.record("s1", tool="read_file", tool_input={"path": "a.txt"},
              output="hello", duration_ms=5, ok=True)
    tt.record("s1", tool="bash", tool_input="ls", output="x", ok=False,
              error="blocked")
    e = tt.read("s1")
    assert len(e) == 2
    assert e[0]["tool"] == "read_file" and "a.txt" in e[0]["input"]
    assert e[0]["duration_ms"] == 5
    assert e[1]["ok"] is False and e[1]["error"] == "blocked"


def test_output_capped(home):
    tt.record("s2", tool="x", output="A" * 5000)
    assert len(tt.read("s2")[0]["output"]) <= 2000


def test_format_summary(home):
    tt.record("s3", tool="read_file", tool_input="p", output="ok", duration_ms=3)
    out = tt.format_summary(tt.read("s3"))
    assert "read_file" in out and "✓" in out


def test_read_empty_session(home):
    assert tt.read("nope") == []


def test_engine_records_trace(home, monkeypatch, tmp_path):
    from delfin.agent import model_capabilities as mc
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: mc.ModelCapabilities(
        model="m", provider="ollama", context_window=32768, supports_tools=True))
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine(repo_dir=str(tmp_path), backend="api",
                      provider="ollama", model="x")
    eng._trace_pending.append(("read_file", '{"path":"a"}', time.monotonic()))
    eng._record_tool_trace("read_file", "content here", ok=True)
    entries = tt.read(eng.trace_session())
    assert entries and entries[0]["tool"] == "read_file"
    assert "content here" in entries[0]["output"]


def test_bug_report_bundles_trace(home, tmp_path):
    tt.record("sess9", tool="read_file", tool_input={"path": "notes.txt"},
              output="the secret is GIRAFFE", ok=True)
    from delfin.agent.bug_report import write_bug_report
    rd = write_bug_report(chat_messages=[], session_id="sess9",
                          trace_session="sess9", description="x",
                          archive_dir=tmp_path / "arch")
    payload = json.loads((rd / "report.json").read_text(encoding="utf-8"))
    assert payload["tool_trace"] and payload["tool_trace"][0]["tool"] == "read_file"
    assert (rd / "tool_trace.jsonl").is_file()
    md = (rd / "report.md").read_text(encoding="utf-8")
    assert "Tool trace" in md and "read_file" in md
