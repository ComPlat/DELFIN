"""Controls for compactions in the round report (order LF, phase 4).

Whether and how often a session was compacted had to be dug out of the
session files by hand. Wave 5 (26.9.2026) gap 3.

Every summary compaction archives the pre-compaction transcript as one
line in ~/.delfin/transcript_archive/<sid>.jsonl (session_store.
archive_pre_compaction_transcript); the count is the line count. That
is the number reported -- in-place trims (elided store) are NOT
compaction events and do not appear there.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import round_report as rr
from delfin.agent import turn_metrics


NOW = 1_800_000_000.0


def _write_trace(base, session, entries):
    p = base / f"{session}.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e) + "\n")


def _te(tool, ts):
    return {"ts": ts, "tool": tool, "input": "", "output": "",
            "duration_ms": 10, "ok": True, "error": ""}


def _write_compactions(archive_dir, session, n):
    """n compaction events, one archive line each (the shape
    archive_pre_compaction_transcript writes)."""
    p = archive_dir / f"{session}.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a", encoding="utf-8") as f:
        for i in range(n):
            f.write(json.dumps({
                "compacted_at": NOW - 60 * (i + 1),
                "n_messages": 10,
                "info": {"messages_compacted": 10, "tokens_before": 1},
                "messages": [],
            }) + "\n")


@pytest.fixture
def state(tmp_path, monkeypatch):
    traces = tmp_path / "tool_traces"
    archive = tmp_path / "transcript_archive"
    metrics = tmp_path / "turn_metrics"
    monkeypatch.setattr(turn_metrics, "_DIR", metrics)

    _write_trace(traces, "nacht-a", [_te("read_file", NOW - 100)])
    _write_trace(traces, "nacht-b", [_te("read_file", NOW - 90)])
    _write_trace(traces, "alt-1", [_te("read_file", NOW - 10 * 86_400)])

    _write_compactions(archive, "nacht-a", 3)   # three compactions
    _write_compactions(archive, "nacht-b", 0)   # none: file is empty
    _write_compactions(archive, "alt-1", 2)     # before the cut: excluded
    return {"traces": traces, "archive": archive}


def _collect(state, **kw):
    return rr.collect(since_s=NOW - 3_600, trace_root=state["traces"],
                      **kw)


class TestCompactionsCollected:
    def test_counted_per_session(self, state):
        data = _collect(state, archive_root=state["archive"])
        by_id = {s["session_id"]: s for s in data["sessions"]}
        assert by_id["nacht-a"]["compactions"] == 3
        assert by_id["nacht-b"]["compactions"] == 0

    def test_zero_without_an_archive(self, state):
        data = _collect(state)      # no archive_root: dir does not exist
        by_id = {s["session_id"]: s for s in data["sessions"]}
        assert by_id["nacht-a"]["compactions"] == 0
        assert by_id["nacht-b"]["compactions"] == 0

    def test_summed_in_totals(self, state):
        data = _collect(state, archive_root=state["archive"])
        assert data["totals"]["compactions"] == 3

    def test_render_shows_compactions_per_session_and_total(self, state):
        data = _collect(state, archive_root=state["archive"])
        text = rr.render_text(data)
        assert "compactions: 3" in text
        assert "compactions" in text.split("TOTAL", 1)[1]

    def test_cli_collect_passes_the_real_archive_dir(self, state,
                                                     monkeypatch, capsys):
        """The public path: cmd_report --since must hand the REAL
        transcript archive to collect, not just the fabricated stores
        of the tests. Red before the wiring existed."""
        import argparse
        from delfin.agent import cli
        monkeypatch.setattr(rr, "collect",
                            lambda **kw: (_ for _ in ()).throw(
                                AssertionError(
                                    f"collect kwargs: {sorted(kw)}")))
        args = argparse.Namespace(since="90m", session="", name="",
                                  json=False)
        with pytest.raises(AssertionError) as exc:
            cli.cmd_report(args)
        assert "archive_root" in str(exc.value)
