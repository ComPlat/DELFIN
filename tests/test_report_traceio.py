"""Tests for delfin.agent.report_traceio (fabricated traces in tmp_path)."""

from __future__ import annotations

import json

from delfin.agent import report_traceio


def _make_store(tmp_path, n_files=3, entries_per_file=4):
    """Fabricate a trace store; returns (root, {session: n_entries})."""
    counts = {}
    for i in range(n_files):
        sid = f"session{i:02d}"
        lines = [
            json.dumps({"ts": 1.0, "tool": "read_file", "input": "",
                        "output": "", "duration_ms": 1, "ok": True,
                        "error": ""})
            for _ in range(entries_per_file + i)
        ]
        # one malformed line per file: readers must skip, not crash
        lines.append("{not json")
        (tmp_path / f"{sid}.jsonl").write_text("\n".join(lines) + "\n",
                                               encoding="utf-8")
        counts[sid] = entries_per_file + i
    (tmp_path / "ignored.txt").write_text("x", encoding="utf-8")
    return tmp_path, counts


def test_collect_counts_and_timings(tmp_path):
    root, counts = _make_store(tmp_path)
    data = report_traceio.collect(root=root, repeats=3)
    assert data["n_files"] == 3
    assert data["n_entries"] == sum(counts.values())   # malformed skipped
    assert data["repeats"] == 3
    assert data["total_mb"] > 0
    for key in ("module_read", "plain_walk"):
        t = data[key]
        assert t["runs"] == 3.0
        assert t["min_ms"] <= t["median_ms"] <= t["max_ms"]
        assert t["spread_ms"] == t["max_ms"] - t["min_ms"]


def test_collect_both_readers_agree_on_entries(tmp_path):
    """Whatever the timing, the two ways must see the same store."""
    root, _ = _make_store(tmp_path, n_files=2, entries_per_file=5)
    data = report_traceio.collect(root=root, repeats=2)
    # module_read populated n_entries; plain walk must match
    plain = report_traceio._plain_walk(root)
    assert data["n_entries"] == sum(plain) == 11  # 5 + 6 per-file entries


def test_collect_empty_and_missing_root(tmp_path):
    data = report_traceio.collect(root=tmp_path / "nope")
    assert data["n_files"] == 0
    assert data["module_read"] == {}
    assert data["plain_walk"] == {}
    text = report_traceio.format_text(data)
    assert "nothing to measure" in text
    empty = tmp_path / "empty"
    empty.mkdir()
    data2 = report_traceio.collect(root=empty)
    assert data2["n_files"] == 0


def test_format_text_never_contains_trace_content(tmp_path):
    """Sizes, counts and timings only -- no tool names, inputs or outputs."""
    root, _ = _make_store(tmp_path)
    # plant recognizable content
    secret = "TOPSECRET_MARKER_CONTENT"
    (root / "s3.jsonl").write_text(
        json.dumps({"tool": "bash", "input": secret, "output": secret,
                    "ok": True}) + "\n", encoding="utf-8")
    text = report_traceio.format_text(report_traceio.collect(root=root))
    assert secret not in text
    assert "bash" not in text
    for needle in ("median", "spread", "module (root=)",
                   "plain open()+loads"):
        assert needle in text


def test_format_text_ratio_line(tmp_path):
    root, _ = _make_store(tmp_path, n_files=2, entries_per_file=3)
    text = report_traceio.format_text(
        report_traceio.collect(root=root, repeats=2))
    assert "module/plain ratio:" in text


def test_main_ok_and_missing(tmp_path, capsys):
    root, _ = _make_store(tmp_path, n_files=1, entries_per_file=1)
    assert report_traceio.main([str(root)]) == 0
    out = capsys.readouterr().out
    assert "Tool-trace read timing" in out
    assert report_traceio.main([]) == 0            # real store, never raises
    assert report_traceio.main([str(tmp_path / "nope")]) == 0
    assert "nothing to measure" in capsys.readouterr().out


def test_plain_walk_skips_bad_lines_and_non_jsonl(tmp_path):
    root, _ = _make_store(tmp_path, n_files=1, entries_per_file=2)
    counts = report_traceio._plain_walk(root)
    assert counts == [2]          # 2 good + 1 malformed -> 2; .txt ignored
