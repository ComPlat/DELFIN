"""Tests for delfin.agent.report_tools with fabricated trace stores."""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import report_tools


def _catalogue(*names: str) -> list[dict]:
    return [
        {"type": "function",
         "function": {"name": n, "description": n, "parameters": {}}}
        for n in names
    ]


def _write_trace(root, session: str, entries: list[dict]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    p = root / f"{session}.jsonl"
    lines = []
    for e in entries:
        # allow fabricating malformed lines via a raw string
        if isinstance(e, str):
            lines.append(e)
        else:
            lines.append(json.dumps(e))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _entry(tool: str, *, ok=True, ts=None, **extra) -> dict:
    e = {"ts": ts if ts is not None else time.time(), "tool": tool,
         "ok": ok, "input": "", "output": ""}
    e.update(extra)
    return e


CAT = _catalogue("bash", "read_file", "search_docs", "never_touched")


def test_counts_calls_failures_and_share(tmp_path):
    _write_trace(tmp_path, "s1", [
        _entry("bash"), _entry("bash"), _entry("bash", ok=False),
        _entry("read_file"),
        _entry("search_docs", ok=False), _entry("search_docs", ok=False),
    ])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert list(data["tools"]) == ["bash", "search_docs", "read_file"]
    assert data["tools"]["bash"] == {
        "calls": 3, "failures": 1, "failure_share": pytest.approx(1 / 3)}
    assert data["tools"]["search_docs"]["failure_share"] == pytest.approx(1.0)
    assert data["total_calls"] == 6
    assert data["sessions"] == 1


def test_never_used_comes_from_advertised_catalogue(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("bash")])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["never_used"] == ["never_touched", "read_file", "search_docs"]
    assert data["n_advertised"] == 4
    assert data["n_never_used"] == 3
    assert data["n_used_advertised"] == 1


def test_default_catalogue_is_the_client_advertised_set(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("bash")])
    data = report_tools.collect(root=tmp_path)
    # Not a typed-out list: derives from api_client's own surface.
    assert data["n_advertised"] > 10
    assert "bash" not in data["never_used"]
    assert "search_docs" in report_tools._advertised_names()
    assert isinstance(data["never_used"], list)


def test_mcp_prefix_stripped_for_comparison(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("mcp__delfin-ops__list_tools")])
    data = report_tools.collect(root=tmp_path, catalogue=CAT + _catalogue("list_tools"))
    assert "list_tools" in data["tools"]
    assert "list_tools" not in data["never_used"]


def test_called_but_not_advertised_is_reported(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("ghost_tool")])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["called_not_advertised"] == ["ghost_tool"]


def test_since_s_window_filters_old_entries(tmp_path):
    now = time.time()
    _write_trace(tmp_path, "s1", [
        _entry("bash", ts=now - 999_999),
        _entry("read_file", ts=now - 10),
    ])
    data = report_tools.collect(since_s=3600.0, root=tmp_path, catalogue=CAT)
    assert list(data["tools"]) == ["read_file"]
    assert data["window_s"] == 3600.0
    # Without a window both count.
    data_all = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data_all["total_calls"] == 2


def test_empty_and_malformed_store(tmp_path):
    _write_trace(tmp_path, "s1", ["{not json", json.dumps({"tool": "bash"}),
                                   json.dumps({"ok": True})])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["tools"]["bash"]["calls"] == 1  # missing ts counts (no window)
    assert data["total_calls"] == 1


def test_empty_store_yields_all_advertised_unused(tmp_path):
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["total_calls"] == 0
    assert data["n_never_used"] == 4
    text = report_tools.format_text(data)
    assert "No tool calls recorded." in text
    assert "- never_touched" in text


def test_format_text_orders_most_used_first_and_shows_failures():
    text = report_tools.format_text({
        "tools": {
            "bash": {"calls": 10, "failures": 0, "failure_share": 0.0},
            "grep_file": {"calls": 4, "failures": 2, "failure_share": 0.5},
        },
        "never_used": ["never_touched"], "n_advertised": 3,
        "n_never_used": 1, "total_calls": 14, "sessions": 2,
        "window_s": None,
    })
    assert text.index("bash") < text.index("grep_file")
    assert "50.0%" in text
    assert "0.0%" in text
    assert "never called: 1 of 3 (33.3%)" in text
    assert "- never_touched" in text


def test_main_prints_and_returns_zero(tmp_path, capsys):
    _write_trace(tmp_path, "s1", [_entry("bash")])
    rc = report_tools.main(["--trace-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Tool usage" in out
    assert "never called" in out


def test_real_store_never_raises():
    # The default store: whatever is on this machine, collect must not raise.
    data = report_tools.collect()
    assert isinstance(data["tools"], dict)
    assert data["n_advertised"] >= 50
