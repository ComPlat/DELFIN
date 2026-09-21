"""Session B tests: render_markdown / render_terminal / write_session_report.

The renderers are pure functions over the shared SessionReport contract
(defined by session A); tests build reports by hand and assert sections.
The write hook is best-effort and must never raise.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent.session_report import (
    SessionReport,
    render_markdown,
    render_terminal,
    write_session_report,
)


def _sample_report() -> SessionReport:
    now = time.time()
    return SessionReport(
        session_id="sess-123",
        model="kit.glm-5.3",
        started_at=now - 3665.0,
        ended_at=now,
        tool_calls=[
            {"name": "read_file", "count": 5, "ok": 5, "failed": 0},
            {"name": "edit_file", "count": 3, "ok": 2, "failed": 1},
        ],
        files_changed=[
            {"path": "delfin/agent/session_report.py", "change": "modified"},
            {"path": "tests/test_render.py", "change": "created"},
        ],
        commands_run=["grep -n def delfin/agent/cli.py", "pytest -q"],
        tests_run=[
            {"target": "tests/test_session_report.py", "status": "passed",
             "passed": 5, "failed": 0},
            {"target": "tests/test_other.py", "status": "failed",
             "passed": 2, "failed": 1},
        ],
        denials=[{"kind": "deny_pattern", "detail": "rm -rf blocked"}],
        cost_usd=4.5,
        input_tokens=12000,
        output_tokens=3456,
    )


def _empty_report() -> SessionReport:
    return SessionReport(session_id="")


class TestRenderMarkdown:
    def test_header_contains_session_model_duration(self):
        md = render_markdown(_sample_report())
        assert "# Session Report: sess-123" in md
        assert "**Model:** kit.glm-5.3" in md
        assert "**Duration:** 1h 1m 5s" in md

    def test_all_sections_present(self):
        md = render_markdown(_sample_report())
        for section in ("## Tool Calls", "## Files Changed",
                        "## Commands Run", "## Tests Run",
                        "## Denials", "## Cost & Tokens"):
            assert section in md

    def test_tool_table_rows(self):
        md = render_markdown(_sample_report())
        assert "| Tool | Calls | OK | Failed |" in md
        assert "| read_file | 5 | 5 | 0 |" in md
        assert "| edit_file | 3 | 2 | 1 |" in md

    def test_files_table_rows(self):
        md = render_markdown(_sample_report())
        assert "| delfin/agent/session_report.py | modified |" in md
        assert "| tests/test_render.py | created |" in md

    def test_commands_in_code_block(self):
        md = render_markdown(_sample_report())
        assert "```\ngrep -n def delfin/agent/cli.py\npytest -q\n```" in md

    def test_tests_table_rows(self):
        md = render_markdown(_sample_report())
        assert "| tests/test_session_report.py | passed | 5 | 0 |" in md
        assert "| tests/test_other.py | failed | 2 | 1 |" in md

    def test_denials_listed(self):
        md = render_markdown(_sample_report())
        assert "**deny_pattern**: rm -rf blocked" in md

    def test_cost_and_tokens(self):
        md = render_markdown(_sample_report())
        assert "**Cost:** $4.50" in md
        assert "**Input tokens:** 12,000" in md
        assert "**Output tokens:** 3,456" in md

    def test_empty_report_renders_none_markers(self):
        md = render_markdown(_empty_report())
        assert "(unknown session)" in md
        # Every section still exists, each with a _(none)_ marker.
        assert md.count("_(none)_") == 5  # tools, files, commands, tests, denials
        assert "**Model:** -" in md
        assert "**Duration:** -" in md

    def test_pipe_in_cell_is_escaped(self):
        report = SessionReport(
            session_id="s", files_changed=[{"path": "a|b", "change": "modified"}])
        md = render_markdown(report)
        assert "a\\|b" in md


class TestRenderTerminal:
    def test_is_plain_ascii(self):
        text = render_terminal(_sample_report())
        text.encode("ascii")  # raises if any colour / non-ascii crept in

    def test_three_compact_lines(self):
        lines = render_terminal(_sample_report()).splitlines()
        assert len(lines) == 3

    def test_summary_line_contents(self):
        lines = render_terminal(_sample_report()).splitlines()
        assert lines[0].startswith("Session sess-123")
        assert "model kit.glm-5.3" in lines[0]
        assert "duration 1h 1m 5s" in lines[0]
        assert "Tool calls: 8 (1 failed)" in lines[1]
        assert "files changed: 2" in lines[1]
        assert "tests: 1 passed, 1 failed" in lines[1]
        assert "denials: 1" in lines[1]
        assert "Cost $4.50" in lines[2]
        assert "tokens in 12,000 out 3,456" in lines[2]

    def test_empty_report_does_not_crash(self):
        lines = render_terminal(_empty_report()).splitlines()
        assert lines[0].startswith("Session (unknown)")
        assert "duration -" in lines[0]
        assert "Tool calls: 0 (0 failed)" in lines[1]


class TestWriteSessionReport:
    def test_writes_markdown_file(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "delfin.agent.session_report.collect_session_report",
            lambda sid: _sample_report())
        monkeypatch.setattr(
            "delfin.agent.session_report._report_dir",
            lambda: tmp_path)
        path = write_session_report("sess-123")
        assert path is not None and path.exists()
        text = path.read_text(encoding="utf-8")
        assert "# Session Report: sess-123" in text
        assert "## Cost & Tokens" in text
        assert not list(tmp_path.glob("*.tmp"))  # no leftover temp file

    def test_overwrite_is_latest_wins(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "delfin.agent.session_report.collect_session_report",
            _sample_report)
        monkeypatch.setattr(
            "delfin.agent.session_report._report_dir",
            lambda: tmp_path)
        write_session_report("sess-123")
        other = _sample_report()
        other.model = "changed"
        monkeypatch.setattr(
            "delfin.agent.session_report.collect_session_report",
            lambda sid: other)
        write_session_report("sess-123")
        assert "changed" in (tmp_path / "sess-123.md").read_text()

    def test_empty_session_id_returns_none(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "delfin.agent.session_report._report_dir", lambda: tmp_path)
        assert write_session_report("") is None
        assert write_session_report("   ") is None

    def test_never_raises_on_collector_error(self, tmp_path, monkeypatch):
        def boom(sid):
            raise RuntimeError("collector exploded")
        monkeypatch.setattr(
            "delfin.agent.session_report.collect_session_report", boom)
        monkeypatch.setattr(
            "delfin.agent.session_report._report_dir", lambda: tmp_path)
        assert write_session_report("sess-123") is None

    def test_unsafe_session_id_cannot_traverse(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "delfin.agent.session_report.collect_session_report",
            lambda sid: _sample_report())
        monkeypatch.setattr(
            "delfin.agent.session_report._report_dir", lambda: tmp_path)
        path = write_session_report("../../etc/passwd")
        assert path is not None
        assert path.parent == tmp_path
        assert ".." not in path.name
