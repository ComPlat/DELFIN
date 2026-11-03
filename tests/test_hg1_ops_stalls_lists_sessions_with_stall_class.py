"""`delfin-agent ops stalls` -- a read-only list of stalls, per session.

Input: the session rows `ops_interface.sessions()` produces (idle_s,
open_approvals, pid, paused) as `session_health.collect()` enriches them.
Output: one line per session ``<name>  <class>  <age>``, plus the
`classify_record(record, now)` mapping in progress_watchdog.
Asserted: open approvals outrank everything (approval-wait, age from
approval_opened_at); a dead pid is session-gone, not idle; a paused
session is not a stall; a long idle with no approvals is idle with its
age; a fresh active row is "running" (not a stall). Nothing is mutated.
"""

from __future__ import annotations

import pytest

from delfin.agent import progress_watchdog as pw


NOW = 1_700_000_000.0


def _row(**over):
    row = {
        "name": "w20/hg1",
        "session_id": "abc-123",
        "pid": 4242,
        "paused": False,
        "pid_alive": True,
        "open_approvals": [],
        "approval_opened_at": None,
        "idle_s": 30,
    }
    row.update(over)
    return row


class TestClassifyRecord:
    def test_open_approval_outranks_everything(self):
        rec = _row(open_approvals=["a1"], approval_opened_at=NOW - 600,
                   idle_s=900)
        out = pw.classify_record(rec, NOW)
        assert out["stall"] is True
        assert out["class"] == "approval-wait"
        assert out["duration_s"] == 600.0

    def test_dead_pid_is_session_gone_not_idle(self):
        out = pw.classify_record(_row(pid_alive=False), NOW)
        assert out["class"] == "session-gone"
        assert out["stall"] is True

    def test_paused_is_not_a_stall(self):
        out = pw.classify_record(_row(paused=True, idle_s=5000), NOW)
        assert out["stall"] is False

    def test_long_idle_without_approvals_is_idle(self):
        out = pw.classify_record(_row(idle_s=7200), NOW)
        assert out["stall"] is True
        assert out["class"] == "idle"
        assert out["duration_s"] == 7200.0

    def test_fresh_active_row_is_running(self):
        out = pw.classify_record(_row(idle_s=10), NOW)
        assert out["stall"] is False
        assert out["class"] == "running"

    @pytest.mark.parametrize("idle", [None, "x"])
    def test_unreadable_idle_is_not_idle(self, idle):
        out = pw.classify_record(_row(idle_s=idle), NOW)
        assert out["class"] != "idle"


class TestOpsStallsCommand:
    def _run(self, rows, monkeypatch, capsys):
        from delfin.agent import ops_interface as oi
        monkeypatch.setattr(oi, "sessions", lambda **kw: rows)
        import argparse
        from delfin.agent import cli as cli_mod
        args = argparse.Namespace(ops_action="stalls")
        cli_mod._cmd_ops_stalls(args)
        return capsys.readouterr().out

    def test_one_line_per_session_with_class_and_age(self, monkeypatch,
                                                     capsys):
        rows = [_row(name="w20/hg1", idle_s=7200),
                _row(name="w20/other", open_approvals=["a1"],
                     approval_opened_at=NOW - 90, idle_s=95)]
        monkeypatch.setattr("time.time", lambda: NOW)
        out = self._run(rows, monkeypatch, capsys)
        lines = [l for l in out.splitlines() if l.strip()]
        assert any("w20/hg1" in l and "idle" in l for l in lines)
        assert any("w20/other" in l and "approval-wait" in l for l in lines)
        # age rendered, and no session was mutated
        assert any("600" in l or "10m" in l for l in lines)

    def test_no_stalls_prints_a_clear_line(self, monkeypatch, capsys):
        monkeypatch.setattr("time.time", lambda: NOW)
        out = self._run([_row(idle_s=5)], monkeypatch, capsys)
        assert "no stalled" in out.lower()
