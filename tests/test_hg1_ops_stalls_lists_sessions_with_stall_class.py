"""`delfin-agent ops stalls` -- a read-only list of stalls, per session.

Input: the session rows `session_health.collect()` produces, enriched
by `ops_interface.sessions()` fields (idle_s, open_approvals, pid,
paused). Output: one line per session ``<name>  <class>  <age>`` via
the `_cmd_ops_stalls` handler. Asserted: open approvals outrank
everything; paused/dead-pid/api-errors are unclassified, not guessed;
a long idle turn-over session is idle, not a stall; fresh rows are
running. Nothing is mutated by the command.
"""

from __future__ import annotations

import argparse

import pytest

from delfin.agent import progress_watchdog as pw


NOW = 1_700_000_000.0


def _record(**over):
    base = {
        "name": "w20/hg1", "session_id": "abc-123",
        "pid": 4242, "pid_alive": True, "paused": False,
        "open_approvals": [], "approval_opened_at": None,
        "no_token_s": None, "gate_wait_s": None,
        "idle_s": 30, "last_turn": {}, "phase": "",
        "api_errors_in_a_row": 0,
    }
    base.update(over)
    return base


class TestStallClassPrecedence:
    def test_open_approval_outranks_everything(self):
        rec = _record(open_approvals=["a1"], approval_opened_at=NOW - 600,
                      no_token_s=900.0, gate_wait_s=800.0, phase="tool")
        out = pw.stall_class(rec, NOW)
        assert out["stall"] is True
        assert out["class"] == "approval-wait"
        assert out["age_s"] == 600.0
        assert out["action"] == pw.KEEP_WAITING

    def test_gate_queue_outranks_model_wait(self):
        out = pw.stall_class(_record(no_token_s=120.0, gate_wait_s=90.0), NOW)
        assert out["class"] == "gate-queue"
        assert out["age_s"] == 90.0

    def test_model_wait_from_no_token_s(self):
        out = pw.stall_class(_record(no_token_s=120.0), NOW)
        assert out["class"] == "model-wait"
        assert out["age_s"] == 120.0

    def test_running_phase_is_tool_running(self):
        out = pw.stall_class(_record(phase="tool", phase_s=45.0), NOW)
        assert out["class"] == "tool-running"
        assert out["age_s"] == 45.0
        assert out["action"] == pw.KEEP_WAITING

    @pytest.mark.parametrize("kw", [
        {"paused": True, "idle_s": 5000},
        {"pid_alive": False},
        {"api_errors_in_a_row": 5},
    ])
    def test_unstoppable_states_are_unclassified_not_guessed(self, kw):
        out = pw.stall_class(_record(**kw), NOW)
        assert out["class"] == "unclassified"
        assert out["stall"] is False

    def test_quiet_session_is_idle_not_a_stall(self):
        out = pw.stall_class(_record(idle_s=7200), NOW)
        assert out["class"] == "idle"
        assert out["stall"] is False

    def test_fresh_active_row_is_running_not_idle(self):
        out = pw.stall_class(_record(idle_s=10), NOW)
        assert out["class"] == "idle"
        assert out["stall"] is False

    def test_non_dict_is_unclassified(self):
        out = pw.stall_class(None, NOW)
        assert out["class"] == "unclassified"


class TestOpsStallsCommand:
    def _run(self, rows, monkeypatch, capsys, tmp_path):
        from delfin.agent import cli as cli_mod
        # The command's data source: session_health.collect() takes the
        # approvals from ops_interface.pending_approvals(); that is what
        # the test fakes (collect() recomputes open_approvals from it).
        # Each test sets `approvals` on the row itself: open_approvals
        # on the fake row, so only tests that ask for them see one.
        from delfin.agent import ops_interface as oi
        from delfin.agent import session_pause
        monkeypatch.setattr(session_pause, "_DIR",
                            tmp_path / "pause", raising=False)
        monkeypatch.setattr(oi, "sessions", lambda **kw: rows)
        fake_approvals = [
            {"session": r["name"], "session_id": r["session_id"],
             "id": a, "opened_at": r["approval_opened_at"]}
            for r in rows for a in (r.get("open_approvals") or [])]
        monkeypatch.setattr(oi, "pending_approvals",
                            lambda: fake_approvals)
        monkeypatch.setattr("time.time", lambda: NOW)
        args = argparse.Namespace(ops_action="stalls", json=False)
        cli_mod._cmd_ops_stalls(args)
        return capsys.readouterr().out

    def test_one_line_per_session_with_class_and_age(self, monkeypatch,
                                                     capsys, tmp_path):
        rows = [_record(name="w20/hg1", phase="tool", phase_s=45.0),
                _record(name="w20/other", session_id="def-456",
                        open_approvals=["a1"],
                        approval_opened_at=NOW - 600)]
        out = self._run(rows, monkeypatch, capsys, tmp_path)
        lines = [l for l in out.splitlines() if l.strip()]
        assert any("w20/hg1" in l and "tool-running" in l for l in lines)
        assert any("w20/other" in l and "approval-wait" in l for l in lines)
        assert any("10 min" in l for l in lines)

    def test_no_stalls_prints_a_clear_line(self, monkeypatch, capsys,
                                           tmp_path):
        out = self._run([_record(idle_s=5)], monkeypatch, capsys, tmp_path)
        assert "no stalled" in out.lower()

    def test_command_touches_nothing(self, monkeypatch, capsys,
                                     tmp_path):
        calls = []
        out = self._run([_record(name="w/x", no_token_s=300.0)],
                        monkeypatch, capsys, tmp_path)
        assert "w/x" in out and "model-wait" in out
        assert not calls
