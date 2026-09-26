"""Every budget breach reaches the audit log (phase 5).

A breach the round report cannot count is a breach that did not
happen, observably. record_breach() writes one `process_budget` audit
record per breach — session, profile, limit, value, ceiling — and the
guards call it on their own, so foreground, background, hook and test
breaches all land there without each caller doing its own bookkeeping.
"""
import json
import time

import pytest

from delfin.agent import audit_log, bash_jobs, process_budget


def _fanout(n: int) -> str:
    return f"for i in $(seq 1 {n}); do sleep 15 & done; sleep 15"


def _shrink(monkeypatch, name, max_processes):
    monkeypatch.setattr(
        process_budget, "limits_for_profile",
        lambda n: process_budget.BudgetLimits(
            max_processes=max_processes if n == name else 4,
            max_rss_mb=1e9, max_cpu_seconds=1e9))


def test_record_breach_writes_a_process_budget_event(tmp_path):
    breach = process_budget.Breach(
        limit="max_processes", value=12, ceiling=4,
        message="Process budget exceeded: max_processes ...")
    process_budget.record_breach(
        breach, profile="background", session_id="sess-1",
        log_path=tmp_path / "audit.log")
    lines = (tmp_path / "audit.log").read_text().splitlines()
    assert len(lines) == 1
    rec = json.loads(lines[0])
    assert rec["tool"] == "process_budget"
    assert rec["decision"] == "block"
    assert rec["session_id"] == "sess-1"
    assert rec["profile"] == "background"
    assert rec["limit"] == "max_processes"
    assert rec["value"] == 12
    assert rec["ceiling"] == 4


def test_record_breach_never_raises(tmp_path):
    breach = process_budget.Breach(
        limit="max_rss_mb", value=99, ceiling=8, message="x")
    # A log path that cannot be written must not propagate.
    process_budget.record_breach(
        breach, profile="hook", log_path=tmp_path / "no" / "dir" / "x.log")
    # And no exception — audit must not break the guard.


def test_guard_poll_once_records_the_breach(monkeypatch, tmp_path):
    """A breach found by poll_once lands in the audit log on its own —
    the foreground/hook/tests call path."""
    _shrink(monkeypatch, "hook", 4)
    monkeypatch.setattr(audit_log, "_default_log_path",
                        lambda: tmp_path / "audit.log")
    import subprocess as sp
    proc = sp.Popen(
        ["/bin/bash", "-c", _fanout(12)], cwd=str(tmp_path),
        stdout=sp.DEVNULL, stderr=sp.DEVNULL,
        stdin=sp.DEVNULL, start_new_session=True)
    log = tmp_path / "audit.log"
    try:
        guard = process_budget.BudgetGuard(
            proc.pid, poll_s=0.2, profile="hook")
        deadline = time.monotonic() + 10
        breach = None
        while time.monotonic() < deadline:
            breach = guard.poll_once()
            if breach is not None:
                break
            time.sleep(0.1)
        assert breach is not None
        recs = audit_log.read_last_n(10, log_path=log)
        assert any(r.get("tool") == "process_budget"
                   and r.get("profile") == "hook"
                   and r.get("limit") == "max_processes" for r in recs), recs
    finally:
        if proc.poll() is None:
            proc.kill()
        proc.wait()


def test_background_breach_reaches_the_audit_log(monkeypatch, tmp_path):
    """The bash_background path: a breach recorded on the job also
    writes the audit event — through the public registry start."""
    _shrink(monkeypatch, "background", 4)
    monkeypatch.setattr(audit_log, "_default_log_path",
                        lambda: tmp_path / "audit.log")
    reg = bash_jobs._Registry()
    job = reg.start(_fanout(12), cwd=str(tmp_path),
                    timeout_s=60, workspace=str(tmp_path),
                    session_id="sess-bg")
    deadline = time.monotonic() + 25
    while time.monotonic() < deadline and job.poll() is None:
        time.sleep(0.2)
    assert job.poll() is not None
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and job.budget_breach is None:
        time.sleep(0.2)
    assert job.budget_breach is not None
    recs = audit_log.read_last_n(10, log_path=tmp_path / "audit.log")
    assert any(r.get("tool") == "process_budget"
               and r.get("profile") == "background"
               and r.get("session_id") == "sess-bg"
               for r in recs), recs
