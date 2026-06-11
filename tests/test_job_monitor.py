"""Tests for the proactive SLURM job monitor.

Contract: opt-in (disabled by default → daemon refuses to run), the
watch loop is LLM-free, diagnosis is a single read-only turn that can be
disabled separately, findings round-trip for the dashboard banner, and
the PID lock prevents double daemons.  No SLURM and no LLM needed —
squeue/sacct and the engine are injectable.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import job_monitor as jm


_ON = {"agent": {"job_monitor": {"enabled": True}}}
_OFF = {"agent": {"job_monitor": {"enabled": False}}}


# ---------------------------------------------------------------------------
# Opt-in: default OFF, daemon refuses while disabled
# ---------------------------------------------------------------------------

def test_default_settings_disabled():
    from delfin.user_settings import DEFAULT_SETTINGS
    cfg = DEFAULT_SETTINGS["agent"]["job_monitor"]
    assert cfg["enabled"] is False          # token safety: explicit opt-in
    assert cfg["auto_diagnose"] is True
    assert cfg["interval_s"] >= 60
    # dedicated diagnosis provider/model configurable, empty = agent default
    assert cfg["provider"] == "" and cfg["model"] == ""


def test_run_loop_exits_immediately_when_disabled():
    assert jm.run_loop(settings=_OFF, max_iterations=99) == 0


def test_monitor_settings_reads_dedicated_model():
    cfg = jm.monitor_settings({"agent": {"job_monitor": {
        "enabled": True, "provider": "kit", "model": "azure.gpt-5-nano"}}})
    assert cfg["provider"] == "kit"
    assert cfg["model"] == "azure.gpt-5-nano"


# ---------------------------------------------------------------------------
# Watch list + LLM-free failure detection
# ---------------------------------------------------------------------------

def _fake_run(squeue: str = "", sacct: str = ""):
    def run(cmd):
        if cmd[0] == "squeue":
            return squeue
        if cmd[0] == "sacct":
            return sacct
        return ""
    return run


def test_watch_add_remove_roundtrip(tmp_path):
    wp = tmp_path / "watched.json"
    jm.add_watch("111", "/calc/a", path=wp)
    jm.add_watch("222", path=wp)
    assert set(jm.load_watched(wp)["jobs"]) == {"111", "222"}
    jm.remove_watch("111", path=wp)
    assert set(jm.load_watched(wp)["jobs"]) == {"222"}


def test_failed_job_produces_one_finding(tmp_path):
    wp = tmp_path / "watched.json"
    jm.add_watch("4976064", str(tmp_path), path=wp)
    run = _fake_run(sacct="4976064  FAILED\n")
    findings = jm.check_once(path=wp, run_fn=run)
    assert len(findings) == 1
    assert findings[0].job_id == "4976064"
    assert findings[0].state == "FAILED"
    # No duplicate alarm on the next poll (state already recorded).
    assert jm.check_once(path=wp, run_fn=run) == []


def test_running_and_completed_jobs_do_not_alert(tmp_path):
    wp = tmp_path / "watched.json"
    jm.add_watch("1", path=wp)
    jm.add_watch("2", path=wp)
    run = _fake_run(squeue="1 RUNNING\n", sacct="2  COMPLETED\n")
    assert jm.check_once(path=wp, run_fn=run) == []


def test_error_signature_scan_finds_jeromes_venv_case(tmp_path):
    err = tmp_path / "delfin_4976064.err"
    err.write_text(
        "slurm_script: line 287: /scratch/.../delfin_venv/bin/activate: "
        "No such file or directory\n"
    )
    labels = jm.scan_error_signatures(str(tmp_path))
    assert "venv-activation-failed" in labels


def test_signature_scan_handles_missing_folder():
    assert jm.scan_error_signatures("/does/not/exist") == []


# ---------------------------------------------------------------------------
# Findings log round-trip (dashboard banner source)
# ---------------------------------------------------------------------------

def test_findings_roundtrip_and_since_filter(tmp_path):
    fp = tmp_path / "findings.jsonl"
    f1 = jm.Finding("1", "/a", "FAILED", ts=100.0)
    f2 = jm.Finding("2", "/b", "TIMEOUT", ts=200.0)
    jm.record_finding(f1, path=fp)
    jm.record_finding(f2, path=fp)
    assert len(jm.load_findings(path=fp)) == 2
    newer = jm.load_findings(since=150.0, path=fp)
    assert [r["job_id"] for r in newer] == ["2"]


# ---------------------------------------------------------------------------
# Diagnosis: injectable engine, read-only, separately disableable
# ---------------------------------------------------------------------------

def test_auto_diagnose_off_spends_zero_tokens():
    calls = []
    f = jm.diagnose_finding(
        jm.Finding("1", "/x", "FAILED"),
        settings={"agent": {"job_monitor": {"enabled": True,
                                             "auto_diagnose": False}}},
        engine_factory=lambda *a, **k: calls.append(1),
    )
    assert calls == []                      # engine never built
    assert "no LLM" in f.summary


def test_diagnosis_uses_engine_and_records_summary(monkeypatch):
    class _FakeEngine:
        def stream_response(self, *, user_message, on_token, **kw):
            assert "READ-ONLY" in user_message
            on_token("Ursache: venv-Tarball nicht entpackt. Fix: …")

    saved = {}
    monkeypatch.setattr(
        "delfin.agent.session_store.save_session",
        lambda sid, **kw: saved.update({"sid": sid, **kw}),
    )
    f = jm.diagnose_finding(
        jm.Finding("4976064", "/calc/test", "FAILED",
                   signatures=["venv-activation-failed"]),
        settings=_ON,
        engine_factory=lambda folder, settings: _FakeEngine(),
    )
    assert "venv-Tarball" in f.summary
    assert f.diagnosis_session.startswith("monitor-4976064-")
    assert saved["sid"] == f.diagnosis_session
    assert "🚨" in saved.get("title", "")


def test_diagnosis_attaches_fix_assessment_zero_tokens(monkeypatch, tmp_path):
    # An OOM kill with a real submit script -> finding carries a prepared
    # one-click resource fix, even with auto_diagnose OFF (LLM-free).
    (tmp_path / "submit.sh").write_text(
        "#!/bin/bash\n#SBATCH --mem=16G\n#SBATCH --time=4:00:00\norca x.inp\n"
    )
    f = jm.diagnose_finding(
        jm.Finding("555", str(tmp_path), "FAILED",
                   signatures=["out-of-memory"]),
        settings={"agent": {"job_monitor": {"enabled": True,
                                             "auto_diagnose": False}}},
        engine_factory=lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("engine must not be built when auto_diagnose off")),
    )
    assert f.fix.get("fix_class") == "slurm_resource"
    assert f.fix.get("one_click") is True
    assert "--mem=24G" in f.fix["proposal"]["new_line"]


def test_diagnosis_failure_is_contained(monkeypatch):
    # Isolate the session store — without this the test writes a REAL
    # "🚨 Job 9 failed" session into the user's ~/.delfin/agent_sessions
    # (which confused a live user once).
    monkeypatch.setattr(
        "delfin.agent.session_store.save_session",
        lambda sid, **kw: None,
    )

    def _boom(folder, settings):
        raise RuntimeError("no api key")
    f = jm.diagnose_finding(jm.Finding("9", "/x", "FAILED"),
                            settings=_ON, engine_factory=_boom)
    assert "diagnosis failed" in f.summary   # never raises


# ---------------------------------------------------------------------------
# PID lock: single instance
# ---------------------------------------------------------------------------

def test_pid_lock_blocks_second_instance(tmp_path):
    pp = tmp_path / "monitor.pid"
    assert jm.acquire_pid_lock(pp) is True          # we own it
    assert pp.read_text().strip() == str(os.getpid())
    assert jm.acquire_pid_lock(pp) is False         # same live pid → blocked
    jm.release_pid_lock(pp)
    assert not pp.exists()


def test_stale_pid_is_reclaimed(tmp_path):
    pp = tmp_path / "monitor.pid"
    pp.write_text("999999999")                       # definitely dead
    assert jm.acquire_pid_lock(pp) is True
    jm.release_pid_lock(pp)
