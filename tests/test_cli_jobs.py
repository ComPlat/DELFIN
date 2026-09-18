"""The cluster's jobs, reachable from the terminal and from a session.

``job_monitor`` already knew everything about SLURM jobs — the watch
lists, the three-state scheduler query, the daemon with its PID lock —
but only the dashboard had a window into it. The terminal had no
subcommand and a session had no slash command (background bash jobs had
``/bash``; SLURM jobs had nothing).

These tests pin the MISSING SURFACE, not a second mechanism: everything
funnels through ``query_job_states_detailed`` and the existing watch
files. No squeue binary is needed — every scheduler touch is injected,
and a host without squeue must say so in one line, not a traceback.
"""

from __future__ import annotations

import json
import subprocess

import pytest

from delfin.agent import cli_jobs


# ---------------------------------------------------------------------------
# collect: one row per watched job, running/queued and recently finished
# ---------------------------------------------------------------------------

def _no_scheduler(cmd):
    """run_fn that answers as a host with no squeue/sacct at all."""
    raise FileNotFoundError(cmd[0])


def test_rows_from_real_squeue_output(tmp_path):
    """Driven against real squeue/sacct output shapes, injected."""
    agent_ws = tmp_path / "ws"
    agent_ws.mkdir()
    # shared watch list (dashboard-written) with one running job
    shared = tmp_path / "shared_watched.json"
    shared.write_text(json.dumps({"jobs": {
        "4976064": {"folder": "/tmp/run1", "last_state": "RUNNING"},
    }}), encoding="utf-8")
    # agent watch list with one queued job and one finished
    agent_file = agent_ws / ".delfin" / "agent_watched_jobs.json"
    agent_file.parent.mkdir(parents=True)
    agent_file.write_text(json.dumps({"jobs": {
        "111": {"kind": "slurm", "description": "prod run",
                "added_at": 0.0, "last_state": ""},
        "222": {"kind": "slurm", "description": "short opt",
                "added_at": 0.0, "last_state": "COMPLETED"},
    }}), encoding="utf-8")

    def run_fn(cmd):
        if cmd[0] == "squeue":
            # real shape: -h -o "%i %T"; 222 left the queue already
            return "111 PENDING\n"
        if cmd[0] == "sacct":
            return "222 COMPLETED\n"
        return None

    rows = cli_jobs.collect_job_rows(
        run_fn=run_fn, shared_watched_path=shared,
        agent_workspaces=[str(agent_ws)])
    by_id = {r["job_id"]: r for r in rows}
    # 4976064 is only on the shared list, which the injected scheduler
    # does not know — "not known" must be said, not rendered as quiet.
    assert by_id["4976064"]["state"] == "UNKNOWN"
    assert by_id["111"]["state"] == "PENDING"
    assert by_id["222"]["state"] == "COMPLETED"
    assert by_id["222"]["workspace"].endswith("ws")
    assert by_id["4976064"]["workspace"] == ""  # shared list has no workspace


def test_running_first_then_queued_then_finished(tmp_path):
    shared = tmp_path / "shared_watched.json"
    shared.write_text(json.dumps({"jobs": {
        "3": {"folder": "", "last_state": "COMPLETED"},
        "1": {"folder": "", "last_state": "RUNNING"},
        "2": {"folder": "", "last_state": "PENDING"},
    }}), encoding="utf-8")

    def run_fn(cmd):
        if cmd[0] == "squeue":
            return "1 RUNNING\n2 PENDING\n3 COMPLETED\n"
        if cmd[0] == "sacct":
            return "3 COMPLETED\n"
        return None

    rows = cli_jobs.collect_job_rows(
        run_fn=run_fn, shared_watched_path=shared)
    assert [r["job_id"] for r in rows] == ["1", "2", "3"]


def test_unknown_job_id_is_said_not_hidden(tmp_path):
    """Aged-out ids must appear as unknown — a quiet table is a lie."""
    shared = tmp_path / "shared_watched.json"
    shared.write_text(json.dumps({"jobs": {
        "999": {"folder": "", "last_state": "RUNNING"},
    }}), encoding="utf-8")

    def run_fn(cmd):
        if cmd[0] == "squeue":
            return None       # squeue -j fails once an id left the queue
        if cmd[0] == "sacct":
            return ""         # accounting does not know it either
        return None

    rows = cli_jobs.collect_job_rows(
        run_fn=run_fn, shared_watched_path=shared)
    assert len(rows) == 1
    assert rows[0]["state"] == "UNKNOWN"


def test_bash_and_ci_watches_are_skipped(tmp_path):
    """This surface is the cluster's jobs; /bash already owns bash jobs."""
    ws = tmp_path / "ws"
    agent_file = ws / ".delfin" / "agent_watched_jobs.json"
    agent_file.parent.mkdir(parents=True)
    agent_file.write_text(json.dumps({"jobs": {
        "a1b2c3d4": {"kind": "bash", "description": "xtb scan",
                     "added_at": 0.0},
        "ci:owner/repo@abc123": {"kind": "ci", "description": "CI",
                                 "added_at": 0.0},
    }}), encoding="utf-8")
    rows = cli_jobs.collect_job_rows(
        run_fn=_no_scheduler, shared_watched_path=tmp_path / "none.json",
        agent_workspaces=[str(ws)])
    assert rows == []


def test_no_squeue_at_all_is_one_line_not_a_traceback(tmp_path):
    """Universal: a host without squeue says so and carries on."""
    shared = tmp_path / "shared_watched.json"
    shared.write_text(json.dumps({"jobs": {
        "5": {"folder": "/tmp/x", "last_state": "RUNNING"},
    }}), encoding="utf-8")
    # collect raises the named condition; the caller prints one line.
    with pytest.raises(cli_jobs.SchedulerUnavailable):
        cli_jobs.collect_job_rows(
            run_fn=_no_scheduler, shared_watched_path=shared)
    text = cli_jobs.scheduler_note(1, FileNotFoundError("squeue"))
    assert "squeue" in text    # the host says what it lacks
    assert "Traceback" not in text
    # and an empty watch list on the same host is empty, not degraded
    assert cli_jobs.collect_job_rows(
        run_fn=_no_scheduler,
        shared_watched_path=tmp_path / "none.json") == []


def test_empty_lists_render_a_hint():
    text = cli_jobs.render_jobs([])
    assert "no jobs" in text.lower()


def test_render_shows_id_state_elapsed_and_workspace(tmp_path):
    shared = tmp_path / "shared_watched.json"
    shared.write_text(json.dumps({"jobs": {
        "10": {"folder": "/tmp/run", "last_state": "RUNNING",
               "added_at": 0.0},
    }}), encoding="utf-8")

    def run_fn(cmd):
        if cmd[0] == "squeue":
            return "10 RUNNING\n"
        return "" if cmd[0] == "sacct" else None

    rows = cli_jobs.collect_job_rows(
        run_fn=run_fn, shared_watched_path=shared, now=lambda: 120.0)
    text = cli_jobs.render_jobs(rows)
    assert "10" in text and "RUNNING" in text
    assert "2m" in text          # elapsed since added_at=0
    assert "/tmp/run" in text    # the workspace each job runs in


# ---------------------------------------------------------------------------
# watch: daemon on/off/status, and the disabled-default honesty
# ---------------------------------------------------------------------------

def test_watch_status_reports_disabled_default(tmp_path, monkeypatch):
    monkeypatch.setattr(cli_jobs, "_settings_path", lambda: tmp_path / "s.json")
    text = cli_jobs.watch_report(pid_path=tmp_path / "no.pid",
                                 watched_path=tmp_path / "none.json")
    assert "disabled" in text
    assert "agent.job_monitor.enabled" in text


def test_watch_report_running_daemon_and_watched_count(tmp_path, monkeypatch):
    monkeypatch.setattr(cli_jobs, "_settings_path", lambda: tmp_path / "s.json")
    (tmp_path / "s.json").write_text(json.dumps(
        {"agent": {"job_monitor": {"enabled": True}}}), encoding="utf-8")
    me = __import__("os").getpid()
    (tmp_path / "pid").write_text(str(me), encoding="utf-8")
    (tmp_path / "w.json").write_text(json.dumps(
        {"jobs": {"1": {}, "2": {}}}), encoding="utf-8")
    text = cli_jobs.watch_report(pid_path=tmp_path / "pid",
                                 watched_path=tmp_path / "w.json")
    assert "running" in text and str(me) in text
    assert "2" in text


def test_watch_on_persists_the_setting(tmp_path, monkeypatch):
    monkeypatch.setattr(cli_jobs, "_settings_path", lambda: tmp_path / "s.json")
    out = cli_jobs.watch_set("on", pid_path=tmp_path / "no.pid")
    assert "enabled" in out
    saved = json.loads((tmp_path / "s.json").read_text(encoding="utf-8"))
    assert saved["agent"]["job_monitor"]["enabled"] is True
    # off again
    cli_jobs.watch_set("off", pid_path=tmp_path / "no.pid")
    saved = json.loads((tmp_path / "s.json").read_text(encoding="utf-8"))
    assert saved["agent"]["job_monitor"]["enabled"] is False


def test_watch_on_off_do_not_launch_or_kill_a_daemon(tmp_path, monkeypatch):
    """Setting the flag must not spawn anything — only the flag."""
    monkeypatch.setattr(cli_jobs, "_settings_path", lambda: tmp_path / "s.json")
    spawned = []
    monkeypatch.setattr(subprocess, "Popen",
                        lambda *a, **k: spawned.append(a))
    cli_jobs.watch_set("on", pid_path=tmp_path / "no.pid")
    cli_jobs.watch_set("off", pid_path=tmp_path / "no.pid")
    assert spawned == []


# ---------------------------------------------------------------------------
# wiring: the subcommand and the slash command exist
# ---------------------------------------------------------------------------

def test_cli_registers_jobs_subcommands():
    from delfin.agent import cli as agent_cli
    parser = agent_cli.build_parser()
    assert "jobs" in agent_cli._subcommand_names(parser)


def test_cli_help_lists_jobs():
    from delfin.agent import help_gen, repl_commands
    text = help_gen.generate_help(repl_commands.palette_rows())
    assert "/jobs" in text


def test_repl_has_a_jobs_command():
    from delfin.agent import repl_commands
    assert "/jobs" in repl_commands.BUILTINS


def test_dashboard_jobs_prefix_still_routes_to_builtin():
    from delfin.dashboard import tab_agent
    assert tab_agent.slash_command_routes_to_builtin("/jobs")
    assert tab_agent.slash_command_routes_to_builtin("/jobs watch status")
