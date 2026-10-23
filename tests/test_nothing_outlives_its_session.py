"""Nothing an agent session started outlives it.

A dashboard server started on 2026-09-11 was still running on 2026-09-16:
its launcher had gone, no window had ever opened a kernel, and the server
only stopped itself when a kernel ended. The agent's background shells and
MCP servers run in process groups of their own and did not end with the
kernel either. The rule now: without a kept session everything ends with
the last window; with one, once it is no longer kept.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import bash_jobs as bj
from delfin.dashboard import resume_server as R


def _wait_ended(job, seconds=10.0):
    deadline = time.time() + seconds
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.1)
    return job.poll() is not None


def test_a_server_nobody_opens_ends_by_itself(monkeypatch):
    monkeypatch.delenv(R.STAY_UP_ENV, raising=False)
    assert (f"--ServerApp.shutdown_no_activity_timeout="
            f"{int(R.NEVER_CONNECTED_SECONDS)}") in R.cull_config_args()
    monkeypatch.setenv(R.STAY_UP_ENV, "1")
    assert not any("shutdown_no_activity_timeout" in a
                   for a in R.cull_config_args())


def test_closing_a_session_stops_its_own_shells_only(tmp_path):
    registry = bj.get_registry()
    mine = registry.start(command="sleep 60", cwd=str(tmp_path),
                          workspace=str(tmp_path), session_id="S")
    theirs = registry.start(command="sleep 60", cwd=str(tmp_path),
                            workspace=str(tmp_path), session_id="T")
    try:
        assert registry.stop_running(session_id="S") == [mine.job_id]
        assert _wait_ended(mine)
        assert theirs.poll() is None
    finally:
        registry.kill(theirs.job_id)
        registry.kill(mine.job_id)


def test_a_closed_tab_ends_the_shells_it_started(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    import atexit

    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    registered = []
    monkeypatch.setattr(atexit, "register",
                        lambda fn, *a, **kw: registered.append((fn, a, kw)))
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    _tab, refs = tab_agent.create_tab(ctx)

    refs["state"]["active_session_id"] = "S"
    job = bj.get_registry().start(command="sleep 60", cwd=str(tmp_path),
                                  workspace=str(tmp_path), session_id="S")
    try:
        refs["shutdown"]()
        assert _wait_ended(job)
        refs["shutdown"]()                      # a second close is harmless
    finally:
        bj.get_registry().kill(job.job_id)

    # The kernel ending closes the tab and stops what the process started.
    assert any(fn is refs["shutdown"] and kw == {"save": False}
               for fn, _a, kw in registered)


def test_the_kernel_ending_stops_every_shell_and_mcp_server(tmp_path,
                                                            monkeypatch):
    from delfin.agent import mcp_client
    from delfin.dashboard import tab_agent

    reset = []
    monkeypatch.setattr(mcp_client, "reset_registry",
                        lambda workspace=None: reset.append(workspace))
    job = bj.get_registry().start(command="sleep 60", cwd=str(tmp_path),
                                  workspace=str(tmp_path))
    try:
        tab_agent._stop_what_this_process_started()
        assert _wait_ended(job)
        assert reset == [None], "every MCP registry is stopped"
    finally:
        bj.get_registry().kill(job.job_id)
