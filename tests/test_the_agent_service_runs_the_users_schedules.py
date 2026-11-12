"""The agent service runs the user's schedules, each run a session.

* A due schedule of this host becomes a session in its folder, given its
  task as the first prompt; how the run ended, its answer, cost and
  session are recorded, and the schedule is free for its next run.
* A run that fails is recorded as failed; one held by an approval nobody
  answers is reported as waiting, then stopped at its time ceiling.
* A schedule of another host is not run here.
* The service stays up while this host has a schedule to run, and goes
  after its idle time once there is none.
* A scheduled session never starts in bypass, whatever the settings say,
  and carries its cost and time ceilings to the worker.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from delfin.agent import schedules
from delfin.agent.service import rundir
from tui_service import ServiceProcess


def _service(tmp_path, *, idle_exit: float = 120.0) -> ServiceProcess:
    svc = ServiceProcess(tmp_path, start=False, idle_exit=idle_exit)
    svc.env.update(DELFIN_SCHEDULES_FAKE="0.003", DELFIN_SCHEDULES_POLL_S="0.2",
                   DELFIN_SCHEDULES_MINUTE_S="2")
    return svc


@pytest.fixture
def store(monkeypatch):
    """Point this process at the service's store: set per test below."""
    def at(svc: ServiceProcess) -> Path:
        from delfin.agent import trust_store
        path = svc.scratch / ".delfin" / "user_schedules"
        monkeypatch.setattr(schedules, "_DIR", path)
        # The service's account key: it checks the confirmations made here.
        monkeypatch.setattr(trust_store, "store_dir",
                            lambda: svc.scratch / ".delfin" / "trust_store")
        trust_store.forget_cached_keys()
        return path
    return at


def _add(svc, text, **kw):
    kw.setdefault("confirm", True)
    return schedules.add(text, workspace=str(svc.workspace()), created_by="test", **kw)


def _wait_run(entry_id: int, statuses: set[str], timeout: float = 60) -> dict:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        for row in schedules.runs(limit=500):
            if row["id"] == entry_id and row["status"] in statuses:
                return row
        time.sleep(0.1)
    raise AssertionError(f"no {statuses} run of #{entry_id}: {schedules.runs()}")


def test_a_due_schedule_runs_as_a_session_and_is_recorded(tmp_path, store):
    svc = _service(tmp_path)
    store(svc)
    try:
        entry = _add(svc, "alle 1h Lies die neuen Ergebnisse")
        svc.start()
        schedules.request_run(entry["id"])
        run = _wait_run(entry["id"], {"done"})
        assert run["session"] and run["text"]
        after = schedules.get(entry["id"])
        assert after["running"] == "" and after["last_status"] == "done"
        assert after["last_session"] == run["session"]
        with svc.client() as c:
            listed = c.request("sessions.list")["sessions"]
        row = next(s for s in listed if s["sessionId"] == run["session"])
        assert f"#{entry['id']}" in (row["title"] or "")
        assert Path(row["workspace"]).resolve() == svc.workspace().resolve()
    finally:
        svc.stop()


def test_a_failing_run_is_recorded_as_failed(tmp_path, store):
    svc = _service(tmp_path)
    store(svc)
    try:
        entry = _add(svc, "alle 1h fail please")
        svc.start()
        schedules.request_run(entry["id"])
        _wait_run(entry["id"], {"failed"})
        assert schedules.get(entry["id"])["running"] == ""
    finally:
        svc.stop()


def test_a_run_held_by_an_approval_is_reported_then_stopped(tmp_path, store):
    svc = _service(tmp_path)
    store(svc)
    try:
        entry = _add(svc, "alle 1h bash it", max_minutes=1)
        svc.start()
        schedules.request_run(entry["id"])
        waiting = _wait_run(entry["id"], {"waiting"})
        assert waiting["session"]
        _wait_run(entry["id"], {"timeout"})
        assert schedules.get(entry["id"])["running"] == ""
    finally:
        svc.stop()


def test_another_hosts_schedule_is_not_run_here(tmp_path, store):
    svc = _service(tmp_path)
    store(svc)
    try:
        entry = _add(svc, "alle 1h Lies", host="another-node")
        svc.start()
        schedules.request_run(entry["id"])
        time.sleep(2.0)
        assert schedules.runs() == []
        assert schedules.get(entry["id"])["run_requested"] is True
    finally:
        svc.stop()


def test_the_service_stays_up_for_a_schedule_and_goes_without_one(tmp_path, store):
    svc = _service(tmp_path, idle_exit=1.0)
    store(svc)
    try:
        entry = _add(svc, "täglich 07:30 Lies")
        svc.start()
        time.sleep(4.0)
        assert svc.proc.poll() is None, svc.log()
        schedules.delete(entry["id"])
        deadline = time.monotonic() + 30
        while svc.proc.poll() is None and time.monotonic() < deadline:
            time.sleep(0.2)
        assert svc.proc.poll() is not None
    finally:
        svc.stop()


def test_a_scheduled_session_is_never_bypass_and_carries_its_ceilings(tmp_path):
    from delfin.agent.service import schedule_runner
    entry = {"id": 4, "task": "x", "when": "alle 1h", "workspace": str(tmp_path),
             "budget_run_usd": 0.5, "max_minutes": 20}
    bypass = {"agent": {"scheduler": {"permission_mode": "bypassPermissions",
                                      "unattended": True}}}
    params = schedule_runner.session_params(entry, bypass)
    assert params["permissionMode"] in ("plan", "default", "acceptEdits")
    assert params["unattended"] is True and params["attach"] is False
    assert params["maxBudgetUsd"] == 0.5 and params["maxRunSeconds"] == 20 * 60
    assert params["cwd"] == str(tmp_path)
    assert schedule_runner.session_params(entry, {})["permissionMode"] == "default"


def test_the_worker_takes_the_ceilings_onto_its_session(tmp_path):
    from delfin.agent.service import worker
    args = worker.parse_args(["--cwd", str(tmp_path), "--unattended",
                              "--max-budget-usd", "0.5", "--max-run-seconds", "1200"])
    ns = worker.chat_namespace(args, tmp_path)
    assert ns.unattended is True
    assert ns.max_budget_usd == 0.5 and ns.max_run_seconds == 1200


def test_the_host_name_is_the_services(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "s")
    ws = tmp_path / "w"
    ws.mkdir()
    entry = schedules.add("alle 1h x", workspace=str(ws), created_by="t", confirm=True)
    assert entry["host"] == rundir.hostname()
