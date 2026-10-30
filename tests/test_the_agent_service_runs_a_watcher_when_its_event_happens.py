"""The agent service looks at the watchers of its host and runs them when
their event happens.

* A file that appears starts the watcher's task as a session, the event
  named in its prompt; it fires once.
* A watcher without a task only reports the event (no session).
* A watcher that cannot work here (a link out of its folder) is paused,
  and that is reported.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import pytest

from delfin.agent import schedules
from tui_service import ServiceProcess


@pytest.fixture
def svc(tmp_path, monkeypatch):
    from delfin.agent import trust_store
    service = ServiceProcess(tmp_path, start=False)
    service.env.update(DELFIN_SCHEDULES_FAKE="0.003", DELFIN_SCHEDULES_POLL_S="0.2",
                       DELFIN_SCHEDULES_MINUTE_S="2")
    monkeypatch.setattr(schedules, "_DIR", service.scratch / ".delfin" / "user_schedules")
    monkeypatch.setattr(trust_store, "store_dir",
                        lambda: service.scratch / ".delfin" / "trust_store")
    trust_store.forget_cached_keys()
    yield service
    service.stop()


def _wait_run(entry_id: int, statuses: set[str], timeout: float = 60) -> dict:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        for row in schedules.runs(limit=500):
            if row["id"] == entry_id and row["status"] in statuses:
                return row
        time.sleep(0.1)
    raise AssertionError(f"no {statuses} run of #{entry_id}: {schedules.runs()}")


def _fast(entry_id: int) -> None:
    """Look now rather than in a minute."""
    data = __import__("json").loads(schedules._file().read_text())
    data["entries"][str(entry_id)]["next_fire_at"] = 0
    schedules._file().write_text(__import__("json").dumps(data))


def test_a_file_that_appears_runs_the_task_once(svc):
    ws = svc.workspace()
    entry = schedules.add("wenn datei result.txt erscheint Lies das Ergebnis",
                          workspace=str(ws), created_by="test", confirm=True)
    svc.start()
    _fast(entry["id"])
    time.sleep(1.0)
    assert schedules.runs() == []
    (ws / "result.txt").write_text("42")
    _fast(entry["id"])
    run = _wait_run(entry["id"], {"done"})
    assert run["session"]
    after = schedules.get(entry["id"])
    assert after["watch_done"] is True and after["trigger_note"] == ""
    with svc.client() as c:
        events = c.request("session.replay", {"sessionId": run["session"]})["events"]
    said = [e.get("text", "") for e in events if e.get("type") == "user.message"]
    assert any("result.txt appeared" in s for s in said), said


def test_a_watcher_without_a_task_only_reports(svc):
    ws = svc.workspace()
    entry = schedules.add("wenn datei done.flag erscheint", workspace=str(ws),
                          created_by="test", confirm=True)
    (ws / "done.flag").write_text("")
    svc.start()
    _fast(entry["id"])
    run = _wait_run(entry["id"], {"triggered"})
    assert run["session"] == "" and "done.flag appeared" in run["text"]


def test_a_watcher_that_cannot_work_here_is_paused(svc, tmp_path):
    ws = svc.workspace()
    entry = schedules.add("wenn datei link erscheint", workspace=str(ws),
                          created_by="test", confirm=True)
    os.symlink(tmp_path, ws / "link")
    svc.start()
    _fast(entry["id"])
    run = _wait_run(entry["id"], {"failed"})
    assert "outside" in run["text"]
    assert schedules.get(entry["id"])["paused"] is True
    assert Path(ws / "link").is_symlink()
