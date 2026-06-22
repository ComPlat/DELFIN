"""Tests for delfin.agent.scheduler."""

from __future__ import annotations

import json
import tempfile
import time
from pathlib import Path

import pytest

from delfin.agent import scheduler as S
from delfin.agent.api_client import _DocToolExecutor


@pytest.fixture
def fresh_path():
    with tempfile.NamedTemporaryFile(
        suffix=".json", delete=False,
    ) as f:
        p = Path(f.name)
    p.unlink(missing_ok=True)
    yield p
    p.unlink(missing_ok=True)


def test_schedule_once_creates_entry(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_once(delay_seconds=300, prompt="check build")
    assert ent.kind == "once"
    assert ent.delay_seconds == 300
    assert ent.next_fire_at > time.time()


def test_interval_rejects_too_short(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    with pytest.raises(ValueError):
        sch.schedule_interval(every_seconds=10, prompt="x")


def test_persistence_round_trip(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    sch.schedule_once(delay_seconds=120, prompt="A")
    sch.schedule_interval(every_seconds=600, prompt="B")
    # New scheduler reads same path
    sch2 = S.Scheduler(path=fresh_path)
    kinds = sorted(e.kind for e in sch2.list_entries())
    assert kinds == ["interval", "once"]


def test_delete_removes_entry(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_once(delay_seconds=60, prompt="x")
    assert sch.delete(ent.id) is True
    assert sch.delete(ent.id) is False
    assert sch.list_entries() == []


def test_tick_fires_due_entries(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    fired: list = []
    sch.set_fire_callback(lambda e: fired.append(e.prompt))
    ent = sch.schedule_once(delay_seconds=60, prompt="ready?")
    # Manually back-date to make it due
    ent.next_fire_at = time.time() - 1
    n = sch.tick()
    assert n == 1
    assert fired == ["ready?"]
    # 'once' entries self-delete after firing
    assert sch.list_entries() == []


def test_interval_reschedules_after_fire(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    sch.set_fire_callback(lambda _e: None)
    ent = sch.schedule_interval(every_seconds=120, prompt="poll")
    ent.next_fire_at = time.time() - 1
    sch.tick()
    assert sch.list_entries() and sch.list_entries()[0].fire_count == 1
    new_next = sch.list_entries()[0].next_fire_at
    assert new_next > time.time() + 100   # ~120s ahead


def test_no_callback_no_fire(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_once(delay_seconds=60, prompt="x")
    ent.next_fire_at = time.time() - 1
    n = sch.tick()
    # No callback bound => no fire
    assert n == 0
    # But entry stays
    assert len(sch.list_entries()) == 1


def test_callback_exception_does_not_crash(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    sch.set_fire_callback(lambda _e: (_ for _ in ()).throw(RuntimeError("boom")))
    ent = sch.schedule_once(delay_seconds=60, prompt="x")
    ent.next_fire_at = time.time() - 1
    # Tick must complete without raising
    sch.tick()


def test_tool_dispatch_schedule_wakeup(fresh_path, monkeypatch):
    monkeypatch.setattr(S, "_DEFAULT_PATH", fresh_path)
    S.reset_scheduler()
    out = _DocToolExecutor().execute(
        "schedule_wakeup",
        {"delay_seconds": 120, "prompt": "later", "reason": "test"},
        permissions=None,
    )
    payload = json.loads(out)
    assert payload["status"] == "ok"
    assert payload["id"]


def test_tool_dispatch_cron_create_and_list(fresh_path, monkeypatch):
    monkeypatch.setattr(S, "_DEFAULT_PATH", fresh_path)
    S.reset_scheduler()
    create_out = _DocToolExecutor().execute(
        "cron_create",
        {"every_seconds": 600, "prompt": "poll", "reason": "test"},
        permissions=None,
    )
    eid = json.loads(create_out)["id"]
    list_out = _DocToolExecutor().execute(
        "cron_list", {}, permissions=None,
    )
    listing = json.loads(list_out)
    assert any(e["id"] == eid for e in listing["entries"])

    delete_out = _DocToolExecutor().execute(
        "cron_delete", {"entry_id": eid}, permissions=None,
    )
    assert json.loads(delete_out)["status"] == "ok"


def test_tool_rejects_too_short_interval(fresh_path, monkeypatch):
    monkeypatch.setattr(S, "_DEFAULT_PATH", fresh_path)
    S.reset_scheduler()
    out = _DocToolExecutor().execute(
        "cron_create",
        {"every_seconds": 5, "prompt": "x"},
        permissions=None,
    )
    payload = json.loads(out)
    assert "error" in payload
