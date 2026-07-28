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


def test_tick_accepts_explicit_callback_argument(fresh_path):
    """The headless daemon passes its callback per tick — no binding."""
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_once(delay_seconds=60, prompt="x")
    ent.next_fire_at = time.time() - 1
    fired: list = []
    n = sch.tick(fire_callback=lambda e: fired.append(e.id))
    assert n == 1
    assert fired == [ent.id]


def test_workspace_and_budget_persist(fresh_path, tmp_path):
    sch = S.Scheduler(path=fresh_path)
    sch.schedule_once(
        delay_seconds=60, prompt="x",
        workspace=str(tmp_path), budget_usd=2.5,
    )
    ent = S.Scheduler(path=fresh_path).list_entries()[0]
    assert ent.workspace == str(tmp_path)
    assert ent.budget_usd == 2.5


def test_workspace_defaults_to_cwd(fresh_path):
    import os
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_interval(every_seconds=120, prompt="x")
    assert ent.workspace == os.getcwd()


def test_disabled_entry_is_never_fired(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_once(delay_seconds=60, prompt="x")
    ent.next_fire_at = time.time() - 1
    ent.disabled = True
    fired: list = []
    assert sch.tick(fire_callback=lambda e: fired.append(e.id)) == 0
    assert fired == []
    assert len(sch.list_entries()) == 1      # kept, just inert


def test_failures_count_and_disable_after_three(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_once(delay_seconds=60, prompt="x")
    ent.next_fire_at = time.time() - 1

    def _boom(_e):
        raise RuntimeError("kaputt")

    for expected in (1, 2, 3):
        assert sch.tick(fire_callback=_boom) == 0
        assert sch.list_entries()[0].fail_count == expected
    ent2 = sch.list_entries()[0]
    assert ent2.disabled is True
    assert "consecutive failures" in ent2.disabled_reason
    assert "kaputt" in ent2.disabled_reason
    # Both the flag and the count survive a reload.
    reloaded = S.Scheduler(path=fresh_path).list_entries()[0]
    assert reloaded.disabled is True
    assert reloaded.fail_count == 3


def test_success_resets_consecutive_failure_count(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_interval(every_seconds=120, prompt="x")
    ent.next_fire_at = time.time() - 1

    def _boom(_e):
        raise RuntimeError("kaputt")

    sch.tick(fire_callback=_boom)
    sch.tick(fire_callback=_boom)
    assert sch.list_entries()[0].fail_count == 2
    assert sch.tick(fire_callback=lambda _e: None) == 1
    assert sch.list_entries()[0].fail_count == 0


def test_disable_entry_exception_disables_immediately(fresh_path):
    sch = S.Scheduler(path=fresh_path)
    ent = sch.schedule_interval(every_seconds=120, prompt="x")
    ent.next_fire_at = time.time() - 1

    def _refuse(_e):
        raise S.DisableEntry("workspace no longer exists: /gone")

    assert sch.tick(fire_callback=_refuse) == 0
    ent2 = sch.list_entries()[0]
    assert ent2.disabled is True
    assert ent2.disabled_reason == "workspace no longer exists: /gone"
    assert ent2.fail_count == 0              # not an ordinary failure
    assert S.Scheduler(path=fresh_path).list_entries()[0].disabled is True


def test_legacy_entry_without_new_fields_still_loads(fresh_path):
    """cron.json written before the daemon fields existed must load."""
    legacy = {
        "entries": [{
            "id": "abc123", "kind": "interval", "prompt": "old",
            "reason": "", "delay_seconds": 0, "every_seconds": 600,
            "created_at": time.time(), "next_fire_at": time.time() + 600,
            "last_fired_at": 0.0, "fire_count": 0,
        }],
    }
    fresh_path.write_text(json.dumps(legacy), encoding="utf-8")
    ent = S.Scheduler(path=fresh_path).list_entries()[0]
    assert ent.id == "abc123"
    assert ent.workspace == ""
    assert ent.disabled is False
    assert ent.fail_count == 0


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
