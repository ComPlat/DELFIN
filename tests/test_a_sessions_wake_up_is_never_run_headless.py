"""A wake-up a dashboard session scheduled is never run by the headless daemon.

The daemon ticks with its own fire callback, and that override fired every
due entry -- including one a session had scheduled for itself ("check in 30
minutes whether Session A answered"). With the dashboard gone, that became
an agent turn nobody watched. Entries without an owner (the CLI, older
builds) are still the daemon's.
"""
import time

import pytest

from delfin.agent import scheduler as S


@pytest.fixture
def sch(tmp_path):
    return S.Scheduler(path=tmp_path / "cron.json")


def _due(sch, **kw):
    ent = sch.schedule_once(delay_seconds=60, prompt="check", **kw)
    ent.next_fire_at = time.time() - 1
    sch._save()
    return ent


def test_the_daemon_leaves_a_sessions_wake_up_alone(sch):
    ent = _due(sch, session_id="A")
    fired = []
    sch.tick(fire_callback=lambda e: fired.append(e.id))
    assert fired == []
    assert not sch._entries[ent.id].disabled, "still due for its session"


def test_the_daemon_still_runs_an_unowned_wake_up(sch):
    ent = _due(sch)
    fired = []
    sch.tick(fire_callback=lambda e: fired.append(e.id))
    assert fired == [ent.id]
