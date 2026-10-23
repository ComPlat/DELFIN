"""A wake-up that is long overdue must not fire when the daemon returns.

Observed, not theorised. Starting the daemon after registering the
nightly benchmark fired an unrelated one-shot entry that had been due for
85 minutes — its workspace was a `/tmp` folder from a test run, its
prompt was "check the job", and firing it spent a real agent turn on
nothing before it could be stopped by hand.

`is_due` asks only "is now past the fire time". For a recurring entry
that is right: a missed interval should run at the next opportunity. For
a one-shot it is wrong in a way that costs money. "Wake me in two hours
to check the calculation" is a statement about a moment. Honouring it
three days later, on a machine whose state has moved on, is not a late
version of what was asked for — it is a different request nobody made.

An interval entry keeps firing late on purpose: that IS the schedule.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import scheduler as S


def _scheduler(tmp_path):
    return S.Scheduler(path=tmp_path / "cron.json")


def _fired(sched):
    seen = []
    sched.tick(fire_callback=lambda e: seen.append(e.id))
    return seen


# ---------------------------------------------------------------------------
# One-shot entries
# ---------------------------------------------------------------------------

def test_a_one_shot_that_is_due_now_fires(tmp_path):
    sched = _scheduler(tmp_path)
    ent = sched.schedule_once(delay_seconds=1, prompt="check the job")
    ent.next_fire_at = time.time() - 5
    assert _fired(sched) == [ent.id]


def test_a_one_shot_overdue_by_hours_does_not_fire(tmp_path):
    """The case that spent a turn on a deleted /tmp workspace."""
    sched = _scheduler(tmp_path)
    ent = sched.schedule_once(delay_seconds=7200, prompt="check the job")
    ent.next_fire_at = time.time() - S._STALE_ONCE_GRACE_S - 60
    assert _fired(sched) == []


def test_a_dropped_one_shot_says_why(tmp_path):
    sched = _scheduler(tmp_path)
    ent = sched.schedule_once(delay_seconds=7200, prompt="check the job")
    ent.next_fire_at = time.time() - S._STALE_ONCE_GRACE_S - 60
    sched.tick(fire_callback=lambda e: None)
    stored = [e for e in sched.list_entries() if e.id == ent.id]
    assert stored, "a stale entry was deleted instead of explained"
    assert stored[0].disabled
    assert "overdue" in stored[0].disabled_reason.lower()


def test_a_stale_entry_is_not_re_evaluated_every_poll(tmp_path):
    """Disabled means disabled -- not a fire on the next pass."""
    sched = _scheduler(tmp_path)
    ent = sched.schedule_once(delay_seconds=10, prompt="check the job")
    ent.next_fire_at = time.time() - S._STALE_ONCE_GRACE_S - 60
    sched.tick(fire_callback=lambda e: None)
    assert _fired(sched) == []


def test_the_grace_window_is_generous_enough_for_a_reboot(tmp_path):
    """A laptop closed over lunch must still get its wake-up."""
    assert S._STALE_ONCE_GRACE_S >= 3600
    sched = _scheduler(tmp_path)
    ent = sched.schedule_once(delay_seconds=60, prompt="check the job")
    ent.next_fire_at = time.time() - (S._STALE_ONCE_GRACE_S - 300)
    assert _fired(sched) == [ent.id]


# ---------------------------------------------------------------------------
# Recurring entries are the opposite case
# ---------------------------------------------------------------------------

def test_an_interval_entry_still_fires_after_a_long_gap(tmp_path):
    """Missing a night is not a reason to skip the nightly run."""
    sched = _scheduler(tmp_path)
    ent = sched.schedule_interval(every_seconds=86400, prompt="[bench-nightly]")
    ent.next_fire_at = time.time() - 5 * 86400
    assert _fired(sched) == [ent.id]


def test_an_interval_entry_does_not_replay_every_missed_run(tmp_path):
    """It fires once and is rescheduled forward, not once per missed slot."""
    sched = _scheduler(tmp_path)
    ent = sched.schedule_interval(every_seconds=86400, prompt="[bench-nightly]")
    ent.next_fire_at = time.time() - 5 * 86400
    _fired(sched)
    assert _fired(sched) == []
