"""A wake-up a terminal or TUI session scheduled fires in that session.

Field report (TUI): after a turn ended, its "Wake-up" row stayed on screen
("due") and nothing ever came. A wake-up a session schedules
(schedule_wakeup, owned by that session) was fired only by the dashboard's
listener (tab_agent add_fire_listener); the scheduler thread skips an
owned entry with no listener, and the terminal's idle prompt never looked
-- so in the TUI and the REPL it waited forever and its row never left.

Output, each with its control:

* Scheduler.take_due(owns) hands over the due entries a session owns and
  removes a one-shot (rescheduling an interval); another session's entry
  and one not due yet stay;
* the idle prompt's wake path (TerminalAgent._wake_text, which the service
  worker's agent shares) turns this session's due wake-up into the next
  turn, once; the row it showed is gone from the background list.
"""
from __future__ import annotations

import time
from types import SimpleNamespace

import pytest

from delfin.agent import background_view, job_wake, scheduler
from delfin.agent import repl as R


@pytest.fixture
def sch(tmp_path, monkeypatch):
    s = scheduler.Scheduler(path=tmp_path / "schedule.json")
    monkeypatch.setattr(scheduler, "get_scheduler", lambda path=None: s)
    return s


def _due(entry):
    entry.next_fire_at = time.time() - 1


def test_take_due_hands_over_only_this_sessions_due_entries(sch):
    mine = sch.schedule_once(delay_seconds=60, prompt="check the run", session_id="s1")
    other = sch.schedule_once(delay_seconds=60, prompt="theirs", session_id="s2")
    later = sch.schedule_once(delay_seconds=600, prompt="later", session_id="s1")
    loop = sch.schedule_interval(every_seconds=60, prompt="every minute", session_id="s1")
    for e in (mine, other, loop):
        _due(e)
    got = sch.take_due(lambda owner: owner == "s1")
    assert sorted(e.prompt for e in got) == ["check the run", "every minute"]
    ids = {e.id for e in sch.list_entries()}
    assert mine.id not in ids and other.id in ids and later.id in ids and loop.id in ids
    assert sch.take_due(lambda owner: owner == "s1") == []


@pytest.fixture
def agent(monkeypatch):
    a = R.TerminalAgent.__new__(R.TerminalAgent)
    a.engine = SimpleNamespace(session_id="s1")
    monkeypatch.setattr(job_wake, "wake_enabled", lambda *a, **k: True)
    monkeypatch.setattr(job_wake, "finished_shells", lambda seen, **k: [])
    return a


def test_the_idle_prompt_runs_the_sessions_wake_up_once(agent, sch):
    e = sch.schedule_once(delay_seconds=60, prompt="look at job 42 again", reason="job 42", session_id="s1")
    _due(e)
    text = agent._wake_text("")
    assert "look at job 42 again" in text
    agent._wake_last_look = 0.0
    assert "look at job 42 again" not in agent._wake_text("")
    rows = background_view.collect("/tmp", session_id="s1")["wakeups"]
    assert not [r for r in rows if r["id"] == e.id]


def test_nothing_is_taken_while_something_is_typed(agent, sch):
    e = sch.schedule_once(delay_seconds=60, prompt="look again", session_id="s1")
    _due(e)
    assert agent._wake_text("half a mess") == ""
    assert e.id in {x.id for x in sch.list_entries()}
