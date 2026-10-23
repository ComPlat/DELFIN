"""Background work belongs to the session that started it.

Several agent sessions run side by side in one dashboard, and two of them
can work in the same workspace. A wake-up went to whichever session was
built last, and a finished job to whichever session asked first.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import background_view as bgv
from delfin.agent import job_monitor as jm
from delfin.agent import scheduler as S

_SHA_A, _SHA_B = "269429eb", "d507039a"


@pytest.fixture
def sch(tmp_path, monkeypatch):
    scheduler = S.Scheduler(path=tmp_path / "cron.json")
    monkeypatch.setattr(S, "_GLOBAL", scheduler)
    return scheduler


@pytest.fixture
def ws(tmp_path, monkeypatch, sch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    return tmp_path


def _green(url):
    return {"workflow_runs": [
        {"name": "CI", "status": "completed", "conclusion": "success",
         "head_sha": sha + "0" * 32, "html_url": "u", "jobs_url": "j"}
        for sha in (_SHA_A, _SHA_B)]}


# ---------------------------------------------------------------------------
# Wake-ups
# ---------------------------------------------------------------------------

def test_a_wake_up_goes_to_the_session_that_scheduled_it(sch):
    fired_a, fired_b = [], []
    sch.add_fire_listener("a", fired_a.append, lambda owner: owner == "A")
    sch.add_fire_listener("b", fired_b.append, lambda owner: owner == "B")
    ent = sch.schedule_once(delay_seconds=60, prompt="check", session_id="A")
    ent.next_fire_at = time.time() - 1

    assert sch.tick() == 1
    assert [e.id for e in fired_a] == [ent.id] and fired_b == []


def test_a_wake_up_of_a_session_not_open_here_waits_for_it(sch):
    fired = []
    sch.add_fire_listener("a", fired.append, lambda owner: owner == "A")
    ent = sch.schedule_once(delay_seconds=60, prompt="check", session_id="C")
    ent.next_fire_at = time.time() - 1

    assert sch.tick() == 0 and fired == []
    assert ent.id in {e.id for e in sch.list_entries()}


def test_a_wake_up_nobody_owns_wakes_no_session(sch):
    """Scheduled before wake-ups had an owner: a new conversation was woken
    by "wait for the CI result of 8feb1765", which an old one had left."""
    first, last, daemon = [], [], []
    sch.add_fire_listener("a", first.append, lambda owner: owner == "A")
    sch.add_fire_listener("b", last.append, lambda owner: owner == "B")
    ent = sch.schedule_once(delay_seconds=60, prompt="check")
    ent.next_fire_at = time.time() - 1

    assert sch.tick() == 0 and first == [] and last == []
    assert sch.tick(fire_callback=daemon.append) == 1, "the daemon still runs it"


def test_the_owner_is_kept_beside_the_schedule_not_in_it(sch):
    """A build from before this field passes every key of an entry to
    ScheduleEntry: an unknown key made it drop the entry, and then delete
    it as removed elsewhere. The home directory is shared by every login
    node, and not all of them run the same build."""
    ent = sch.schedule_once(delay_seconds=60, prompt="check", session_id="A")
    old_keys = {"id", "kind", "prompt", "reason", "delay_seconds",
                "every_seconds", "created_at", "next_fire_at",
                "last_fired_at", "fire_count", "workspace", "budget_usd",
                "fail_count", "disabled", "disabled_reason"}
    (stored,) = json.loads(sch.path.read_text(encoding="utf-8"))["entries"]
    assert set(stored) <= old_keys
    assert sch.owner_of(ent.id) == "A"


def test_an_entry_from_a_newer_build_is_read_not_dropped(tmp_path):
    path = tmp_path / "cron.json"
    path.write_text(json.dumps({"entries": [{
        "id": "x1", "kind": "once", "prompt": "p",
        "next_fire_at": time.time() + 600, "added_by_a_newer_build": True}]}),
        encoding="utf-8")
    assert [e.id for e in S.Scheduler(path=path).list_entries()] == ["x1"]


# ---------------------------------------------------------------------------
# Watched jobs
# ---------------------------------------------------------------------------

def test_a_finished_job_is_reported_to_its_own_session(ws):
    jid = f"ci:ComPlat/DELFIN@{_SHA_B}"
    jm.register_agent_job(ws, jid, "B's push", extra={"session_id": "B"})

    assert jm.check_agent_jobs(ws, fetch_fn=_green, session_id="A") == []
    done = jm.check_agent_jobs(ws, fetch_fn=_green, session_id="B")
    assert [d["job_id"] for d in done] == [jid]


def test_a_watch_nobody_owns_is_no_sessions(ws):
    jm.register_agent_job(ws, f"ci:ComPlat/DELFIN@{_SHA_A}", "older watch")
    assert jm.check_agent_jobs(ws, fetch_fn=_green, session_id="A") == []
    assert len(jm.check_agent_jobs(ws, fetch_fn=_green)) == 1, "the daemon's"


def test_the_background_panel_shows_its_own_sessions_work(ws, sch):
    jm.register_agent_job(ws, "ci:ComPlat/DELFIN@0badc0de", "left by an old one")
    jm.register_agent_job(ws, f"ci:ComPlat/DELFIN@{_SHA_A}", "A's push",
                          extra={"session_id": "A"})
    jm.register_agent_job(ws, f"ci:ComPlat/DELFIN@{_SHA_B}", "B's push",
                          extra={"session_id": "B"})
    sch.schedule_once(delay_seconds=600, prompt="p", reason="B's wake-up",
                      session_id="B")

    view = bgv.collect(ws, session_id="A")
    assert [w["label"] for w in view["watches"]] == ["A's push"]
    assert view["wakeups"] == []
    everything = bgv.collect(ws)
    assert len(everything["watches"]) == 3 and len(everything["wakeups"]) == 1
