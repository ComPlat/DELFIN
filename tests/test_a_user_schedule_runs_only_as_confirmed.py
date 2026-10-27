"""A user schedule runs only as the user confirmed it, on its own host.

* What the user writes is confirmed by writing it; a proposal waits for
  the user's ``ok``; a change to a confirmed schedule -- by any route,
  the file edited by hand included -- withdraws it until confirmed again.
* Due means confirmed, not paused, not already running, on this host, and
  past its time (or asked to run now). Its next time is set when it is
  taken, so a slow start never runs it twice; missed runs are made up
  once within the catch-up window and skipped (with a record) beyond it.
* The day's cost ceiling stops further runs that day, with a record.
* The store is the user's alone: 0700/0600, masked in every cage and
  refused to the agent's file tools.
"""

from __future__ import annotations

import fnmatch
import json
import os
import stat

import pytest

from delfin.agent import schedules


@pytest.fixture
def store(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "user_schedules")
    monkeypatch.setattr(schedules, "_hostname", lambda: "here")
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


T0 = 1_800_000_000.0          # a fixed instant; never this machine's clock


def _add(ws, text="alle 1h Jobs prüfen", **kw):
    kw.setdefault("confirm", True)
    return schedules.add(text, workspace=str(ws), created_by="test", now=T0,
                         tz="UTC", **kw)


def test_what_the_user_writes_runs_at_its_time_once(store):
    entry = _add(store)
    assert entry["confirmed"] and entry["host"] == "here"
    assert schedules.take_due("here", now=T0 + 3599) == []
    due = schedules.take_due("here", now=T0 + 3600)
    assert [d["id"] for d in due] == [entry["id"]]
    # Running: not taken again, whatever the clock says.
    assert schedules.take_due("here", now=T0 + 99999) == []
    schedules.finished(due[0], status="done", text="ok", cost=0.1, now=T0 + 3700)
    again = schedules.get(entry["id"])
    assert again["running"] == "" and again["next_fire_at"] == T0 + 3600 + 3600


def test_a_proposal_waits_for_the_users_ok(store):
    entry = _add(store, confirm=False)
    assert not entry["confirmed"]
    assert schedules.take_due("here", now=T0 + 10 * 3600) == []
    # Confirmed after its first time went by: it looks ahead from then.
    done = schedules.confirm(entry["id"], now=T0 + 10 * 3600,
                             expect=schedules.fingerprint(schedules.get(entry["id"])))
    assert done["next_fire_at"] == T0 + 11 * 3600
    assert [d["id"] for d in schedules.take_due("here", now=T0 + 11 * 3600)] == [entry["id"]]


def test_a_change_after_the_confirmation_withdraws_it(store):
    entry = _add(store)
    shown = schedules.fingerprint(schedules.get(entry["id"]))
    path = schedules._file()
    data = json.loads(path.read_text())
    data["entries"][str(entry["id"])]["task"] = "rm -rf ~"
    path.write_text(json.dumps(data))
    assert not schedules.get(entry["id"])["confirmed"]
    assert schedules.take_due("here", now=T0 + 3600) == []
    # Confirming what was shown before the change is refused.
    with pytest.raises(schedules.ScheduleError):
        schedules.confirm(entry["id"], expect=shown)
    # Asked to run anyway: refused, and the refusal is recorded.
    schedules.request_run(entry["id"])
    assert schedules.take_due("here", now=T0 + 3600) == []
    assert schedules.runs()[-1]["status"] == "unconfirmed"


def test_only_its_own_host_runs_it(store):
    entry = _add(store, host="other-node")
    assert schedules.take_due("here", now=T0 + 10 * 3600) == []
    assert [d["id"] for d in schedules.take_due("other-node", now=T0 + 3600)] == [entry["id"]]
    assert schedules.active_here("here") is False
    assert schedules.active_here("other-node") is True


def test_paused_does_not_run_and_resuming_looks_ahead(store):
    entry = _add(store)
    schedules.pause(entry["id"])
    assert schedules.take_due("here", now=T0 + 5 * 3600) == []
    resumed = schedules.resume(entry["id"], now=T0 + 5 * 3600)
    assert resumed["next_fire_at"] == T0 + 6 * 3600


def test_run_now_runs_before_its_time(store):
    entry = _add(store)
    schedules.request_run(entry["id"])
    assert [d["id"] for d in schedules.take_due("here", now=T0 + 1)] == [entry["id"]]


def test_a_missed_run_is_made_up_once_or_skipped_with_a_record(store):
    entry = _add(store, text="täglich 08:00 Bericht")
    first = schedules.get(entry["id"])["next_fire_at"]
    # A little late: made up once, and the next is tomorrow, not now again.
    due = schedules.take_due("here", now=first + 3600)
    assert [d["id"] for d in due] == [entry["id"]]
    schedules.finished(due[0], status="done", now=first + 3700)
    assert schedules.get(entry["id"])["next_fire_at"] == first + 86400
    # Very late (the host was off for days): skipped, said, rescheduled.
    late = first + 86400 + schedules.CATCH_UP_S + 60
    assert schedules.take_due("here", now=late) == []
    assert schedules.runs()[-1]["status"] == "missed"
    assert schedules.get(entry["id"])["next_fire_at"] > late


def test_the_days_ceiling_stops_more_runs_that_day(store):
    _add(store, budget_run_usd=1.0, budget_day_usd=2.0)
    due = schedules.take_due("here", now=T0 + 3600)
    schedules.finished(due[0], status="done", cost=2.0, now=T0 + 3650)
    assert schedules.take_due("here", now=T0 + 7200) == []
    assert schedules.runs()[-1]["status"] == "budget"
    with pytest.raises(schedules.ScheduleError):
        _add(store, budget_run_usd=5.0, budget_day_usd=2.0)


def test_runs_are_numbered_for_whoever_reports_them(store):
    entry = _add(store)
    due = schedules.take_due("here", now=T0 + 3600)
    one = schedules.finished(due[0], status="failed", text="boom", now=T0 + 3601)
    assert schedules.runs(after_seq=one["seq"]) == []
    assert schedules.last_seq() == one["seq"]
    assert schedules.runs()[-1]["text"] == "boom"
    assert schedules.runs()[-1]["id"] == entry["id"]


def test_a_service_that_starts_releases_runs_of_the_one_before(store):
    entry = _add(store)
    schedules.take_due("here", now=T0 + 3600)
    assert schedules.get(entry["id"])["running"]
    schedules.release_all("here")
    assert schedules.get(entry["id"])["running"] == ""


@pytest.mark.parametrize("text, code", [
    ("alle 1h", "task"), ("Jobs prüfen", "when"), ("alle 1m x", "interval"),
])
def test_what_is_refused_says_which_part(store, text, code):
    with pytest.raises(schedules.ScheduleError) as info:
        _add(store, text=text)
    assert info.value.code == code


def test_the_folder_must_exist(store, tmp_path):
    with pytest.raises(schedules.ScheduleError) as info:
        schedules.add("alle 1h x", workspace=str(tmp_path / "gone"),
                      created_by="t", confirm=True)
    assert info.value.code == "workspace"


def test_the_store_is_the_users_alone(store):
    _add(store)
    d = schedules._DIR
    assert stat.S_IMODE(d.stat().st_mode) == 0o700
    for f in (schedules._file(),):
        assert stat.S_IMODE(f.stat().st_mode) == 0o600


def test_the_store_is_hidden_from_agent_commands_and_file_tools():
    from delfin.agent import api_client, boundaries, sandbox, state_paths
    assert ".delfin/user_schedules" in sandbox._HOME_SECRET_DIRS
    assert ".delfin/user_schedules" in boundaries.MASKED_HOME_ENTRIES
    for path in ("user_schedules/schedules.json", ".delfin/user_schedules/runs.jsonl"):
        assert any(fnmatch.fnmatch(path, g) for g in api_client._DEFAULT_PATH_DENY_GLOBS), path
    assert ("delfin.agent.schedules", "_DIR", "user_schedules") in state_paths.USER_STATE_SINKS


def test_the_command_line_refuses_inside_an_agent_session(monkeypatch, store, capsys):
    from delfin.agent import approval_guard, schedule_commands
    monkeypatch.setenv(approval_guard.ENV_MARKER, "1")
    code = schedule_commands.cli_main(["add", "alle", "1h", "x", "--cwd", str(store)])
    assert code != 0
    assert "agent session" in capsys.readouterr().err
    assert schedules.entries() == []
    assert os.environ.get(approval_guard.ENV_MARKER) == "1"
