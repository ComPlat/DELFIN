"""Watchers: a schedule that runs when something happens, not at a time.

Written like a schedule, the event first, then the task (none: only a
report)::

    wenn job 123456 fertig Fasse die Ergebnisse zusammen
    when pid 4711 exits
    wenn datei out/result.xyz erscheint Werte sie aus
    when file calc.log changes every 2 min
    wenn ausgabe "squeue -u max" sich ändert

* The check is cheap and runs no model: the agent service looks at the
  job queue, a process, a file's existence or stamp, or a command's output.
  Only when the event happens does the task run (as a schedule's run, the
  event named in its prompt); without a task the event is only reported.
* "fertig", "endet", "erscheint" fire once; "ändert sich" fires on every
  change, at most once per 5 minutes when a task runs, and the first look
  only learns the state.
* A command is one of a fixed set of reading commands (squeue, sacct, ls,
  tail, ...), run without a shell, with a clean environment and a time
  limit; a path -- a watched file or a command's file argument -- must lie
  inside the schedule's folder. So a watcher, from whichever surface, can
  neither change anything nor read past its folder.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time

import pytest

from delfin.agent import schedules, watch_checks


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "user_schedules")
    monkeypatch.setattr(schedules, "_hostname", lambda: "here")
    folder = tmp_path / "ws"
    folder.mkdir()
    return folder


@pytest.mark.parametrize("text, what, target, event, rest", [
    ("wenn job 123456 fertig Fasse zusammen", "job", "123456", "done", "Fasse zusammen"),
    ("when job 77_3 finishes summarise", "job", "77_3", "done", "summarise"),
    ("nach job 9 Bericht", "job", "9", "done", "Bericht"),
    ("when pid 4711 exits", "pid", "4711", "done", ""),
    ("wenn prozess 4711 endet: los", "pid", "4711", "done", "los"),
    ("wenn datei out/result.xyz erscheint Werte aus", "file", "out/result.xyz", "appears", "Werte aus"),
    ("when file out.txt appears", "file", "out.txt", "appears", ""),
    ("sobald datei a.log existiert x", "file", "a.log", "appears", "x"),
    ("wenn datei calc.log sich ändert x", "file", "calc.log", "changes", "x"),
    ("when file calc.log changes every 2 min x", "file", "calc.log", "changes", "x"),
    ('wenn ausgabe "squeue -u max" sich ändert x', "output", "squeue -u max", "changes", "x"),
    ("when output `ls runs` changes", "output", "ls runs", "changes", ""),
])
def test_the_event_is_read_as_written(text, what, target, event, rest):
    spec, after = schedules.split_when(text)
    assert (spec["kind"], spec["what"], spec["target"], spec["event"], after) == \
        ("watch", what, target, event, rest)


def test_the_check_interval_is_bounded():
    spec, _ = schedules.split_when("when file a changes every 2 min")
    assert spec["every"] == 120
    spec, _ = schedules.split_when("wenn datei a erscheint")
    assert spec["every"] == schedules.WATCH_EVERY_S
    with pytest.raises(schedules.ScheduleError):
        schedules.split_when("when file a changes every 5s")


@pytest.mark.parametrize("text, code", [
    ('wenn ausgabe "rm -rf runs" sich ändert', "command"),
    ('wenn ausgabe "python -c 1" sich ändert', "command"),
    ('wenn ausgabe "grep -f /etc/shadow x" sich ändert', "command"),
    ('wenn ausgabe "tail /etc/passwd" sich ändert', "path"),
    ('wenn ausgabe "sacct -o /tmp/x" sich ändert', "path"),
    ("wenn datei /etc/passwd erscheint", "path"),
    ("wenn datei ../outside.txt erscheint", "path"),
    ("wenn job abc fertig", "when"),
    ("wenn datei erscheint", "when"),
])
def test_what_a_watcher_may_not_do_is_refused(ws, text, code):
    with pytest.raises(schedules.ScheduleError) as info:
        schedules.add(text, workspace=str(ws), created_by="t", confirm=True)
    assert info.value.code == code


def test_a_path_that_leaves_the_folder_by_a_link_is_refused_at_the_check(ws, tmp_path):
    (tmp_path / "secret").write_text("x")
    entry = schedules.add("wenn datei link erscheint", workspace=str(ws),
                          created_by="t", confirm=True)
    os.symlink(tmp_path / "secret", ws / "link")
    with pytest.raises(schedules.ScheduleError):
        watch_checks.check(schedules.parse_when(entry["when"]), str(ws), None)


# -- the checks ----------------------------------------------------------------

def _spec(text):
    return schedules.parse_when(text)


def test_a_file_that_appears_fires_once(ws):
    spec = _spec("wenn datei out/r.txt erscheint")
    state, fired, _ = watch_checks.check(spec, str(ws), None)
    assert not fired
    (ws / "out").mkdir()
    (ws / "out" / "r.txt").write_text("1")
    state, fired, note = watch_checks.check(spec, str(ws), state)
    assert fired and "out/r.txt" in note


def test_a_change_fires_after_the_first_look(ws):
    spec = _spec("wenn datei calc.log sich ändert")
    (ws / "calc.log").write_text("a")
    state, fired, _ = watch_checks.check(spec, str(ws), None)
    assert not fired
    state, fired, _ = watch_checks.check(spec, str(ws), state)
    assert not fired
    time.sleep(0.02)
    (ws / "calc.log").write_text("ab")
    state, fired, _ = watch_checks.check(spec, str(ws), state)
    assert fired


def test_a_process_that_ends_fires(ws):
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        spec = _spec(f"when pid {proc.pid} exits")
        state, fired, _ = watch_checks.check(spec, str(ws), None)
        assert not fired
    finally:
        proc.kill()
        proc.wait()
    _, fired, note = watch_checks.check(spec, str(ws), state)
    assert fired and str(proc.pid) in note


def _fake_slurm(tmp_path, monkeypatch, queued: bool):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    (bin_dir / "squeue").write_text(
        "#!/bin/sh\n" + ("echo RUNNING\n" if queued else "") + "exit 0\n")
    (bin_dir / "sacct").write_text("#!/bin/sh\necho COMPLETED\n")
    for name in ("squeue", "sacct"):
        (bin_dir / name).chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")


def test_a_job_that_leaves_the_queue_fires_with_its_end_state(ws, tmp_path, monkeypatch):
    spec = _spec("wenn job 123 fertig")
    _fake_slurm(tmp_path, monkeypatch, queued=True)
    state, fired, _ = watch_checks.check(spec, str(ws), None)
    assert not fired and state == "RUNNING"
    _fake_slurm(tmp_path, monkeypatch, queued=False)
    _, fired, note = watch_checks.check(spec, str(ws), state)
    assert fired and "COMPLETED" in note


def test_a_command_output_change_fires_and_the_command_gets_no_secrets(ws, monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "sk-should-not-pass")
    spec = _spec('wenn ausgabe "ls" sich ändert')
    state, fired, _ = watch_checks.check(spec, str(ws), None)
    assert not fired
    (ws / "new.txt").write_text("x")
    state, fired, note = watch_checks.check(spec, str(ws), state)
    assert fired and "new.txt" in note
    assert "OPENAI_API_KEY" not in watch_checks._command_env()


# -- in the store ----------------------------------------------------------------

def test_a_watcher_never_fires_by_the_clock(ws):
    entry = schedules.add("wenn datei a.txt erscheint Werte aus", workspace=str(ws),
                          created_by="t", confirm=True, now=1000.0)
    assert schedules.take_due("here", now=10 ** 10) == []
    assert [e["id"] for e in schedules.watches_to_check("here", now=1000.0 + 61)] == [entry["id"]]
    # Looked at: not again before its interval.
    assert schedules.watches_to_check("here", now=1000.0 + 62) == []


def test_an_event_with_a_task_starts_a_run_naming_it(ws):
    entry = schedules.add("wenn datei a.txt erscheint Werte aus", workspace=str(ws),
                          created_by="t", confirm=True, now=1000.0)
    schedules.watch_seen(entry["id"], state="x", fired=True, note="a.txt appeared", now=2000.0)
    due = schedules.take_due("here", now=2000.0)
    assert [d["id"] for d in due] == [entry["id"]]
    assert due[0]["trigger_note"] == "a.txt appeared"
    # "erscheint" fires once: it is not looked at again.
    assert schedules.get(entry["id"])["watch_done"] is True
    assert schedules.watches_to_check("here", now=10 ** 10) == []


def test_an_event_without_a_task_is_only_reported(ws):
    entry = schedules.add("wenn datei a.txt sich ändert", workspace=str(ws),
                          created_by="t", confirm=True, now=1000.0)
    assert entry["task"] == ""
    schedules.watch_seen(entry["id"], state="y", fired=True, note="a.txt changed", now=2000.0)
    assert schedules.take_due("here", now=2000.0) == []
    run = schedules.runs()[-1]
    assert run["status"] == "triggered" and run["text"] == "a.txt changed"
    # A change repeats: still watched.
    assert not schedules.get(entry["id"])["watch_done"]


def test_a_change_with_a_task_waits_five_minutes_before_the_next_look(ws):
    entry = schedules.add("wenn datei a.txt sich ändert Werte aus", workspace=str(ws),
                          created_by="t", confirm=True, now=1000.0)
    schedules.watch_seen(entry["id"], state="y", fired=True, note="changed", now=2000.0)
    assert schedules.get(entry["id"])["next_fire_at"] >= 2000.0 + schedules.MIN_INTERVAL_S


def test_an_unconfirmed_watcher_is_not_looked_at(ws):
    schedules.add("wenn datei a.txt erscheint x", workspace=str(ws),
                  created_by="agent", confirm=False, now=1000.0)
    assert schedules.watches_to_check("here", now=10 ** 10) == []


@pytest.mark.parametrize("command", [
    "du --files0-from=/etc/x", "wc --files0-from=list", "ls --directory=/etc",
    "grep --include=/etc/x pat", "find . -fprint out",
])
def test_an_option_cannot_carry_a_path_out_or_read_a_list(ws, command):
    with pytest.raises(schedules.ScheduleError) as info:
        schedules.add(f'wenn ausgabe "{command}" sich ändert', workspace=str(ws),
                      created_by="t", confirm=True)
    assert info.value.code in ("command", "path")
