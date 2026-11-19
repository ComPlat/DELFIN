"""A disabled schedule has a way back, and a damaged cron.json is kept.

Two defects of delfin/agent/scheduler.py, both measured on the shipped
code before the change these tests pin:

1. ``Scheduler._disable`` is reached from four places (three consecutive
   failures, a one-shot past its grace, ``DisableEntry`` from the fire
   callback, the emergency stop) and nothing in ``delfin/`` set
   ``disabled`` back to False. The only documented repair was editing
   ``~/.delfin/cron.json`` by hand.

2. ``_read_file_or_none`` distinguishes "absent" from "unreadable", but
   ``_save`` read through the lenient ``_read_file``, which answers ``{}``
   for a corrupt file. One ``add()`` from a fresh process therefore
   replaced the file with the single entry just created.

The emergency-stop case is the one refusal that matters: a stop is a
person saying "nothing of mine runs", and a repair command must not undo
it. The test builds that reason through ``stop_all.request`` itself
rather than writing the string, so the two cannot drift apart.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent import scheduler as S


@pytest.fixture
def cron(tmp_path) -> Path:
    return tmp_path / ".delfin" / "cron.json"


def _disabled_entry(sch: S.Scheduler, reason: str, *,
                    workspace: str = "") -> S.ScheduleEntry:
    ent = sch.schedule_interval(every_seconds=3600, prompt="watch the run",
                                reason="nightly watch",
                                workspace=workspace or str(Path.cwd()))
    ent.fail_count = 3
    sch._disable(ent, reason)
    sch._save()
    return ent


# --- the way back ------------------------------------------------------

def test_enable_clears_the_disable_the_reason_and_the_failure_count(
        cron, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    sch = S.Scheduler(path=cron)
    ent = _disabled_entry(sch, "disabled after 3 consecutive failures "
                               "(last: boom)", workspace=str(ws))

    changed, message = sch.enable(ent.id)

    assert changed is True
    assert ent.id in message
    again = S.Scheduler(path=cron).list_entries()
    assert len(again) == 1
    assert again[0].disabled is False
    assert again[0].disabled_reason == ""
    assert again[0].fail_count == 0


def test_an_enabled_entry_fires_again(cron, tmp_path):
    """The point of the repair: the entry runs, it is not merely a flag."""
    ws = tmp_path / "ws"
    ws.mkdir()
    sch = S.Scheduler(path=cron)
    ent = _disabled_entry(sch, "disabled after 3 consecutive failures "
                               "(last: boom)", workspace=str(ws))
    ent.next_fire_at = time.time() - 1
    assert sch.tick(lambda e: None) == 0        # disabled: never fired

    sch.enable(ent.id)
    assert sch.tick(lambda e: None) == 1


def test_enabling_an_enabled_entry_is_a_no_op_not_an_error(cron, tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_interval(every_seconds=3600, prompt="x",
                                workspace=str(ws))

    changed, message = sch.enable(ent.id)

    assert changed is False
    assert message                              # it says what it found
    assert sch.list_entries()[0].disabled is False


def test_enable_names_an_id_it_does_not_have(cron):
    changed, message = S.Scheduler(path=cron).enable("deadbeef00")
    assert changed is False
    assert "deadbeef00" in message


def test_enable_refuses_while_the_workspace_is_gone(cron, tmp_path):
    """Re-enabling an entry whose workspace is gone only re-disables it on
    the next tick, so the refusal is the honest answer."""
    ws = tmp_path / "gone"
    ws.mkdir()
    sch = S.Scheduler(path=cron)
    ent = _disabled_entry(sch, f"workspace {ws} no longer exists",
                          workspace=str(ws))
    ws.rmdir()

    changed, message = sch.enable(ent.id)

    assert changed is False
    assert str(ws) in message
    assert S.Scheduler(path=cron).list_entries()[0].disabled is True


def test_enable_refuses_an_entry_the_emergency_stop_disabled(
        cron, tmp_path, monkeypatch):
    """The reason comes from stop_all itself, not from a copy of its text."""
    from delfin.agent import stop_all

    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.setattr(S, "_DEFAULT_PATH", cron)
    monkeypatch.setattr(stop_all, "_PATH", tmp_path / ".delfin" / "stop.json")
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_interval(every_seconds=3600, prompt="x",
                                workspace=str(ws))

    summary = stop_all.request("test", end_here=False)
    assert summary["schedules_disabled"] == 1

    reread = S.Scheduler(path=cron)
    changed, message = reread.enable(ent.id)

    assert changed is False
    assert "emergency stop" in message.lower()
    assert S.Scheduler(path=cron).list_entries()[0].disabled is True


def test_the_emergency_stop_marker_matches_what_stop_all_writes(
        cron, tmp_path, monkeypatch):
    """A unit check on the predicate, so a reworded reason fails here and
    not only in the integration test above."""
    from delfin.agent import stop_all

    monkeypatch.setattr(S, "_DEFAULT_PATH", cron)
    monkeypatch.setattr(stop_all, "_PATH", tmp_path / ".delfin" / "stop.json")
    S.Scheduler(path=cron).schedule_interval(every_seconds=3600, prompt="x",
                                             workspace=str(tmp_path))
    stop_all.request("test", end_here=False)

    reason = S.Scheduler(path=cron).list_entries()[0].disabled_reason
    assert S.disabled_by_emergency_stop(reason) is True
    assert S.disabled_by_emergency_stop(
        "disabled after 3 consecutive failures (last: boom)") is False
    assert S.disabled_by_emergency_stop("") is False


# --- the store a corruption must survive -------------------------------

_CORRUPT = '{"entries": [{"id": "aaa1", "kind": "inter'


def test_a_fresh_process_refuses_to_write_over_a_corrupt_store(cron):
    cron.parent.mkdir(parents=True, exist_ok=True)
    cron.write_text(_CORRUPT, encoding="utf-8")

    sch = S.Scheduler(path=cron)                # loads nothing
    assert sch.list_entries() == []
    note = sch._save()

    assert cron.read_text(encoding="utf-8") == _CORRUPT
    assert "not changed" in note.lower()
    assert not list(cron.parent.glob("cron.json.unreadable-*"))


def test_a_corrupt_store_is_moved_aside_before_it_is_written(cron, tmp_path):
    good = S.Scheduler(path=cron)
    kept = good.schedule_interval(every_seconds=3600, prompt="keep me",
                                  workspace=str(tmp_path))
    cron.write_text(_CORRUPT, encoding="utf-8")

    good.schedule_once(delay_seconds=60, prompt="new",
                       workspace=str(tmp_path))
    note = good.take_store_note()

    aside = list(cron.parent.glob("cron.json.unreadable-*"))
    assert len(aside) == 1
    assert aside[0].read_text(encoding="utf-8") == _CORRUPT
    assert aside[0].name in note
    ids = {e.id for e in S.Scheduler(path=cron).list_entries()}
    assert kept.id in ids                       # memory survived the move
    assert json.loads(cron.read_text(encoding="utf-8"))["entries"]


def test_a_second_corruption_does_not_overwrite_the_first_rescue(
        cron, tmp_path):
    """os.replace overwrites its destination without a word, so the dated
    name is reserved with O_EXCL. Two corruptions inside one second
    otherwise leave only the second file."""
    sch = S.Scheduler(path=cron)
    sch.schedule_interval(every_seconds=3600, prompt="a",
                          workspace=str(tmp_path))
    first = '{"entries": [{"id": "first corruption'
    second = '{"entries": [{"id": "second corruption'
    for text in (first, second):
        cron.write_text(text, encoding="utf-8")
        sch._save()

    aside = sorted(p.read_text(encoding="utf-8")
                   for p in cron.parent.glob("cron.json.unreadable-*"))
    assert aside == sorted([first, second])


def test_an_intact_store_is_never_quarantined(cron, tmp_path):
    """Control: none of the above touches an ordinary schedule."""
    sch = S.Scheduler(path=cron)
    sch.schedule_interval(every_seconds=3600, prompt="a",
                          workspace=str(tmp_path))
    sch.schedule_once(delay_seconds=60, prompt="b", workspace=str(tmp_path))

    assert sch.take_store_note() == ""
    assert not list(cron.parent.glob("cron.json.unreadable-*"))
    assert len(S.Scheduler(path=cron).list_entries()) == 2


def test_store_problem_tells_a_reader_the_file_is_damaged(cron):
    cron.parent.mkdir(parents=True, exist_ok=True)
    assert S.Scheduler(path=cron).store_problem() == ""
    cron.write_text(_CORRUPT, encoding="utf-8")
    assert str(cron) in S.Scheduler(path=cron).store_problem()


# --- the route that is not a hand-edited file ---------------------------

def test_the_cli_enables_an_entry_and_says_so(cron, tmp_path, capsys,
                                              monkeypatch):
    from delfin.agent import cli

    ws = tmp_path / "ws"
    ws.mkdir()
    monkeypatch.setattr(S, "_DEFAULT_PATH", cron)
    sch = S.Scheduler(path=cron)
    ent = _disabled_entry(sch, "disabled after 3 consecutive failures "
                               "(last: boom)", workspace=str(ws))

    rc = cli.main(["scheduler", "enable", ent.id])

    out = capsys.readouterr().out
    assert rc == 0, out
    assert ent.id in out
    assert S.Scheduler(path=cron).list_entries()[0].disabled is False


def test_the_cli_refuses_an_emergency_stop_disable_with_a_non_zero_exit(
        cron, tmp_path, capsys, monkeypatch):
    from delfin.agent import cli, stop_all

    monkeypatch.setattr(S, "_DEFAULT_PATH", cron)
    monkeypatch.setattr(stop_all, "_PATH", tmp_path / ".delfin" / "stop.json")
    ent = S.Scheduler(path=cron).schedule_interval(
        every_seconds=3600, prompt="x", workspace=str(tmp_path))
    stop_all.request("test", end_here=False)

    rc = cli.main(["scheduler", "enable", ent.id])

    captured = capsys.readouterr()
    assert rc != 0
    assert "emergency stop" in (captured.out + captured.err).lower()


def test_the_cli_status_marks_a_disabled_entry_and_names_the_repair(
        cron, tmp_path, capsys, monkeypatch):
    from delfin.agent import cli

    monkeypatch.setattr(S, "_DEFAULT_PATH", cron)
    sch = S.Scheduler(path=cron)
    ent = _disabled_entry(sch, "disabled after 3 consecutive failures "
                               "(last: boom)", workspace=str(tmp_path))

    assert cli.main(["scheduler", "status"]) == 0

    out = capsys.readouterr().out
    assert "DISABLED" in out
    assert "consecutive failures" in out
    assert f"scheduler enable {ent.id}" in out


def test_enable_over_a_damaged_file_says_so_instead_of_not_found(cron):
    """Every id looks absent in a file that cannot be read, and "no such
    entry" sends the person looking for an id that is in it."""
    cron.parent.mkdir(parents=True, exist_ok=True)
    cron.write_text(_CORRUPT, encoding="utf-8")

    changed, message = S.Scheduler(path=cron).enable("aaa1")

    assert changed is False
    assert str(cron) in message
    assert "cannot be read" in message


def test_the_cli_warns_when_a_write_moved_the_damaged_file_aside(
        cron, tmp_path, capsys, monkeypatch):
    """`scheduler add-bench` is the one CLI path that writes an entry; the
    note from that write must reach the person, not only the file."""
    from delfin.agent import cli

    monkeypatch.setattr(S, "_DEFAULT_PATH", cron)
    cron.parent.mkdir(parents=True, exist_ok=True)
    cron.write_text(_CORRUPT, encoding="utf-8")

    rc = cli.main(["scheduler", "add-bench", "--model", "m", "--every", "24h"])

    captured = capsys.readouterr()
    assert rc == 0
    aside = list(cron.parent.glob("cron.json.unreadable-*"))
    assert len(aside) == 1
    assert aside[0].read_text(encoding="utf-8") == _CORRUPT
    assert aside[0].name in captured.err
    assert len(S.Scheduler(path=cron).list_entries()) == 1


def test_the_agent_facing_list_shows_the_disable_and_its_reason():
    """cron_list answered id/kind/next_fire_at and never `disabled`, so an
    entry that had stopped read as one that was simply not due yet."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client._DocToolExecutor)
    block = src.split('if name == "cron_list"', 1)
    assert len(block) == 2, "cron_list moved; find it and re-pin this"
    body = block[1].split('if name == "cron_delete"', 1)[0]
    assert '"disabled"' in body
    assert '"disabled_reason"' in body
