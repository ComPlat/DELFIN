"""The phone hears what a watcher saw, and lists a fired one as done.

* A watcher's report (no task) arrives as "🔔 #N <event> · ausgelöst"
  with what was seen; a watcher that fired once is listed as done.
* The session and the shell list it the same way; every new phone text
  has an English entry; the help names the watchers.
"""

from __future__ import annotations

import asyncio

from delfin.agent import schedules
from delfin.agent.service.bridge import Bridge
from delfin.messenger import gateway, lang
from delfin.messenger.gateway import Rooms


def test_a_report_says_what_was_seen():
    text = Rooms.run_text({"id": 4, "when": "wenn datei out.txt erscheint",
                           "status": "triggered", "text": "out.txt appeared",
                           "task": "", "costUsd": 0})
    assert text.startswith("🔔 #4 wenn datei out.txt erscheint · ausgelöst")
    assert "out.txt appeared" in text


def test_a_fired_watcher_is_listed_as_done(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "s")
    ws = tmp_path / "w"
    ws.mkdir()
    entry = schedules.add("wenn datei a erscheint x", workspace=str(ws),
                          created_by="t", confirm=True)
    schedules.watch_seen(entry["id"], state="there", fired=True, note="a appeared")
    view = Bridge(allowed_uid=1, service_path=tmp_path / "none.sock")._schedule_view(
        schedules.get(entry["id"]))
    assert view["state"] == "done" and view["next"] == "done"
    from delfin.agent import schedule_commands
    monkeypatch.setattr(schedule_commands, "_service_line", lambda: "")
    assert "done" in schedule_commands.listing(host=entry["host"])


def test_the_new_phone_texts_have_english_entries():
    texts = [gateway.SCHEDULE_HELP, gateway._RUN_WORDS["triggered"][1],
             gateway._STATE_WORDS["done"]]
    assert [x[:50] for x in texts if x not in lang.EN] == []
    assert "wenn" in gateway.SCHEDULE_HELP and "when" in lang.EN[gateway.SCHEDULE_HELP]


def test_a_watcher_is_made_from_the_phone_like_a_schedule(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "s")
    root = tmp_path / "phone"
    root.mkdir()
    monkeypatch.setenv("DELFIN_MESSENGER_WORKSPACES", str(root))
    b = Bridge(allowed_uid=1, service_path=tmp_path / "none.sock")
    made = asyncio.run(b._run_command({"cmd": "schedule", "action": "new",
                                       "text": "wenn datei fertig.txt erscheint"}))
    assert made["ok"] and made["confirmed"] is False
    refused = asyncio.run(b._run_command({"cmd": "schedule", "action": "new",
                                          "text": 'wenn ausgabe "rm -rf x" sich ändert'}))
    assert refused.get("error") and refused["code"] == "command"
    assert refused["error"] in lang.EN
