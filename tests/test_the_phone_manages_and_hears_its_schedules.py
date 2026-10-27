"""The paired phone manages the user's schedules and hears how runs end.

Bridge (``schedules``, ``schedule``):
* a schedule made from the phone is a proposal until the phone's ``ok``,
  which must name the check it was shown -- an ``ok`` without it, or for
  content changed since, shows the schedule again and confirms nothing;
* made in the schedules chat it runs in the host's phone folder
  (``DELFIN_MESSENGER_WORKSPACES``/zeitplaene), made in a session's chat
  in that session's folder -- and only for a session enabled for the
  messenger; the phone never names a folder;
* pause, resume, delete and run-now act by number; run-now only for this
  host's schedules;
* the listing and the runs leave the bridge redacted and cut.

Gateway:
* ``zeitpläne`` opens the chat "<machine> · Zeitpläne"; there ``neu``,
  ``ok N``, ``pause N``, ``weiter N``, ``löschen N``, ``jetzt N``,
  ``pläne`` work; ``plan <when> <task>`` in a session's chat proposes one
  in that session's folder;
* each run is reported there once, as the schedule's ``notify`` says,
  and not again after a restart.
"""

from __future__ import annotations

import asyncio

import pytest

from delfin.agent import schedules
from delfin.agent.service.bridge import COMMANDS, Bridge
from delfin.messenger import activation, store
from delfin.messenger.gateway import Rooms


# ---------------------------------------------------------------------------
# The bridge
# ---------------------------------------------------------------------------

@pytest.fixture
def bridge(tmp_path, monkeypatch):
    root = tmp_path / "phone"
    root.mkdir()
    sess = tmp_path / "sess-folder"
    sess.mkdir()
    monkeypatch.setenv("DELFIN_MESSENGER_WORKSPACES", str(root))
    monkeypatch.setattr(activation, "_DIR", tmp_path / "activation")
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "user_schedules")
    activation.configure("s1")
    b = Bridge(allowed_uid=1, service_path=tmp_path / "none.sock")

    async def request(method, params=None):
        if method == "sessions.list":
            return {"sessions": [{"sessionId": "s1", "title": "t", "workspace": str(sess)},
                                 {"sessionId": "s2", "title": "u", "workspace": str(sess)}]}
        raise AssertionError(method)

    monkeypatch.setattr(b, "_request", request)
    return b, root, sess


def _cmd(b, **params):
    return asyncio.run(b._run_command({"cmd": "schedule", **params}))


def test_the_commands_are_in_the_bridges_fixed_set():
    assert "schedules" in COMMANDS and "schedule" in COMMANDS
    assert set(COMMANDS["schedule"]) == {"action", "id", "text", "session", "check"}


def test_made_on_the_phone_it_waits_for_ok_with_its_check(bridge):
    b, root, _ = bridge
    made = _cmd(b, action="new", text="täglich 07:30 Fasse die Jobs zusammen")
    assert made["ok"] and made["confirmed"] is False
    entry = schedules.get(made["id"])
    assert entry["created_by"] == "phone" and not entry["confirmed"]
    assert entry["workspace"] == str((root / "zeitplaene").resolve())
    # ok without the check: shown again, nothing confirmed.
    again = _cmd(b, action="ok", id=str(made["id"]))
    assert again.get("confirmationRequired") and again["check"] == made["check"]
    assert not schedules.get(made["id"])["confirmed"]
    done = _cmd(b, action="ok", id=str(made["id"]), check=made["check"])
    assert done["ok"] and schedules.get(made["id"])["confirmed"]


def test_an_ok_for_changed_content_is_refused(bridge):
    b, _, _ = bridge
    made = _cmd(b, action="new", text="alle 2h x")
    import json
    data = json.loads(schedules._file().read_text())
    data["entries"][str(made["id"])]["task"] = "something else"
    schedules._file().write_text(json.dumps(data))
    reply = _cmd(b, action="ok", id=str(made["id"]), check=made["check"])
    assert not reply.get("ok") and not schedules.get(made["id"])["confirmed"]


def test_made_in_a_session_chat_it_runs_in_that_sessions_folder(bridge):
    b, _, sess = bridge
    made = _cmd(b, action="new", session="s1", text="mo-fr 18:00 Prüfe")
    assert made["ok"]
    assert schedules.get(made["id"])["workspace"] == str(sess.resolve())
    refused = _cmd(b, action="new", session="s2", text="mo-fr 18:00 Prüfe")
    assert refused.get("error") and len(schedules.entries()) == 1


def test_without_a_phone_folder_the_schedules_chat_cannot_make_one(bridge, monkeypatch):
    b, _, _ = bridge
    monkeypatch.delenv("DELFIN_MESSENGER_WORKSPACES")
    assert _cmd(b, action="new", text="alle 2h x").get("error")
    assert schedules.entries() == []


def test_a_wrong_time_is_answered_in_words(bridge):
    b, _, _ = bridge
    reply = _cmd(b, action="new", text="irgendwann x")
    assert "täglich 07:30" in reply["error"]


def test_pause_resume_run_delete_by_number(bridge):
    b, _, _ = bridge
    n = str(_cmd(b, action="new", text="alle 2h x")["id"])
    assert _cmd(b, action="pause", id=n)["ok"] and schedules.get(n)["paused"]
    assert _cmd(b, action="resume", id=n)["ok"] and not schedules.get(n)["paused"]
    assert _cmd(b, action="run", id=n).get("error")          # not confirmed yet
    check = _cmd(b, action="ok", id=n)["check"]
    _cmd(b, action="ok", id=n, check=check)
    assert _cmd(b, action="run", id=n)["ok"] and schedules.get(n)["run_requested"]
    assert _cmd(b, action="delete", id=n)["ok"] and schedules.entries() == []
    assert _cmd(b, action="pause", id="77").get("error")


def test_the_listing_and_the_runs_leave_redacted(bridge):
    b, _, _ = bridge
    secret = "sk-ant-api03-" + "A" * 40
    n = _cmd(b, action="new", text=f"alle 2h use {secret}")["id"]
    entry = schedules.get(n)
    schedules.record_run(entry, status="done", text=f"result {secret}")
    reply = asyncio.run(b._run_command({"cmd": "schedules"}))
    assert reply["ok"] and reply["schedules"][0]["id"] == n
    assert reply["runs"][0]["status"] == "done"
    assert secret not in str(reply)
    later = asyncio.run(b._run_command({"cmd": "schedules", "after": str(reply["last"])}))
    assert later["runs"] == []


# ---------------------------------------------------------------------------
# The gateway
# ---------------------------------------------------------------------------

class Daemon:
    def __init__(self):
        self.sent, self.groups, self.invited = [], [], []

    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
        return {"chatItems": []}

    async def create_group(self, title):
        self.groups.append(title)
        return "#" + str(500 + len(self.groups))

    async def invite(self, chat, owner):
        self.invited.append((chat, owner))


class FakeBridge:
    def __init__(self):
        self.calls = []
        self.runs = []
        self.shown = {}

    def ask(self, command):
        self.calls.append(command)
        if command["cmd"] == "sessions":
            return {"sessions": [{"id": "s1", "name": "t", "state": "idle",
                                  "messengerEnabled": True}]}
        if command["cmd"] == "session":
            return {"ok": True, "title": "t", "state": "idle", "messengerEnabled": True}
        if command["cmd"] == "schedules":
            after = int(command.get("after") or 0)
            return {"ok": True, "host": "pi", "schedules": [
                {"id": 3, "when": "täglich 07:30", "tz": "Europe/Berlin",
                 "next": "Mi 28.10. 07:30", "task": "Jobs", "state": "active",
                 "host": "pi", "lastStatus": ""}],
                "runs": [r for r in self.runs if r["seq"] > after],
                "last": max([r["seq"] for r in self.runs] or [0])}
        if command["cmd"] == "schedule":
            action = command["action"]
            if action == "new":
                return {"ok": True, "id": 3, "check": "abc", "confirmed": False,
                        "when": "täglich 07:30", "tz": "Europe/Berlin",
                        "next": "Mi 28.10. 07:30", "task": "Jobs", "folder": "zeitplaene",
                        "runUsd": 1.0, "dayUsd": 2.0, "minutes": 30}
            if action == "ok" and command.get("check") != "abc":
                return {"ok": False, "confirmationRequired": True, "id": 3, "check": "abc",
                        "when": "täglich 07:30", "tz": "Europe/Berlin",
                        "next": "Mi 28.10. 07:30", "task": "Jobs", "folder": "zeitplaene",
                        "runUsd": 1.0, "dayUsd": 2.0, "minutes": 30}
            return {"ok": True, "id": 3, "next": "Mi 28.10. 07:30"}
        raise AssertionError(command)


@pytest.fixture
def phone(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "m"))
    monkeypatch.setenv("DELFIN_MESSENGER_LABEL", "Pi")
    store.save("7")
    store.put_room("#10", "s1", "t", owner="7", ready=True)
    daemon, bridge = Daemon(), FakeBridge()
    return Rooms(daemon, bridge), daemon, bridge


def _say(rooms, chat, text, mid):
    asyncio.run(rooms.handle({"chat": chat, "text": text, "sender": "7", "id": mid}))


def test_zeitplaene_opens_the_schedules_chat(phone):
    rooms, daemon, _ = phone
    _say(rooms, "@7", "zeitpläne", "1")
    assert daemon.groups == ["DELFIN-Pi · Zeitpläne"]
    chat = store.schedules_room()
    assert chat and daemon.invited == [(chat, "7")]
    # Asked again: the same chat, no second group.
    _say(rooms, "@7", "zeitpläne", "2")
    assert len(daemon.groups) == 1


def test_the_schedules_chat_makes_confirms_and_manages(phone):
    rooms, daemon, bridge = phone
    _say(rooms, "@7", "zeitpläne", "1")
    chat = store.schedules_room()
    _say(rooms, chat, "neu täglich 07:30 Jobs", "2")
    assert bridge.calls[-1] == {"cmd": "schedule", "action": "new", "text": "täglich 07:30 Jobs"}
    assert "ok 3" in daemon.sent[-1][1]
    _say(rooms, chat, "ok 3", "3")
    assert bridge.calls[-1] == {"cmd": "schedule", "action": "ok", "id": "3", "check": "abc"}
    for word, action in (("pause 3", "pause"), ("weiter 3", "resume"),
                         ("jetzt 3", "run"), ("löschen 3", "delete")):
        _say(rooms, chat, word, word)
        assert bridge.calls[-1]["action"] == action and bridge.calls[-1]["id"] == "3"
    _say(rooms, chat, "pläne", "list")
    assert bridge.calls[-1]["cmd"] == "schedules" and "täglich 07:30" in daemon.sent[-1][1]


def test_an_ok_the_gateway_did_not_show_shows_it_first(phone):
    rooms, daemon, bridge = phone
    _say(rooms, "@7", "zeitpläne", "1")
    chat = store.schedules_room()
    _say(rooms, chat, "ok 3", "2")
    assert bridge.calls[-1].get("check") is None
    assert "ok 3" in daemon.sent[-1][1] and "Jobs" in daemon.sent[-1][1]
    _say(rooms, chat, "ok 3", "3")
    assert bridge.calls[-1].get("check") == "abc"


def test_plan_in_a_session_chat_proposes_in_its_folder(phone):
    rooms, daemon, bridge = phone
    _say(rooms, "#10", "plan mo-fr 18:00 Prüfe", "1")
    assert {"cmd": "schedule", "action": "new", "session": "s1",
            "text": "mo-fr 18:00 Prüfe"} in bridge.calls
    assert "Zeitpläne" in daemon.sent[-1][1]


def test_each_run_is_reported_once_as_its_schedule_says(phone):
    rooms, daemon, bridge = phone
    _say(rooms, "@7", "zeitpläne", "1")
    chat = store.schedules_room()
    bridge.runs = [{"seq": 1, "id": 3, "status": "done", "text": "alles gut",
                    "costUsd": 0.04, "task": "Jobs", "when": "täglich 07:30",
                    "notify": "always"}]
    asyncio.run(rooms.schedule_news())
    # The first look only learns where the record stands: no backlog.
    assert not any("alles gut" in text for _, text in daemon.sent)
    bridge.runs.append({"seq": 2, "id": 3, "status": "done", "text": "zweiter Lauf",
                        "costUsd": 0.0, "task": "Jobs", "when": "täglich 07:30",
                        "notify": "always"})
    bridge.runs.append({"seq": 3, "id": 3, "status": "done", "text": "still",
                        "costUsd": 0.0, "task": "Jobs", "when": "täglich 07:30",
                        "notify": "errors"})
    bridge.runs.append({"seq": 4, "id": 3, "status": "failed", "text": "kaputt",
                        "costUsd": 0.0, "task": "Jobs", "when": "täglich 07:30",
                        "notify": "errors"})
    asyncio.run(rooms.schedule_news())
    texts = [text for c, text in daemon.sent if c == chat]
    assert any("zweiter Lauf" in t for t in texts)
    assert not any("still" in t for t in texts)
    assert any("kaputt" in t for t in texts)
    count = len(daemon.sent)
    # Another gateway (a restart) reports nothing twice.
    Rooms(daemon, bridge)
    asyncio.run(Rooms(daemon, bridge).schedule_news())
    assert len(daemon.sent) == count


def test_every_schedule_text_the_phone_can_get_has_an_english_entry():
    from delfin.agent.service import bridge as bridge_mod
    from delfin.messenger import gateway, lang
    texts = [gateway.SCHEDULE_HELP, *bridge_mod.SCHEDULE_ERRORS.values(),
             *(word for _, word in gateway._RUN_WORDS.values()),
             *gateway._STATE_WORDS.values()]
    assert [x[:50] for x in texts if x not in lang.EN] == []
