"""The messenger is a SimpleX bot the way SimpleX defines one.

Official sources (simplex-chat repository, bots/README.md and
bots/api/COMMANDS.md, TYPES.md): a bot profile has ``peerType: "bot"``; its
commands are ``preferences.commands`` (ChatBotCommand: command | menu),
shown by the apps as a menu (APIUpdateProfile ``/_profile <userId> <json>``);
per-contact overrides go through APISetContactPrefs
(``/_set prefs @<contactId> <json>``); a group's menu is
``groupPreferences.commands`` of its profile (APIUpdateGroupProfile); a
reply quotes with ``quotedItemId`` in the ComposedMessage (APISendMessages);
delivery receipts are a user setting (``/_set receipts contacts|groups
<userId> on``, chat core parser); members are added with role ``member``
(APIAddMember), which cannot add others.

Output, through the real Daemon command builder (a fake socket) and the
real gateway (the rig), each with its control.
"""
from __future__ import annotations

import asyncio
import json

import pytest

from delfin.messenger import render, simplex, store


class _Socket:
    def __init__(self, daemon, replies):
        self.daemon, self.replies, self.sent = daemon, list(replies), []

    async def send(self, raw):
        payload = json.loads(raw)
        self.sent.append(payload["cmd"])
        reply = self.replies.pop(0) if self.replies else {"type": "cmdOk"}
        self.daemon._pending[payload["corrId"]].set_result(reply)


def _daemon(replies=()):
    daemon = simplex.Daemon(port=1)
    daemon._ws = _Socket(daemon, replies)
    return daemon


USER = {"type": "activeUser", "user": {"userId": 1, "profile": {
    "displayName": "DELFIN-Pi", "fullName": "", "preferences": {"calls": {"allow": "no"}}}}}


def _json_of(cmd, skip):
    return json.loads(cmd.split(" ", skip)[skip])


def test_the_profile_becomes_a_bot_with_a_command_menu():
    daemon = _daemon([USER, USER, {"type": "userProfileUpdated"}])
    asyncio.run(daemon.set_bot_profile(render.menu("en")))
    cmd = daemon._ws.sent[-1]
    assert cmd.startswith("/_profile 1 ")
    profile = _json_of(cmd, 2)
    assert profile["peerType"] == "bot" and profile["displayName"] == "DELFIN-Pi"
    assert profile["preferences"]["calls"] == {"allow": "no"}, "other preferences are kept"
    keywords = {c["keyword"] for c in profile["preferences"]["commands"] if c["type"] == "command"}
    assert {"sessions", "get", "new", "status", "help"} <= keywords


def test_a_contact_gets_the_menu_in_its_language():
    daemon = _daemon()
    asyncio.run(daemon.set_contact_menu("5", render.menu("de")))
    cmd = daemon._ws.sent[-1]
    assert cmd.startswith("/_set prefs @5 ")
    keywords = {c["keyword"] for c in _json_of(cmd, 3)["commands"] if c["type"] == "command"}
    assert {"sitzungen", "hol", "neu", "hilfe"} <= keywords


def test_a_session_chat_keeps_its_menu_when_renamed():
    daemon = _daemon()
    asyncio.run(daemon.rename_group("#9", "Butan", commands=render.group_menu("en")))
    profile = _json_of(daemon._ws.sent[-1], 2)
    assert profile["displayName"] == "Butan"
    assert {c["keyword"] for c in profile["groupPreferences"]["commands"]} >= {"status", "stop"}


def test_a_reply_quotes_the_message_it_answers():
    daemon = _daemon([{"type": "newChatItems", "chatItems": []}])
    asyncio.run(daemon.send_chat("@5", "ok", quote="77"))
    composed = _json_of(daemon._ws.sent[-1], 3)
    assert composed[0]["quotedItemId"] == 77
    control = _daemon([{"type": "newChatItems", "chatItems": []}])
    asyncio.run(control.send_chat("@5", "ok"))
    assert "quotedItemId" not in _json_of(control._ws.sent[-1], 3)[0]


def test_receipts_are_switched_on_best_effort():
    refused = {"type": "chatCmdError", "chatError": {"type": "error", "errorType": {"type": "commandError"}}}
    daemon = _daemon([USER, refused, refused])
    asyncio.run(daemon.enable_receipts())                 # a refusal never raises
    assert any(c.startswith("/_set receipts contacts 1 on") for c in daemon._ws.sent)


def test_members_are_added_as_members():
    daemon = _daemon()
    asyncio.run(daemon.invite("#9", "5"))
    assert daemon._ws.sent[-1] == "/_add #9 5 member"


@pytest.mark.parametrize("lang", ["en", "de"])
def test_every_menu_command_is_one_the_gateway_understands(lang):
    from delfin.messenger.gateway import Rooms
    samples = {"<N>": "1", "<NAME>": "Probe", "<en|de|auto>": "auto"}
    for entry in render.menu(lang) + render.group_menu(lang):
        text = "/" + entry["keyword"] + (" " + samples.get(entry.get("params") or "", "x")
                                         if entry.get("params") else "")
        assert Rooms.understands(text), text


@pytest.fixture
def acct(tmp_path, monkeypatch):
    from _messenger_people_rig import Account
    monkeypatch.delenv("DELFIN_MESSENGER_LANG", raising=False)
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    return account


def test_the_gateway_sets_the_bot_up_and_quotes_its_replies(acct):
    asyncio.run(acct.rooms.bot_setup())
    assert acct.daemon.profile_commands, "the bot profile got its menu"
    acct.say("21", "#a1", "status", mid="4711")
    assert acct.daemon.quotes[-1] == ("#a1", "4711")


def test_a_language_change_changes_the_contact_menu(acct):
    acct.say("21", "@21", "sprache de")
    cid, commands = acct.daemon.contact_menus[-1]
    assert cid == "21" and any(c.get("keyword") == "hol" for c in commands)
