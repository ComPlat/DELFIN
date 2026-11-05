"""A phone reaches only its own person's sessions, through the real path.

Input: three people on one account -- the owner, anna, ben (people mode) --
each with an agent service of their own behind the REAL bridge, the REAL
gateway in front of it, phones bound in the gateway's store
(tests/_messenger_people_rig.py). Before this, one paired contact could
reach every session of the one service.

Output, each attack with its control:

* anna's phone in anna's chat reaches anna's session (control);
* anna's phone writing in ben's group (forged membership) is refused and
  audited; nothing reaches ben;
* anna's phone naming ben's session (``hol``, ``chat``) reaches nothing;
* an approval code of ben's request, sent from anna's room, is refused by
  the EXECUTOR (the bridge's person context): ben's service sees no
  permission.reply. A bridge command forged with person=anna and ben's
  session id is refused there too; the control, person=ben, applies;
* a revoked phone's late reply to a request it was shown frees nothing;
* anna's forum reaches anna's sessions only; ``@NAME`` addresses one;
* ``neu NAME`` starts the session in anna's own service and folder;
* an unknown or unverified contact in a group is ignored.
"""
from __future__ import annotations

import asyncio

import pytest

from _messenger_people_rig import CONFIRM, Account
from delfin.agent import security_events
from delfin.messenger import store


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    account = Account(tmp_path, monkeypatch)
    store.save("7")
    store.save("21", person="anna", label="anna phone")
    store.save("31", person="ben", label="ben phone")
    account.enable("anna", "s-anna")
    account.enable("ben", "s-ben")
    account.enable("", "s-own")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    store.put_room("#af", None, "Forum", person="anna", ready=True)
    store.put_room("#b1", "s-ben", "Benzol", person="ben", ready=True)
    store.put_room("#o1", "s-own", "Eigene", person="", ready=True)
    yield account
    security_events.clear()


def _prompts(acct, person):
    return [p for p in acct.services[person].calls("session.prompt")
            if not str(p.get("text", "")).startswith("/sync")]


def test_a_phone_reaches_its_persons_session(acct):
    acct.say("21", "#a1", "rechne das")
    assert [p["sessionId"] for p in _prompts(acct, "anna")] == ["s-anna"]
    assert not _prompts(acct, "ben") and not _prompts(acct, "")
    # every bridge command of anna's phone carried anna
    assert all(c.get("person") == "anna" for c in acct.client.commands if c["cmd"] != "pairing")


def test_a_phone_in_another_persons_group_is_refused_and_audited(acct):
    acct.say("21", "#b1", "rechne das")
    assert not _prompts(acct, "ben") and not _prompts(acct, "anna")
    assert not acct.client.commands
    assert "messenger_cross_person_refused" in [e.kind for e in security_events.recent()]
    # the control: ben's phone in ben's group
    acct.say("31", "#b1", "rechne das")
    assert [p["sessionId"] for p in _prompts(acct, "ben")] == ["s-ben"]


@pytest.mark.parametrize("text", ["hol s-ben", "hol Benzol", "chat s-ben"])
def test_a_phone_cannot_name_another_persons_session(acct, text):
    acct.say("21", "@21", text)
    assert not acct.services["ben"].asked
    assert acct.bridge._person_bridge("ben")._record("s-ben").get("by") == "session"
    assert not any(chat == "#b1" for chat, _ in acct.daemon.invited)


def _ben_code(acct):
    acct.services["ben"].pending["s-ben"] = [{**CONFIRM, "requestId": "rb-1"}]
    reply = asyncio.run(acct.bridge._run_command(
        {"cmd": "session", "person": "ben", "session": "s-ben", "action": "pending"}))
    return reply["items"][0]["id"]


def test_an_approval_of_another_persons_request_is_refused_by_the_executor(acct):
    code = _ben_code(acct)
    acct.say("21", "#a1", f"/approve {code}")
    assert not acct.services["ben"].calls("permission.reply")
    assert acct.services["ben"].pending["s-ben"]
    # forged straight at the bridge: anna's context, ben's session and code
    forged = asyncio.run(acct.bridge._run_command(
        {"cmd": "session", "person": "anna", "session": "s-ben", "action": "approve", "id": code}))
    assert forged.get("error") and not acct.services["ben"].calls("permission.reply")
    # anna's context, anna's session, ben's code
    forged = asyncio.run(acct.bridge._run_command(
        {"cmd": "session", "person": "anna", "session": "s-anna", "action": "approve", "id": code}))
    assert forged.get("error") and not acct.services["anna"].calls("permission.reply")
    # the control: ben's phone in ben's room
    acct.say("31", "#b1", f"/approve {code}")
    assert acct.services["ben"].calls("permission.reply") == [
        {"sessionId": "s-ben", "requestId": "rb-1", "choice": "y", "persist": False}]


def test_numbered_approvals_stay_within_the_person(acct):
    # a scheduled session of ben's without a chat: its approvals are numbered
    acct.services["ben"].sessions.append({"sessionId": "s-ben-cron", "title": "nachts",
                                          "workspace": "/w/n", "state": "waiting"})
    acct.services["ben"].pending["s-ben-cron"] = [{**CONFIRM, "requestId": "rb-2"}]
    asyncio.run(_approvals_news_for(acct, "ben"))
    assert any("rm -r results" in t for t in acct.daemon.to("@31"))
    assert not acct.daemon.to("@21") and not acct.daemon.to("@7")
    acct.say("21", "@21", "ok 1")                # anna answers ben's number
    assert not acct.services["ben"].calls("permission.reply")
    acct.say("31", "@31", "ok 1")                # the control
    assert acct.services["ben"].calls("permission.reply")


async def _approvals_news_for(acct, person):
    from delfin.messenger.gateway import _PERSON
    token = _PERSON.set(person)
    try:
        await acct.rooms.approvals_news()
    finally:
        _PERSON.reset(token)


def test_a_revoked_phones_late_reply_frees_nothing(acct):
    acct.services["ben"].pending["s-ben"] = [{**CONFIRM, "requestId": "rb-3"}]
    acct.say("31", "#b1", "pending")
    shown = [k for k in acct.rooms._reply_targets if k[0] == "#b1"]
    assert shown
    store.revoke_contact("31", by="test")
    acct.say("31", "#b1", "ja", quote=shown[0][1])
    assert not acct.services["ben"].calls("permission.reply")


def test_a_persons_forum_reaches_only_their_sessions(acct):
    acct.say("21", "#af", "Zwischenstand bitte")
    assert [p["sessionId"] for p in _prompts(acct, "anna")] == ["s-anna"]
    assert not _prompts(acct, "ben") and not _prompts(acct, "")


def test_at_name_in_the_forum_addresses_one_session(acct):
    acct.enable("anna", "s-anna2")
    store.put_room("#a2", "s-anna2", "Ethan", person="anna", ready=True)
    acct.say("21", "#af", "@Ethan nur du")
    assert [(p["sessionId"], p["text"].split("\n")[0]) for p in _prompts(acct, "anna")] == [
        ("s-anna2", "nur du")]
    acct.say("21", "#af", "@Benzol hallo")       # ben's session: not anna's to address
    assert not _prompts(acct, "ben")
    assert "Keine aktivierte Sitzung" in acct.daemon.to("#af")[-1]


def test_neu_starts_the_session_in_the_persons_own_service_and_folder(acct):
    acct.say("21", "@21", "neu Propan")
    created = acct.services["anna"].calls("session.create")
    assert created and created[0]["cwd"] == str(acct.home / "people" / "anna" / "DELFIN-Handy" / "Propan")
    assert not acct.services[""].calls("session.create") and not acct.services["ben"].calls("session.create")


@pytest.mark.parametrize("sender", ["99", "25"])
def test_unknown_and_unverified_contacts_are_ignored_in_groups(acct, sender):
    if sender == "25":
        store.bind("25", person="anna", label="x", request="r" * 32, code="123456")
    acct.say(sender, "#a1", "rechne das")
    assert not _prompts(acct, "anna") and not acct.daemon.to("#a1")


def test_every_message_in_a_session_chat_names_the_session(acct):
    acct.say("21", "#a1", "status")
    assert acct.daemon.to("#a1")[-1].startswith("[Butan")


def test_a_schedules_command_from_a_person_is_refused(acct):
    acct.say("21", "@21", "zeitpläne")
    assert "Kontoinhaber" in acct.daemon.to("@21")[-1]
    reply = asyncio.run(acct.bridge._run_command({"cmd": "schedules", "person": "anna"}))
    assert reply.get("error")
