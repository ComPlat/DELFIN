"""Each phone is answered in its person's language: English unless they
write German, or ask for it.

Rules (delfin/messenger/language.py):

* a new contact starts in DELFIN_MESSENGER_LANG, English when unset;
* two German messages in a row switch it to German, two English ones back;
  one stray message never flips it;
* ``language en|de|auto`` / ``sprache en|de|auto`` fixes it (or frees it);
* the gateway's texts (help, notices, cards) follow it; what the agent
  wrote is never translated; one person's language never touches another's;
* the setting survives a gateway restart.

Output, through the real gateway and bridge (the rig), each with its
control.
"""
from __future__ import annotations

import asyncio

import pytest

from _messenger_people_rig import CONFIRM, Account
from delfin.messenger import language, lang, store
from delfin.messenger.gateway import HELP, Rooms, _PERSON

HANDED = "Message delivered."


@pytest.fixture
def acct(tmp_path, monkeypatch):
    monkeypatch.delenv("DELFIN_MESSENGER_LANG", raising=False)
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    store.save("31", person="ben", label="ben phone")
    account.enable("anna", "s-anna")
    account.enable("ben", "s-ben")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    store.put_room("#b1", "s-ben", "Benzol", person="ben", ready=True)
    return account


def last(acct, chat):
    return acct.daemon.to(chat)[-1]


def test_a_new_phone_is_answered_in_english(acct):
    acct.say("21", "@21", "help")
    assert last(acct, "@21") == HELP


def test_the_installation_default_names_the_language_of_a_new_phone(acct, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")
    acct.say("21", "@21", "help")
    assert last(acct, "@21") == lang.DE[HELP]


def test_two_german_messages_switch_and_two_english_switch_back(acct):
    acct.say("21", "#a1", "bitte rechne das noch einmal")
    assert last(acct, "#a1") == HANDED, "one German message never flips it"
    acct.say("21", "#a1", "und dann zeig mir die Energie")
    assert last(acct, "#a1") == lang.DE[HANDED]
    acct.say("21", "#a1", "please show the energy again")
    assert last(acct, "#a1") == lang.DE[HANDED]
    acct.say("21", "#a1", "and then compute the barrier")
    assert last(acct, "#a1") == HANDED


def test_neutral_messages_never_switch(acct):
    for text in ("42", "Butan", "ok", "status"):
        acct.say("21", "#a1", text)
    assert language.of_contact("21") == "en"


def test_an_explicit_choice_holds_until_auto(acct):
    acct.say("21", "@21", "sprache de")
    assert "Deutsch" in last(acct, "@21")
    for text in ("please compute it", "and show it to me"):
        acct.say("21", "#a1", text)
    assert last(acct, "#a1") == lang.DE[HANDED]
    acct.say("21", "@21", "language auto")
    assert language.of_contact("21") == "de"
    acct.say("21", "@21", "language en")
    assert "English" in last(acct, "@21")


def test_the_cards_follow_the_persons_language(acct):
    language.set_contact("21", "de")
    acct.services["anna"].pending["s-anna"] = [{**CONFIRM, "requestId": "r1"}]

    async def go():
        token = _PERSON.set("anna")
        try:
            await acct.rooms.refresh()
        finally:
            _PERSON.reset(token)
    asyncio.run(go())
    assert "Antworte mit *ja* oder *nein*." in last(acct, "#a1")


def test_what_the_agent_wrote_is_never_translated(acct):
    language.set_contact("21", "de")
    row = acct.services["anna"].sessions[0]

    async def go():
        token = _PERSON.set("anna")
        try:
            await acct.rooms.refresh()
        finally:
            _PERSON.reset(token)
    asyncio.run(go())
    row.update(state="idle", answer="The barrier is 12.3 kcal/mol.")
    asyncio.run(go())
    assert last(acct, "#a1") == "The barrier is 12.3 kcal/mol."


def test_one_persons_language_never_touches_anothers(acct):
    language.set_contact("21", "de")
    acct.say("31", "@31", "help")
    assert last(acct, "@31") == HELP


def test_the_setting_survives_a_restart(acct):
    acct.say("21", "@21", "sprache de")
    acct.rooms = Rooms(acct.daemon, acct.client)
    acct.say("21", "@21", "help")
    assert last(acct, "@21") == lang.DE[HELP]
