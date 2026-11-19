"""A phone's yes or no is read the way people write it, and never guessed.

Field report: "Jo mach" and "ok" as replies to an approval card were
refused; a plain "yes" without the reply gesture was refused even with one
approval open, and never reached the agent when none was.

Output, through the real gateway and bridge (the rig), each with its
control:

* lang.yes_no reads 40+ spellings (German and English, emoji, "please",
  punctuation) and returns None for anything else ("maybe", "ja nein");
* a reply to a card answers it with any of those; anything else is asked
  "Yes or no?" and decides nothing;
* a plain yes/no without the reply gesture answers the approval only when
  exactly one approval card was sent to this chat in the last ten minutes,
  it is still open, and it is at least 2 s old; with two cards, or a card
  that just arrived, it decides nothing and asks for the reply gesture;
* with no approval open, a plain "yes" is a message to the agent.
"""
from __future__ import annotations

import asyncio

import pytest

from _messenger_people_rig import CONFIRM, Account
from delfin.messenger import lang, store
from delfin.messenger.gateway import _PERSON

YES = ["yes", "Yes!", "y", "ja", "Ja.", "j", "jo", "jap", "jep", "yep", "yeah", "yup",
       "ok", "OK", "okay", "okey", "klar", "sicher", "gerne", "gern", "mach", "mach das",
       "Jo, mach das!", "los", "go", "do it", "passt", "genau", "freigeben", "approve",
       "erlauben", "allow", "sure", "of course", "natürlich", "jawohl", "ja bitte",
       "yes please", "👍", "✅", "✔️"]
NO = ["no", "No!", "n", "nein", "Nein.", "ne", "nö", "nee", "nope", "nah", "stop", "stopp",
      "abbrechen", "cancel", "deny", "ablehnen", "lieber nicht", "auf keinen fall",
      "no way", "nein danke", "no thanks", "verweigern", "refuse", "👎", "❌"]
NEITHER = ["maybe", "vielleicht", "ja nein", "yes no", "", "rm -r results", "ja aber erst morgen",
           "why?", "2"]


@pytest.mark.parametrize("text", YES)
def test_every_yes_is_read_as_yes(text):
    assert lang.yes_no(text) == "y"


@pytest.mark.parametrize("text", NO)
def test_every_no_is_read_as_no(text):
    assert lang.yes_no(text) == "n"


@pytest.mark.parametrize("text", NEITHER)
def test_anything_else_is_neither(text):
    assert lang.yes_no(text) is None


def test_the_table_is_large():
    assert len(YES) + len(NO) >= 40


@pytest.fixture
def acct(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    account.service = account.services["anna"]
    return account


def _open(acct, rid="r1"):
    acct.service.pending.setdefault("s-anna", []).append({**CONFIRM, "requestId": rid})


def _cards(acct, age=5.0):
    async def go():
        token = _PERSON.set("anna")
        try:
            await acct.rooms.refresh()
        finally:
            _PERSON.reset(token)
    asyncio.run(go())
    for value in acct.rooms._reply_targets.values():
        value["_sent_at"] -= age
    log = getattr(acct.rooms, "_cards_sent", {})
    for chat, entries in log.items():
        log[chat] = type(entries)((when - age, code) for when, code in entries)
    return [m for (c, m) in acct.rooms._reply_targets if c == "#a1"]


def _replies(acct):
    return acct.service.calls("permission.reply")


def test_a_reply_to_the_card_takes_a_lenient_yes(acct):
    _open(acct)
    mid, = _cards(acct)
    acct.say("21", "#a1", "Jo, mach das!", quote=mid)
    assert [r["choice"] for r in _replies(acct)] == ["y"]


def test_a_reply_to_the_card_takes_a_lenient_no(acct):
    _open(acct)
    mid, = _cards(acct)
    acct.say("21", "#a1", "nö", quote=mid)
    assert [r["choice"] for r in _replies(acct)] == ["n"]


def test_an_unclear_reply_is_asked_yes_or_no_and_decides_nothing(acct):
    _open(acct)
    mid, = _cards(acct)
    acct.say("21", "#a1", "vielleicht", quote=mid)
    assert not _replies(acct)
    assert acct.daemon.to("#a1")[-1] == "Yes or no?"


def test_a_plain_yes_answers_the_one_open_approval(acct):
    _open(acct)
    _cards(acct)
    acct.say("21", "#a1", "ok")
    assert [r["choice"] for r in _replies(acct)] == ["y"]


def test_a_plain_yes_with_two_cards_decides_nothing(acct):
    _open(acct, "r1")
    _cards(acct)
    acct.service.pending["s-anna"] = []
    acct.say("21", "#a1", "pending")          # r1 answered elsewhere
    _open(acct, "r2")
    _cards(acct)
    acct.say("21", "#a1", "yes")
    assert not _replies(acct)
    assert "reply" in acct.daemon.to("#a1")[-1].lower()


def test_a_plain_yes_right_after_the_card_decides_nothing(acct):
    _open(acct)
    _cards(acct, age=0.0)
    acct.say("21", "#a1", "yes")
    assert not _replies(acct)


def test_with_nothing_open_a_plain_yes_is_a_message_to_the_agent(acct):
    acct.say("21", "#a1", "yes")
    assert [p["text"] for p in acct.service.calls("session.prompt")] and not _replies(acct)
    assert acct.service.calls("session.prompt")[-1]["text"].startswith("yes")
