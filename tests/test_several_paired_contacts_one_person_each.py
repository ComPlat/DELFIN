"""Several paired contacts, each with the person it belongs to.

Input: one gateway process, its pairing store (delfin/messenger/store.py
format 2) with TWO verified contacts of TWO people (owner ``""`` and
``anna``), each contact's own room in its person's context.

Output:

* every verified contact of every person is served; each message is routed
  by its SENDER's person (store.verified_person), never by the process's
  stored owner contact;
* an unknown contact is still refused exactly as today: no reply, nothing
  dispatched, result refused:unknown in the eventlog;
* a person's sessions stay reachable only through that person's contacts
  (cross-person routing refused, as today).
"""
from __future__ import annotations

import pytest

from _messenger_people_rig import Account
from delfin.agent import security_events
from delfin.messenger import eventlog, store


@pytest.fixture(autouse=True)
def _german_phone(monkeypatch):
    """These tests read the German phone texts (lang.DE); English is the
    default (test_the_phone_speaks_the_persons_language)."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    account = Account(tmp_path, monkeypatch)
    # Two contacts of two people in ONE process's store: the head's and
    # anna's. The head is a person like the others here (part 4).
    store.save("7", person="", label="head phone")
    store.save("21", person="anna", label="anna phone")
    account.enable("", "s-own")
    account.enable("anna", "s-anna")
    store.put_room("#o1", "s-own", "Eigene", person="", ready=True)
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    yield account
    security_events.clear()


def _prompts(acct, person):
    return [p for p in acct.services[person].calls("session.prompt")
            if not str(p.get("text", "")).startswith("/sync")]


def test_two_contacts_one_process_each_reaches_its_own_persons_session(acct):
    acct.say("7", "#o1", "rechne das")
    acct.say("21", "#a1", "rechne das")
    assert [p["sessionId"] for p in _prompts(acct, "")] == ["s-own"]
    assert [p["sessionId"] for p in _prompts(acct, "anna")] == ["s-anna"]
    # each bridge command carried its own person ("" is never set: absent
    # means the owner, gateway.py `if person:` skips a falsy person)
    assert all(c.get("person", "") == "" for c in acct.client.commands
               if c["cmd"] == "session" and c.get("session") == "s-own")
    assert all(c.get("person") == "anna" for c in acct.client.commands
               if c["cmd"] == "session" and c.get("session") == "s-anna")


def test_contact_a_never_continues_bs_session(acct):
    """The whole point: person A's contact cannot drive person B's session,
    even by naming it or writing in B's room (part 2, pinned as attack)."""
    acct.say("21", "#o1", "rechne das")          # anna writes in the head's room
    assert not _prompts(acct, "") and not _prompts(acct, "anna")
    assert "messenger_cross_person_refused" in [e.kind for e in security_events.recent()]
    # the control: the head's contact in its own room still works
    acct.say("7", "#o1", "rechne das")
    assert [p["sessionId"] for p in _prompts(acct, "")] == ["s-own"]


def test_unknown_contact_refused_exactly_as_today(acct, tmp_path, monkeypatch):
    """Same refusal as the single-contact gateway: no reply, no dispatch,
    result refused:unknown. The refusal text says no more than today."""
    events = []
    monkeypatch.setattr(eventlog, "emit",
                        lambda direction, kind, **kw: events.append((direction, kw)))
    before = len(acct.daemon.sent)
    acct.say("99", "#o1", "rechne das")
    assert len(acct.daemon.sent) == before       # nothing answered
    assert not _prompts(acct, "") and not _prompts(acct, "anna")
    assert ("in", {"room": "#o1", "result": "refused:unknown"}) in \
        [(k, {x: kw.get(x) for x in ("room", "result")}) for k, kw in events]


def test_head_contact_is_a_person_like_the_others(acct):
    """Part 4: the head's contact needs no special path; it is routed like
    anna's, by its own person ("" the owner), and reaches only its own
    room and sessions."""
    acct.say("7", "#a1", "rechne das")           # head writes in anna's room
    assert not _prompts(acct, "anna") and not _prompts(acct, "")
    assert "messenger_cross_person_refused" in [e.kind for e in security_events.recent()]
    acct.say("7", "#o1", "status")
    assert any(c.get("session") == "s-own" for c in acct.client.commands)
