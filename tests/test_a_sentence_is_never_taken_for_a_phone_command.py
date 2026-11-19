"""A sentence that starts with a command word goes to the agent, not to the command.

Review finding 4: ``get the energies from step 3``, ``hol die Energien aus
Schritt 3``, ``phones are fine, continue``, ``help me fix ...`` were taken
for the phone commands get, phones and help and never reached the session.

Output, through the real gateway and bridge (tests/_messenger_people_rig.py):

* the commands match only as a whole short message: ``sessions``,
  ``hilfe``/``help``, ``handys``/``phones``, ``weg``/``off``,
  ``entkoppeln``/``unpair``, and ``hol``/``get`` with exactly one token;
* in a session's chat every longer line is a prompt for that session, in
  German and English alike (ten sentences); ``hol NAME`` that names no
  session falls through to the session as a prompt;
* controls: ``hol 2``, ``sessions`` and ``hilfe`` still work.
"""
from __future__ import annotations

import pytest

from _messenger_people_rig import Account
from delfin.messenger import store

SENTENCES = [
    "get the energies from step 3",
    "hol die Energien aus Schritt 3",
    "phones are fine, continue",
    "help me fix the input file",
    "Hilfe, das Skript bricht bei Schritt 2 ab",
    "sessions of the solvent model take longer, check it",
    "weg mit den alten Dateien in scratch",
    "off by one in the loop over atoms",
    "hol mir bitte die Ergebnisse",
    "entkoppeln der Schwingungen bitte prüfen",
]


@pytest.fixture
def acct(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    account.enable("anna", "s-anna")
    store.put_room("#a1", "s-anna", "Butan", person="anna", ready=True)
    return account


def _prompts(acct):
    return [p["text"].split("\n\n[Messenger")[0] for p in acct.services["anna"].calls("session.prompt")
            if p["text"] != "/sync status"]


@pytest.mark.parametrize("sentence", SENTENCES)
def test_a_sentence_reaches_the_session(acct, sentence):
    acct.say("21", "#a1", sentence)
    assert _prompts(acct) == [sentence]
    assert store.contact("21") is not None                 # not unpaired


def test_a_get_that_names_no_session_falls_through_to_the_agent(acct):
    acct.say("21", "#a1", "hol Energien")
    assert _prompts(acct) == ["hol Energien"]


def test_the_short_commands_still_work(acct):
    acct.say("21", "@21", "sessions")
    assert "1. Butan" in acct.daemon.to("@21")[-1]
    acct.say("21", "@21", "hol 2")
    assert "Ethan" in acct.daemon.to("@21")[-1]
    acct.say("21", "@21", "hilfe")
    assert "get N" in acct.daemon.to("@21")[-1]          # the help, in English by default
    assert not _prompts(acct)
