"""``hol`` stays, bounded: own tree, a per-person switch, announced, not from a new phone.

Review finding 9 (design decision): pulling a running session to the phone
stays easy, within these bounds -- each through the real gateway and
bridge (tests/_messenger_people_rig.py), each with its control:

* only a session of the person that runs in that person's own workspace
  tree (``~/people/<name>``; the owner: the account home outside
  ``~/people``);
* a per-person setting, default on: ``/sync policy pull off`` at a window
  stops every pull of that person's phones, ``/sync policy pull on``
  allows it again;
* every pull is announced in the session and to the person's other phones;
* a phone verified less than 10 minutes ago cannot pull.
"""
from __future__ import annotations

import time

import pytest

from _messenger_people_rig import Account
from delfin.messenger import activation, store


@pytest.fixture(autouse=True)
def _german_phone(monkeypatch):
    """These tests read the German phone texts (lang.DE); English is the
    default (test_the_phone_speaks_the_persons_language)."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")

OLD = time.time() - 3600


@pytest.fixture
def acct(tmp_path, monkeypatch):
    account = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone", verified_at=OLD)
    store.save("22", person="anna", label="anna tablet", verified_at=OLD)
    return account


def _on(acct, sid):
    return activation.enabled(sid, directory=acct.activation_dir("anna"))


def test_a_session_outside_the_persons_tree_is_not_pulled(acct, tmp_path):
    outside = tmp_path / "elsewhere" / "proj"
    outside.mkdir(parents=True)
    acct.services["anna"].sessions.append({"sessionId": "s-out", "title": "Draussen",
                                           "workspace": str(outside), "state": "idle"})
    acct.say("21", "@21", "hol Draussen")
    assert not _on(acct, "s-out")
    assert "Ordner" in acct.daemon.to("@21")[-1]
    acct.say("21", "@21", "hol Butan")                      # control
    assert _on(acct, "s-anna")


def test_the_person_can_switch_pulling_off_and_on_at_a_window(acct):
    out = activation.configure("s-anna", "policy pull off", directory=acct.activation_dir("anna"))
    assert "off" in out and not activation.pull_allowed(directory=acct.activation_dir("anna"))
    acct.say("21", "@21", "hol Butan")
    assert not _on(acct, "s-anna")
    activation.configure("s-anna", "policy pull on", directory=acct.activation_dir("anna"))
    acct.say("21", "@21", "hol Butan")
    assert _on(acct, "s-anna")


def test_every_pull_is_announced_to_the_other_phones_and_the_session(acct):
    acct.say("21", "@21", "hol Butan")
    assert any("Butan" in t and "anna phone" in t for t in acct.daemon.to("@22"))
    assert {"sessionId": "s-anna", "text": "/sync status"} in acct.services["anna"].calls("session.prompt")


def test_a_freshly_verified_phone_cannot_pull(acct):
    store.save("23", person="anna", label="new phone", verified_at=time.time() - 60)
    acct.say("23", "@23", "hol Butan")
    assert not _on(acct, "s-anna")
    assert "10 Minuten" in acct.daemon.to("@23")[-1]
    acct.say("21", "@21", "hol Butan")                      # control: an older phone
    assert _on(acct, "s-anna")
