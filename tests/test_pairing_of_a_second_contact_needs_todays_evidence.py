"""Pairing a second contact keeps today's evidence, and the single-user
path keeps working.

Input: the pairing store (delfin/messenger/store.py format 2) in a
temporary state dir; for the full second-person scenario, the REAL
gateway and pairing over the people rig (tests/_messenger_people_rig.py).
Before: exactly one contact paired, as format 1 had.

Output:

* a SECOND contact is bound only through the same evidence a first one
  needs: a one-time invitation of that person's, consumed once, then the
  6-digit code confirmed in time by the session that started the pairing;
  per-person limits (MAX_PHONES, invitations per window) apply to the
  second person exactly as to the first;
* an invitation can be issued only for a real person of the account --
  the head's pairing is never a shortcut for anyone else;
* today's single-user path is unchanged: the owner's first verified
  contact is paired_contact() as before, and it stays so beside anna's.
"""
from __future__ import annotations

import time

import pytest

from _messenger_people_rig import Account
from conftest import truthful_peer
from delfin.agent import security_events
from delfin.messenger import pairing, store
from delfin.messenger.people import OWNER


@pytest.fixture(autouse=True)
def _german_phone(monkeypatch):
    """These tests read the German phone texts (lang.DE); English is the
    default (test_the_phone_speaks_the_persons_language)."""
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "de")


@pytest.fixture
def state(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "messenger"))
    store._REFUSED_SEEN.clear()
    yield tmp_path / "messenger"
    store._REFUSED_SEEN.clear()


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    account = Account(tmp_path, monkeypatch)
    yield account
    security_events.clear()


def _pair(acct, person="anna", requester="session:s-anna"):
    rid = pairing.begin(requester, directory=acct.pairing_dir(person))
    acct.round()
    return rid, pairing.read(rid, requester, directory=acct.pairing_dir(person))


def _code(acct, cid):
    text = acct.daemon.to("@" + cid)[0]
    return text.split("Code ")[1].split(".")[0]


def test_a_second_contact_needs_the_same_invitation_and_code(state):
    """The head's phone is paired; anna's second contact is bound ONLY by
    anna's own one-time invitation plus her code -- nothing about being
    second makes it easier."""
    store.save("7")                              # the head's phone, today's path
    store.add_invitation("41", person="anna", request="r" * 32)
    # no code yet: bound but without rights
    store.bind("21", person="anna", label="anna phone", request="r" * 32,
               code="123456")
    assert store.verified_person("21") is None
    assert store.open_invitation("anna", "r" * 32) is not None
    # a wrong request code consumes nothing
    got, rec = store.check_code("21", "000000")
    assert got == "wrong" and store.verified_person("21") is None
    # the right code verifies; the invitation was the evidence, and it is
    # used once (a replay is refused)
    conn, _ = store.open_invitation("anna", "r" * 32)
    why = store.consume_invitation(conn, now=time.time())[1]
    assert why == "ok"
    assert store.consume_invitation(conn, now=time.time())[1] == "replayed"
    got, rec = store.check_code("21", "123456")
    assert got == "verified" and store.verified_person("21") == "anna"
    # and today's single-user answer is untouched by the second contact
    assert store.paired_contact() == "7"
    assert store.verified_person("7") == OWNER
    assert store.persons() == [OWNER, "anna"]


def test_an_invitation_for_a_malformed_person_is_refused(state):
    """The head's verified contact does not authorise anything for anyone
    else: an invitation is issued only for the owner or a well-formed
    person name (people.valid); the account-existence check happens when
    the request is collected (people.home). The head's pairing is never a
    format-level shortcut."""
    store.save("7")
    for bad in ("Ghost", "../anna", "-x", None):
        with pytest.raises(store.PairingRefused):
            store.add_invitation("41", person=bad, request="r" * 32)


def test_the_second_person_hits_the_same_limits(state):
    store.save("7")
    for i in range(pairing.MAX_PHONES):
        store.save(str(100 + i), person="anna")
    with pytest.raises(store.PairingRefused):
        store.add_invitation("50", person="anna", request="r" * 32)
    # the control: another person is not affected
    store.add_invitation("51", person="ben", request="s" * 32)
    assert store.open_invitation("ben", "s" * 32) is not None


def test_single_user_path_unchanged(state):
    """Today's installation: one contact, no people mode in the store --
    paired_contact(), verified_person and persons() answer as before."""
    store.save("7", label="phone")
    assert store.paired_contact() == "7"
    assert store.verified_person("7") == OWNER
    assert store.phones(OWNER, verified_only=True) == [("7", store.contact("7"))]
    assert store.persons() == [OWNER]


def test_a_second_person_pairs_through_the_real_path_while_the_head_is_paired(acct, truthful_peer):
    """Part 3, the whole scenario: the head's phone is already verified,
    then anna pairs her own phone through the SAME path (her /pair in her
    own directory, her one-time invitation, her code confirmed in her
    session). No shortcut, no sharing of evidence; the head stays the
    owner contact, anna is a person beside it."""
    store.save("7", label="head phone")
    rid, res = _pair(acct, "anna")
    assert res["state"] == "offered"
    acct.connect("21", "901", name="Anna Pixel")
    rec = store.contact("21")
    assert rec["person"] == "anna" and rec["state"] == "unverified"
    code = _code(acct, "21")
    # anna's session confirms with the code; the head's session confirms nothing
    assert pairing.confirm("session:s-other", code,
                           directory=acct.pairing_dir("anna")) is None
    assert pairing.confirm("session:s-anna", code,
                           directory=acct.pairing_dir("anna")) is not None
    acct.round()
    assert store.verified_person("21") == "anna"
    # the head's contact is untouched, still the owner's
    assert store.paired_contact() == "7"
    assert store.verified_person("7") == OWNER
    assert store.persons() == [OWNER, "anna"]
