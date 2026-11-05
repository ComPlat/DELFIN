"""The gateway's pairing store binds each phone to exactly one person.

Input: delfin/messenger/store.py, format 2 (contacts -> person, rooms per
person, one-time invitations, an HMAC under store.key). Output, each with
its control:

* a format-1 store (one ``contact_id``) is read and migrated without
  re-pairing: its contact becomes the owner's verified phone, its rooms the
  owner's;
* a tampered mapping, a corrupt file and a format-1 file written after the
  key exists (a downgrade) refuse the WHOLE store (nobody is paired, writes
  raise) and record one ``messenger_store_refused`` event;
* an invitation is used once; a second use is ``replayed``, a late one
  ``expired``; two processes racing for one invitation bind once;
* a contact bound to one person is never rebound in place;
* at most MAX_PHONES phones and INVITES_PER_WINDOW invitations per person;
* every file is 0600 in a 0700 directory.
"""
from __future__ import annotations

import json
import os
import stat
import threading

import pytest

from delfin.agent import security_events
from delfin.messenger import pairing, store


@pytest.fixture
def state(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "pairing"))
    monkeypatch.setattr(store, "_REFUSED_SEEN", set())
    security_events.clear()
    yield tmp_path / "pairing"
    security_events.clear()


def _write_raw(path, data):
    path.mkdir(parents=True, exist_ok=True)
    (path / "pair.json").write_text(json.dumps(data))


def test_a_format_1_store_migrates_to_the_owner_without_re_pairing(state):
    _write_raw(state, {"contact_id": "7",
                       "rooms": {"#10": {"session_id": "s1", "title": "t", "ready": True}},
                       "schedules_room": "#40", "schedule_seq": 3})
    # read as it is: the owner's verified phone, the owner's room
    assert store.paired_contact() == "7"
    assert store.verified_person("7") == ""
    assert store.rooms()["#10"]["person"] == ""
    assert store.schedules_room() == "#40" and store.schedule_seq() == 3
    # the first write makes it format 2 and keeps everything
    store.put_room("#20", "s2", "u", owner="7", ready=True)
    data = json.loads((state / "pair.json").read_text())
    assert data["version"] == 2 and "mac" in data
    assert data["contacts"]["7"]["state"] == "verified"
    assert set(store.rooms()) == {"#10", "#20"}
    assert not store.refused()


def test_an_untampered_store_is_read_control(state):
    store.save("7")
    store.save("8", person="anna")
    assert store.verified_person("8") == "anna" and not store.refused()


def test_a_tampered_mapping_refuses_the_whole_store_and_is_audited(state):
    store.save("7")
    store.save("8", person="anna")
    data = json.loads((state / "pair.json").read_text())
    data["contacts"]["8"]["person"] = ""          # anna's phone made the owner's
    (state / "pair.json").write_text(json.dumps(data))
    assert store.refused()
    assert store.verified_person("8") is None and store.verified_person("7") is None
    assert store.paired_contact() is None and store.persons() == []
    with pytest.raises(store.StoreRefused):
        store.save("9")
    kinds = [e.kind for e in security_events.recent()]
    assert kinds.count("messenger_store_refused") == 1


@pytest.mark.parametrize("damage", ["not json", "downgrade", "no key"])
def test_a_damaged_store_fails_closed(state, damage):
    store.save("7")
    if damage == "not json":
        (state / "pair.json").write_text("{nope")
    elif damage == "downgrade":
        (state / "pair.json").write_text(json.dumps({"contact_id": "666"}))
    else:
        (state / "store.key").unlink()
    assert store.refused()
    assert store.paired_contact() is None and store.contact("666") is None


def test_an_invitation_is_used_once(state):
    store.add_invitation("41", person="anna", request="r" * 32, now=1000.0)
    inv, why = store.consume_invitation("41", now=1001.0)
    assert why == "ok" and inv["person"] == "anna"
    again, why = store.consume_invitation("41", now=1002.0)
    assert again is None and why == "replayed"


def test_an_invitation_past_its_time_is_expired(state):
    store.add_invitation("42", person="anna", request="r" * 32, now=1000.0)
    inv, why = store.consume_invitation("42", now=1000.0 + pairing.INVITATION_TTL_S + 1)
    assert inv is None and why == "expired"
    # the control: in time it works
    store.add_invitation("43", person="anna", request="s" * 32, now=5000.0)
    assert store.consume_invitation("43", now=5001.0)[1] == "ok"


def test_the_admin_invitation_is_not_the_gateways(state):
    store.add_invitation("44", person="", request="admin", by="admin")
    assert store.consume_invitation("44", by="gateway") == (None, "other")
    assert store.consume_invitation("44", by="admin")[1] == "ok"
    # used by the admin: still "other" for the gateway, not a replay
    assert store.consume_invitation("44", by="gateway") == (None, "other")


def test_two_processes_racing_for_one_invitation_bind_once(state):
    store.add_invitation("45", person="anna", request="r" * 32)
    results, barrier = [], threading.Barrier(8)

    def race():
        barrier.wait()
        results.append(store.consume_invitation("45")[1])

    threads = [threading.Thread(target=race) for _ in range(8)]
    for th in threads:
        th.start()
    for th in threads:
        th.join()
    assert sorted(results) == ["ok"] + ["replayed"] * 7


def test_a_bound_contact_is_never_rebound(state):
    store.bind("9", person="anna", label="phone", request="r" * 32, code="123456")
    with pytest.raises(store.PairingRefused):
        store.bind("9", person="ben", label="phone", request="s" * 32, code="654321")
    with pytest.raises(store.PairingRefused):
        store.save("9", person="ben")
    assert store.contact("9")["person"] == "anna"
    # rebinding goes through a revocation, which is recorded
    entry = store.revoke_contact("9", by="test")
    assert entry["person"] == "anna" and store.uncleaned()[0]["contact"] == "9"
    store.bind("9", person="ben", label="phone", request="s" * 32, code="654321")
    assert store.contact("9")["person"] == "ben"


def test_a_phone_is_verified_only_by_its_code(state):
    store.bind("9", person="anna", label="phone", request="r" * 32, code="123456")
    assert store.verified_person("9") is None             # unverified: no rights
    assert store.check_code("9", "000000")[0] == "wrong"
    assert store.check_code("9", "123 456")[0] == "verified"
    assert store.verified_person("9") == "anna"
    assert "code" not in store.contact("9")               # the hash is gone too


def test_three_wrong_codes_remove_the_phone(state):
    store.bind("9", person="anna", label="phone", request="r" * 32, code="123456")
    assert [store.check_code("9", "1")[0] for _ in range(3)] == ["wrong", "wrong", "failed"]
    assert store.contact("9") is None
    assert store.check_code("9", "123456")[0] == "none"


def test_a_late_code_is_refused(state):
    store.bind("9", person="anna", label="phone", request="r" * 32, code="123456", now=1000.0)
    late = 1000.0 + pairing.VERIFY_TTL_S + 1
    assert store.check_code("9", "123456", now=late)[0] == "expired"
    assert store.contact("9") is None


def test_phones_and_invitations_per_person_are_limited(state):
    for i in range(pairing.MAX_PHONES):
        store.save(str(100 + i), person="anna")
    with pytest.raises(store.PairingRefused) as limit:
        store.add_invitation("50", person="anna", request="r" * 32)
    assert limit.value.reason == "limit"
    # the control: another person is not affected; then the rate limit
    for i in range(pairing.INVITES_PER_WINDOW):
        store.add_invitation(str(60 + i), person="ben", request=f"{i:032d}", now=2000.0 + i)
    with pytest.raises(store.PairingRefused) as rate:
        store.add_invitation("70", person="ben", request="x" * 32, now=2010.0)
    assert rate.value.reason == "rate"
    store.add_invitation("71", person="ben", request="y" * 32,
                         now=2000.0 + pairing.INVITE_WINDOW_S + 5)


def test_a_room_belongs_to_one_person(state):
    store.save("7")
    store.save("8", person="anna")
    store.put_room("#10", "s1", "t", person="anna", ready=True)
    with pytest.raises(ValueError):
        store.put_room("#10", "s1", "t", person="", ready=True)        # another person
    with pytest.raises(ValueError):
        store.put_room("#11", "s2", "t", owner="7", person="anna")     # owner's phone, anna's room
    with pytest.raises(ValueError):
        store.configure_room("#10", owner="7", forum=False)
    assert store.room_of_session("s1", "anna") == "#10"
    assert store.room_of_session("s1", "") is None


def test_the_store_files_are_private(state):
    store.save("7")
    assert stat.S_IMODE(os.stat(state).st_mode) == 0o700
    for name in ("pair.json", "store.key"):
        assert stat.S_IMODE(os.stat(state / name).st_mode) == 0o600, name
