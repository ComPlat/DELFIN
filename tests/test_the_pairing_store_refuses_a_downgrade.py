"""The pairing store refuses a downgrade to format 1 and a damaged key cleanly.

Review finding 7: deleting ``store.key`` and writing a format-1 file
(``{"contact_id": ...}``) made that contact the owner's verified phone;
a key file of the wrong length made every write raise TypeError.

Output, each with its control:

* once the store has been format 2 (its ``format`` marker), a format-1
  file is refused -- with or without the key -- and nobody is paired;
* a key file of the wrong length refuses the store; a write raises
  StoreRefused, not TypeError;
* a genuine format-1 store (never migrated) is still migrated (control).
"""
from __future__ import annotations

import json

import pytest

from delfin.messenger import store


@pytest.fixture
def state(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "pairing"))
    monkeypatch.setattr(store, "_REFUSED_SEEN", set())
    return tmp_path / "pairing"


@pytest.mark.parametrize("key", ["deleted", "kept"])
def test_format_1_after_format_2_is_refused(state, key):
    store.save("7")
    if key == "deleted":
        (state / "store.key").unlink()
    (state / "pair.json").write_text(json.dumps({"contact_id": "666"}))
    assert store.refused()
    assert store.paired_contact() is None and store.verified_person("666") is None
    with pytest.raises(store.StoreRefused):
        store.save("8")


def test_a_key_of_the_wrong_length_refuses_cleanly(state):
    store.save("7")
    (state / "store.key").write_bytes(b"short")
    assert store.refused() and store.paired_contact() is None
    with pytest.raises(store.StoreRefused):
        store.save("8")


def test_a_genuine_format_1_store_is_still_migrated_control(state):
    state.mkdir(parents=True)
    (state / "pair.json").write_text(json.dumps({"contact_id": "7"}))
    assert store.paired_contact() == "7"
    store.save("8", person="anna")
    assert store.verified_person("7") == "" and (state / "format").read_text().strip() == "2"
