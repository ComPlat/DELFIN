"""A credential store that cannot be read is moved aside, never replaced.

Measured on 2026-11-17, before this behaviour existed: a store truncated
the way a full disk or a kill mid-write truncates it, plus ONE
`set_credential` call, left a file holding only the key just written. The
two keys it had held were gone and there was no copy anywhere.

The cause was one shared read: `_read_store` answered "there is nothing
stored" for both "the file is absent" and "the file cannot be read", and
every writer built its new store from that answer. Reading leniently is
right for a LOOKUP -- a damaged store must not stop the agent from using
a key out of the environment -- and wrong for a WRITE.

So `_read_store(strict=True)` raises `CredentialStoreUnreadable`, and:

* `set_credential` quarantines the unreadable file (a rename, so the bytes
  are not rewritten and the 0600 mode travels with it), stores the key,
  and then reports -- the key is not lost and the person does not type it
  twice;
* `delete_credential` refuses and leaves the file exactly as it is, since
  a store written from nothing drops every other name;
* `migrate_to_keyring` refuses rather than reporting every name "absent",
  which would read as "there was nothing to move";
* `list_credentials` reports the damage instead of an empty store, because
  "no credential" and "the file is broken" lead a person to opposite
  actions.
"""

from __future__ import annotations

import json
import os
import stat

import pytest

from delfin.agent import credentials as C

TRUNCATED = '{"KIT_TOOLBOX_API_KEY": "STAND-IN-KIT", "ANTHRO'


@pytest.fixture
def store(tmp_path, monkeypatch):
    """A private store holding two stand-in values, at a path of our own."""
    for name in ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY",
                 "ANTHROPIC_API_KEY"):
        monkeypatch.delenv(name, raising=False)
    d = tmp_path / ".delfin"
    d.mkdir()
    os.chmod(d, 0o700)
    p = d / "credentials.json"
    p.write_text(json.dumps({"KIT_TOOLBOX_API_KEY": "STAND-IN-KIT",
                             "ANTHROPIC_API_KEY": "STAND-IN-ANT"}) + "\n")
    os.chmod(p, 0o600)
    return p


def _aside(path):
    return sorted(q for q in path.parent.iterdir()
                  if "unreadable" in q.name)


def test_a_lookup_still_answers_from_a_damaged_store(store):
    """Lenient on purpose: the agent must not be stopped by a broken file."""
    store.write_text(TRUNCATED)
    assert C.load_credential("KIT_TOOLBOX_API_KEY", path=store) == ""


def test_setting_a_key_moves_the_damaged_store_aside_and_keeps_it(store):
    store.write_text(TRUNCATED)
    os.chmod(store, 0o600)

    with pytest.raises(C.CredentialStoreUnreadable) as exc:
        C.set_credential("OPENAI_API_KEY", "STAND-IN-NEW", path=store,
                         store="file")

    # The key asked for is stored: the report is not a failed write.
    assert C.load_credential("OPENAI_API_KEY", path=store) == "STAND-IN-NEW"
    moved = _aside(store)
    assert len(moved) == 1, "the unreadable store must be kept, not deleted"
    assert "STAND-IN-KIT" in moved[0].read_text(), (
        "the quarantined file must still hold what was in it")
    assert stat.S_IMODE(moved[0].stat().st_mode) == 0o600, (
        "a rename keeps the mode; the secrets stay as private as before")
    # The message names the file, never a value.
    assert "STAND-IN-KIT" not in str(exc.value)
    assert "unreadable" in str(exc.value)


def test_a_second_damage_does_not_overwrite_the_first_rescue(store):
    """One `.bak` slot is how the original is lost. The name is dated."""
    store.write_text(TRUNCATED)
    with pytest.raises(C.CredentialStoreUnreadable):
        C.set_credential("OPENAI_API_KEY", "STAND-IN-1", path=store,
                         store="file")
    store.write_text('{"OTHER_KEY": "STAND-IN-OTHER", "BROK')
    with pytest.raises(C.CredentialStoreUnreadable):
        C.set_credential("OPENAI_API_KEY", "STAND-IN-2", path=store,
                         store="file")
    moved = _aside(store)
    assert len(moved) == 2, "both damaged stores are kept"
    kept = "".join(q.read_text() for q in moved)
    assert "STAND-IN-KIT" in kept and "STAND-IN-OTHER" in kept


def test_deleting_from_a_damaged_store_refuses_and_changes_nothing(store):
    store.write_text(TRUNCATED)
    before = store.read_text()
    with pytest.raises(C.CredentialStoreUnreadable):
        C.delete_credential("KIT_TOOLBOX_API_KEY", path=store)
    assert store.read_text() == before, "the file must be untouched"
    assert not _aside(store), "a delete does not quarantine on its own"


def test_migrating_a_damaged_store_refuses_instead_of_saying_absent(
        store, monkeypatch):
    monkeypatch.setattr(C, "configured_store", lambda: "auto")
    monkeypatch.setattr(C, "_usable_backend", lambda: None)
    store.write_text(TRUNCATED)
    with pytest.raises(C.CredentialStoreUnreadable):
        C.migrate_to_keyring(path=store)


def test_the_listing_says_the_store_is_damaged_not_that_it_is_empty(store):
    store.write_text(TRUNCATED)
    items = C.list_credentials(path=store)
    assert items.get("", {}).get("source") == "unreadable"
    assert "detail" in items[""], "a person needs the reason"
    assert items[""]["value"] == "", "no part of a value is reported"


def test_an_intact_store_is_untouched_by_any_of_this(store):
    """The control: nothing above may change the normal path."""
    assert C.set_credential("OPENAI_API_KEY", "STAND-IN-NEW", path=store,
                            store="file") is True
    kept = json.loads(store.read_text())
    assert sorted(kept) == ["ANTHROPIC_API_KEY", "KIT_TOOLBOX_API_KEY",
                            "OPENAI_API_KEY"], (
        "the other keys must survive an ordinary write")
    assert not _aside(store), "nothing is quarantined when the store parses"
    assert C.delete_credential("ANTHROPIC_API_KEY", path=store) is True
    assert sorted(json.loads(store.read_text())) == [
        "KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY"]
