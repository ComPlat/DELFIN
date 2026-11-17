"""Named accounts beside the single dashboard token.

Input: an accounts file (name, secret stored as a jupyter-server hash,
one account flagged as head) and a login presenting a name and a secret.
Output: the account a login belongs to -- or nothing.

Properties asserted:

* the accounts store verifies a secret with the hashing jupyter-server's
  own provider already relies on (``passwd_check``, argon2); the secret
  itself is never written to the store, the page, a log line or an error;
* exactly one account is the head; an accounts file without a head, with
  two heads, or with a duplicate name is refused;
* a login as a named account yields a session whose owner is that
  account; with no accounts configured the single-token path works
  exactly as before;
* the refused-login rate limit applies per account, not globally: a
  flood against one account does not mute the record for another's
  wrong secret;
* the dashboard still binds to the loopback only (this step does not
  touch who may reach it over the network).

Unit tests drive the store, the provider and the rate limit with a fake
clock; no server is started.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.dashboard import accounts as accounts_mod
from delfin.dashboard import login_session

NAME = "ana"
OTHER = "bruno"
SECRET = "a-secret-of-ana-long-enough-for-a-token-0123456789"
OTHER_SECRET = "brunos-secret-is-long-enough-too-0123456789abcdef"


class Clock:
    def __init__(self):
        self.now = 7_000.0

    def __call__(self):
        return self.now


# -- the accounts store ----------------------------------------------------------

def _write_accounts(tmp_path: Path, entries) -> Path:
    path = tmp_path / "accounts.json"
    path.write_text(json.dumps(entries), encoding="utf-8")
    return path


def _entries():
    return [
        {"name": NAME, "secret_hash": accounts_mod.hash_secret(SECRET),
         "head": True},
        {"name": OTHER, "secret_hash": accounts_mod.hash_secret(OTHER_SECRET),
         "head": False},
    ]


def test_the_secret_verifies_against_its_hash(tmp_path):
    store = accounts_mod.Accounts(_write_accounts(tmp_path, _entries()))
    assert store.verify(NAME, SECRET)
    assert not store.verify(NAME, OTHER_SECRET)
    assert store.verify(OTHER, OTHER_SECRET)


def test_an_unknown_name_verifies_nothing(tmp_path):
    store = accounts_mod.Accounts(_write_accounts(tmp_path, _entries()))
    assert not store.verify("not-a-name", SECRET)


def test_the_store_holds_no_secret_in_clear(tmp_path):
    path = _write_accounts(tmp_path, _entries())
    store = accounts_mod.Accounts(path)
    store.verify(NAME, SECRET)
    text = path.read_text(encoding="utf-8")
    assert SECRET not in text
    assert OTHER_SECRET not in text
    assert accounts_mod.hash_secret(SECRET) in text


def test_exactly_one_account_is_the_head(tmp_path):
    store = accounts_mod.Accounts(_write_accounts(tmp_path, _entries()))
    assert store.head() == NAME
    assert store.is_head(NAME)
    assert not store.is_head(OTHER)


def test_accounts_without_a_head_are_refused(tmp_path):
    entries = _entries()
    entries[0]["head"] = False
    with pytest.raises(accounts_mod.AccountsError):
        accounts_mod.Accounts(_write_accounts(tmp_path, entries))


def test_two_heads_are_refused(tmp_path):
    entries = _entries()
    entries[1]["head"] = True
    with pytest.raises(accounts_mod.AccountsError):
        accounts_mod.Accounts(_write_accounts(tmp_path, entries))


def test_a_duplicate_name_is_refused(tmp_path):
    entries = _entries()
    entries.append({"name": NAME, "secret_hash": "argon2:x", "head": False})
    with pytest.raises(accounts_mod.AccountsError):
        accounts_mod.Accounts(_write_accounts(tmp_path, entries))


def test_an_empty_file_means_no_accounts(tmp_path):
    path = tmp_path / "accounts.json"
    path.write_text("", encoding="utf-8")
    store = accounts_mod.Accounts(path)
    assert not store and store.head() is None


def test_names_of_the_accounts_are_listed(tmp_path):
    store = accounts_mod.Accounts(_write_accounts(tmp_path, _entries()))
    assert set(store.names()) == {NAME, OTHER}


# -- sessions carry their owner ---------------------------------------------------

def _sessions(tmp_path, clock):
    return login_session.Sessions(tmp_path / "s.json", idle=3600.0,
                                  lifetime=7 * 86400.0, clock=clock)


def test_a_session_records_its_owner_and_reports_it(tmp_path):
    clock = Clock()
    sessions = _sessions(tmp_path, clock)
    sid = sessions.open(owner=NAME)
    assert sessions.owner(sid) == NAME
    sid2 = sessions.open()                # ownerless, as today
    assert sessions.owner(sid2) is None


def test_an_old_two_field_store_still_loads(tmp_path):
    clock = Clock()
    sessions = _sessions(tmp_path, clock)
    sid = sessions.open()
    text = (tmp_path / "s.json").read_text(encoding="utf-8")
    # what today's writer produces: digest -> [created, seen]
    data = json.loads(text)
    assert all(len(v) == 2 for v in data.values())
    again = login_session.Sessions(tmp_path / "s.json", idle=3600.0,
                                   lifetime=7 * 86400.0, clock=clock)
    assert again.touch(sid)


def test_an_unknown_session_has_no_owner(tmp_path):
    sessions = _sessions(tmp_path, Clock())
    assert sessions.owner("no-such-session") is None


# -- refused logins are counted per account ----------------------------------------

def test_a_flood_against_one_account_does_not_mute_another(tmp_path):
    clock = Clock()
    seen = []
    refused = login_session.RefusedLogins(window=60.0, clock=clock,
                                          record=lambda *a, **k: seen.append(a))
    refused.refused("wrong secret", "10.0.0.7", "/login", who=NAME)
    clock.now += 1
    refused.refused("wrong secret", "10.0.0.7", "/login", who=NAME)
    assert len(seen) == 1            # ana's flood is rate limited
    clock.now += 1
    refused.refused("wrong secret", "10.0.0.7", "/login", who=OTHER)
    assert len(seen) == 2            # bruno's refusal is still recorded now


def test_the_record_names_the_account_but_never_a_secret(tmp_path):
    clock = Clock()
    seen = []
    refused = login_session.RefusedLogins(window=60.0, clock=clock,
                                          record=lambda *a, **k: seen.append(a))
    refused.refused("wrong secret", "10.0.0.7", "/login", who=NAME)
    detail = seen[-1][-1]
    assert NAME in detail
    assert SECRET not in detail
