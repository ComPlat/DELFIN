"""A per-person credential store beside the shared one.

The Raspberry Pi is to serve several people, each bringing their own API key.
Before this wave the store is single-tenant: one ``~/.delfin/credentials.json``
belongs to the one OS user every service runs as, so one store serves
everybody. This step adds an optional ``person`` to the store API:

- ``load_credential(name, person="")``: with a person, look in that person's
  store first and fall back to the shared store; with no person the lookup
  order is exactly today's (env > keyring > file > "").
- ``set_credential`` / ``delete_credential`` / ``list_credentials`` take the
  same optional ``person``; ``list_credentials`` for a person never reveals
  that another person holds a credential of the same name.
- One person's store is never read as another's through this module's own
  API: a wrong person name gets the shared fallback or nothing.
- A missing credential still returns ``""`` — a value that names neither a
  credential nor a person (the module's docstring already forbids callers
  from revealing the name; this keeps it true).

The module enforces no filesystem boundary here: who may read which store is
the egress proxy and the workspace roots in later steps. This file pins the
API contract only.
"""

from __future__ import annotations

import json
import sys
import types

import pytest

from delfin.agent import credentials as cred


@pytest.fixture
def store_path(tmp_path, monkeypatch):
    p = tmp_path / "credentials.json"
    monkeypatch.setattr(cred, "_DEFAULT_PATH", p)
    return p


# ---------------------------------------------------------------------------
# task 1 — no-person lookup order stays exactly today's, all four sources
# ---------------------------------------------------------------------------


def _install_fake_keyring(monkeypatch, store_path):
    """A fake keyring backend so the keyring source can be driven without a
    real OS keyring. Modeled on the backend fakes in
    test_credentials_use_the_os_keyring.py."""
    try:
        import keyring.backend  # noqa: F401
        kr = sys.modules["keyring"]
    except ImportError:
        kr = _make_keyring_module(monkeypatch)
        keyring_core = kr
    else:
        keyring_core = sys.modules["keyring"]

    class Fake(keyring_core.backend.KeyringBackend):
        priority = 10

        def __init__(self):
            self.items = {}

        def get_password(self, service, username):
            return self.items.get((service, username))

        def set_password(self, service, username, password):
            self.items[(service, username)] = password

        def delete_password(self, service, username):
            if (service, username) not in self.items:
                raise keyring_core.errors.PasswordDeleteError(username)
            del self.items[(service, username)]

    Fake.__name__ = Fake.__qualname__ = "FakeSecretService"
    inst = Fake()
    monkeypatch.setattr(keyring_core, "get_keyring", lambda: inst, raising=False)
    monkeypatch.setattr(cred, "_keyring_factory", cred._default_keyring_factory)
    cred.reset_keyring_state()
    return inst


def _make_keyring_module(monkeypatch):
    pkg = types.ModuleType("keyring")
    backend = types.ModuleType("keyring.backend")
    errors = types.ModuleType("keyring.errors")

    class KeyringBackend:
        priority = 1

        def get_password(self, service, username):
            raise NotImplementedError

        def set_password(self, service, username, password):
            raise NotImplementedError

        def delete_password(self, service, username):
            raise NotImplementedError

    class PasswordDeleteError(Exception):
        pass

    backend.KeyringBackend = KeyringBackend
    errors.PasswordDeleteError = PasswordDeleteError
    pkg.backend = backend
    pkg.errors = errors
    pkg.get_keyring = lambda: None
    monkeypatch.setitem(sys.modules, "keyring", pkg)
    monkeypatch.setitem(sys.modules, "keyring.backend", backend)
    monkeypatch.setitem(sys.modules, "keyring.errors", errors)
    return pkg


def test_no_person_order_is_env_then_keyring_then_file_then_empty(
        store_path, monkeypatch):
    """With no person the lookup order is EXACTLY today's: env > keyring >
    file > ''. Driving all four sources pins it."""
    kr = _install_fake_keyring(monkeypatch, store_path)
    monkeypatch.delenv("M4_KEY", raising=False)

    # 1. env wins over everything.
    monkeypatch.setenv("M4_KEY", "from-env")
    assert load("M4_KEY", store_path, kr, person="") == "from-env"

    # 2. keyring wins over the file.
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.reset_keyring_state()
    put_keyring(kr, "M4_KEY", "from-keyring")
    put_file(store_path, "M4_KEY", "from-file")
    assert load("M4_KEY", store_path, kr, person="") == "from-keyring"

    # 3. file when only it holds the value.
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.reset_keyring_state()
    del_keyring(kr, "M4_KEY", store_path)
    put_file(store_path, "M4_KEY", "from-file")
    assert load("M4_KEY", store_path, kr, person="") == "from-file"

    # 4. '' when nothing holds it.
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.reset_keyring_state()
    del_keyring(kr, "M4_KEY", store_path)
    clear_file(store_path)
    assert load("M4_KEY", store_path, kr, person="") == ""


def put_file(store_path, name, value):
    data = {}
    if store_path.exists():
        data = json.loads(store_path.read_text(encoding="utf-8"))
    data[name] = value
    store_path.write_text(json.dumps(data), encoding="utf-8")


def clear_file(store_path):
    store_path.write_text("{}", encoding="utf-8")


def put_keyring(kr, name, value):
    kr.set_password("delfin", name, value)


def del_keyring(kr, name, store_path):
    from pathlib import Path
    idx = Path(store_path.with_name(store_path.stem + ".keyring.json"))
    names = []
    if idx.exists():
        names = json.loads(idx.read_text(encoding="utf-8")).get("names", [])
    if name not in names:
        idx.write_text(json.dumps({"names": [name]}), encoding="utf-8")
    try:
        kr.delete_password("delfin", name)
    except Exception:
        pass


def load(name, store_path, kr, person=""):
    return cred.load_credential(name, path=store_path, person=person)


def test_single_user_default_unchanged(store_path, monkeypatch):
    """Whether person is implemented of not, the no-person path must behave
    as before: an env var wins, else the file, else ''."""
    cred.set_credential("M4_KEY", "file-value", path=store_path)
    monkeypatch.setenv("M4_KEY", "env-value")
    assert cred.load_credential("M4_KEY", path=store_path) == "env-value"
    monkeypatch.delenv("M4_KEY", raising=False)
    assert cred.load_credential("M4_KEY", path=store_path) == "file-value"
    assert cred.load_credential("NEVER", path=store_path) == ""


# ---------------------------------------------------------------------------
# task 1 — a person's store first, then shared fallback
# ---------------------------------------------------------------------------


def test_person_store_read_before_shared(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.set_credential("M4_KEY", "shared-value", path=store_path)
    cred.set_credential("M4_KEY", "person-value", path=store_path, person="alice")
    # The person's own store wins over the shared fallback.
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") \
        == "person-value"
    # The shared store is untouched by the person write.
    assert cred.load_credential("M4_KEY", path=store_path) == "shared-value"


def test_person_falls_back_to_shared(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    # Only the shared store holds it; the person store does not exist.
    cred.set_credential("M4_KEY", "shared-value", path=store_path)
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") \
        == "shared-value"


def test_person_env_still_wins(store_path, monkeypatch):
    cred.set_credential("M4_KEY", "shared-value", path=store_path)
    cred.set_credential("M4_KEY", "person-value", path=store_path, person="alice")
    monkeypatch.setenv("M4_KEY", "env-value")
    # env is the session-level override; it wins even for a person.
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") \
        == "env-value"


# ---------------------------------------------------------------------------
# task 2 — set / delete / list take a person; list never reveals another's
# ---------------------------------------------------------------------------


def test_set_person_store_is_separate_from_shared(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.set_credential("M4_KEY", "alice-value", path=store_path, person="alice")
    # visible to alice, not through the shared/no-person view
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") \
        == "alice-value"
    assert cred.load_credential("M4_KEY", path=store_path) == ""


def test_two_people_do_not_see_each_other(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.set_credential("M4_KEY", "alice-value", path=store_path,
                        person="alice")
    cred.set_credential("M4_KEY", "bob-value", path=store_path, person="bob")
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") \
        == "alice-value"
    assert cred.load_credential("M4_KEY", path=store_path, person="bob") \
        == "bob-value"


def test_delete_person_only_deletes_that_person(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.set_credential("M4_KEY", "shared-value", path=store_path)
    cred.set_credential("M4_KEY", "alice-value", path=store_path,
                        person="alice")
    assert cred.delete_credential("M4_KEY", path=store_path,
                                  person="alice") is True
    # Alice no longer holds it; the shared value survives.
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") \
        == "shared-value"
    assert cred.load_credential("M4_KEY", path=store_path) == "shared-value"


def test_delete_person_missing_returns_false(store_path):
    assert cred.delete_credential("NEVER", path=store_path,
                                  person="alice") is False


def test_list_credentials_for_a_person_shows_their_own_store(
        store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    cred.set_credential("M4_KEY", "alice-value", path=store_path,
                        person="alice")
    items = cred.list_credentials(path=store_path, person="alice")
    assert "M4_KEY" in items
    assert items["M4_KEY"]["source"] == "file"


def test_list_person_does_not_reveal_another_person_same_name(
        store_path, monkeypatch):
    """list_credentials(person="bob") must not reveal that alice holds a
    credential called M4_KEY."""
    monkeypatch.delenv("M4_KEY", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    cred.set_credential("M4_KEY", "alice-value", path=store_path,
                        person="alice")
    items = cred.list_credentials(path=store_path, person="bob")
    assert "M4_KEY" not in items


def test_list_person_does_not_reveal_shared_same_name(store_path, monkeypatch):
    """A person's own list is scoped to their store + the session env; a
    shared-store value with the same name is not surfaced as theirs."""
    monkeypatch.delenv("M4_KEY", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    cred.set_credential("M4_KEY", "shared-value", path=store_path)
    items = cred.list_credentials(path=store_path, person="bob")
    assert "M4_KEY" not in items


# ---------------------------------------------------------------------------
# task 3 — a wrong person name gets shared fallback or nothing
# ---------------------------------------------------------------------------


def test_wrong_person_gets_shared_fallback(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    cred.set_credential("M4_KEY", "shared-value", path=store_path)
    # 'nobody' is not a person with a store; the shared value is the fallback.
    assert cred.load_credential("M4_KEY", path=store_path, person="nobody") \
        == "shared-value"


def test_wrong_person_gets_nothing_not_anothers_value(store_path, monkeypatch):
    monkeypatch.delenv("M4_KEY", raising=False)
    # Only alice holds M4_KEY; naming the wrong person must yield ''
    # (the shared store has nothing), never alice's value.
    cred.set_credential("M4_KEY", "alice-value", path=store_path,
                        person="alice")
    assert cred.load_credential("M4_KEY", path=store_path, person="wrong") == ""


# ---------------------------------------------------------------------------
# task 4 — a missing credential yields a refusal that names no name/person
# ---------------------------------------------------------------------------


def test_missing_credential_yields_even_when_no_person_store(store_path):
    """With a person named but no store of theirs and no shared value, the
    result is exactly '' — a refusal that carries neither a credential name
    nor a person name, and never raises."""
    assert cred.load_credential("M4_KEY", path=store_path, person="alice") == ""


def test_missing_credential_yields_no_person_in_no_person_call(store_path):
    assert cred.load_credential("M4_KEY", path=store_path) == ""
