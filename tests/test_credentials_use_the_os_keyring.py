"""API keys go into the operating system's credential store when one is usable.

Input: the setting ``agent.credentials.store`` (auto | keyring | file, read
from the user's own settings file only), the ``keyring`` package when it is
installed, and the backend it selects. Output: where ``set_credential``
writes, where ``load_credential`` reads, what ``credentials list`` shows and
what ``credentials migrate`` moves.

Every backend here is a fake: a ``keyring.backend.KeyringBackend`` subclass
handed to the real factory (``keyring.get_keyring``). Where the package is
not installed an in-memory stand-in module takes its place, so the same
cases run on every host. HOME and every store live under tmp_path; no test
reaches the real keyring, the real ``~/.delfin`` or the real settings file.
"""
from __future__ import annotations

import json
import subprocess
import sys
import threading
import time
import types

import pytest

from delfin.agent import credentials as cred

SECRET = "sk-test-0123456789abcdefXYZ"
OTHER = "sk-other-9876543210fedcbaQRS"


# ---------------------------------------------------------------------------
# Fixtures: a private HOME, a settings file, a keyring module with a fake
# ---------------------------------------------------------------------------


def _keyring_module(monkeypatch):
    """The real ``keyring`` package where installed, else a stand-in."""
    try:
        import keyring  # noqa: F401
        import keyring.backend  # noqa: F401
        return sys.modules["keyring"]
    except ImportError:
        pass
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


@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    (h / ".delfin").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(h))
    store = h / ".delfin" / "credentials.json"
    monkeypatch.setattr(cred, "_DEFAULT_PATH", store)
    monkeypatch.setattr(cred, "HEADLESS_TIMEOUT_S", 0.3)
    monkeypatch.setattr(cred, "_interactive", lambda: False)
    cred.reset_keyring_state()
    yield h
    cred.reset_keyring_state()


def _set_mode(home, mode):
    (home / ".delfin_settings.json").write_text(
        json.dumps({"agent": {"credentials": {"store": mode}}}), encoding="utf-8")
    cred.reset_keyring_state()


@pytest.fixture
def keyring_with(monkeypatch, home):
    """Install a fake backend behind the real factory. Returns a function
    taking the backend class; the instance is returned."""
    kr = _keyring_module(monkeypatch)
    releases: list[threading.Event] = []

    def install(cls):
        inst = cls()
        monkeypatch.setattr(kr, "get_keyring", lambda: inst, raising=False)
        monkeypatch.setattr(cred, "_keyring_factory", cred._default_keyring_factory)
        cred.reset_keyring_state()
        return inst

    install.kr = kr
    install.releases = releases
    yield install
    for ev in releases:
        ev.set()


def _fake_class(kr, *, priority=5, name="FakeSecretService"):
    class Fake(kr.backend.KeyringBackend):
        def __init__(self):
            self.items: dict = {}
            self.calls: list = []

        def get_password(self, service, username):
            self.calls.append(("get", service, username))
            return self.items.get((service, username))

        def set_password(self, service, username, password):
            self.calls.append(("set", service, username))
            self.items[(service, username)] = password

        def delete_password(self, service, username):
            self.calls.append(("delete", service, username))
            if (service, username) not in self.items:
                raise kr.errors.PasswordDeleteError(username)
            del self.items[(service, username)]

    Fake.priority = priority
    Fake.__name__ = Fake.__qualname__ = name
    return Fake


def _secret_on_disk(root, secret=SECRET) -> list:
    hits = []
    for p in root.rglob("*"):
        if p.is_file():
            try:
                if secret in p.read_text(encoding="utf-8", errors="replace"):
                    hits.append(p)
            except OSError:
                pass
    return hits


# ---------------------------------------------------------------------------
# The setting
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw,expected", [
    ("auto", "auto"), ("keyring", "keyring"), ("file", "file"),
    ("FILE", "file"), ("wat", "auto"), ("", "auto"),
])
def test_the_setting_is_read_from_the_user_settings_file(home, raw, expected):
    _set_mode(home, raw)
    assert cred.configured_store() == expected


def test_no_setting_means_auto(home):
    assert cred.configured_store() == "auto"


def test_a_workspace_settings_file_does_not_choose_the_store(home, tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "settings.json").write_text(
        json.dumps({"agent": {"credentials": {"store": "file"}}}), encoding="utf-8")
    monkeypatch.chdir(ws)
    _set_mode(home, "keyring")
    assert cred.configured_store() == "keyring"


# ---------------------------------------------------------------------------
# keyring absent: the file, as before
# ---------------------------------------------------------------------------


def test_without_the_package_auto_writes_the_file(home, monkeypatch):
    monkeypatch.setitem(sys.modules, "keyring", None)
    monkeypatch.setattr(cred, "_keyring_factory", cred._default_keyring_factory)
    cred.reset_keyring_state()
    status = cred.keyring_status()
    assert status["usable"] is False
    assert "not installed" in status["reason"]
    assert cred.set_credential("KIT_TOOLBOX_API_KEY", SECRET) is True
    stored = json.loads((home / ".delfin" / "credentials.json").read_text())
    assert stored == {"KIT_TOOLBOX_API_KEY": SECRET}
    assert cred.load_credential("KIT_TOOLBOX_API_KEY") == SECRET
    assert cred.stored_in("KIT_TOOLBOX_API_KEY") == ["file"]


def test_without_the_package_an_explicit_keyring_setting_refuses_to_write(home, monkeypatch):
    monkeypatch.setitem(sys.modules, "keyring", None)
    monkeypatch.setattr(cred, "_keyring_factory", cred._default_keyring_factory)
    _set_mode(home, "keyring")
    with pytest.raises(cred.CredentialStoreUnavailable):
        cred.set_credential("KIT_TOOLBOX_API_KEY", SECRET)
    assert _secret_on_disk(home) == []


# ---------------------------------------------------------------------------
# A secure backend: store / load / delete
# ---------------------------------------------------------------------------


def test_auto_stores_in_the_keyring_and_writes_no_secret_to_disk(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    assert cred.keyring_status()["usable"] is True
    assert cred.set_credential("KIT_TOOLBOX_API_KEY", SECRET) is True
    assert fake.items[(cred.KEYRING_SERVICE, "KIT_TOOLBOX_API_KEY")] == SECRET
    assert _secret_on_disk(home) == []
    assert cred.stored_in("KIT_TOOLBOX_API_KEY") == ["keyring"]
    # A fresh process state reads it back through the index, not a cache.
    cred.reset_keyring_state()
    assert cred.load_credential("KIT_TOOLBOX_API_KEY") == SECRET


def test_setting_the_same_value_twice_changes_nothing(keyring_with, home):
    keyring_with(_fake_class(keyring_with.kr))
    assert cred.set_credential("X_KEY", SECRET) is True
    assert cred.set_credential("X_KEY", SECRET) is False


def test_the_environment_still_wins_over_the_keyring(keyring_with, home, monkeypatch):
    keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("KIT_TOOLBOX_API_KEY", SECRET)
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", OTHER)
    assert cred.load_credential("KIT_TOOLBOX_API_KEY") == OTHER


def test_delete_removes_the_entry_from_every_store(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("X_KEY", OTHER, store="file")
    cred.set_credential("X_KEY", SECRET)
    assert sorted(cred.stored_in("X_KEY")) == ["file", "keyring"]
    assert cred.delete_credential("X_KEY") is True
    assert fake.items == {}
    assert cred.stored_in("X_KEY") == []
    assert cred.load_credential("X_KEY") == ""


def test_a_name_never_put_in_the_keyring_does_not_touch_it(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    _set_mode(home, "file")
    cred.set_credential("X_KEY", SECRET)
    fake.calls.clear()
    assert cred.load_credential("X_KEY") == SECRET
    assert fake.calls == []


def test_file_mode_keeps_the_file_even_with_a_keyring(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    _set_mode(home, "file")
    cred.set_credential("X_KEY", SECRET)
    assert fake.items == {}
    assert cred.stored_in("X_KEY") == ["file"]


def test_an_existing_file_entry_keeps_working_after_the_switch(keyring_with, home):
    keyring_with(_fake_class(keyring_with.kr))
    _set_mode(home, "file")
    cred.set_credential("OLD_KEY", SECRET)
    _set_mode(home, "auto")
    assert cred.load_credential("OLD_KEY") == SECRET


def test_the_keyring_copy_wins_over_an_older_file_copy(keyring_with, home):
    keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("X_KEY", OTHER, store="file")
    cred.set_credential("X_KEY", SECRET)
    assert cred.load_credential("X_KEY") == SECRET


# ---------------------------------------------------------------------------
# Backends that are not a secure store, or fail: the file
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("priority,name", [
    (0, "Keyring"),               # keyring.backends.fail.Keyring
    (-1, "Keyring"),              # keyring.backends.null.Keyring
    (0.5, "PlaintextKeyring"),    # keyrings.alt: a file in clear text
    (5, "PlaintextKeyring"),
])
def test_an_insecure_backend_is_not_used(keyring_with, home, priority, name):
    fake = keyring_with(_fake_class(keyring_with.kr, priority=priority, name=name))
    assert cred.keyring_status()["usable"] is False
    cred.set_credential("X_KEY", SECRET)
    assert fake.items == {}
    assert cred.stored_in("X_KEY") == ["file"]


def test_a_backend_that_fails_on_write_falls_back_to_the_file_under_auto(keyring_with, home):
    Fake = _fake_class(keyring_with.kr)

    class Locked(Fake):
        def set_password(self, service, username, password):
            raise RuntimeError("collection is locked")

    keyring_with(Locked)
    assert cred.set_credential("X_KEY", SECRET) is True
    assert cred.stored_in("X_KEY") == ["file"]
    assert cred.load_credential("X_KEY") == SECRET


def test_a_read_back_that_differs_is_not_reported_as_stored(keyring_with, home):
    Fake = _fake_class(keyring_with.kr)

    class Mangles(Fake):
        def get_password(self, service, username):
            v = super().get_password(service, username)
            return v[:-1] if v else v

    keyring_with(Mangles)
    _set_mode(home, "keyring")
    with pytest.raises(cred.CredentialStoreUnavailable):
        cred.set_credential("X_KEY", SECRET)
    assert cred.stored_in("X_KEY") == []


# ---------------------------------------------------------------------------
# Headless: a backend that blocks (an unlock prompt nobody answers)
# ---------------------------------------------------------------------------


def _blocking_class(keyring_with):
    Fake = _fake_class(keyring_with.kr)
    gate = threading.Event()
    keyring_with.releases.append(gate)

    class Blocks(Fake):
        def get_password(self, service, username):
            gate.wait(30)
            return super().get_password(service, username)

        def set_password(self, service, username, password):
            gate.wait(30)
            return super().set_password(service, username, password)

    return Blocks


def test_a_blocking_backend_times_out_to_the_file_in_headless_mode(keyring_with, home):
    fake = keyring_with(_blocking_class(keyring_with))
    t0 = time.monotonic()
    assert cred.set_credential("X_KEY", SECRET) is True
    first = time.monotonic() - t0
    assert first < 5.0
    assert cred.stored_in("X_KEY") == ["file"]
    assert fake.items == {}
    status = cred.keyring_status()
    assert status["usable"] is False and "timed out" in status["reason"]
    # The verdict is kept: the next call does not wait again.
    t1 = time.monotonic()
    assert cred.load_credential("X_KEY") == SECRET
    cred.set_credential("Y_KEY", OTHER)
    assert time.monotonic() - t1 < 0.25


def test_an_interactive_session_waits_longer_for_an_unlock(home, monkeypatch):
    monkeypatch.setattr(cred, "_interactive", lambda: True)
    assert cred._keyring_timeout() == cred.INTERACTIVE_TIMEOUT_S
    monkeypatch.setattr(cred, "_interactive", lambda: False)
    assert cred._keyring_timeout() == cred.HEADLESS_TIMEOUT_S
    assert cred.HEADLESS_TIMEOUT_S < cred.INTERACTIVE_TIMEOUT_S


def test_a_terminal_without_a_display_is_headless_on_linux(monkeypatch):
    real = cred._interactive
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)
    assert real() is False


# ---------------------------------------------------------------------------
# Migration
# ---------------------------------------------------------------------------


def test_migrate_moves_file_entries_after_a_verified_read_back(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("A_KEY", SECRET, store="file")
    cred.set_credential("B_KEY", OTHER, store="file")
    rows = {r["name"]: r for r in cred.migrate_to_keyring()}
    assert rows["A_KEY"]["action"] == "migrated"
    assert rows["B_KEY"]["action"] == "migrated"
    assert fake.items[(cred.KEYRING_SERVICE, "A_KEY")] == SECRET
    assert json.loads((home / ".delfin" / "credentials.json").read_text()) == {}
    assert _secret_on_disk(home) == [] and _secret_on_disk(home, OTHER) == []
    assert cred.load_credential("A_KEY") == SECRET


def test_migrate_keeps_the_file_entry_when_the_read_back_differs(keyring_with, home):
    Fake = _fake_class(keyring_with.kr)

    class Mangles(Fake):
        def get_password(self, service, username):
            v = super().get_password(service, username)
            return (v + "!") if v else v

    keyring_with(Mangles)
    cred.set_credential("A_KEY", SECRET, store="file")
    rows = cred.migrate_to_keyring()
    assert rows[0]["action"] == "failed"
    assert json.loads((home / ".delfin" / "credentials.json").read_text()) == {"A_KEY": SECRET}
    assert cred.load_credential("A_KEY") == SECRET


def test_migrate_does_not_overwrite_a_different_keyring_value(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("A_KEY", OTHER)            # into the keyring
    cred.set_credential("A_KEY", SECRET, store="file")
    rows = cred.migrate_to_keyring()
    assert rows[0]["action"] == "differs"
    assert fake.items[(cred.KEYRING_SERVICE, "A_KEY")] == OTHER
    assert json.loads((home / ".delfin" / "credentials.json").read_text()) == {"A_KEY": SECRET}


def test_migrate_without_a_usable_keyring_changes_nothing(home, monkeypatch):
    monkeypatch.setitem(sys.modules, "keyring", None)
    monkeypatch.setattr(cred, "_keyring_factory", cred._default_keyring_factory)
    cred.reset_keyring_state()
    cred.set_credential("A_KEY", SECRET)
    with pytest.raises(cred.CredentialStoreUnavailable):
        cred.migrate_to_keyring()
    assert json.loads((home / ".delfin" / "credentials.json").read_text()) == {"A_KEY": SECRET}


def test_migrate_refuses_when_the_setting_says_file(keyring_with, home):
    fake = keyring_with(_fake_class(keyring_with.kr))
    _set_mode(home, "file")
    cred.set_credential("A_KEY", SECRET)
    with pytest.raises(cred.CredentialStoreUnavailable):
        cred.migrate_to_keyring()
    assert fake.items == {}


def test_nothing_migrates_without_being_asked(keyring_with, home):
    keyring_with(_fake_class(keyring_with.kr))
    _set_mode(home, "file")
    cred.set_credential("A_KEY", SECRET)
    _set_mode(home, "auto")
    cred.load_credential("A_KEY")
    cred.list_credentials()
    cred.set_credential("B_KEY", OTHER)
    assert json.loads((home / ".delfin" / "credentials.json").read_text()) == {"A_KEY": SECRET}


# ---------------------------------------------------------------------------
# list: masked, with the store
# ---------------------------------------------------------------------------


def test_list_names_the_store_and_masks(keyring_with, home):
    keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("A_KEY", SECRET)
    cred.set_credential("B_KEY", OTHER, store="file")
    items = cred.list_credentials()
    assert items["A_KEY"]["source"] == "keyring"
    assert items["B_KEY"]["source"] == "file"
    assert items["A_KEY"]["value"] == cred.mask(SECRET)
    assert SECRET not in json.dumps(items) and OTHER not in json.dumps(items)


def test_cli_list_and_migrate(keyring_with, home, capsys):
    keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("A_KEY", SECRET)
    cred.set_credential("B_KEY", OTHER, store="file")
    from delfin.agent import cli as agent_cli
    assert agent_cli.main(["credentials", "list"]) == 0
    out = capsys.readouterr().out
    assert "[keyring]" in out and "[file]" in out
    assert SECRET not in out and OTHER not in out
    assert agent_cli.main(["credentials", "migrate"]) == 0
    out = capsys.readouterr().out
    assert "B_KEY" in out and OTHER not in out
    assert cred.stored_in("B_KEY") == ["keyring"]


def test_cli_migrate_says_why_when_no_keyring(home, monkeypatch, capsys):
    monkeypatch.setitem(sys.modules, "keyring", None)
    monkeypatch.setattr(cred, "_keyring_factory", cred._default_keyring_factory)
    cred.reset_keyring_state()
    cred.set_credential("A_KEY", SECRET)
    from delfin.agent import cli as agent_cli
    assert agent_cli.main(["credentials", "migrate"]) == 1
    assert "not installed" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Neither store reaches the agent's tools
# ---------------------------------------------------------------------------


def test_a_keyring_value_is_a_known_secret_for_scrub_and_egress(keyring_with, home):
    keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("A_KEY", SECRET)
    cred.reset_keyring_state()          # a later process: nothing cached
    assert SECRET in cred.held_secret_values()
    from delfin.agent.output_guard import scrub_secrets
    from delfin.agent.web_tools import outbound_secret_reason
    assert SECRET not in scrub_secrets(f"value={SECRET} end")
    assert outbound_secret_reason(f"https://example.org/?q={SECRET}")


def test_no_tool_route_returns_a_stored_secret(keyring_with, home, tmp_path):
    """A copy of the key in the workspace, read by every reading tool,
    reaches the model's context redacted; the keyring CLIs are refused."""
    from delfin.agent.api_client import (
        KitToolPermissions, _doc_executor, _redact_tool_result, _DOC_TOOLS_OPENAI)
    keyring_with(_fake_class(keyring_with.kr))
    cred.set_credential("A_KEY", SECRET)
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "leak.txt").write_text(f"token: {SECRET}\n", encoding="utf-8")
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions",
                               confirm_callback=None)
    for name, args in (
        ("read_file", {"path": "leak.txt"}),
        ("grep_file", {"path": "leak.txt", "pattern": "token"}),
        ("bash", {"command": "cat leak.txt", "description": "read"}),
    ):
        raw = _doc_executor.execute(name, args, perms)
        assert SECRET not in _redact_tool_result(raw), name
    for cmd in (
        "secret-tool lookup service delfin username A_KEY",
        "security find-generic-password -s delfin -a A_KEY -w",
        "keyring get delfin A_KEY",
        "python3 -m keyring get delfin A_KEY",
        "python3 -c \"import keyring; print(keyring.get_password('delfin','A_KEY'))\"",
    ):
        out = json.loads(_doc_executor.execute(
            "bash", {"command": cmd, "description": "x"}, perms))
        assert "error" in out and "deny" in out["error"].lower(), cmd
    # The stores themselves are deny paths.
    for p in (cred.credentials_path(), cred.index_path()):
        assert perms.matches_path_deny(str(p)), p
    # No tool reads a credential: none is named for it.
    names = [t["function"]["name"] for t in _DOC_TOOLS_OPENAI]
    assert not [n for n in names
                if any(w in n.lower() for w in ("credential", "keyring", "secret"))]


def test_a_child_process_reaches_no_keyring(tmp_path):
    from conftest import child_env
    env = child_env(tmp_path)
    assert env.get("PYTHON_KEYRING_BACKEND") == "keyring.backends.null.Keyring"
    code = ("from delfin.agent import credentials as c; "
            "print(c.keyring_status()['usable'])")
    out = subprocess.run([sys.executable, "-c", code], env=env,
                         capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().endswith("False")


def test_the_suite_default_reaches_no_keyring():
    cred.reset_keyring_state()
    assert cred.keyring_status()["usable"] is False
