"""The account key and the pins beside it are out of every agent route.

Input: a home under ``tmp_path`` with ``~/.delfin/trust_store``; the key
created by ``trust_store.derived_key``. Output: the key file's mode, the
derived keys, and what the real executor (``_doc_executor.execute``)
returns for every read route an agent has, with the shell isolation off
and on auto, with and without ``~/.delfin`` granted.

Semantics:

* the key is 32 bytes from ``secrets``, 0600 in a 0700 directory, created
  once and the same on every later read;
* a key per purpose, derived from the one account key;
* a process with the agent-session marker gets no key and creates none;
* a file in the key's place that is a link, group- or world-readable or
  of the wrong size is refused and recorded (``trust_store_refused``);
* read_file, grep_file, list_files and the shell's ``cat`` do not return
  the key's bytes; every cage masks the directory (bubblewrap tmpfs after
  the binds, Landlock hide outside a write root).
"""
from __future__ import annotations

import os
import stat

import pytest

from delfin.agent import api_client as A
from delfin.agent import security_events as SE

MARKER = "DELFIN_AGENT_SESSION"


@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.delenv(MARKER, raising=False)
    from delfin.agent import trust_store
    monkeypatch.setattr(trust_store, "store_dir",
                        lambda: h / ".delfin" / "trust_store")
    monkeypatch.setattr(trust_store, "_REFUSAL_NOTED", False)
    trust_store.forget_cached_keys()
    SE.clear()
    yield h
    trust_store.forget_cached_keys()
    SE.clear()


def test_the_key_is_private_and_stable(home):
    from delfin.agent import trust_store as T
    k1 = T.account_key()
    path = home / ".delfin" / "trust_store" / "account.key"
    assert k1 is not None and len(k1) == 32
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert stat.S_IMODE(path.parent.stat().st_mode) == 0o700
    T.forget_cached_keys()
    assert T.account_key() == k1


def test_each_purpose_has_its_own_key(home):
    from delfin.agent import trust_store as T
    a = T.derived_key("session_messages")
    b = T.derived_key("audit_log")
    assert a and b and a != b and a != T.account_key()
    tag = T.mac("audit_log", b"x")
    assert T.check_mac("audit_log", b"x", tag) is True
    assert T.check_mac("audit_log", b"y", tag) is False
    assert T.check_mac("session_messages", b"x", tag) is False


def test_an_agent_session_gets_no_key_and_makes_none(home, monkeypatch):
    from delfin.agent import trust_store as T
    monkeypatch.setenv(MARKER, "123")
    assert T.account_key() is None
    assert T.derived_key("session_messages") is None
    assert T.mac("session_messages", b"x") == ""
    assert not (home / ".delfin" / "trust_store" / "account.key").exists()
    assert T.write_record("mcp_pins.json", {"a": 1}) is False


@pytest.mark.parametrize("shape", ["world-readable", "symlink", "short"])
def test_a_key_file_that_is_not_private_is_refused(home, tmp_path, shape):
    from delfin.agent import trust_store as T
    d = home / ".delfin" / "trust_store"
    d.mkdir(parents=True, mode=0o700)
    path = d / "account.key"
    if shape == "symlink":
        target = tmp_path / "elsewhere.key"
        target.write_bytes(b"k" * 32)
        os.chmod(target, 0o600)
        os.symlink(target, path)
    else:
        path.write_bytes(b"k" * (32 if shape == "world-readable" else 5))
        os.chmod(path, 0o644 if shape == "world-readable" else 0o600)
    assert T.account_key() is None
    assert T.derived_key("audit_log") is None
    assert "trust_store_refused" in [e.kind for e in SE.recent(50)]


def test_a_record_is_written_whole_and_owner_only(home):
    from delfin.agent import trust_store as T
    assert T.write_record("mcp_pins.json", {"servers": {}}) is True
    path = home / ".delfin" / "trust_store" / "mcp_pins.json"
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert T.read_record("mcp_pins.json") == {"servers": {}}
    assert T.read_record("account.key") is None
    assert T.write_record("account.key", {}) is False


@pytest.mark.parametrize("granted", ["none", "the .delfin dir"])
@pytest.mark.parametrize("isolation", ["off", "auto"])
def test_no_agent_read_route_reaches_the_key(home, tmp_path, monkeypatch,
                                             granted, isolation):
    from delfin.agent import trust_store as T
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", isolation)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    key = T.account_key()
    path = home / ".delfin" / "trust_store" / "account.key"
    # The key as text a reader would print: its hex and its raw bytes.
    marks = (key.hex(), key.decode("latin-1"))
    ws = tmp_path / "ws"
    ws.mkdir()
    extra = (home / ".delfin",) if granted != "none" else ()
    perms = A.KitToolPermissions(workspace=ws.resolve(),
                                 mode="bypassPermissions",
                                 extra_workspace_dirs=extra)
    routes = {
        "read_file": ("read_file", {"path": str(path)}),
        "grep_file": ("grep_file", {"path": str(path), "pattern": "."}),
        "list_files": ("list_files", {"path": str(path.parent)}),
        "bash cat": ("bash", {"command": f"cat {path}"}),
        "bash xxd via ~": ("bash", {"command":
                                    "od -An -tx1 ~/.delfin/trust_store/"
                                    "account.key | tr -d ' \\n'"}),
    }
    leaked = []
    for label, (tool, args) in routes.items():
        out = str(A._doc_executor.execute(tool, args, perms))
        if any(m and m in out for m in marks):
            leaked.append(label)
    assert not leaked, (granted, isolation, leaked)


def test_the_cli_backend_denies_the_store():
    rules = A._cli_secret_deny_rules()
    assert "Read(~/.delfin/trust_store/**)" in rules
    assert "Edit(~/.delfin/trust_store/**)" in rules


def test_every_cage_masks_the_store(home, tmp_path, monkeypatch):
    """Built as on a host with bubblewrap, on any host: the directory is a
    tmpfs AFTER the binds, and Landlock hides it outside a write root."""
    from delfin.agent import boundaries, sandbox
    assert ".delfin/trust_store" in sandbox._HOME_SECRET_DIRS
    assert boundaries.MASKED_HOME_ENTRIES == tuple(sandbox._HOME_SECRET_DIRS)
    monkeypatch.setattr(A.shutil, "which",
                        lambda name, *a, **k: f"/usr/bin/{name}")
    ws = tmp_path / "ws"
    ws.mkdir()
    (home / ".delfin").mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(),
                                 mode="bypassPermissions",
                                 extra_workspace_dirs=(home / ".delfin",))
    store = str((home / ".delfin" / "trust_store").resolve())
    argv = A._bash_isolation_argv("true", ws, perms, mode="bwrap")
    pairs = list(zip(argv, argv[1:]))
    masks = [i for i, (a, b) in enumerate(pairs)
             if a == "--tmpfs" and b == store]
    binds = [i for i, (a, b) in enumerate(pairs)
             if a in ("--bind", "--bind-try")
             and b == str((home / ".delfin").resolve())]
    assert masks and binds and max(masks) > max(binds), argv
    ws_only = A.KitToolPermissions(workspace=ws.resolve(),
                                   mode="bypassPermissions")
    ll = A._landlock_argv(["true"], ws_only)
    assert store in [b for a, b in zip(ll, ll[1:]) if a == "--hide"], ll


def test_the_refusal_kind_has_a_label():
    assert "trust_store_refused" in SE.known_kinds()


def test_the_store_is_a_redirected_user_state_resolver():
    from delfin.agent import state_paths
    assert ("delfin.agent.trust_store", "store_dir", "trust_store") in \
        state_paths.USER_STATE_RESOLVERS


def test_the_suite_keeps_the_store_out_of_the_tests_workspace(tmp_path):
    """The suite's store lies beside tmp_path, as ~/.delfin lies beside a
    workspace: inside a write root its Landlock hide would refuse every
    new file at that root's top (tests/conftest.py)."""
    from delfin.agent import trust_store as T
    store = T.store_dir().resolve()
    assert tmp_path.resolve() not in store.parents, store
