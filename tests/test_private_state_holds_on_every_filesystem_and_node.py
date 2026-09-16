"""Private state stays private on every filesystem and across login nodes.

Audit 2026-09-16 (host dependencies):
- the credential store's temporary file was created under the umask and
  chmodded afterwards, world-readable for a moment on a default umask and
  for good on a filesystem that ignores modes;
- a session lock taken on another login node was judged by its pid, which
  names nothing on this node, and broken;
- the parent-death signal loaded glibc by name and did nothing on musl;
- MCP servers and hooks inherited the model provider keys.
"""
import json
import os
import stat

import pytest

from delfin.agent import credentials as CR
from delfin.agent import hooks as HK
from delfin.agent import mcp_client as MC
from delfin.agent import session_store as SS


def test_the_key_store_is_created_private(tmp_path, monkeypatch):
    monkeypatch.setattr(os, "umask", os.umask)      # untouched, just explicit
    old = os.umask(0o022)
    try:
        path = tmp_path / "creds" / "credentials.json"
        CR.set_credential("KIT_TOOLBOX_API_KEY", "value-1234567890", path=path)
    finally:
        os.umask(old)
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert json.loads(path.read_text())["KIT_TOOLBOX_API_KEY"] == "value-1234567890"


def test_a_filesystem_that_cannot_keep_it_private_gets_no_key(tmp_path, monkeypatch):
    real_fstat = os.fstat

    def ignores_modes(fd):
        st = real_fstat(fd)
        fields = list(st)
        fields[stat.ST_MODE] = stat.S_IFREG | 0o777
        return os.stat_result(fields)

    monkeypatch.setattr(os, "fstat", ignores_modes)
    path = tmp_path / "credentials.json"
    with pytest.raises(CR.CredentialStoreNotPrivate):
        CR.set_credential("KIT_TOOLBOX_API_KEY", "secret-value-123456", path=path)
    assert not path.exists()
    assert not path.with_suffix(".json.tmp").exists()


def test_a_symlink_planted_as_the_temporary_file_is_not_followed(tmp_path):
    victim = tmp_path / "victim.txt"
    victim.write_text("untouched\n")
    path = tmp_path / "credentials.json"
    path.with_suffix(".json.tmp").symlink_to(victim)
    CR.set_credential("OPENAI_API_KEY", "value-abcdefghijk", path=path)
    assert victim.read_text() == "untouched\n"


def test_a_lock_taken_on_another_login_node_holds(tmp_path, monkeypatch):
    monkeypatch.setattr(SS, "_ensure_dir", lambda: tmp_path)
    monkeypatch.setattr(SS, "_lock_path", lambda sid: tmp_path / f"{sid}.lock")
    import time
    (tmp_path / "s1.lock").write_text(json.dumps(
        {"pid": 999999999, "ts": time.time(), "host": "another-login-node"}))
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("s1")
    (tmp_path / "s1.lock").write_text(json.dumps(
        {"pid": os.getpid(), "ts": time.time(), "host": "another-login-node"}))
    with pytest.raises(SS.SessionLockedError):
        SS.acquire_session_lock("s1")          # same pid number, other node
    SS.release_session_lock("s1")
    assert (tmp_path / "s1.lock").exists()     # not ours to release


def test_a_dead_local_holder_still_breaks_the_lock(tmp_path, monkeypatch):
    monkeypatch.setattr(SS, "_ensure_dir", lambda: tmp_path)
    import time
    (tmp_path / "s2.lock").write_text(json.dumps(
        {"pid": 999999999, "ts": time.time(), "host": SS._this_host()}))
    SS.acquire_session_lock("s2")
    assert json.loads((tmp_path / "s2.lock").read_text())["pid"] == os.getpid()


def test_the_parent_death_signal_does_not_name_glibc():
    import inspect
    from delfin.agent import lifeline
    assert 'CDLL("libc.so.6"' not in inspect.getsource(lifeline.parent_death_signal)


def test_hooks_and_mcp_servers_do_not_get_the_provider_keys(monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "kit-provider-key")
    monkeypatch.setenv("GITHUB_TOKEN", "a-servers-own-token")
    env = HK._build_env("PreToolUse", tool_name="bash")
    assert "KIT_TOOLBOX_API_KEY" not in env and env["GITHUB_TOKEN"] == "a-servers-own-token"
    import inspect
    src = inspect.getsource(MC)
    assert "if k not in _DELFIN_PROVIDER_KEYS}" in src
