"""The remote layer end to end against a real sshd: pin, log in, run, capture, transfer.

Opt-in: ``DELFIN_TEST_SSHD=<host>`` (an sshd that accepts this user) and
``DELFIN_TEST_SSHD_KEY=<private key>`` (key authentication, no password).
Input: ctun's config in a home under ``tmp_path`` (its default place, so the
registry accepts it), the host key pinned through ctun's own API, DELFIN's
``remote`` settings as a dict. Output: what ``delfin.remote.connection``
returns for a command and a round trip of files; the transfers stay under
``tmp_path`` (the "cluster" is this machine).

Measured on 2026-10-20 against the local OpenSSH sshd with an ed25519 key.
"""

from __future__ import annotations

import filecmp
import getpass
import os
import shutil
import tempfile
from pathlib import Path

import pytest

from cluster_tunnel import api as ctun_api
from delfin.remote import connection, errors

HOST = os.environ.get("DELFIN_TEST_SSHD", "")
KEY = os.environ.get("DELFIN_TEST_SSHD_KEY", "")

pytestmark = pytest.mark.skipif(
    not (HOST and KEY) or not shutil.which("ssh") or not shutil.which("rsync"),
    reason="set DELFIN_TEST_SSHD=<host> and DELFIN_TEST_SSHD_KEY=<key> to run against an sshd")


@pytest.fixture
def cluster(monkeypatch):
    # A short root: ssh's control socket path must fit AF_UNIX's ~100 bytes.
    root = Path(tempfile.mkdtemp(prefix="dsshd"))
    home = root / "h"
    ctun_dir = home / ".config" / "cluster-tunnel"
    ctun_dir.mkdir(parents=True)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(home / ".config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(home / ".cache"))
    monkeypatch.delenv("CTUN_CONFIG", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    config = ctun_dir / "config.yaml"
    config.write_text(
        f"clusters:\n  local:\n    host: {HOST}\n    user: {getpass.getuser()}\n"
        f"    identity_file: {Path(KEY).expanduser()}\n    requires_password: false\n",
        encoding="utf-8")
    settings = {"remote": {"clusters": {"local": {"work_dir": str(root / "work")}}}}
    try:
        yield config, settings, root
    finally:
        try:
            ctun_api.logout("local", config_path=config)
        finally:
            shutil.rmtree(root, ignore_errors=True)


def test_an_unpinned_sshd_is_refused_and_a_pinned_one_runs_and_transfers(cluster):
    config, settings, root = cluster
    with pytest.raises(errors.ClusterNotUsable):
        connection.connect("local", settings=settings)

    # Pin the key the server presents, through ctun as `hostkey --pin` does.
    offered = ctun_api.host_key("local", config_path=config).offered
    ctun_api.pin_host_key("local", offered.fingerprint, config_path=config)

    conn = connection.connect("local", settings=settings)
    assert conn.is_live() is False
    with pytest.raises(errors.LoginRequired):
        conn.run(["true"])
    ctun_api.login("local", config_path=config, timeout=60)
    assert conn.is_live() is True

    result = conn.run(["sh", "-c", 'printf "%s|" "$1" "$2"; exit 3', "sh", "a b", "$HOME"])
    assert (result.exit_code, result.stdout) == (3, "a b|$HOME|")

    src = root / "up"
    (src / "sub").mkdir(parents=True)
    (src / "x.txt").write_text("over the tunnel\n", encoding="utf-8")
    (src / "sub" / "y.bin").write_bytes(bytes(range(256)) * 8)
    remote_dir = root / "remote"
    remote_dir.mkdir()
    assert conn.upload(src, f"{remote_dir}/").ok
    back = root / "back"
    assert conn.download(f"{remote_dir}/up", back).ok
    got = back / "up" if (back / "up").is_dir() else back
    cmp = filecmp.dircmp(src, got)
    assert not (cmp.diff_files or cmp.left_only or cmp.right_only)
    assert filecmp.cmp(src / "sub" / "y.bin", got / "sub" / "y.bin", shallow=False)
