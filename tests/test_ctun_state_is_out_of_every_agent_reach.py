"""cluster-tunnel's config, sessions and control sockets are out of every agent route.

Input: a home under ``tmp_path`` with ctun's config
(``~/.config/cluster-tunnel/config.yaml``), a session file and a listening
control socket under ``~/.cache/cluster-tunnel``. Output: what the real
executor (``_doc_executor.execute``) returns for each read and write route an
agent has, with the shell isolation off and on auto, with and without
``~/.config`` granted; whether a caged command ever connects to the socket;
and the lists every cage masks from.

Why: the control socket of a logged-in master is the user's whole identity on
a cluster -- a command that reaches it acts there with no login at all. The
config holds the ssh options and budget scripts DELFIN's own process runs
with, so a line an agent wrote into it (a ProxyCommand) would run outside
every cage; the session file holds the budget limit. delfin.remote alone
drives ctun, from DELFIN's own process.

Semantics:

* read_file, grep_file, list_files and the shell's ``cat`` return none of the
  config's or the session's bytes, and write_file and the shell change
  neither file;
* a command that finds the socket on its own (walking the home, never
  spelling the path) does not connect to it;
* every cage masks ctun's directories at their real place (bubblewrap tmpfs
  after the binds, Landlock hide, the process cage's doors), the CLI backend
  denies them, and a checkout of the cluster-tunnel repository is not
  mistaken for them.
"""

from __future__ import annotations

import shutil
import socket
import tempfile
import threading
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import sandbox

MARK = "ctun-state-marker-6c1f0e"
PROBE = """\
import os, socket
for root, dirs, files in os.walk(os.path.expanduser("~")):
    for name in files + dirs:
        path = os.path.join(root, name)
        try:
            s = socket.socket(socket.AF_UNIX)
            s.settimeout(2)
            s.connect(path)
            print("REACHED", path)
        except OSError:
            pass
"""


@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(h / ".config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(h / ".cache"))
    monkeypatch.delenv("CTUN_CONFIG", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    config = h / ".config" / "cluster-tunnel" / "config.yaml"
    config.parent.mkdir(parents=True)
    config.write_text(f"# {MARK}\nclusters:\n  kcist:\n    host: login.example\n",
                      encoding="utf-8")
    sessions = h / ".cache" / "cluster-tunnel" / "sessions"
    sessions.mkdir(parents=True)
    (sessions / "kcist.json").write_text(f'{{"limit": 100, "note": "{MARK}"}}',
                                         encoding="utf-8")
    return h


def _perms(tmp_path, home, granted):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    extra = (home / ".config",) if granted != "none" else ()
    return A.KitToolPermissions(workspace=ws.resolve(), mode="bypassPermissions",
                                extra_workspace_dirs=extra), ws


@pytest.mark.parametrize("granted", ["none", "the .config dir"])
@pytest.mark.parametrize("isolation", ["off", "auto"])
def test_no_agent_read_route_reaches_ctuns_files(home, tmp_path, monkeypatch,
                                                  granted, isolation):
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", isolation)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    perms, _ = _perms(tmp_path, home, granted)
    config = home / ".config" / "cluster-tunnel" / "config.yaml"
    session = home / ".cache" / "cluster-tunnel" / "sessions" / "kcist.json"
    routes = {
        "read_file config": ("read_file", {"path": str(config)}),
        "read_file session": ("read_file", {"path": str(session)}),
        "grep_file": ("grep_file", {"path": str(config), "pattern": "."}),
        "list_files": ("list_files", {"path": str(session.parent)}),
        "bash cat": ("bash", {"command": f"cat {config}"}),
        "bash cat via ~": ("bash", {"command": "cat ~/.config/cluster-tunnel/config.yaml"}),
        "bash grep -r ~": ("bash", {"command": f"grep -r {MARK[:12]} ~/.cache"}),
    }
    leaked = [label for label, (tool, args) in routes.items()
              if MARK in str(A._doc_executor.execute(tool, args, perms))]
    assert not leaked, (granted, isolation, leaked)


@pytest.mark.parametrize("isolation", ["off", "auto"])
def test_no_agent_write_route_changes_ctuns_config(home, tmp_path, monkeypatch, isolation):
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", isolation)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    perms, _ = _perms(tmp_path, home, "none")
    config = home / ".config" / "cluster-tunnel" / "config.yaml"
    before = config.read_bytes()
    planted = "defaults:\n  ssh_config: /tmp/evil\n"
    A._doc_executor.execute("write_file", {"path": str(config), "content": planted}, perms)
    A._doc_executor.execute(
        "bash", {"command": "printf 'x\\n' >> ~/.config/cluster-tunnel/config.yaml"}, perms)
    assert config.read_bytes() == before


@pytest.fixture()
def short_home(monkeypatch):
    """A home whose socket path fits AF_UNIX's ~100 bytes (pytest's tmp_path
    does not)."""
    root = Path(tempfile.mkdtemp(prefix="dsock"))
    h = root / "home"
    (h / ".config" / "cluster-tunnel").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(h / ".config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(h / ".cache"))
    monkeypatch.delenv("CTUN_CONFIG", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    yield h
    shutil.rmtree(root, ignore_errors=True)


@pytest.mark.parametrize("isolation", ["off", "auto"])
def test_a_command_that_finds_the_control_socket_cannot_connect(short_home, monkeypatch,
                                                                 isolation):
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", isolation)
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    home = short_home
    sockets = home / ".cache" / "cluster-tunnel" / "sockets"
    sockets.mkdir(parents=True)
    path = sockets / "kcist"
    assert len(str(path)) < 100, path
    server = socket.socket(socket.AF_UNIX)
    server.bind(str(path))
    server.listen(4)
    server.settimeout(0.2)
    reached: list[bool] = []
    stop = threading.Event()

    def accept():
        while not stop.is_set():
            try:
                conn, _ = server.accept()
                reached.append(True)
                conn.close()
            except OSError:
                continue

    t = threading.Thread(target=accept, daemon=True)
    t.start()
    try:
        perms, ws = _perms(home.parent, home, "none")
        (ws / "probe.py").write_text(PROBE, encoding="utf-8")
        out = str(A._doc_executor.execute("bash", {"command": "python3 probe.py"}, perms))
    finally:
        stop.set()
        t.join(timeout=5)
        server.close()
    assert not reached and "REACHED" not in out, (isolation, out[-400:])


def test_every_cage_masks_ctuns_directories_where_they_really_are(tmp_path, monkeypatch):
    """Redirected by XDG, as on a system that sets it: the real places count."""
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg-config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg-cache"))
    config_dir = (tmp_path / "xdg-config" / "cluster-tunnel").resolve()
    cache_dir = (tmp_path / "xdg-cache" / "cluster-tunnel").resolve()
    masked = sandbox.runtime_masked_dirs()
    assert str(config_dir) in masked and str(cache_dir) in masked
    assert config_dir.is_dir() and cache_dir.is_dir(), "created so nothing is planted first"
    assert {str(config_dir), str(cache_dir)} <= set(A._home_secret_paths())
    assert {str(config_dir), str(cache_dir)} <= set(A._session_doors())
    opts = A._process_cage_options()
    tmpfs = {opts[i + 1] for i, o in enumerate(opts) if o == "--tmpfs"}
    assert {str(config_dir), str(cache_dir)} <= tmpfs


def test_the_default_locations_are_in_the_sandboxs_home_list():
    assert ".config/cluster-tunnel" in sandbox._HOME_SECRET_DIRS
    assert ".cache/cluster-tunnel" in sandbox._HOME_SECRET_DIRS


def test_the_cli_backend_denies_ctuns_directories():
    rules = A._cli_secret_deny_rules()
    for d in ("~/.config/cluster-tunnel/**", "~/.cache/cluster-tunnel/**"):
        assert f"Read({d})" in rules and f"Edit({d})" in rules


def test_a_checkout_of_the_cluster_tunnel_repository_stays_ordinary(tmp_path):
    perms = A.KitToolPermissions(workspace=tmp_path.resolve(), mode="default")
    assert not perms.matches_path_deny("/home/u/project/cluster-tunnel/README.md")
    assert not perms.matches_path_deny("cluster-tunnel/cluster_tunnel/api.py")
    assert perms.matches_path_deny("/home/u/.config/cluster-tunnel/config.yaml")
    assert perms.matches_path_deny("/home/u/.cache/cluster-tunnel/sockets/kcist")
    assert perms.matches_path_deny("/home/u/.config/cluster-tunnel")
