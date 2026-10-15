"""The dashboard approval runner shows the confinement it applies.

Input: a ``SandboxConfig`` and the host's isolation probes. Output: the
state the dashboard prints before running an approved command
(``mode=..., fs=..., net=...``) and the argv the runner starts.

Measured 2026-10-15 on a host with Landlock ABI 4 and the socket guard:
in ``allowlist`` mode with ``allow_network=False`` the dashboard printed
``net=off`` while the command connected to a TCP listener on 127.0.0.1;
the Landlock helper was started without ``--net`` (its default is
"open"), and ``DELFIN_AGENT_SANDBOX_NETWORK`` acted only on bubblewrap and
firejail. On a host without Landlock the same mode ran ``bash -c`` with
no confinement at all, also shown as ``net=off``.

The behavioural test holds on every host: whatever the runner reports,
the command's connect has to agree with it.
"""
from __future__ import annotations

import re
import socket
import threading
from pathlib import Path

import pytest

from delfin.agent import sandbox as SB


def _shown_network(cfg) -> str:
    """The network state the dashboard prints for *cfg*."""
    if hasattr(SB, "confinement"):
        return SB.confinement(cfg).network
    # The expression the dashboard used before the truth function existed.
    return "on" if cfg.allow_network else "off"


@pytest.fixture
def listener():
    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.bind(("127.0.0.1", 0))
    srv.listen(8)
    srv.settimeout(20)
    stop = threading.Event()

    def accept():
        while not stop.is_set():
            try:
                conn, _ = srv.accept()
                conn.close()
            except OSError:
                return

    t = threading.Thread(target=accept, daemon=True)
    t.start()
    yield srv.getsockname()[1]
    stop.set()
    srv.close()


def _connect_cmd(port: int) -> str:
    return ("python3 -c \"import socket; "
            f"socket.create_connection(('127.0.0.1', {port}), timeout=5); "
            "print('CONNECTED')\"")


def test_a_network_shown_as_off_is_off(tmp_path, listener, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    repo = tmp_path / "repo"
    repo.mkdir()
    cfg = SB.SandboxConfig("allowlist", False, 60)
    shown = _shown_network(cfg)
    res = SB.run_agent_command(_connect_cmd(listener), repo, config=cfg)
    assert not res.blocked, res.block_reason
    connected = "CONNECTED" in res.stdout
    if shown in ("off", "off-tcp"):
        assert not connected, (
            f"the dashboard shows net={shown} and the command reached "
            f"127.0.0.1:{listener}")
    else:
        assert shown == "open", shown


# ---------------------------------------------------------------------------
# One state per host, from the same probes the runner uses


def _host(monkeypatch, *, landlock: bool, abi: int, guard: bool,
          filt: bool = True):
    from delfin.agent import api_client as A
    from delfin.agent import landlock_exec as LL
    from delfin.agent import socket_guard as SG
    monkeypatch.setattr(A, "_landlock_functional", lambda: landlock)
    monkeypatch.setattr(LL, "abi_version", lambda: abi)
    monkeypatch.setattr(SG, "available", lambda: guard)
    monkeypatch.setattr(SG, "filter_supported", lambda: filt)


def test_allowlist_with_the_socket_guard_is_off(monkeypatch):
    _host(monkeypatch, landlock=True, abi=4, guard=True)
    c = SB.confinement(SB.SandboxConfig("allowlist", False, 60))
    assert (c.filesystem, c.network, c.landlock_net) == ("landlock", "off", "none")
    assert c.notice == ""


def test_allowlist_with_landlock_tcp_only_says_udp_is_open(monkeypatch):
    _host(monkeypatch, landlock=True, abi=4, guard=False)
    c = SB.confinement(SB.SandboxConfig("allowlist", False, 60))
    assert (c.network, c.landlock_net) == ("off-tcp", "none")
    assert "UDP" in c.notice


def test_allowlist_without_a_network_layer_is_open(monkeypatch):
    _host(monkeypatch, landlock=True, abi=3, guard=False)
    c = SB.confinement(SB.SandboxConfig("allowlist", False, 60))
    assert (c.filesystem, c.network, c.landlock_net) == ("landlock", "open", "open")
    assert "NOT restricted" in c.notice


def test_allowlist_without_landlock_is_unconfined_and_says_so(monkeypatch):
    _host(monkeypatch, landlock=False, abi=0, guard=True)
    c = SB.confinement(SB.SandboxConfig("allowlist", False, 60))
    assert (c.filesystem, c.network, c.landlock_net) == ("none", "open", "")
    assert "unconfined" in c.notice
    assert "fs=none" in c.label() and "net=open" in c.label()


def test_allowlist_with_the_network_opted_in_is_on(monkeypatch):
    _host(monkeypatch, landlock=True, abi=4, guard=True)
    c = SB.confinement(SB.SandboxConfig("allowlist", True, 60))
    assert (c.network, c.landlock_net, c.notice) == ("on", "open", "")


@pytest.mark.parametrize("mode,fs", [("bwrap", "bwrap"), ("firejail", "firejail")])
def test_a_namespace_mode_is_off(monkeypatch, mode, fs):
    _host(monkeypatch, landlock=False, abi=0, guard=False)
    c = SB.confinement(SB.SandboxConfig(mode, False, 60))
    assert (c.filesystem, c.network, c.notice) == (fs, "off", "")
    assert SB.confinement(SB.SandboxConfig(mode, True, 60)).network == "on"


def test_the_off_mode_is_open_and_says_so(monkeypatch):
    _host(monkeypatch, landlock=True, abi=4, guard=True)
    c = SB.confinement(SB.SandboxConfig("off", False, 60))
    assert (c.filesystem, c.network) == ("none", "open")
    assert c.notice


@pytest.mark.parametrize("landlock,abi,guard,filt", [
    (True, 4, True, True), (True, 4, False, True), (True, 4, False, False),
    (True, 3, False, True), (True, 2, True, True), (False, 0, False, False),
])
@pytest.mark.parametrize("net", [False, True])
def test_the_argv_applies_what_is_shown(monkeypatch, tmp_path, landlock, abi,
                                        guard, filt, net):
    _host(monkeypatch, landlock=landlock, abi=abi, guard=guard, filt=filt)
    cfg = SB.SandboxConfig("allowlist", net, 60)
    c = SB.confinement(cfg)
    argv = SB._unsandboxed_argv("true", tmp_path, "allowlist", allow_network=net)
    if c.filesystem == "none":
        assert argv == ["bash", "-c", "true"]
        return
    head = argv[:argv.index("--")]
    assert head[head.index("--net") + 1] == c.landlock_net
    # A restricted network is applied strictly: if the helper cannot hold
    # it, the command is refused instead of running with the network open.
    restricted = c.network in ("off", "off-tcp")
    assert head[head.index("--strict") + 1] == ("1" if restricted else "0")


def test_the_dashboard_prints_the_truth_function():
    src = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    start = src.index("# --- Bash commands: route through delfin.agent.sandbox")
    window = src[start:start + 3000]
    assert "confinement(" in window
    assert not re.search(r"cfg\.allow_network\s+else\s+'off'", window)
