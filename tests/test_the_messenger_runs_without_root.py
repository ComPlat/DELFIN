"""The messenger without root: ``delfin-agent messenger setup --user``.

For hosts without root (cluster login nodes). The same gateway, bridge and
pinned simplex-chat daemon, run as the user (delfin/messenger/rootless.py):

* systemd --user units where available, else supervised processes with pid
  files (``stop --user`` ends them, by pid AND start time);
* low priority (nice 10, idle I/O class), no whisper, an idle exit;
* the daemon's API on 127.0.0.1 at a random port. Without nftables any
  local user can connect to it: setup says so and runs only with
  ``--accept-loopback-risk``; it refuses an unpinned simplex-chat build and
  an agent session;
* the gateway and the agent service are one user now, so the bridge admits
  its own uid only for the gateway's recorded process (pid and start time
  from the pid file), never for a process carrying the agent marker.

Every test uses fakes for systemctl, the spawner and the daemon binary; no
daemon, no network.
"""
from __future__ import annotations

import hashlib
import os
import socket
import subprocess
import sys
import time

import pytest

from conftest import child_env
from delfin.messenger import rootless


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "data"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    from delfin.agent.approval_guard import ENV_MARKER
    monkeypatch.delenv(ENV_MARKER, raising=False)
    binary = tmp_path / "simplex-chat"
    binary.write_bytes(b"fake simplex-chat build")
    monkeypatch.setattr(rootless, "PINNED", rootless.PINNED | {hashlib.sha256(b"fake simplex-chat build").hexdigest()})
    return tmp_path, binary


def test_the_plan_runs_low_and_local_without_whisper(home):
    _, binary = home
    p = rootless.plan(simplex_bin=str(binary), port=34567, bot_name="DELFIN-login", idle_exit_s=3600)
    daemon = p["processes"]["daemon"]
    assert "--chat-server-port" in daemon and "34567" in daemon
    assert p["env"]["DELFIN_MESSENGER_PORT"] == "34567"
    assert p["env"]["DELFIN_MESSENGER_ROOTLESS"] == "1"
    assert not any("WHISPER" in k for k in p["env"])
    assert "--idle-exit" in p["processes"]["gateway"]
    for text in p["units"].values():
        assert "Nice=10" in text and "IOSchedulingClass=idle" in text
    assert str(rootless.state_dir()) in " ".join(daemon)


def test_the_port_is_random_and_on_loopback():
    ports = {rootless.free_port() for _ in range(5)}
    assert len(ports) > 1 and all(1024 < p < 65536 for p in ports)


def _setup(home, *extra, systemd=False):
    _, binary = home
    calls, spawned, out = [], [], []

    def run(argv, **kw):
        calls.append(list(argv))
        ok = systemd or argv[:3] != ["systemctl", "--user", "is-system-running"]
        return subprocess.CompletedProcess(argv, 0 if ok else 1, "running\n" if ok else "", "")

    def spawn(name, argv, env, log):
        spawned.append((name, list(argv), dict(env)))
        return 100000 + len(spawned)

    code = rootless.setup(["--simplex-chat", str(binary), "--no-pair", *extra],
                          out=out.append, run=run, spawn=spawn)
    return code, calls, spawned, "\n".join(out)


def test_setup_says_the_loopback_risk_and_waits_for_consent(home):
    code, calls, spawned, out = _setup(home)
    assert code == 2 and not spawned
    assert "other users of this host" in out and "--accept-loopback-risk" in out
    assert "login node" in out


def test_setup_refuses_an_agent_session(home, monkeypatch):
    from delfin.agent.approval_guard import ENV_MARKER
    monkeypatch.setenv(ENV_MARKER, "1")
    code, _, spawned, out = _setup(home, "--accept-loopback-risk")
    assert code == 2 and not spawned and "agent session" in out


def test_setup_refuses_an_unpinned_build(home, monkeypatch):
    monkeypatch.setattr(rootless, "PINNED", frozenset({"0" * 64}))
    code, _, spawned, out = _setup(home, "--accept-loopback-risk")
    assert code == 2 and not spawned and "pinned" in out


def test_with_systemd_user_it_writes_and_starts_user_units(home):
    tmp, _ = home
    code, calls, spawned, out = _setup(home, "--accept-loopback-risk", systemd=True)
    assert code == 0, out
    units = sorted(p.name for p in (tmp / "config" / "systemd" / "user").glob("*.service"))
    assert units == sorted(rootless.UNITS)
    assert ["systemctl", "--user", "daemon-reload"] in calls
    assert any(c[:4] == ["systemctl", "--user", "enable", "--now"] for c in calls)
    assert not spawned


def test_without_systemd_it_supervises_three_low_priority_processes(home):
    code, calls, spawned, out = _setup(home, "--accept-loopback-risk")
    assert code == 0, out
    assert [name for name, _, _ in spawned] == ["daemon", "bridge", "gateway"]
    for _, argv, env in spawned:
        assert argv[:3] == ["nice", "-n", "10"]
        assert env["DELFIN_MESSENGER_ROOTLESS"] == "1"


def test_stop_ends_only_the_recorded_processes(home, monkeypatch):
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"], env=child_env(home[0]))
    try:
        run_dir = rootless.state_dir() / "run"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "gateway.pid").write_text(f"{proc.pid} {rootless.start_ticks(proc.pid)}\n")
        (run_dir / "bridge.pid").write_text(f"{proc.pid} 1\n")          # a reused pid: other start
        killed = []
        assert rootless.stop(out=lambda *_: None, kill=lambda pid, sig: killed.append(pid)) == 0
        assert killed == [proc.pid]
    finally:
        proc.kill()
        proc.wait(5)


_CLIENT = r"""
import socket, sys, time
s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
s.connect(sys.argv[1])
print("connected", flush=True)
time.sleep(20)
"""


def _peer(tmp_path, marker, record):
    path = str(tmp_path / "b.sock")
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    server.bind(path)
    server.listen(1)
    env = child_env(tmp_path)
    from delfin.agent.approval_guard import ENV_MARKER
    if marker:
        env[ENV_MARKER] = "1"
    proc = subprocess.Popen([sys.executable, "-c", _CLIENT, path], env=env, stdout=subprocess.PIPE, text=True)
    conn, _ = server.accept()
    assert proc.stdout.readline().strip() == "connected"
    pidfile = tmp_path / "gateway.pid"
    ticks = rootless.start_ticks(proc.pid)
    pidfile.write_text(f"{proc.pid if record != 'other' else proc.pid + 1} {ticks}\n")
    return proc, conn, server, pidfile


@pytest.mark.parametrize("marker,record,ok", [(False, "gateway", True), (True, "gateway", False),
                                             (False, "other", False)])
def test_the_bridge_admits_its_own_uid_only_for_the_gateway(tmp_path, monkeypatch, marker, record, ok):
    from delfin.agent.service import bridge
    # Where the peer's pid and start time cannot be read (no SO_PEERCRED and
    # /proc), nobody is admitted: fail closed.
    ok = ok and sys.platform.startswith("linux")
    proc, conn, server, pidfile = _peer(tmp_path, marker, record)
    try:
        monkeypatch.setenv("DELFIN_MESSENGER_ROOTLESS", "1")
        monkeypatch.setenv("DELFIN_MESSENGER_GATEWAY_PIDFILE", str(pidfile))
        assert bridge.rootless_peer_ok(conn) is ok
        monkeypatch.delenv("DELFIN_MESSENGER_ROOTLESS")
        assert bridge.rootless_peer_ok(conn) is False, "only in rootless mode"
    finally:
        proc.kill()
        proc.wait(5)
        conn.close()
        server.close()


def test_the_pins_are_the_installers_pins():
    script = (rootless.Path(rootless.__file__).resolve().parents[2] / "scripts" / "pi" / "pi_messenger.sh").read_text()
    import re
    assert set(re.findall(r'SIMPLEX_SHA256_[0-9_a-z]+="([0-9a-f]{64})"', script)) == set(rootless.PINNED)
    assert f'SIMPLEX_VERSION="{rootless.SIMPLEX_VERSION}"' in script


def test_the_gateway_exits_when_idle(tmp_path, monkeypatch):
    import asyncio
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path / "gw"))
    from delfin.messenger.gateway import Rooms

    class Quiet:
        async def subscribe(self):
            return {}

        async def events(self):
            while True:
                await asyncio.sleep(3600)
                yield {}

    class Bridge:
        def ask(self, command):
            return {"sessions": []}

    started = time.monotonic()
    asyncio.run(asyncio.wait_for(Rooms(Quiet(), Bridge()).run(idle_exit_s=0.3, idle_poll_s=0.05), 10))
    assert time.monotonic() - started < 5
