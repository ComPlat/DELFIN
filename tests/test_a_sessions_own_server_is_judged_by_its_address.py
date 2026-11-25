"""A session's own server on the loopback is judged by its ADDRESS, not only
by its port number.

Input: ``egress_proxy._own_server`` and ``EgressProxy.settle`` for a
loopback request, where the session's own background job listens on one
loopback address and something else listens on a DIFFERENT address of the
same port number.

Before: ``_own_ports`` returned the port NUMBERS the session's jobs listen
on, and ``_own_server`` answered True for every loopback address carrying
one of those numbers -- ``_loopback_addr`` was consulted only afterwards.
A job of the session listening on ``[::1]:P`` therefore opened
``127.0.0.1:P`` to whatever was there, a host service included.

Measured 2026-11-25 (``scratch/x2probe_portshare.py``): two listeners on one
port number in one network namespace is possible (``::1:P`` and
``127.0.0.1:P`` are different endpoints), and
``socket_guard.loopback_listeners`` is address-aware while the port set was
not -- so the guard refused the wrong address and the proxy admitted it.

After: the rule is ``socket_guard.owns_loopback`` at the address the proxy
actually connects to, for the inodes the session's jobs hold. The rule the
package asks for keeps its two sides: the session's own listener is
reached, a listener that is not the session's stays refused.

The process and socket tree is RECORDED (``socket_guard._PROC``,
``bash_jobs._PROC``), so every test here runs on every machine.
"""
from __future__ import annotations

import ipaddress
import sys

import pytest

from delfin.agent import bash_jobs as BJ
from delfin.agent import egress_proxy as EP
from delfin.agent import socket_guard as SG

LOOPBACK = ipaddress.ip_address("127.0.0.1")
LOOPBACK6 = ipaddress.ip_address("::1")

#: The port both listeners share: the session's job on ``::1``, a host
#: service on ``127.0.0.1``.
SHARED_PORT = 9100
OWN_INODE = "11"
FOREIGN_INODE = "22"


def _hex_address(ip: str, port: int) -> str:
    """A ``local_address`` field as the kernel of THIS machine prints it:
    every 32-bit word in host byte order."""
    raw = ipaddress.ip_address(ip).packed
    if sys.byteorder == "little":
        raw = b"".join(raw[i:i + 4][::-1] for i in range(0, len(raw), 4))
    return f"{raw.hex().upper()}:{port:04X}"


class ProcTree:
    """A recorded /proc: processes with their parent and socket fds, and
    the kernel's two TCP tables."""

    def __init__(self, root):
        self.root = root
        (root / "net").mkdir(parents=True, exist_ok=True)
        self.rows = {"tcp": [], "tcp6": []}
        self._write()

    def process(self, pid: int, ppid: int = 1, socket_inodes=()):
        d = self.root / str(pid)
        d.mkdir(exist_ok=True)
        (d / "stat").write_text(f"{pid} (python) S {ppid} " + "0 " * 20 + "\n")
        fd = d / "fd"
        fd.mkdir(exist_ok=True)
        for i, inode in enumerate(socket_inodes):
            (fd / str(i)).symlink_to(f"socket:[{inode}]")
        return self

    def listener(self, ip: str, port: int, inode: str):
        row = (f"   0: {_hex_address(ip, port)} 00000000:0000 0A "
               f"00000000:00000000 00:00000000 00000000  1000  0 {inode} 1\n")
        self.rows["tcp6" if ":" in ip else "tcp"].append(row)
        self._write()
        return self

    def _write(self):
        for table, rows in self.rows.items():
            (self.root / "net" / table).write_text(
                "  sl  local_address rem_address   st\n" + "".join(rows))


@pytest.fixture
def proc(tmp_path, monkeypatch):
    """The recorded tree, read by both modules that consult /proc."""
    tree = ProcTree(tmp_path / "proc")
    monkeypatch.setattr(SG, "_PROC", str(tree.root))
    monkeypatch.setattr(BJ, "_PROC", str(tree.root))
    return tree


class _Perms:
    task_session_id = "sess-own"
    mode = "bypassPermissions"
    confirm_callback = None


def _own_job_pids(pid: int):
    """The seam BOTH ``own_listening_sockets`` and ``own_listening_ports``
    read: with the session's job pid named, each derives its answer from the
    recorded /proc through its own code, so a test never patches the rule it
    is measuring."""
    return lambda sid: {pid} if sid == "sess-own" else set()


@pytest.fixture
def session_listening_on_ipv6(proc, monkeypatch):
    """The session's job listens on ``[::1]:PORT``; a service that is NOT
    the session's listens on ``127.0.0.1:PORT`` -- one port number, two
    endpoints, in one network namespace."""
    proc.process(600, socket_inodes=[OWN_INODE])
    proc.process(700, socket_inodes=[FOREIGN_INODE])
    proc.listener("::1", SHARED_PORT, OWN_INODE)
    proc.listener("127.0.0.1", SHARED_PORT, FOREIGN_INODE)
    monkeypatch.setattr(BJ, "_own_job_pids", _own_job_pids(600))
    return proc


def _proxy():
    return EP.EgressProxy(allowed_domains=(), allowed_ports=(80, 443))


def test_the_own_address_is_still_the_sessions_own(session_listening_on_ipv6):
    """The control on both sides: at the address it listens on, the job's
    server is the session's, before and after the rule is tightened."""
    assert SG.loopback_listeners(LOOPBACK6, SHARED_PORT) == {OWN_INODE}
    assert SG.owns_loopback(LOOPBACK6, SHARED_PORT, {OWN_INODE})
    assert EP._own_server(_Perms(), "::1", SHARED_PORT) is True


def test_the_same_port_number_on_another_address_is_somebody_elses(
        session_listening_on_ipv6):
    """The defect: the port number alone said "ours" for an address the
    session does not listen on."""
    # The address-aware rule already answers correctly ...
    assert SG.loopback_listeners(LOOPBACK, SHARED_PORT) == {FOREIGN_INODE}
    assert not SG.owns_loopback(LOOPBACK, SHARED_PORT, {OWN_INODE})
    # ... and the proxy must answer the way the guard does.
    assert EP._own_server(_Perms(), "127.0.0.1", SHARED_PORT) is False


def test_a_post_to_that_address_does_not_pass_as_the_sessions_own(
        session_listening_on_ipv6):
    """The write side of the same defect: with the port-number rule a POST
    there passed ``settle`` as the session's own server. It is refused now,
    and by the rule that owns the reason: the port gate (the session's port
    is not otherwise an allowed one) before the upload rule could speak."""
    p = _proxy()
    try:
        why = p.settle("127.0.0.1", SHARED_PORT, method="POST", body_bytes=10,
                       perms=_Perms())
        assert why, "a POST to a listener that is not the session's passed"
        assert "port 9100 is not allowed" in why, why
        # The session's own server on its own address is reached, not asked.
        assert p.settle("::1", SHARED_PORT, method="POST", body_bytes=10,
                        perms=_Perms()) is None
    finally:
        p.close()


def test_a_job_listening_on_the_address_still_passes_a_post(proc, monkeypatch):
    """The behaviour the by-port rule existed for, unchanged: the agent's
    POST to the server its own job started is not an upload."""
    proc.process(600, socket_inodes=[OWN_INODE])
    proc.listener("127.0.0.1", SHARED_PORT, OWN_INODE)
    monkeypatch.setattr(BJ, "_own_job_pids", _own_job_pids(600))
    p = _proxy()
    try:
        assert p.settle("127.0.0.1", SHARED_PORT, method="POST", body_bytes=10,
                        perms=_Perms()) is None
        # Another session's port is not this session's.
        assert p.settle("127.0.0.1", SHARED_PORT, method="POST", body_bytes=10)
    finally:
        p.close()


def test_another_session_does_not_inherit_the_listener(proc, monkeypatch):
    """The inodes are the session's own: a session that holds none gets no
    loopback server, on the address or elsewhere."""
    proc.process(600, socket_inodes=[OWN_INODE])
    proc.listener("127.0.0.1", SHARED_PORT, OWN_INODE)
    monkeypatch.setattr(BJ, "_own_job_pids", _own_job_pids(600))
    assert EP._own_server(_Perms(), "127.0.0.1", SHARED_PORT) is True
    # Another session names no job of its own: nothing is its own.
    monkeypatch.setattr(BJ, "_own_job_pids", lambda sid: set())
    assert EP._own_server(_Perms(), "127.0.0.1", SHARED_PORT) is False
