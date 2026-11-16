"""A sandboxed command reaches the servers it started itself on the
loopback, and no other server there.

Measured 2026-11-15 in an unattended swarm run: 118 refusals of the egress
proxy ("a POST request to 127.0.0.1 is an upload, and nobody can be
asked"), each a test posting to the scripted model endpoint it served to
itself; the agent could not run the test its task asked for. The loopback
of a shared host also carries other users' services and the user's own
unsandboxed ones, so it is not opened as a whole: a connect passes when
every socket listening there is the command's own -- held by its process
tree, or a listener of the session's background jobs (``--own-listener``).
The socket guard judges a direct connect that way, the proxy a request by
the pid of the helper that forwarded it (``egress_proxy._own_server``).

The rule is checked here against a RECORDED process and socket tree
(``socket_guard._PROC``, ``bash_jobs._PROC``), so every test in this file
runs on every machine. A recorded tree also states cases a live host does
not hand out on demand: a second listener on the same port, and a
listener of a process that is not readable. The tests that need the real
seccomp guard and Landlock -- a child process under the guard, reaching
its own server and no other -- live with the other guard tests in
tests/test_a_sandboxed_command_reaches_the_network_through_its_proxy.py.
"""
from __future__ import annotations

import ipaddress
import sys

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ
from delfin.agent import egress_proxy as EP
from delfin.agent import socket_guard as SG

LOOPBACK = ipaddress.ip_address("127.0.0.1")


def _hex_address(ip: str, port: int) -> str:
    """A ``local_address`` field as the kernel of THIS machine prints it:
    every 32-bit word in host byte order."""
    raw = ipaddress.ip_address(ip).packed
    if sys.byteorder == "little":
        raw = b"".join(raw[i:i + 4][::-1] for i in range(0, len(raw), 4))
    return f"{raw.hex().upper()}:{port:04X}"


def _row(ip: str, port: int, inode: str, state: str = "0A") -> str:
    return (f"   0: {_hex_address(ip, port)} 00000000:0000 {state} "
            f"00000000:00000000 00:00000000 00000000  1000  0 {inode} 1\n")


class ProcTree:
    """A recorded /proc: processes with their parent and their socket fds,
    and the kernel's two TCP tables."""

    def __init__(self, root):
        self.root = root
        (root / "net").mkdir(parents=True, exist_ok=True)
        self.rows = {"tcp": [], "tcp6": []}
        self._write()

    def process(self, pid: int, ppid: int = 1, socket_inodes=(), readable=True):
        d = self.root / str(pid)
        d.mkdir(exist_ok=True)
        (d / "stat").write_text(f"{pid} (python) S {ppid} " + "0 " * 20 + "\n")
        fd = d / "fd"
        fd.mkdir(exist_ok=True)
        if not readable:
            # A process of another user: the directory cannot be listed.
            fd.chmod(0o000)
            return self
        for i, inode in enumerate(socket_inodes):
            (fd / str(i)).symlink_to(f"socket:[{inode}]")
        return self

    def listener(self, ip: str, port: int, inode: str):
        self.rows["tcp6" if ":" in ip else "tcp"].append(_row(ip, port, inode))
        self._write()
        return self

    def connection(self, ip: str, port: int, inode: str):
        """Not a listener: state 01 = ESTABLISHED."""
        table = "tcp6" if ":" in ip else "tcp"
        self.rows[table].append(_row(ip, port, inode, state="01"))
        self._write()
        return self

    def _write(self):
        for table, rows in self.rows.items():
            (self.root / "net" / table).write_text(
                "  sl  local_address rem_address   st\n" + "".join(rows))


@pytest.fixture
def proc(tmp_path, monkeypatch):
    tree = ProcTree(tmp_path / "proc")
    monkeypatch.setattr(SG, "_PROC", str(tree.root))
    monkeypatch.setattr(BJ, "_PROC", str(tree.root))
    yield tree
    for d in tree.root.glob("*/fd"):
        d.chmod(0o755)


# -- the address the kernel printed -----------------------------------------

def test_a_proc_net_address_is_read_in_host_order():
    for ip, port in (("127.0.0.1", 8080), ("::1", 80), ("0.0.0.0", 443)):
        assert SG._proc_net_address(_hex_address(ip, port)) == (
            ipaddress.ip_address(ip), port)
    if sys.byteorder == "little":
        assert SG._proc_net_address("0100007F:1F90") == (LOOPBACK, 8080)
        assert SG._proc_net_address(
            "00000000000000000000000001000000:0050") == (
                ipaddress.ip_address("::1"), 80)


# -- the ownership rule ------------------------------------------------------

def test_only_a_listener_the_caller_holds_is_its_own(proc):
    proc.process(100, socket_inodes=["11"]).listener("127.0.0.1", 9000, "11")
    assert SG.loopback_listeners(LOOPBACK, 9000) == {"11"}
    assert SG.owns_loopback(LOOPBACK, 9000, (), 100)
    # Named as a background job's listener, without the tree.
    assert SG.owns_loopback(LOOPBACK, 9000, {"11"})
    # Below a process that holds nothing it is somebody else's, and the
    # root's own sockets do not count when the root is excluded.
    proc.process(200)
    assert not SG.owns_loopback(LOOPBACK, 9000, (), 200)
    assert not SG.owns_loopback(LOOPBACK, 9000, (), 100, include_root=False)
    # A child of the root holds it: that is the command's too.
    proc.process(101, ppid=100, socket_inodes=["11"])
    assert SG.owns_loopback(LOOPBACK, 9000, (), 100, include_root=False)


def test_nothing_listening_and_a_foreign_address_are_not_owned(proc):
    proc.process(100, socket_inodes=["11"])
    # Nothing listens on 9000: a connect there reaches nobody.
    assert not SG.owns_loopback(LOOPBACK, 9000, (), 100)
    # An established connection on the port is not a listener.
    proc.connection("127.0.0.1", 9000, "11")
    assert SG.loopback_listeners(LOOPBACK, 9000) == set()
    assert not SG.owns_loopback(LOOPBACK, 9000, (), 100)
    # Not the loopback at all.
    proc.listener("10.0.0.5", 9000, "11")
    assert not SG.owns_loopback(ipaddress.ip_address("10.0.0.5"), 9000, (), 100)


def test_a_second_listener_on_the_port_makes_it_somebody_elses(proc):
    """What a live host does not hand out on demand: the command's own
    server and another program both reachable at that address."""
    proc.process(100, socket_inodes=["11"])
    proc.process(300, socket_inodes=["22"], readable=False)
    proc.listener("127.0.0.1", 9000, "11").listener("0.0.0.0", 9000, "22")
    assert SG.loopback_listeners(LOOPBACK, 9000) == {"11", "22"}
    assert not SG.owns_loopback(LOOPBACK, 9000, (), 100)
    # Both named: every listener there is accounted for.
    assert SG.owns_loopback(LOOPBACK, 9000, {"11", "22"}, 100)


def test_a_wildcard_and_an_ipv6_listener_are_reachable_there(proc):
    proc.process(100, socket_inodes=["11", "12"])
    proc.listener("0.0.0.0", 9100, "11")
    assert SG.loopback_listeners(LOOPBACK, 9100) == {"11"}
    proc.listener("::1", 9200, "12")
    assert SG.loopback_listeners(ipaddress.ip_address("::1"), 9200) == {"12"}
    assert SG.owns_loopback(ipaddress.ip_address("::1"), 9200, (), 100)


# -- the proxy ---------------------------------------------------------------

class _Perms:
    task_session_id = "sess-own"
    mode = "bypassPermissions"
    confirm_callback = None


def test_a_post_to_a_background_jobs_server_is_not_an_upload(monkeypatch):
    """Before, ``decide`` let the session's own port through and ``settle``
    then refused every write to it with nobody to ask."""
    monkeypatch.setattr(BJ, "own_listening_ports",
                        lambda sid: {8393} if sid == "sess-own" else set())
    p = EP.EgressProxy(allowed_domains=(), allowed_ports=(80, 443))
    try:
        assert p.settle("127.0.0.1", 8393, method="POST", body_bytes=10,
                        perms=_Perms()) is None
        assert p.settle("127.0.0.1", 8394, method="POST", body_bytes=10,
                        perms=_Perms())
        # The process token (no session) gets no loopback.
        assert p.settle("127.0.0.1", 8393, method="POST", body_bytes=10)
    finally:
        p.close()


def test_the_proxy_asks_the_pid_that_forwarded_the_request(proc, monkeypatch):
    monkeypatch.setattr(BJ, "own_listening_ports", lambda sid: set())
    proc.process(100, socket_inodes=["11"]).listener("127.0.0.1", 9300, "11")
    # The helper forwards for the command below it; its own sockets, such
    # as the forwarder, do not make a port the command's.
    proc.process(101, ppid=100, socket_inodes=["11"])
    assert EP._own_server(_Perms(), "127.0.0.1", 9300, peer_pid=100)
    assert not EP._own_server(_Perms(), "127.0.0.1", 9300, peer_pid=200)
    p = EP.EgressProxy(allowed_domains=(), allowed_ports=(80, 443))
    try:
        assert p.decide("127.0.0.1", 9300, perms=_Perms(),
                        peer_pid=100)[0] == "allow"
    finally:
        p.close()


# -- the bash tool -----------------------------------------------------------

def test_the_bash_tool_names_the_sessions_job_listeners(monkeypatch):
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(BJ, "own_listening_sockets",
                        lambda sid: {"4242"} if sid == "sess-own" else set())
    args, _dir = A._sandbox_network_args(_Perms(), allow_slurm=False)
    assert args[args.index("--own-listener") + 1] == "4242", args


def test_a_jobs_listener_inodes_are_read_from_its_processes(proc):
    proc.process(400, socket_inodes=["31", "32"])
    proc.process(401, ppid=400, socket_inodes=["33"])
    proc.listener("127.0.0.1", 9400, "31").listener("127.0.0.1", 9401, "33")
    proc.connection("127.0.0.1", 9999, "32")
    assert BJ._listening_of({400, 401}) == {"31": 9400, "33": 9401}
    assert BJ._listening_ports_of({400, 401}) == {9400, 9401}
