"""The datagram (UDP) endpoints a sandboxed MPI command may use.

Measured in the JC1/JC2/JC3 live runs (14.-15.11.2026, SLURM jobs): OpenMPI
inside the agent's command cage failed its interface discovery before any
connect() -- "opal_ifinit: socket() failed with errno=13". The seccomp
filter of a restricted network refuses every non-stream INET socket
(socket_guard._program, block_inet_datagrams), because a UDP datagram can
leave to the real network without ever calling connect(); the guard's
connect-supervision never sees it. The ownership rule of the loopback
(#742) judges TCP listeners only (/proc/net/tcp{,6}); a UDP endpoint is
judged from /proc/net/udp{,6}, whose bound endpoint is a socket the
process tree holds.

The rule checked here, on a RECORDED process and socket tree so every
test runs on every machine: a datagram INET socket may be created when
the network mode is restricted ONLY for loopback addresses (127.0.0.0/8
and ::1 -- an MPI rank pair talks over the loopback it was given), and
the send/recv endpoints of the same command's process tree. A UDP socket
bound to a non-loopback address is refused at creation, so nothing can
leave past the proxy. The TCP ownership rule (owns_loopback) is
untouched.
"""

from __future__ import annotations

import ipaddress

import pytest

from delfin.agent import socket_guard as SG


LOOPBACK = ipaddress.ip_address("127.0.0.1")


def _hex_address(ip: str, port: int) -> str:
    raw = ipaddress.ip_address(ip).packed
    if SG.sys.byteorder == "little":
        raw = b"".join(raw[i:i + 4][::-1] for i in range(0, len(raw), 4))
    return f"{raw.hex().upper()}:{port:04X}"


def _row(ip: str, port: int, inode: str) -> str:
    return (f"   0: {_hex_address(ip, port)} 00000000:0000 07 "
            f"00000000:00000000 00:00000000 00000000  1000  0 {inode} 1\n")


class UdpTree:
    """A recorded /proc with UDP tables and processes with socket fds."""

    def __init__(self, root):
        self.root = root
        (root / "net").mkdir(parents=True, exist_ok=True)
        self.rows = {"udp": [], "udp6": []}

    def process(self, pid, ppid=1, socket_inodes=(), readable=True):
        d = self.root / str(pid)
        d.mkdir(exist_ok=True)
        (d / "stat").write_text(f"{pid} (mpirun) S {ppid} " + "0 " * 20 + "\n")
        fd = d / "fd"
        fd.mkdir(exist_ok=True)
        if not readable:
            fd.chmod(0o000)
            return self
        for i, inode in enumerate(socket_inodes):
            (fd / str(i)).symlink_to(f"socket:[{inode}]")
        return self

    def bound(self, ip: str, port: int, inode: str):
        self.rows["udp6" if ":" in ip else "udp"].append(_row(ip, port, inode))
        return self

    def write(self):
        for table, rows in self.rows.items():
            (self.root / "net" / table).write_text(
                "  sl  local_address rem_address   st\n" + "".join(rows))


@pytest.fixture
def udp_tree(tmp_path, monkeypatch):
    tree = UdpTree(tmp_path / "proc")
    monkeypatch.setattr(SG, "_PROC", str(tree.root))
    yield tree
    for d in tree.root.glob("*/fd"):
        d.chmod(0o755)


def test_a_udp_endpoint_of_the_command_tree_is_its_own(udp_tree):
    """An MPI rank bound to the loopback, the socket held by a process of
    the command tree: it may use its endpoint."""
    udp_tree.process(100, socket_inodes=["11"])
    udp_tree.process(101, ppid=100, socket_inodes=["12"])
    udp_tree.bound("127.0.0.1", 4001, "12").write()
    assert SG.udp_datagram_ends(LOOPBACK, 4001) == {"12"}
    assert SG.owns_udp(LOOPBACK, 4001, (), 100)


def test_a_foreign_or_unheld_udp_endpoint_is_not_its_own(udp_tree):
    """Another user's rank on the same port (its fds unreadable), and a
    bound endpoint whose holder is not in the tree: both stay refused."""
    udp_tree.process(100)
    udp_tree.process(300, socket_inodes=(), readable=False)
    udp_tree.bound("127.0.0.1", 4002, "33").write()
    assert SG.udp_datagram_ends(LOOPBACK, 4002) == {"33"}
    assert not SG.owns_udp(LOOPBACK, 4002, (), 100)
    # Named inodes (a session's background jobs) count as own, as for TCP.
    assert SG.owns_udp(LOOPBACK, 4002, {"33"}, 0)
    # Nothing bound there: a connect/refusal question, not an owner.
    assert not SG.owns_udp(LOOPBACK, 4003, (), 100)


def test_a_non_loopback_udp_endpoint_is_never_owned(udp_tree):
    """A datagram endpoint bound to a real interface address is what the
    filter must refuse at creation: it is how a datagram leaves past the
    proxy."""
    udp_tree.process(100, socket_inodes=["11"])
    udp_tree.bound("192.168.1.10", 4004, "11").write()
    assert not SG.owns_udp(ipaddress.ip_address("192.168.1.10"), 4004, (), 100)


def test_established_and_tcp_rows_are_not_udp(udp_tree):
    """The UDP tables' state 07 is the bound entry; a TCP table row or an
    ESTABLISHED-style row must not be mistaken for a UDP endpoint."""
    udp_tree.process(100, socket_inodes=["11"])
    (udp_tree.root / "net" / "udp").write_text(
        "  sl  local_address rem_address   st\n"
        + _row("127.0.0.1", 4005, "41").replace(" 07 ", " 01 ")
        + _row("127.0.0.1", 4006, "42"))
    assert SG.udp_datagram_ends(LOOPBACK, 4005) == set()
    assert SG.udp_datagram_ends(LOOPBACK, 4006) == {"42"}
