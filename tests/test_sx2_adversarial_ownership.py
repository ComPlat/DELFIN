"""Adversarial tests for the SX2 loopback ownership rule (reviewer nacht-s2).

Attacks the ownership functions from the angles the builder named and
named himself: the /proc read race, the parent-chain trust, the
non-loopback and error paths, port-shape edge cases, and the boundary
between own-listener inodes and the process tree.

These tests never touch a real /proc: they point SG._PROC at a tmp
tree, as the builder's own test does. On a real machine nothing here
runs privileged.
"""

from __future__ import annotations

import ipaddress
import os
from pathlib import Path

import pytest

from delfin.agent import socket_guard as SG

LOOPBACK = ipaddress.ip_address("127.0.0.1")


def _hex_address(ip: str, port: int) -> str:
    raw = ipaddress.ip_address(ip).packed
    if ":" in ip:
        # the kernel prints each 32-bit word of an IPv6 address in host
        # byte order, so reverse each 4-byte group (little-endian host)
        words = b"".join(raw[i:i + 4][::-1] for i in range(0, 16, 4))
        return words.hex().upper() + f":{port:04X}"
    word = int.from_bytes(raw, "little")
    return f"{word:08X}:{port:04X}"


def _row(ip: str, port: int, inode: str, state: str = "07") -> str:
    return (f"   0: {_hex_address(ip, port)} 00000000:0000 {state} "
            f"00000000:0000 00:00000000 00000000     0        0 "
            f"{inode} 1 0000000000000000 100 0 0 10 0\n")


class ProcTree:
    def __init__(self, root: Path):
        self.root = root
        self.rows = {"udp": [], "udp6": [], "tcp": [], "tcp6": []}

    def process(self, pid, ppid=1, socket_inodes=(), readable=True):
        d = self.root / str(pid)
        d.mkdir(parents=True, exist_ok=True)
        stat = f"{pid} (proc) S {ppid} 1 1 0 0 0 0 0 0 0\n"
        (d / "stat").write_text(stat)
        fd = d / "fd"
        fd.mkdir(exist_ok=True)
        if not readable:
            fd.chmod(0o000)
            return self
        for i, inode in enumerate(socket_inodes):
            (fd / str(i)).symlink_to(f"socket:[{inode}]")
        return self

    def bound(self, ip: str, port: int, inode: str, table: str = None):
        if table is None:
            table = "udp6" if ":" in ip else "udp"
        self.rows[table].append(_row(ip, port, inode))
        return self

    def write(self):
        for table, rows in self.rows.items():
            path = self.root / "net" / table
            if rows:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(
                    "  sl  local_address rem_address   st\n" + "".join(rows))


@pytest.fixture
def tree(tmp_path, monkeypatch):
    t = ProcTree(tmp_path / "proc")
    monkeypatch.setattr(SG, "_PROC", str(t.root))
    yield t
    for d in t.root.glob("*/fd"):
        try:
            d.chmod(0o755)
        except OSError:
            pass


def test_proc_net_address_column_shapes(tree):
    """The address parser must not silently accept a malformed column:
    a too-short or non-hex field must be a ValueError, not a crash or a
    wrong address."""
    malformed = [
        "ZZZZZZZZ:1F90",
        "0100007F",
        "0100007F:ZZZZ",
        "",
        "0100007F:1F90:extra",
    ]
    for column in malformed:
        with pytest.raises(ValueError):
            SG._proc_net_address(column)


def test_wildcard_listener_counts_as_reachable(tree):
    """A listener bound to 0.0.0.0 (wildcard) is reachable from 127.0.0.1
    and must be counted in the reach set -- both for TCP and UDP."""
    tree.process(100, socket_inodes=["77"])
    tree.bound("0.0.0.0", 5001, "77").write()
    assert SG.loopback_listeners(LOOPBACK, 5001) == set()  # no tcp table rows written
    tree.rows["tcp"] = [_row("0.0.0.0", 5001, "77", state="0A")]
    tree.write()
    assert SG.loopback_listeners(LOOPBACK, 5001) == {"77"}
    assert not SG.owns_loopback(LOOPBACK, 5001, (), 0)
    assert SG.owns_loopback(LOOPBACK, 5001, {"77"}, 0)
    tree.rows["udp"] = [_row("0.0.0.0", 5001, "77")]
    tree.write()
    assert SG.udp_datagram_ends(LOOPBACK, 5001) == {"77"}


def test_ipv6_mapped_listener_on_loopback_port(tree):
    """An IPv6 socket bound to ::FFFF:127.0.0.1 (mapped) must be seen as
    the loopback address, and ::1 counts as loopback too."""
    tree.process(100, socket_inodes=["88"])
    tree.rows["tcp6"] = [
        _row("::FFFF:127.0.0.1", 5002, "88", state="0A"),
        _row("::1", 5003, "88", state="0A"),
    ]
    tree.write()
    assert SG.loopback_listeners(LOOPBACK, 5002) == {"88"}
    assert SG.loopback_listeners(ipaddress.ip_address("::1"), 5003) == {"88"}
    assert SG.owns_loopback(LOOPBACK, 5002, {"88"}, 0)
    assert SG.owns_loopback(ipaddress.ip_address("::1"), 5003, {"88"}, 0)


def test_parent_chain_trust(tree):
    """The parent closure trusts /proc stat. A process that LIES about
    its parent (points at a pid outside the tree that owns the socket)
    gains nothing: the socket inode is not in the tree's set."""
    tree.process(500, ppid=1, socket_inodes=["91"])
    tree.process(600, ppid=500, socket_inodes=())
    # a fake tree-root with a descendant claiming a foreign parent
    tree.process(700, ppid=600, socket_inodes=["92"])
    tree.bound("127.0.0.1", 5004, "91").write()
    assert not SG.owns_udp(LOOPBACK, 5004, (), 600)
    # but from the true root (500) it is owned
    assert SG.owns_udp(LOOPBACK, 5004, (), 500)
    # a descendant's own socket is in the tree from 600 as well
    tree.bound("127.0.0.1", 5005, "92").write()
    assert SG.owns_udp(LOOPBACK, 5005, (), 600)


def test_descendant_with_hidden_fds_is_not_owned(tree):
    """A descendant whose fds are unreadable (another user's rank inside
    the tree is impossible, but a hostile setuid child is not): its
    sockets cannot be shown, so the ownership check must FAIL CLOSED."""
    tree.process(800, socket_inodes=())
    tree.process(810, ppid=800, socket_inodes=["95"], readable=False)
    tree.bound("127.0.0.1", 5006, "95").write()
    assert not SG.owns_udp(LOOPBACK, 5006, (), 800)
    assert not SG.owns_loopback(LOOPBACK, 5006, (), 800)


def test_include_root_semantics(tree):
    """include_root=False excludes the root pid's own sockets: only the
    descendants' sockets count as own."""
    tree.process(900, socket_inodes=["101"])
    tree.process(910, ppid=900, socket_inodes=["102"])
    tree.bound("127.0.0.1", 5007, "101").write()
    assert SG.owns_udp(LOOPBACK, 5007, (), 900)
    assert not SG.owns_udp(LOOPBACK, 5007, (), 900, include_root=False)
    tree.bound("127.0.0.1", 5008, "102").write()
    assert SG.owns_udp(LOOPBACK, 5008, (), 900, include_root=False)


def test_mixed_own_and_foreign_endpoints_fail_closed(tree):
    """One endpoint mine, one foreign (bound to the same port by an
    unreadable process): the reach set is not a subset of own -> refuse.
    This is the shared-host case: another user's rank on 127.0.0.1."""
    tree.process(1000, socket_inodes=["111"])
    tree.process(1100, ppid=1, socket_inodes=(), readable=True)
    tree.bound("127.0.0.1", 5009, "111")
    tree.bound("127.0.0.1", 5009, "999")
    tree.write()
    assert SG.udp_datagram_ends(LOOPBACK, 5009) == {"111", "999"}
    assert not SG.owns_udp(LOOPBACK, 5009, (), 1000)


def test_port_zero_and_huge_port(tree):
    """Port 0 (ephemeral query) and a huge port number must not crash
    and must behave like any other port."""
    tree.process(1200, socket_inodes=["121"])
    tree.bound("127.0.0.1", 0, "121")
    tree.bound("127.0.0.1", 65535, "121")
    tree.write()
    assert SG.owns_udp(LOOPBACK, 0, (), 1200)
    assert SG.owns_udp(LOOPBACK, 65535, (), 1200)
    assert not SG.owns_udp(LOOPBACK, 0, (), 1200, include_root=False)


def test_owns_functions_never_raise(tree, monkeypatch):
    """The public ownership functions are used inside a seccomp
    supervisor path: they must return False, never raise, on any
    internal error."""
    monkeypatch.setattr(SG, "_PROC", "/nonexistent-proc-root-9x7y")
    assert not SG.owns_loopback(LOOPBACK, 5010, (), 1234)
    assert not SG.owns_udp(LOOPBACK, 5010, (), 1234)
    assert SG.loopback_listeners(LOOPBACK, 5010) == set()
    assert SG.udp_datagram_ends(LOOPBACK, 5010) == set()
    assert SG.socket_inodes_below(1234) == set()


def test_own_inodes_mixed_types(tree):
    """own_inodes given as ints and strings mixed: the set comparison
    must normalise to strings (the Supervisor keeps digit strings)."""
    tree.process(1300, socket_inodes=())
    tree.bound("127.0.0.1", 5011, "131").write()
    assert SG.owns_udp(LOOPBACK, 5011, (131,), 0)
    assert SG.owns_udp(LOOPBACK, 5011, ("131",), 0)
    assert SG.owns_udp(LOOPBACK, 5011, [131, "141"], 0)


def test_non_digit_inode_strings_in_own_set(tree):
    """Garbage in the own set must not match a real inode and must not
    crash (the Supervisor filters non-digit strings, but a caller may
    not)."""
    tree.process(1400, socket_inodes=())
    tree.bound("127.0.0.1", 5012, "141").write()
    assert not SG.owns_udp(LOOPBACK, 5012, ("socket:[141]",), 0)
    assert SG.owns_udp(LOOPBACK, 5012, ("141",), 0)


def test_stat_parsing_survives_brackets_and_spaces(tree):
    """A process name with ')' and spaces inside must not break the
    parent-chain parse (rsplit on the last ')'): with a well-formed
    kernel stat line the child is recognised; with a line the parser
    cannot read the child is SKIPPED, so the check fails closed rather
    than adopting a guessed parent."""
    d = tree.root / "1500"
    d.mkdir(parents=True, exist_ok=True)
    # well-formed kernel stat line: after the comm bracket the fields are
    # state then ppid, exactly what rsplit(")",1)[1].split()[1] reads
    (d / "stat").write_text("1500 (comm name) S 1300 1 1 0 0 0 0 0 0 0\n")
    (d / "fd").mkdir(exist_ok=True)
    os.symlink("socket:[151]", str(d / "fd" / "0"))
    tree.process(1300, socket_inodes=())
    tree.bound("127.0.0.1", 5013, "151").write()
    assert SG.owns_udp(LOOPBACK, 5013, (), 1300)
    # an unparseable stat line (extra field shifts state/ppid): the
    # parser must SKIP the process, not adopt a wrong parent
    (d / "stat").write_text(
        "1500 (a process) with ) brackets S 1300 1 1 0 0 0 0 0 0 0\n")
    tree.bound("127.0.0.1", 5014, "152").write()
    assert not SG.owns_udp(LOOPBACK, 5014, (), 1300)
