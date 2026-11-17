"""A sandboxed command cannot reach the user's sessions through a socket.

Landlock holds a command to its folders, but it does not govern connecting
to a Unix socket. Through one a command reaches what runs OUTSIDE every
sandbox with the user's rights: a tmux or screen server (``send-keys`` types
into the user's terminal), the SSH agent, the systemd user manager
(``systemd-run --user``), a desktop session, the Docker daemon. Bubblewrap
hides those sockets behind a mount namespace; without namespaces the files
cannot be hidden, so the connection itself is judged.

Standalone (stdlib only), used by ``landlock_exec`` in the process it
starts the command from. The command runs under a seccomp filter that
hands every ``connect()`` to this supervisor (``SECCOMP_RET_USER_NOTIF``).
The supervisor copies the address out of the command's memory, decides,
and performs the connect ITSELF on a duplicate of the command's socket
(``pidfd_getfd``) with the address it copied -- so a thread that rewrites
the address after the check changes nothing. Unix sockets are allowed only
beneath an allow-list (the workspace, the private temp directory, the
cluster's authentication and name services); abstract Unix sockets are
refused; other address families are connected as asked. The filter also
refuses datagram Unix sockets (``sendto`` carries its own address) and
io_uring (which connects without the system call).

Under the proxy a TCP connect to the loopback is refused unless every
socket listening there is the command's own: held by its process tree or
named as a listener of the session's own background jobs
(``Supervisor._own_loopback``). A test that serves a fake endpoint to
itself works; the user's notebook or model server and other users'
services on the same host stay out of reach.

If the supervisor ends, the listener closes and every later connect fails:
the guard fails closed. Linux x86_64 and aarch64; ``available()`` says
whether this kernel and architecture support it.
"""
from __future__ import annotations

import ctypes
import errno
import ipaddress
import os
import platform
import select
import socket
import struct
import sys
import threading

_ARCHES = {
    # machine: (AUDIT_ARCH, seccomp, connect, socket, io_uring_setup,
    #           pidfd_open, pidfd_getfd)
    # The datagram surface (bind, sendto, sendmsg, sendmmsg) lives in
    # _ARCH_DGRAM: _arch() keeps its original 7-tuple contract, which a
    # main-line test unpacks by position.
    "x86_64": (0xC000003E, 317, 42, 41, 425, 434, 438),
    "aarch64": (0xC00000B7, 277, 203, 198, 425, 434, 438),
}

#: machine: (bind, sendto, sendmsg, sendmmsg)
# From the kernel headers (operator live-run #3, recvfrom(45) was
# refused and sendmmsg missed with the old rows): x86_64
# asm/unistd_64.h -- bind 49, sendto 44, sendmsg 46, sendmmsg 307;
# aarch64 asm-generic -- bind 200, sendto 206, sendmsg 211, sendmmsg
# 269. Guarded by tests/test_sx2_syscall_numbers.py (header table +
# name mapping), never hand-copied again.
_ARCH_DGRAM = {
    "x86_64": (49, 44, 46, 307),
    "aarch64": (200, 206, 211, 269),
}


def _arch_dgram():
    return _ARCH_DGRAM.get(platform.machine())

_SECCOMP_SET_MODE_FILTER = 1
_SECCOMP_GET_NOTIF_SIZES = 3
_SECCOMP_FILTER_FLAG_NEW_LISTENER = 1 << 3
_RET_ALLOW = 0x7FFF0000
_RET_ERRNO = 0x00050000
_RET_USER_NOTIF = 0x7FC00000
_NOTIF_RECV = 0xC0502100
_NOTIF_SEND = 0xC0182101
_NOTIF_ID_VALID = 0x40082102
_NOTIF_ID_VALID_OLD = 0x80082102
#: SECCOMP_USER_NOTIF_FLAG_CONTINUE: the kernel carries out the notified
#: syscall itself after the supervisor answers 0 with this flag.
_NOTIF_FLAG_CONTINUE = 1

_libc = ctypes.CDLL(None, use_errno=True)
_libc.syscall.restype = ctypes.c_long
_libc.ioctl.restype = ctypes.c_int
_libc.connect.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
_libc.connect.restype = ctypes.c_int

#: Beneath these a Unix socket may be reached from any sandbox: the cluster
#: authentication (munge, for sbatch/squeue) and the name services (sssd,
#: nscd), without which user names do not resolve.
DEFAULT_ALLOWED = (
    "/run/munge", "/var/run/munge", "/var/lib/sss/pipes", "/run/nscd",
    "/var/run/nscd",
)

#: systemd-resolved's sockets answer host-name lookups only, and a lookup
#: carries its name to the DNS: ``<data>.attacker.example`` reaches the
#: attacker's name server. A command behind the proxy needs no lookup (the
#: proxy resolves), and its inet datagrams are refused, so this is offered
#: only to a command with a direct route (SLURM controller, license
#: server, SSH remote) on a host whose nsswitch resolves names through
#: systemd-resolved, where the program could not reach that route
#: otherwise (``name_lookup_sockets``). Measured 2026-11-13: allowed for
#: every command before, a caged command connected to
#: /run/systemd/resolve/io.systemd.Resolve.
NAME_LOOKUP_SOCKETS = ("/run/systemd/resolve",)


def name_lookup_sockets(direct_routes: bool, net_mode: str,
                        nsswitch: str = "/etc/nsswitch.conf") -> tuple:
    """``NAME_LOOKUP_SOCKETS`` when a command needs them, else ().

    Needed: an open network (the command can send DNS anyway), or direct
    routes on a host whose ``hosts:`` line in ``nsswitch`` names
    ``resolve``. Otherwise a lookup through them is only a way to carry
    a name out."""
    if net_mode == "open":
        return NAME_LOOKUP_SOCKETS
    if not direct_routes:
        return ()
    try:
        with open(nsswitch, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                words = line.split("#", 1)[0].split()
                if words[:1] == ["hosts:"] and "resolve" in words[1:]:
                    return NAME_LOOKUP_SOCKETS
    except OSError:
        pass
    return ()


#: Metadata services that hand out cloud credentials; link-local ranges are
#: refused as a whole besides these.
_METADATA = {ipaddress.ip_address(a) for a in (
    "169.254.169.254", "100.100.100.200", "168.63.129.16", "fd00:ec2::254")}


#: Where the process and socket tables are read. A test points it at a
#: recorded tree, so the ownership rule is checked on every machine.
_PROC = "/proc"


def _proc_net_address(field: str):
    """``(ip_address, port)`` of a ``local_address`` in /proc/net/tcp{,6}.

    The kernel prints each 32-bit word of the address in host byte order."""
    hexaddr, hexport = field.split(":")
    raw = bytes.fromhex(hexaddr)
    if sys.byteorder == "little":
        raw = b"".join(raw[i:i + 4][::-1] for i in range(0, len(raw), 4))
    return ipaddress.ip_address(raw), int(hexport, 16)


def loopback_listeners(ip, port: int) -> set:
    """Inodes of the listening TCP sockets a connect to ``ip:port`` can
    reach: bound to that address or to a wildcard, IPv4 and IPv6 (an IPv6
    socket that also accepts IPv4 shows its address as mapped). Read from
    /proc/net, which shows every user's sockets in this network namespace;
    empty when nothing listens there or the tables cannot be read."""
    out: set = set()
    for table in (f"{_PROC}/net/tcp", f"{_PROC}/net/tcp6"):
        try:
            with open(table, encoding="latin-1") as fh:
                next(fh, None)
                for line in fh:
                    cols = line.split()
                    # local_address, state 0A = LISTEN, inode at index 9
                    if len(cols) < 10 or cols[3] != "0A":
                        continue
                    try:
                        addr, lport = _proc_net_address(cols[1])
                    except ValueError:
                        continue
                    if lport != port:
                        continue
                    if getattr(addr, "ipv4_mapped", None) is not None:
                        addr = addr.ipv4_mapped
                    if addr == ip or addr.is_unspecified:
                        out.add(cols[9])
        except OSError:
            continue
    return out


def socket_inodes_below(root: int, include_root: bool = True) -> set:
    """Inodes of the sockets ``root`` and every process below it hold
    (/proc); without ``include_root``, those of the processes below it only.
    A process of another user is not readable and adds nothing."""
    tree = {root}
    parent_of: dict = {}
    try:
        names = os.listdir(_PROC)
    except OSError:
        names = []
    for name in names:
        if not name.isdigit():
            continue
        try:
            with open(f"{_PROC}/{name}/stat", encoding="latin-1") as fh:
                stat = fh.read()
            parent_of[int(name)] = int(stat.rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            continue
    grew = True
    while grew:
        grew = False
        for pid, ppid in parent_of.items():
            if ppid in tree and pid not in tree:
                tree.add(pid)
                grew = True
    if not include_root:
        tree.discard(root)
    inodes: set = set()
    for pid in tree:
        try:
            fds = os.listdir(f"{_PROC}/{pid}/fd")
        except OSError:
            continue
        for fd in fds:
            try:
                link = os.readlink(f"{_PROC}/{pid}/fd/{fd}")
            except OSError:
                continue
            if link.startswith("socket:["):
                inodes.add(link[8:-1])
    return inodes


def owns_loopback(ip, port: int, own_inodes=(), root_pid: int = 0,
                  include_root: bool = True) -> bool:
    """Whether every socket a connect to loopback ``ip:port`` can reach is
    the caller's own: one of ``own_inodes`` (the listening sockets of the
    session's background jobs) or held by ``root_pid``'s process tree.

    On a shared host the loopback carries other users' services and the
    user's own unsandboxed ones (a notebook server, a model server, the
    dashboard); those are in neither and stay refused. Judged from /proc
    when asked, so a port the server gave up and another program took is no
    longer the caller's. False when nothing listens there, for a non-loopback
    address and on any error."""
    try:
        if not ip.is_loopback:
            return False
        reach = loopback_listeners(ip, port)
        if not reach:
            return False
        own = {str(i) for i in own_inodes}
        if reach <= own:
            return True
        return bool(root_pid) and reach <= (
            own | socket_inodes_below(int(root_pid), include_root))
    except Exception:
        return False


def _udp_datagram_rows():
    """The bound endpoint rows of the kernel's UDP tables.

    State 07 is UDP's bound entry (an unconnected socket's resting state,
    which is how OpenMPI holds its discovery and rank-to-rank sockets)."""
    out: set = set()
    for table in (f"{_PROC}/net/udp", f"{_PROC}/net/udp6"):
        try:
            with open(table, encoding="latin-1") as fh:
                next(fh, None)
                for line in fh:
                    cols = line.split()
                    # local_address, state 07 = the bound endpoint,
                    # inode at index 9
                    if len(cols) < 10 or cols[3] != "07":
                        continue
                    try:
                        addr, lport = _proc_net_address(cols[1])
                    except ValueError:
                        continue
                    if getattr(addr, "ipv4_mapped", None) is not None:
                        addr = addr.ipv4_mapped
                    out.add((addr, lport, cols[9]))
        except OSError:
            continue
    return out


def udp_datagram_ends(ip, port: int) -> set:
    """Inodes of the bound UDP endpoints a datagram to ``ip:port`` can
    reach: bound to that address or to a wildcard, IPv4 and IPv6 (an IPv6
    socket that also serves IPv4 shows its address as mapped). Read from
    /proc/net, which shows every user's sockets in this network namespace;
    empty when nothing is bound there or the tables cannot be read."""
    return {(inode) for (addr, lport, inode) in _udp_datagram_rows()
            if lport == port
            and (addr == ip or addr.is_unspecified)}


def owns_udp(ip, port: int, own_inodes=(), root_pid: int = 0,
             include_root: bool = True) -> bool:
    """Whether every UDP endpoint a datagram to ``ip:port`` can reach is
    the caller's own: one of ``own_inodes`` (the session's background
    jobs' sockets) or held by ``root_pid``'s process tree.

    Same rule as :func:`owns_loopback`, on the UDP tables: on a shared
    host the loopback carries other users' ranks and the user's own
    unsandboxed senders; those are in neither and stay refused. False
    when nothing is bound there, for a non-loopback address (a datagram
    bound to a real interface is what the filter must refuse at
    creation) and on any error."""
    try:
        if not ip.is_loopback:
            return False
        reach = udp_datagram_ends(ip, port)
        if not reach:
            return False
        own = {str(i) for i in own_inodes}
        if reach <= own:
            return True
        return bool(root_pid) and reach <= (
            own | socket_inodes_below(int(root_pid), include_root))
    except Exception:
        return False


def _arch():
    return _ARCHES.get(platform.machine())


class _SockFilter(ctypes.Structure):
    _fields_ = [("code", ctypes.c_uint16), ("jt", ctypes.c_uint8),
                ("jf", ctypes.c_uint8), ("k", ctypes.c_uint32)]


class _SockFprog(ctypes.Structure):
    _fields_ = [("len", ctypes.c_uint16), ("filter", ctypes.POINTER(_SockFilter))]


def available() -> bool:
    """Whether the kernel offers user notification and this architecture is
    supported. Does not change the process."""
    a = _arch()
    if a is None or not sys.platform.startswith("linux") or sys.byteorder != "little":
        return False
    sizes = (ctypes.c_uint16 * 3)()
    rc = _libc.syscall(ctypes.c_long(a[1]), ctypes.c_uint(_SECCOMP_GET_NOTIF_SIZES),
                       ctypes.c_uint(0), ctypes.byref(sizes))
    return rc == 0 and sizes[0] >= 80 and sizes[1] >= 24


def _assemble(items) -> list[tuple[int, int, int, int]]:
    """Resolve labelled jumps. An item is ``("label", name)``, a statement
    ``(code, k)``, or a jump ``(code, k, true_label, false_label)`` where a
    label of None means the next instruction."""
    labels, count = {}, 0
    for item in items:
        if item[0] == "label":
            labels[item[1]] = count
        else:
            count += 1
    out, pc = [], 0
    for item in items:
        if item[0] == "label":
            continue
        if len(item) == 4:
            code, k, t, f = item
            jt = 0 if t is None else labels[t] - pc - 1
            jf = 0 if f is None else labels[f] - pc - 1
            if not (0 <= jt <= 255 and 0 <= jf <= 255):
                raise ValueError("jump out of range")
            out.append((code, jt, jf, k))
        else:
            out.append((item[0], 0, 0, item[1]))
        pc += 1
    return out


_LD, _JEQ, _JSET, _RET, _AND = 0x20, 0x15, 0x45, 0x06, 0x54


def _head(audit: int):
    deny = _RET_ERRNO | errno.EPERM
    items = [(_LD, 4), (_JEQ, audit, None, "deny_all"), (_LD, 0)]
    if audit == 0xC000003E:
        # x32 system calls carry bit 30; refuse them all.
        items.append((_JSET, 0x40000000, "deny_all", None))
    return items, [("label", "deny_all"), (_RET, deny)]


def _program(block_inet_datagrams: bool = False):
    """connect() to the supervisor, no io_uring, no datagram Unix socket,
    and with a restricted network the datagram surface is judged, not
    blocked: socket(AF_INET/6, SOCK_DGRAM) is notified (creation sends
    nothing; OpenMPI reads interfaces over it), and the address-bearing
    datagram sends (sendto/sendmsg/sendmmsg) are notified too -- a UDP
    datagram otherwise leaves without connect(), past the supervisor.
    A foreign architecture (32-bit entry points use other numbers)
    runs nothing."""
    (audit, _sc, nr_connect, nr_socket, nr_uring, _po, _pg) = _arch()
    (nr_bind, nr_sendto, nr_sendmsg, nr_sendmmsg) = _arch_dgram()
    items, tail = _head(audit)
    items += [
        (_JEQ, nr_connect, "notify", None),
        (_JEQ, nr_uring, "enosys", None),
    ]
    if block_inet_datagrams:
        # The address-bearing datagram sends are judged BEFORE the
        # socket() comparison: sendto(nr) never matches nr_socket, so a
        # JEQ after socket()'s false edge was unreachable -- and with a
        # false edge to allow it was an open escape (reviewer walk).
        # sendto is judged only when the destination pointer (args[4],
        # seccomp_data args[1] of sendto) is non-NULL: a NULL
        # destination uses the connected peer, which connect() already
        # judged -- the supervisor stays out of the no-destination hot
        # path. sendmsg/sendmmsg carry the destination inside the
        # msghdr, invisible to BPF: they are refused in-kernel (the
        # measured OpenMPI path needs socket+ioctl and TCP connect
        # only; adding them back needs the pidfd_getfd handoff first).
        items += [(_JEQ, nr_sendto, None, "after_sends"),
                  (_LD, 40),                          # args[1]: dest ptr
                  (_JEQ, 0, "allow", "notify")]       # NULL -> peer judged, allow;
        # dest non-NULL -> notified and judged by _datagram_send.
        items += [("label", "after_sends"),
                  # sendmsg/sendmmsg: destination inside the msghdr,
                  # invisible to BPF -- refused in-kernel, type-agnostic
                  # (a datagram could leave past the proxy otherwise,
                  # reviewer walk on 9ab3816e9).
                  (_JEQ, nr_sendmsg, "eacces", None),
                  (_JEQ, nr_sendmmsg, "eacces", None)]
    # socket(): its false edge falls THROUGH to the family dispatch, it
    # must not terminate in allow (that would skip the sends above).
    items += [(_JEQ, nr_socket, None, None)]
    items += [
        (_LD, 16),
        (_JEQ, socket.AF_UNIX, "unix", None),
    ]
    if block_inet_datagrams:
        items += [(_JEQ, socket.AF_INET, "inet", None),
                  (_JEQ, socket.AF_INET6, "inet", None)]
    items += [(_RET, _RET_ALLOW),
              ("label", "unix"), (_LD, 24), (_AND, 0xF),
              (_JEQ, socket.SOCK_DGRAM, "eacces", "allow")]
    if block_inet_datagrams:
        # Datagram creation is notified; anything else (streams, the
        # proxy's TCP path) falls through to the allow label.
        items += [("label", "inet"), (_LD, 24), (_AND, 0xF),
                  (_JEQ, socket.SOCK_DGRAM, "notify", None)]
    items += [("label", "allow"), (_RET, _RET_ALLOW),
              ("label", "eacces"), (_RET, _RET_ERRNO | errno.EACCES),
              ("label", "enosys"), (_RET, _RET_ERRNO | errno.ENOSYS),
              ("label", "notify"), (_RET, _RET_USER_NOTIF)] + tail
    return _assemble(items)


def _block_unix_program(block_inet_datagrams: bool = False):
    """Where the kernel cannot hand connects to a supervisor: no Unix socket
    at all (and no io_uring), and with a restricted network no inet socket
    at all either."""
    audit, _sc, _nc, nr_socket, nr_uring, _po, _pg, *_x = _arch()
    items, tail = _head(audit)
    items += [(_JEQ, nr_uring, "enosys", None),
              (_JEQ, nr_socket, None, "allow"),
              (_LD, 16),
              (_JEQ, socket.AF_UNIX, "eacces", None)]
    if block_inet_datagrams:
        items += [(_JEQ, socket.AF_INET, "eacces", None),
                  (_JEQ, socket.AF_INET6, "eacces", None)]
    items += [("label", "allow"), (_RET, _RET_ALLOW),
              ("label", "eacces"), (_RET, _RET_ERRNO | errno.EACCES),
              ("label", "enosys"), (_RET, _RET_ERRNO | errno.ENOSYS)] + tail
    return _assemble(items)


def filter_supported() -> bool:
    return (_arch() is not None and sys.platform.startswith("linux")
            and sys.byteorder == "little")


def _install(insns, flags: int) -> int:
    a = _arch()
    arr = (_SockFilter * len(insns))(*[_SockFilter(*i) for i in insns])
    prog = _SockFprog(len(insns), arr)
    fd = _libc.syscall(ctypes.c_long(a[1]), ctypes.c_uint(_SECCOMP_SET_MODE_FILTER),
                       ctypes.c_uint(flags), ctypes.byref(prog))
    if fd < 0:
        raise OSError(ctypes.get_errno(), "seccomp filter: " + os.strerror(ctypes.get_errno()))
    return int(fd)


def install_block_unix_filter(block_inet_datagrams: bool = False) -> None:
    _install(_block_unix_program(block_inet_datagrams), 0)


def install_filter(block_inet_datagrams: bool = False) -> int:
    """In the process that will run the command: install the filter and
    return the listener descriptor. Needs no_new_privs set already."""
    return _install(_program(block_inet_datagrams), _SECCOMP_FILTER_FLAG_NEW_LISTENER)


def parse_notification(raw: bytes) -> tuple[int, int, int, tuple]:
    """(id, pid, syscall number, args) of a struct seccomp_notif.

    Layout: u64 id @0, u32 pid @8, u32 flags @12, then struct seccomp_data
    @16 -- s32 nr @16, u32 arch @20, u64 instruction_pointer @24, u64
    args[6] @32. The number was read at offset 24 (the instruction
    pointer's low half), so no socket() or sendto() was ever recognised
    and each fell into the connect() path, which answered EINVAL (live run
    2026-11-17: OpenMPI's opal_ifinit in the agent's cage)."""
    notif_id, pid, _flags = struct.unpack_from("<QII", raw, 0)
    (nr,) = struct.unpack_from("<i", raw, 16)
    args = struct.unpack_from("<6Q", raw, 32)
    return int(notif_id), int(pid), int(nr), tuple(int(a) for a in args)


class Supervisor:
    """Answers the command's connect() calls until the command ends."""

    def __init__(self, listener: int, allowed: list[str], net_mode: str = "open",
                 allowed_tcp=(), denied=(), own_listeners=()):
        self.listener = listener
        # The listening sockets of the session's own background jobs (their
        # inodes, read when the command started); with the sockets of the
        # command's own process tree they are the servers it may reach on
        # the loopback (_own_loopback).
        self.own_listeners = {str(i) for i in own_listeners if str(i).isdigit()}
        self.root_pid = 0
        self.allowed = [os.path.realpath(p).rstrip("/") for p in allowed if p]
        self.denied = [os.path.realpath(p).rstrip("/") for p in denied if p]
        self.net_mode = net_mode if net_mode in ("open", "proxy", "none") else "none"
        self.allowed_tcp = set()
        for entry in allowed_tcp:
            host, _, port = str(entry).rpartition(":")
            try:
                self.allowed_tcp.add((ipaddress.ip_address(host.strip("[]")), int(port)))
            except ValueError:
                continue
        a = _arch()
        self._nr_pidfd_open, self._nr_pidfd_getfd = a[5], a[6]
        self.refused: list[str] = []

    # -- policy -------------------------------------------------------------
    def inet_refusal(self, ip, port: int) -> str:
        """Why a connect to ``ip:port`` is refused, or ""."""
        if getattr(ip, "ipv4_mapped", None) is not None:
            ip = ip.ipv4_mapped
        if ip.is_link_local or ip in _METADATA:
            return f"cloud metadata / link-local address {ip}"
        if self.net_mode == "open":
            return ""
        # An entry with port 0 admits every port of that address (a license
        # server whose vendor daemon listens on a port of its own).
        if self.net_mode == "proxy" and ((ip, port) in self.allowed_tcp
                                         or (ip, 0) in self.allowed_tcp):
            return ""
        if self.net_mode == "proxy" and ip.is_loopback and self._own_loopback(ip, port):
            return ""
        return (f"direct connection to {ip}:{port}; the sandbox reaches the "
                "network through its proxy only")

    def _own_loopback(self, ip, port: int) -> bool:
        """A server of the command's own tree or one of the session's own
        background jobs (``owns_loopback``). The supervisor connects right
        after the check (_handle); a socket replaced in between is not seen.

        Rejected: allowing the whole loopback under the proxy (reaches what
        runs outside every sandbox with the user's rights); a port list
        fixed when the command starts (a port given up and taken by another
        program stays open)."""
        return owns_loopback(ip, port, self.own_listeners, self.root_pid)

    def unix_path_allowed(self, path: str) -> bool:
        real = os.path.realpath(path)

        def beneath(prefixes):
            return any(real == p or real.startswith(p + "/") or p == ""
                       for p in prefixes)

        if self.denied and beneath(self.denied):
            return False
        return beneath(self.allowed)

    # -- plumbing -----------------------------------------------------------
    def _send(self, notif_id: int, val: int, err: int, flags: int = 0) -> None:
        # struct seccomp_notif_resp: id, val, error, flags. CONTINUE goes in
        # flags with val and error both 0 -- in the error field it made the
        # kernel answer the command's socket() with EINVAL (live run 17.11.).
        resp = struct.pack("<QqiI", notif_id, val, err, flags)
        buf = ctypes.create_string_buffer(resp, 24)
        _libc.ioctl(self.listener, ctypes.c_ulong(_NOTIF_SEND), buf)

    def _valid(self, notif_id: int) -> bool:
        idb = ctypes.c_uint64(notif_id)
        for req in (_NOTIF_ID_VALID, _NOTIF_ID_VALID_OLD):
            if _libc.ioctl(self.listener, ctypes.c_ulong(req), ctypes.byref(idb)) == 0:
                return True
        return False

    @staticmethod
    def _thread_group(tid: int) -> int:
        """The process a thread belongs to. A notification names the calling
        THREAD, and pidfd_open accepts only a process: a connect() from a
        resolver thread (curl, most HTTP clients) failed with EBADF."""
        try:
            with open(f"/proc/{tid}/status", encoding="ascii", errors="replace") as fh:
                for line in fh:
                    if line.startswith("Tgid:"):
                        return int(line.split()[1])
        except (OSError, ValueError, IndexError):
            pass
        return tid

    def _dup_fd(self, pid: int, fd: int) -> int:
        pidfd = _libc.syscall(ctypes.c_long(self._nr_pidfd_open),
                              ctypes.c_int(self._thread_group(pid)), ctypes.c_uint(0))
        if pidfd < 0:
            return -1
        try:
            return int(_libc.syscall(ctypes.c_long(self._nr_pidfd_getfd), ctypes.c_int(pidfd),
                                     ctypes.c_int(fd), ctypes.c_uint(0)))
        finally:
            os.close(pidfd)

    def _datagram_socket(self, notif_id: int, pid: int, args: tuple,
                         nr: int) -> None:
        """socket() under a restricted network: creating a datagram socket
        is the command's own act (creation sends nothing; OpenMPI reads
        its interfaces over one). Allowed with
        SECCOMP_USER_NOTIF_FLAG_CONTINUE, so the kernel carries out the
        real call -- no address rewrite, no copied descriptor.

        The escape is judged where the address appears: sendto/
        sendmsg/sendmmsg with a destination go to _datagram_send."""
        self._send(notif_id, 0, 0, _NOTIF_FLAG_CONTINUE)

    def _datagram_send(self, notif_id: int, pid: int, args: tuple) -> None:
        """sendto() with a destination address: the destination is judged
        by the UDP ownership rule (owns_udp -- the command tree's own
        ranks, or a session background job's socket, loopback only). The
        filter notifies sendto only when the destination pointer is
        non-NULL (a NULL destination uses the connected peer, which
        connect() already judged, and is allowed in-kernel);
        sendmsg/sendmmsg are refused in-kernel (the destination is
        inside the msghdr, invisible to BPF). Everything else is refused
        with EACCES, fail closed."""
        try:
            addr_ptr, addr_len = int(args[1]), int(args[2]) & 0xFFFFFFFF
            if addr_len < 2:
                # no destination reached the judge: the filter only
                # notifies sendto with a non-NULL pointer
                self._send(notif_id, -1, -errno.EACCES)
                return
            with open(f"/proc/{pid}/mem", "rb", buffering=0) as mem:
                mem.seek(addr_ptr)
                raw = mem.read(addr_len)
            if len(raw) != addr_len or not self._valid(notif_id):
                raise ValueError("unreadable destination address")
            family = struct.unpack_from("<H", raw)[0]
            if family == socket.AF_INET and len(raw) >= 8:
                port = struct.unpack_from(">H", raw, 2)[0]
                ip = ipaddress.ip_address(raw[4:8])
            elif family == socket.AF_INET6 and len(raw) >= 24:
                port = struct.unpack_from(">H", raw, 2)[0]
                ip = ipaddress.ip_address(raw[8:24])
            else:
                raise ValueError(f"address family {family}")
            if owns_udp(ip, port, self.own_listeners,
                        self.root_pid, include_root=True):
                self._send(notif_id, 0, 0, _NOTIF_FLAG_CONTINUE)
                return
            self.refused.append(f"send {ip}:{port}")
            self._send(notif_id, -1, -errno.EACCES)
        except Exception:
            try:
                self._send(notif_id, -1, -errno.EACCES)
            except Exception:
                pass

    def _handle(self, notif_id: int, pid: int, args: tuple, nr: int = -1) -> None:
        # The datagram surface, by the architecture's own numbers: socket()
        # is the command's own creation; sendto with a destination is
        # judged by the UDP ownership rule. bind() is not notified (the
        # filter allows it in-kernel; the escape is judged at sendto).
        # sendmsg/sendmmsg never arrive: the filter refuses them
        # in-kernel (the destination sits inside the msghdr, invisible
        # to BPF).
        (audit, _sc, nr_connect, nr_socket, nr_uring, _po, _pg) = _arch()
        (_nr_bind, nr_sendto, _nsm, _nsmm) = _arch_dgram()
        if nr == nr_socket:
            self._datagram_socket(notif_id, pid, args, nr)
            return
        if nr == nr_sendto:
            self._datagram_send(notif_id, pid, args)
            return
        err = errno.EACCES
        val = -1
        try:
            fd, addr_ptr, addr_len = int(args[0]), int(args[1]), int(args[2]) & 0xFFFFFFFF
            if addr_len < 2 or addr_len > 256:
                self._send(notif_id, -1, -errno.EINVAL)
                return
            with open(f"/proc/{pid}/mem", "rb", buffering=0) as mem:
                mem.seek(addr_ptr)
                raw = mem.read(addr_len)
            if len(raw) != addr_len or not self._valid(notif_id):
                self._send(notif_id, -1, -errno.EACCES)
                return
            family = struct.unpack_from("<H", raw)[0]
            if family in (socket.AF_INET, socket.AF_INET6):
                if family == socket.AF_INET and len(raw) >= 8:
                    port = struct.unpack_from(">H", raw, 2)[0]
                    ip = ipaddress.ip_address(raw[4:8])
                elif family == socket.AF_INET6 and len(raw) >= 24:
                    port = struct.unpack_from(">H", raw, 2)[0]
                    ip = ipaddress.ip_address(raw[8:24])
                else:
                    self._send(notif_id, -1, -errno.EINVAL)
                    return
                why = self.inet_refusal(ip, port)
                if why:
                    self.refused.append(why)
                    self._send(notif_id, -1, -errno.EACCES)
                    return
            elif family not in (socket.AF_UNIX, socket.AF_UNSPEC) and self.net_mode != "open":
                self.refused.append(f"address family {family}")
                self._send(notif_id, -1, -errno.EACCES)
                return
            if family == socket.AF_UNIX:
                name = raw[2:]
                if not name or name[0] == 0:
                    self.refused.append("abstract unix socket")
                    self._send(notif_id, -1, -errno.EACCES)
                    return
                path = name.split(b"\0", 1)[0].decode("utf-8", "surrogateescape")
                if not os.path.isabs(path):
                    path = os.path.join(os.readlink(f"/proc/{pid}/cwd"), path)
                real = os.path.realpath(path)
                if not self.unix_path_allowed(real):
                    self.refused.append(real)
                    self._send(notif_id, -1, -errno.EACCES)
                    return
                encoded = real.encode("utf-8", "surrogateescape")
                if len(encoded) >= 108:
                    self._send(notif_id, -1, -errno.ENAMETOOLONG)
                    return
                raw = struct.pack("<H", socket.AF_UNIX) + encoded + b"\0"
            dup = self._dup_fd(pid, fd)
            if dup < 0:
                self._send(notif_id, -1, -errno.EBADF)
                return
            try:
                rc = _libc.connect(dup, raw, len(raw))
                if rc == 0:
                    val, err = 0, 0
                else:
                    err = ctypes.get_errno()
            finally:
                os.close(dup)
            self._send(notif_id, val, -err if err else 0)
        except Exception:
            try:
                self._send(notif_id, -1, -errno.EACCES)
            except Exception:
                pass

    def serve(self, child_pid: int) -> int:
        """Answer notifications until ``child_pid`` exits; return its status."""
        self.root_pid = int(child_pid)
        poller = select.poll()
        poller.register(self.listener, select.POLLIN)
        while True:
            try:
                wpid, status = os.waitpid(child_pid, os.WNOHANG)
            except ChildProcessError:
                return 0
            if wpid == child_pid:
                return status
            for _fd, ev in poller.poll(100):
                if ev & (select.POLLHUP | select.POLLERR | select.POLLNVAL):
                    continue
                buf = ctypes.create_string_buffer(80)
                if _libc.ioctl(self.listener, ctypes.c_ulong(_NOTIF_RECV), buf) != 0:
                    continue
                notif_id, pid, nr, args = parse_notification(buf.raw)
                threading.Thread(target=self._handle,
                                 args=(notif_id, pid, args, nr),
                                 daemon=True).start()
