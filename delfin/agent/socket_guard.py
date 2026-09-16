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
    # machine: (AUDIT_ARCH, seccomp, connect, socket, io_uring_setup, pidfd_open, pidfd_getfd)
    "x86_64": (0xC000003E, 317, 42, 41, 425, 434, 438),
    "aarch64": (0xC00000B7, 277, 203, 198, 425, 434, 438),
}

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

_libc = ctypes.CDLL(None, use_errno=True)
_libc.syscall.restype = ctypes.c_long
_libc.ioctl.restype = ctypes.c_int
_libc.connect.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_uint32]
_libc.connect.restype = ctypes.c_int

#: Beneath these a Unix socket may be reached from any sandbox: the cluster
#: authentication (munge, for sbatch/squeue) and the name services (sssd,
#: nscd, systemd-resolved), without which user names do not resolve.
DEFAULT_ALLOWED = (
    "/run/munge", "/var/run/munge", "/var/lib/sss/pipes", "/run/nscd",
    "/var/run/nscd", "/run/systemd/resolve",
)


#: Metadata services that hand out cloud credentials; link-local ranges are
#: refused as a whole besides these.
_METADATA = {ipaddress.ip_address(a) for a in (
    "169.254.169.254", "100.100.100.200", "168.63.129.16", "fd00:ec2::254")}


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
    and with a restricted network no non-stream inet socket (UDP would
    leave without connect(), past the supervisor). A foreign architecture
    (32-bit entry points use other numbers) runs nothing."""
    audit, _sc, nr_connect, nr_socket, nr_uring, _po, _pg = _arch()
    items, tail = _head(audit)
    items += [
        (_JEQ, nr_connect, "notify", None),
        (_JEQ, nr_uring, "enosys", None),
        (_JEQ, nr_socket, None, "allow"),
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
        items += [("label", "inet"), (_LD, 24), (_AND, 0xF),
                  (_JEQ, socket.SOCK_STREAM, "allow", "eacces")]
    items += [("label", "allow"), (_RET, _RET_ALLOW),
              ("label", "eacces"), (_RET, _RET_ERRNO | errno.EACCES),
              ("label", "enosys"), (_RET, _RET_ERRNO | errno.ENOSYS),
              ("label", "notify"), (_RET, _RET_USER_NOTIF)] + tail
    return _assemble(items)


def _block_unix_program(block_inet_datagrams: bool = False):
    """Where the kernel cannot hand connects to a supervisor: no Unix socket
    at all (and no io_uring), and with a restricted network no inet socket
    at all either."""
    audit, _sc, _nc, nr_socket, nr_uring, _po, _pg = _arch()
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


class Supervisor:
    """Answers the command's connect() calls until the command ends."""

    def __init__(self, listener: int, allowed: list[str], net_mode: str = "open",
                 allowed_tcp=(), denied=()):
        self.listener = listener
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
        if self.net_mode == "proxy" and (ip, port) in self.allowed_tcp:
            return ""
        return (f"direct connection to {ip}:{port}; the sandbox reaches the "
                "network through its proxy only")

    def unix_path_allowed(self, path: str) -> bool:
        real = os.path.realpath(path)

        def beneath(prefixes):
            return any(real == p or real.startswith(p + "/") or p == ""
                       for p in prefixes)

        if self.denied and beneath(self.denied):
            return False
        return beneath(self.allowed)

    # -- plumbing -----------------------------------------------------------
    def _send(self, notif_id: int, val: int, err: int) -> None:
        resp = struct.pack("<QqiI", notif_id, val, err, 0)
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

    def _handle(self, notif_id: int, pid: int, args: tuple) -> None:
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
                notif_id, pid, _flags = struct.unpack_from("<QII", buf.raw)
                args = struct.unpack_from("<6Q", buf.raw, 32)
                threading.Thread(target=self._handle, args=(notif_id, pid, args),
                                 daemon=True).start()
