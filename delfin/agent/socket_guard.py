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


def _program() -> list[tuple[int, int, int, int]]:
    audit, _sc, nr_connect, nr_socket, nr_uring, _po, _pg = _arch()
    LD, JEQ, JSET, RET, AND = 0x20, 0x15, 0x45, 0x06, 0x54
    deny = _RET_ERRNO | errno.EPERM
    prog: list[tuple[int, int, int, int]] = []
    # 0 arch must match, or nothing runs (32-bit entry points bypass numbers)
    prog.append((LD, 0, 0, 4))
    prog.append((JEQ, 1, 0, audit))
    prog.append((RET, 0, 0, deny))
    prog.append((LD, 0, 0, 0))
    if audit == 0xC000003E:
        # x32 system calls carry bit 30; refuse them all.
        prog.append((JSET, 0, 1, 0x40000000))
        prog.append((RET, 0, 0, deny))
    prog.append((JEQ, 0, 1, nr_connect))
    prog.append((RET, 0, 0, _RET_USER_NOTIF))
    prog.append((JEQ, 0, 1, nr_uring))
    prog.append((RET, 0, 0, _RET_ERRNO | errno.ENOSYS))
    prog.append((JEQ, 0, 7, nr_socket))         # not socket(): ALLOW below
    prog.append((LD, 0, 0, 16))                 # args[0] low: domain
    prog.append((JEQ, 0, 5, socket.AF_UNIX))    # not AF_UNIX: ALLOW below
    prog.append((LD, 0, 0, 24))                 # args[1] low: type
    prog.append((AND, 0, 0, 0xF))
    prog.append((JEQ, 1, 0, socket.SOCK_DGRAM))
    prog.append((RET, 0, 0, _RET_ALLOW))
    prog.append((RET, 0, 0, _RET_ERRNO | errno.EACCES))
    prog.append((RET, 0, 0, _RET_ALLOW))
    return prog


def _block_unix_program() -> list[tuple[int, int, int, int]]:
    """Where the kernel cannot hand connects to a supervisor: no Unix socket
    at all (and no io_uring), everything else as usual."""
    audit, _sc, _nc, nr_socket, nr_uring, _po, _pg = _arch()
    LD, JEQ, JSET, RET = 0x20, 0x15, 0x45, 0x06
    deny = _RET_ERRNO | errno.EPERM
    prog = [(LD, 0, 0, 4), (JEQ, 1, 0, audit), (RET, 0, 0, deny), (LD, 0, 0, 0)]
    if audit == 0xC000003E:
        prog += [(JSET, 0, 1, 0x40000000), (RET, 0, 0, deny)]
    prog += [(JEQ, 0, 1, nr_uring), (RET, 0, 0, _RET_ERRNO | errno.ENOSYS),
             (JEQ, 0, 3, nr_socket), (LD, 0, 0, 16),
             (JEQ, 0, 1, socket.AF_UNIX), (RET, 0, 0, _RET_ERRNO | errno.EACCES),
             (RET, 0, 0, _RET_ALLOW)]
    return prog


def filter_supported() -> bool:
    return (_arch() is not None and sys.platform.startswith("linux")
            and sys.byteorder == "little")


def install_block_unix_filter() -> None:
    a = _arch()
    insns = _block_unix_program()
    arr = (_SockFilter * len(insns))(*[_SockFilter(*i) for i in insns])
    prog = _SockFprog(len(insns), arr)
    rc = _libc.syscall(ctypes.c_long(a[1]), ctypes.c_uint(_SECCOMP_SET_MODE_FILTER),
                       ctypes.c_uint(0), ctypes.byref(prog))
    if rc != 0:
        raise OSError(ctypes.get_errno(), "seccomp filter: " + os.strerror(ctypes.get_errno()))


def install_filter() -> int:
    """In the process that will run the command: install the filter and
    return the listener descriptor. Needs no_new_privs set already."""
    a = _arch()
    insns = _program()
    arr = (_SockFilter * len(insns))(*[_SockFilter(*i) for i in insns])
    prog = _SockFprog(len(insns), arr)
    fd = _libc.syscall(ctypes.c_long(a[1]), ctypes.c_uint(_SECCOMP_SET_MODE_FILTER),
                       ctypes.c_uint(_SECCOMP_FILTER_FLAG_NEW_LISTENER),
                       ctypes.byref(prog))
    if fd < 0:
        raise OSError(ctypes.get_errno(), "seccomp filter: " + os.strerror(ctypes.get_errno()))
    return int(fd)


class Supervisor:
    """Answers the command's connect() calls until the command ends."""

    def __init__(self, listener: int, allowed: list[str]):
        self.listener = listener
        self.allowed = [os.path.realpath(p).rstrip("/") for p in allowed if p]
        a = _arch()
        self._nr_pidfd_open, self._nr_pidfd_getfd = a[5], a[6]
        self.refused: list[str] = []

    # -- policy -------------------------------------------------------------
    def unix_path_allowed(self, path: str) -> bool:
        real = os.path.realpath(path)
        return any(real == p or real.startswith(p + "/") for p in self.allowed)

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

    def _dup_fd(self, pid: int, fd: int) -> int:
        pidfd = _libc.syscall(ctypes.c_long(self._nr_pidfd_open), ctypes.c_int(pid),
                              ctypes.c_uint(0))
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
