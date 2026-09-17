"""Run a command under a Landlock filesystem policy, then become it.

Filesystem isolation for the agent's shell uses bubblewrap where
unprivileged user namespaces work. Many hosts refuse those (HPC login
nodes, hardened distributions, CI runners) while their kernel has
Landlock (Linux 5.13+), which any process may apply to itself. This
helper is the second way to the same promise:

* writes only beneath the workspace roots, a private temp directory and
  ``/dev``;
* reads everywhere except the credential locations it is told to hide.

It is executed BY PATH with ``python -I`` in a fresh single-threaded
process, applies the policy to itself and ``exec``s the command, which
inherits the policy and cannot drop it. Standalone on purpose: no DELFIN
import, so nothing in the workspace or the environment can change what
runs before the policy holds. Any failure refuses the command (exit 126)
instead of running it unconfined.

    python -I landlock_exec.py --write DIR ... --hide PATH ... \
        [--tmpdir DIR] [--allow-socket DIR ...] [--strict 1] [--fs 0] \
        [--net open|proxy|none] [--proxy-socket PATH] [--allow-tcp IP:PORT ...] \
        [--deny-socket DIR ...] \
        -- /bin/bash -c CMD

Unix sockets are guarded too (``socket_guard``): the command reaches them
only beneath the write roots, the private temp directory and the cluster's
authentication and name services.
    python -I landlock_exec.py --probe
"""
from __future__ import annotations

import ctypes
import os
import stat
import sys

_SYS_CREATE_RULESET = 444
_SYS_ADD_RULE = 445
_SYS_RESTRICT_SELF = 446
_CREATE_RULESET_VERSION = 1
_RULE_PATH_BENEATH = 1
_PR_SET_NO_NEW_PRIVS = 38

_EXECUTE = 1 << 0
_WRITE_FILE = 1 << 1
_READ_FILE = 1 << 2
_READ_DIR = 1 << 3
_REMOVE_DIR = 1 << 4
_REMOVE_FILE = 1 << 5
_MAKE_CHAR = 1 << 6
_MAKE_DIR = 1 << 7
_MAKE_REG = 1 << 8
_MAKE_SOCK = 1 << 9
_MAKE_FIFO = 1 << 10
_MAKE_BLOCK = 1 << 11
_MAKE_SYM = 1 << 12
_REFER = 1 << 13
_TRUNCATE = 1 << 14

_READ = _READ_FILE | _READ_DIR
_FILE_RIGHTS = _EXECUTE | _WRITE_FILE | _READ_FILE | _TRUNCATE


class _RulesetAttr(ctypes.Structure):
    _fields_ = [("handled_access_fs", ctypes.c_uint64)]


class _PathBeneathAttr(ctypes.Structure):
    _pack_ = 1
    _fields_ = [("allowed_access", ctypes.c_uint64),
                ("parent_fd", ctypes.c_int32)]


_libc = ctypes.CDLL(None, use_errno=True)
_libc.syscall.restype = ctypes.c_long


def abi_version() -> int:
    """The kernel's Landlock ABI version, or 0 where there is none."""
    try:
        v = _libc.syscall(ctypes.c_long(_SYS_CREATE_RULESET), None,
                          ctypes.c_size_t(0),
                          ctypes.c_uint32(_CREATE_RULESET_VERSION))
    except Exception:
        return 0
    return int(v) if v > 0 else 0


def _write_rights(abi: int) -> int:
    rights = (_WRITE_FILE | _REMOVE_DIR | _REMOVE_FILE | _MAKE_CHAR
              | _MAKE_DIR | _MAKE_REG | _MAKE_SOCK | _MAKE_FIFO
              | _MAKE_BLOCK | _MAKE_SYM)
    if abi >= 2:
        rights |= _REFER
    if abi >= 3:
        rights |= _TRUNCATE
    return rights


class _Refused(Exception):
    pass


def _add_rule(ruleset_fd: int, path: str, access: int) -> None:
    try:
        fd = os.open(path, os.O_PATH | os.O_CLOEXEC)
    except OSError:
        return                      # gone or unreadable: nothing to grant
    try:
        if not stat.S_ISDIR(os.fstat(fd).st_mode):
            access &= _FILE_RIGHTS
        if not access:
            return
        attr = _PathBeneathAttr(allowed_access=access, parent_fd=fd)
        rc = _libc.syscall(ctypes.c_long(_SYS_ADD_RULE),
                           ctypes.c_int(ruleset_fd),
                           ctypes.c_int(_RULE_PATH_BENEATH),
                           ctypes.byref(attr), ctypes.c_uint32(0))
        if rc != 0:
            raise _Refused(f"landlock_add_rule({path}) failed: "
                           f"{os.strerror(ctypes.get_errno())}")
    finally:
        os.close(fd)


def _grant_except(ruleset_fd: int, directory: str, hidden: list[str],
                  access: int, *, list_dirs: bool) -> None:
    """Grant ``access`` beneath ``directory`` except beneath ``hidden``.

    Landlock only allows, and a rule covers everything beneath its path,
    so a directory that holds a hidden path is granted entry by entry and
    only the branch towards the hidden path is walked. With ``list_dirs``
    the directories on that branch stay listable (``ls ~`` works); the
    hidden content still needs a right it does not get. Symlinks are not
    granted: access through one is judged at its target."""
    prefix = directory.rstrip("/") + "/"
    if not any(h == directory or h.startswith(prefix) for h in hidden):
        _add_rule(ruleset_fd, directory, access)
        return
    if directory in hidden:
        return
    if list_dirs:
        _add_rule(ruleset_fd, directory, _READ_DIR)
    try:
        entries = list(os.scandir(directory))
    except OSError:
        return
    for entry in entries:
        full = prefix + entry.name
        if full in hidden:
            continue
        try:
            if entry.is_symlink():
                continue
            if entry.is_dir(follow_symlinks=False):
                _grant_except(ruleset_fd, full, hidden, access,
                              list_dirs=list_dirs)
            else:
                _add_rule(ruleset_fd, full, access & _FILE_RIGHTS)
        except OSError:
            continue


def apply(write: list[str], hide: list[str], read_only: list[str] = ()) -> None:
    """Restrict this process. Raises ``_Refused`` when it cannot.

    ``read_only`` paths stay readable (from the top-level read grant) but are
    carved out of every write grant, so a directory inside a write root can
    be readable and executable without being writable -- git's ``hooks``."""
    abi = abi_version()
    if abi < 1:
        raise _Refused("Landlock is not available on this kernel")
    write_rights = _write_rights(abi)
    handled = _READ | write_rights
    attr = _RulesetAttr(handled_access_fs=handled)
    ruleset_fd = _libc.syscall(ctypes.c_long(_SYS_CREATE_RULESET),
                               ctypes.byref(attr),
                               ctypes.c_size_t(ctypes.sizeof(attr)),
                               ctypes.c_uint32(0))
    if ruleset_fd < 0:
        raise _Refused("landlock_create_ruleset failed: "
                       f"{os.strerror(ctypes.get_errno())}")
    try:
        hidden = sorted({os.path.realpath(h) for h in hide if h})
        no_write = sorted(set(hidden) | {os.path.realpath(r) for r in read_only if r})
        _grant_except(ruleset_fd, "/", hidden, _READ, list_dirs=True)
        # Write roots get the write rights only; reading comes from the
        # walk above, so a hidden path inside a workspace stays hidden,
        # and it is not writable either. A read_only path is skipped by the
        # write walk but not by the read grant: readable, not writable.
        for path in write:
            if path:
                _grant_except(ruleset_fd, os.path.realpath(path), no_write,
                              write_rights, list_dirs=False)
        if _libc.prctl(_PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0) != 0:
            raise _Refused("prctl(PR_SET_NO_NEW_PRIVS) failed")
        rc = _libc.syscall(ctypes.c_long(_SYS_RESTRICT_SELF),
                           ctypes.c_int(ruleset_fd), ctypes.c_uint32(0))
        if rc != 0:
            raise _Refused("landlock_restrict_self failed: "
                           f"{os.strerror(ctypes.get_errno())}")
    finally:
        os.close(ruleset_fd)


def _load_socket_guard():
    """socket_guard.py beside this file, loaded by path: under ``python -I``
    the script's directory is not on sys.path, which is the point."""
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "socket_guard.py")
    spec = importlib.util.spec_from_file_location("_delfin_socket_guard", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _exit_code(status: int) -> int:
    code = os.waitstatus_to_exitcode(status)
    return 128 - code if code < 0 else code


class _Forwarder:
    """A loopback port inside the sandbox that leads to the egress proxy's
    Unix socket. Runs in the supervising process, outside the filter."""

    def __init__(self, unix_path: str):
        import socket as _socket
        self._socket = _socket
        self.unix_path = unix_path
        self.server = _socket.socket(_socket.AF_INET, _socket.SOCK_STREAM)
        self.server.bind(("127.0.0.1", 0))
        self.server.listen(64)
        self.port = self.server.getsockname()[1]

    def start(self) -> None:
        import threading
        threading.Thread(target=self._accept, daemon=True).start()

    def _accept(self) -> None:
        import threading
        while True:
            try:
                conn, _ = self.server.accept()
            except OSError:
                return
            threading.Thread(target=self._pipe, args=(conn,), daemon=True).start()

    def _pipe(self, conn) -> None:
        import select as _select
        up = self._socket.socket(self._socket.AF_UNIX, self._socket.SOCK_STREAM)
        try:
            up.connect(self.unix_path)
            pair = [conn, up]
            while True:
                readable, _, _ = _select.select(pair, [], [], 600)
                if not readable:
                    return
                for s in readable:
                    data = s.recv(65536)
                    if not data:
                        return
                    (up if s is conn else conn).sendall(data)
        except OSError:
            return
        finally:
            for s in (conn, up):
                try:
                    s.close()
                except OSError:
                    pass


_PROXY_VARS = ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY")
_UNSET_VARS = ("no_proxy", "NO_PROXY", "all_proxy", "ALL_PROXY",
               "SSH_AUTH_SOCK", "SSH_AGENT_PID", "TMUX", "TMUX_PANE", "STY",
               "DBUS_SESSION_BUS_ADDRESS")


def _run_guarded(write: list[str], hide: list[str], allowed: list[str],
                 command: list[str], strict: bool, *, use_fs: bool = True,
                 net_mode: str = "open", proxy_socket: str = "",
                 allow_tcp: list[str] | None = None,
                 deny_socket: list[str] | None = None,
                 read_only: list[str] | None = None) -> int:
    """Apply the filesystem policy (unless ``use_fs`` is off, inside
    bubblewrap), the socket guard and the network mode, and run the command.

    With user notification the command runs as a child under a filter that
    hands every connect() to the supervisor in this process. Without it,
    Unix sockets are refused outright and, with a restricted network, inet
    sockets too; where no filter can be installed at all a strict run
    refuses the command."""
    allow_tcp = list(allow_tcp or [])
    restricted = net_mode in ("proxy", "none")
    try:
        sg = _load_socket_guard()
    except Exception:
        sg = None
    guard = sg is not None and sg.available()
    forwarder = None
    env_updates: dict[str, str] = {}
    if net_mode == "proxy" and proxy_socket and guard:
        try:
            token_path = os.path.join(os.path.dirname(proxy_socket), "token")
            with open(token_path, encoding="utf-8") as fh:
                token = fh.read().strip()
            forwarder = _Forwarder(proxy_socket)
            url = f"http://delfin:{token}@127.0.0.1:{forwarder.port}"
            env_updates = {k: url for k in _PROXY_VARS}
            allow_tcp.append(f"127.0.0.1:{forwarder.port}")
        except OSError as exc:
            if strict:
                print(f"refused: the sandbox's network proxy is not reachable "
                      f"({exc}); the command was not run.", file=sys.stderr)
                return 126
            forwarder = None
    if not guard:
        if restricted and strict:
            print("refused: this host cannot restrict the command's network "
                  "(no seccomp user notification); the command was not run.",
                  file=sys.stderr)
            return 126
        try:
            if use_fs:
                apply(write, hide, read_only or [])
            elif _libc.prctl(_PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0) != 0:
                raise _Refused("prctl(PR_SET_NO_NEW_PRIVS) failed")
        except _Refused as exc:
            print(f"refused: filesystem isolation could not be applied ({exc}); "
                  "the command was not run.", file=sys.stderr)
            return 126
        try:
            if sg is None or not sg.filter_supported():
                raise OSError("no seccomp filter for this platform")
            sg.install_block_unix_filter(block_inet_datagrams=False)
        except Exception as exc:
            if strict:
                print("refused: the command could not be kept from the user's "
                      f"sessions ({exc}); the command was not run.", file=sys.stderr)
                return 126
        for key in _UNSET_VARS:
            os.environ.pop(key, None)
        try:
            os.execvp(command[0], command)
        except OSError as exc:
            print(f"landlock_exec: cannot run {command[0]}: {exc}", file=sys.stderr)
            return 127
        return 127

    import socket as _socket
    parent_end, child_end = _socket.socketpair()
    pid = os.fork()
    if pid == 0:                                    # the command
        parent_end.close()
        try:
            if forwarder is not None:
                forwarder.server.close()
            if use_fs:
                apply(write, hide, read_only or [])
            elif _libc.prctl(_PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0) != 0:
                raise _Refused("prctl(PR_SET_NO_NEW_PRIVS) failed")
            listener = sg.install_filter(block_inet_datagrams=restricted)
            _socket.send_fds(child_end, [b"L"], [listener])
            os.close(listener)                      # the command must not hold it
            child_end.recv(1)
            child_end.close()
            for key in _UNSET_VARS:
                os.environ.pop(key, None)
            os.environ.update(env_updates)
            os.execvp(command[0], command)
        except BaseException as exc:                # noqa: BLE001 - child must exit
            try:
                os.write(2, f"refused: sandbox could not be applied ({exc}); "
                            "the command was not run.\n".encode())
            finally:
                os._exit(126)
    child_end.close()
    try:
        _msg, fds, _flags, _addr = _socket.recv_fds(parent_end, 16, 1)
    except OSError:
        fds = []
    if not fds:
        _pid, status = os.waitpid(pid, 0)
        return _exit_code(status)
    listener = fds[0]
    parent_end.send(b"k")
    parent_end.close()
    try:
        _libc.prctl(4, 0, 0, 0, 0)                  # PR_SET_DUMPABLE: not readable
    except Exception:
        pass
    if forwarder is not None:
        forwarder.start()
    status = sg.Supervisor(listener, allowed, net_mode, allow_tcp,
                           denied=deny_socket or ()).serve(pid)
    return _exit_code(status)


def main(argv: list[str]) -> int:
    if argv[:1] == ["--probe"]:
        print(abi_version())
        return 0
    write: list[str] = []
    hide: list[str] = []
    allow_socket: list[str] = []
    allow_tcp: list[str] = []
    deny_socket: list[str] = []
    read_only: list[str] = []
    strict = False
    use_fs = True
    net_mode = "open"
    proxy_socket = ""
    tmpdir = ""
    i = 0
    while i < len(argv) and argv[i] != "--":
        flag = argv[i]
        value = argv[i + 1] if i + 1 < len(argv) else ""
        if flag == "--write":
            write.append(value)
        elif flag == "--hide":
            hide.append(value)
        elif flag == "--tmpdir":
            tmpdir = value
        elif flag == "--allow-socket":
            allow_socket.append(value)
        elif flag == "--allow-tcp":
            allow_tcp.append(value)
        elif flag == "--deny-socket":
            deny_socket.append(value)
        elif flag == "--read-only":
            read_only.append(value)
        elif flag == "--strict":
            strict = value == "1"
        elif flag == "--fs":
            use_fs = value != "0"
        elif flag == "--net":
            net_mode = value if value in ("open", "proxy", "none") else "none"
        elif flag == "--proxy-socket":
            proxy_socket = value
        else:
            print(f"landlock_exec: unknown option {flag}", file=sys.stderr)
            return 126
        i += 2
    command = argv[i + 1:]
    if not command:
        print("landlock_exec: no command", file=sys.stderr)
        return 126
    if tmpdir:
        write.append(tmpdir)
        for key in ("TMPDIR", "TMP", "TEMP"):
            os.environ[key] = tmpdir
    roots = [w for w in write if w]
    write.append("/dev")
    try:
        sg_defaults = list(_load_socket_guard().DEFAULT_ALLOWED)
    except Exception:
        sg_defaults = []
    return _run_guarded(write, hide, sg_defaults + roots + allow_socket,
                        command, strict, use_fs=use_fs, net_mode=net_mode,
                        proxy_socket=proxy_socket, allow_tcp=allow_tcp,
                        deny_socket=deny_socket, read_only=read_only)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
