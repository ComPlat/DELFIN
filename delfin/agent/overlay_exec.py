"""Mount the workspace overlay, then become the cage.

Run BY PATH with ``python -I`` in a fresh single-threaded process, the way
``landlock_exec`` is: no DELFIN import, so nothing in the workspace or the
environment changes what runs before the overlay holds.

    python -I -S overlay_exec.py --status FILE \
        --layer ROOT UPPER WORK [--layer ...] -- bwrap ARGS...

Steps, each a refusal (exit ``SETUP_FAILED``, one line on stderr starting
with ``MARK``) if it fails -- the command is never run without the
overlay it was promised. FILE (opened before any mount, close-on-exec)
receives ``ok`` once every mount holds: the caller reads it, not the exit
code or stderr, which the command itself could imitate.

1. a user and a mount namespace of its own (``unshare``), mapping the
   caller's uid and gid to 0 inside;
2. for each ``--layer``: an overlay over ROOT itself
   (``lowerdir=ROOT,upperdir=UPPER,workdir=WORK``), so every write the
   command makes under ROOT lands in UPPER and ROOT stays as it was;
3. the root's real trash (``ROOT/.delfin-trash``, opened before the
   overlay covered it) bound read-only over the merged view: the command
   reads the trash, cannot change it, and never looks up UPPER or WORK
   through the overlay (both lie inside the trash; the kernel answers a
   lookup of a layer through its own overlay with ELOOP); the trash's
   ``.layers`` directory is covered by an empty read-only tmpfs, so a
   ``find`` over the root meets neither the layers nor the kernel's
   mode-000 work entries;
4. ``exec`` of the bubblewrap argv given after ``--``, with
   ``--unshare-user --uid UID --gid GID`` added, so the command runs under
   the caller's own ids again.

``--probe DIR`` mounts a throwaway overlay in DIR and exits 0 when it
worked.
"""
from __future__ import annotations

import ctypes
import os
import sys

SETUP_FAILED = 97
MARK = "delfin-overlay:"
TRASH = ".delfin-trash"
LAYERS = ".layers"

_CLONE_NEWNS = 0x00020000
_CLONE_NEWUSER = 0x10000000
_MS_RDONLY = 1
_MS_NOSUID = 2
_MS_NODEV = 4
_MS_REMOUNT = 32
_MS_BIND = 4096
_MS_REC = 16384
_MS_PRIVATE = 1 << 18

_libc = ctypes.CDLL(None, use_errno=True)
_libc.mount.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
                        ctypes.c_ulong, ctypes.c_char_p]
_libc.unshare.argtypes = [ctypes.c_int]


def _fail(why: str) -> None:
    sys.stderr.write(f"{MARK} {why}\n")
    sys.stderr.flush()
    os._exit(SETUP_FAILED)


def _mount(src, target, fstype, flags, data=None) -> None:
    rc = _libc.mount(None if src is None else os.fsencode(src),
                     os.fsencode(target),
                     None if fstype is None else fstype.encode(),
                     flags, None if data is None else os.fsencode(data))
    if rc != 0:
        err = ctypes.get_errno()
        raise OSError(err, f"mount {fstype or 'bind'} on {target}: "
                           f"{os.strerror(err)}")


def _write(path: str, text: str) -> None:
    fd = os.open(path, os.O_WRONLY)
    try:
        os.write(fd, text.encode())
    finally:
        os.close(fd)


def _enter_namespace() -> tuple[int, int]:
    uid, gid = os.getuid(), os.getgid()
    if _libc.unshare(_CLONE_NEWUSER | _CLONE_NEWNS) != 0:
        err = ctypes.get_errno()
        raise OSError(err, f"unshare: {os.strerror(err)}")
    _write("/proc/self/setgroups", "deny")
    _write("/proc/self/uid_map", f"0 {uid} 1")
    _write("/proc/self/gid_map", f"0 {gid} 1")
    _mount(None, "/", None, _MS_REC | _MS_PRIVATE)
    return uid, gid


#: Characters the overlay option string cannot carry in a path.
UNSAFE = (",", ":", "\\", "\n")


def _overlay(root: str, upper: str, work: str) -> None:
    for p in (root, upper, work):
        if any(c in p for c in UNSAFE):
            raise OSError(0, f"path not expressible as an overlay option: {p!r}")
    trash = os.path.join(root, TRASH)
    fd = os.open(trash, os.O_PATH | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        _mount("overlay", root, "overlay", 0,
               f"lowerdir={root},upperdir={upper},workdir={work}")
        _mount(f"/proc/self/fd/{fd}", trash, None, _MS_BIND | _MS_REC)
        layers = os.path.join(trash, LAYERS)
        if os.path.isdir(layers):
            _mount("tmpfs", layers, "tmpfs",
                   _MS_RDONLY | _MS_NOSUID | _MS_NODEV, "size=4k,mode=0555")
        _mount(None, trash, None,
               _MS_BIND | _MS_REMOUNT | _MS_RDONLY | _MS_NOSUID | _MS_NODEV)
    finally:
        os.close(fd)


def _probe(base: str) -> int:
    dirs = {k: os.path.join(base, k) for k in ("lower", "upper", "work")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join(dirs["lower"], TRASH), exist_ok=True)
    try:
        _enter_namespace()
        _overlay(dirs["lower"], dirs["upper"], dirs["work"])
        with open(os.path.join(dirs["lower"], "probe"), "w") as fh:
            fh.write("1")
    except OSError as exc:
        _fail(str(exc))
    return 0


def main(argv: list[str]) -> int:
    if len(argv) == 2 and argv[0] == "--probe":
        return _probe(argv[1])
    layers: list[tuple[str, str, str]] = []
    status = ""
    i = 0
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--layer" and i + 3 < len(argv):
            layers.append((argv[i + 1], argv[i + 2], argv[i + 3]))
            i += 4
            continue
        if argv[i] == "--status" and i + 1 < len(argv):
            status = argv[i + 1]
            i += 2
            continue
        _fail(f"unknown argument {argv[i]!r}")
    cmd = argv[i + 1:]
    if not layers or not status or not cmd \
            or os.path.basename(cmd[0]) != "bwrap":
        _fail("needs --status, at least one --layer and a bwrap argv "
              "after --")
    try:
        sfd = os.open(status, os.O_WRONLY | os.O_CREAT | os.O_TRUNC
                      | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600)
    except OSError as exc:
        _fail(f"status file: {exc}")
    try:
        uid, gid = _enter_namespace()
        for root, upper, work in layers:
            _overlay(root, upper, work)
    except OSError as exc:
        _fail(str(exc))
    cmd = [cmd[0], "--unshare-user", "--uid", str(uid), "--gid", str(gid),
           *cmd[1:]]
    os.write(sfd, b"ok")
    try:
        os.execvp(cmd[0], cmd)
    except OSError as exc:
        os.ftruncate(sfd, 0)
        _fail(f"exec {cmd[0]}: {exc}")
    return SETUP_FAILED


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
