"""Every file change a bypass shell command makes, reversible, whatever its size.

Input: the session's permissions and one foreground ``bash`` command in
``bypassPermissions`` whose cage is bubblewrap. Output: the command runs
over an overlay of the workspace roots; afterwards :func:`apply` moves
what it wrote into the real roots and returns one :class:`Applied` per
changed path, which the executor journals (``change_journal``,
``overlay: True`` records) so ``undo_changes`` restores it.

Mechanism, per command:

1. :func:`prepare` creates, for each eligible root R,
   ``R/.delfin-trash/.layers/<id>/{u,w}`` (upper and work directory, on
   R's own filesystem) and records the start time.
2. :func:`wrap` puts ``overlay_exec`` in front of the bubblewrap argv: in
   a user and mount namespace of its own it mounts
   ``overlay(lowerdir=R, upperdir=u, workdir=w)`` over R itself, binds the
   real trash read-only over the merged view (the layers hidden beneath an
   empty mount), and execs bubblewrap, whose binds of R now resolve
   through the overlay. R is not written while the command runs: every
   create, change, rename and delete lands in u.
3. :func:`apply`, once the command has ended (and with it the namespace
   and the mount), walks u. A new path is renamed into R. A changed or
   deleted one (whiteout, opaque directory) has its original renamed into
   ``R/.delfin-trash/<UTC date>/cmd-<id>/<relative path>`` first, then the
   new version renamed in. Every move is a rename on one filesystem: no
   byte of a file is copied, a 50 MiB file costs what a 5-byte one does.
   A whole new or deleted directory is one rename and one journal record.
4. the layer directory is removed (its work dir holds mode-000 entries
   the kernel created; they are made writable first).

Decisions, each with the reason:

* Apply after every command, not at the turn's end: the person's
  dashboard and the next command see the change at once; a file tool
  writing the same path later in the turn would otherwise be hidden under
  a stale upper copy; an undo has one command's changes to take back.
* A path the person (or any other process) changed while the command ran
  (its ctime is after the start) is not overwritten: the command's version
  is placed beside it as ``<name>.delfin-conflict-<id>`` and the note says
  so; a deletion of such a path is not applied. Both versions are kept.
* Background commands (``bash_background``) do not run over the overlay:
  their output must appear while they run (a calculation's log, a
  watcher's file), and an overlay shows nothing until it is applied. They
  keep the shell change journal and the trash refusal (``unrestorable``).
* ``.git`` inside a root is in the overlay like everything else. Created
  files under ``.git/objects`` are applied without a journal record
  (content-addressed and append-only: a creation there replaces nothing);
  refs, the index, packed-refs, HEAD and every deletion (gc) go through
  the trash and the journal like any file, so a reset or a branch delete
  is undone with the files. A git directory outside the roots (a
  worktree's common dir) is not overlaid.
* Disk: originals move into the trash instead of being freed, kept
  ``TRASH_DAYS`` days. A root whose filesystem has less than
  ``MIN_FREE_BYTES`` free (or ``MIN_FREE_FRACTION``) is not overlaid for
  that command and the note says so. A command that rewrites part of a
  large existing file causes the kernel to copy it up whole first.
* The journal: overlay records keep no pre-image in the undo store (the
  trash holds it), so they do not count toward its 500-entry cap; they
  have a cap of their own (``change_journal.MAX_OVERLAY_ENTRIES``).
  ``undo_changes`` renames the original back and the replaced version
  into the trash, so an undo is itself undoable.

Fallback: a root on which the overlay cannot be mounted (the probe
failed: no unprivileged user namespaces, a filesystem without user
xattrs such as NFS), a root with another filesystem mounted beneath it
(an overlay does not cross mounts, the command would see an empty
directory), a path the mount option string cannot carry, too little free
disk -- that root is left out and keeps today's behaviour, the shell
change journal plus the trash refusal for what nothing could restore. If
the mount fails at run time the command has not run; the executor checks
the refusal and runs it the old way. A command is never run unconfined
for want of an overlay.

Not seen: hard links between files the command copied up are split;
extended attributes and a directory's own mtime are not carried back;
a change made under the real root by another process while the command
runs is visible to the command or not (the kernel leaves that open), and
is never overwritten.

Rejected alternatives: snapshotting every file before every command (the
shell change journal: bounded at 1 MiB per file and 64 MiB in all,
because it copies bytes); reflink or hard-link copies of the tree (ZFS
has no reflink here; a hard link does not protect an in-place write);
bubblewrap's own ``--overlay`` (0.11+; this server has 0.9, and the
upper layer must be on the root's filesystem for the apply to be a
rename); keeping the upper layer for the whole turn (see above);
putting the layers outside the root (only the root's filesystem gives a
rename, and the session has no other writable place on it).
"""

from __future__ import annotations

import errno
import hashlib
import itertools
import os
import shutil
import stat
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from . import overlay_exec as _ox
from .unrestorable import TRASH, TRASH_DAYS, prepare_trash

LAYERS = ".layers"
MIN_FREE_BYTES = 1 << 30
MIN_FREE_FRACTION = 0.02
#: A directory's fingerprint (the undo guard of a whole tree) is taken up
#: to this many entries; beyond, the directory's inode guards it.
MAX_TREE_ENTRIES = 20_000
#: Created or changed files up to this size get a sha256 post-image hash
#: (the undo guard and the self-written detection); larger ones are
#: guarded by (size, mtime, inode).
MAX_HASH_BYTES = 64 * 1024 * 1024
#: Environment switch: "off" never overlays (an operator's escape hatch).
ENV = "DELFIN_WORKSPACE_OVERLAY"

_LOCK = threading.Lock()
_PROBED: dict[int, str] = {}          # st_dev -> "" (works) or the reason
_RUNTIME_FAILURE: str = ""            # set when a mount failed at run time
_COUNTER = itertools.count(1)


@dataclass
class Layer:
    root: str
    base: str            # R/.delfin-trash/.layers/<id>
    upper: str
    work: str
    trash: str           # R/.delfin-trash/<date>/cmd-<id> (made on demand)


@dataclass
class Plan:
    cmd_id: str
    start_ns: int
    layers: list[Layer] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def covers(self, path) -> bool:
        return any(_under(str(path), l.root) for l in self.layers)


@dataclass
class Applied:
    path: str
    kind: str                 # created | modified | deleted
    entry: str                # file | dir | symlink | other (the new one)
    root: str
    trash: str = ""           # where the original is now ("" if none)
    post_hash: str = ""
    post_sig: list = field(default_factory=list)
    post_tree: str = ""
    journal: bool = True


@dataclass
class Outcome:
    changes: list[Applied] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    kept: list[str] = field(default_factory=list)   # layers left unapplied


# ---------------------------------------------------------------------------
# Eligibility
# ---------------------------------------------------------------------------

def _under(path: str, root: str) -> bool:
    root = root.rstrip(os.sep)
    return path == root or path.startswith(root + os.sep)


def _mount_points() -> list[str]:
    out: list[str] = []
    try:
        with open("/proc/self/mountinfo", encoding="utf-8",
                  errors="replace") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) > 4:
                    out.append(_unescape(parts[4]))
    except OSError:
        pass
    return out


def _unescape(text: str) -> str:
    return (text.replace("\\040", " ").replace("\\011", "\t")
            .replace("\\012", "\n").replace("\\134", "\\"))


def _probe(root: str) -> str:
    """"" when an overlay mounts on ``root``'s filesystem, else why not.
    Cached per filesystem; costs one helper process the first time."""
    try:
        dev = os.stat(root).st_dev
    except OSError as exc:
        return f"the root cannot be read ({exc.strerror})"
    with _LOCK:
        if dev in _PROBED:
            return _PROBED[dev]
    reason = ""
    base = os.path.join(root, TRASH, LAYERS, f"probe-{os.getpid()}")
    try:
        _trash_day(root)
        os.makedirs(base, exist_ok=True)
        proc = subprocess.run(
            [sys.executable, "-I", "-S", os.path.abspath(_ox.__file__),
             "--probe", base], capture_output=True, text=True, timeout=20)
        if proc.returncode != 0:
            reason = ((proc.stderr or "").strip().replace(_ox.MARK, "")
                      .strip() or f"probe exited {proc.returncode}")
            reason = f"an overlay cannot be mounted here ({reason[:160]})"
    except Exception as exc:
        reason = f"the overlay probe failed ({exc})"
    finally:
        _remove_tree(base)
    with _LOCK:
        _PROBED[dev] = reason
    return reason


_DAYS: dict[tuple[str, str], str] = {}


def _trash_day(root: str) -> str:
    """``unrestorable.prepare_trash`` once per root and UTC day (it starts
    a version-control process to keep the trash out of the index, and
    sweeps old days); the day folder after that."""
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    with _LOCK:
        day = _DAYS.get((root, today))
    if day and os.path.isdir(day):
        return day
    day = str(prepare_trash(Path(root)))
    with _LOCK:
        _DAYS[(root, today)] = day
    return day


def _root_reason(root: str, mounts: list[str]) -> str:
    """"" when ``root`` can be overlaid, else why not."""
    if any(c in root for c in _ox.UNSAFE):
        return "its path holds a character the overlay mount cannot carry"
    for m in mounts:
        if m != root and _under(m, root) and not _under(
                m, os.path.join(root, TRASH)):
            return (f"another filesystem is mounted beneath it ({m}); an "
                    "overlay does not reach into it")
    try:
        st = os.statvfs(root)
        free, total = st.f_bavail * st.f_frsize, st.f_blocks * st.f_frsize
        if free < MIN_FREE_BYTES or (total and free < total * MIN_FREE_FRACTION):
            return (f"too little free disk for the trash ({free >> 20} MiB "
                    "free)")
    except OSError:
        pass
    return _probe(root)


def enabled() -> bool:
    return (sys.platform.startswith("linux")
            and os.environ.get(ENV, "").strip().lower() != "off"
            and not _RUNTIME_FAILURE)


def eligible_roots(roots, skip=(), exclude=()) -> tuple[list[str], dict]:
    """(roots to overlay, {root: reason} for the ones left out).

    Only the outermost of nested roots is mounted: a root inside it is
    reached through its overlay. ``skip``: read-only and calc roots, which
    are never overlaid; ``exclude``: roots left out silently (the
    session's scratch)."""
    def norm(items):
        out = []
        for r in items or ():
            try:
                p = os.path.realpath(str(r))
            except (OSError, ValueError):
                continue
            if p not in out:
                out.append(p)
        return out
    skip_n, excl_n = norm(skip), norm(exclude)
    cand = [r for r in norm(roots) if os.path.isdir(r)
            and r not in excl_n and not any(_under(r, s) for s in skip_n)]
    outer = [r for r in cand
             if not any(o != r and _under(r, o) for o in cand)]
    mounts = _mount_points()
    ok, left = [], {}
    for r in outer:
        why = _root_reason(r, mounts)
        if why:
            left[r] = why
        else:
            ok.append(r)
    return ok, left


# ---------------------------------------------------------------------------
# Prepare / wrap
# ---------------------------------------------------------------------------

def prepare(roots: list[str]) -> Optional[Plan]:
    """The layer directories for one command on each root; None when no
    root could be prepared."""
    now = datetime.now(timezone.utc)
    cmd_id = f"{now.strftime('%H%M%S')}-{os.getpid()}-{next(_COUNTER)}"
    plan = Plan(cmd_id=cmd_id, start_ns=0)
    for root in roots:
        try:
            day = _trash_day(root)
            base = os.path.join(root, TRASH, LAYERS, cmd_id)
            upper, work = os.path.join(base, "u"), os.path.join(base, "w")
            os.makedirs(upper)
            os.makedirs(work)
            os.chmod(os.path.join(root, TRASH, LAYERS), 0o700)
            _sweep_layers(os.path.join(root, TRASH, LAYERS))
            plan.layers.append(Layer(root=root, base=base, upper=upper,
                                     work=work,
                                     trash=os.path.join(str(day),
                                                        f"cmd-{cmd_id}")))
        except OSError as exc:
            plan.notes.append(f"{root}: not overlaid ({exc.strerror or exc})")
    if not plan.layers:
        return None
    # The start: a path whose ctime is later was changed by someone else
    # while the command ran. Taken from the filesystem's own clock (the
    # ctime of a file just touched), not from time.time().
    plan.start_ns = _fs_now(plan.layers[0].base)
    return plan


def _fs_now(base: str) -> int:
    mark = os.path.join(base, "start")
    try:
        with open(mark, "w"):
            pass
        return os.stat(mark).st_ctime_ns
    except OSError:
        return time.time_ns()
    finally:
        try:
            os.unlink(mark)
        except OSError:
            pass


def can_wrap(argv: list[str]) -> bool:
    """The argv is the full bubblewrap cage (the "/" bound read-only),
    which is what the overlay is mounted beneath."""
    return list(argv[:4]) == ["bwrap", "--ro-bind", "/", "/"]


def _status_file(plan: Plan) -> str:
    return os.path.join(plan.layers[0].base, "status")


def wrap(argv: list[str], plan: Plan) -> list[str]:
    out = [sys.executable, "-I", "-S", os.path.abspath(_ox.__file__),
           "--status", _status_file(plan)]
    for l in plan.layers:
        out += ["--layer", l.root, l.upper, l.work]
    return out + ["--", *argv]


def setup_failed(plan: Plan) -> bool:
    """The helper did not report every mount in place: the command has not
    run. Read from the status file the helper opened before any mount, not
    from the exit code or stderr, which the command could imitate."""
    try:
        with open(_status_file(plan), "rb") as fh:
            return fh.read(8) != b"ok"
    except OSError:
        return True


def note_runtime_failure(stderr: str) -> str:
    """Remember for this process that the mount fails; the reason."""
    global _RUNTIME_FAILURE
    line = next((l for l in (stderr or "").splitlines() if _ox.MARK in l), "")
    _RUNTIME_FAILURE = line.replace(_ox.MARK, "").strip() or "setup failed"
    return _RUNTIME_FAILURE


def discard(plan: Optional[Plan]) -> None:
    for l in (plan.layers if plan else ()):
        _remove_tree(l.base)


# ---------------------------------------------------------------------------
# Apply
# ---------------------------------------------------------------------------

def _is_whiteout(st) -> bool:
    return stat.S_ISCHR(st.st_mode) and st.st_rdev == 0


def _opaque(path: str) -> bool:
    for name in ("user.overlay.opaque", "trusted.overlay.opaque"):
        try:
            if os.getxattr(path, name, follow_symlinks=False) == b"y":
                return True
        except OSError:
            continue
    return False


def _strip_xattrs(path: str) -> None:
    try:
        names = os.listxattr(path, follow_symlinks=False)
    except OSError:
        return
    for n in names:
        if n.startswith(("user.overlay.", "trusted.overlay.")):
            try:
                os.removexattr(path, n, follow_symlinks=False)
            except OSError:
                pass


def _entry(st) -> str:
    if stat.S_ISLNK(st.st_mode):
        return "symlink"
    if stat.S_ISDIR(st.st_mode):
        return "dir"
    if stat.S_ISREG(st.st_mode):
        return "file"
    return "other"


def _lstat(path: str):
    try:
        return os.lstat(path)
    except OSError:
        return None


def file_hash(path: str, size: int) -> str:
    if size > MAX_HASH_BYTES:
        return ""
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return ""
    return h.hexdigest()


def tree_print(path: str) -> str:
    """sha256 over (relative path, type, size, mtime) of every entry under
    ``path``; "" beyond ``MAX_TREE_ENTRIES``."""
    h = hashlib.sha256()
    n = 0
    for dirpath, dirs, files in os.walk(path, followlinks=False):
        dirs.sort()
        for name in sorted(dirs + files):
            n += 1
            if n > MAX_TREE_ENTRIES:
                return ""
            p = os.path.join(dirpath, name)
            st = _lstat(p)
            if st is None:
                continue
            h.update(f"{os.path.relpath(p, path)}\0{st.st_mode}\0"
                     f"{st.st_size}\0{st.st_mtime_ns}\n".encode(
                         "utf-8", "surrogateescape"))
    return h.hexdigest()


def post_state(path: str) -> tuple[str, str, list, str]:
    """(entry, post_hash, post_sig, post_tree) of what is at ``path`` now."""
    st = _lstat(path)
    if st is None:
        return "", "", [], ""
    kind = _entry(st)
    sig = [st.st_size, st.st_mtime_ns, st.st_ino]
    if kind == "file":
        return kind, file_hash(path, st.st_size), sig, ""
    if kind == "symlink":
        try:
            return kind, hashlib.sha256(os.fsencode(os.readlink(path))
                                        ).hexdigest(), sig, ""
        except OSError:
            return kind, "", sig, ""
    if kind == "dir":
        return kind, "", sig, tree_print(path)
    return kind, "", sig, ""


def same_content(a: str, b: str, size: int) -> bool:
    try:
        with open(a, "rb") as fa, open(b, "rb") as fb:
            while True:
                x, y = fa.read(1 << 20), fb.read(1 << 20)
                if x != y:
                    return False
                if not x:
                    return True
    except OSError:
        return False


class _Applier:
    def __init__(self, plan: Plan, layer: Layer, out: Outcome):
        self.plan, self.layer, self.out = plan, layer, out
        self.failed: list[str] = []

    # Trash: the original keeps its relative path under cmd-<id>/.
    def to_trash(self, real: str, rel: str) -> str:
        dest = os.path.join(self.layer.trash, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        _no_link(os.path.dirname(dest))
        _no_link(os.path.dirname(real))
        os.rename(real, dest)
        return dest

    def newer(self, st) -> bool:
        return st is not None and st.st_ctime_ns > self.plan.start_ns

    def record(self, real: str, kind: str, trash: str = "",
               journal: bool = True) -> None:
        if not journal:
            entry, ph, sig, tree = "file", "", [], ""
        elif kind == "deleted":
            st = _lstat(trash) if trash else None
            entry, ph, sig, tree = (_entry(st) if st else ""), "", [], ""
        else:
            entry, ph, sig, tree = post_state(real)
        self.out.changes.append(Applied(
            path=real, kind=kind, entry=entry or "file", root=self.layer.root,
            trash=trash, post_hash=ph, post_sig=sig, post_tree=tree,
            journal=journal))

    def place(self, src: str, real: str, rel: str, lst) -> None:
        """Move the command's ``src`` to ``real``; the original (``lst``,
        None when absent) to the trash first."""
        if lst is not None and self.newer(lst):
            side = f"{real}.delfin-conflict-{self.plan.cmd_id}"
            _no_link(os.path.dirname(side))
            os.rename(src, side)
            _strip_tree_xattrs(side)
            self.out.notes.append(
                f"{real} was changed by someone else while the command ran; "
                f"it was kept, and the command's version is beside it as "
                f"{os.path.basename(side)}")
            self.record(side, "created")
            return
        trash = self.to_trash(real, rel) if lst is not None else ""
        _no_link(os.path.dirname(real))
        os.rename(src, real)
        _strip_tree_xattrs(real)
        objects = rel.startswith(".git" + os.sep + "objects" + os.sep)
        self.record(real, "modified" if lst is not None else "created",
                    trash=trash, journal=not (objects and lst is None))

    def remove(self, real: str, rel: str, lst) -> None:
        if lst is None:
            return
        if self.newer(lst):
            self.out.notes.append(
                f"{real} was changed by someone else while the command ran; "
                "the command deleted it, and that deletion was not applied")
            return
        trash = self.to_trash(real, rel)
        self.record(real, "deleted", trash=trash)

    def walk(self, udir: str, rdir: str, rel: str, opaque: bool) -> None:
        try:
            names = sorted(os.listdir(udir))
        except OSError as exc:
            raise OSError(exc.errno, f"cannot list {udir}: {exc.strerror}")
        if opaque:
            # Everything the command's directory does not hold was
            # removed by it (rm -r d; mkdir d).
            try:
                lower = sorted(os.listdir(rdir))
            except OSError:
                lower = []
            for name in lower:
                if name in names:
                    continue
                if not rel and name == TRASH:
                    continue
                try:
                    self.remove(os.path.join(rdir, name),
                                os.path.join(rel, name),
                                _lstat(os.path.join(rdir, name)))
                except OSError as exc:
                    self.failed.append(f"{os.path.join(rdir, name)}: "
                                       f"{exc.strerror or exc}")
        for name in names:
            r = os.path.join(rel, name) if rel else name
            if not rel and name == TRASH:
                continue            # the trash is never the command's to change
            try:
                self.entry(udir, rdir, name, r)
            except OSError as exc:
                self.failed.append(f"{os.path.join(rdir, name)}: "
                                   f"{exc.strerror or exc}")

    def entry(self, udir: str, rdir: str, name: str, r: str) -> None:
        src, real = os.path.join(udir, name), os.path.join(rdir, name)
        ust, lst = _lstat(src), _lstat(real)
        if ust is None:
            return
        if _is_whiteout(ust):
            self.remove(real, r, lst)
            return
        if stat.S_ISDIR(ust.st_mode):
            if lst is not None and stat.S_ISDIR(lst.st_mode):
                if stat.S_IMODE(ust.st_mode) != stat.S_IMODE(lst.st_mode):
                    try:
                        os.chmod(real, stat.S_IMODE(ust.st_mode))
                    except OSError:
                        pass
                self.walk(src, real, r, _opaque(src))
                return
            self.place(src, real, r, lst)
            return
        if (lst is not None and stat.S_ISREG(ust.st_mode)
                and stat.S_ISREG(lst.st_mode) and not self.newer(lst)
                and ust.st_size == lst.st_size
                and same_content(src, real, ust.st_size)):
            # Copied up, bytes unchanged (opened for writing, touched,
            # chmod-ed): the metadata only, in place.
            try:
                if stat.S_IMODE(ust.st_mode) != stat.S_IMODE(lst.st_mode):
                    os.chmod(real, stat.S_IMODE(ust.st_mode))
                os.utime(real, ns=(ust.st_atime_ns, ust.st_mtime_ns))
            except OSError:
                pass
            return
        self.place(src, real, r, lst)


def _no_link(directory: str) -> None:
    """Refuse to move anything through a directory reached over a
    symbolic link: the roots are resolved, so a resolved path that differs
    means a link was put in the way after the command ran."""
    if os.path.realpath(directory) != directory:
        raise OSError(errno.ELOOP,
                      f"{directory} is reached through a symbolic link")


def _strip_tree_xattrs(path: str) -> None:
    _strip_xattrs(path)
    if os.path.isdir(path) and not os.path.islink(path):
        n = 0
        for dirpath, dirs, files in os.walk(path):
            for name in dirs + files:
                n += 1
                if n > MAX_TREE_ENTRIES:
                    return
                _strip_xattrs(os.path.join(dirpath, name))


def apply(plan: Optional[Plan]) -> Outcome:
    """Move the command's changes into the real roots. Never raises; a
    layer that could not be applied whole is kept and named."""
    out = Outcome()
    if plan is None:
        return out
    for layer in plan.layers:
        app = _Applier(plan, layer, out)
        try:
            app.walk(layer.upper, layer.root, "", False)
        except Exception as exc:
            app.failed.append(str(exc))
        if app.failed:
            out.kept.append(layer.base)
            out.notes.append(
                f"the command's changes under {layer.root} could not all be "
                f"applied ({'; '.join(app.failed[:3])}); what was not is "
                f"kept in {layer.upper}")
            continue
        _remove_tree(layer.base)
        try:
            os.rmdir(layer.trash)          # nothing was trashed
        except OSError:
            pass
    return out


# ---------------------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------------------

def _remove_tree(path: str) -> None:
    """rmtree that first makes the kernel's mode-000 work entries
    writable. Never raises."""
    if not os.path.lexists(path):
        return
    for dirpath, dirs, _files in os.walk(path):
        for d in dirs:
            p = os.path.join(dirpath, d)
            try:
                if not os.path.islink(p):
                    os.chmod(p, 0o700)
            except OSError:
                pass
    try:
        os.chmod(path, 0o700)
    except OSError:
        pass
    shutil.rmtree(path, ignore_errors=True)


def _sweep_layers(layers_dir: str, days: int = TRASH_DAYS) -> None:
    """Layers left by a run that ended before applying them are kept as
    long as the trash keeps a file, then removed."""
    cutoff = time.time() - days * 86400
    try:
        for e in os.scandir(layers_dir):
            try:
                if e.stat(follow_symlinks=False).st_mtime < cutoff:
                    _remove_tree(e.path)
            except OSError:
                continue
    except OSError:
        pass


def undo_trash_dir(root: str, seq: int) -> str:
    """Where an undo puts what it replaces: the trash, like a command."""
    day = _trash_day(root)
    return os.path.join(str(day), f"undo-{time.strftime('%H%M%S')}-{seq}")


__all__ = ["Plan", "Layer", "Applied", "Outcome", "enabled",
           "eligible_roots", "prepare", "can_wrap", "wrap", "apply",
           "discard", "setup_failed", "note_runtime_failure", "post_state",
           "undo_trash_dir", "LAYERS", "ENV"]
