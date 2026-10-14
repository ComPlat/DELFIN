"""What a shell command changed in the writable roots.

Input: the writable roots of a session (``perms.all_workspace_roots()``).
Output: per command, the files created, modified and deleted under them,
each with the pre-image held from before the command where one could be
held. The executor journals every one of them (``change_journal``), so a
write that bypassed write_file -- an interpreter, a script, a compiler --
is undoable, reaches the completion check and the self-written
detection like any file-tool write.

Mechanism: :func:`before` walks the roots and records ``(size, mtime_ns,
inode, ctime_ns)`` per regular file, and holds the bytes of files up to
``MAX_PRE_IMAGE_BYTES``, newest first, up to ``MAX_HELD_BYTES`` in all.
Held bytes are cached across commands by the same signature, so a file
is read once and again only after it changed: in steady state a command
reads what the previous one changed, and nothing else. :func:`after`
walks again and compares signatures; only files whose signature moved
are read, and one whose bytes came out identical is not a change.

Bounds, each said rather than silently applied:

* ``.git``, ``__pycache__`` and ``.delfin`` subdirectories are not
  walked (git internals, bytecode, the agent's own state); a directory
  with more than ``MAX_ENTRIES_PER_DIR`` entries is not walked either.
* a root with more than ``MAX_FILES_PER_ROOT`` files, or whose walk
  takes longer than ``MAX_WALK_SECONDS``, is not tracked; the note says
  so once per session and root.
* a changed file whose pre-image was not held (over the per-file cap,
  or the held budget was spent) is journalled WITHOUT undo, with the
  reason, rather than not at all.
* a created file over ``MAX_HASH_BYTES`` is journalled without undo:
  its post-image hash, which guards the undo, is not computed.

Not seen: directories (an empty directory left behind by an undone
creation stays), symlinks (not followed, not recorded), a change that
leaves size, mtime, inode and ctime all equal, and writes outside the
roots (which the cage refuses where it runs). A change another process
makes during the command is attributed to the command; one the session
journalled itself during the command (seq after the start) is not
journalled twice.

Rejected alternatives: reading every file before every command (2 000
files: the full read each time, where the cache reads only what moved);
inotify (Linux only, no pre-image); an overlay mount in the cage (needs
bubblewrap, and the journal has to hold where bubblewrap does not run).
"""

from __future__ import annotations

import os
import threading
import time
from dataclasses import dataclass, field
from typing import Optional

MAX_FILES_PER_ROOT = 20_000
MAX_ENTRIES_PER_DIR = 5_000
MAX_WALK_SECONDS = 2.0
MAX_PRE_IMAGE_BYTES = 1024 * 1024        # per file, held before the run
MAX_HELD_BYTES = 64 * 1024 * 1024        # all held pre-images together
MAX_HASH_BYTES = 64 * 1024 * 1024        # post-image hashed for the guard

_SKIP_DIR_NAMES = frozenset({".git", "__pycache__", ".delfin"})

_Sig = tuple  # (size, mtime_ns, ino, ctime_ns)

_LOCK = threading.Lock()
# path -> (sig, bytes): pre-images read for an earlier snapshot.
_HELD: dict[str, tuple[_Sig, bytes]] = {}
_HELD_TOTAL = 0
# (session, root, reason) already said.
_SAID: set[tuple[str, str, str]] = set()


@dataclass
class Snapshot:
    """The state of the roots before one command."""
    roots: list[str]
    files: dict[str, _Sig] = field(default_factory=dict)
    held: dict[str, bytes] = field(default_factory=dict)
    # Roots not tracked for this command: root -> reason.
    untracked: dict[str, str] = field(default_factory=dict)
    start_seq: int = 0
    session_id: str = ""


@dataclass
class Change:
    path: str
    kind: str                     # "created" | "modified" | "deleted"
    pre: Optional[bytes] = None   # None: created, or not held
    reason: str = ""              # why there is no undo ("" if there is)


def _walk(root: str, other_roots: tuple[str, ...], deadline: float,
          skipped: dict[str, str]) -> tuple[dict[str, _Sig], str]:
    """(files, reason) under one root; reason non-empty = not tracked.
    A directory over ``MAX_ENTRIES_PER_DIR`` goes into ``skipped``."""
    files: dict[str, _Sig] = {}
    stack = [root]
    while stack:
        d = stack.pop()
        try:
            with os.scandir(d) as it:
                entries = list(it)
        except OSError:
            continue
        if len(entries) > MAX_ENTRIES_PER_DIR and d != root:
            skipped[d] = (f"too many files to track (more than "
                          f"{MAX_ENTRIES_PER_DIR} entries in one directory)")
            continue
        for e in entries:
            try:
                if e.is_dir(follow_symlinks=False):
                    if e.name in _SKIP_DIR_NAMES or e.path in other_roots:
                        continue
                    stack.append(e.path)
                    continue
                if not e.is_file(follow_symlinks=False):
                    continue
                st = e.stat(follow_symlinks=False)
            except OSError:
                continue
            files[e.path] = (st.st_size, st.st_mtime_ns, st.st_ino,
                             st.st_ctime_ns)
            if len(files) > MAX_FILES_PER_ROOT:
                return {}, (f"too many files to track (more than "
                            f"{MAX_FILES_PER_ROOT})")
        if time.monotonic() > deadline:
            return {}, (f"too many files to track (the walk took longer "
                        f"than {MAX_WALK_SECONDS:g} s)")
    return files, ""


def _norm_roots(roots) -> list[str]:
    out: list[str] = []
    for r in roots or ():
        try:
            p = os.path.realpath(str(r))
        except (OSError, ValueError):
            continue
        if os.path.isdir(p) and p not in out:
            out.append(p)
    return out


def _read(path: str, limit: int) -> Optional[bytes]:
    try:
        with open(path, "rb") as fh:
            data = fh.read(limit + 1)
    except OSError:
        return None
    return data if len(data) <= limit else None


def before(roots, *, session_id: str = "", start_seq: int = 0) -> Snapshot:
    """Snapshot the roots before a command. Never raises."""
    global _HELD_TOTAL
    snap = Snapshot(roots=_norm_roots(roots), start_seq=int(start_seq or 0),
                    session_id=str(session_id or ""))
    try:
        all_roots = tuple(snap.roots)
        for root in snap.roots:
            others = tuple(r for r in all_roots if r != root)
            files, reason = _walk(root, others,
                                  time.monotonic() + MAX_WALK_SECONDS,
                                  snap.untracked)
            if reason:
                snap.untracked[root] = reason
                continue
            snap.files.update(files)

        # Pre-images, newest first: a command most often rewrites what
        # was written last.
        order = sorted(
            (p for p, s in snap.files.items() if s[0] <= MAX_PRE_IMAGE_BYTES),
            key=lambda p: snap.files[p][1], reverse=True)
        budget = MAX_HELD_BYTES
        with _LOCK:
            cached = dict(_HELD)
        for p in order:
            sig = snap.files[p]
            if sig[0] > budget:
                continue
            hit = cached.get(p)
            if hit is not None and hit[0] == sig:
                data: Optional[bytes] = hit[1]
            else:
                data = _read(p, MAX_PRE_IMAGE_BYTES)
                if data is None or len(data) != sig[0]:
                    continue
            snap.held[p] = data
            budget -= len(data)
        with _LOCK:
            prefixes = tuple(r.rstrip(os.sep) + os.sep for r in snap.roots)
            for p in [p for p in _HELD if p.startswith(prefixes)]:
                _HELD_TOTAL -= len(_HELD.pop(p)[1])
            if _HELD_TOTAL > MAX_HELD_BYTES:
                _HELD.clear()
                _HELD_TOTAL = 0
            for p, data in snap.held.items():
                _HELD[p] = (snap.files[p], data)
                _HELD_TOTAL += len(data)
    except Exception:
        snap.files.clear()
        snap.held.clear()
        for r in snap.roots:
            snap.untracked.setdefault(r, "the snapshot failed")
    return snap


def after(snap: Snapshot) -> list[Change]:
    """Files under the tracked roots whose content moved since ``snap``.

    A root that was tracked before but is over a bound now (the command
    created many files) yields nothing and is added to
    ``snap.untracked``: a partial walk would report the rest as deleted.
    Never raises."""
    changes: list[Change] = []
    try:
        all_roots = tuple(snap.roots)
        now: dict[str, _Sig] = {}
        for root in snap.roots:
            if root in snap.untracked:
                continue
            others = tuple(r for r in all_roots if r != root)
            files, reason = _walk(root, others,
                                  time.monotonic() + MAX_WALK_SECONDS,
                                  snap.untracked)
            if reason:
                snap.untracked[root] = reason
                continue
            now.update(files)
        tracked = tuple(r.rstrip(os.sep) + os.sep for r in snap.roots
                        if r not in snap.untracked)
        if not tracked:
            return []
        # A directory skipped before OR after is unseen on one side: its
        # files would read as created or deleted.
        unseen = tuple(d.rstrip(os.sep) + os.sep for d in snap.untracked)

        for p in sorted(set(snap.files) | set(now)):
            if not p.startswith(tracked) or (unseen and p.startswith(unseen)):
                continue
            old, new = snap.files.get(p), now.get(p)
            if old == new:
                continue
            held = snap.held.get(p)
            if old is None:
                if new[0] > MAX_HASH_BYTES:
                    changes.append(Change(p, "created", reason=(
                        f"too large to guard an undo ({new[0]} bytes > "
                        f"{MAX_HASH_BYTES})")))
                else:
                    changes.append(Change(p, "created"))
                continue
            if new is None:
                changes.append(Change(p, "deleted", pre=held, reason="" if
                                      held is not None else _no_pre(old)))
                continue
            if held is None:
                if old[:3] == new[:3]:
                    continue        # ctime only: metadata, not content
                changes.append(Change(p, "modified", reason=_no_pre(old)))
                continue
            post = _read(p, MAX_PRE_IMAGE_BYTES) if new[0] == len(held) \
                else None
            if post is not None and post == held:
                continue
            changes.append(Change(p, "modified", pre=held))
    except Exception:
        return changes
    return changes


def _no_pre(sig: _Sig) -> str:
    if sig[0] > MAX_PRE_IMAGE_BYTES:
        return (f"too large to keep a copy of before the command "
                f"({sig[0]} bytes > {MAX_PRE_IMAGE_BYTES})")
    return ("no copy was kept before the command (the "
            f"{MAX_HELD_BYTES // (1024 * 1024)} MB held for undo was spent)")


def take_notes(snap: Snapshot) -> list[str]:
    """One line per untracked root, each said once per session."""
    out: list[str] = []
    with _LOCK:
        for root, reason in snap.untracked.items():
            key = (snap.session_id, root, reason.split(" (")[0])
            if key in _SAID:
                continue
            _SAID.add(key)
            out.append(
                f"changes under {root} are not journalled: {reason}; "
                "undo_changes cannot take back what a shell command writes "
                "there")
    return out
