"""Ensure every bind target of a cage argv exists before bwrap runs it.

The cage argv builders (``api_client``: ``_bash_isolation_argv``,
``_test_run_argv``) bind or mask absolute paths. A mount target under the
read-only ``--ro-bind / /`` root must exist on the host when bwrap sets the
mount up: for a plain mount bwrap creates the missing leaf INSIDE the new
root, and under a read-only parent that mkdir fails with
``bwrap: Can't mkdir parents for <path>: Read-only file system`` -- the
command dies before it runs. Whether a target exists depends on what the
session did earlier (first approval code, first HOME refresh), which made
the failure look random (measured, wave 31).

Rule: before an argv is handed to bwrap, every target that a mask or bind
names and that the session may create is made to exist, with the owner's
mode -- 0700 for anything under ``~/.delfin``. Nothing is degraded to a
``-try`` mount: a skipped ``--tmpfs-try`` mask would re-expose the very
path the mask hides. A target the session must not create keeps its plain
mount only when it exists; the builders already skip such targets at build
time (``sandbox`` helpers check ``os.path.isdir/isfile``), and this module
re-checks them at the exec boundary -- a target that vanished between build
and exec is re-created where creation is right.

Public entry points:

* :func:`ensure_argv_bind_targets` -- walk a finished argv, create what is
  missing. Called by the argv builders right before they return; never
  raises (a cage that cannot be prepared is not made better by a traceback).

Tested by ``tests/test_bw1_cage_binds.py`` (unit, runs everywhere) and
``tests/test_bw1_missing_bind_targets.py`` (real bwrap, operator run)."""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["ensure_argv_bind_targets", "BIND_FLAGS", "MASK_FLAGS"]

#: Mount flags whose first operand (the HOST source / target) must exist
#: when bwrap sets the mount up.
BIND_FLAGS = ("--tmpfs", "--ro-bind", "--bind", "--dev-bind", "--ro-bind-try",
              "--bind-try")

#: Mount flags that MASK a path (hide it). A mask is never skipped and
#: never weakens: its target is created, never degraded.
MASK_FLAGS = ("--tmpfs",)

#: First path components under HOME that are private to the owner (the
#: home secret doors the cage masks: .ssh, .gnupg, ... -- anything a mask
#: or door target sits under directly). A target created here gets the
#: owner-only mode.
_MODE_0700_ROOTS = (".delfin", ".ssh", ".aws", ".gnupg", ".kube",
                    ".docker", ".azure", ".gcp", ".config", ".screen",
                    ".anthropic", ".openai", ".claude", ".local")


def _owner_mode(target: str) -> int:
    """Mode a created target gets: 0700 for a private home area (under
    ~/.delfin, a secret door of its own like ~/.ssh, or a path directly
    under HOME -- a door target _process_cage_options masks), else 0755."""
    try:
        home = str(Path.home().resolve())
    except OSError:
        home = ""
    t = str(Path(target).resolve()) if os.path.isabs(target) else target
    if home and t.startswith(home.rstrip("/") + "/"):
        rel = t[len(home.rstrip("/")) + 1:]
        first = rel.split("/")[0]
        if first in _MODE_0700_ROOTS:
            return 0o700
        if "/" not in rel:
            return 0o700  # a leaf directly under HOME is a private door
    return 0o755


def _mkdir_parents(target: str, mode: int) -> bool:
    """Create the target's missing parents and the target itself. True when
    the target exists afterwards (it did before, or was just made)."""
    try:
        p = Path(target)
        if p.exists():
            return True
        p.mkdir(parents=True, exist_ok=True, mode=mode)
        # mkdir's mode is masked by umask; state the mode explicitly.
        os.chmod(p, mode)
        return True
    except OSError:
        return False


def ensure_argv_bind_targets(argv: list[str]) -> list[str]:
    """Make every createable bind target in ``argv`` exist; return ``argv``.

    Walks the mount flags and, for each target that does not exist, creates
    it with the owner's mode (0700 under ``~/.delfin``). Never raises and
    never changes the argv: a target that cannot be created is left to
    bwrap's own refusal (which names the path and is diagnosable), and a
    mask is never skipped. Safe to call on any argv: only mount-flag
    operands are read, the command tail is untouched."""
    for i, tok in enumerate(argv):
        if tok == "--":
            break  # the command tail is never a mount
        if tok in BIND_FLAGS:
            _ensure_flag_target(argv, i, tok)
    return argv


def _ensure_flag_target(argv: list[str], i: int, flag: str) -> None:
    """Create the target(s) one mount flag names, if missing.

    ``--tmpfs <target>``: one operand. ``--*-bind <src> <dst>``: the
    destination is the in-cage path -- it must exist under the new root and
    is covered by the read-only "/" bind, so the HOST path to create is the
    same string (the builders bind paths at themselves). ``--dev-bind`` and
    ``--ro-bind /dev/null <file>`` follow the same rule."""
    rest = argv[i + 1:]
    if not rest:
        return
    if flag == "--tmpfs":
        if os.path.isabs(rest[0]):
            _mkdir_parents(rest[0], _owner_mode(rest[0]))
        return
    src = rest[0] if os.path.isabs(rest[0]) else None
    dst = rest[1] if len(rest) > 1 and os.path.isabs(rest[1]) else None
    if src and dst and _is_filelike(src):
        # A file-type mask (--ro-bind /dev/null <file>): the destination
        # must be a regular file too, or bwrap refuses the type mismatch.
        _ensure_file(dst)
        _mkdir_parents(src, _owner_mode(src))
        return
    for t in (a for a in rest[:2] if os.path.isabs(a)):
        if not _mkdir_parents(t, _owner_mode(t)):
            # Not createable (no permission, rofs): leave it to bwrap's
            # own refusal -- never degrade the mount to a -try form.
            continue


def _is_filelike(path: str) -> bool:
    """A regular file or a device node: the bind source is file-type, so
    its destination must be a file (bwrap refuses a dir for /dev/null)."""
    try:
        if os.path.isfile(path):
            return True
        return Path(path).is_char_device()
    except OSError:
        return False


def _ensure_file(path: str) -> bool:
    """Create ``path``'s parents and an empty regular file at the leaf."""
    try:
        p = Path(path)
        if p.is_file():
            return True
        if p.is_dir():
            return True  # a dir over a file mask still hides: tmpfs-like
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch(exist_ok=True)
        os.chmod(p, _owner_mode(str(p.parent)))
        return True
    except OSError:
        return False
