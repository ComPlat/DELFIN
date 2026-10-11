"""patch_targets -- the files a git/patch command would WRITE.

``git apply``, ``git am``, ``patch`` and ``git checkout -- paths`` /
``git restore`` do not name their write targets on the command line the way
``write_file`` and ``tee`` do: ``git apply changes.diff`` writes whatever the
patch's headers name, and the write gate (``api_client._bash_write_targets``)
saw no target at all -- so ``write_allow_globs`` and the workspace boundary
did not hold there, and only a wave-level deny pattern kept those commands
from writing anywhere the account could reach.

This module parses the git/patch command line and, when the patch is a named
file, the patch itself, into the absolute paths the command would write. The
caller then pushes each one through the same per-path write gate (the
executor's ``_gate_write_path``) as any other write.

A command whose write set cannot be read from its text -- a patch piped on
stdin (``git show … | git apply``, ``git apply -``, ``patch -p1 < f.diff``),
or a mailbox (``git am``) -- returns :data:`UNKNOWABLE`. The caller turns
that into a refusal: "we could not check what it would write, so it may not
run." Fail-closed: an undecidable command is refused, never let through.

The parsing is deliberately conservative -- every path a patch might touch is
listed, so the gate can refuse it. A command that lists one more target than
reality only over-refuses a legit operation; missing one would let a write
escape the gate, which is the failure this module exists to stop.
"""

from __future__ import annotations

import os
import re
import shlex
from pathlib import Path


class _Unknowable:
    """Sentinel for 'the write targets cannot be read from the command text'."""

    _instance = None

    def __new__(cls) -> "_Unknowable":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:  # pragma: no cover - repr only
        return "UNKNOWABLE"

    def __bool__(self) -> bool:
        # A sentinel is truthy-neutral; callers test identity, not truthiness.
        return True


#: Returned by :func:`targets` when a command writes but its write set cannot
#: be read from the command text (patch on stdin, pipe, mailbox). The caller
#: must refuse such a command.
UNKNOWABLE = _Unknowable()


#: git subcommands whose working-tree writes are parsed by this module. The
#: set is only documentation: each is dispatched explicitly in
#: :func:`_targets_git`.
_WORKTREE_WRITE_SUBCOMMANDS = frozenset({
    "apply", "am", "checkout", "restore",
})


def _resolve(path: str, cwd: Path) -> str:
    """Resolve a possibly-relative path against *cwd*, as an absolute string."""
    p = os.path.expanduser(path)
    if os.path.isabs(p):
        return p
    return str((cwd / p))


def targets(command: str, cwd: str | os.PathLike | None = None) -> object:
    """The files a git/patch command in *command* would write, absolutized.

    Returns one of:
      * :data:`UNKNOWABLE` -- the command writes but its targets cannot be
        read from the text (patch from stdin / a pipe / ``git am``). Refuse.
      * a list of absolute path strings -- one per file the command would
        add, change or delete. The caller passes each through the write gate.
      * ``[]`` -- the command is not a git/patch working-tree write, or its
        only write is a ``--check``/validation (nothing lands). Nothing to
        add; the caller falls through to its existing behaviour.

    ``cwd`` defaults to the current working directory; relative target paths
    are resolved against it. Never raises; an unparseable command yields
    ``[]`` (the caller's other safeguards still hold).
    """
    base = Path(cwd) if cwd is not None else Path.cwd()
    return _targets_for_command(command, base)


# --- dispatch --------------------------------------------------------------

def _targets_for_command(command: str, base: Path) -> object:
    """Top-level: split the command into shell segments and collect the write
    targets of any git/patch segment. ``[]`` when no segment writes a file.
    A segment returning :data:`UNKNOWABLE` makes the whole command unknowable
    (the caller refuses it): a pipe whose writer reads stdin cannot be told
    what it would write by listing the other segment's targets."""
    collected: list[str] = []
    for segment in _split_segments(command):
        result = _targets_for_segment(segment, base)
        if result is UNKNOWABLE:
            return UNKNOWABLE
        collected.extend(result)
    return _dedup(collected)


def _split_segments(command: str) -> list[str]:
    """Split a shell command on statement separators/``&&``/``||``/newline,
    losing nothing. ``|`` is intentionally NOT a split point: a pipeline
    ``git show x | git apply`` must stay one segment so the writer half is
    seen as reading stdin (unknowable), not as a standalone command."""
    return [s for s in re.split(r"[;&]{1,2}|\n", command) if s.strip()]


def _targets_for_segment(segment: str, base: Path) -> object:
    """Write targets of ONE shell segment. If the segment is a pipeline
    (``A | git apply``), the patch writer reads its patch from stdin, which
    cannot be read from the command text -> UNKNOWABLE (refuse)."""
    if "|" in segment:
        return _targets_pipeline(segment, base)
    try:
        argv = shlex.split(segment, comments=True)
    except ValueError:
        argv = segment.split()
    # Skip leading env-var prefixes (FOO=bar cmd ...).
    while argv and re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", argv[0]):
        argv = argv[1:]
    if not argv:
        return []
    name = Path(argv[0]).name
    rest = argv[1:]

    if name == "git":
        return _targets_git(rest, base)
    if name == "patch":
        return _targets_patch(rest, base)
    if name == "git-apply":
        # `git-apply` is the core built as a standalone binary.
        return _targets_git_apply(rest, base)
    return []


def _targets_pipeline(segment: str, base: Path) -> object:
    """A pipeline's write targets. Any part that is a patch writer WITHOUT a
    named patch file (``git apply``, ``git am``, bare ``patch``) reads its
    patch from the pipe -> UNKNOWABLE. Otherwise the parts are independent and
    collected."""
    collected: list[str] = []
    for part in segment.split("|"):
        part = part.strip()
        if not part:
            continue
        result = _targets_for_segment(part, base)
        if result is UNKNOWABLE:
            return UNKNOWABLE
        collected.extend(result)
    return _dedup(collected)


def _dedup(paths: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _git_subcommand(rest: list[str]) -> str:
    """The git subcommand in *rest* (after ``git``), skipping global options
    that carry a value (``-C dir``, ``-c k=v``, ``--git-dir=``…)."""
    i = 0
    value_opts = {"-C", "-c", "--git-dir", "--work-tree", "--namespace",
                  "--exec-path", "--super-prefix"}
    while i < len(rest):
        arg = rest[i]
        if not arg.startswith("-"):
            return arg
        if "=" in arg:
            i += 1
            continue
        if arg in value_opts:
            i += 2
        else:
            i += 1
    return ""


def _targets_git(rest: list[str], base: Path) -> object:
    """Write targets of a ``git <sub> …`` command. ``rest`` is the argument
    list after ``git`` with any leading env-var prefix already stripped."""
    try:
        sub = _git_subcommand(rest)
    except Exception:
        return []
    if sub == "apply":
        return _targets_git_apply(_after_subcommand(rest, sub), base)
    if sub == "am":
        # A mailbox: the patch lives inside the mbox body, not on the command
        # line. Reading it robustly is not worth the risk of a missed target;
        # refuse the command outright.
        return UNKNOWABLE
    if sub == "checkout":
        return _targets_git_checkout(_after_subcommand(rest, sub), base)
    if sub == "restore":
        return _targets_git_restore(_after_subcommand(rest, sub), base)
    return []


def _after_subcommand(rest: list[str], sub: str) -> list[str]:
    """The argument list after the *sub* command token in *rest*."""
    try:
        idx = rest.index(sub)
    except ValueError:
        return []
    return rest[idx + 1:]


def _targets_git_apply(rest: list[str], base: Path) -> object:
    """git apply [-pN] [--directory=prefix] [-R] [--include=…] [patch].

    The files it writes are read from the patch file named on the command
    line. A patch on stdin (``git apply`` with no file, ``git apply -``, or a
    pipe) is UNKNOWABLE -> the caller refuses. ``--check`` validates without
    writing -> nothing to add."""
    strip = 1
    prefix = None
    check_only = False
    patch_file = None
    i = 0
    while i < len(rest):
        a = rest[i]
        if a == "--check":
            check_only = True
        elif a == "-R" or a == "--reverse":
            pass            # a reverse apply still touches the same a/ b/ files
        elif a.startswith("--directory="):
            prefix = a.split("=", 1)[1]
        elif a == "--directory" and i + 1 < len(rest):
            i += 1
            prefix = rest[i]
        elif a == "-p" and i + 1 < len(rest):
            i += 1
            try:
                strip = int(rest[i])
            except ValueError:
                strip = 1
        elif a.startswith("-p"):
            try:
                strip = int(a[2:])
            except ValueError:
                strip = 1
        elif a.startswith("--include=") or a.startswith("--exclude="):
            pass            # includes/excludes select, they don't name new files
        elif a in ("--include", "--exclude") and i + 1 < len(rest):
            i += 1
        elif not a.startswith("-") and patch_file is None:
            patch_file = a      # the patch file (only one is applied)
        i += 1
    if check_only:
        return []
    if patch_file is None:
        # Patch on stdin (positionalless, `-`, or piped in). Not knowable.
        return UNKNOWABLE
    result = _targets_from_patch(patch_file, strip=strip, prefix=prefix)
    if result is UNKNOWABLE:
        return UNKNOWABLE
    return [_resolve(p, base) for p in result]


def _targets_git_checkout(rest: list[str], base: Path) -> object:
    """git checkout [branch] -- <paths>  /  git checkout <path> ...

    The write targets are the explicit pathspecs. ``--`` separates a branch
    from the paths; without ``--`` the traversal skips git flags and their
    values (``-b <name>``, ``-B <name>``, ``--orphan <name>``, …) and takes
    the remaining operands as paths (path-mode checkout). A branch name over-
    gated as a path is harmless: branches are not writable files."""
    return [_resolve(p, base) for p in _operands_skipping_flags(
        rest, {"-b", "-B", "--branch", "--orphan", "--guess", "-q",
               "--quiet", "-t", "--track", "--no-track", "-m", "--merge",
               "-p", "--patch", "--ours", "--theirs", "-f", "--force"})]


def _targets_git_restore(rest: list[str], base: Path) -> object:
    """git restore [--source=<ref>] [--staged] -- <paths> / git restore <paths>.

    A restore overwrites the listed pathspecs in the working tree. ``--source
    <ref>`` carries a value, so its operand is consumed, not treated as a
    path."""
    return [_resolve(p, base) for p in _operands_skipping_flags(
        rest, {"--staged", "--worktree", "--source", "-s", "-S", "-p",
               "--patch", "-q", "--quiet", "--porcelain", "-W", "--ours",
               "--theirs", "-m", "--merge"})]


def _operands_skipping_flags(rest: list[str], value_flags: set[str]) -> list[str]:
    """The non-option operands of *rest*, skipping every flag and, for flags
    in *value_flags* that take a separate value, that value too. ``--`` stops
    flag interpretation and DISCARDS every operand before it (for checkout,
    that leading operand is the branch/commit name, not a writable path);
    everything after ``--`` is a path. Without ``--``, a leading operand is a
    path only if it is not the value of a preceding value-flag."""
    ops: list[str] = []
    is_value = False
    seen_dd = False
    for a in rest:
        if a == "--":
            seen_dd = True
            ops = []                 # drop branch/option operands seen so far
            is_value = False
            continue
        if seen_dd:
            ops.append(a)
            continue
        if is_value:
            is_value = False         # consume the value of the previous flag
            continue
        if a.startswith("-"):
            if "=" in a:
                continue              # --flag=value: value bundled, already consumed
            if a in value_flags:
                is_value = True       # the NEXT operand is its value
            continue
        ops.append(a)
    return ops


def _targets_patch(rest: list[str], base: Path) -> object:
    """GNU ``patch [-pN] [-d DIR] [-i PATCHFILE] [orig [patchfile]]``.

    Writes the files named by the patch headers. The patch is readable only
    when supplied as a named argument: a ``-i file`` or the second positional
    (patchfile). Otherwise the patch is on stdin -> UNKNOWABLE."""
    strip = 1
    directory = None
    patch_file = None
    original = None
    i = 0
    args = list(rest)
    while i < len(args):
        a = args[i]
        if a == "-p" and i + 1 < len(args):
            i += 1
            try:
                strip = int(args[i])
            except ValueError:
                strip = 1
        elif a.startswith("-p"):
            try:
                strip = int(a[2:])
            except ValueError:
                strip = 1
        elif a in ("-d", "--directory") and i + 1 < len(args):
            i += 1
            directory = args[i]
        elif a.startswith("-d") and len(a) > 2:
            directory = a[2:]
        elif a in ("-i", "--input") and i + 1 < len(args):
            i += 1
            patch_file = args[i]
        elif a.startswith("-i") and len(a) > 2:
            patch_file = a[2:]
        elif a == "--dry-run":
            return []
        elif not a.startswith("-"):
            if original is None:
                original = a
            else:
                patch_file = a    # second positional is the patchfile
        i += 1
    if patch_file is None:
        # Patch read from stdin (no -i, no patchfile positional).
        return UNKNOWABLE
    result = _targets_from_gnupatch(
        patch_file, strip=strip, directory=directory)
    if result is UNKNOWABLE:
        return UNKNOWABLE
    targets = list(result)
    # GNU patch additionally writes `orig` when a raw patch has no headers.
    if original is not None and not targets:
        targets.append(_apply_prefix_strip(original, strip, directory))
    return [_resolve(p, base) for p in targets]


def _apply_prefix_strip(path: str, strip: int,
                        prefix: str | None) -> str:
    """Apply ``-p<strip>`` (drop *strip* leading components) and an optional
    ``--directory``/``-d`` prefix to a bare header/operand path."""
    parts = path.rstrip("/").split("/")
    stripped = "/".join(parts[strip:]) if strip > 0 and len(parts) > strip \
        else path
    return f"{prefix.rstrip('/')}/{stripped.lstrip('/')}" \
        if prefix else stripped


# --- patch parsing --------------------------------------------------------
#
# The two file-write formats that matter:
#   * git apply:  ``diff --git a/X b/Y`` headers; the file touched is the
#     ``b/`` name for a modify/add and the ``a/`` name for a delete/rename.
#   * GNU patch:  ``*** X``/``--- X`` (before) and ``+++ Y`` (after) headers,
#     plus an ``Index: X`` line for svn-style patches.
# Both are collected with a path STRIP (-pN) and PREFIX (--directory) applied.
# Every path a hunk could touch is listed; over-listing only over-refuses.

_GIT_DIFF_HEADER = re.compile(r"^diff --git a/(.*) b/(.*?)\s*$")
_UNIFIED_HEADER = re.compile(r"^[+-]{3}\s+(.+?)(?:\t.*)?\s*$")
_GNUPATCH_AFTER = re.compile(r"^\+\+\+\s+(.+?)(?:\t.*)?\s*$")
_INDEX_LINE = re.compile(r"^Index:\s+(.*?)(?:\t.*)?\s*$")


def _read_patch(patch_file: str) -> list[str] | None:
    """The lines of a patch file, or None if it cannot be read (the caller
    must then refuse: 'could not check what it would write, so it may not
    run'). Reading is for the filesystem only; the texts are not executed."""
    try:
        with open(os.path.expanduser(patch_file), "r", encoding="utf-8",
                  errors="replace") as fh:
            return fh.readlines()
    except (OSError, ValueError):
        return None


def _strip_component(path: str, strip: int) -> str:
    """Remove the first *strip* path components from a header path. Git apply
    default ``-p1`` removes the single ``a/``/``b/`` prefix it writes, so the
    header path is passed WITH that marker and strip drops it."""
    if strip <= 0:
        return path
    parts = path.split("/")
    if len(parts) > strip:
        return "/".join(parts[strip:])
    return path


def _targets_from_patch(patch_file: str, *, strip: int,
                        prefix: str | None) -> object:
    """The paths a ``git apply`` of *patch_file* would add/change/delete.

    Returns UNKNOWABLE when the patch file cannot be read (so we cannot name
    what it would write), else the de-duplicated path list with ``-p``/``--directory``
    applied. The ``a/``/``b/`` marker git writes is kept as the pattable first
    component, so ``-p1`` drops it and ``-p0`` keeps it. Both names are listed:
    a delete removes the ``a/`` name, a modify/add writes the ``b/`` name, and
    a rename touches both -- gating them all is the safe direction."""
    lines = _read_patch(patch_file)
    if lines is None:
        return UNKNOWABLE
    targets: list[str] = []
    for line in lines:
        m = _GIT_DIFF_HEADER.match(line)
        if m:
            a, b = m.group(1), m.group(2)
            targets.append(_strip_component(f"a/{a}", strip))
            if b and b != "/dev/null":
                targets.append(_strip_component(f"b/{b}", strip))
            continue
        m = _UNIFIED_HEADER.match(line)
        if m and "+" in line[:3]:
            p = m.group(1)
            if p and p != "/dev/null":
                targets.append(_strip_component(p, strip))
    if prefix:
        targets = [
            f"{prefix.rstrip('/')}/{t.lstrip('/')}" for t in targets]
    return _dedup([t for t in targets if t])


def _targets_from_gnupatch(patch_file: str, *, strip: int,
                           directory: str | None) -> object:
    """The paths a GNU ``patch`` of *patch_file* would write, from the
    ``+++`` / ``*`` / ``Index:`` headers. UNKNOWABLE when unreadable."""
    lines = _read_patch(patch_file)
    if lines is None:
        return UNKNOWABLE
    targets: list[str] = []
    for line in lines:
        m = _GNUPATCH_AFTER.match(line)
        if m:
            p = m.group(1)
            if p and p != "/dev/null":
                targets.append(_strip_component(p, strip))
            continue
        m = _INDEX_LINE.match(line)
        if m:
            p = m.group(1)
            if p:
                targets.append(_strip_component(p, strip))
    if directory:
        targets = [
            f"{directory.rstrip('/')}/{t.lstrip('/')}" for t in targets]
    return _dedup([t for t in targets if t])
