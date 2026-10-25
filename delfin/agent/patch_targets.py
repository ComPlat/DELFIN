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


#: git subcommands that rewrite work-tree files without naming them. With
#: the tree moved by --work-tree / GIT_WORK_TREE / core.worktree the tree
#: itself is the target; ``-C dir`` is named by the caller's destination
#: options already.
_TREE_WIDE_WRITES = frozenset({
    "reset", "stash", "merge", "pull", "rebase", "cherry-pick", "revert",
    "switch", "clean", "rm", "mv", "read-tree", "checkout-index",
    "sparse-checkout", "submodule", "bisect",
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
    argv, stdin = _drop_redirections(argv)
    env: dict[str, str] = {}
    argv, base, xargs = _unwrap(argv, base, env)
    if not argv:
        return []
    name = Path(argv[0]).name
    rest = argv[1:]
    if xargs and name in ("git", "patch", "git-apply") and (
            name != "git" or _git_subcommand(rest) in _WORKTREE_WRITE_SUBCOMMANDS):
        # xargs appends operands read from stdin: the patch files or
        # pathspecs are not in the text.
        return UNKNOWABLE

    if name == "git":
        return _targets_git(rest, base, stdin=stdin,
                            env_tree=env.get("GIT_WORK_TREE"))
    if name == "patch":
        return _targets_patch(rest, base, stdin=stdin)
    if name == "git-apply":
        # `git-apply` is the core built as a standalone binary.
        return _targets_git_apply(rest, base, stdin=stdin)
    return []


#: A redirection word as shlex leaves it: ``2>&1``, ``2>/dev/null``, ``>``,
#: ``<a.diff``, ``<<EOF``. Group 1 is the operator, group 2 an attached
#: operand.
_REDIRECTION = re.compile(r"^(?:\d*|&)(<<<|<<-?|>>|>&|>\||<&|<>|>|<)(.*)$")


def _drop_redirections(argv: list[str]) -> tuple[list[str], object]:
    """*argv* without its redirection words, and what stdin is.

    Returns ``(argv, stdin)``: ``stdin`` is the file a ``< file`` names,
    False for a here-document or here-string (a patch in the command text
    that this module does not read), else None. Without this ``git apply
    a.diff 2>&1`` read ``2>`` as a second patch file and was refused."""
    out: list[str] = []
    stdin: object = None
    i = 0
    while i < len(argv):
        m = _REDIRECTION.match(argv[i])
        if not m:
            out.append(argv[i])
            i += 1
            continue
        fd = argv[i][:m.start(1)]
        op, operand = m.group(1), m.group(2)
        if not operand and i + 1 < len(argv):
            operand = argv[i + 1]
            i += 1
        if fd in ("", "0"):
            if op.startswith("<<"):
                stdin = False
            elif op in ("<", "<>"):
                stdin = operand
        i += 1
    return out, stdin


#: Commands that run the rest of their argument list as a command, with the
#: options of each that take a separate value.
_WRAPPERS: dict[str, frozenset] = {
    "env": frozenset({"-u", "--unset", "-C", "--chdir", "-S",
                      "--split-string"}),
    "command": frozenset(),
    "builtin": frozenset(),
    "exec": frozenset({"-a"}),
    "nice": frozenset({"-n", "--adjustment"}),
    "nohup": frozenset(),
    "time": frozenset({"-f", "--format", "-o", "--output"}),
    "timeout": frozenset({"-s", "--signal", "-k", "--kill-after"}),
    "stdbuf": frozenset({"-i", "-o", "-e", "--input", "--output",
                         "--error"}),
    "ionice": frozenset({"-c", "--class", "-n", "--classdata", "-p",
                         "-P", "-u"}),
    "chrt": frozenset(),
    "taskset": frozenset(),
    "sudo": frozenset({"-u", "--user", "-g", "--group", "-h", "--host",
                       "-p", "--prompt", "-C", "--close-from", "-r",
                       "--role", "-t", "--type", "-U", "--other-user",
                       "-T", "--command-timeout", "-D", "--chdir"}),
    "doas": frozenset({"-u", "-C"}),
    "xargs": frozenset({"-a", "--arg-file", "-d", "--delimiter", "-E",
                        "-e", "-I", "-i", "-L", "-l", "-n", "--max-args",
                        "-P", "--max-procs", "-s", "--max-chars"}),
}


def _unwrap(argv: list[str], base: Path,
            env: dict[str, str]) -> tuple[list[str], Path, bool]:
    """Strip ``VAR=value`` words and command wrappers (env, nice, timeout,
    sudo, xargs ... -- :data:`_WRAPPERS`) from the front of *argv*.

    Returns ``(argv, base, xargs)``: the wrapped command, the directory it
    runs in (``env -C dir``, ``sudo -D dir``), and whether xargs feeds it
    operands. Assignments land in *env*. A wrapper's first non-option word
    is the command, except timeout's DURATION and taskset/chrt's leading
    number, which are skipped."""
    xargs = False
    while argv:
        word = argv[0]
        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", word):
            key, _, value = word.partition("=")
            env[key] = value
            argv = argv[1:]
            continue
        name = Path(word).name
        if name not in _WRAPPERS:
            break
        xargs = xargs or name == "xargs"
        takes = _WRAPPERS[name]
        i = 1
        while i < len(argv):
            a = argv[i]
            if a == "--":
                i += 1
                break
            if not a.startswith("-") or a == "-":
                if name == "env" and re.match(
                        r"^[A-Za-z_][A-Za-z0-9_]*=", a):
                    key, _, value = a.partition("=")
                    env[key] = value
                    i += 1
                    continue
                break
            key, eq, value = a.partition("=")
            if key in takes and not eq:
                value = argv[i + 1] if i + 1 < len(argv) else ""
                i += 1
            if key in ("-C", "--chdir", "-D") and name in ("env", "sudo") \
                    and value:
                base = Path(_resolve(value, base))
            i += 1
        if name == "timeout" and i < len(argv):
            i += 1               # DURATION
        if name in ("taskset", "chrt") and i < len(argv):
            i += 1               # mask / priority
        argv = argv[i:]
    return argv, base, xargs


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


def _git_tree(rest: list[str], base: Path,
              env_tree: str | None) -> tuple[Path, bool]:
    """The directory ``git``'s path operands resolve against.

    Every leading ``-C dir`` moves the cwd, each relative to the last
    (``-C a -C b`` is ``a/b``). The work tree is then, first match wins,
    ``--work-tree[=]dir``, ``GIT_WORK_TREE`` (*env_tree*) or ``-c
    core.worktree=dir``, relative to that cwd; without one it is the cwd.
    ``--git-dir`` moves only ``.git`` and is ignored. A ``core.worktree``
    set in a config file is not in the text and is not seen. The second
    value is True when one of the three named the tree."""
    cwd = base
    tree = None
    config_tree = None
    i = 0
    while i < len(rest):
        arg = rest[i]
        if not arg.startswith("-"):
            break                    # the subcommand: globals end here
        nxt = rest[i + 1] if i + 1 < len(rest) else ""
        if arg == "-C":
            if nxt:
                cwd = Path(_resolve(nxt, cwd))
            i += 2
            continue
        if arg == "--work-tree":
            tree = nxt
            i += 2
            continue
        if arg.startswith("--work-tree="):
            tree = arg.split("=", 1)[1]
        elif arg == "-c":
            key, _, value = nxt.partition("=")
            if key.lower() == "core.worktree":
                config_tree = value
            i += 2
            continue
        elif arg in ("--git-dir", "--namespace", "--super-prefix",
                     "--config-env"):
            i += 2
            continue
        i += 1
    chosen = tree or env_tree or config_tree
    return (Path(_resolve(chosen, cwd)), True) if chosen else (cwd, False)


def _targets_git(rest: list[str], base: Path, *, stdin: object = None,
                 env_tree: str | None = None) -> object:
    """Write targets of a ``git <sub> …`` command. ``rest`` is the argument
    list after ``git``; *stdin* is the file a ``< file`` names (see
    :func:`_drop_redirections`); *env_tree* is a ``GIT_WORK_TREE`` set on
    the command line."""
    try:
        sub = _git_subcommand(rest)
    except Exception:
        return []
    eff_base, moved = _git_tree(rest, base, env_tree)
    after = _after_subcommand(rest, sub)
    if sub == "apply":
        return _targets_git_apply(after, eff_base, stdin=stdin)
    if sub == "am":
        # A mailbox: the patch lives inside the mbox body, not on the command
        # line. Reading it robustly is not worth the risk of a missed target;
        # refuse the command outright.
        return UNKNOWABLE
    if sub == "checkout":
        named = _targets_git_checkout(after, eff_base)
        # A branch switch names no file but rewrites the moved tree.
        return [str(eff_base)] if moved and named == [] else named
    if sub == "restore":
        return _targets_git_restore(after, eff_base)
    if moved and sub in _TREE_WIDE_WRITES:
        # The files are not in the text; the tree they are in is.
        return [str(eff_base)]
    return []


def _after_subcommand(rest: list[str], sub: str) -> list[str]:
    """The argument list after the *sub* command token in *rest*."""
    try:
        idx = rest.index(sub)
    except ValueError:
        return []
    return rest[idx + 1:]


def _targets_git_apply(rest: list[str], base: Path, *,
                       stdin: object = None) -> object:
    """git apply [-pN] [--directory=prefix] [-R] [--include=…] [patch ...].

    The files it writes are read from EVERY patch file named on the command
    line (``git apply a.diff b.diff`` applies both). With no file, or ``-``,
    the patch is stdin: the file a ``< file`` redirection names is read, and
    a pipe, here-document or terminal is UNKNOWABLE -> the caller refuses.

    ``--3way`` and ``--index`` apply the patch's hunks to the index AND
    write the working tree, so the patch's files are named. ``--cached``
    stages to the index alone and writes no working-tree file, so it names
    nothing (``[]``). ``--check`` validates without writing -> nothing."""
    strip = 1
    prefix = None
    check_only = False
    cached_only = False
    patch_files: list[str] = []

    i = 0
    while i < len(rest):
        a = rest[i]
        if a == "--check":
            check_only = True
        elif a == "--cached":
            cached_only = True  # index only; nothing written to the work tree
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
        elif a.startswith("-p") and not a.startswith("--"):
            try:
                strip = int(a[2:])
            except ValueError:
                strip = 1
        elif a in ("--include", "--exclude", "-C", "--whitespace") \
                and i + 1 < len(rest):
            i += 1              # the option's value, not a patch file
        elif a == "-":
            patch_files.append(a)
        elif not a.startswith("-"):
            patch_files.append(a)   # every positional is a patch file
        i += 1
    if check_only or cached_only:
        return []
    files = [f for f in patch_files if f != "-"]
    if "-" in patch_files or not files:
        if not isinstance(stdin, str):
            # Patch on stdin from a pipe, here-document or terminal.
            return UNKNOWABLE
        files.append(stdin)
    found: list[str] = []
    for pf in files:
        result = _targets_from_patch(
            _resolve(pf, base), strip=strip, prefix=prefix)
        if result is UNKNOWABLE:
            return UNKNOWABLE
        found.extend(result)
    return _dedup([_resolve(p, base) for p in found])


def _targets_git_checkout(rest: list[str], base: Path) -> object:
    """git checkout [<tree-ish>] [--] <paths>  /  git checkout <branch>.

    After ``--`` every operand is a path. Without it, with two or more
    operands git reads the first as a tree-ish and the rest as paths;
    a single operand is a path when it is path-shaped (``/``, glob,
    pathspec magic) or names an existing file under *base*, else a branch
    -- ``git checkout main`` writes no file ``main``. With ``-b``/``-B``/
    ``--orphan``/``--detach`` the operand is a start point, never a path.
    ``--pathspec-from-file`` names its paths outside the text ->
    UNKNOWABLE. Of the flags only ``-b``, ``-B``, ``--orphan`` and
    ``--conflict``/``-U`` (when not attached) take a separate value;
    ``-f``, ``-p``, ``--ours`` and the rest stand alone."""
    if any(a.startswith("--pathspec-from-file") for a in rest):
        return UNKNOWABLE
    value_flags = {"-b", "-B", "--orphan", "--conflict", "-U"}
    ops = _operands_skipping_flags(rest, value_flags)
    if "--" in rest:
        return _dedup([_resolve(p, base) for p in ops])
    if any(a in ("-b", "-B", "--orphan", "--detach") for a in rest):
        return []
    if len(ops) >= 2:
        paths = ops[1:]
        if _is_path_shaped(ops[0]) or (base / ops[0]).exists():
            paths = ops
    else:
        # GT1: a single operand that is path-shaped but names an EXISTING
        # LOCAL BRANCH of the repository is a branch switch, not a
        # pathspec. Branch names may contain "/" (``agent/s3-ds1b``), so
        # ``_is_path_shaped`` alone cannot separate the two readings --
        # wave-21 measured it: the write gate refused
        # ``git checkout agent/s3-ds1b`` as a write to
        # ``<worktree>/agent/s3-ds1b``. The branch reading wins ONLY when
        # git's own ref store confirms the branch (``git show-ref --verify
        # refs/heads/<name>``) AND no file of that name exists under
        # *base* -- the file reading wins the ambiguity (git checks the
        # index and the tree before the refs), so nothing here is looser
        # than before: every operand that is not a confirmed local branch
        # is classified exactly as it was. Fail-closed: a git probe that
        # errors or times out leaves the operand a path, the direction
        # that refuses.
        paths = [p for p in ops if _is_path_shaped(p) or (base / p).exists()]
        if paths and len(ops) == 1 and _is_path_shaped(ops[0]) \
                and not (base / ops[0]).exists() \
                and _branch_exists_locally(ops[0], base):
            paths = []
    return _dedup([_resolve(p, base) for p in paths])


def _targets_git_restore(rest: list[str], base: Path) -> object:
    """git restore [-s <tree>] [-S] [-W] [--] <paths>.

    A restore writes the work tree unless ``--staged``/``-S`` is given
    without ``--worktree``/``-W`` (index only -> ``[]``). Short options
    bundle (``-SW``); ``-s``/``--source`` is the one that takes a separate
    value. ``--pathspec-from-file`` -> UNKNOWABLE."""
    if any(a.startswith("--pathspec-from-file") for a in rest):
        return UNKNOWABLE
    staged = worktree = False
    for a in rest:
        if a == "--":
            break
        if a == "--staged":
            staged = True
        elif a == "--worktree":
            worktree = True
        elif a.startswith("-") and not a.startswith("--"):
            letters = a[1:]
            if "s" in letters:
                letters = letters[:letters.index("s")]   # -sHEAD: value
            staged = staged or "S" in letters
            worktree = worktree or "W" in letters
    if staged and not worktree:
        return []
    return _dedup([_resolve(p, base) for p in _operands_skipping_flags(
        rest, {"-s", "--source", "--conflict"})])


def _operands_skipping_flags(rest: list[str], value_flags: set[str]) -> list[str]:
    """The non-option operands of *rest*, skipping every flag and, for flags
    in *value_flags* that take a separate value, that value too. ``--``
    stops flag interpretation and DISCARDS every operand before it (for
    checkout, that leading operand is the branch/commit name, not a
    writable path); everything after ``--`` is a path."""
    ops: list[str] = []
    is_value = False
    seen_dd = False
    for a in rest:
        if a == "--" and not seen_dd:
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


def _is_path_shaped(p: str) -> bool:
    """True when *p* looks like a git pathspec rather than a branch name: it
    carries a ``/``, a glob (``*?[]``) or a pathspec magic (``:(top)`` …)."""
    return ("/" in p or ":" in p or any(c in p for c in "*?[]"))


def _branch_exists_locally(name: str, base: Path) -> bool:
    """GT1: does ``name`` resolve to an EXISTING LOCAL branch of the
    repository *base* belongs to?

    Verified against git's own ref store, never against the tree, and
    ALWAYS in the repository of *base* (wave-22 operator finding: a
    probe in the process cwd answered for the wrong repo and made a
    path read as a branch). The call goes through
    :func:`delfin.agent.safe_git.run` -- every git command DELFIN runs
    must; a direct ``subprocess`` git call breaks the safe_git guard
    test. Fails closed -- any refusal, error, timeout or non-zero
    return means "no branch confirmed", and the caller keeps the
    operand classified as a path. Reads nothing outside the
    repository; the only I/O is the probe itself, which costs a
    measured ~5 ms and runs only for a single path-shaped checkout
    operand that is not an existing file."""
    try:
        from delfin.agent import safe_git
        probe = safe_git.run(
            ["show-ref", "--verify", "--quiet", "refs/heads/" + name],
            base, timeout=5)
        return probe.returncode == 0
    except Exception:
        return False


#: GNU patch options that take a separate value, by short and long name.
_PATCH_VALUE_OPTS = {
    "-p": "strip", "--strip": "strip",
    "-d": "directory", "--directory": "directory",
    "-i": "input", "--input": "input",
    "-o": "output", "--output": "output",
    "-r": "reject", "--reject-file": "reject",
    "-B": "prefix", "--prefix": "prefix",
    "-Y": "basename_prefix", "--basename-prefix": "basename_prefix",
    "-z": "suffix", "--suffix": "suffix",
    "-D": "", "--ifdef": "", "-F": "", "--fuzz": "",
    "-V": "", "--version-control": "", "-g": "", "--get": "",
    "--quoting-style": "",
}


def _targets_patch(rest: list[str], base: Path, *,
                   stdin: object = None) -> object:
    """GNU ``patch [options] [orig [patchfile]]``.

    Writes the files named by the patch headers, plus the side files its
    options name: ``-r``/``--reject-file`` (the reject file), ``-B``/
    ``--prefix`` (prefix + each target), ``-Y``/``--basename-prefix`` and
    ``-z``/``--suffix`` (backups beside each target). ``-o`` is named by
    the caller's destination options. The patch is ``-i``/``--input``, the
    second positional, or the file a ``< file`` names; a pipe or terminal
    on stdin is UNKNOWABLE. ``-d``/``--directory`` applies before anything
    is read."""
    opts: dict[str, str] = {}
    positionals: list[str] = []
    i = 0
    args = list(rest)
    while i < len(args):
        a = args[i]
        if a == "--dry-run":
            return []
        if a == "--":
            positionals.extend(args[i + 1:])
            break
        if a.startswith("--"):
            key, eq, value = a.partition("=")
            if key in _PATCH_VALUE_OPTS:
                if not eq and i + 1 < len(args):
                    i += 1
                    value = args[i]
                if _PATCH_VALUE_OPTS[key]:
                    opts[_PATCH_VALUE_OPTS[key]] = value
        elif a.startswith("-") and len(a) > 1:
            key = a[:2]
            if key in _PATCH_VALUE_OPTS:
                value = a[2:]
                if not value and i + 1 < len(args):
                    i += 1
                    value = args[i]
                if _PATCH_VALUE_OPTS[key]:
                    opts[_PATCH_VALUE_OPTS[key]] = value
        else:
            positionals.append(a)
        i += 1
    try:
        strip = int(opts.get("strip", "1"))
    except ValueError:
        strip = 1
    original = positionals[0] if positionals else None
    patch_file = opts.get("input")
    if patch_file is None and len(positionals) >= 2:
        patch_file = positionals[1]
    if patch_file is None or patch_file == "-":
        if not isinstance(stdin, str):
            # Patch read from a pipe, here-document or terminal.
            return UNKNOWABLE
        patch_file = stdin
    directory = opts.get("directory")
    # GNU patch changes into -d DIR before it reads anything.
    work = Path(_resolve(directory, base)) if directory else base
    result = _targets_from_gnupatch(
        _resolve(patch_file, work), strip=strip, directory=None)
    if result is UNKNOWABLE:
        return UNKNOWABLE
    found = list(result)
    if original is not None:
        # An explicit file operand is the one file patched.
        found.append(original)
    found = [_resolve(p, work) for p in found]
    side: list[str] = []
    if opts.get("reject"):
        side.append(_resolve(opts["reject"], work))
    for t in found:
        if opts.get("prefix"):
            side.append(_resolve(opts["prefix"] + os.path.relpath(t, work)
                                 if t.startswith(str(work) + os.sep)
                                 else opts["prefix"] + t, work))
        if opts.get("basename_prefix"):
            side.append(os.path.join(os.path.dirname(t),
                                     opts["basename_prefix"]
                                     + os.path.basename(t)))
        if opts.get("suffix"):
            side.append(t + opts["suffix"])
    return _dedup(found + side)


# --- patch parsing --------------------------------------------------------
#
# The two file-write formats that matter:
#   * git apply:  ``diff --git a/X b/Y`` headers, ``rename to``/``copy to``
#     lines, and a traditional ``--- X`` / ``+++ Y`` pair (``+++ /dev/null``
#     is a delete of X).
#   * GNU patch:  ``--- X``/``+++ Y`` (unified), ``*** X``/``--- Y``
#     (context), and an ``Index: X`` line; patch picks among those names,
#     so every one is listed.
# Both are collected with a path STRIP (-pN) and PREFIX (--directory) applied.
# Every path a hunk could touch is listed; over-listing only over-refuses.

_GIT_DIFF_HEADER = re.compile(r"^diff --git a/(.*) b/(.*?)\s*$")
_GIT_NAME_LINE = re.compile(r"^(?:rename|copy) (?:from|to) (.+?)\s*$")
_OLD_HEADER = re.compile(r"^---\s+(.+?)(?:\t.*)?\s*$")
_NEW_HEADER = re.compile(r"^\+\+\+\s+(.+?)(?:\t.*)?\s*$")
_CONTEXT_OLD = re.compile(r"^\*\*\*\s+(.+?)(?:\t.*)?\s*$")
_INDEX_LINE = re.compile(r"^Index:\s+(.*?)(?:\t.*)?\s*$")
#: A context-diff hunk range (``*** 1,3 ****``/``--- 1,3 ----``), not a name.
_HUNK_RANGE = re.compile(r"^\d+(?:,\d+)?\s+(?:\*{4}|-{4})$")

#: Largest patch file read. A bigger one is UNKNOWABLE (refused), so the
#: gate's own read is bounded in memory and time.
_MAX_PATCH_BYTES = 16 * 1024 * 1024


def _read_patch(patch_file: str) -> list[str] | None:
    """The lines of a patch file, or None if it cannot be read (the caller
    must then refuse: 'could not check what it would write, so it may not
    run').

    Only a regular file of at most :data:`_MAX_PATCH_BYTES` is read. It is
    opened non-blocking and checked with ``fstat`` on the open descriptor:
    a FIFO named as the patch blocked this read until a writer appeared,
    and ``/dev/zero`` never ends. Reading is for the filesystem only; the
    text is not executed."""
    import stat
    try:
        fd = os.open(os.path.expanduser(patch_file),
                     os.O_RDONLY | os.O_NONBLOCK | getattr(os, "O_NOCTTY", 0))
    except (OSError, ValueError):
        return None
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode) or st.st_size > _MAX_PATCH_BYTES:
            return None
        chunks: list[bytes] = []
        left = _MAX_PATCH_BYTES + 1
        while left > 0:
            chunk = os.read(fd, min(left, 1 << 20))
            if not chunk:
                break
            chunks.append(chunk)
            left -= len(chunk)
        if left <= 0:
            return None
    except OSError:
        return None
    finally:
        os.close(fd)
    return b"".join(chunks).decode("utf-8", errors="replace").splitlines(True)


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


def _header_names(lines: list[str]) -> tuple[list[str], list[str]]:
    """``(header names, rename/copy names)`` of a patch, ``/dev/null``
    dropped. Header names are every ``diff --git`` side (with its ``a/``/
    ``b/`` marker), ``+++`` name, ``---`` name followed by a ``+++`` line,
    ``***``/``---`` context pair and ``Index:`` name; rename/copy names
    are git's full paths, which ``-p`` does not strip."""
    names: list[str] = []
    full: list[str] = []
    for k, line in enumerate(lines):
        nxt = lines[k + 1] if k + 1 < len(lines) else ""
        m = _GIT_DIFF_HEADER.match(line)
        if m:
            names += [f"a/{m.group(1)}", f"b/{m.group(2)}"]
            continue
        m = _GIT_NAME_LINE.match(line)
        if m:
            full.append(m.group(1))
            continue
        m = _NEW_HEADER.match(line)
        if m:
            names.append(m.group(1))
            continue
        m = _OLD_HEADER.match(line)
        if m and (nxt.startswith("+++ ") or nxt.startswith("+++\t")
                  or lines[k - 1:k] and _CONTEXT_OLD.match(lines[k - 1])):
            names.append(m.group(1))
            continue
        m = _CONTEXT_OLD.match(line)
        if m and _OLD_HEADER.match(nxt):
            names.append(m.group(1))
            continue
        m = _INDEX_LINE.match(line)
        if m:
            names.append(m.group(1))
    keep = [n for n in names
            if n and n != "/dev/null" and not _HUNK_RANGE.match(n)]
    return keep, [n for n in full if n and n != "/dev/null"]


def _targets_from_patch(patch_file: str, *, strip: int,
                        prefix: str | None) -> object:
    """The paths a ``git apply`` of *patch_file* would add/change/delete.

    Returns UNKNOWABLE when the patch file cannot be read, else the
    de-duplicated path list with ``-p``/``--directory`` applied. Both
    names of every header are listed -- a delete removes the old name, a
    modify/add writes the new one, a rename touches both -- and the
    ``rename``/``copy`` lines, which git applies as written."""
    lines = _read_patch(patch_file)
    if lines is None:
        return UNKNOWABLE
    names, full = _header_names(lines)
    found = [_strip_component(n, strip) for n in names] + full
    if prefix:
        found = [f"{prefix.rstrip('/')}/{t.lstrip('/')}" for t in found]
    return _dedup([t for t in found if t])


def _targets_from_gnupatch(patch_file: str, *, strip: int,
                           directory: str | None) -> object:
    """The paths a GNU ``patch`` of *patch_file* would write, from every
    header name (see :func:`_header_names`). UNKNOWABLE when unreadable."""
    lines = _read_patch(patch_file)
    if lines is None:
        return UNKNOWABLE
    names, full = _header_names(lines)
    found = [_strip_component(n, strip) for n in names + full]
    if directory:
        found = [f"{directory.rstrip('/')}/{t.lstrip('/')}" for t in found]
    return _dedup([t for t in found if t])
