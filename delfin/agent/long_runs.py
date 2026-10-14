"""Long foreground runs: a whole test suite, and the same command respelled.

Input: a shell command line (or run_tests' ``target`` and ``pytest_args``)
and the directory it runs in. Output:

* ``strip_command_prefixes(cmd)`` -- the command without the words that
  change how long it is waited for but not what it does (``timeout N``,
  ``nice``, ``nohup``, ``env``/``NAME=value`` assignments, ``set -o
  pipefail;``, a leading ``cd X &&``) and without the output trimming
  after it (``2>&1``, ``| tail -N``, ``| head -N``, ``; echo $?``). Two
  spellings with the same result are the same command for a repeat check.
* ``suite_files_in_command(cmd, cwd)`` / ``suite_files_in_run(...)`` --
  how many test files a pytest run covers when it is a run over whole
  directories with no selection (``-k``, ``-m``, ``--lf``, a file or a
  node id); 0 for anything else. At or above ``LARGE_SUITE_MIN_FILES``
  the run is long enough that it is handed to a background job instead
  of holding the turn.

Measured on the shared bug archive, 2026-09-16 .. 10-14: 38 whole-suite
invocations in 8 sessions, at least 144 minutes of blocking tool time; one
session re-ran the same suite five times with ``timeout`` 540/580/560/550,
each cut at the 600 s ceiling.

Rejected: a wall-clock estimate from earlier runs (no record exists for a
first run, which is the case that blocked); counting collected tests
(collection itself took minutes on the suites that blocked). A file count
is bounded (the walk stops at the threshold) and costs < 50 ms here.
"""

from __future__ import annotations

import os
import re
import shlex
from pathlib import Path

__all__ = [
    "LARGE_SUITE_MIN_FILES", "strip_command_prefixes",
    "suite_files_in_command", "suite_files_in_run",
]

#: Test files at or above which a whole-directory run is long. A project
#: with a handful of test files runs them in seconds and keeps the
#: foreground; DELFIN's own suite has over two thousand.
LARGE_SUITE_MIN_FILES = 40

#: Directory entries the count may visit before it gives up (answers 0).
_WALK_BUDGET = 20_000

_SKIP_DIRS = frozenset({
    ".git", ".hg", ".svn", "node_modules", "__pycache__", ".venv", "venv",
    "env", ".tox", ".nox", ".mypy_cache", ".pytest_cache", "build", "dist",
    "site-packages", ".eggs",
})

_OPERATORS = frozenset({"&&", "||", ";", "|", "&", "|&", ";;"})
_ASSIGNMENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")
_DURATION = re.compile(r"^\d+(?:\.\d+)?[smhd]?$")

#: pytest options that take the next word as their value.
_VALUE_OPTIONS = frozenset({
    "-p", "-c", "-o", "-W", "-n", "-r", "--deselect", "--ignore",
    "--ignore-glob", "--tb", "--rootdir", "--confcutdir", "--basetemp",
    "--junitxml", "--junit-xml", "--maxfail", "--durations", "--timeout",
    "--cov", "--cov-report", "--log-level", "--import-mode", "--dist",
    "--color", "--capture", "--report-log", "--json-report-file",
})
#: Options that select a part of the suite, or run none of it.
_SELECTING = frozenset({
    "-k", "-m", "--lf", "--last-failed", "--co", "--collect-only",
    "--fixtures", "--markers", "--version", "-h", "--help", "--setup-plan",
    "--setup-only", "--trace-config",
})


def _tokens(cmd: str) -> list[str]:
    text = str(cmd or "").replace("\n", " ; ")
    try:
        lex = shlex.shlex(text, posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        return list(lex)
    except ValueError:
        return text.split()


def _segments(cmd: str) -> list[list[str]]:
    """Simple commands of ``cmd``, operators removed, redirections
    dropped with their targets."""
    out: list[list[str]] = []
    seg: list[str] = []
    toks = _tokens(cmd)
    i = 0
    while i < len(toks):
        tok = toks[i]
        if tok in _OPERATORS:
            if seg:
                out.append(seg)
            seg = []
            out.append([tok])
            i += 1
            continue
        if tok and set(tok) <= set("<>&"):
            # A redirection: drop it and its target; a file descriptor
            # number written before it (`2>&1`) went in as its own token.
            if seg and seg[-1].isdigit():
                seg.pop()
            i += 2
            continue
        seg.append(tok)
        i += 1
    if seg:
        out.append(seg)
    return out


def _strip_words(words: list[str]) -> list[str]:
    """``words`` without the leading wrappers that only change waiting."""
    i, n = 0, len(words)
    while i < n:
        w = words[i]
        if _ASSIGNMENT.match(w):
            i += 1
        elif w == "env":
            i += 1
            while i < n and (words[i].startswith("-") or _ASSIGNMENT.match(words[i])):
                i += 2 if words[i] in ("-u", "--unset", "-C", "--chdir") else 1
        elif w == "timeout":
            i += 1
            while i < n and words[i].startswith("-"):
                i += 2 if words[i] in ("-s", "-k", "--signal", "--kill-after") else 1
            if i < n and _DURATION.match(words[i]):
                i += 1
        elif w == "nice":
            i += 1
            if i < n and words[i] == "-n":
                i += 2
            elif i < n and re.match(r"^(?:-\d+|-n\d+|--adjustment=\S+)$", words[i]):
                i += 1
        elif w in ("ionice", "stdbuf"):
            i += 1
            while i < n and words[i].startswith("-"):
                i += 2 if words[i] in ("-c", "-n", "-p", "-i", "-o", "-e") else 1
        elif w in ("nohup", "time", "command", "exec", "builtin"):
            i += 1
        else:
            break
    return words[i:]


def _is_noise(words: list[str]) -> bool:
    """A segment that changes nothing about what a command computes."""
    if not words:
        return True
    head = words[0]
    if head in ("cd", "pushd", "popd", "set", "true", ":"):
        return True
    if head in ("tail", "head"):
        return True
    if head == "echo" and any(re.search(r"\$\?|PIPESTATUS|\brc\b|EXIT", w)
                              for w in words[1:]):
        return True
    return False


def strip_command_prefixes(cmd: str) -> str:
    """The command a repeat is compared on (see the module docstring).

    Whitespace is collapsed; operators are kept between the segments that
    remain, as ``;``. An empty result falls back to the collapsed text, so
    a command made only of noise still has a key of its own."""
    parts: list[str] = []
    for seg in _segments(cmd):
        if len(seg) == 1 and seg[0] in _OPERATORS:
            continue
        words = _strip_words(seg)
        if _is_noise(words):
            continue
        parts.append(" ".join(words))
    key = " ; ".join(parts).strip()
    return key or " ".join(str(cmd or "").split())


def _pytest_args(words: list[str]):
    """The arguments after a pytest invocation in ``words``, or None."""
    if not words:
        return None
    head = os.path.basename(words[0])
    if head in ("pytest", "py.test"):
        return words[1:]
    if re.match(r"^python(?:\d+(?:\.\d+)?)?$", head):
        i = 1
        while i < len(words) and words[i].startswith("-") and words[i] != "-m":
            i += 1
        if (i + 1 < len(words) and words[i] == "-m"
                and words[i + 1] in ("pytest", "py.test")):
            return words[i + 2:]
    return None


def _positionals(args: list[str]):
    """Positional paths of a pytest argument list, or None when an option
    selects part of the suite."""
    out: list[str] = []
    i = 0
    while i < len(args):
        a = args[i]
        name = a.split("=", 1)[0]
        if name in _SELECTING or (a[:2] in ("-k", "-m") and len(a) > 2
                                  and not a.startswith("--")):
            return None
        if a in _VALUE_OPTIONS:
            i += 2
            continue
        if a.startswith("-") or a.isdigit():
            i += 1
            continue
        out.append(a)
        i += 1
    return out


def _count_test_files(root: Path, limit: int) -> int:
    """Test files under ``root``, counting stops at ``limit``. A walk
    that exceeds its entry budget answers 0 (not long, not refused)."""
    seen = 0
    count = 0
    stack = [root]
    while stack:
        d = stack.pop()
        try:
            entries = list(os.scandir(d))
        except OSError:
            continue
        for e in entries:
            seen += 1
            if seen > _WALK_BUDGET:
                return 0
            try:
                if e.is_dir(follow_symlinks=False):
                    if e.name not in _SKIP_DIRS and not e.name.startswith("."):
                        stack.append(Path(e.path))
                elif e.name.endswith(".py") and (
                        e.name.startswith("test_") or e.name.endswith("_test.py")):
                    count += 1
                    if count >= limit:
                        return count
            except OSError:
                continue
    return count


def _files_covered(paths: list[str], cwd) -> int:
    """Test files a run over ``paths`` covers; 0 when any path is a file
    or a node id (a chosen part, not the suite)."""
    base = Path(cwd)
    if not paths:
        tests = base / "tests"
        roots = [tests if tests.is_dir() else base]
    else:
        roots = []
        for p in paths:
            if "::" in p:
                return 0
            cand = Path(p) if os.path.isabs(p) else base / p
            if not cand.is_dir():
                return 0
            roots.append(cand)
    total = 0
    for r in roots:
        total += _count_test_files(r, LARGE_SUITE_MIN_FILES - total)
        if total >= LARGE_SUITE_MIN_FILES:
            break
    return total


def suite_files_in_command(cmd: str, cwd) -> int:
    """Test files covered by a whole-directory pytest run in ``cmd``
    when that is at least ``LARGE_SUITE_MIN_FILES``; else 0. Never raises."""
    try:
        for seg in _segments(cmd):
            args = _pytest_args(_strip_words(seg))
            if args is None:
                continue
            paths = _positionals(args)
            if paths is None:
                continue
            n = _files_covered(paths, cwd)
            if n >= LARGE_SUITE_MIN_FILES:
                return n
    except Exception:
        return 0
    return 0


def suite_files_in_run(target: str, pytest_args, cwd) -> int:
    """The same for a run_tests call (``target`` plus ``pytest_args``)."""
    try:
        args = [str(a) for a in (pytest_args or [])]
        paths = _positionals(args)
        if paths is None:
            return 0
        if str(target or "").strip():
            paths = [str(target).strip()] + paths
        n = _files_covered(paths, cwd)
        return n if n >= LARGE_SUITE_MIN_FILES else 0
    except Exception:
        return 0
