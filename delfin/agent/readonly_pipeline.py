"""Recognise a whole shell pipeline/sequence as safely read-only.

The gate auto-allows a compound only when EVERY shell segment is
individually auto-allowed. That is correct but coarse: a harmless
read-only form such as ``grep x f 2>&1 | tail -5`` or ``git log ...;
git show ...`` falls to the confirm gate because a full pipeline is not
a single match. This module answers the question "is this whole command
read-only" for those forms, conservatively:

* every segment (a pipeline ``a | b``, a sequence ``a ; b`` / ``a && b``
  / ``a || b``) must be a known read-only command;
* no redirection to a FILE (``>``, ``>>``, ``2>``) -- ``2>&1`` (merging
  stderr into stdout) is not a file write and stays allowed;
* no command substitution (``$( ... )``, backticks);
* no interpreter or command-runner as the thing that is run (python,
  bash, sh, perl, xargs, eval, env, ...).

The read-only command set is deliberately the intrinsically read-only
core. A command is admitted only when, once the shell-level write
channels are closed, no flag or argument can make it write or run a
program: ``sed``/``awk`` are left OUT (both can write or run code), as
are ``find`` (``-delete``/``-exec``) and ``hostname`` (a positional name
sets the host). Two admitted commands carry a write in a FLAG and are
guarded, not trusted: ``sort -o/--output`` writes a file and
``date -s/--set`` sets the clock, so a line holding that flag is refused.
``git branch``/``git tag`` are admitted only in LIST form because the
same subcommand also carries every write (``-d`` delete, ``-m`` move,
``-a`` annotate, a bare name creates) -- refused otherwise. Conservative
refusals are the direction this module is allowed to err in.
"""

from __future__ import annotations

import re

#: Commands that cannot write to a file or start a program on their own,
#: so they are read-only for any flags/arguments once shell redirection
#: to a file and command substitution are banned. ``env`` is excluded (it
#: runs a program), as are sed/awk/find for the reasons in the module
#: docstring, and ``xargs`` (runs its argument as a command).
#: ``hostname`` is excluded: with a positional argument it SETS the host.
#: ``sort``/``date`` are flagged below (``sort -o`` writes a file;
#: ``date -s`` sets the clock).
_INTRINSIC_READONLY = frozenset({
    "ls", "cat", "head", "tail", "wc", "grep", "egrep", "fgrep",
    "sort", "uniq", "cut", "tr", "file", "stat", "basename", "dirname",
    "realpath", "readlink", "echo", "printf", "seq", "pwd", "which",
    "type", "ps", "pgrep", "date", "whoami", "uname", "id",
    "nl", "fold", "column", "uuidgen",
})

#: Flags that make a specific intrinsic command write (a file or a system
#: side effect), so the command is NOT read-only when one is present.
#: ``sort -o f`` / ``sort --output f`` / ``sort --output=f`` write a file;
#: ``date -s`` / ``date --set`` set the clock. A write in the flag, like
#: a write in a git subcommand, must never yield a read-only True.
_WRITE_FLAG = {
    "sort": frozenset({"-o", "--output"}),
    "date": frozenset({"-s", "--set"}),
}

#: ``git`` subcommands that only read. Anything else (add, commit, push,
#: checkout, reset, merge, ...) is a write or a history change and is
#: refused.
_READONLY_GIT = frozenset({
    "log", "show", "status", "diff", "rev-parse", "blame", "describe",
    "branch", "tag", "ls-files", "ls-tree", "name-rev", "shortlog",
    "count-objects", "reflog",
})

#: ``git`` subcommands that are read-only only in LIST form. ``branch`` and
#: ``tag`` default to a listing query, but the same subcommand ALSO carries
#: every write: ``-d``/``-D`` delete, ``-m``/``--move`` rename, ``-a``/
#: ``-A``/``--annotate`` create/annotate, and a bare positional name
#: creates one. So these are safe only when the rest of the line holds
#: plain flags and none of them writes.
_GIT_LIST_ONLY = frozenset({"branch", "tag"})

#: These flags (and the ones below) make a ``git branch``/``git tag`` line
#: write or history-changing, never read-only.
_GIT_WRITE_FLAGS = frozenset({
    "-d", "-D", "-m", "-M", "--move", "--delete", "-a", "-A", "-f",
    "--force", "-e", "--edit", "--force-create", "--orphan",
    "--rename", "--track",
})

#: A command that RUNS another command (or code), so it is never read-only
#: even when its own name is on no command list here. Additional guard on
#: top of the command whitelist.
_RUNNERS = frozenset({
    "python", "python2", "python3", "bash", "sh", "zsh", "dash", "ksh",
    "perl", "ruby", "php", "node", "env", "xargs", "exec", "eval",
    "source", "nohup", "timeout", "nice", "ssh", "sudo",
})

# Shell option-setting prefix that the operator routinely prepends to a
# gate invocation. It changes no file and starts no program; a command
# that begins with it is judged on what follows.
_SHELL_OPT_PREFIX = re.compile(
    r"^\s*set\s+(?:(?:-[a-zA-Z]+\s*)|(?:-o\s+\w+\s*))*;\s*")

# A shell word making a FILE redirection. ``2>&1`` is the stream merge and
# is allowed; ``>``/``>>``/``2>``/``2>>`` writing a file are refused.
_FILEREDIRECT = re.compile(r"(?<!\d)(?:[12]?|&)>{1,2}(?!&1\b)")

# Heredocs / here-strings, command substitution and process substitution
# (``<( python … )`` runs a program into a pipe).
_HEREDOC = re.compile(r"<<|<<<")
_CMDSUBST = re.compile(r"\$\(|`")
_PROCSUBST = re.compile(r"<\(|>\(|\$<")

_GIT_WORD = re.compile(r"(?:\S*/)?git\b")


def _strip_quotes(s: str) -> str:
    """Remove enclosing quotes/backslashes so ``git 'log'`` reads as git."""
    out = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "\\" and i + 1 < n:
            out.append(s[i + 1])
            i += 2
            continue
        if c in ("'", '"'):
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def _split_segments(cmd: str) -> list[str]:
    """Split on unquoted ``|``, ``;``, ``&&``, ``||`` and newlines.

    A slower, standalone sibling of the gate's splitter -- this module
    must not depend on api_client internals. Operators inside single or
    double quotes (a grep pattern ``'a|b'``, ``git log -p ';'``) stay
    intact. Empty segments are dropped.
    """
    segs: list[str] = []
    buf: list[str] = []
    q: str | None = None
    i, n = 0, len(cmd)
    while i < n:
        c = cmd[i]
        if c == "\\" and q != "'" and i + 1 < n:
            buf.append(cmd[i:i + 2])
            i += 2
            continue
        if q is not None:
            buf.append(c)
            if c == q:
                q = None
            i += 1
            continue
        if c in ("'", '"'):
            q = c
            buf.append(c)
            i += 1
            continue
        if cmd[i:i + 2] in ("||", "&&", ";;"):
            if "".join(buf).strip():
                segs.append("".join(buf))
            buf = []
            i += 2
            continue
        if c in ("|", ";", "\n"):
            if "".join(buf).strip():
                segs.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    if "".join(buf).strip():
        segs.append("".join(buf))
    return segs


def _command_of(segment: str) -> str:
    """The command word of a segment (leading env/whitespace discarded)."""
    words = _strip_quotes(segment).split()
    for w in words:
        if not w:
            continue
        if "=" in w and not w.startswith("=") and _is_env_assignment(w):
            continue
        return w
    return ""


def _is_env_assignment(word: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=[^;|&\s]*", word))


def _git_subcommand(segment: str) -> str:
    """The first non-option subcommand after ``git`` ('' when none).

    Options that carry a VALUE (``-C dir``, ``--git-dir=...``) are
    stepped over whole so their value is not mistaken for the
    subcommand.
    """
    stripped = _strip_quotes(segment)
    m = _GIT_WORD.search(stripped)
    if not m:
        return ""
    rest = stripped[m.end():].split()
    i = 0
    while i < len(rest):
        w = rest[i]
        if w in ("-C", "--git-dir", "--work-tree", "--exec-path"):
            i += 2
            continue
        if w.startswith("-"):
            i += 1
            continue
        return w
    return ""


def _git_line_safe(segment: str) -> bool:
    """Read-only git? A read-only subcommand, and for the ``branch``/``tag``
    subcommands (which ALSO carry every write) only the LIST form -- the
    line after the subcommand holds plain flags and none of them writes.
    ``git branch -d x``, ``git tag v1.0``, ``git branch -m a b`` are all
    refused because a True here skips the confirm dialog.
    """
    text = _strip_quotes(segment)
    m = _GIT_WORD.search(text)
    if not m:
        return False
    rest = text[m.end():].split()
    i = 0
    while i < len(rest):
        w = rest[i]
        if w in ("-C", "--git-dir", "--work-tree", "--exec-path"):
            i += 2
            continue
        if w.startswith("-"):
            i += 1
            continue
        sub = w
        break
    else:
        return False
    if sub not in _READONLY_GIT:
        return False
    args = rest[i + 1:]
    if sub not in _GIT_LIST_ONLY:
        return True
    for a in args:
        if not a.startswith("-"):
            return False
        if a in _GIT_WRITE_FLAGS:
            return False
    return True


def _segment_is_safe(segment: str) -> bool:
    """True when one pipeline/sequence segment is a read-only command.

    Only the COMMAND position matters for "is it a runner": a runner as
    an argument (``grep python file``) is inert text, not a launch. The
    command is whitelisted to intrinsic read-only commands or read-only
    ``git`` subcommands; everything else -- sed, awk, find, xargs, any
    interpreter -- is refused.
    """
    text = _strip_quotes(segment)
    if _FILEREDIRECT.search(text):
        return False
    if _HEREDOC.search(text):
        return False
    if _CMDSUBST.search(segment) or _PROCSUBST.search(segment):
        return False
    cmd = _command_of(segment)
    if not cmd:
        return False
    if cmd == "git":
        return _git_line_safe(segment)
    if cmd in _RUNNERS:
        return False
    if cmd not in _INTRINSIC_READONLY:
        return False
    banned = _WRITE_FLAG.get(cmd)
    if banned:
        for tok in _strip_quotes(segment).split():
            if tok in banned or _flag_base(tok) in banned:
                return False
    return True


def _flag_base(tok: str) -> str:
    """``--output=file`` -> ``--output`` so a value form matches its flag."""
    return tok.split("=", 1)[0]


def is_safe(cmd: str) -> bool:
    """True only when *cmd* is a pipeline/sequence of read-only commands.

    Conservative: any segment that is not unambiguously read-only makes
    the whole command unsafe (False). An empty or whitespace-only string
    is unsafe -- there is nothing safe to approve about it.
    """
    if not isinstance(cmd, str) or not cmd.strip():
        return False
    cmd2 = _SHELL_OPT_PREFIX.sub("", cmd)
    segments = _split_segments(cmd2)
    if not segments:
        return False
    return all(_segment_is_safe(s) for s in segments)
