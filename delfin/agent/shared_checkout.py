"""Refuse moving the branch of the shared checkout.

Input: a shell command and the directory it runs in. Output of
:func:`refusal`: "" when the command moves no branch of a shared checkout,
else the refusal text with the routes that work.

A *shared checkout* is the checkout this process imports DELFIN from
(``Path(delfin.__file__).parents[1]``) and the one the running dashboard
names in ``DELFIN_SHARED_CHECKOUT``. Everyone works from main there: a
``git switch`` / ``git checkout <branch>`` in it moves every person and
agent session on it at once, and a session that stays open keeps running
the modules it imported before the move. Testing a branch is a preview
(the dashboard header's ⎇ menu and ▶ Test: its own worktree and server);
work on a branch happens in a worktree (``enter_worktree``).

A *branch move* is a git invocation (``git [-C dir] ...``, after ``cd dir``
in the same command line) of ``switch`` (any form but ``--help``), or of
``checkout`` with ``-b``/``-B``/``--orphan``/``--detach``, or with exactly
one operand, no ``--``, that names a commit (``git rev-parse --verify``),
``-`` included. ``git checkout -- <file>`` and ``git checkout <ref> --
<file>`` restore files and are not moves.

Rejected: matching the command text only (``git checkout README.md`` is a
restore, ``git checkout feat`` a move -- only the repository can tell);
refusing every git write in the shared checkout (a commit on main by a
maintainer, a fetch, a status are unaffected by this rule).
"""

from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path
from typing import Iterable, Optional

from . import safe_git

ENV = "DELFIN_SHARED_CHECKOUT"
_SEPARATORS = {";", "&&", "||", "|", "&", "\n"}


def _git(cwd, *args: str) -> subprocess.CompletedProcess:
    return safe_git.run(args, cwd, timeout=15)


def _toplevel(path) -> Optional[Path]:
    r = _git(path, "rev-parse", "--show-toplevel")
    if r.returncode != 0 or not r.stdout.strip():
        return None
    try:
        return Path(r.stdout.strip()).resolve()
    except OSError:
        return None


def shared_checkouts(env=None) -> list[Path]:
    env = os.environ if env is None else env
    out: list[Path] = []
    try:
        import delfin
        code = Path(delfin.__file__).resolve().parent.parent
        if (code / ".git").exists():
            out.append(code)
    except Exception:
        pass
    named = str(env.get(ENV) or "").strip()
    if named:
        top = _toplevel(named)
        if top is not None:
            out.append(top)
    return out


def _segments(command: str) -> list[list[str]]:
    try:
        lex = shlex.shlex((command or "").replace("\n", " ; "), posix=True,
                          punctuation_chars=";&|")
        lex.whitespace_split = True
        toks = list(lex)
    except ValueError:
        return []
    segs: list[list[str]] = [[]]
    for t in toks:
        if t in _SEPARATORS or set(t) <= set(";&|"):
            segs.append([])
        else:
            segs[-1].append(t)
    return [s for s in segs if s]


def _moves(words: list[str], cwd: Path) -> Optional[tuple[Path, str]]:
    """``(directory, sub-command text)`` when *words* is a branch move."""
    while words and "=" in words[0] and words[0].split("=", 1)[0].isidentifier():
        words = words[1:]
    if not words or os.path.basename(words[0]) != "git":
        return None
    i, where = 1, cwd
    while i < len(words) and words[i].startswith("-"):
        w = words[i]
        if w == "-C" and i + 1 < len(words):
            p = Path(os.path.expanduser(words[i + 1]))
            where = p if p.is_absolute() else where / p
            i += 2
            continue
        if w in ("-c", "--git-dir", "--work-tree", "--namespace") and i + 1 < len(words):
            i += 2
            continue
        i += 1
    if i >= len(words):
        return None
    sub, rest = words[i], words[i + 1:]
    if sub == "switch":
        if any(a in ("-h", "--help") for a in rest):
            return None
        return where, " ".join(words[i:])
    if sub != "checkout":
        return None
    if "--" in rest:
        return None
    if any(a in ("-b", "-B", "--orphan", "--detach") or a.startswith("--orphan=")
           for a in rest):
        return where, " ".join(words[i:])
    operands = [a for a in rest if not a.startswith("-") or a == "-"]
    if len(operands) != 1:
        return None
    op = operands[0]
    if op == "-":
        return where, " ".join(words[i:])
    if _git(where, "rev-parse", "--verify", "--quiet",
            f"{op}^{{commit}}").returncode == 0:
        return where, " ".join(words[i:])
    return None


def refusal(command: str, cwd, *, shared: Optional[Iterable] = None) -> str:
    """The refusal for *command* run in *cwd*, or ""."""
    if "git" not in (command or ""):
        return ""
    try:
        shared_tops = [Path(p).resolve() for p in
                       (shared if shared is not None else shared_checkouts())]
    except OSError:
        return ""
    if not shared_tops:
        return ""
    here = Path(cwd or ".")
    for words in _segments(command):
        if words and words[0] == "cd" and len(words) >= 2:
            p = Path(os.path.expanduser(words[1]))
            here = p if p.is_absolute() else here / p
            continue
        hit = _moves(words, here)
        if hit is None:
            continue
        where, text = hit
        top = _toplevel(where)
        if top is None or top not in shared_tops:
            continue
        return (
            f"blocked: `git {text}` would move the branch of the shared "
            f"checkout {top}, which everyone works from and which stays on "
            "main (other sessions and the dashboard run its code). "
            "HINT: for a change and its pull request, "
            "enter_worktree(base_ref=\"origin/main\"), commit there, "
            "`git push -u origin HEAD:<user>/<topic>`, then "
            "`gh pr create --head <user>/<topic>`; to run or test an existing "
            "branch, enter_worktree(base_ref=\"<branch>\"); to only create "
            "one, `git branch <name> origin/main` (no switch). For the user "
            "to try a branch in the dashboard: pick it in the header's ⎇ "
            "menu; ▶ Test opens it in its own worktree and server.")
    return ""
