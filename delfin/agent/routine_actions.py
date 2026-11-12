"""Routine actions that run without a dialog.

Goal of the person (wave 20): agents as free as possible and as safe as
possible -- legitimate normal actions are allowed, as for a careful
human developer. Measured by the wave, the operator answered by hand,
again and again, for four shapes that are entirely ordinary in one's
own tree: aborting a merge, starting a new local branch, moving one of
one's own files inside the write scope, and reading one journal page.
Each shape here is a small, explicit condition checked BY CODE, not by
the model: every argument is matched against the shape, and anything
that does not fit exactly is not routine and is judged as today.

Nothing here turns a refusal into an allow. ``routine_verdict`` is
called from ``action_classifier.screen`` AFTER the user's "never" rules
and the injection guard, so a user who said "never" keeps precedence and
hostile content still asks. The path gates (write scope, sandbox, deny
list) run independently of this verdict -- the classifier only answers
"the model need not be asked about this shape"; an out-of-scope ``git
mv`` still stops at the write gate.

Deliberately NOT routine (operator, stage 1d): under bypass mode curl,
wget, pip install, ssh to other nodes and ``git commit -am`` already
run in wave trees, so none of those shapes may become routine. A bash
command naming them is never matched here; every one of them is caught
by the exact-shape matching (a git shape only matches ``git`` with one
of the three sub-commands, and no web/curl/wget/pip/ssh shape exists).

A list of plain reads (``grep x f | head``, ``git status; git log -3``)
is the fifth shape (2026-11-11): every part must pass the plain-read
floor on its own; only 2>&1 and 2>/dev/null may redirect. Red tests in
tests/test_a_list_of_reads_is_routine.py.

Look-alikes to stay judged (all red-tested in
tests/test_fr1_routine_shapes.py): shell operators anywhere in the
command, force flags (-f/-B/-C/-k), paths after ``checkout``, wildcard
mv, non-http(s) URLs, a POST body or extra arguments on web_fetch.
"""

from __future__ import annotations

import re
import shlex
from collections import namedtuple

from . import action_classifier as _ac

__all__ = ["ROUTINE_KINDS", "routine_verdict"]

Verdict = _ac.Verdict

ROUTINE_KINDS = ("merge_abort", "new_local_branch", "existing_local_branch",
                 "git_mv", "web_fetch_get", "compound_read")


# Shell operators or nested execution anywhere in a command: the moment
# one appears the command is more than the one action being judged.
# A newline separates commands exactly like ``;`` does.
_ROUTINE_FORBIDDEN_CHARS = (">", ">>", ">|", ";", "&&", "||", "|", "<",
                            "<<", "\n", "\r", "$", "`", "&")


def _plain_tokens(command: str):
    """The command's tokens, or None when anything forbidden is in it."""
    cmd = (command or "").strip()
    if not cmd:
        return None
    for ch in _ROUTINE_FORBIDDEN_CHARS:
        if ch in cmd:
            return None
    return cmd.split()


# compound_read: the list operators that join two plain reads, and the
# two redirects that only move or drop the error stream. Every other
# operator (a file redirect, ``<``, ``&``, ``|&``, ``&>``) ends the shape.
_LIST_OPERATORS = ("&&", "||", "|", ";")
_STDERR_REDIRECTS = ("2>&1", "2>/dev/null")
_OPERATOR_CHARS = frozenset(";&|<>")
# find runs a program per match with these; inside a list the ``;`` that
# ends them is quoted and so no longer visible to the plain-read floor.
_FIND_RUNNERS = frozenset({"-exec", "-execdir", "-ok", "-okdir"})
# ``sed -n N,Mp`` after a pipe reads stdin; the floor's sed form wants a
# file operand.
_SED_STDIN_RANGE = re.compile(r"\d+(?:,\d+)?p")


def _list_parts(command: str):
    """The parts of a ``| && || ;`` list, or None.

    Input: one bash command. Output: the text of each part, with every
    QUOTED or backslash-escaped operator character replaced by ``_`` so
    the plain-read floor reads it as data (``grep 'a|b'``), not as an
    operator; an option keeps its prefix (``--output=a;b`` becomes
    ``--output=a_b``). The two stderr redirects are removed. None for an
    unbalanced quote, an empty part, or any other unquoted operator.
    Quote state is tracked here, not by shlex: shlex returns a quoted
    ``'|'`` as the same token as an unquoted ``|``, and splitting there
    would hand ``cat --output=x`` to the floor as a part of its own
    (``git log '|' cat --output=x`` writes x).
    """
    def masked(ch):
        return "_" if ch in _OPERATOR_CHARS else ch
    parts, cur, quote, i, n = [], [], None, 0, len(command)
    while i < n:
        ch = command[i]
        if quote:
            if ch == quote:
                quote = None
                cur.append(ch)
            elif ch == "\\" and quote == '"' and i + 1 < n:
                cur.append(ch + masked(command[i + 1]))
                i += 1
            else:
                cur.append(masked(ch))
            i += 1
            continue
        if ch in "'\"":
            quote = ch
            cur.append(ch)
            i += 1
            continue
        if ch == "\\":
            if i + 1 >= n:
                return None
            cur.append(ch + masked(command[i + 1]))
            i += 2
            continue
        if ch in _OPERATOR_CHARS:
            redirect = next((r for r in _STDERR_REDIRECTS
                             if command.startswith(r[1:], i)), None)
            end = i + len(redirect) - 1 if redirect else i
            if (redirect and "".join(cur).endswith((" 2", "\t2"))
                    and (end == n or command[end] in " \t|;&")):
                cur.pop()
                i = end
                continue
            op = next((o for o in _LIST_OPERATORS
                       if command.startswith(o, i)), None)
            if op is None or command.startswith("|&", i):
                return None
            parts.append("".join(cur))
            cur = []
            i += len(op)
            continue
        cur.append(ch)
        i += 1
    if quote:
        return None
    parts.append("".join(cur))
    parts = [p.strip() for p in parts]
    return None if any(not p for p in parts) else parts


def _compound_read_routine(command: str) -> Verdict | None:
    """A ``| && || ;`` list whose every part is a plain read.

    Each part alone passes ``action_classifier._bash_is_plain_read`` and
    never reaches the model; joining reads adds no write and no program
    beyond the parts. Checked per part: the floor, no find runner, and
    two additions -- ``set -o pipefail`` and ``sed -n N,Mp`` reading the
    pipe. ``$`` and backticks end the shape anywhere, quoted or not
    (expansion inside double quotes; a single-quoted ``$`` is rare
    enough to stay judged). Measured: 37 % of the bash commands judged
    by a model on 2026-11-11 had this shape.
    """
    cmd = (command or "").strip()
    if not cmd or any(ch in cmd for ch in "$`\n\r"):
        return None
    parts = _list_parts(cmd)
    if not parts:
        return None
    for n, part in enumerate(parts):
        try:
            words = shlex.split(part)
        except ValueError:
            return None
        if _FIND_RUNNERS & set(words):
            return None
        if words == ["set", "-o", "pipefail"]:
            continue
        if (n > 0 and len(words) == 3 and words[:2] == ["sed", "-n"]
                and _SED_STDIN_RANGE.fullmatch(words[2])):
            continue
        if not _ac._bash_is_plain_read(part):
            return None
    return Verdict(
        "allow", "routine: a list of plain reads; each part is a read "
        "on its own.")


def _git_routine(tokens: list, sub: str,
                 workspace: str = "") -> Verdict | None:
    """Judge the three routine git shapes. tokens[0] is 'git'."""
    if len(tokens) < 2 or tokens[1] != sub:
        return None
    args = tokens[2:]
    if sub == "merge":
        # Exactly ``git merge --abort`` and nothing else.
        if args == ["--abort"]:
            return Verdict(
                "allow", "routine: git merge --abort in the own worktree.")
        return None
    if sub == "checkout" or sub == "switch":
        # A NEW local branch, no force, no paths, no other flags:
        # checkout -b NAME [start] | switch -c NAME [start].
        flag = "-b" if sub == "checkout" else "-c"
        forced = "-B" if sub == "checkout" else "-C"
        if not args or args[0] != flag:
            # GT1: a plain switch to an EXISTING local branch of the
            # session's own worktree is the same tree-hygiene class:
            # exactly one operand, no flags, no "--", nothing detach-
            # shaped, and the operand must resolve to a local branch
            # (verified against git's own ref store, never against a
            # path in the tree). Anything else is judged as today.
            v = _existing_branch_routine(sub, args, workspace)
            if v is not None:
                return v
            return None
        # Positional-only beyond the flag: any other dash-token could
        # carry a destructive mode (--discard-changes, --track, ...) --
        # reject the shape, not a list of spellings.
        if any(a.startswith("-") for a in args[1:]):
            return None
        if "--" in args or len(args) not in (2, 3):
            return None
        return Verdict(
            "allow",
            f"routine: {sub} {flag} starts a new local branch in the "
            "own worktree.")
    if sub == "mv":
        # Exactly ``git mv SRC DST`` (or one --forceness-free flag path is
        # not allowed at all): two positionals, no flags, no wildcards.
        if len(args) != 2:
            return None
        if any(a.startswith("-") for a in args):
            return None
        if any(ch in a for a in args for ch in "*?"):
            return None
        return Verdict(
            "allow", "routine: git mv of an own file inside the scope; "
            "the write gate still holds the path.")
    return None


def _existing_branch_routine(sub: str, args: list,
                             workspace: str = "") -> Verdict | None:
    """GT1: ``git checkout <branch>`` / ``git switch <branch>`` for a
    branch that EXISTS locally.

    The shape is exact: one operand, no flags, no ``--``, nothing that
    carries a detach/force/orphan mode. The operand must resolve to a
    LOCAL branch in the session's own repository (``git show-ref
    --verify refs/heads/<name>``) and must NOT be an existing file --
    if the name also names a file in the tree, the command could be a
    path checkout, and that stays judged (fail-safe direction). Never
    raises: a git probe that fails means "not routine", which is where
    the command already was.
    """
    if len(args) != 1 or args[0].startswith("-") or args[0] == "--":
        return None
    name = args[0]
    # A path-shaped operand could be a pathspec; a bare name with no
    # separator never is, but an EXISTING file wins the ambiguity
    # against the branch reading (git checks the index and the tree
    # before the refs). Probe the tree first, fail-closed.
    # Asked in the session's workspace, where the command runs: the
    # process's cwd is another tree in the dashboard, and a branch there
    # said nothing about a same-named file here.
    if not workspace:
        return None
    try:
        import os
        if os.path.lexists(os.path.join(workspace, name)):
            return None
        # Every git call DELFIN runs goes through safe_git (the guard
        # test pins that). The session's own worktree is the process
        # cwd here -- the routine shape is explicitly "the session's
        # OWN worktree" -- so the repository to ask is Path.cwd().
        from delfin.agent import safe_git
        probe = safe_git.run(
            ["show-ref", "--verify", "--quiet", "refs/heads/" + name],
            workspace, timeout=5)
        if probe.returncode != 0:
            return None
    except Exception:
        return None
    return Verdict(
        "allow",
        f"routine: {sub} to the existing local branch '{name}' of the "
        "own worktree; no paths, no force, no detach.")


def _url_routine(action: dict) -> Verdict | None:
    """Judge a native web_fetch read: http(s) URL, timeout_s at most,
    nothing else in the arguments."""
    if str(action.get("tool") or "") != "web_fetch":
        return None
    url = str(action.get("url") or "").strip()
    if not (url.startswith("http://") or url.startswith("https://")):
        return None
    # No userinfo: user@host in the authority is a credential-carrying
    # URL, not an ordinary literature fetch -- reject the shape.
    if "@" in url.split("//", 1)[-1].split("/", 1)[0]:
        return None
    # Nothing that can carry data out: a hostile page asks the model to
    # fetch ``https://its-host/?d=<file contents>``. A DOI, an article or a
    # manual page needs no query; any other URL is judged as before, not
    # refused. A fragment stays: the client never sends it to the server.
    if "?" in url.split("#", 1)[0] or len(url) > 200:
        return None
    # The action dict always carries "tool"; only url/timeout_s may be
    # in the ARGUMENTS on top of it.
    allowed = {"url", "timeout_s"}
    if set(action) - allowed - {"tool"}:
        return None
    timeout = action.get("timeout_s")
    if timeout is not None:
        if not isinstance(timeout, (int, float)) or isinstance(timeout, bool):
            return None
        if timeout <= 0 or timeout > 60:
            return None
    return Verdict(
        "allow", "routine: read-only literature fetch through DELFIN's "
        "web tools; the untrusted fence applies to the page text.")


def routine_verdict(action: dict, context: dict) -> Verdict | None:
    """A code-checked Verdict("allow") for a routine shape, else None.

    Called from screen() after the user's "never" rules and the
    injection guard. None is not an opinion: the action is judged
    exactly as today (the model prompt, or the plain-read floor).
    Only native bash and native web_fetch are judged here; an MCP
    variant of web_fetch is a different surface and stays with the
    model.
    """
    tool = str(action.get("tool") or "")
    if tool in ("bash", "bash_background"):
        compound = _compound_read_routine(str(action.get("command") or ""))
        if compound is not None:
            return compound
        tokens = _plain_tokens(str(action.get("command") or ""))
        # A command with a cwd of its own runs elsewhere than the
        # workspace the shapes are checked in.
        cwd = str(action.get("cwd") or "").strip()
        workspace = ("" if cwd not in ("", ".")
                     else str(context.get("workspace") or ""))
        if tokens and tokens[0] == "git":
            for sub in ("merge", "checkout", "switch", "mv"):
                v = _git_routine(tokens, sub, workspace)
                if v is not None:
                    return v
        return None
    if tool == "web_fetch":
        return _url_routine(action)
    return None
