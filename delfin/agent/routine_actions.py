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

Look-alikes to stay judged (all red-tested in
tests/test_fr1_routine_shapes.py): shell operators anywhere in the
command, force flags (-f/-B/-C/-k), paths after ``checkout``, wildcard
mv, non-http(s) URLs, a POST body or extra arguments on web_fetch.
"""

from __future__ import annotations

from collections import namedtuple

from . import action_classifier as _ac

__all__ = ["ROUTINE_KINDS", "routine_verdict"]

Verdict = _ac.Verdict

ROUTINE_KINDS = ("merge_abort", "new_local_branch", "git_mv", "web_fetch_get")


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


def _git_routine(tokens: list, sub: str) -> Verdict | None:
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
        tokens = _plain_tokens(str(action.get("command") or ""))
        if tokens and tokens[0] == "git":
            for sub in ("merge", "checkout", "switch", "mv"):
                v = _git_routine(tokens, sub)
                if v is not None:
                    return v
        return None
    if tool == "web_fetch":
        return _url_routine(action)
    return None
