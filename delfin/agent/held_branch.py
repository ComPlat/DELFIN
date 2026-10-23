"""Run a branch that another worktree holds as a detached checkout.

git lets one branch name be checked out in one worktree at a time; a
second ``git switch main`` fails with "fatal: 'main' is already used by
worktree at '<path>'" although the code is the same commit. This module
checks out that commit with a detached HEAD instead and says where the
name is held.

Input: a checkout path and a branch name. Output: a :class:`MoveResult`
(``ok``, ``message``, ``sha``, ``held_at``, ``target``). Semantics:

* tracked uncommitted changes refuse the move (nothing is stashed, reset
  or carried along); untracked files are left as git leaves them;
* the target is the local branch, or its remote-tracking ref when
  ``remote=True`` and that ref is strictly ahead of the local branch;
* :func:`update_detached` moves a checkout that was detached at a branch
  forward to that branch's current commit, fast-forward only;
* the branch a checkout was detached at is read from the last HEAD
  reflog entry (``checkout: moving from X to main``), so no marker file is
  written and a detach made in a terminal is recognised the same way.

Rejected: a per-worktree marker file (a second source of truth that a
terminal ``git checkout`` would leave stale); ``--ignore-other-worktrees``
(two worktrees on one branch name, where a commit in one silently moves
the other's branch under its index).
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

#: git >= 2.42 says "already used by worktree at", older "already checked out at".
HELD_RE = re.compile(
    r"'(?P<branch>[^']+)' is already (?:used by worktree|checked out) at "
    r"'(?P<path>[^']+)'")

NOT_CLEAN = ("working tree is not clean. Please commit, stash, or discard "
             "local changes first.")


@dataclass(frozen=True)
class MoveResult:
    ok: bool
    message: str
    sha: str = ""
    held_at: str = ""
    target: str = ""


def _git(repo, *args: str, timeout: float = 30.0) -> subprocess.CompletedProcess:
    from . import safe_git
    return safe_git.run(args, repo, timeout=timeout)


def _out(repo, *args: str) -> str:
    r = _git(repo, *args)
    return r.stdout.strip() if r.returncode == 0 else ""


def parse_held_error(text: str) -> Optional[tuple[str, str]]:
    """``(branch, holder path)`` from git's refusal, or None."""
    m = HELD_RE.search(text or "")
    return (m.group("branch"), m.group("path")) if m else None


def _same_path(a: str, b: str) -> bool:
    try:
        return Path(a).resolve() == Path(b).resolve()
    except OSError:
        return a == b


def held_elsewhere(repo, branch: str) -> Optional[str]:
    """Path of another worktree that has *branch* checked out, or None."""
    top = _out(repo, "rev-parse", "--show-toplevel")
    raw = _out(repo, "worktree", "list", "--porcelain")
    path = ""
    for line in raw.splitlines():
        if line.startswith("worktree "):
            path = line[len("worktree "):]
        elif line == f"branch refs/heads/{branch}" and path:
            if not (top and _same_path(path, top)):
                return path
    return None


def _short(repo, ref: str = "HEAD") -> str:
    return _out(repo, "rev-parse", "--short", ref)


def _is_ancestor(repo, a: str, b: str) -> bool:
    return _git(repo, "merge-base", "--is-ancestor", a, b).returncode == 0


def _tracked_changes(repo) -> str:
    return _out(repo, "status", "--porcelain", "--untracked-files=no")


def _upstream(repo, branch: str) -> str:
    return _out(repo, "rev-parse", "--abbrev-ref", "--symbolic-full-name",
                f"{branch}@{{upstream}}")


def _fetch(repo, upstream: str) -> None:
    remote, _, rbranch = upstream.partition("/")
    if remote and rbranch:
        _git(repo, "fetch", "--quiet", remote, rbranch, timeout=60.0)


def detached_from(repo) -> Optional[str]:
    """Branch a detached checkout was moved to, or None.

    None when HEAD is on a branch, or when the last HEAD movement was not
    a checkout of a branch (a commit, a reset, a tag)."""
    if _out(repo, "branch", "--show-current"):
        return None
    if _git(repo, "symbolic-ref", "-q", "HEAD").returncode == 0:
        return None
    # Newest first: commits made on the detached HEAD since are skipped;
    # any other movement (reset, rebase, ...) means it is not "at" a branch.
    m = None
    for subject in _out(repo, "reflog", "-20", "--format=%gs",
                        "HEAD").splitlines():
        m = re.match(r"checkout: moving from \S+ to (\S+)$", subject)
        if m:
            break
        if not subject.startswith("commit"):
            return None
    if not m:
        return None
    name = m.group(1)
    if _git(repo, "show-ref", "--verify", "--quiet",
            f"refs/heads/{name}").returncode == 0:
        return name
    if _git(repo, "show-ref", "--verify", "--quiet",
            f"refs/remotes/{name}").returncode == 0 and "/" in name:
        local = name.split("/", 1)[1]
        if _upstream(repo, local) == name:
            return local
    return None


def checkout_label(repo) -> str:
    """The branch name, or ``main @ 6bc0297 (detached, ...)`` when detached.

    Empty when *repo* is not a git checkout."""
    cur = _out(repo, "branch", "--show-current")
    if cur:
        return cur
    short = _short(repo)
    if not short:
        return ""
    name = detached_from(repo)
    if not name:
        return f"HEAD @ {short} (detached)"
    notes = []
    if held_elsewhere(repo, name):
        notes.append("held by another checkout")
    tip = _out(repo, "rev-parse", name)
    up = _upstream(repo, name)
    newest = up if up and _is_ancestor(repo, tip, up) else name
    if not _is_ancestor(repo, "HEAD", newest):
        notes.append(f"with commits {name} does not have")
    elif _out(repo, "rev-parse", newest) != _out(repo, "rev-parse", "HEAD"):
        notes.append(f"{name} has moved on")
    return f"{name} @ {short} (detached" + "".join(
        ", " + n for n in notes) + ")"


def detach_at_branch(repo, branch: str, *, remote: bool = False,
                     fast_forward_only: bool = False) -> MoveResult:
    """Check out *branch*'s commit with a detached HEAD.

    ``remote=True`` fetches the branch's upstream and takes it when it is
    ahead of the local branch. ``fast_forward_only`` refuses a target that
    HEAD is not an ancestor of."""
    refusal = can_detach(repo, branch)
    if refusal:
        return MoveResult(False, refusal)
    target = branch
    up = _upstream(repo, branch)
    if up:
        _fetch(repo, up)
        if remote and _is_ancestor(repo, branch, up) and (
                _out(repo, "rev-parse", branch) != _out(repo, "rev-parse", up)):
            target = up
    if fast_forward_only and not _is_ancestor(repo, "HEAD", target):
        return MoveResult(
            False,
            f"Update aborted: this checkout has commits that {target} does "
            f"not contain, so moving to {target} is not a fast-forward. "
            "Name them first with git branch <name> (no switch needed).",
            target=target)
    before = _out(repo, "rev-parse", "HEAD")
    r = _git(repo, "checkout", "--quiet", "--detach", target)
    if r.returncode != 0:
        return MoveResult(False, (r.stderr or r.stdout).strip()
                          or f"git checkout --detach {target} failed",
                          target=target)
    return detached_result(
        repo, branch, target,
        unchanged=fast_forward_only
        and before == _out(repo, "rev-parse", "HEAD"))


def detached_result(repo, branch: str, target: str, *,
                    unchanged: bool = False) -> MoveResult:
    """The success report after HEAD was detached at *target*."""
    sha = _short(repo)
    holder = held_elsewhere(repo, branch) or ""
    where = f"{branch} is checked out in {holder}; " if holder else ""
    state = "already runs" if unchanged else "now runs"
    src = f" ({target})" if target != branch else ""
    msg = f"{where}this checkout {state} the same commit {sha}{src} detached."
    if holder:
        msg += (f" To take the name {branch} here, release it there first "
                "(git checkout --detach in that checkout).")
    return MoveResult(True, msg, sha=sha, held_at=holder, target=target)


def can_detach(repo, branch: str) -> Optional[str]:
    """None when *branch* can be detached at here, else the refusal."""
    if _git(repo, "show-ref", "--verify", "--quiet",
            f"refs/heads/{branch}").returncode != 0:
        return f"Branch not available locally: {branch}"
    if _tracked_changes(repo):
        return f"Switch to {branch} aborted: {NOT_CLEAN}"
    return None


def update_detached(repo) -> Optional[MoveResult]:
    """Fast-forward a checkout detached at a branch to its current commit.

    None when HEAD is on a branch or was not detached at one."""
    name = detached_from(repo)
    if not name:
        return None
    return detach_at_branch(repo, name, remote=True, fast_forward_only=True)


def plain_branch_move(command: str) -> Optional[str]:
    """The branch of ``git checkout <b>`` / ``git switch <b>``, else None.

    Only that exact two-word form: any option, path or second command
    makes it something else, which is left to git's own answer."""
    try:
        words = shlex.split(command or "")
    except ValueError:
        return None
    if (len(words) == 3 and words[0] == "git"
            and words[1] in ("checkout", "switch")
            and not words[2].startswith("-")):
        return words[2]
    return None
