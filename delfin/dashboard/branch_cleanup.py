"""Clean up: merged branches and stale worktrees, deleted on request, undoable.

Input: a checkout, an optional ``gh``, ``stale_days`` and the paths of
worktrees that must not be touched (running previews). Output of
:func:`candidates`: ``(candidates, kept)``.

A local **branch** is a candidate when it is fully merged into
``origin/main``, or when its pull request is merged or closed (a squash
merge leaves its commits outside ``origin/main``; ``unmerged`` counts
them). Never a candidate: ``main``, the branch the checkout is on, a
branch checked out in a worktree that is kept. A branch held only by a
worktree that is itself a candidate carries ``after=<that path>`` and is
deleted after it, in the same run.

A **worktree** (never the main one, never the checkout asked about, never
a protected path) is a candidate when its directory is gone, or when it
is clean and its branch is merged / has a merged or closed pull request /
lost its upstream, or it is detached at a commit ``origin/main`` contains,
or nothing changed in it for ``stale_days``. Kept, with the reason:
uncommitted or untracked changes, commits on no branch, a lock, a live
agent session (its sidecar, asked through
``agent_sessions._worktree_is_orphaned``), running background jobs.

:func:`delete` re-checks every rule at the moment of deletion -- the sha
must be the one listed, the branch must not be ``main``, current or
checked out anywhere -- and needs ``confirm_unmerged=True`` for a branch
with commits not in ``origin/main``. It runs ``git worktree remove``
(without ``--force``: git refuses a tree that is not clean) and
``git branch -D``; nothing is pushed and no remote-tracking ref is
touched. Every deletion is appended to ``~/.delfin/branch_cleanup.jsonl``
with the sha and the command that restores it.

Rejected: ``git worktree remove --force`` (deletes uncommitted work that
no log can bring back) and deleting remote branches (other people's
checkouts track them; the remote is not this list's business).
"""

from __future__ import annotations

import json
import os
import shlex
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from delfin.dashboard import branch_list as BL

DEFAULT_STALE_DAYS = 14


@dataclass
class Candidate:
    kind: str                 # branch | worktree
    name: str                 # branch name, or the worktree path
    sha: str
    reason: str
    unmerged: int = 0         # commits not in origin/main
    path: str = ""            # worktree path
    branch: str = ""          # the worktree's branch, "" when detached
    after: str = ""           # branch: the worktree removed first
    upstream: str = ""
    prunable: bool = False

    @property
    def key(self) -> str:
        return f"{self.kind}:{self.name}"


@dataclass(frozen=True)
class Kept:
    kind: str
    name: str
    why: str


@dataclass
class Report:
    deleted: list = field(default_factory=list)
    skipped: list = field(default_factory=list)   # (name, why)
    undo: list = field(default_factory=list)
    log_path: str = ""


def _log_path() -> Path:
    """Resolved per call, so a redirected home moves it."""
    return Path.home() / ".delfin" / "branch_cleanup.jsonl"


def _same(a, b) -> bool:
    try:
        return Path(a).resolve() == Path(b).resolve()
    except OSError:
        return str(a) == str(b)


def worktrees(repo) -> list[dict]:
    """``git worktree list --porcelain`` as dicts; the first is the main one."""
    raw = BL._out(repo, "worktree", "list", "--porcelain")
    out: list[dict] = []
    cur: dict = {}
    for line in raw.splitlines() + [""]:
        if not line.strip():
            if cur.get("worktree"):
                out.append(cur)
            cur = {}
            continue
        key, _, value = line.partition(" ")
        cur[key] = value if value else True
    for row in out:
        row["branch"] = str(row.get("branch") or "").replace("refs/heads/", "")
    return out


def _held_by(repo) -> dict[str, list[str]]:
    held: dict[str, list[str]] = {}
    for row in worktrees(repo):
        if row.get("branch"):
            held.setdefault(row["branch"], []).append(row["worktree"])
    return held


def _toplevel(repo) -> str:
    return BL._out(repo, "rev-parse", "--show-toplevel") or str(repo)


def _upstream(repo, branch: str) -> tuple[str, bool]:
    """(upstream name, gone?)"""
    raw = BL._out(repo, "for-each-ref",
                  "--format=%(upstream:short)|%(upstream:track)",
                  f"refs/heads/{branch}")
    up, _, track = raw.partition("|")
    return up, "gone" in track


def _changes(path: str) -> str:
    return BL._out(path, "status", "--porcelain")


def _commits_on_no_branch(path: str) -> int:
    raw = BL._out(path, "rev-list", "--count", "HEAD", "--not",
                  "--branches", "--remotes", "--tags")
    try:
        return int(raw)
    except ValueError:
        return 0


def _last_activity(path: str) -> float:
    times = []
    raw = BL._out(path, "log", "-1", "--format=%ct", "HEAD")
    try:
        times.append(float(raw))
    except ValueError:
        pass
    gitdir = BL._out(path, "rev-parse", "--path-format=absolute", "--git-dir")
    if gitdir:
        for name in ("index", "HEAD"):
            try:
                times.append(os.stat(os.path.join(gitdir, name)).st_mtime)
            except OSError:
                pass
    return max(times) if times else 0.0


def _session_hold(path: str) -> str:
    """Why an agent session worktree must stay, or ""."""
    try:
        from delfin.dashboard import agent_sessions as AS
        side = AS.read_worktree_sidecar(path)
        if not side:
            return ""
        if not _same(side.get("path") or "", path):
            return ""
        why = AS._worktree_is_orphaned(
            side, saved_workspaces=AS._saved_session_workspaces())
        return f"agent session worktree: {why}" if why else ""
    except Exception:
        return ""


def _jobs_hold(path: str) -> str:
    try:
        from delfin.agent.worktree import jobs_holding_worktree
        jobs = jobs_holding_worktree(path)
    except Exception:
        jobs = []
    return f"{len(jobs)} background job(s) running in it" if jobs else ""


def _pr_note(prs: dict, branch: str):
    pr = prs.get(branch)
    if pr is not None and pr.state in ("MERGED", "CLOSED"):
        return pr
    return None


def candidates(repo, *, gh: Optional[str] = "auto",
               stale_days: float = DEFAULT_STALE_DAYS,
               now: Optional[float] = None,
               protect_paths: Iterable[str] = (),
               refresh: bool = False) -> tuple[list[Candidate], list[Kept]]:
    """What Clean up offers to delete, and what it keeps (with why)."""
    now = time.time() if now is None else float(now)
    base = BL.base_ref(repo)
    prs = BL.pull_requests(repo, gh, refresh=refresh) or {}
    merged = set(BL.merged_local_branches(repo, base))
    current = BL.current_branch(repo)
    top = _toplevel(repo)
    protect = [str(p) for p in protect_paths]
    cands: list[Candidate] = []
    kept: list[Kept] = []
    wt_cands: dict[str, Candidate] = {}

    rows = worktrees(repo)
    for i, row in enumerate(rows):
        path = row["worktree"]
        if i == 0 or _same(path, top):
            continue
        branch = row.get("branch") or ""
        sha = str(row.get("HEAD") or "")
        if any(_same(path, p) for p in protect):
            kept.append(Kept("worktree", path, "a running branch preview"))
            continue
        if row.get("locked"):
            kept.append(Kept("worktree", path, "locked (git worktree lock)"))
            continue
        if row.get("prunable") or not os.path.isdir(path):
            c = Candidate("worktree", path, sha, "its directory is gone",
                          path=path, branch=branch, prunable=True)
            cands.append(c)
            wt_cands[path] = c
            continue
        hold = (_session_hold(path) or _jobs_hold(path))
        if hold:
            kept.append(Kept("worktree", path, hold))
            continue
        if _changes(path):
            kept.append(Kept("worktree", path, "uncommitted or untracked changes"))
            continue
        if not branch and _commits_on_no_branch(path):
            kept.append(Kept("worktree", path,
                             "has commits that are on no branch"))
            continue
        reason = ""
        if branch:
            pr = _pr_note(prs, branch)
            _up, gone = _upstream(repo, branch)
            if branch in merged:
                reason = f"branch {branch} is merged into {base}"
            elif pr is not None:
                reason = f"branch {branch}: PR #{pr.number} {pr.state.lower()}"
            elif gone:
                reason = f"branch {branch}: its upstream is gone"
        elif sha and BL._git(repo, "merge-base", "--is-ancestor", sha,
                             base).returncode == 0:
            reason = f"detached at a commit {base} contains"
        if not reason:
            idle = now - _last_activity(path)
            if idle > float(stale_days) * 86400:
                reason = f"no change for {int(idle // 86400)} days"
        if not reason:
            kept.append(Kept("worktree", path, "in use (recent, not merged)"))
            continue
        c = Candidate("worktree", path, sha, reason, path=path, branch=branch)
        cands.append(c)
        wt_cands[path] = c

    held = _held_by(repo)
    refs: dict[str, tuple[str, str]] = {}
    for line in BL._out(repo, "for-each-ref",
                        "--format=%(refname:short)|%(objectname)|%(upstream:short)",
                        "refs/heads").splitlines():
        name, _, rest = line.partition("|")
        sha, _, up = rest.partition("|")
        refs[name.strip()] = (sha.strip(), up.strip())
    for name in sorted(refs):
        if not name or name == "main" or name == current:
            continue
        pr = _pr_note(prs, name)
        if name in merged:
            reason = f"merged into {base}"
        elif pr is not None:
            reason = f"PR #{pr.number} {pr.state.lower()}"
        else:
            continue
        holders = held.get(name, [])
        after = ""
        if holders:
            removable = [p for p in holders if p in wt_cands]
            if len(removable) != len(holders) or len(holders) != 1:
                continue
            after = removable[0]
        sha, up = refs[name]
        unmerged = 0 if name in merged else (BL.ahead_of_base(repo, name, base) or 0)
        cands.append(Candidate("branch", name, sha, reason, unmerged=unmerged,
                               after=after, upstream=up))
    return cands, kept


def _branch_refusal(repo, name: str) -> str:
    if name == "main":
        return "main is never deleted"
    if name == BL.current_branch(repo):
        return "it is the branch this checkout is on"
    for row in worktrees(repo):
        if row.get("branch") == name:
            return f"it is checked out in {row['worktree']}"
    if BL._git(repo, "show-ref", "--verify", "--quiet",
               f"refs/heads/{name}").returncode != 0:
        return "it is not a local branch"
    return ""


def _undo_branch(repo, name: str, sha: str, upstream: str) -> str:
    q = shlex.quote
    cmd = f"git -C {q(str(repo))} branch {q(name)} {sha}"
    if upstream:
        cmd += f" && git -C {q(str(repo))} branch --set-upstream-to={q(upstream)} {q(name)}"
    return cmd


def _undo_worktree(repo, c: Candidate) -> str:
    q = shlex.quote
    if c.branch:
        return f"git -C {q(str(repo))} worktree add {q(c.path)} {q(c.branch)}"
    return f"git -C {q(str(repo))} worktree add --detach {q(c.path)} {c.sha}"


def _log(entry: dict) -> str:
    try:
        from delfin.agent.state_paths import ensure_dir, open_append, secure_file
        path = _log_path()
        ensure_dir(path.parent)
        with open_append(path) as fh:
            fh.write(json.dumps(entry) + "\n")
        secure_file(path)
        return str(path)
    except Exception:
        return ""


def _log_writable() -> bool:
    try:
        from delfin.agent.state_paths import ensure_dir, open_append
        path = _log_path()
        ensure_dir(path.parent)
        with open_append(path):
            pass
        return True
    except Exception:
        return False


def delete(repo, selected: Iterable[Candidate], *,
           confirm_unmerged: bool = False) -> Report:
    """Delete *selected* (see the module docstring). Worktrees go first."""
    report = Report()
    selected = list(selected)
    if selected and not _log_writable():
        report.skipped = [(c.name, f"the undo log {_log_path()} cannot be "
                           "written, so nothing is deleted") for c in selected]
        return report
    top = _toplevel(repo)
    removed_paths: set[str] = set()
    order = ([c for c in selected if c.kind == "worktree"]
             + [c for c in selected if c.kind == "branch"])
    for c in order:
        if c.kind == "worktree":
            rows = {r["worktree"]: r for r in worktrees(repo)}
            row = next((r for p, r in rows.items() if _same(p, c.path)), None)
            if row is None:
                report.skipped.append((c.name, "no longer a worktree"))
                continue
            if _same(c.path, top) or _same(c.path, next(iter(rows), "")):
                report.skipped.append((c.name, "the main checkout is never removed"))
                continue
            if c.prunable or not os.path.isdir(c.path):
                r = BL._git(repo, "worktree", "prune", "--expire=now")
            else:
                if str(row.get("HEAD") or "") != c.sha:
                    report.skipped.append((c.name, "it moved since the listing"))
                    continue
                if _changes(c.path):
                    report.skipped.append((c.name, "it has uncommitted changes now"))
                    continue
                r = BL._git(repo, "worktree", "remove", c.path)
            if r.returncode != 0:
                report.skipped.append((c.name, (r.stderr or r.stdout).strip()
                                       or "git worktree remove failed"))
                continue
            removed_paths.add(c.path)
            undo = _undo_worktree(repo, c)
            entry = {"at": time.time(), "repo": str(repo), "kind": "worktree",
                     "name": c.path, "path": c.path, "branch": c.branch,
                     "sha": c.sha, "reason": c.reason, "undo": undo}
            report.log_path = _log(entry) or report.log_path
            report.deleted.append(entry)
            report.undo.append(undo)
            continue
        # a branch
        if c.after and c.after not in removed_paths:
            report.skipped.append((c.name, f"its worktree {c.after} was not removed"))
            continue
        why = _branch_refusal(repo, c.name)
        if why:
            report.skipped.append((c.name, why))
            continue
        sha = BL._out(repo, "rev-parse", f"refs/heads/{c.name}")
        if sha != c.sha:
            report.skipped.append((c.name, "it moved since the listing"))
            continue
        base = BL.base_ref(repo)
        unmerged = BL.ahead_of_base(repo, c.name, base) or 0
        if unmerged and not confirm_unmerged:
            report.skipped.append((c.name, (
                f"has {unmerged} commit(s) not in {base}; needs the second "
                "confirmation")))
            continue
        up, _gone = _upstream(repo, c.name)
        r = BL._git(repo, "branch", "-D", c.name)
        if r.returncode != 0:
            report.skipped.append((c.name, (r.stderr or r.stdout).strip()))
            continue
        undo = _undo_branch(repo, c.name, sha, up)
        entry = {"at": time.time(), "repo": str(repo), "kind": "branch",
                 "name": c.name, "sha": sha, "upstream": up,
                 "unmerged": unmerged, "reason": c.reason, "undo": undo}
        report.log_path = _log(entry) or report.log_path
        report.deleted.append(entry)
        report.undo.append(undo)
        try:
            BL.forget(repo, c.name)
        except Exception:
            pass
    return report


def delete_branch(repo, name: str, *, confirm_unmerged: bool = False
                  ) -> tuple[bool, str, Report]:
    """Delete one branch picked in the list, under the rules of :func:`delete`."""
    why = _branch_refusal(repo, name)
    if why:
        return False, f"{name} was not deleted: {why}.", Report()
    sha = BL._out(repo, "rev-parse", f"refs/heads/{name}")
    c = Candidate("branch", name, sha, "picked in the branch list",
                  unmerged=BL.ahead_of_base(repo, name) or 0)
    report = delete(repo, [c], confirm_unmerged=confirm_unmerged)
    if report.deleted:
        return True, f"Deleted {name} ({sha[:10]}). Undo: {report.undo[0]}", report
    why = report.skipped[0][1] if report.skipped else "not deleted"
    return False, f"{name} was not deleted: {why}.", report


def unmerged_count(repo, name: str) -> int:
    return BL.ahead_of_base(repo, name) or 0
