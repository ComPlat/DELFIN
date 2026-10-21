"""Compose a pull-request body from a session's journal and its commits.

Report 2.7: code handover between sessions is a merge, not a drop — one
session's code was missing for another because it lived on a session branch,
and there was no "session branch -> PR -> main -> others pull main" path.
``gh pr create`` is gated under the push grant (#126); the dashboard PR
button runs the ``/pr`` skill, which until now wrote the title and body by
hand from ``git log --oneline origin/main..HEAD`` and never consulted the
session journal. This module gives that skill a body that actually captures
what the session did and why.

``body(workspace, base)`` turns two data sources — one session's tasks in
the workspace task list (``<workspace>/.delfin/session_tasks.json``) and the
non-merge commits ``base..HEAD`` — into a PR body with what, why, how tested
(the ``Control:`` lines), and open points (blocked tasks). The task loader
and the git reader are injectable, so tests drive it with a fake repository.
The identity (home/account/host) the redaction guards protect can be passed
in; without it the process identity is used.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Optional

#: GitHub refuses a pull-request body over 65536 characters; the margin
#: leaves room for the line that says the text was cut.
BODY_LIMIT = 60000

#: Per-commit body share before the total is checked. A branch with a few
#: long commit messages keeps all of them; a long branch keeps every
#: subject line and the head of each body.
_COMMIT_BODY_LIMIT = 2000

#: ``Co-Authored-By:`` / ``Signed-off-by:`` and similar git trailers name
#: people or tools, not the change; they are not copied into the body.
_TRAILER_RE = re.compile(r"(?im)^[ \t]*[A-Za-z][A-Za-z-]*-by:[^\n]*\n?")

# ---------------------------------------------------------------------------
# Injectable seams (defaults read the real workspace; tests pass fakes)
# ---------------------------------------------------------------------------

#: A task record as the journal keeps it: subject, description, status,
#: blocked_reason, timestamps.
TaskRecord = dict[str, Any]

#: One commit we hand to the body: abbreviated hash, subject, full body.
CommitRecord = dict[str, str]

#: loader(workspace: str) -> list[TaskRecord]
TaskLoader = Callable[[str], list[TaskRecord]]

#: git_log(base: str, head: str) -> list[CommitRecord]
GitLogRunner = Callable[[str, str], list[CommitRecord]]

#: What identity to protect: home (str), account (str), host (str).
Identity = dict[str, str]


def _default_task_loader(workspace: str,
                         session_id: "str | None" = None) -> list[TaskRecord]:
    """The non-deleted tasks of one session in the workspace journal.

    The store is shared by every session that ran in the workspace (one
    checkout held 253 tasks from 147 sessions). With ``session_id`` only
    that session's tasks are returned; without one, only tasks that carry
    no session id (the CLI backend mints none), never another session's.
    """
    from . import agent_tasks

    store = agent_tasks.get_store(Path(workspace))
    sid = str(session_id or "").strip()
    if sid:
        return store.list(include_deleted=False, session_id=sid)
    return [t for t in store.list(include_deleted=False)
            if not str(t.get("session_id", "") or "").strip()]


def _default_git_log(base: str, head: str,
                     cwd: "str | None" = None) -> list[CommitRecord]:
    """``git log --no-merges base..head`` as hash, subject and body.

    Merge commits are left out: on a branch that took main in, their
    subjects ("Merge remote-tracking branch ...") describe no change.
    No prompt and no stdin (a credential or pager question fails instead
    of hanging), a timeout, and safe_git, so a checkout config cannot
    make the read run ``gpg.program``, a textconv or any other program. A failing git
    (unknown base, not a repository) raises ``RuntimeError``: an empty
    list would be published as "no commits".
    """
    import os

    sep = "\x1e"  # record separator, practically absent from commit text
    fmt = f"%h%x1f%s%x1f%b%{sep}"
    # A ref is never an option: ``git log --output=<file>`` writes a file,
    # so a name that starts with a dash is not handed to git at all.
    if str(base).startswith("-") or str(head).startswith("-"):
        return []
    from . import safe_git
    proc = safe_git.run(
        ["log", "--no-merges", f"--format={fmt}", "--end-of-options",
         f"{base}..{head}"], cwd or os.getcwd(), timeout=30)
    if proc.returncode != 0:
        err = (proc.stderr or "").strip().splitlines()
        raise RuntimeError(
            f"git log {base}..{head} failed: {err[0] if err else proc.returncode}")
    commits: list[CommitRecord] = []
    for line in proc.stdout.split(sep):
        parts = line.split("\x1f")
        if len(parts) >= 3:
            commits.append({
                "hash": parts[0].strip(),
                "subject": parts[1].strip(),
                "body": parts[2].strip(),
            })
    return commits


# ---------------------------------------------------------------------------
# Redaction — reuse the shared output guards, never reinvent them
# ---------------------------------------------------------------------------

def _redact(text: str, identity: Identity) -> str:
    """A public-safe version of *text*: secrets out, identity paths out.

    Reuses ``delfin.agent.output_guard.scrub_secrets`` for token-shaped and
    known secret values, and ``home_paths_in`` for ``/home/<acct>`` /
    ``/pfs/<store>/home/<group>/<acct>`` path lookalikes. ``identity`` is
    passed in (never read here), so tests pass a fixture identity and the
    production caller can pass none without the helper probing the user's
    home directory.
    """
    from . import output_guard

    safe = output_guard.scrub_secrets(text)
    hits = output_guard.home_paths_in(
        safe, home=identity.get("home", ""),
        account=identity.get("account", ""), host=identity.get("host", ""),
    )
    for hit in sorted(hits, key=len, reverse=True):
        safe = safe.replace(hit, "[redacted]")
    return safe


def _goal_from(tasks: list[TaskRecord]) -> str:
    """One-sentence goal: the in_progress task, else the newest open one."""
    in_progress = [t for t in tasks if str(t.get("status")) == "in_progress"]
    open_tasks = [t for t in tasks
                  if str(t.get("status")) in ("in_progress", "pending")]
    pick = in_progress or open_tasks
    if not pick:
        return ""
    # Newest first by updated_at (or created_at) when there are several.
    def _ts(t):
        return t.get("updated_at") or t.get("created_at") or ""
    pick = sorted(pick, key=_ts, reverse=True)
    head = pick[0]
    subject = str(head.get("subject", "")).strip()
    description = str(head.get("description", "")).strip()
    if description:
        return f"{subject} — {description}".strip(" .")
    return subject


def _open_points(tasks: list[TaskRecord]) -> list[str]:
    """Blocked tasks and their reason — the open points the reviewer needs."""
    points = []
    blocked = [t for t in tasks if str(t.get("status")) == "blocked"]
    # Newest first, same ordering as _goal_from.
    def _ts(t):
        return t.get("updated_at") or t.get("created_at") or ""
    for t in sorted(blocked, key=_ts, reverse=True):
        subject = str(t.get("subject", "")).strip()
        reason = str(t.get("blocked_reason", "")).strip()
        line = subject if subject else "a blocked task"
        if reason:
            line = f"{line} — {reason}"
        points.append(line)
    return points


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def body(
    workspace: str,
    base: str,
    *,
    head: Optional[str] = None,
    task_loader: Optional[TaskLoader] = None,
    git_log: Optional[GitLogRunner] = None,
    identity: Optional[Identity] = None,
    session_id: Optional[str] = None,
) -> str:
    """Compose a PR body from the journal and ``git log base..HEAD``.

    ``workspace`` is the checkout the session ran in (its task journal lives
    under ``<workspace>/.delfin/session_tasks.json``). ``base`` is the branch
    the PR opens into. ``head`` defaults to the current branch.
    ``session_id`` selects this session's tasks in the shared journal
    (``task_list`` shows it); without it only unowned tasks are read.
    ``identity`` defaults to the home, account and host of the process
    (the same values the commit-message gate protects).

    Returns the full body text, ready for ``gh pr create --body``: a short
    goal line, what changed (one line per commit), why and how it was tested
    (the commit bodies, incl. ``Control:`` lines), and open points (blocked
    tasks). Commit trailers are dropped and the text is bounded by
    ``BODY_LIMIT``. The assembled text passes the shared account-name /
    secret guards.
    """
    load_tasks = task_loader or (
        lambda ws: _default_task_loader(ws, session_id))
    # The commits of THIS workspace, not of wherever the process stands.
    run_log = git_log or (
        lambda b, h: _default_git_log(b, h, cwd=str(workspace)))
    if identity is None:
        from .api_client import _identity_for_commit_texts
        identity = _identity_for_commit_texts()
    who = dict(identity)

    tasks = [dict(t) for t in (load_tasks(workspace) or [])]
    head_ref = head or "HEAD"
    commits = _redact_commits(run_log(base, head_ref), who)

    goal = _goal_from(tasks)
    open_points = _open_points(tasks)

    # Each piece is redacted, and the assembled text once more: the join
    # can put a name next to a path and form a match neither piece had.
    lines = []
    lines.append("# Pull request from a session branch")
    lines.append("")
    if goal:
        lines.append("## Goal")
        lines.append(goal)
        lines.append("")
    lines.append("## What changed")
    if commits:
        for c in commits:
            lines.append(f"- {c['hash']} {c['subject']}")
    else:
        lines.append("No commits between the base and this branch.")
    lines.append("")
    lines.append("## Why / how tested")
    room = BODY_LIMIT - sum(len(x) + 1 for x in lines) - 2000
    for c in commits:
        text = c.get("body") or ""
        if len(text) > _COMMIT_BODY_LIMIT:
            text = text[:_COMMIT_BODY_LIMIT].rstrip() + " [...]"
        if not text:
            continue
        entry = f"{c['hash']}: {text}"
        if len(entry) + 1 > room:
            lines.append("[further commit messages cut: see git log]")
            break
        lines.append(entry)
        room -= len(entry) + 1
    lines.append("")
    lines.append("## Open points")
    if open_points:
        for p in open_points:
            lines.append(f"- {p}")
    else:
        lines.append("None.")
    lines.append("")

    text = _redact("\n".join(lines), who)
    if len(text) > BODY_LIMIT:
        text = text[:BODY_LIMIT - 40].rstrip() + "\n\n[body cut at the length limit]\n"
    return text


def _redact_commits(
    commits: list[CommitRecord], identity: Identity,
) -> list[CommitRecord]:
    """The commits with any identity-bearing text scrubbed from subject/body."""
    out = []
    for c in commits:
        out.append({
            "hash": str(c.get("hash", "")),
            "subject": _redact(str(c.get("subject", "")), identity),
            "body": _redact(
                _TRAILER_RE.sub("", str(c.get("body", ""))).strip(), identity),
        })
    return out
