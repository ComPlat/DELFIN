"""Read-only inventory of git worktrees attached to a DELFIN repo.

Surfacing this in the dashboard makes it obvious when the agent is
running inside an isolated worktree (vs. the main checkout) and lets
users see what other branches are checked out somewhere else on disk.
We intentionally don't add/remove worktrees from the UI — those are
power-user actions and the user should drive them via ``git worktree``.
"""
from __future__ import annotations

import html as _html
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class WorktreeEntry:
    """One row from ``git worktree list --porcelain``."""

    path: str          # absolute worktree path
    head: str          # commit SHA at HEAD
    branch: str        # branch name (or "" for detached / bare)
    is_main: bool      # True if this is the main repo worktree
    is_bare: bool      # True if this is a bare repo


def _run_git(repo: Path, *args: str, timeout: float = 3.0) -> str:
    """Run a git command and return stdout (empty string on failure)."""
    try:
        out = subprocess.run(
            ["git", *args], cwd=str(repo),
            capture_output=True, text=True, timeout=timeout,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        return ""
    if out.returncode != 0:
        return ""
    return out.stdout


def list_worktrees(repo: Path | str) -> list[WorktreeEntry]:
    """Return all worktrees attached to ``repo``.

    Empty list if git is missing, the path is not a repo, or the listing
    fails for any reason.
    """
    raw = _run_git(Path(repo), "worktree", "list", "--porcelain")
    if not raw.strip():
        return []
    out: list[WorktreeEntry] = []
    main_path = ""
    cur: dict[str, str] = {}

    def _flush() -> None:
        if not cur.get("worktree"):
            return
        path = cur["worktree"]
        nonlocal main_path
        if not main_path:
            main_path = path
        out.append(WorktreeEntry(
            path=path,
            head=cur.get("HEAD", ""),
            branch=cur.get("branch", "").replace("refs/heads/", ""),
            is_main=path == main_path,
            is_bare="bare" in cur,
        ))

    for raw_line in raw.splitlines():
        line = raw_line.rstrip()
        if not line:
            _flush()
            cur = {}
            continue
        if " " in line:
            key, _, value = line.partition(" ")
            cur[key] = value
        else:
            cur[line] = "1"
    _flush()
    return out


def render_worktrees_html(items: list[WorktreeEntry]) -> str:
    """Compact HTML inventory.  Empty input → helpful hint."""
    if not items:
        return (
            '<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
            'font-weight:700;">Git Worktrees</h4>'
            '<p style="margin:0;color:#9ca3af;font-size:11px;">'
            'No worktrees detected. Add one with '
            '<code>git worktree add &lt;path&gt; &lt;branch&gt;</code>.</p>'
        )

    rows: list[str] = []
    for entry in items:
        if entry.is_bare:
            badge = ('<span style="background:#9ca3af;color:#fff;font-size:10px;'
                     'padding:1px 6px;border-radius:8px;font-weight:700;">BARE</span>')
        elif entry.is_main:
            badge = ('<span style="background:#3b82f6;color:#fff;font-size:10px;'
                     'padding:1px 6px;border-radius:8px;font-weight:700;">MAIN</span>')
        else:
            badge = ('<span style="background:#f59e0b;color:#fff;font-size:10px;'
                     'padding:1px 6px;border-radius:8px;font-weight:700;">'
                     'WORKTREE</span>')
        branch = entry.branch or "(detached)"
        head = entry.head[:8] if entry.head else ""
        rows.append(
            '<li style="padding:6px 8px;border-bottom:1px solid #f3f4f6;'
            'font-size:11px;display:flex;gap:10px;align-items:center;">'
            f'{badge}'
            f'<code style="color:#111827;font-weight:600;">{_html.escape(branch)}</code>'
            f'<code style="color:#9ca3af;">{_html.escape(head)}</code>'
            f'<span style="color:#374151;flex:1;">'
            f'{_html.escape(entry.path)}</span>'
            '</li>'
        )

    return (
        f'<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
        f'font-weight:700;">Git Worktrees ({len(items)})</h4>'
        f'<div style="border:1px solid #e5e7eb;border-radius:8px;'
        f'background:#fff;padding:6px 4px;">'
        f'<ul style="list-style:none;margin:0;padding:0;">'
        + "".join(rows) +
        '</ul></div>'
    )


__all__ = ["WorktreeEntry", "list_worktrees", "render_worktrees_html"]
