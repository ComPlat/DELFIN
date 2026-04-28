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


def _classify_worktree(entry: WorktreeEntry) -> str:
    """Bucket a worktree path into one of four groups.

    - ``main``       — the repo's primary worktree(s). Always shown.
    - ``commit_sweep`` — created by the agent_workspace/commit_sweep tool.
    - ``bisect``     — created by ad-hoc bisect helpers. Path patterns:
                       /tmp/bvs_*, /tmp/vf_*, /tmp/vsf_*, /tmp/vs51_*,
                       /tmp/sweep_wt_*, /tmp/sd_*, /tmp/sw_*, /tmp/r20*,
                       /tmp/feschiff_*, /tmp/main_wt, /tmp/wt_pre_*.
    - ``other``      — anything else (user-created worktrees, etc.).
    """
    if entry.is_main or entry.is_bare:
        return "main"
    p = entry.path or ""
    if "/agent_workspace/commit_sweep/" in p:
        return "commit_sweep"
    bisect_prefixes = (
        "/tmp/bvs_", "/tmp/vf_", "/tmp/vsf_", "/tmp/vs51_",
        "/tmp/sweep_wt_", "/tmp/sd_", "/tmp/sw_", "/tmp/r20",
        "/tmp/feschiff_", "/tmp/main_wt", "/tmp/wt_pre_",
        "/tmp/delfin_bisect", "/tmp/delfin-main",
    )
    if any(p.startswith(pref) for pref in bisect_prefixes):
        return "bisect"
    return "other"


def _group_worktrees(
    items: list[WorktreeEntry],
) -> dict[str, list[WorktreeEntry]]:
    """Split a worktree list into ``{main, commit_sweep, bisect, other}``."""
    groups: dict[str, list[WorktreeEntry]] = {
        "main": [], "commit_sweep": [], "bisect": [], "other": [],
    }
    for entry in items:
        groups[_classify_worktree(entry)].append(entry)
    return groups


def _render_worktree_row(entry: WorktreeEntry) -> str:
    """Render one worktree as a compact list item."""
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
    return (
        '<li style="padding:6px 8px;border-bottom:1px solid #f3f4f6;'
        'font-size:11px;display:flex;gap:10px;align-items:center;">'
        f'{badge}'
        f'<code style="color:#111827;font-weight:600;">{_html.escape(branch)}</code>'
        f'<code style="color:#9ca3af;">{_html.escape(head)}</code>'
        f'<span style="color:#374151;flex:1;">'
        f'{_html.escape(entry.path)}</span>'
        '</li>'
    )


def _render_group(label: str, items: list[WorktreeEntry], *, open_default: bool) -> str:
    """Render one group as a <details> element. Empty groups are hidden."""
    n = len(items)
    if n == 0:
        return ""
    rows = "".join(_render_worktree_row(e) for e in items)
    open_attr = " open" if open_default else ""
    summary_color = "#111827" if open_default else "#374151"
    return (
        f'<details{open_attr} style="border:1px solid #e5e7eb;border-radius:8px;'
        f'background:#fff;margin-bottom:8px;">'
        f'<summary style="padding:8px 12px;cursor:pointer;font-size:12px;'
        f'font-weight:600;color:{summary_color};">'
        f'{_html.escape(label)} <span style="color:#6b7280;font-weight:500;">'
        f'({n})</span></summary>'
        f'<ul style="list-style:none;margin:0;padding:6px 4px;">{rows}</ul>'
        f'</details>'
    )


def render_worktrees_html(items: list[WorktreeEntry]) -> str:
    """Grouped + collapsible HTML inventory.

    - Main worktrees expanded by default (always small, always relevant).
    - External-tool buckets (commit_sweep, bisect) collapsed because
      they can each balloon to hundreds of entries from sweep runs.
    - "Other" expanded when present (user-created worktrees we didn't
      classify — likely worth a glance).
    - Empty input → helpful hint.
    """
    if not items:
        return (
            '<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
            'font-weight:700;">Git Worktrees</h4>'
            '<p style="margin:0;color:#9ca3af;font-size:11px;">'
            'No worktrees detected. Add one with '
            '<code>git worktree add &lt;path&gt; &lt;branch&gt;</code>.</p>'
        )

    groups = _group_worktrees(items)
    n_main = len(groups["main"])
    n_sweep = len(groups["commit_sweep"])
    n_bisect = len(groups["bisect"])
    n_other = len(groups["other"])

    # Header summary so the user sees the breakdown without expanding
    summary_parts = [f"{n_main} main"]
    if n_sweep:
        summary_parts.append(f"{n_sweep} commit_sweep")
    if n_bisect:
        summary_parts.append(f"{n_bisect} bisect")
    if n_other:
        summary_parts.append(f"{n_other} other")
    breakdown = " · ".join(summary_parts)

    sections = "".join([
        _render_group("Main worktrees", groups["main"], open_default=True),
        _render_group("External: commit_sweep", groups["commit_sweep"], open_default=False),
        _render_group("External: bisect", groups["bisect"], open_default=False),
        _render_group("External: other", groups["other"], open_default=True),
    ])

    return (
        f'<h4 style="margin:12px 0 6px 0;color:#111827;font-size:13px;'
        f'font-weight:700;">Git Worktrees ({len(items)}) '
        f'<span style="color:#6b7280;font-weight:500;font-size:11px;">'
        f'— {_html.escape(breakdown)}</span></h4>'
        f'{sections}'
    )


__all__ = [
    "WorktreeEntry", "list_worktrees", "render_worktrees_html",
    "_classify_worktree", "_group_worktrees",
]
