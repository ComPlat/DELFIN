"""Git worktree isolation (Claude-Code-compatible enter/exit).

Spawning a sub-task inside a temporary ``git worktree`` lets the
agent edit, build, and test against a side-branch without
disturbing the user's working tree. On exit:

  - if the agent made no commits, the worktree is removed and the
    branch is deleted;
  - if commits exist, the worktree's path and branch are returned
    so the user can review / merge.

This module exposes a small synchronous API plus a context manager
for tests::

    info = enter_worktree(repo_dir)
    ...
    exit_worktree(info, keep_if_changed=True)

    with worktree_session(repo_dir) as info:
        ...

Hard requirements: ``git`` on PATH, the source path is a git repo.
Failures raise ``WorktreeError`` rather than crashing the agent.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional


class WorktreeError(RuntimeError):
    pass


@dataclass
class WorktreeInfo:
    repo_dir: Path        # source repo
    path: Path            # worktree directory (under tmp by default)
    branch: str
    base_ref: str
    created_at: float
    cleaned_up: bool = False
    final_path: Optional[Path] = None    # set if kept after exit
    had_changes: bool = False


def _run_git(repo: Path, *args: str) -> str:
    try:
        result = subprocess.run(
            ["git", *args], cwd=str(repo),
            capture_output=True, text=True, timeout=30,
        )
    except FileNotFoundError as exc:
        raise WorktreeError("git is not installed") from exc
    if result.returncode != 0:
        raise WorktreeError(
            f"git {' '.join(args)} failed: {result.stderr.strip()}"
        )
    return result.stdout


def _is_git_repo(path: Path) -> bool:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--git-dir"], cwd=str(path),
            capture_output=True, text=True, timeout=10,
        )
        return out.returncode == 0
    except (FileNotFoundError, subprocess.SubprocessError):
        return False


def _current_head(repo: Path) -> str:
    return _run_git(repo, "rev-parse", "HEAD").strip()


def enter_worktree(
    repo_dir: Path | str,
    *,
    branch_prefix: str = "agent",
    parent: Path | str | None = None,
) -> WorktreeInfo:
    """Create a fresh worktree for the agent. Returns WorktreeInfo."""
    repo = Path(repo_dir).resolve()
    if not _is_git_repo(repo):
        raise WorktreeError(f"not a git repo: {repo}")

    base_ref = _current_head(repo)
    parent_dir = Path(parent) if parent else Path(tempfile.gettempdir())
    parent_dir.mkdir(parents=True, exist_ok=True)
    suffix = uuid.uuid4().hex[:8]
    branch_name = f"{branch_prefix}/{suffix}"
    wt_path = parent_dir / f"delfin-wt-{suffix}"
    if wt_path.exists():
        raise WorktreeError(f"worktree path collision: {wt_path}")
    _run_git(
        repo, "worktree", "add", "--quiet",
        "-b", branch_name, str(wt_path), base_ref,
    )
    return WorktreeInfo(
        repo_dir=repo,
        path=wt_path,
        branch=branch_name,
        base_ref=base_ref,
        created_at=time.time(),
    )


def _has_changes(info: WorktreeInfo) -> bool:
    """True if the worktree differs from base_ref (commits or dirty files)."""
    try:
        head_now = _run_git(info.path, "rev-parse", "HEAD").strip()
        if head_now != info.base_ref:
            return True
        status = _run_git(info.path, "status", "--porcelain").strip()
        return bool(status)
    except WorktreeError:
        return False


def exit_worktree(
    info: WorktreeInfo,
    *,
    keep_if_changed: bool = True,
) -> WorktreeInfo:
    """Tear down the worktree. Keep it iff there are changes."""
    if info.cleaned_up:
        return info
    info.had_changes = _has_changes(info)
    if info.had_changes and keep_if_changed:
        info.final_path = info.path
        info.cleaned_up = True
        return info
    # Otherwise remove the worktree and (if no commits exist on it)
    # delete the branch so we don't leave clutter behind.
    try:
        _run_git(info.repo_dir, "worktree", "remove", "--force", str(info.path))
    except WorktreeError:
        # Best-effort cleanup if git refuses for some reason
        if info.path.exists():
            shutil.rmtree(info.path, ignore_errors=True)
    if not info.had_changes:
        try:
            _run_git(info.repo_dir, "branch", "-D", info.branch)
        except WorktreeError:
            pass
    info.cleaned_up = True
    return info


@contextmanager
def worktree_session(
    repo_dir: Path | str,
    *,
    branch_prefix: str = "agent",
    keep_if_changed: bool = True,
) -> Iterator[WorktreeInfo]:
    info = enter_worktree(repo_dir, branch_prefix=branch_prefix)
    try:
        yield info
    finally:
        exit_worktree(info, keep_if_changed=keep_if_changed)


__all__ = [
    "WorktreeError",
    "WorktreeInfo",
    "enter_worktree",
    "exit_worktree",
    "worktree_session",
]
