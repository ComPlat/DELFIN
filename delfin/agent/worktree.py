"""Git worktree isolation (.delfin-native enter/exit).

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
from dataclasses import dataclass, field
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
    # Background jobs that were still running inside the tree when a
    # teardown was attempted; non-empty means the tree was NOT removed.
    held_by_jobs: list = field(default_factory=list)


def jobs_holding_worktree(path: Path | str) -> list:
    """Background jobs still running inside ``path``.

    ``bash_background`` registers a job under the workspace it was started
    in, and for an isolated sub-agent that workspace IS the worktree —
    including the ``.delfin/bash_jobs.json`` that is the only record of the
    job. The child is started in its own session, so ``git worktree remove
    --force`` does not stop it: it keeps running with its working directory
    deleted underneath it, its completion event destroyed with the registry
    file, and ``bash_status`` answering "unknown job_id", which reads as a
    typo rather than as a live process holding a shared node.

    Best-effort: a registry that cannot be read reports nothing, so a
    teardown behaves exactly as it did before this check existed.
    """
    try:
        from .bash_jobs import running_jobs_for_workspace
        return running_jobs_for_workspace(path) or []
    except Exception:
        return []


def _held_message(path: Path | str, jobs: list) -> str:
    ids = ", ".join(str((j or {}).get("job_id") or "?") for j in jobs[:8])
    return (f"the worktree was NOT removed: {len(jobs)} background job(s) "
            f"started in it are still running ({ids}). Removing it would "
            f"delete their working directory and their completion records "
            f"while the processes keep running. End them with "
            f"bash_kill(job_id) or wait for them; the tree stays at {path}.")


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
    """Tear down the worktree. Keep it iff there are changes — or iff
    background jobs are still running inside it (see
    :func:`jobs_holding_worktree`); that check lives here rather than in
    one caller so every teardown path is covered by it."""
    if info.cleaned_up:
        return info
    info.held_by_jobs = jobs_holding_worktree(info.path)
    if info.held_by_jobs:
        info.had_changes = _has_changes(info)
        info.final_path = info.path
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


@dataclass
class MergeResult:
    ok: bool                 # operation completed without error
    applied: bool            # changes were actually written to the target tree
    files: list              # changed paths (relative to the worktree)
    message: str             # human/agent-facing summary


def _changed_files(wt: Path, base_ref: str) -> list:
    try:
        out = _run_git(wt, "diff", "--cached", "--name-only", base_ref)
    except WorktreeError:
        return []
    return [ln for ln in out.splitlines() if ln.strip()]


def _git_apply(repo: Path, patch: str, *, check: bool) -> tuple[bool, str]:
    """Run `git apply` reading the patch from stdin. Returns (ok, stderr)."""
    args = ["git", "apply", "--whitespace=nowarn"]
    if check:
        args.append("--check")
    args.append("-")
    try:
        proc = subprocess.run(
            args, cwd=str(repo), input=patch,
            capture_output=True, text=True, timeout=60,
        )
    except FileNotFoundError as exc:
        raise WorktreeError("git is not installed") from exc
    return proc.returncode == 0, proc.stderr.strip()


def _cleanup_merged_worktree(info: WorktreeInfo) -> None:
    """Remove a worktree + its throwaway ``agent/*`` branch after a successful
    merge (the changes now live in the target tree, so both are redundant).
    Best-effort — never raises; the worktree must be removed BEFORE the branch
    so git isn't holding the branch checked out. Prevents parallel-writer
    fan-outs from leaking ``/tmp`` worktrees and orphan branches.

    A merge does not stop the jobs the worktree is running, so the same
    guard applies here: the changes are already in the target tree, and the
    tree can be removed once the processes using it are done."""
    wt = info.final_path or info.path
    info.held_by_jobs = jobs_holding_worktree(wt)
    if info.held_by_jobs:
        info.final_path = Path(wt)
        return
    try:
        _run_git(info.repo_dir, "worktree", "remove", "--force", str(wt))
    except WorktreeError:
        if wt.exists():
            shutil.rmtree(wt, ignore_errors=True)
    try:
        _run_git(info.repo_dir, "branch", "-D", info.branch)
    except WorktreeError:
        pass
    info.cleaned_up = True
    info.final_path = None


def merge_worktree(
    info: WorktreeInfo,
    *,
    target_repo: Path | str | None = None,
    cleanup: bool = True,
) -> MergeResult:
    """Bring a worktree's changes (vs its ``base_ref``) into the target repo's
    working tree — but ONLY if they apply cleanly. Never forces, never leaves a
    half-applied tree: on any conflict the target is left untouched and the
    worktree is kept for manual merge.

    Completes the fan-out-writers → review (``diff_summary``) → merge flow.
    ``target_repo`` defaults to ``info.repo_dir`` (the source). Changes land in
    the working tree *uncommitted* so the parent can review and commit them.
    """
    repo = Path(target_repo).resolve() if target_repo else info.repo_dir
    wt = info.final_path or info.path
    if not _is_git_repo(repo):
        raise WorktreeError(f"not a git repo: {repo}")
    if not wt.exists():
        raise WorktreeError(f"worktree path missing: {wt}")
    # Stage everything (including untracked) so the patch captures the full
    # worktree state, then diff against the shared branch point.
    _run_git(wt, "add", "-A")
    patch = _run_git(wt, "diff", "--cached", "--binary", info.base_ref)
    if not patch.strip():
        return MergeResult(True, False, [], "no changes to merge")
    files = _changed_files(wt, info.base_ref)
    ok, err = _git_apply(repo, patch, check=True)
    if not ok:
        return MergeResult(
            False, False, files,
            f"changes do not apply cleanly to {repo} (conflicts) — target tree "
            f"left untouched; worktree kept at {wt} for manual merge. "
            f"git: {err}"[:600],
        )
    applied_ok, err2 = _git_apply(repo, patch, check=False)
    if not applied_ok:  # passed --check but failed for real → nothing partial trusted
        return MergeResult(
            False, False, files,
            f"merge failed mid-apply into {repo}: {err2}"[:600],
        )
    # Success: the changes now live in the target working tree, so the
    # isolated worktree + its throwaway branch have served their purpose.
    if cleanup:
        _cleanup_merged_worktree(info)
    message = (f"merged {len(files)} file(s) into the working tree of {repo}; "
               f"review with `git -C {repo} diff` and commit.")
    if info.held_by_jobs:
        message += " " + _held_message(info.final_path or info.path,
                                       info.held_by_jobs)
    return MergeResult(True, True, files, message)


def diff_summary(info: WorktreeInfo, *, max_chars: int = 2000) -> str:
    """Concise review of what a worktree changed vs its base: the ``--stat`` of
    tracked changes plus any new untracked files. Empty when the worktree is
    gone or unchanged. Lets the parent agent SEE what an isolated (e.g.
    parallel-writer) subagent actually did, not just that it changed something."""
    path = info.final_path or info.path
    try:
        stat = _run_git(path, "diff", "--stat", info.base_ref)
        untracked = _run_git(path, "ls-files", "--others", "--exclude-standard")
    except WorktreeError:
        return ""
    parts: list[str] = []
    if stat.strip():
        parts.append(stat.strip())
    if untracked.strip():
        parts.append("Untracked (new):\n" + "\n".join(
            f"  {f}" for f in untracked.strip().splitlines()))
    return "\n".join(parts).strip()[:max_chars]


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
    "MergeResult",
    "enter_worktree",
    "exit_worktree",
    "jobs_holding_worktree",
    "merge_worktree",
    "diff_summary",
    "worktree_session",
]
