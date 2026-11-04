"""A session in its own worktree and branch of the workspace's project.

Several sessions of one person (dashboard, terminal, messenger) in the
same agent workspace edit the same files when they share the folder. A
session started with ``worktree=True`` (``session.create``) runs in a
worktree of the workspace's repository instead, on a branch of its own,
and its work comes back by :func:`merge` when the user asks.

Input / output:

* :func:`start` -- ``(workspace, init_project=False)`` -> ``{"path",
  "branch", "base_ref", "repo"}``. A workspace that is not a git
  repository raises :class:`NeedsLocalProject` unless ``init_project`` is
  set, which is the user's consent given at the session start (the
  window asks; ``local_project.consent_question`` is the text); then
  ``local_project.init_local_project`` runs first, with its refusals
  (home, system paths, inside another repository).
* :func:`merge` -- ``(worktree_path)`` -> ``{"status": "ok"|"conflict"|
  "error", "files", "message", "target"}``: the worktree's committed and
  uncommitted changes since its branch point, applied to the main
  checkout's working tree, uncommitted, for review -- through the
  ``worktree_merge`` executor (target containment, the per-file write
  gate: protected files and files outside the write scope are refused,
  nobody being there to confirm them). A conflict leaves the main
  checkout untouched and the worktree kept; a clean merge removes the
  worktree and its branch.
* :func:`listing` -- the session worktrees of a repository.

Where the worktrees live: ``<repo>/.delfin/worktrees/<name>``, the layout
the dashboard's "Own worktree" option has used since it existed (moved
here from ``delfin/dashboard/agent_sessions.py``, which now calls
:func:`create`, so the dashboard and the service make the same tree), kept
out of ``git status`` by an entry in the clone's ``info/exclude``, with the
sidecar ``.delfin/session_worktree.json`` inside the tree that
``agent_sessions.release_session_worktree`` reads. Inside the person's own
workspace: per person, covered by the home directory's backup, not
readable through /tmp. Rejected: a sibling ``<repo>.worktrees`` directory
(a second layout beside the one the dashboard, ``ops_status`` and the
git-errors skill already know).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

from . import local_project as _lp
from . import worktree as _wt

BRANCH_PREFIX = "session"


class NeedsLocalProject(RuntimeError):
    """The workspace is not a git repository and no consent was given."""

    def __init__(self, workspace: Path, question: str):
        super().__init__(
            f"{workspace} is not a git repository. A session in its own "
            "worktree needs one: start it again with the user's consent to "
            "make the folder a local project (initProject).")
        self.workspace = workspace
        self.question = question


#: Beside the worktree's own git metadata, not in a central registry: the
#: record and the thing it describes are then removed by the same act, so
#: there is no list that can outlive its entries.
WORKTREE_SIDECAR = ".delfin/session_worktree.json"


def exclude_locally(root: str, common_dir: str, pattern: str = ".delfin/") -> None:
    """Keep ``pattern`` out of ``git status`` in this clone only: an entry
    in its info/exclude, nothing committed. Skipped when already ignored."""
    try:
        probe = _wt.safe_git.run(
            ["check-ignore", "-q", ".delfin/worktrees/probe"], root, timeout=5)
        if probe.returncode == 0:
            return
        exclude = Path(common_dir) / "info" / "exclude"
        exclude.parent.mkdir(parents=True, exist_ok=True)
        text = exclude.read_text(encoding="utf-8") if exclude.exists() else ""
        if pattern not in text.splitlines():
            sep = "" if not text or text.endswith("\n") else "\n"
            exclude.write_text(f"{text}{sep}{pattern}\n", encoding="utf-8")
    except Exception:
        pass


def write_worktree_sidecar(info) -> None:
    """Record repo, branch and owner for a session worktree. Never raises.

    pid and host are recorded so a later reader can tell a tree whose
    session is gone from one that is still being worked in, without
    judging another node's pid."""
    import socket
    from . import proc_identity as _proc_identity
    try:
        path = Path(info.path) / WORKTREE_SIDECAR
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({
            "repo_dir": str(info.repo_dir),
            "path": str(info.path),
            "branch": str(info.branch),
            "base_ref": str(info.base_ref),
            "created_at": float(info.created_at),
            "pid": os.getpid(),
            # With the pid: a number handed out again later names a
            # different process (proc_identity).
            "proc_start": _proc_identity.process_start(os.getpid()),
            "host": (socket.gethostname() or "").strip(),
        }), encoding="utf-8")
    except Exception:
        pass


def create(workspace) -> tuple[str, "_wt.WorktreeInfo"]:
    """A fresh worktree of the repository *workspace* is in, on a
    ``session/<hex>`` branch, under ``<repo>/.delfin/worktrees``. Returns
    the path at the same place inside the tree as *workspace* is inside
    the repository, and the WorktreeInfo. Raises ValueError outside a
    repository, WorktreeError when git refuses."""
    from . import session_presence as _presence
    repo = _presence.repository_of(str(workspace))
    if not repo["root"]:
        raise ValueError(f"{workspace} is not in a git repository")
    exclude_locally(repo["root"], repo["common_dir"])
    info = _wt.enter_worktree(repo["root"], branch_prefix=BRANCH_PREFIX,
                              parent=Path(repo["root"]) / ".delfin" / "worktrees")
    write_worktree_sidecar(info)
    try:
        inside = Path(workspace).resolve().relative_to(Path(repo["root"]))
    except ValueError:
        inside = Path(".")
    return str((info.path / inside).resolve()), info


def start(workspace, *, init_project: bool = False) -> dict:
    """A worktree + branch of *workspace*'s repository for one session."""
    ws = Path(workspace).expanduser().resolve()
    if not _wt._is_git_repo(ws):
        why = _lp.refusal(ws)
        if why:
            raise _wt.WorktreeError(
                f"{ws} is not a git repository and cannot become a local "
                f"project: {why}")
        if not init_project:
            raise NeedsLocalProject(ws, _lp.consent_question(_lp.plan(ws)))
        try:
            _lp.init_local_project(ws)
        except _lp.LocalProjectError as exc:
            raise _wt.WorktreeError(
                f"making {ws} a local project failed and was rolled back: "
                f"{exc}") from None
        try:
            from . import session_presence as _presence
            _presence._git_cache.pop(str(ws), None)
        except Exception:
            pass
    path, info = create(ws)
    return {"path": path, "worktree": str(info.path), "branch": info.branch,
            "base_ref": info.base_ref, "repo": str(info.repo_dir)}


def listing(workspace) -> list[dict]:
    """The session worktrees (branch ``session/*``) of *workspace*'s
    repository: ``[{"path", "branch", "head"}]``; [] outside a repository."""
    ws = Path(workspace).expanduser()
    r = _wt.safe_git.run(["worktree", "list", "--porcelain"], ws, timeout=15)
    if r.returncode != 0:
        return []
    out, cur = [], {}
    for line in r.stdout.splitlines() + [""]:
        if not line.strip():
            if cur.get("branch", "").startswith(f"{BRANCH_PREFIX}/"):
                out.append(cur)
            cur = {}
        elif line.startswith("worktree "):
            cur["path"] = line[len("worktree "):]
        elif line.startswith("HEAD "):
            cur["head"] = line[len("HEAD "):]
        elif line.startswith("branch "):
            cur["branch"] = line[len("branch "):].removeprefix("refs/heads/")
    return out


def merge(worktree_path, *, target_dir: Optional[str] = None) -> dict:
    """Merge a session worktree back into its repository's main checkout.

    Runs the ``worktree_merge`` executor with a permission object rooted at
    the main checkout, mode ``acceptEdits`` and no one to ask: the user's
    request is the approval of the merge, and every file the write gate
    would ask about (the agent's safety layer, git configuration files,
    stored calculations) is refused rather than written unseen."""
    from .api_client import KitToolPermissions, _DocToolExecutor
    wt_path = Path(worktree_path).expanduser().resolve()
    r = _wt.safe_git.run(["worktree", "list", "--porcelain"], wt_path,
                         timeout=15)
    if r.returncode != 0:
        return {"status": "error", "files": [], "target": "",
                "message": f"not a worktree: {wt_path}"}
    main = next((line[len("worktree "):] for line in r.stdout.splitlines()
                 if line.startswith("worktree ")), "")
    target = Path(target_dir).expanduser().resolve() if target_dir \
        else Path(main)
    perms = KitToolPermissions(workspace=target, mode="acceptEdits")
    perms.task_session_id = "session-worktree-merge"
    raw = _DocToolExecutor()._execute_worktree_merge(
        {"path": str(wt_path), "target_dir": str(target)}, perms,
        framework_worktree=True)
    try:
        res = json.loads(raw)
    except ValueError:
        res = {"error": raw}
    if res.get("error"):
        return {"status": "error", "files": [], "target": str(target),
                "message": str(res["error"])}
    return {"status": str(res.get("status") or ""),
            "files": list(res.get("files") or []),
            "target": str(target),
            "message": str(res.get("message") or "")}


__all__ = ["BRANCH_PREFIX", "WORKTREE_SIDECAR", "NeedsLocalProject", "create",
           "start", "merge", "listing", "exclude_locally",
           "write_worktree_sidecar"]
