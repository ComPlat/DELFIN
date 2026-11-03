"""Branch operations for every surface: the dashboard's logic, UI-free.

The dashboard's branch widgets (``delfin/dashboard/branch_panel.py``,
``branch_list.py``, ``branch_preview.py``, ``branch_cleanup.py``) hold
their logic in modules with no widget imports already. This module is the
agent side's single door to them: it imports that logic -- it does not
copy it -- and answers the questions the Pi's surfaces ask in plain
values:

* :func:`rows` -- the dropdown's rows of one view (``mine``, ``all``,
  ``hidden``, ``merged``, ``remote``), as dicts. The dashboard calls
  ``branch_list.all_rows`` and ``view_rows`` directly; both answer the
  same question here, so there is one implementation of "what branches
  exist and what each is".
* :func:`pull_main` -- fast-forward-only update of ``main``. Extracted
  from the dashboard's header (``delfin/dashboard/__init__.py``), which
  now calls this back; its lock helpers moved to
  ``delfin/dashboard/git_write.py`` so the dashboard and this module take
  the same lock.
* :func:`test_branch` -- a worktree at the branch (the preview's
  ``ensure_worktree``; the worktree never holds the shared checkout's
  branch when the checkout itself is on it), the branch's changed test
  files run there (capped wall clock, ``nice``), pass/fail counts and
  the first failure. The worktree is removed afterwards unless kept.
* :func:`cleanup_candidates` / :func:`cleanup_delete` -- what
  ``branch_cleanup`` offers and deletes. The deletion is
  ``branch_cleanup.delete`` itself, with its rules re-read at deletion
  time; a branch with commits not in the base needs the second
  confirmation, and nothing is deleted without ``confirm=True``.

Nothing here switches the shared checkout: a pull only fast-forwards the
branch the checkout is on when that is ``main``; a test runs in a second
worktree; a deletion refuses the checked-out branch.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Callable, Optional

#: Views every surface offers, in this order.
VIEWS = ("mine", "all", "hidden", "merged", "remote")

#: The wall-clock cap of one test run in a test worktree. A branch's own
#: test files should not take a work day; a run cut off at this says so.
TEST_TIMEOUT_S = 900

#: The most test files one run gets: a branch that touches more has its
#: runs capped by time anyway, and one pytest call per file would not.
MAX_TEST_FILES = 20


def _repo_of(repo) -> Path:
    return Path(repo).expanduser().resolve()


def _is_git(repo) -> bool:
    from delfin.dashboard import branch_list as BL
    return BL._git(_repo_of(repo), "rev-parse", "--git-dir").returncode == 0


def _not_git(text: str = "not a git checkout") -> dict:
    return {"ok": False, "text": text}


# -- rows (the dashboard's branch_list, imported, not copied) -----------------

def rows(repo, view: str = "mine", text: str = "", *, gh="auto",
         refresh: bool = False):
    """The rows of *view* (one of :data:`VIEWS`), narrowed by *text*, as
    plain dicts. A non-git checkout is not an empty list -- that would
    read as "no branches" -- it is the refusal dict
    ``{"ok": False, "text": ..., "rows": []}``."""
    from delfin.dashboard import branch_list as BL
    if not _is_git(repo):
        return {"ok": False, "text": "not a git checkout", "rows": []}
    all_rows = BL.all_rows(_repo_of(repo), gh, refresh=refresh)
    picked = BL.view_rows(all_rows, view if view in VIEWS else "mine", text)
    return {
        "ok": True,
        "rows": [{
        "ref": r.ref, "kind": r.kind, "sha": r.sha, "date": r.date,
        "author": r.author, "merged": r.merged, "pinned": r.pinned,
        "hidden": r.hidden, "mine": r.mine, "ahead": r.ahead,
        "pr": r.pr, "ci": r.ci, "current": r.current,
    } for r in picked]}


# -- pull (extracted from the dashboard's header handler) ---------------------

def pull_main(repo, *, read=None, write=None) -> dict:
    """Fast-forward-only update of ``main`` of the checkout the service
    runs in. The shared checkout is never switched, never merged:
    off-``main`` and a diverged ``main`` are refused with a reason.

    ``read``/``write`` stand in for the dashboard's capture/write runners
    in tests; the default is :mod:`delfin.dashboard.git_write`."""
    from delfin.dashboard import git_write as GW
    read = read or (lambda r, *a: GW.run_git_capture(r, *a))
    write = write or GW.run_git_write
    repo = _repo_of(repo)
    if not _is_git(repo):
        return _not_git()
    try:
        with GW.dashboard_git_write_lock(repo, "git pull"):
            if GW.ensure_git_index_lock_clear(repo, "Pull"):
                pass  # a stale lock removed; the pull below decides the rest
            from delfin.agent import held_branch as HB
            moved = HB.update_detached(repo)
            if moved is not None:
                return {"ok": True, "text": moved.message.splitlines()[0][:160]}
            current = (read(repo, "branch", "--show-current").stdout
                       or "").strip()
            if current != "main":
                return {"ok": False, "text": (
                    "Pull only updates main; the checkout is on "
                    f"{current or 'a detached commit'}. Commit or move its "
                    "changes, then run git switch main in a terminal. To "
                    "try a branch, use Test.")}
            before = (read(repo, "rev-parse", "HEAD").stdout or "").strip()
            result = write(repo, "pull", "--ff-only", "origin", "main")
            if result.returncode != 0:
                return {"ok": False, "text": (
                    (result.stdout or result.stderr or
                     "git pull failed").strip()[:400] or
                    "Couldn't update main (fast-forward only).")}
            n = (read(repo, "rev-list", "--count",
                      f"{before}..HEAD").stdout or "").strip()
            return {"ok": True, "text": (
                "main is up to date." if n in ("", "0")
                else f"Updated main: {n} new commit(s).")}
    except RuntimeError as exc:
        return {"ok": False, "text": str(exc)}


# -- test (a worktree, the branch's changed tests) ----------------------------

def _changed_test_files(worktree: Path, merge_merge_base_sha: str) -> list[str]:
    """Test files the branch changes against its merge base with the base
    ref -- plus the tests of files a branch may not touch directly is not
    guessed here: only what git names."""
    from delfin.dashboard import branch_list as BL
    base = BL.base_ref(worktree)
    if not base:
        base = "HEAD~1"
    diff = BL._git(worktree, "diff", "--name-only", "--diff-filter=d",
                   f"{merge_merge_base_sha}", "--", "tests/")
    files = [line.strip() for line in (diff.stdout or "").splitlines()
             if line.strip().endswith(".py")]
    return files[:MAX_TEST_FILES]


def _merge_base_sha(repo) -> str:
    from delfin.dashboard import branch_list as BL
    base = BL.base_ref(repo)
    sha = BL._out(repo, "rev-parse", "--verify", f"{base}") if base else ""
    return sha or (BL._out(repo, "rev-parse", "HEAD") or "")


def run_tests(worktree, files, timeout_s=TEST_TIMEOUT_S) -> dict:
    """The default runner: pytest over *files* in *worktree* (nice, capped).

    Answers ``{"passed": int, "failed": int, "first_failure": str,
    "status": str}``. A capped run says ``status == "timeout"`` and
    reports nothing about the suite."""
    from delfin.agent import test_runner
    target = " ".join(files) if files else ""

    def nice(argv, report_dir):
        # nice, not renice: the branch's tests must not slow the service
        # or a running calculation on the same machine.
        nice_bin = "nice" if os.path.exists("/usr/bin/nice") else "nice"
        return [nice_bin, "-n", "10"] + list(argv)

    result = test_runner.run_tests(
        worktree, target=target, timeout_s=timeout_s, wrap=nice)
    summary = result.get("summary") or {}
    failures = result.get("failures") or []
    first = ""
    if failures:
        row = failures[0]
        first = f"{row.get('nodeid', '')}: {(row.get('message') or '')[:300]}"
    return {
        "status": result.get("status", "error"),
        "passed": int(summary.get("passed", 0) or 0),
        "failed": int(summary.get("failed", 0) or 0)
        + int(summary.get("errors", 0) or 0),
        "first_failure": first,
        "note": result.get("note", ""),
    }


def test_branch(repo, ref: str, *, keep: bool = False,
                runner: Optional[Callable] = None) -> dict:
    """A worktree at *ref*, the branch's changed test files run there.

    Nothing switches the shared checkout; the worktree is removed again
    unless *keep*. The runner stands in for :func:`run_tests` in tests."""
    from delfin.dashboard import branch_preview as BP
    from delfin.dashboard import branch_list as BL
    runner = runner or run_tests
    repo = _repo_of(repo)
    if not _is_git(repo):
        return _not_git("not a git checkout: nothing to test")
    target, why = BP.resolve_target(repo, ref)
    if target is None:
        return {"ok": False, "text": why}
    ok, message, worktree, attached = BP.ensure_worktree(repo, target)
    if not ok:
        return {"ok": False, "text": message}
    files = _changed_test_files(worktree, _merge_base_sha(repo))
    result = runner(worktree, tuple(files), TEST_TIMEOUT_S)
    text_parts = [f"{target.ref}: {result.get('passed', 0)} passed, "
                  f"{result.get('failed', 0)} failed"]
    if result.get("status") == "timeout":
        text_parts.append(f"cut off at {TEST_TIMEOUT_S}s")
    if result.get("first_failure"):
        text_parts.append(f"first failure: {result['first_failure']}")
    if result.get("note"):
        text_parts.append(str(result["note"]))
    if not files and result.get("passed", 0) == 0 and \
            result.get("failed", 0) == 0:
        text_parts.append("no test files changed on this branch")
    kept = bool(keep)
    if not kept:
        removed = remove_test_worktree(repo, worktree)
        if not removed:
            text_parts.append("the worktree could not be removed; "
                              f"it stays at {worktree}")
            kept = True
    return {"ok": True, "text": "; ".join(text_parts),
            "passed": result.get("passed", 0),
            "failed": result.get("failed", 0),
            "first_failure": result.get("first_failure", ""),
            "worktree": str(worktree), "kept": kept}


def remove_test_worktree(repo, worktree) -> bool:
    """Remove a test worktree created by :func:`test_branch`. Uncommitted
    changes keep it (the worktree may hold work); the branch is never
    touched. The caches a test run leaves (``__pycache__``,
    ``.pytest_cache``) are not work: they are removed first, so a run
    that changed nothing does not keep the worktree alive."""
    from delfin.dashboard import branch_preview as BP
    path = Path(worktree)
    if not path.is_dir():
        return True
    import shutil
    for cache in path.rglob("__pycache__"):
        shutil.rmtree(cache, ignore_errors=True)
    for cache in (".pytest_cache", ".ruff_cache", ".mypy_cache"):
        shutil.rmtree(path / cache, ignore_errors=True)
    if BP.changes(path):
        return False
    from delfin.dashboard import branch_list as BL
    r = BL._git(_repo_of(repo), "worktree", "remove", str(path))
    if r.returncode != 0:
        r = BL._git(_repo_of(repo), "worktree", "remove", "--force",
                    str(path))
    if r.returncode == 0:
        BL._git(_repo_of(repo), "worktree", "prune", "--expire=now")
        return True
    return False


# -- cleanup (branch_cleanup.delete itself; the rules live there) -------------

def cleanup_candidates(repo, *, include_unmerged: bool = False,
                       refresh: bool = False):
    """What may be deleted (``candidates``) and what is kept (``kept``),
    as dicts. ``branch_cleanup.candidates`` answers the same question for
    the dashboard; this is the same call. A non-git workspace answers
    ``None, None`` -- the caller carries the refusal (an empty list would
    read as "nothing to clean up")."""
    from delfin.dashboard import branch_cleanup as BC
    if not _is_git(repo):
        return None, None
    cands, kept = BC.candidates(_repo_of(repo), refresh=refresh,
                                include_unmerged=include_unmerged)
    return ([_candidate_dict(c) for c in cands],
            [{"name": k.name, "reason": k.reason} for k in kept])


def _candidate_dict(c) -> dict:
    return {"kind": c.kind, "name": c.name, "sha": c.sha, "reason": c.reason,
            "path": c.path, "branch": c.branch, "prunable": c.prunable,
            "unmerged": c.unmerged, "after": c.after,
            "upstream": c.upstream}


def cleanup_delete(repo, names, *, confirm: bool = False,
                   confirm_unmerged: Optional[bool] = None) -> dict:
    """Delete *names* through ``branch_cleanup.delete`` -- the dashboard's
    own deletion, with its rules re-read at deletion time (refusals,
    moved-branch check, second confirmation for unmerged branches).

    Nothing is deleted without ``confirm=True``: a listing may be answered
    from an earlier read, a deletion is a deletion."""
    from delfin.dashboard import branch_cleanup as BC
    repo = _repo_of(repo)
    if not _is_git(repo):
        return _not_git()
    names = [n for n in names if str(n).strip()]
    if not names:
        return {"ok": False, "text": "nothing named to delete"}
    if not confirm:
        return {"ok": False, "text": (
            "Not deleted: deleting needs the confirmation "
            "(cleanup ok <number>, or confirm=true). Nothing was changed.")}
    refusals = BC.delete_refusals(repo, names)
    candidates_by_name = {}
    for name in names:
        why = refusals.get(name)
        if why:
            return {"ok": False, "text": f"{name} was not deleted: {why}."}
        sha = BC.BL._out(repo, "rev-parse", f"refs/heads/{name}") \
            if BC.BL.ref_kind(repo, name) == "local" else ""
        if not sha:
            return {"ok": False,
                    "text": f"{name} was not deleted: it is not a local "
                            "branch."}
        candidates_by_name[name] = BC.Candidate(
            "branch", name, sha, "picked in a branch list",
            unmerged=BC.unmerged_count(repo, name))
    # The second confirmation rides the same boolean the dashboard's
    # dialog sets after the user has seen the count: one word (ok) after
    # the listing names it on the phone.
    second = (not confirm_unmerged) if confirm_unmerged is not None \
        else confirm
    report = BC.delete(repo, list(candidates_by_name.values()),
                       confirm_unmerged=second)
    if report.deleted:
        lines = [f"Deleted {e['name']} ({e['sha'][:10]})."
                 for e in report.deleted]
        if report.undo:
            lines.append(f"Undo: {report.undo[-1]}")
        for name, why in report.skipped:
            lines.append(f"Kept {name}: {why}.")
        return {"ok": True, "text": "\n".join(lines),
                "deleted": [e["name"] for e in report.deleted],
                "skipped": [s[0] for s in report.skipped]}
    lines = [f"{name} was not deleted: {why}."
             for name, why in report.skipped] or ["not deleted"]
    return {"ok": False, "text": "\n".join(lines)}


__all__ = ["VIEWS", "TEST_TIMEOUT_S", "MAX_TEST_FILES", "rows", "pull_main",
           "test_branch", "run_tests", "remove_test_worktree",
           "cleanup_candidates", "cleanup_delete"]
