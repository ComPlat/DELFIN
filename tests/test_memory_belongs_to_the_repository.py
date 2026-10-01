"""The per-project memory key is the repository, not the working directory.

Input: a workspace path. Output: the slug naming its memory store.

Measured on 2026-09-24 in ~/.delfin/projects:

    throwaway keys   538 directories, 358 notes
    real projects      4 directories,  13 notes

96% of everything the agent had written sat under keys derived from paths
that no longer exist: session worktrees under ``.delfin/worktrees/``,
``/tmp`` scratch trees, probe directories. A worktree of a repository is
the same repository; keying on the resolved working directory gave each
one a private store that was orphaned when the directory went.

Consequence for machinery that already exists and could not act:
near-duplicate merging cannot merge across stores, BM25 selection cannot
rank what it does not load, and the 90-day disuse decay never fires on a
store nothing recalls from.

Semantics:
  - inside a git repository: the main worktree's root, so every linked
    worktree shares one store
  - outside a repository: the resolved path, unchanged
  - unreadable git, missing binary, timeout: the resolved path, unchanged

The fallback direction is deliberate. A wrong answer that SPLITS a store
costs recall; a wrong answer that MERGES two projects puts one project's
notes in another's prompt.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.agent.memory_store import _project_slug


def _repo(path):
    path.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "init", "-q"], cwd=str(path), check=True)
    for k, v in (("user.email", "a@b.c"), ("user.name", "t")):
        subprocess.run(["git", "config", k, v], cwd=str(path), check=True)
    (path / "f.txt").write_text("x\n", encoding="utf-8")
    subprocess.run(["git", "add", "-A"], cwd=str(path), check=True)
    subprocess.run(["git", "commit", "-q", "-m", "x"], cwd=str(path),
                   check=True)
    return path


def test_a_linked_worktree_shares_the_repository_store(tmp_path):
    repo = _repo(tmp_path / "repo")
    wt = tmp_path / "wt"
    subprocess.run(["git", "worktree", "add", "-q", str(wt), "-b", "side"],
                   cwd=str(repo), check=True)
    assert _project_slug(wt) == _project_slug(repo)


def test_a_worktree_inside_the_repository_shares_it_too(tmp_path):
    """The shape DELFIN itself creates: .delfin/worktrees/<name>."""
    repo = _repo(tmp_path / "repo")
    wt = repo / ".delfin" / "worktrees" / "delfin-wt-475b4d87"
    wt.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "worktree", "add", "-q", str(wt), "-b", "s2"],
                   cwd=str(repo), check=True)
    assert _project_slug(wt) == _project_slug(repo)


def test_a_subdirectory_shares_it(tmp_path):
    repo = _repo(tmp_path / "repo")
    sub = repo / "delfin" / "agent"
    sub.mkdir(parents=True)
    assert _project_slug(sub) == _project_slug(repo)


def test_two_different_repositories_stay_apart(tmp_path):
    """The failure that costs more: one project's notes in another's
    prompt."""
    a = _repo(tmp_path / "a")
    b = _repo(tmp_path / "b")
    assert _project_slug(a) != _project_slug(b)


def test_a_directory_outside_any_repository_keeps_its_path(tmp_path):
    plain = tmp_path / "plain"
    plain.mkdir()
    assert _project_slug(plain) == "-" + str(plain.resolve()).replace(
        "/", "-").lstrip("-")


def test_the_slug_is_still_a_single_path_segment(tmp_path):
    repo = _repo(tmp_path / "repo")
    slug = _project_slug(repo)
    assert "/" not in slug
    assert slug.startswith("-")


def test_it_falls_back_instead_of_raising(tmp_path, monkeypatch):
    """No git binary, a timeout, a corrupt repository: answer the path.
    Splitting a store costs recall; merging two projects costs more."""
    import delfin.agent.memory_store as M

    def _boom(*a, **k):
        raise OSError("no git here")

    monkeypatch.setattr(M.subprocess, "run", _boom)
    plain = tmp_path / "anything"
    plain.mkdir()
    assert _project_slug(plain) == "-" + str(plain.resolve()).replace(
        "/", "-").lstrip("-")


def test_the_answer_is_cached_per_path(tmp_path):
    """It is asked on every memory read and write; one subprocess per
    path, not per call."""
    repo = _repo(tmp_path / "repo")
    calls = {"n": 0}
    import delfin.agent.memory_store as M

    real = M.subprocess.run

    def _count(*a, **k):
        calls["n"] += 1
        return real(*a, **k)

    M._repo_root_cache.clear()
    M.subprocess.run = _count
    try:
        for _ in range(5):
            _project_slug(repo)
    finally:
        M.subprocess.run = real
    assert calls["n"] <= 1, f"{calls['n']} git calls for five lookups"
