"""A session's own worktree was created and never released, by any path.

``session_worktree`` returned a path string and dropped the WorktreeInfo
on the floor, so no caller COULD release the tree: neither
``close_session`` nor the tab's shutdown had a worktree step, and there is
no ``exit_worktree`` caller anywhere in the dashboard. Every "Own
worktree" session therefore left a checkout under
``<repo>/.delfin/worktrees/`` and an orphan ``session/<hex>`` branch
behind, whether it exited cleanly or not.

The decision about whether a tree is spare is NOT restated here. It is
``worktree.exit_worktree(keep_if_changed=True)``, which already keeps a
tree that has changes and one that background jobs are still running
inside -- the same rule the subagent teardown path uses. A second copy of
"is this tree spare" would be a second answer that can disagree with it.

Separately: the working-directory picker is a plain filesystem listing, so
a worktree belonging to a RUNNING session was offered and typeable.
Starting there makes ``repository_of()`` resolve the worktree as the root
and nests a worktree inside a worktree. The hint existed; nothing refused.

Universal: a repository is made in ``tmp_path`` with ``git init`` and one
commit, and presence is stubbed. No host paths, no assumption about the
user's home. git is shadowed nowhere -- it is a hard prerequisite of these
paths, and the suite already requires it.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from delfin.dashboard import agent_sessions as AS


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=str(repo), check=True,
                   capture_output=True)


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "t@example.invalid")
    _git(root, "config", "user.name", "t")
    (root / "a.txt").write_text("one\n", encoding="utf-8")
    _git(root, "add", "a.txt")
    _git(root, "commit", "-qm", "first")
    return root


def _worktrees(repo) -> set[str]:
    out = subprocess.run(["git", "worktree", "list", "--porcelain"],
                         cwd=str(repo), capture_output=True, text=True)
    return {line.split(" ", 1)[1] for line in out.stdout.splitlines()
            if line.startswith("worktree ")}


def _branches(repo) -> set[str]:
    out = subprocess.run(["git", "branch", "--format=%(refname:short)"],
                         cwd=str(repo), capture_output=True, text=True)
    return {b.strip() for b in out.stdout.splitlines() if b.strip()}


# ---------------------------------------------------------------------------
# The record that makes a release possible
# ---------------------------------------------------------------------------

def test_a_new_worktree_records_what_a_release_needs(repo):
    ws = AS.session_worktree(str(repo))
    side = AS.read_worktree_sidecar(ws)
    assert side is not None, "nothing recorded; no caller could release it"
    assert Path(side["repo_dir"]).resolve() == repo.resolve()
    assert side["branch"].startswith("session/")
    assert Path(side["path"]).is_dir()
    assert side["pid"] and side["host"]


def test_the_record_is_found_from_a_subdirectory(repo):
    """session_worktree returns the path at the same place INSIDE the tree
    as the original workspace was inside the repository, so the reader has
    to walk upwards."""
    ws = Path(AS.session_worktree(str(repo)))
    deep = ws / "a" / "b"
    deep.mkdir(parents=True)
    assert AS.read_worktree_sidecar(deep) is not None


def test_a_plain_directory_has_no_record(tmp_path):
    assert AS.read_worktree_sidecar(tmp_path) is None


def test_a_corrupt_record_does_not_raise(repo):
    ws = Path(AS.session_worktree(str(repo)))
    (ws / AS._WORKTREE_SIDECAR).write_text("{not json", encoding="utf-8")
    assert AS.read_worktree_sidecar(ws) is None
    out = AS.release_session_worktree(ws)
    assert out["released"] is False
    assert out["kept"] == ""


# ---------------------------------------------------------------------------
# Releasing
# ---------------------------------------------------------------------------

def test_a_clean_worktree_is_released_with_its_branch(repo):
    ws = Path(AS.session_worktree(str(repo)))
    branch = AS.read_worktree_sidecar(ws)["branch"]
    assert str(ws) in _worktrees(repo)
    assert branch in _branches(repo)

    out = AS.release_session_worktree(ws)
    assert out == {"released": True, "kept": ""}
    assert not ws.exists()
    assert str(ws) not in _worktrees(repo)
    assert branch not in _branches(repo), "an orphan session branch was left"


def test_uncommitted_work_is_never_thrown_away(repo):
    ws = Path(AS.session_worktree(str(repo)))
    (ws / "unsaved.txt").write_text("do not lose me\n", encoding="utf-8")

    out = AS.release_session_worktree(ws)
    assert out["released"] is False
    assert "uncommitted" in out["kept"]
    assert ws.is_dir()
    assert (ws / "unsaved.txt").read_text() == "do not lose me\n"


def test_a_modified_tracked_file_counts_as_work(repo):
    ws = Path(AS.session_worktree(str(repo)))
    (ws / "a.txt").write_text("changed\n", encoding="utf-8")
    out = AS.release_session_worktree(ws)
    assert out["released"] is False
    assert ws.is_dir()


def test_a_tree_a_job_is_running_in_is_kept(repo, monkeypatch):
    """The decision is exit_worktree's; this asserts that this path is
    covered by it rather than deciding for itself."""
    from delfin.agent import worktree as WT

    ws = Path(AS.session_worktree(str(repo)))
    monkeypatch.setattr(WT, "jobs_holding_worktree",
                        lambda path: [{"id": "job-1"}])
    out = AS.release_session_worktree(ws)
    assert out["released"] is False
    assert "jobs" in out["kept"]
    assert ws.is_dir()


def test_releasing_twice_is_not_an_error(repo):
    ws = Path(AS.session_worktree(str(repo)))
    assert AS.release_session_worktree(ws)["released"] is True
    out = AS.release_session_worktree(ws)
    assert out["released"] is False
    assert out["kept"] == ""


def test_nothing_to_release_is_not_a_failure(tmp_path):
    assert AS.release_session_worktree(tmp_path) == {"released": False,
                                                     "kept": ""}
    assert AS.release_session_worktree("") == {"released": False, "kept": ""}


# ---------------------------------------------------------------------------
# A directory a live session works in is refused
# ---------------------------------------------------------------------------

def _presence(monkeypatch, rows):
    from delfin.agent import session_presence as P
    monkeypatch.setattr(P, "open_sessions", lambda **kw: rows)


def test_a_directory_a_live_session_holds_is_named(tmp_path, monkeypatch):
    d = tmp_path / "shared"
    d.mkdir()
    _presence(monkeypatch, [{"key": "abcd1234", "title": "Session A",
                             "workspace": str(d), "host": "node7"}])
    label = AS._live_session_in(d)
    assert "Session A" in label
    assert "node7" in label


def test_a_directory_nobody_holds_is_free(tmp_path, monkeypatch):
    d = tmp_path / "free"
    d.mkdir()
    _presence(monkeypatch, [{"key": "abcd1234", "title": "Session A",
                             "workspace": str(tmp_path / "elsewhere")}])
    assert AS._live_session_in(d) == ""


def test_a_record_from_another_host_is_trusted_as_live(tmp_path, monkeypatch):
    """A pid from another login node names nothing here, so it cannot be
    checked and must not be assumed dead."""
    d = tmp_path / "shared"
    d.mkdir()
    _presence(monkeypatch, [{"key": "k", "title": "", "workspace": str(d),
                             "host": "some-other-node"}])
    assert AS._live_session_in(d) != ""


def test_an_unreadable_registry_does_not_block_a_start(tmp_path, monkeypatch):
    from delfin.agent import session_presence as P

    def boom(**kw):
        raise OSError("registry unreadable")

    monkeypatch.setattr(P, "open_sessions", boom)
    assert AS._live_session_in(tmp_path) == ""


def test_the_start_path_consults_it():
    import inspect

    src = inspect.getsource(AS)
    start = src.index("def _on_start")
    body = src[start:src.index("def _on_resume", start)]
    assert "_live_session_in" in body
    assert "Own" in body, "the refusal must name the way out"
    # Ticking "own worktree" IS the answer to a shared directory, so the
    # refusal must not fire then.
    assert "_fresh" in body


def test_closing_a_session_releases_its_worktree():
    import inspect

    src = inspect.getsource(AS)
    start = src.index("def close_session")
    body = src[start:start + 2000]
    assert "release_session_worktree" in body
    i_shutdown = body.index('refs"].get("shutdown")')
    assert i_shutdown < body.index("release_session_worktree"), (
        "a job holding the tree must be asked to stop before the decision")
