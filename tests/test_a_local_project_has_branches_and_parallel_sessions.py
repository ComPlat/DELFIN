"""A local project (no remote) works like a hosted one, minus the remote.

* The branch views list its local branches; Pull is unavailable with the
  reason, and nothing raises.
* Several sessions of one person in one workspace run in their own
  worktrees and branches (session.create with ``worktree``), commit there,
  and are merged back on request; a conflict is reported, not written.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.agent import local_project as LP
from delfin.agent import session_worktree as SW
from delfin.dashboard import branch_cleanup as BC
from delfin.dashboard import branch_list as BL


def _git(repo, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t",
                           "-c", "user.email=t@t", *args], check=True,
                          capture_output=True, text=True).stdout


def _local_project(tmp_path) -> Path:
    ws = tmp_path / "agent_workspace"
    ws.mkdir()
    (ws / "shared.py").write_text("A = 1\nB = 2\n")
    LP.init_local_project(ws)
    return ws


# -- branch views without origin -------------------------------------------------

def test_the_branch_views_list_local_branches_without_origin(tmp_path, monkeypatch):
    monkeypatch.setenv(BL.GH_ENV, "off")
    ws = _local_project(tmp_path)
    _git(ws, "branch", "feature")
    _git(ws, "switch", "-q", "feature")
    (ws / "f.py").write_text("F = 1\n")
    _git(ws, "add", "f.py")
    _git(ws, "commit", "-qm", "feature")
    _git(ws, "switch", "-q", "main")
    _git(ws, "branch", "done")                       # merged into main
    assert BL.base_ref(ws) == "main"
    rows = BL.all_rows(ws, None)
    assert {r.ref: r.kind for r in rows} == {"main": "main", "feature": "local",
                                            "done": "local"}
    assert [r.ref for r in BL.view_rows(rows, "mine")] == ["main", "feature"]
    assert [r.ref for r in BL.view_rows(rows, "merged")] == ["done"]
    assert BL.view_rows(rows, "remote") == []
    BL.fill_ahead(ws, rows)
    assert next(r for r in rows if r.ref == "feature").ahead == 1
    cands, kept = BC.candidates(ws, gh=None)
    assert "done" in [c.name for c in cands]


def test_pull_is_unavailable_with_the_reason_without_origin(tmp_path):
    ws = _local_project(tmp_path)
    why = BL.pull_unavailable(ws)
    assert "no remote 'origin'" in why and "git remote add" in why
    plain = tmp_path / "plain"
    plain.mkdir()
    assert "not a git checkout" in BL.pull_unavailable(plain)
    _git(ws, "remote", "add", "origin", str(tmp_path / "elsewhere.git"))
    assert BL.pull_unavailable(ws) == ""


# -- sessions in their own worktrees --------------------------------------------

def test_a_session_worktree_needs_consent_in_a_plain_folder(tmp_path):
    ws = tmp_path / "agent_workspace"
    ws.mkdir()
    (ws / "a.txt").write_text("a\n")
    with pytest.raises(SW.NeedsLocalProject) as exc:
        SW.start(ws)
    assert "LOCAL git repository" in exc.value.question
    assert not (ws / ".git").exists()
    got = SW.start(ws, init_project=True)
    # The dashboard's own layout: inside the person's workspace, ignored.
    assert Path(got["path"]).parent == ws.resolve() / ".delfin" / "worktrees"
    assert got["branch"].startswith("session/")
    assert _git(ws, "status", "--porcelain").strip() == ""


def test_two_sessions_commit_in_their_worktrees_and_merge_back(tmp_path):
    ws = _local_project(tmp_path)
    one, two = SW.start(ws), SW.start(ws)
    assert one["path"] != two["path"] and one["branch"] != two["branch"]
    (Path(one["path"]) / "one.py").write_text("ONE = 1\n")
    _git(one["path"], "add", "-A")
    _git(one["path"], "commit", "-qm", "session one")
    (Path(two["path"]) / "shared.py").write_text("A = 1\nB = 3\n")
    _git(two["path"], "commit", "-qam", "session two")
    # Neither touched the main checkout or the other's files.
    assert not (ws / "one.py").exists()
    assert not (Path(two["path"]) / "one.py").exists()
    assert {w["branch"] for w in SW.listing(ws)} == {one["branch"], two["branch"]}
    first = SW.merge(one["path"])
    second = SW.merge(two["path"])
    assert (first["status"], second["status"]) == ("ok", "ok"), (first, second)
    assert (ws / "one.py").read_text() == "ONE = 1\n"
    assert (ws / "shared.py").read_text() == "A = 1\nB = 3\n"
    assert SW.listing(ws) == []                      # merged trees removed


def test_a_conflicting_session_is_reported_and_kept(tmp_path):
    ws = _local_project(tmp_path)
    one, two = SW.start(ws), SW.start(ws)
    for got, value in ((one, "10"), (two, "20")):
        (Path(got["path"]) / "shared.py").write_text(f"A = {value}\nB = 2\n")
        _git(got["path"], "commit", "-qam", value)
    assert SW.merge(one["path"])["status"] == "ok"
    res = SW.merge(two["path"])
    assert res["status"] == "conflict" and "untouched" in res["message"]
    assert (ws / "shared.py").read_text() == "A = 10\nB = 2\n"
    assert Path(two["path"]).is_dir()


def test_the_service_starts_and_merges_worktree_sessions(tmp_path):
    from delfin.agent.service import protocol
    from delfin.agent.service.client import ServiceError
    from tui_service import ServiceProcess

    with ServiceProcess(tmp_path) as svc:
        c = svc.client()
        ws = svc.workspace("agent_workspace")
        (ws / "shared.py").write_text("A = 1\n")
        params = {"cwd": str(ws), "fake": True, "fakeDelay": 0.003,
                  "worktree": True}
        with pytest.raises(ServiceError) as exc:
            c.request("session.create", params, timeout=60)
        assert exc.value.code == protocol.NEEDS_LOCAL_PROJECT
        assert "LOCAL git repository" in exc.value.message
        assert not (ws / ".git").exists()
        a = c.request("session.create", {**params, "initProject": True},
                      timeout=60)
        b = c.request("session.create", params, timeout=60)
        paths = [a["worktree"]["path"], b["worktree"]["path"]]
        assert paths[0] != paths[1]
        assert a["session"]["workspace"] == paths[0]
        listed = c.request("session.worktrees", {"cwd": str(ws)})["worktrees"]
        assert sorted(w["path"] for w in listed) == sorted(paths)
        assert all(w["sessions"] for w in listed)
        # Merging under a running session is refused: a clean merge
        # removes the tree the session works in.
        with pytest.raises(ServiceError):
            c.request("session.merge", {"path": paths[0]})
        (Path(paths[0]) / "a.py").write_text("A_FILE = 1\n")
        _git(paths[0], "add", "-A")
        _git(paths[0], "commit", "-qm", "a")
        (Path(paths[1]) / "b.py").write_text("B_FILE = 1\n")
        _git(paths[1], "add", "-A")
        _git(paths[1], "commit", "-qm", "b")
        for got in (a, b):
            c.request("session.shutdown", {"sessionId": got["sessionId"]},
                      timeout=60)
        for path in paths:
            for _ in range(200):
                try:
                    res = c.request("session.merge", {"path": path}, timeout=60)
                    break
                except ServiceError as err:
                    assert "work in" in err.message
                    import time
                    time.sleep(0.05)
            assert res["status"] == "ok", res
        assert (ws / "a.py").exists() and (ws / "b.py").exists()
