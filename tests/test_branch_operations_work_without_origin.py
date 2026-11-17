"""branch_ops on a local project (a repository with no origin).

The surfaces' branch logic (/branches in the terminal, the messenger's
branch words, the service's branch.* methods) answers on such a
repository: the rows are its local branches, Pull is unavailable with
the reason and runs no git, a cleanup listing works, nothing raises.
"""

from __future__ import annotations

import subprocess

from delfin.agent import branch_ops as BO
from delfin.agent import local_project as LP
from delfin.dashboard import branch_list as BL


def _git(repo, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t",
                           "-c", "user.email=t@t", *args], check=True,
                          capture_output=True, text=True).stdout


def _project(tmp_path):
    ws = tmp_path / "agent_workspace"
    ws.mkdir()
    (ws / "a.py").write_text("A = 1\n")
    LP.init_local_project(ws)
    _git(ws, "branch", "feature")
    _git(ws, "switch", "-q", "feature")
    (ws / "b.py").write_text("B = 1\n")
    _git(ws, "add", "b.py")
    _git(ws, "commit", "-qm", "feature")
    _git(ws, "switch", "-q", "main")
    return ws


def test_rows_list_the_local_branches(tmp_path, monkeypatch):
    monkeypatch.setenv(BL.GH_ENV, "off")
    ws = _project(tmp_path)
    got = BO.rows(ws, "all", gh=None)
    assert got["ok"] and [r["ref"] for r in got["rows"]] == ["main", "feature"]
    assert BO.rows(ws, "remote", gh=None)["rows"] == []


def test_pull_is_unavailable_with_the_reason_and_runs_no_git(tmp_path):
    ws = _project(tmp_path)
    calls = []
    out = BO.pull_main(ws, read=lambda *a: calls.append(a),
                       write=lambda *a: calls.append(a))
    assert out["ok"] is False and out.get("unavailable") is True
    assert "no remote 'origin'" in out["text"]
    assert calls == []


def test_cleanup_answers_without_origin(tmp_path, monkeypatch):
    monkeypatch.setenv(BL.GH_ENV, "off")
    ws = _project(tmp_path)
    _git(ws, "branch", "done")
    cands, kept = BO.cleanup_candidates(ws)
    assert "done" in [c["name"] for c in cands]
