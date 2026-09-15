"""Sessions in one repository stay out of each other's way.

On 2026-09-15 the DELFIN agent and Claude Code edited the same checkout at
the same time. Open sessions now say where they work; a new session in a
repository another one is using gets a worktree of its own, and an agent
is told who else works in its repository.
"""

from __future__ import annotations

import json
import os
import subprocess
import time

import pytest

from delfin.agent import session_presence as P


def _git(cwd, *args):
    subprocess.run(["git", "-C", str(cwd), *args], check=True,
                   capture_output=True, text=True)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._git_cache.clear()
    P._last_written.clear()
    root = tmp_path / "project"
    root.mkdir()
    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "t@example.invalid")
    _git(root, "config", "user.name", "t")
    (root / "a.py").write_text("x = 1\n", encoding="utf-8")
    _git(root, "add", "a.py")
    _git(root, "commit", "-q", "-m", "init")
    return root


def test_a_worktree_is_the_same_repository(repo):
    from delfin.agent import worktree as W
    info = W.enter_worktree(repo, branch_prefix="session",
                            parent=repo / ".delfin" / "worktrees")
    here, there = P.repository_of(str(repo)), P.repository_of(str(info.path))
    assert here["common_dir"] == there["common_dir"]
    assert here["root"] != there["root"]
    assert there["branch"].startswith("session/")


def test_the_sessions_in_a_repository_are_found(repo, tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    P.announce("A", title="refactor", workspace=str(repo))
    P.announce("B", title="docs", workspace=str(repo))
    P.announce("C", title="unrelated", workspace=str(elsewhere))
    found = P.in_same_repository(str(repo), exclude_key="A")
    assert [r["title"] for r in found] == ["docs"]


def test_a_session_that_is_gone_is_not_listed(repo):
    P.announce("A", title="old", workspace=str(repo))
    record = json.loads(P._path("A").read_text(encoding="utf-8"))
    P._path("A").write_text(json.dumps({**record,
                                        "updated_at": time.time() - 3600}),
                            encoding="utf-8")
    assert P.open_sessions() == []

    P._last_written.clear()
    P.announce("B", title="crashed", workspace=str(repo))
    record = json.loads(P._path("B").read_text(encoding="utf-8"))
    dead = subprocess.Popen(["true"])
    dead.wait()
    P._path("B").write_text(json.dumps({**record, "pid": dead.pid}),
                            encoding="utf-8")
    assert P.open_sessions() == []

    P.withdraw("B")
    assert not P._path("B").exists()


def test_the_agent_is_told_who_else_works_in_its_repository(repo, monkeypatch):
    from delfin.agent.engine import AgentEngine

    P.announce("mine", title="this session", workspace=str(repo))
    P.announce("theirs", title="fix the collector", workspace=str(repo))

    class _Perms:
        workspace = str(repo)
        presence_key = "mine"

    monkeypatch.setattr(AgentEngine, "kit_permissions",
                        property(lambda self: _Perms()))
    block = AgentEngine.__new__(AgentEngine)._build_other_sessions_block()
    assert "fix the collector" in block and "this session" not in block
    assert "git stash" in block

    P.withdraw("theirs")
    assert AgentEngine.__new__(AgentEngine)._build_other_sessions_block() == ""


def test_a_new_session_in_a_repository_in_use_gets_its_own_worktree(
        repo, tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.agent import scheduler as S
    from delfin.dashboard import agent_sessions as AS
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(repo))

    def build(ctx):
        state = {"active_session_id": "", "chat_messages": [],
                 "streaming": False, "engine": None}
        return widgets.VBox(), {"state": state, "shutdown": lambda: None}

    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "agent_workspace")
    _widget, refs = AS.create_tab(ctx, build=build)
    form = refs["form"]

    form["new"].click()
    form["workdir"].value = str(repo)
    assert form["own_worktree"].value is True, "another session works here"
    form["start"].click()

    first, second = refs["sessions"]()
    assert first["workspace"] == str(repo)
    assert second["workspace"].startswith(str(repo / ".delfin" / "worktrees"))
    assert P.repository_of(second["workspace"])["branch"].startswith("session/")
    status = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"],
                            capture_output=True, text=True).stdout
    assert ".delfin" not in status, "the worktrees stay out of git status"
    for rec in refs["sessions"]():
        refs["close"](rec["key"])
    assert os.listdir(P._DIR) == [] or all(
        not f.startswith(first["key"]) for f in os.listdir(P._DIR))
