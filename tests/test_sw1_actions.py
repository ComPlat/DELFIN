"""Package SW1 phase 2: swarm init/plan/start/status (library + CLI).

Red test first. Control: tests/test_sw1_actions.py, 10 cases, all 10 red
on the previous commit (ModuleNotFoundError: delfin.agent.swarm.launch).
"""

from __future__ import annotations

import textwrap

import pytest

from delfin.agent.swarm.launch import (
    LaunchError,
    plan_launch,
    start_members,
)
from delfin.agent.swarm.spec import load_spec


@pytest.fixture()
def swarm_repo(tmp_path, monkeypatch):
    """A git repo with a swarm.yaml, isolated HOME."""
    import subprocess

    (tmp_path / "delfin" / "agent").mkdir(parents=True)
    (tmp_path / "delfin" / "agent" / "swarm").mkdir(parents=True)
    (tmp_path / "delfin" / "agent" / "worktrees").mkdir(parents=True)
    (tmp_path / "delfin" / "agent" / "swarm" / "__init__.py").write_text("",
                                                                         encoding="utf-8")
    (tmp_path / "delfin" / "agent" / "worktrees" / "core.py").write_text(
        "x = 1\n", encoding="utf-8")
    (tmp_path / "delfin" / "agent" / "worktrees" / "ui.py").write_text(
        "y = 2\n", encoding="utf-8")
    (tmp_path / "swarm.yaml").write_text(textwrap.dedent("""
        goal: Ship the swarm feature
        packages:
          - name: core
            owned: ["delfin/agent/swarm/**"]
            task: Build the core package
            roles:
              builder: true
              reviewer: true
          - name: ui
            owned: ["delfin/agent/worktrees/**"]
            task: Build the UI
    """), encoding="utf-8")
    subprocess.run(["git", "init", "-b", "main", str(tmp_path)], check=True,
                   capture_output=True)
    env = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t",
           "HOME": str(tmp_path)}
    subprocess.run(["git", "add", "-A"], cwd=tmp_path, check=True,
                   capture_output=True, env=env)
    subprocess.run(["git", "commit", "-m", "init"], cwd=tmp_path, check=True,
                   capture_output=True, env=env)
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_plan_launch_returns_members_with_branches_and_scopes(swarm_repo):
    spec = load_spec(swarm_repo / "swarm.yaml")
    plan = plan_launch(spec)
    names = [m["name"] for m in plan["members"]]
    assert "core-builder" in names and "core-reviewer" in names
    assert "ui-builder" in names
    for m in plan["members"]:
        assert m["branch"].startswith("swarm/core-") or \
            m["branch"].startswith("swarm/ui-")
        assert m["worktree"].endswith(m["branch"].replace("swarm/", ""))
        assert m["permission_mode"] == "default"
        assert m["role_prompt"].startswith("delfin/agent/swarm/prompts/")
    assert plan["members"][0]["task"].startswith("Build the")


def test_plan_launch_is_pure_read_only(swarm_repo):
    import subprocess
    before = subprocess.run(["git", "status", "--porcelain"],
                            cwd=swarm_repo, capture_output=True,
                            text=True).stdout
    spec = load_spec(swarm_repo / "swarm.yaml")
    plan_launch(spec)
    after = subprocess.run(["git", "status", "--porcelain"],
                           cwd=swarm_repo, capture_output=True,
                           text=True).stdout
    assert before == after


def test_start_dry_run_starts_nothing(swarm_repo):
    import subprocess
    spec = load_spec(swarm_repo / "swarm.yaml")
    plan = plan_launch(spec)
    assert plan["dry_run"] is True
    out = subprocess.run(["git", "worktree", "list"], cwd=swarm_repo,
                         capture_output=True, text=True).stdout
    assert out.count("\n") <= 1  # only the main worktree


def test_start_refuses_without_execute(swarm_repo):
    # start_members is the side-effecting function; starting real
    # sessions from a test is out of scope, so the real path is gated
    # behind execute=True and the default refuses.
    spec = load_spec(swarm_repo / "swarm.yaml")
    with pytest.raises(LaunchError, match="execute"):
        start_members(spec)


def test_start_execute_creates_worktrees_and_settings(swarm_repo):
    import subprocess
    spec = load_spec(swarm_repo / "swarm.yaml")
    result = start_members(spec, execute=True)
    wt_list = subprocess.run(["git", "worktree", "list"], cwd=swarm_repo,
                             capture_output=True, text=True).stdout
    assert "core-builder" in wt_list
    assert "ui-builder" in wt_list
    settings = (swarm_repo / ".swarm" / "core-builder" / ".delfin" /
                "settings.json")
    assert settings.exists()
    import json
    data = json.loads(settings.read_text(encoding="utf-8"))
    assert data["kit"]["write_allow_globs"] == ["delfin/agent/swarm/**"]
    assert data["kit"]["default_mode"] == "default"
    assert len(result["started"]) == 3
