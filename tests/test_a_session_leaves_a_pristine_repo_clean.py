"""A session in a pristine repository leaves ``git status`` empty.

Input: a one-commit git repository and an agent engine rooted in it.
Output: ``git status --porcelain`` after a prompt build, and where the
background-job registry lives.

Measured on main with the foreign-repo probe: every one of four pristine
repositories showed ``?? .delfin/`` after the first prompt build -- the
completion-event drain takes a cross-process lock beside the registry,
and the registry was ``<workspace>/.delfin/bash_jobs.json``. The first
thing the agent then read about the repository was an untracked change
it had made itself, and a careless ``git add -A`` commits it.

The registry also lived in the one folder the agent may write and a
cloned repository may ship: a record there names the files whose tail is
put into the model's context.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from delfin.agent import bash_jobs as BJ

from conftest import child_env


def _git(ws: Path, tmp_path: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ws, env=child_env(tmp_path), check=True,
        capture_output=True, text=True).stdout


@pytest.fixture
def repo(tmp_path) -> Path:
    ws = tmp_path / "proj"
    ws.mkdir()
    (ws / "README.md").write_text("x\n", encoding="utf-8")
    _git(ws, tmp_path, "init", "-q", "-b", "main")
    _git(ws, tmp_path, "add", "README.md")
    _git(ws, tmp_path, "-c", "user.name=t", "-c", "user.email=t@example.org",
         "commit", "-q", "-m", "init")
    assert _git(ws, tmp_path, "status", "--porcelain") == ""
    return ws


def test_draining_and_listing_jobs_writes_nothing_into_the_workspace(repo):
    BJ.drain_finished_events(repo)
    BJ.drain_all_finished_events(repo)
    BJ.running_jobs_for_workspace(repo)
    assert not (repo / ".delfin").exists()


def test_the_registry_lives_outside_the_workspace(repo):
    path = BJ._registry_path(repo)
    with pytest.raises(ValueError):
        path.resolve().relative_to(repo.resolve())


def test_a_started_job_is_recorded_outside_the_workspace(repo):
    BJ._persist_job_start(str(repo), {"job_id": "kit_bg_x", "pid": 1,
                                      "started_at": 1.0})
    assert not (repo / ".delfin").exists()
    assert "kit_bg_x" in BJ._load_registry_file(repo)["jobs"]


def test_two_workspaces_keep_two_registries(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    assert BJ._registry_path(a) != BJ._registry_path(b)


def test_a_registry_shipped_in_the_repository_is_not_read(repo):
    shipped = repo / ".delfin" / "bash_jobs.json"
    shipped.parent.mkdir()
    shipped.write_text(json.dumps({"jobs": {"kit_bg_planted": {
        "pid": 1, "stdout_path": "/etc/hostname"}}}), encoding="utf-8")
    assert BJ._load_registry_file(repo)["jobs"] == {}


def test_a_prompt_build_leaves_git_status_empty(repo, tmp_path):
    """The real client and permission layer, as the foreign-repo probe
    builds them; nothing is sent (no request is made)."""
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine(repo_dir=repo, backend="api", provider="kit",
                      api_key="stub-key", model="kit.glm-5.3", mode="solo",
                      permission_mode="default")
    eng.messages = [{"role": "user", "content": "fix the failing test"}]
    eng._build_current_system_prompt(task_text="fix the failing test")
    assert _git(repo, tmp_path, "status", "--porcelain",
                "--untracked-files=all") == ""
