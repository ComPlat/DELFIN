"""pr_helper.body: what reaches a public pull request, and from where.

The body is published to GitHub. These cases pin the inputs it may draw
on (this session's tasks, this branch's own commits) and the bounds it
keeps (the process identity redacted by default, no commit trailers, a
length GitHub accepts).
"""

from __future__ import annotations

import json
import subprocess

import pytest

from conftest import child_env
from delfin.agent import pr_helper

_FAKE_IDENTITY = {"home": "/home/fixtureacct", "account": "fixtureacct",
                  "host": "fixturehost"}


def _store(ws, tasks):
    (ws / ".delfin").mkdir(parents=True, exist_ok=True)
    (ws / ".delfin" / "session_tasks.json").write_text(
        json.dumps({"next_id": len(tasks) + 1, "tasks": tasks}),
        encoding="utf-8")


def _one_commit(body=""):
    return lambda b, h: [{"hash": "abc1234", "subject": "A change",
                          "body": body}]


# --------------------------------------------------------------------------
# The journal: one workspace store holds every session's tasks
# --------------------------------------------------------------------------

def _two_sessions(ws):
    _store(ws, [
        {"id": 1, "subject": "mine in progress", "status": "in_progress",
         "session_id": "s-mine"},
        {"id": 2, "subject": "other in progress", "status": "in_progress",
         "session_id": "s-other", "updated_at": "2099-01-01T00:00:00Z"},
        {"id": 3, "subject": "other blocked", "status": "blocked",
         "blocked_reason": "another session's plan", "session_id": "s-other"},
    ])


def test_another_sessions_tasks_do_not_reach_the_body(tmp_path):
    _two_sessions(tmp_path)
    out = pr_helper.body(tmp_path, "main", git_log=_one_commit(),
                         identity=_FAKE_IDENTITY, session_id="s-mine")
    assert "mine in progress" in out
    assert "other in progress" not in out
    assert "another session's plan" not in out


def test_without_a_session_id_owned_tasks_are_not_read(tmp_path):
    _two_sessions(tmp_path)
    out = pr_helper.body(tmp_path, "main", git_log=_one_commit(),
                         identity=_FAKE_IDENTITY)
    assert "in progress" not in out
    assert "another session's plan" not in out


# --------------------------------------------------------------------------
# Identity: the default protects the account the process runs as
# --------------------------------------------------------------------------

def test_the_default_identity_is_the_process_identity(tmp_path, monkeypatch):
    from delfin.agent import api_client
    monkeypatch.setattr(api_client, "_identity_for_commit_texts",
                        lambda: dict(_FAKE_IDENTITY))
    out = pr_helper.body(
        tmp_path, "main",
        git_log=_one_commit("Measured on fixturehost as fixtureacct."))
    assert "fixturehost" not in out
    assert "fixtureacct" not in out


# --------------------------------------------------------------------------
# Trailers and size
# --------------------------------------------------------------------------

def test_commit_trailers_are_not_copied(tmp_path):
    out = pr_helper.body(
        tmp_path, "main", identity=_FAKE_IDENTITY,
        git_log=_one_commit("Why.\n\nCo-Authored-By: Some Model "
                            "<noreply@example.com>\nSigned-off-by: A <a@b>"))
    assert "Why." in out
    assert "co-authored-by" not in out.lower()
    assert "signed-off-by" not in out.lower()


def test_the_body_fits_the_github_limit(tmp_path):
    many = [{"hash": f"{i:07x}", "subject": f"change {i}",
             "body": "x" * 3000} for i in range(400)]
    out = pr_helper.body(tmp_path, "main", identity=_FAKE_IDENTITY,
                         git_log=lambda b, h: many)
    assert len(out) <= pr_helper.BODY_LIMIT <= 65536
    assert "change 0" in out


# --------------------------------------------------------------------------
# The real git reader
# --------------------------------------------------------------------------

def _git(repo, *args, env):
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@e.invalid",
                    "-c", "commit.gpgsign=false", *args], cwd=repo, env=env,
                   check=True, capture_output=True, stdin=subprocess.DEVNULL)


@pytest.fixture
def repo(tmp_path):
    env = child_env(tmp_path)
    r = tmp_path / "repo"
    r.mkdir()
    _git(r, "init", "-q", "-b", "main", env=env)
    _git(r, "commit", "-q", "--allow-empty", "-m", "base", env=env)
    _git(r, "checkout", "-q", "-b", "side", env=env)
    _git(r, "commit", "-q", "--allow-empty", "-m", "side work", env=env)
    _git(r, "checkout", "-q", "-b", "feature", "main", env=env)
    _git(r, "commit", "-q", "--allow-empty", "-m", "feature work", env=env)
    _git(r, "merge", "-q", "--no-ff", "--no-edit", "side", env=env)
    return r


def test_merge_commits_are_not_listed(repo):
    commits = pr_helper._default_git_log("main", "HEAD", cwd=str(repo))
    subjects = [c["subject"] for c in commits]
    assert "feature work" in subjects and "side work" in subjects
    assert not any(s.startswith("Merge") for s in subjects)


def test_an_unknown_base_is_an_error_not_an_empty_branch(repo):
    with pytest.raises(RuntimeError):
        pr_helper._default_git_log("no-such-base", "HEAD", cwd=str(repo))
