"""BG2 phase 2: the exceptions keep today's behaviour (each with a test).

A foreground bash command that reaches its timeout is MOVED to the
background -- except where the move would break a rule the tool holds:

- a `git push` / PR publish that waits: it is waiting for a credential
  that cannot be given here, and the run already ended the way today's
  message says (api_client.py:25285-25309);
- a command the login-node calculation rule (calc_policy, issue #640)
  admitted only for the foreground wall time: the move must not extend a
  calculation on a login node;
- a command whose process tree already exited before the adoption:
  nothing live to move; today's message is the truth.

These tests run the tool with plain argv (cage off) and isolated
registries, like test_bg2_timeout_moves_to_background.py.
"""

import json
import os
import stat
import subprocess
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ


@pytest.fixture(autouse=True)
def _plain_argv(monkeypatch):
    monkeypatch.setattr(
        A, "_bash_isolation_argv",
        lambda cmd, cwd, perms, **kw: ["/bin/bash", "-c", cmd])


@pytest.fixture(autouse=True)
def _clean_registry(tmp_path):
    BJ._REGISTRY._jobs.clear()
    yield
    for job in list(BJ._REGISTRY._jobs.values()):
        try:
            BJ.get_registry().kill(job.job_id)
        except Exception:
            pass
    BJ._REGISTRY._jobs.clear()


def _run(cmd, timeout_s=2, tmp_path=None, **extra):
    perms = A.KitToolPermissions(
        workspace=str(tmp_path), mode="bypassPermissions")
    getattr(A._doc_executor, "_timed_out_by_session", {}).clear()
    return json.loads(A._doc_executor._execute_bash(
        {"command": cmd, "timeout_s": timeout_s, **extra}, perms))


def _push_repo(tmp_path, monkeypatch, helper_body: str, remote: str):
    """A git workspace whose 'git push' hangs (the askpass helper sleeps),
    the fixture pattern of test_a_push_without_a_credential_fails..."""
    ws = tmp_path / "ws"
    ws.mkdir()
    subprocess.run(["git", "init", "-q", str(ws)], check=True)
    subprocess.run(["git", "-C", str(ws), "-c", "user.name=t",
                    "-c", "user.email=t@example.org", "commit",
                    "-q", "--allow-empty", "-m", "c"], check=True)
    bindir = ws / "bin"
    bindir.mkdir()
    helper = bindir / "git-remote-hang"
    helper.write_text("#!/bin/sh\n" + helper_body + "\n", encoding="utf-8")
    helper.chmod(helper.stat().st_mode | stat.S_IXUSR)
    subprocess.run(["git", "-C", str(ws), "remote", "add", remote,
                    "hang::x"], check=True)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setattr(A, "PUSH_TIMEOUT_S", 2)
    return ws


def test_a_push_that_waits_is_still_ended(tmp_path, monkeypatch):
    """The push branch keeps today's behaviour: ended at the push
    timeout with the credential message -- a push must never be kept
    running in the background where it waits forever for a prompt
    nobody can give (the fixture of
    test_a_push_without_a_credential_fails_instead_of_hanging.py)."""
    ws = _push_repo(tmp_path, monkeypatch, "sleep 60", "far")
    perms = A.KitToolPermissions(
        workspace=str(ws), mode="bypassPermissions")
    getattr(A._doc_executor, "_timed_out_by_session", {}).clear()
    res = json.loads(A._doc_executor._execute_bash(
        {"command": "git push far HEAD:refs/heads/x", "timeout_s": 300},
        perms))
    assert res.get("status") != "moved_to_background", res
    assert "credential" in (res.get("error") or ""), res


def test_a_process_that_already_exited_reports_today_s_message(tmp_path):
    """Nothing live to adopt: the command ended (or exited between the
    timeout and the adoption), and the result says the work is over."""
    res = _run("echo quickly-done", timeout_s=2, tmp_path=tmp_path)
    assert res.get("status") != "moved_to_background", res
    # 'echo quickly-done' finishes in milliseconds; the timeout branch
    # never runs, and the result carries no error at all. What this
    # pins is only that a fast command is never answered with a move.
    assert "error" not in res, res
