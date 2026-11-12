"""The agent's git commit gets the co-author trailer mechanically.

The prompt asked the model to end its commits with
``Co-Authored-By: DELFIN-Agent (<model>)``, and the model often omitted it
(22 of 77 commits on one colleague's PR). The executor now appends the
trailer itself, through git's own trailer mechanism, to every `git commit`
a model turn runs. Never twice (git dedups identical --trailer values),
never with ``agent.commit_coauthor = false``, and never for a command the
user runs from the agent's prompt (`!cmd`), which is not the model's
commit.

The tests drive the REAL executor (_doc_executor.execute), not a stub of
the rewrite, on a throwaway repository.
"""
import os
import subprocess
from pathlib import Path

import pytest

from delfin.agent import api_client as A

TRAILER = ("Co-Authored-By: DELFIN-Agent (kit.glm-5.3-flash) "
           "<330454853+delfin-agent@users.noreply.github.com>")


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A git repository and executor permissions pointing at it."""
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
    monkeypatch.setenv("GIT_CONFIG_SYSTEM", "/dev/null")
    repo = tmp_path / "repo"
    repo.mkdir()
    env = {"GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_SYSTEM": "/dev/null",
           "PATH": os.environ["PATH"]}
    run = lambda *a: subprocess.run(["git", *a], cwd=repo, env=env,
                                    check=True, capture_output=True)
    run("init", "-q", "-b", "main")
    run("-c", "user.email=t@t", "-c", "user.name=t",
        "commit", "-q", "--allow-empty", "-m", "init")
    return repo


def _agent_commit(repo, monkeypatch, command, model="kit.glm-5.3-flash",
                  settings=None):
    """One model bash call through the real executor; last commit message."""
    import delfin.user_settings as us
    # bash_isolation "off" ({...off}, escaped as such) lets the executor run
    # the git path without the host having to install a seccomp filter,
    # which this sandboxed test host cannot (nested seccomp: EBUSY). The
    # trailer rewrite is exercised through the REAL executor either way.
    base = {"agent": {"bash_isolation": "off"}}
    if settings:
        base["agent"].update(settings["agent"])
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: base)
    A.set_commit_trailer_model(model)
    try:
        perms = A.KitToolPermissions(workspace=str(repo), mode="bypassPermissions")
        result = A._doc_executor.execute("bash", {"command": command}, perms)
        assert "error" not in result, result
    finally:
        A.set_commit_trailer_model("")
    msg = subprocess.run(
        ["git", "log", "-1", "--format=%B"], cwd=repo,
        capture_output=True, text=True).stdout
    return result, msg


def test_an_agent_commit_carries_exactly_one_trailer(repo, monkeypatch):
    _, msg = _agent_commit(
        repo, monkeypatch,
        "echo hi > f.txt && git add f.txt && "
        "git -c user.email=t@t -c user.name=t commit -q -m work")
    assert msg.count("Co-Authored-By:") == 1, msg
    assert TRAILER in msg, msg


def test_an_amend_does_not_duplicate(repo, monkeypatch):
    _agent_commit(
        repo, monkeypatch,
        "echo hi > f.txt && git add f.txt && "
        "git -c user.email=t@t -c user.name=t commit -q -m work")
    _, msg = _agent_commit(
        repo, monkeypatch,
        "git -c user.email=t@t -c user.name=t commit -q --amend "
        "--allow-empty -m work2")
    assert msg.count("Co-Authored-By:") == 1, msg
    assert TRAILER in msg, msg


def test_a_commit_that_already_carries_it_stays_single(repo, monkeypatch):
    _, msg = _agent_commit(
        repo, monkeypatch,
        "echo hi > f.txt && git add f.txt && "
        "git -c user.email=t@t -c user.name=t commit -q -m work "
        f"--trailer '{TRAILER}'")
    assert msg.count("Co-Authored-By:") == 1, msg


def test_the_setting_turns_it_off(repo, monkeypatch):
    _, msg = _agent_commit(
        repo, monkeypatch,
        "echo hi > f.txt && git add f.txt && "
        "git -c user.email=t@t -c user.name=t commit -q -m work",
        settings={"agent": {"commit_coauthor": False}})
    assert "Co-Authored-By:" not in msg, msg


def test_a_user_command_from_the_prompt_is_untouched(repo, monkeypatch):
    """`!cmd` runs inside output_for_a_person: the commit is the USER's."""
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {})
    A.set_commit_trailer_model("kit.glm-5.3-flash")
    try:
        perms = A.KitToolPermissions(workspace=str(repo), mode="bypassPermissions")
        with A.output_for_a_person():
            res = A._doc_executor._execute_bash(
                {"command":
                 "echo hi > f.txt && git add f.txt && "
                 "git -c user.email=t@t -c user.name=t commit -q -m mine"},
                perms)
        assert "error" not in res, res
    finally:
        A.set_commit_trailer_model("")
    msg = subprocess.run(["git", "log", "-1", "--format=%B"], cwd=repo,
                         capture_output=True, text=True).stdout
    assert "Co-Authored-By:" not in msg, msg


def test_the_email_comes_from_the_setting(repo, monkeypatch):
    _, msg = _agent_commit(
        repo, monkeypatch,
        "echo hi > f.txt && git add f.txt && "
        "git -c user.email=t@t -c user.name=t commit -q -m work",
        settings={"agent": {
            "commit_coauthor_email":
                "12345678+delfin-agent@users.noreply.github.com"}})
    assert ("Co-Authored-By: DELFIN-Agent (kit.glm-5.3-flash) "
            "<12345678+delfin-agent@users.noreply.github.com>") in msg, msg
