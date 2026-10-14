"""A push over HTTPS uses the user's gh login, outside the cage only.

Field report 2026-10-14: `git push` stayed blocked for a user although
`gh auth login` had succeeded. The readiness probe runs outside the
cage and passed; the push itself ran inside it, where ~/.config/gh is
masked (sandbox._HOME_SECRET_DIRS), so git's credential helper `gh auth
git-credential` found no login and the push failed for want of a
credential the user had just given.

  https://github.com/<o>/<r> remote, gh on PATH   DELFIN's push path:
                                                  outside the cage, hooks
                                                  off, only that scheme,
                                                  every other credential
                                                  helper reset, gh's added
  no gh, or the checkout's own config names http.*,
  https.*, credential.* or url.*                  stays in the cage
  the token                                       reaches git's helper
                                                  outside; nothing in the
                                                  cage sees it
  each refusal / failure                          says its own cause
"""

from __future__ import annotations

import json
import os
import re
import stat
import subprocess

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

TOKEN = "gho_standin0000000000000000000000000000"


def _env():
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


@pytest.fixture
def repo(tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()

    def git(*args, cwd=ws):
        return subprocess.run(["git", "-C", str(cwd), "-c", "user.name=t",
                               "-c", "user.email=t@example.org", *args],
                              check=True, capture_output=True, text=True,
                              env=_env())

    git("init", "-q")
    git("commit", "-q", "--allow-empty", "-m", "c")
    git("remote", "add", "origin", "https://github.com/ComPlat/DELFIN.git")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    gh = bindir / "gh"
    gh.write_text(
        "#!/bin/sh\n# stand-in for `gh auth git-credential get`\n"
        "cat >/dev/null\n"
        f"printf 'username=x-access-token\\npassword={TOKEN}\\n'\n",
        encoding="utf-8")
    gh.chmod(gh.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    return ws, git, gh


def test_an_https_github_push_takes_the_gh_route(repo):
    ws, git, gh = repo
    plan = A._github_push_plan("git push -u origin topic", ws)
    assert plan is not None
    argv = plan["argv"]
    assert plan.get("via") == "gh"
    for setting in ("core.hooksPath=/dev/null", "protocol.allow=never",
                    "protocol.https.allow=always", "credential.helper=",
                    "push.recurseSubmodules=no"):
        assert setting in argv, setting
    helper = [a for a in argv if a.startswith("credential.helper=!")]
    assert len(helper) == 1 and str(gh) in helper[0]
    # The reset comes first, so no other helper is asked.
    assert argv.index("credential.helper=") < argv.index(helper[0])


def test_without_gh_it_stays_in_the_cage(repo, monkeypatch):
    ws, git, gh = repo
    gh.unlink()
    monkeypatch.setattr(A.shutil, "which",
                        lambda name, *a, **k: None if name == "gh" else
                        __import__("shutil").which(name, *a, **k))
    assert A._github_push_plan("git push origin topic", ws) is None


@pytest.mark.parametrize("key,value", [
    ("http.proxy", "http://127.0.0.1:9"),
    ("http.sslVerify", "false"),
    ("credential.helper", "store"),
    ("url.https://evil.example/.insteadOf", "https://github.com/"),
])
def test_a_checkout_that_configures_the_transport_stays_in_the_cage(
        repo, key, value):
    """The checkout's .git/config is writable by an agent command; a proxy
    or helper set there must not receive the user's token."""
    ws, git, gh = repo
    git("config", key, value)
    assert A._github_push_plan("git push origin topic", ws) is None


def test_the_helper_hands_the_token_to_git_outside_the_cage(repo):
    ws, git, gh = repo
    plan = A._github_push_plan("git push origin topic", ws)
    opts = plan["argv"][1:plan["argv"].index("push")]
    done = subprocess.run(
        ["git", *opts, "credential", "fill"], cwd=str(ws),
        input="protocol=https\nhost=github.com\n\n", capture_output=True,
        text=True, env={**_env(), "GIT_TERMINAL_PROMPT": "0"}, timeout=20)
    assert f"password={TOKEN}" in done.stdout


def test_no_caged_command_sees_the_token(repo):
    ws, git, gh = repo
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    out = json.loads(eng._execute_bash({
        "command": "printf 'protocol=https\\nhost=github.com\\n\\n' | "
                   "git credential fill 2>&1; env"}, perms))
    assert TOKEN not in out.get("stdout", "") + out.get("stderr", "")


def test_the_route_pushes_to_the_remote(repo, tmp_path, monkeypatch):
    """End to end through the executor, against a local bare remote: the
    URL check is widened to that remote's file:// URL for the test, and
    the plan admits exactly that scheme."""
    ws, git, gh = repo
    bare = tmp_path / "bare.git"
    subprocess.run(["git", "init", "-q", "--bare", str(bare)], check=True,
                   env=_env())
    git("remote", "set-url", "origin", f"file://{bare}")
    monkeypatch.setattr(A, "_GITHUB_HTTPS_URL_RE",
                        re.compile("^" + re.escape(f"file://{bare}") + "$"))
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    out = json.loads(eng._execute_bash(
        {"command": "git push origin HEAD:refs/heads/topic"}, perms))
    assert out["exit_code"] == 0, out
    refs = subprocess.run(["git", "-C", str(bare), "show-ref"],
                          capture_output=True, text=True, env=_env()).stdout
    assert "refs/heads/topic" in refs


def test_a_failed_gh_route_names_the_missing_login(repo, monkeypatch):
    ws, git, gh = repo
    monkeypatch.setattr(A, "_run_github_push", lambda *a, **k:
                        subprocess.CompletedProcess(
                            ["git"], 128, "",
                            "fatal: could not read Username for "
                            "'https://github.com': terminal prompts disabled\n"))
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    out = json.loads(eng._execute_bash(
        {"command": "git push origin topic"}, perms))
    hint = out.get("hint", "")
    assert hint.startswith("push failed: no credentials reachable")
    assert "gh auth login" in hint


# -- the gate's own refusals say which one they are ---------------------------

@pytest.fixture
def gate(repo, monkeypatch):
    ws, git, gh = repo
    monkeypatch.setattr(A, "_push_readiness_rows", lambda *a, **k: [])

    def _run(cmd, *, role="maintainer", grant=True):
        monkeypatch.setattr(A, "_git_role", lambda: role)
        perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
        perms.push_grants = {"push": 1} if grant else {}
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        return eng._run_permission_gate("bash", {"command": cmd}, perms)
    return _run


def test_a_contributor_push_to_the_default_branch_says_so(gate, repo):
    ws, git, gh = repo
    branch = git("rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
    git("update-ref", f"refs/remotes/origin/{branch}", "HEAD")
    git("symbolic-ref", "refs/remotes/origin/HEAD",
        f"refs/remotes/origin/{branch}")
    out = gate(f"git push origin {branch}", role="contributor")
    assert out and "contributor" in out and "/git role maintainer" in out


def test_an_unrequested_push_says_so(gate):
    out = gate("git push origin topic", grant=False)
    assert out and "has not asked for a push" in out


def test_a_requested_branch_push_passes_the_gate(gate):
    assert gate("git push origin topic") is None
