"""The agent does not switch the branch of the shared checkout.

Input: a command and the directory it runs in; the shared checkout is the
one the dashboard names in ``DELFIN_SHARED_CHECKOUT`` (here a repository
made in ``tmp_path``) or the one DELFIN's code is imported from. Output: a
refusal for every branch move there (``git switch``, ``git checkout
<branch>``, ``-b``, ``--detach``, ``-``, through ``-C`` or after ``cd``)
that names the routes that work (enter_worktree, the header's ▶ Test);
nothing for file restores, other git commands, other repositories, or a
worktree of the same repository. The permission gate refuses it in every
mode, before the auto-allow list (which allows a plain ``git switch``) is
consulted.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent import shared_checkout as SC
from delfin.agent.api_client import KitToolPermissions, _doc_executor

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args):
    return subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=True, capture_output=True,
        text=True).stdout.strip()


@pytest.fixture
def shared(tmp_path):
    d = tmp_path / "shared"
    d.mkdir()
    _git(d, "init", "-q", "-b", "main")
    (d / "README.md").write_text("x\n")
    (d / "sub").mkdir()
    (d / "sub" / "f.txt").write_text("y\n")
    _git(d, "add", "-A")
    _git(d, "commit", "-q", "-m", "init")
    _git(d, "branch", "feat")
    return d


@pytest.mark.parametrize("cmd", [
    "git switch feat",
    "git checkout feat",
    "git checkout -b new-one",
    "git switch -c new-one origin/main",
    "git checkout --detach feat",
    "git checkout -",
    "git fetch origin && git switch feat",
    "cd sub && git checkout feat",
    "git status\ngit switch feat",
])
def test_branch_moves_in_the_shared_checkout_are_refused(shared, cmd):
    out = SC.refusal(cmd, shared, shared=[shared])
    assert out.startswith("blocked:"), cmd
    assert "enter_worktree" in out and "▶ Test" in out


@pytest.mark.parametrize("cmd", [
    "git checkout -- README.md",
    "git checkout README.md",
    "git checkout feat -- README.md",
    "git status",
    "git commit -m 'git switch feat'",
    "git switch --help",
    "echo git switch feat",
])
def test_other_commands_are_not_refused(shared, cmd):
    assert SC.refusal(cmd, shared, shared=[shared]) == ""


def test_the_route_through_minus_c_is_refused_from_elsewhere(shared, tmp_path):
    other = tmp_path / "elsewhere"
    other.mkdir()
    assert SC.refusal(f"git -C {shared} switch feat", other, shared=[shared])


def test_another_repository_and_a_worktree_are_free(shared, tmp_path):
    other = tmp_path / "other"
    other.mkdir()
    _git(other, "init", "-q", "-b", "main")
    (other / "a").write_text("a")
    _git(other, "add", "-A")
    _git(other, "commit", "-q", "-m", "a")
    _git(other, "branch", "feat")
    assert SC.refusal("git switch feat", other, shared=[shared]) == ""
    wt = tmp_path / "wt"
    _git(shared, "worktree", "add", "-q", "--detach", str(wt), "main")
    assert SC.refusal("git switch feat", wt, shared=[shared]) == ""


def test_the_shared_checkout_comes_from_the_environment(shared, monkeypatch):
    monkeypatch.setenv(SC.ENV, str(shared / "sub"))
    assert shared.resolve() in SC.shared_checkouts()


class _Yes:
    def callback(self, *_a, **_k):
        return True


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_the_permission_gate_refuses_it_in_every_mode(shared, monkeypatch, mode):
    monkeypatch.setenv(SC.ENV, str(shared))
    perms = KitToolPermissions(workspace=shared, mode=mode,
                               confirm_callback=_Yes().callback)
    out = _doc_executor._run_permission_gate(
        "bash", {"command": "git switch feat"}, perms)
    assert out and "shared checkout" in out
    assert _git(shared, "branch", "--show-current") == "main"
    # the same command in a session worktree is not this rule's business
    out = _doc_executor._run_permission_gate(
        "bash", {"command": "git status"}, perms)
    assert not out or "shared checkout" not in out


def test_the_cli_backend_runner_refuses_it_too(shared, monkeypatch):
    from delfin.agent import sandbox
    monkeypatch.setenv(SC.ENV, str(shared))
    cfg = sandbox.SandboxConfig(mode="off", allow_network=False, timeout_s=10)
    res = sandbox.run_agent_command("git switch feat", shared, config=cfg)
    assert res.blocked and "shared checkout" in res.block_reason
    assert _git(shared, "branch", "--show-current") == "main"
