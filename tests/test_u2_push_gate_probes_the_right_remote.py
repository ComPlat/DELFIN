"""The push gate's readiness probe, driven with real git against local remotes.

The other U2 tests replace ``api_client._push_readiness_rows`` with fixed
rows; these put the real one back and drive the gate
(``_run_permission_gate``) against bare repositories on disk, so the rows
come from the doctor's actual probes. No network: every remote is a local
path, and nothing is pushed -- the gate only decides.

What is pinned, each measured as a defect of the first build:

* the probe runs in the host process, outside the sandbox the agent's own
  commands run in, and must not execute a command the checkout's
  ``.git/config`` names (``remote.<name>.uploadpack`` ran);
* a push is probed in the directory it runs in (``git -C``, the ``cwd``
  argument) and against the remote it names (``upstream``, a URL), not
  always ``origin`` in the workspace;
* an empty remote (no HEAD yet) is reachable, not a refusal;
* an unset commit identity does not refuse a push, which makes no commit;
* a plain ``git push`` never runs ``gh``;
* and the control: an unreachable remote and an unknown remote name are
  still refused.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import doctor as D
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

_REAL_ROWS = A._push_readiness_rows

pytestmark = pytest.mark.skipif(
    subprocess.run(["git", "--version"], capture_output=True).returncode != 0,
    reason="git is required to build the local remotes")


def _git(*args, cwd):
    subprocess.run(["git", *args], cwd=str(cwd), check=True,
                   capture_output=True)


@pytest.fixture
def host(tmp_path, monkeypatch):
    """A repository with one commit, a bare remote with that commit, and
    the user's git configuration replaced by one holding only an
    identity."""
    cfg = tmp_path / "gitconfig"
    cfg.write_text("[user]\n\tname = t\n\temail = t@t\n")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(cfg))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    monkeypatch.setattr(A, "_push_readiness_rows", _REAL_ROWS)
    monkeypatch.setattr(A, "_git_role", lambda: "maintainer")
    bare = tmp_path / "remote.git"
    _git("init", "-q", "--bare", str(bare), cwd=tmp_path)
    repo = tmp_path / "ws" / "repo"
    repo.mkdir(parents=True)
    _git("init", "-q", cwd=repo)
    _git("commit", "-q", "--allow-empty", "-m", "x", cwd=repo)
    _git("push", "-q", str(bare), "HEAD:refs/heads/base", cwd=repo)
    _git("symbolic-ref", "HEAD", "refs/heads/base", cwd=bare)
    return {"tmp": tmp_path, "repo": repo, "bare": bare, "cfg": cfg}


def _gate(workspace, cmd, **args):
    perms = KitToolPermissions(workspace=Path(workspace),
                               mode="bypassPermissions")
    perms.push_grants = {"push": 1}
    return _DocToolExecutor()._run_permission_gate(
        "bash", dict(args, command=cmd), perms)


def test_the_probe_runs_no_command_the_checkout_config_names(host):
    _git("remote", "add", "origin", str(host["bare"]), cwd=host["repo"])
    marker = host["tmp"] / "uploadpack-ran"
    _git("config", "remote.origin.uploadpack",
         f"touch {marker}; git-upload-pack", cwd=host["repo"])
    out = _gate(host["repo"], "git push origin HEAD:feature")
    assert not marker.exists(), "the gate executed the checkout's uploadpack"
    assert out is None, out


def test_a_push_from_a_subdirectory_is_probed_where_it_runs(host):
    _git("remote", "add", "origin", str(host["bare"]), cwd=host["repo"])
    ws = host["repo"].parent  # not a repository itself
    assert _gate(ws, "git -C repo push origin HEAD:feature") is None
    assert _gate(ws, "git push origin HEAD:feature", cwd="repo") is None
    assert _gate(ws, "cd repo && git push origin HEAD:feature") is None


def test_a_push_to_another_remote_probes_that_remote(host):
    _git("remote", "add", "upstream", str(host["bare"]), cwd=host["repo"])
    assert _gate(host["repo"], "git push upstream HEAD:feature") is None
    assert _gate(host["repo"],
                 f"git push {host['bare']} HEAD:feature") is None


def test_an_empty_remote_is_reachable(host):
    empty = host["tmp"] / "empty.git"
    _git("init", "-q", "--bare", str(empty), cwd=host["tmp"])
    _git("remote", "add", "origin", str(empty), cwd=host["repo"])
    assert _gate(host["repo"], "git push -u origin HEAD:feature") is None


def test_an_unset_identity_does_not_refuse_a_push(host):
    _git("remote", "add", "origin", str(host["bare"]), cwd=host["repo"])
    host["cfg"].write_text("")
    out = _gate(host["repo"], "git push origin HEAD:feature")
    assert out is None, out


def test_a_plain_git_push_never_runs_gh(host, monkeypatch):
    _git("remote", "add", "origin", str(host["bare"]), cwd=host["repo"])
    bin_dir = host["tmp"] / "bin"
    bin_dir.mkdir()
    marker = host["tmp"] / "gh-ran"
    gh = bin_dir / "gh"
    gh.write_text(f"#!/bin/sh\ntouch {marker}\nexit 0\n")
    gh.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    assert _gate(host["repo"], "git push origin HEAD:feature") is None
    assert not marker.exists(), "gh auth status ran for a plain git push"


def test_an_unreachable_or_unknown_remote_is_still_refused(host):
    _git("remote", "add", "origin", str(host["tmp"] / "missing.git"),
         cwd=host["repo"])
    out = _gate(host["repo"], "git push origin HEAD:feature")
    assert out and out.startswith("blocked: git remote"), out
    out = _gate(host["repo"], "git push nosuch HEAD:feature")
    assert out and "no remote named 'nosuch'" in out, out


def test_a_diagnosed_refusal_carries_git_s_own_words(host):
    """The remedy says to read what git printed; the refusal must show it."""
    _git("remote", "add", "origin", str(host["tmp"] / "missing.git"),
         cwd=host["repo"])
    row = [r for r in D._check_push({"workspace": str(host["repo"])})
           if r["check"] == "git remote"][0]
    assert row["status"] == D.WARN
    assert "does not appear to be a git repository" in row["detail"], row
