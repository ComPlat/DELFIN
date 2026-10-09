"""Phase 3 of U2: the command gate asks readiness before a push.

Phase 2 added ``doctor.ready_for_push`` -- one readiness function for
"can this host push": git on PATH, identity set, remote reachable, gh on
PATH, gh logged in. Phase 3 wires it into the command gate: the same gate
that guards ``git push`` / ``gh pr create`` now asks readiness before one
gets through, and refuses with the remedy when the host cannot push --
instead of letting the command fail and diagnosing afterwards.

The check runs twice, because the two moments need different probes:
- A cheap *local* probe before asking the user. If git is missing or no
  identity is set, there is nothing to ask about: refuse immediately with
  the remedy, and never prompt. It must not touch the network (no
  ``git ls-remote``), so a push the user might decline costs nothing.
- The *full* readiness check at the surrender point, right before the
  command would run. By then the user has authorized the push, so the
  network probe (``git ls-remote``) is warranted: refuse with the remedy
  if the remote will not answer.

Tests need no live git/gh and no network: they patch ``shutil.which``
(git/gh presence) and ``subprocess.run`` (config / ls-remote / auth
status) in the doctor module namespace, the doctor's own injection point
(see tests/test_delfin_checks_the_tools_it_requires.py).
"""

from __future__ import annotations

import subprocess as _sp

import pytest

from delfin.agent import api_client as A
from delfin.agent import doctor as D
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _boom(*a, **k):
    raise OSError("no subprocess in a test")


@pytest.fixture(autouse=True)
def _maintainer_role(monkeypatch):
    # These tests are about the capability check, not the role; the role
    # has its own file. A maintainer may push to any branch, so the gate
    # reaches the readiness check on the happy path.
    monkeypatch.setattr(A, "_git_role", lambda: "maintainer")
    # Default: a host that lacks every tool, and a runner that refuses.
    # Individual tests override via *_host.
    monkeypatch.setattr(D.shutil, "which", lambda n: None)
    monkeypatch.setattr(D.subprocess, "run", _boom)


def _perms(tmp_path, mode="bypassPermissions"):
    return KitToolPermissions(workspace=tmp_path, mode=mode)


def _gate(perms, cmd):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


def _result(rc_ok, argv=None, stderr=""):
    if rc_ok is True:
        return _sp.CompletedProcess(argv or [], 0, stdout="", stderr="")
    if isinstance(rc_ok, int):
        return _sp.CompletedProcess(argv or [], rc_ok, stdout="",
                                    stderr=stderr)
    raise _boom()


def _field(field, identity, helper):
    if field == "credential.helper":
        return helper
    if identity is True and field in ("user.name", "user.email"):
        return f"User {field}"
    return 1  # git config --get exits 1 when unset


def _host(monkeypatch, *, git=True, identity=True, remote=True, gh=True,
          authed=True, helper=""):
    """Point the doctor probes at a concrete host. The default host has
    everything and can push; callers turn pieces off to make a broken one."""
    def which(n):
        if n == "git" and git:
            return "/usr/bin/git"
        if n == "gh" and gh:
            return "/usr/bin/gh"
        return None

    def run(argv, **kw):
        argv = list(argv)
        if argv[0] == "git":
            if argv[1:3] == ["config", "--get"]:
                return _result(_field(argv[2], identity, helper), argv)
            if argv[1:3] == ["ls-remote"]:
                return _result(remote if remote is not False else 128, argv,
                               "" if remote else "fatal: could not read from remote repository")
            return _result(remote if remote is not False else 128, argv,
                           "" if remote else "fatal: remote error")
        if argv[0] == "gh":
            return _result(authed, argv,
                           "" if authed else "To get started with GitHub CLI")
        raise _boom()

    monkeypatch.setattr(D.shutil, "which", which)
    monkeypatch.setattr(D.subprocess, "run", run)


def test_host_without_git_is_refused_before_it_prompts(tmp_path, monkeypatch):
    """A host with no git must refuse before the consent dialog, and must
    not reach the full (network) probe either."""
    _host(monkeypatch, git=False)
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "git push origin main")

    assert msg and "git is not installed" in msg
    assert not asked, "a host that cannot push must not prompt"


def test_host_without_identity_is_refused_before_it_prompts(
        tmp_path, monkeypatch):
    _host(monkeypatch, git=True, identity=False, gh=True, authed=True)
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "git push origin main")

    assert msg and "identity" in msg.lower()
    assert not asked


def test_gh_logged_out_is_refused_before_pr_create_prompts(
        tmp_path, monkeypatch):
    _host(monkeypatch, git=True, identity=True, gh=True, authed=False)
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "gh pr create --base main")

    assert msg and "gh" in msg.lower()
    assert not asked


def test_a_ready_host_reaches_the_granted_push_unhindered(
        tmp_path, monkeypatch):
    """Everything present and the user asked: the check must not invent a
    blockage. ``remote=True`` so the surrender probe also passes."""
    _host(monkeypatch)
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "push it", new_request=True)

    assert _gate(perms, "git push origin main") is None


def test_a_ready_remote_without_a_grant_still_refuses_as_before(
        tmp_path, monkeypatch):
    """The capability check must not weaken the existing one-request gate:
    a push the user did not ask for is still refused, capability aside."""
    _host(monkeypatch)
    perms = _perms(tmp_path, mode="default")

    msg = _gate(perms, "git push origin main")

    assert (msg and ("has not asked for a push" in msg or "blocked" in msg))


def test_remote_unreachable_refuses_at_surrender_despite_a_grant(
        tmp_path, monkeypatch):
    """A granted push whose remote will not answer is refused when the
    full check runs, with the remedy -- not sent to git to fail."""
    _host(monkeypatch, git=True, identity=True, remote=False, gh=True,
          authed=True)
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "push it", new_request=True)

    msg = _gate(perms, "git push origin main")

    assert msg and "remote" in msg


def test_pr_create_without_gh_but_awarded_grant_refuses(
        tmp_path, monkeypatch):
    """grant=1 must not let a host without gh publish a pull request."""
    _host(monkeypatch, git=True, identity=True, remote=True, gh=False)
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "open the pr", new_request=True)

    msg = _gate(perms, "gh pr create --base main")

    assert msg and "gh" in msg.lower()


def test_host_without_git_never_reaches_the_network(
        tmp_path, monkeypatch):
    """Regression: the local pre-ask probe must not run git ls-remote."""
    seen_ls_remote = False
    monkeypatch.setattr(D.shutil, "which", lambda n: None)

    def guard(argv, **kw):
        nonlocal seen_ls_remote
        argv = list(argv)
        if argv[0] == "git" and "ls-remote" in argv:
            seen_ls_remote = True  # must stay False
        if argv[0] == "git":
            return _result(1, argv, "no git here")
        if argv[0] == "gh":
            return _result(1, argv, "no gh here")
        raise _boom()

    monkeypatch.setattr(D.subprocess, "run", guard)
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    _gate(perms, "git push origin main")
    assert not asked
    assert not seen_ls_remote


def test_git_push_is_allowed_when_gh_is_absent(monkeypatch, tmp_path):
    """git push publishes over SSH/HTTPS and does not need gh: an absent
    gh must not be refused (the over-blocking the operator attacked)."""
    _host(monkeypatch, git=True, identity=True, remote=True, gh=False)
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg is None, f"gh-absent host must be allowed: {msg}"


def test_granted_git_push_refused_when_remote_unreachable(
        monkeypatch, tmp_path):
    """A pre-set grant still must not let the command run when the remote
    will not answer: the surrender check runs on the grant route too, and
    refuses with the remedy."""
    _host(monkeypatch, git=True, identity=True, remote=False, gh=True)
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg and "remote" in msg.lower()
    # the refusal names a remedy, not just the failure
    assert "--" in msg or "push from" in msg or "reach" in msg.lower()


def test_gh_pr_create_logged_out_refuses_naming_the_login(
        monkeypatch, tmp_path):
    """gh pr create needs gh logged in; the refusal names the login."""
    _host(monkeypatch, git=True, identity=True, remote=True,
          gh=True, authed=False)
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "gh pr create --base main")
    assert msg and "gh auth login" in msg


def test_a_doctor_that_raises_is_a_refusal(monkeypatch, tmp_path):
    """Never let a push through when the readiness discovery itself blows
    up: a doctor that raises is a WARN refusal with a fix, not a pass."""
    def _raises(*a, **k):
        raise RuntimeError("doctor is broken")
    monkeypatch.setattr(D, "_check_git_tooling", _raises)
    _host(monkeypatch, git=True, identity=True, remote=True, gh=True)
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg  # a non-None refusal, never None
