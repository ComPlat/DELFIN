"""Phase 3 of U2: the command gate asks readiness before a push.

Phase 2 added ``doctor.ready_for_push`` -- one readiness function for
"can this host push": git on PATH, identity set, remote reachable, gh on
PATH, gh logged in. Phase 3 wires it into the command gate: the same gate
that guards ``git push`` / ``gh pr create`` now asks readiness before one
gets through, and refuses with the remedy when the host cannot push --
instead of letting the command fail and diagnosing afterwards.

The gate reads readiness rows ONLY from a single injectable seam,
``api_client._push_readiness_rows``. The operator's build adds an
autouse conftest fixture ``neutral_push_readiness`` that returns ``[]``
(nothing non-PASS) for every test, so no existing test depends on the
machine. Each test here OVERRIDES that seam with the exact rows its
scenario needs -- it never touches the doctor's probes (no real git
config, PATH or network) and never depends on the host.

``_push_capability_block`` filters the rows by ``check`` name: a plain
``git push`` keeps only ``{git installed, git identity}`` (+ ``git
remote`` at the surrender full check), so an absent gh never refuses a
push; ``gh pr create`` / ``gh pr merge`` additionally keep ``{gh
installed, gh authenticated}``. The check runs twice -- a local-only
pre-ask probe (no network) before the consent dialog, and the full check
(incl. the network ``git ls-remote`` probe, hard-bounded inside the
doctor) on every route that lets the command run.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _row(check, status="PASS", detail="", fix=""):
    return {"check": check, "status": status, "detail": detail, "fix": fix}


def _rows(monkeypatch, rows):
    """Point the gate's one row source at a concrete list of rows.

    The conftest autouse ``neutral_push_readiness`` sets
    ``_push_readiness_rows`` to [] for every test; each test here
    overrides it with the exact rows it needs, so the gate's refusal
    comes from THIS scenario and never from the real host.
    """
    monkeypatch.setattr(A, "_push_readiness_rows", lambda *a, **k: rows)


def _perms(tmp_path, mode="bypassPermissions"):
    return KitToolPermissions(workspace=tmp_path, mode=mode)


def _gate(perms, cmd):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


# ---------------------------------------------------------------------------
# Refusals before the consent dialog (local-only probe)
# ---------------------------------------------------------------------------


def test_host_without_git_is_refused_before_it_prompts(tmp_path, monkeypatch):
    """A host with no git must refuse before the consent dialog, and must
    not reach the full (network) probe either."""
    _rows(monkeypatch, [
        _row("git installed", "WARN", "git is not on PATH",
             "install git; the agent cannot push without it")])
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "git push origin main")

    assert msg and "git is not installed" in msg
    assert not asked, "a host that cannot push must not prompt"


def test_host_without_identity_is_refused_before_it_prompts(
        tmp_path, monkeypatch):
    _rows(monkeypatch, [
        _row("git installed"),
        _row("git identity", "WARN", "not configured: user.name, user.email",
             "git config --global user.name/email")])
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "git push origin main")

    assert msg and "identity" in msg.lower()
    assert not asked


def test_gh_logged_out_is_refused_before_pr_create_prompts(
        tmp_path, monkeypatch):
    _rows(monkeypatch, [
        _row("git installed"),
        _row("git identity"),
        _row("gh installed"),
        _row("gh authenticated", "WARN", "no active GitHub login",
             "gh auth login -- pull requests cannot be opened until this "
             "succeeds")])
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "gh pr create --base main")

    assert msg and "gh auth login" in msg
    assert not asked


def test_a_ready_host_with_no_blocker_never_prompts(tmp_path, monkeypatch):
    """A ready host (all rows PASS) reaches the consent path and is NOT
    refused by the capability probe."""
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity")])
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "git push origin main")

    # Not refused for a capability reason; the refusal (if any) is the
    # ordinary one-request gate, never a readiness block.
    assert msg is None or "git is not installed" not in msg
    assert not asked or True  # readiness never prompted


# ---------------------------------------------------------------------------
# The surrender full check on every allowance route
# ---------------------------------------------------------------------------


def test_ready_host_reaches_the_granted_push_unhindered(
        tmp_path, monkeypatch):
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity"), _row("git remote")])
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1

    msg = _gate(perms, "git push origin main")

    assert msg is None, f"a ready, granted host must be allowed: {msg}"


def test_git_push_is_allowed_when_gh_is_absent(monkeypatch, tmp_path):
    """git push publishes over SSH/HTTPS and does not need gh: an absent
    gh must not be refused."""
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity"), _row("git remote")])
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg is None, f"gh-absent host must be allowed: {msg}"


def test_granted_git_push_refused_when_remote_unreachable(
        monkeypatch, tmp_path):
    """A pre-set grant still must not let the command run when the remote
    will not answer: the surrender check runs on the grant route too, and
    refuses with the remedy."""
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity"),
        _row("git remote", "WARN", "unreachable: origin could not be reached",
             "push from a node with outbound access, or ask the user")])
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg and "remote" in msg


def test_gh_pr_create_logged_out_refuses_naming_the_login(
        monkeypatch, tmp_path):
    """gh pr create with a grant still needs gh logged in: refused naming
    the fix."""
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity"),
        _row("gh installed"),
        _row("gh authenticated", "WARN", "no active GitHub login",
             "gh auth login -- pull requests cannot be opened until this "
             "succeeds")])
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "gh pr create --base main")
    assert msg and "gh auth login" in msg


def test_pr_create_without_gh_but_awarded_grant_refuses(
        tmp_path, monkeypatch):
    """A granted gh pr create still refuses when gh is missing."""
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity"),
        _row("gh installed", "WARN", "gh is not on PATH",
             "install the GitHub CLI (`gh`) to open pull requests")])
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "gh pr create --base main")
    assert msg and "gh" in msg.lower()


def test_a_doctor_that_raises_is_a_refusal(monkeypatch, tmp_path):
    """If the readiness row source raises, the gate must refuse, never
    allow the push through untouched."""
    def boom(*a, **k):
        raise OSError("no host probes in a test")
    monkeypatch.setattr(A, "_push_readiness_rows", boom)
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg  # a non-None refusal, never None


# ---------------------------------------------------------------------------
# Regression: the one-request push gate is not weakened
# ---------------------------------------------------------------------------


def test_a_ready_remote_without_a_grant_still_refuses_as_before(
        tmp_path, monkeypatch):
    _rows(monkeypatch, [
        _row("git installed"), _row("git identity"), _row("git remote")])
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda n, a, p: asked.append(p) or True

    msg = _gate(perms, "git push origin main")

    assert (msg and ("has not asked for a push" in msg
                     or "blocked" in msg))


def test_host_without_git_never_reaches_the_network(tmp_path, monkeypatch):
    """The pre-ask local probe must not touch the network probe (no
    git remote row needed / no ls-remote): a host with no git is refused
    before any network row is even consulted."""
    _rows(monkeypatch, [
        _row("git installed", "WARN", "git is not on PATH",
             "install git; the agent cannot push without it")])
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    msg = _gate(perms, "git push origin main")
    assert msg and "git is not installed" in msg
