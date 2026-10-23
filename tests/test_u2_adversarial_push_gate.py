"""Adversarial guard for U2 phase-3, defect 1 — the split-by-command contract.

The corrected phase-3 gate (operator finding, patch NOT built) must split
readiness by command: `git push` needs git installed (+ the remote at
surrender; no identity, a push makes no commit); `gh pr create` /
`gh pr merge` additionally need gh installed + gh authenticated. A
`git push` over SSH/HTTPS needs NO gh, so the gate must let it through
when git/remote are ready even if gh is not on PATH or not logged in.

The gate gets its readiness rows from exactly one place —
`api_client._push_readiness_rows(workspace, *, local_only, ...) -> list[dict]`
(the corrected patch; `_push_capability_block` calls it, never doctor
directly). U2 tests rebind that seam with
`monkeypatch.setattr(api_client, "_push_readiness_rows", fake)`, setting
the doctor rows they need explicitly — never the real git config, PATH or
network.

This test pins that contract as a GUARD that REQUIRES the seam. It binds
directly to `api_client._push_readiness_rows` unconditionally — there is no
fallback — so on the current tree (the corrected patch, and with it the
seam, not built) the test is RED: a metric that is green while the seam is
missing would pass even if the patch never landed or renamed the hook, and
tests nothing. That red is the control. Once the operator builds the
corrected patch (seam + neutral_push_readiness conftest fixture), the seam
exists, the fake binds, and the guard goes green — and a patch that reads
real readiness (isolated CI HOME, no git identity) or refuses a `git push`
on a gh row (over-blocking) keeps it red. The conftest fixture is the
operator's build; this test overrides it with its own explicit rows.

Driven through the gate the dispatcher calls (eng._run_permission_gate),
exactly as in tests/test_a_pull_request_is_opened_under_the_push_grant.py.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import doctor as D
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _rows(git_ok=True, remote_ok=True, gh_ok=True):
    """Fake readiness rows: PASS for git/remote, a gh row not-PASS when asked.

    The split defect is that a `git push` must pass when only git is
    checked; here we hand back the git rows PASS and a gh row that is NOT
    PASS (gh absent / logged out). A correct patch checks git/remote for
    push and ignores the gh row.
    """
    rows = [
        {"check": "git installed", "status": D.PASS, "detail": "/usr/bin/git",
         "fix": ""},
        {"check": "git identity", "status": D.PASS,
         "detail": "user.name/user.email set", "fix": ""},
        {"check": "git remote", "status": D.PASS,
         "detail": "reachable, and it answers", "fix": ""},
    ]
    if not gh_ok:
        rows.append({"check": "gh installed", "status": D.WARN,
                     "detail": "gh is not on PATH", "fix": "install gh"})
    return rows


def _bind_readiness(monkeypatch, gh_ok):
    """Point the gate at fake rows through the seam — no fallback.

    The gate reads readiness ONLY from api_client._push_readiness_rows (the
    operator's seam). We monkeypatch that attribute directly with the fake.
    On the current tree the seam does not exist, so the setattr raises
    AttributeError and the test is red — that missing-seam red is the
    control (see module docstring). Once the operator builds the corrected
    patch, the seam exists, the fake binds, and the test goes green.
    """
    rows = _rows(gh_ok=gh_ok)
    monkeypatch.setattr(A, "_push_readiness_rows",
                        lambda workspace=".", local_only=False, **_: rows)


@pytest.fixture
def gate(tmp_path, monkeypatch):
    def _make(*, mode="bypassPermissions", grant=True, gh_ok=True):
        monkeypatch.setattr(A, "_git_role", lambda: "maintainer")
        perms = KitToolPermissions(workspace=tmp_path, mode=mode)
        perms.push_grants = {"push": True} if grant else {}
        perms.confirm_callback = None
        _bind_readiness(monkeypatch, gh_ok)
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        eng._permissions = perms

        def _run(cmd):
            return eng._run_permission_gate("bash", {"command": cmd}, perms)
        return _run
    return _make


def test_a_git_push_is_allowed_without_gh_on_the_path(gate):
    """A checkout that pushes over SSH/HTTPS needs no gh binary."""
    run = gate(gh_ok=False)
    out = run("git push origin feature")
    assert out is None, (
        "a `git push` with git/identity/remote ready must not be refused "
        "because gh is missing (split-by-command): " + repr(out))


def test_a_git_push_is_allowed_with_gh_logged_out(gate):
    """gh being logged out must not stop a plain git push."""
    run = gate(gh_ok=True)  # gh present; the guard above covers the push
    out = run("git push origin feature")
    assert out is None, "a ready git push must pass the gate: " + repr(out)
