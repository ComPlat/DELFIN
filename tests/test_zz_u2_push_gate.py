"""Red probes for U2 phase-3 gate defects (operator findings, patch built).

Drives the gate through the real public path (`_DocToolExecutor` via
`_run_permission_gate`, exactly as tests/test_a_pull_request_is_opened_
under_the_push_grant.py does) and asserts the readiness-informed outcome
the phase-3 patch MUST produce: the gate refuses a push / pr create on
EVERY route (granted, dialog-approved, bypass) when readiness says so.
These bind the operator's seam -- api_client._push_readiness_rows, the ONLY
place the gate gets rows -- and verify the seam path refuses correctly.

Injection: monkeypatch api_client._push_readiness_rows with a fake
returning the standard row shape {"check","status","detail","fix"}.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import doctor as D
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _rows(workspace, *, remote_ok=True, gh_ok=True, raise_also=False):
    """Fake api_client._push_readiness_rows — a list of doctor rows."""
    rows = [
        {"check": "git installed", "status": D.PASS,
         "detail": "/usr/bin/git", "fix": ""},
        {"check": "git identity", "status": D.PASS,
         "detail": "user.name/user.email set", "fix": ""},
    ]
    if not remote_ok:
        rows.append({"check": "git remote", "status": D.WARN,
                     "detail": "unreachable",
                     "fix": "push from a node with outbound access, or ask "
                            "the user"})
    else:
        rows.append({"check": "git remote", "status": D.PASS,
                     "detail": "reachable, and it answers", "fix": ""})
    if gh_ok:
        rows.append({"check": "gh installed", "status": D.PASS,
                     "detail": "/usr/bin/gh", "fix": ""})
        rows.append({"check": "gh authenticated", "status": D.PASS,
                     "detail": "logged in", "fix": ""})
    else:
        rows.append({"check": "gh installed", "status": D.PASS,
                     "detail": "/usr/bin/gh", "fix": ""})
        rows.append({"check": "gh authenticated", "status": D.WARN,
                     "detail": "not logged in", "fix": "! gh auth login"})
    return rows


@pytest.fixture
def gate(tmp_path, monkeypatch):
    def _make(*, mode="bypassPermissions", grant=True, fake_rows=None,
              raiser=False, answer=None):
        monkeypatch.setattr(A, "_git_role", lambda: "maintainer")
        perms = KitToolPermissions(workspace=tmp_path, mode=mode)
        perms.push_grants = {"push": True} if grant else {}
        perms.confirm_callback = None
        if answer is not None:
            perms.confirm_callback = lambda *a, **k: answer
        eng = _DocToolExecutor.__new__(_DocToolExecutor)
        eng._permissions = perms

        if raiser:
            def _boom(*a, **k):
                raise RuntimeError("doctor exploded")
            monkeypatch.setattr(A, "_push_readiness_rows", _boom)
        elif fake_rows is not None:
            monkeypatch.setattr(A, "_push_readiness_rows",
                                lambda workspace=".", local_only=False:
                                fake_rows)

        def _run(cmd):
            return eng._run_permission_gate("bash", {"command": cmd}, perms)
        return _run
    return _make


def test_granted_push_refused_when_remote_unreachable(gate):
    """Operator defect 2: granted bypass push still consults readiness."""
    run = gate(grant=True, fake_rows=_rows(".", remote_ok=False, gh_ok=True))
    out = run("git push origin feature")
    assert isinstance(out, str), "granted push with unreachable remote must be refused"
    assert "outbound" in out or "ask the user" in out


def test_pr_create_refused_when_gh_logged_out(gate):
    """Operator defect 2: gh pr create needs gh logged in, refusing ! gh auth login."""
    run = gate(grant=True, fake_rows=_rows(".", remote_ok=True, gh_ok=False))
    out = run("gh pr create")
    assert isinstance(out, str), "gh pr create with gh logged out must be refused"
    assert "gh auth login" in out


def test_dialog_approved_push_still_checks_readiness(gate):
    """Operator defect 2: an approved-in-the-dialog push consults readiness."""
    run = gate(grant=True, mode="default", answer=True,
               fake_rows=_rows(".", remote_ok=False, gh_ok=True))
    out = run("git push origin feature")
    assert isinstance(out, str), "approved push must still refuse an unreachable remote"


def test_doctor_raising_is_a_refusal(gate):
    """A doctor crash must never let a push through."""
    run = gate(grant=True, raiser=True)
    out = run("git push origin feature")
    assert isinstance(out, str), "a doctor that raises must refuse, not allow"
