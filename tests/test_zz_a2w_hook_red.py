"""Package A2W reviewer (nacht-s12): red test for the ONE hook in
api_client.py -- at the bash decision, after every hard tier and before
the ask, workspace_selfservice.classify rates "auto" shapes to run
WITHOUT a question: only with the cage on, only inside the own
workspace; an escape variant still asks.

State on the hookless main: a rated pip --require-hashes install and a
git switch ask (or dead-end "not on the auto-allow list" when unattended)
-- the gate's caged rule refuses package managers and git by text scan,
so the shapes the classifier rates are unreachable today. That is the
defect this test measures red; the operator-built hook turns it green.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    d.mkdir()
    (d / "requirements.txt").write_text("numpy==1.26.4\n")
    venv = d / ".venv" / "bin"
    venv.mkdir(parents=True)
    (venv / "pip").write_text("#!/bin/sh\necho pip-ran\n")
    return d


def _cage(monkeypatch, kind: str) -> None:
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (kind, "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: kind == "bwrap")
    monkeypatch.setattr(A, "_landlock_functional", lambda: kind == "bwrap")


def _gate(ws, cmd, mode="default", attended=False):
    asked = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=_cb if attended else None)
    admitted = []
    verdict = _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms, caged_admission=admitted)
    return verdict, asked, admitted


# --- the rated shapes must RUN without a question ----------------------------
# xfail until the operator builds the ONE hook (api_client ~21705-21724,
# caged-admission consults workspace_selfservice.classify): on hookless
# main these are NOT admitted -- that red is the control. When the hook
# lands, strict=False keeps them green without any edit here.

@pytest.mark.xfail(reason="hook not built yet: rated shape not admitted "
                          "on hookless main", strict=False)
@pytest.mark.parametrize("attended", [True, False])
def test_pip_require_hashes_runs_unasked_when_caged(ws, monkeypatch,
                                                    attended):
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(
        ws, ".venv/bin/pip install -r requirements.txt --require-hashes",
        attended=attended)
    assert admitted == [True], verdict
    assert asked == []
    assert verdict is None, verdict


@pytest.mark.xfail(reason="hook not built yet: rated shape not admitted "
                          "on hookless main", strict=False)
@pytest.mark.parametrize("attended", [True, False])
def test_git_switch_own_branch_runs_unasked_when_caged(ws, monkeypatch,
                                                       attended):
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git switch feature",
                                     attended=attended)
    assert admitted == [True], verdict
    assert asked == []
    assert verdict is None, verdict


@pytest.mark.xfail(reason="hook not built yet: rated shape not admitted "
                          "on hookless main", strict=False)
def test_tool_version_runs_unasked_when_caged(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "tar --version")
    assert admitted == [True], verdict
    assert asked == []
    assert verdict is None, verdict


# --- the escape variants must still ask / be refused -------------------------

def test_pip_r_outside_still_not_admitted(ws, monkeypatch):
    # Measured (probe through the gate): the gate returns None with no
    # dialog and no admission -- the executor's own write gates catch the
    # outside -r file later. The hook must not change that: classify
    # returns None for this shape, so nothing may be admitted here.
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(
        ws, f".venv/bin/pip install -r {ws.parent / 'evil.txt'} --require-hashes")
    assert admitted == []


def test_git_switch_redirect_outside_is_denied(ws, monkeypatch):
    # Measured: a deny-pattern (">\\s*/etc/") refuses this before any
    # dialog. The hook must keep that: classify returns None for the
    # redirect variant, so the refusal stands.
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(
        ws, "git switch feature > /etc/a2w_escape_marker")
    assert admitted == []
    assert "deny-pattern" in verdict


# --- the hook only adds; it never removes ------------------------------------

def test_cage_off_the_rated_shapes_still_ask(ws, monkeypatch):
    """Sandbox off: the classifier's verdict must not matter -- the old
    gate decides (ask/refuse), nothing runs unasked."""
    _cage(monkeypatch, "none")
    for cmd in (".venv/bin/pip install -r requirements.txt --require-hashes",
                "git switch feature", "tar --version"):
        verdict, asked, admitted = _gate(ws, cmd)
        assert admitted == []


def test_a_foreign_branch_switch_is_not_auto(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git switch origin/main")
    assert admitted == []
