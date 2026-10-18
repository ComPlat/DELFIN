"""Package A2W reviewer (nacht-s12): DORMANT pin -- the workspace
self-service classifier is shipped UNWIRED (operator ruling 2026-10-17).

pip/pipx/git/gh are on _CAGED_UNASKED_EXCLUDED_WORDS (api_client.py
:28446-28449) and stay there on purpose: git switch can run .git/hooks, a
binary wheel can drop a .pth that runs at the next python start. So even
with the person's switch ON, a rated shape is NOT admitted at the bash
gate -- the caged rule's excluded-word tier refuses/asks first, exactly as
on main. This file pins that safety in both switch states; if a later
wiring ever admits one of these shapes unasked, these tests go red and
demand the wiring's own security review.
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


def _enable(monkeypatch, _tmp_path_unused, value):
    """Turn agent.workspace_selfservice.enabled to `value` (True/False/
    "yes"/None = no agent block at all). Patched at the SAME loader the
    (future) hook would read through (delfin.user_settings.load_settings),
    the pattern tests/test_a_caged_command_runs_without_asking.py:290
    uses; no real ~/.delfin is touched."""
    import delfin.user_settings as US

    def _payload(*_a, **_k):
        if value is None:
            return {}
        return {"agent": {"workspace_selfservice": {"enabled": value}}}

    monkeypatch.setattr(US, "load_settings", _payload)


# --- the rated shapes are NOT run unasked, switch ON or OFF ------------------

@pytest.mark.parametrize("enabled", [True, False, None])
def test_pip_shape_is_not_admitted_even_enabled(ws, monkeypatch, enabled):
    _enable(monkeypatch, ws.parent.parent, enabled)
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(
        ws, ".venv/bin/pip install -r requirements.txt --require-hashes"
            " --only-binary=:all:")
    assert admitted == [], verdict


@pytest.mark.parametrize("enabled", [True, False, None])
def test_git_switch_is_not_admitted_even_enabled(ws, monkeypatch, enabled):
    _enable(monkeypatch, ws.parent.parent, enabled)
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git switch feature")
    assert admitted == [], verdict


@pytest.mark.parametrize("enabled", [True, False, None])
def test_git_version_is_not_admitted_even_enabled(ws, monkeypatch, enabled):
    _enable(monkeypatch, ws.parent.parent, enabled)
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git --version")
    assert admitted == [], verdict


# --- the switch state itself is honest (module contract, fail closed) --------

def test_switch_defaults_off(ws):
    import delfin.agent.workspace_selfservice as wss
    assert wss.selfservice_defaults() == {"enabled": False}


# --- escapes stay refused whatever the switch says ---------------------------

def test_pip_r_outside_still_not_admitted(ws, monkeypatch):
    _enable(monkeypatch, ws.parent.parent, True)
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(
        ws, f".venv/bin/pip install -r {ws.parent / 'evil.txt'}"
            " --require-hashes --only-binary=:all:")
    assert admitted == []


def test_git_switch_redirect_outside_is_denied(ws, monkeypatch):
    # Measured: a deny-pattern (">\\s*/etc/") refuses this before any
    # dialog. A later wiring must keep that: classify returns None for
    # the redirect variant, so no admission can ever cover it.
    _enable(monkeypatch, ws.parent.parent, True)
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(
        ws, "git switch feature > /etc/a2w_escape_marker")
    assert admitted == []
    assert "deny-pattern" in verdict


# --- cage off: the classifier's verdict must not matter ----------------------

def test_cage_off_the_rated_shapes_still_ask(ws, monkeypatch):
    """Sandbox off: nothing runs unasked, whatever the switch says."""
    _enable(monkeypatch, ws.parent.parent, True)
    _cage(monkeypatch, "none")
    for cmd in (".venv/bin/pip install -r requirements.txt"
                " --require-hashes --only-binary=:all:",
                "git switch feature", "git --version"):
        verdict, asked, admitted = _gate(ws, cmd)
        assert admitted == []


def test_a_foreign_branch_switch_is_not_auto(ws, monkeypatch):
    _enable(monkeypatch, ws.parent.parent, True)
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git switch origin/main")
    assert admitted == []
