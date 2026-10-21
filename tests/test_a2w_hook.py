"""Package A2W, DORMANT: a pin of today's safety, not of a wiring.

A2W ships workspace_selfservice.classify (the three rated shapes) and
its settings contract (agent.workspace_selfservice.enabled, default
OFF, fail closed) WITHOUT wiring it into the gate. The person's
decision (operator, 2026-10-17): the exclusion deadlock found post-build
(pip/pipx/git/gh are on _CAGED_UNASKED_EXCLUDED_WORDS, and git switch
can run .git/hooks, a wheel can drop a .pth) means wiring needs its own
security review first.

This file pins that the rating exists but changes NOTHING at the gate:
with the switch EXPLICITLY on, each rated shape is STILL not run
unasked in the cage -- it reaches the excluded-word refusal/ask exactly
as on main. If a later change wires the hook, this file must be
rewritten along that decision (with its own review), not silently
passed.
"""

from __future__ import annotations

import pytest

from delfin import user_settings
from delfin.agent import api_client as A
from delfin.agent import workspace_selfservice as WSS
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


def _gate(ws, cmd, monkeypatch, enabled, mode="default", attended=True):
    """`enabled=None` leaves the settings file WITHOUT the block (the
    fail-closed default applies); any other value is written as
    agent.workspace_selfservice.enabled into the REAL settings file under
    the redirected HOME."""
    asked = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    if enabled is not None:
        user_settings.save_settings(
            {"agent": {WSS.SETTINGS_KEY: {"enabled": enabled}}})
    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=_cb if attended else None)
    admitted = []
    verdict = _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms, caged_admission=admitted)
    return verdict, asked, admitted


# --- the rating exists, but is dormant ---------------------------------------

def test_pip_rated_shape_is_classified_but_not_wired(ws, monkeypatch):
    cmd = (".venv/bin/pip install -r requirements.txt --require-hashes "
           "--only-binary=:all:")
    assert WSS.classify(cmd, str(ws)) == "auto", "the rating exists"
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, cmd, monkeypatch, True)
    # DORMANT: even with the switch on, nothing runs unasked.
    assert admitted == [], "the rating must not admit anything unwired"
    assert asked == [], "nothing asked either: the excluded-word rule holds"


def test_git_switch_rated_shape_is_classified_but_not_wired(ws, monkeypatch):
    assert WSS.classify("git switch feature", str(ws)) == "auto"
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git switch feature",
                                     monkeypatch, True)
    assert admitted == [], "the rating must not admit anything unwired"
    assert asked == [], "nothing asked either: the excluded-word rule holds"


def test_tool_version_rated_shape_is_classified_but_not_wired(
        ws, monkeypatch):
    assert WSS.classify("git --version", str(ws)) == "auto"
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, "git --version",
                                     monkeypatch, True)
    assert admitted == [], "the rating must not admit anything unwired"


# --- the settings contract stands on its own ----------------------------------

def test_switch_default_is_off(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    assert user_settings.load_settings().get("agent", {}).get(
        WSS.SETTINGS_KEY) is None, "no block on a fresh settings file"
    assert WSS.settings_from_payload({}) == {"enabled": False}, \
        "fail closed"


@pytest.mark.parametrize("junk", [False, "yes", 1])
def test_switch_non_bool_or_false_is_off(junk):
    assert WSS.settings_from_payload(
        {"agent": {WSS.SETTINGS_KEY: {"enabled": junk}}}) == {"enabled": False}
