"""A request of the agent for more permissions is flagged, asked twice, and
never put right after a no.

Measured 2026-10-23 on the TUI: the user refused a write, and the model
answered with remember_permission_bundle -- a workspace root under /tmp
and ten auto-allow patterns into the user's settings -- behind one
ordinary question. DELFIN is used by people who do not read regexes: what
a yes means has to be said in plain words, above the details.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import kit_settings as ks


@pytest.fixture
def setup(tmp_path, monkeypatch):
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "user.json")
    ws = tmp_path / "ws"
    ws.mkdir()
    asked: list = []
    answers: list = []

    def confirm(name, args, preview):
        asked.append(preview)
        return answers.pop(0) if answers else False

    perms = A.KitToolPermissions(workspace=ws, mode="default",
                                 confirm_callback=confirm)
    perms.turn_token = A._next_turn_token()
    return A._DocToolExecutor(), perms, asked, answers, tmp_path


def _allow(ex, perms, pattern=r"^make\b"):
    return json.loads(ex._execute_remember_permission(
        {"kind": "allow_pattern", "value": pattern, "rationale": "build it"},
        perms))


def test_the_question_says_plainly_what_a_yes_means(setup):
    ex, perms, asked, answers, _ = setup
    answers.extend([True, True])
    assert _allow(ex, perms)["status"] == "persisted"
    first, second = asked
    assert first.startswith("⚠ THE AGENT ASKS FOR MORE RIGHTS — PERMANENTLY")
    assert "WITHOUT asking you" in first and "every later session" in first
    assert "wrote the reason below itself" in first
    assert second.startswith("⚠ CONFIRM AGAIN — THIS IS PERMANENT")


def test_a_no_to_the_second_question_writes_nothing(setup):
    ex, perms, asked, answers, tmp_path = setup
    answers.extend([True, False])
    assert _allow(ex, perms)["status"] == "denied"
    assert len(asked) == 2
    assert r"^make\b" not in perms.bash_auto_allow_patterns
    assert not (tmp_path / "user.json").exists()


def test_after_a_no_in_the_turn_it_is_not_asked(setup):
    """The pattern that was measured: refused, then asked for more."""
    ex, perms, asked, answers, _ = setup
    assert A._ask_confirm(perms, "write_file", {"path": "x"}, "write x") is False
    asked.clear()
    answers.extend([True, True])
    out = _allow(ex, perms)
    assert "refused" in out["error"] and "said no" in out["error"]
    assert asked == []


def test_the_next_turn_asks_again(setup):
    ex, perms, asked, answers, _ = setup
    A._ask_confirm(perms, "write_file", {"path": "x"}, "write x")
    perms.turn_token = A._next_turn_token()
    asked.clear()
    answers.extend([True, True])
    assert _allow(ex, perms)["status"] == "persisted"


def test_a_bundle_is_flagged_and_asked_twice(setup):
    ex, perms, asked, answers, tmp_path = setup
    project = tmp_path / "proj"
    project.mkdir()
    answers.extend([True, False])
    out = json.loads(ex._execute_remember_permission_bundle(
        {"profile": "project_dev", "directory": str(project),
         "scope": "user", "rationale": "develop"}, perms))
    assert out["status"] == "denied"
    assert asked[0].startswith("⚠ THE AGENT ASKS FOR MORE RIGHTS")
    assert f"WRITE in {project.resolve()}" in asked[0]
    assert len(asked) == 2
    assert project.resolve() not in perms.all_workspace_roots()


def test_a_deny_rule_narrows_and_is_asked_once(setup):
    ex, perms, asked, answers, _ = setup
    answers.append(True)
    out = json.loads(ex._execute_remember_permission(
        {"kind": "deny_pattern", "value": r"^curl\b", "rationale": "no net"},
        perms))
    assert out["status"] == "persisted"
    assert len(asked) == 1 and not asked[0].startswith("⚠")
