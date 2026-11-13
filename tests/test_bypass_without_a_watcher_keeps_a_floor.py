"""Bypass with no action classifier at all keeps the high-risk floor.

Input: commands the permission rules allow under ``bypassPermissions``,
through the real bash gate, with the classifier switched off
(``_classifier_settings`` None) or with no runner attached.

Output: the high-risk class (``_bypass_failure_asks``: deletes, publishing,
uploads, git commands that drop uncommitted work) is refused with nobody
to ask -- with the route that works named -- and asked with someone to
ask, exactly as when the classifier is on and gives no verdict. Before,
an absent classifier returned the rule result and the class ran on the
deny list alone. Workspace work and ref moves the reflog can undo still
run unasked; a push or pull request the user asked for in this turn (the
push grant) is the user's own verdict and is not held.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A


class Dialog:
    last_timed_out = False

    def __init__(self, answer=True):
        self.answer, self.shown = answer, []

    def ask(self, name, args, preview):
        self.shown.append(preview)
        return self.answer


@pytest.fixture(params=["off", "no-runner"])
def no_watcher(request, monkeypatch):
    if request.param == "off":
        monkeypatch.setattr(A, "_classifier_settings", lambda: None)
    else:
        monkeypatch.setattr(A, "_classifier_settings",
                            lambda: {"model": "m", "rules": [], "broad": False})
    return request.param


def _perms(tmp_path, dialog=None, *, runner=None):
    ws = tmp_path / "wt"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode="bypassPermissions")
    perms.classifier_runner = runner
    if dialog is not None:
        perms.confirm_callback = dialog.ask
    return perms


def _second(perms, cmd):
    return A._DocToolExecutor()._classifier_second_opinion(
        perms, "bash", {"command": cmd}, "bypassPermissions")


def _gate(perms, cmd):
    return A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


@pytest.mark.parametrize("cmd", [
    "rm -f results.csv",
    "git clean -n",
    "git checkout -- notes.md",
    "git restore notes.md",
    "rsync -a out/ host:/data/",
])
def test_the_class_is_refused_unattended_with_a_way_on(tmp_path, no_watcher, cmd):
    out = _second(_perms(tmp_path), cmd)
    assert out is not None, cmd
    assert "classifier" in out and "report" in out, out


def test_a_delete_refusal_names_moving_it_aside(tmp_path, no_watcher):
    out = _second(_perms(tmp_path), "rm -f results.csv")
    assert "mv" in out


@pytest.mark.parametrize("cmd", ["rm -f results.csv", "git restore notes.md"])
def test_the_class_is_asked_attended(tmp_path, no_watcher, cmd):
    dialog = Dialog(True)
    assert _second(_perms(tmp_path, dialog), cmd) is None
    assert len(dialog.shown) == 1 and cmd in dialog.shown[0]


@pytest.mark.parametrize("cmd", [
    "python3 build.py",
    "git status",
    "git commit -m 'remove the rm -rf note'",
    "git commit --amend -m 'better message'",
    "git rebase main",
    "git branch -f topic HEAD~1",
])
def test_workspace_work_and_ref_moves_still_run(tmp_path, no_watcher, cmd):
    dialog = Dialog(False)
    assert _second(_perms(tmp_path, dialog), cmd) is None, cmd
    assert _second(_perms(tmp_path), cmd) is None, cmd
    assert dialog.shown == []


def test_a_granted_push_is_the_users_own_verdict(tmp_path, no_watcher):
    perms = _perms(tmp_path)
    perms.push_grants["push"] = 1
    assert _second(perms, "git push origin topic") is None
    # Not a way past the floor for anything else in the same line.
    assert _second(perms, "git push origin topic && rm -f results.csv")


def test_an_ungranted_push_is_held(tmp_path, no_watcher):
    assert _second(_perms(tmp_path), "git push origin topic") is not None


def test_other_modes_are_unchanged(tmp_path, no_watcher):
    perms = _perms(tmp_path)
    perms.mode = "acceptEdits"
    assert A._DocToolExecutor()._classifier_second_opinion(
        perms, "bash", {"command": "rm -f results.csv"}, "acceptEdits") is None


def test_through_the_gate_unattended(tmp_path, no_watcher):
    out = _gate(_perms(tmp_path), "rm -f results.csv")
    assert out is not None and "classifier" in out
