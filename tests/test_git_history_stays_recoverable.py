"""Git history stays recoverable when nobody reads the command.

Input: shell commands through the real bash gate
(``_DocToolExecutor._run_permission_gate``) under ``bypassPermissions``,
unattended (no confirm callback) and attended, with the action classifier
off, plus the pure predicate ``reversibility.git_history_refusal``.

Output: what removes the way back -- a force push in any spelling
(``--force-with-lease``, a ``+refspec``, ``--mirror``, ``--prune``, a
short-option bundle), expiring or deleting the reflog, an immediate prune,
switching the reflog off for one command or for good, moving a tag (tags
have no reflog) -- is refused in bypass and wherever nobody can be asked,
with the route that works named; asked in an attended session of another
mode. Ref moves the reflog can undo (``branch -f``, a soft or mixed reset,
``rebase``, ``commit --amend``, ``checkout -B``, ``switch -C``,
``update-ref`` of a branch) are not touched by this rule.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import reversibility as R


@pytest.fixture(autouse=True)
def _classifier_off(monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)


class Dialog:
    last_timed_out = False

    def __init__(self, answer=True):
        self.answer, self.shown = answer, []

    def ask(self, name, args, preview):
        self.shown.append(preview)
        return self.answer


def _perms(tmp_path, mode="bypassPermissions", dialog=None):
    ws = tmp_path / "wt"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode=mode)
    if dialog is not None:
        perms.confirm_callback = dialog.ask
    return perms


def _gate(perms, cmd):
    return A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


UNRECOVERABLE = [
    "git push --force-with-lease origin topic",
    "git push --force-with-lease=topic:abc123 origin topic",
    "git push --force-if-includes --force-with-lease origin topic",
    "git push origin +topic",
    "git push origin +HEAD:refs/heads/main",
    "git push --mirror backup",
    "git push --prune origin 'refs/heads/*:refs/heads/*'",
    "git push -uf origin topic",
    "git push -d origin topic",
    "git push --forc origin topic",
    "git -C repo push --force-with-lease",
    "bash -c 'git push --force-with-lease origin topic'",
    "git reflog expire --expire=now --all",
    "git reflog expire --expire-unreachable=now --all",
    "git reflog delete HEAD@{1}",
    "git gc --prune=now",
    "git gc --aggressive --prune=all",
    "git prune",
    "git prune --expire=now",
    "git -c core.logAllRefUpdates=false commit --amend -m x",
    "git -c core.logallrefupdates=0 rebase main",
    "git -c gc.reflogExpire=now gc",
    "git -c gc.reflogExpireUnreachable=now gc",
    "git -c gc.pruneExpire=now gc",
    "git --config-env=core.logAllRefUpdates=X commit --amend",
    "git config core.logAllRefUpdates false",
    "git config --local core.logAllRefUpdates false",
    "git config set gc.reflogExpire now",
    "git config remote.origin.mirror true",
    "git config remote.origin.push '+refs/heads/*:refs/heads/*'",
    "GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=core.logAllRefUpdates "
    "GIT_CONFIG_VALUE_0=false git commit --amend -m x",
    "GIT_CONFIG_PARAMETERS=\"'core.logallrefupdates'='false'\" git rebase main",
    "git tag -f v1.0",
    "git tag -fa v1.0 -m release",
    "git update-ref refs/tags/v1.0 HEAD",
    "git update-ref --stdin < moves.txt",
]

RECOVERABLE = [
    "git branch -f topic HEAD~1",
    "git reset --soft HEAD~1",
    "git reset HEAD~1",
    "git rebase main",
    "git commit --amend -m 'better message'",
    "git checkout -B topic origin/main",
    "git switch -C topic",
    "git update-ref refs/heads/topic HEAD~1",
    "git push origin topic",
    "git push -u origin HEAD:refs/heads/new-name",
    "git push --no-force-with-lease origin topic",
    "git push -o ci.skip origin topic",
    "git gc",
    "git gc --prune",
    "git gc --prune=2.weeks.ago",
    "git prune --dry-run",
    "git prune --expire=2.weeks.ago",
    "git reflog",
    "git reflog show topic",
    "git config --get core.logAllRefUpdates",
    "git config core.logAllRefUpdates true",
    "git -c gc.reflogExpire=never gc",
    "git tag v1.1",
    "git commit -m 'never run git push --mirror or git reflog expire'",
    "grep -rn 'reflog expire' docs",
]


@pytest.mark.parametrize("cmd", UNRECOVERABLE)
def test_the_predicate_names_what_removes_the_way_back(cmd):
    assert R.git_history_refusal(cmd), cmd


@pytest.mark.parametrize("cmd", RECOVERABLE)
def test_the_predicate_leaves_recoverable_moves_alone(cmd):
    assert R.git_history_refusal(cmd) == "", cmd


@pytest.mark.parametrize("cmd", UNRECOVERABLE)
def test_bypass_unattended_refuses_with_the_alternative(tmp_path, cmd):
    out = _gate(_perms(tmp_path), cmd)
    assert out is not None, cmd
    low = out.lower()
    # The deny list answers first for `--force-if-includes` (it reads as
    # --force), naming its own route: a branch of your own.
    assert ("new branch" in low or "branch of your own" in low
            or "reflog" in low or "new tag" in low), out


def test_a_force_push_refusal_names_a_new_branch(tmp_path):
    out = _gate(_perms(tmp_path), "git push --force-with-lease origin topic")
    assert "HEAD:refs/heads/<new-name>" in out


def test_a_reflog_refusal_names_leaving_gc_to_git(tmp_path):
    out = _gate(_perms(tmp_path), "git reflog expire --expire=now --all")
    assert "git gc" in out and "reflog" in out


def test_bypass_attended_refuses_too_without_a_dialog(tmp_path):
    dialog = Dialog()
    out = _gate(_perms(tmp_path, dialog=dialog),
                "git reflog expire --expire=now --all")
    assert out is not None
    assert dialog.shown == []


def test_default_mode_attended_asks_once(tmp_path):
    # `git reflog ...` is on the auto-allow list, so before this rule the
    # expire ran unasked in default mode.
    dialog = Dialog(False)
    out = _gate(_perms(tmp_path, mode="default", dialog=dialog),
                "git reflog expire --expire=now --all")
    assert out is not None
    assert len(dialog.shown) == 1 and "reflog" in dialog.shown[0]


def test_default_mode_with_nobody_to_ask_refuses(tmp_path):
    out = _gate(_perms(tmp_path, mode="default"),
                "git reflog expire --expire=now --all")
    assert out is not None and "reflog" in out


@pytest.mark.parametrize("cmd", [
    "git branch -f topic HEAD~1",
    "git commit --amend -m 'better message'",
    "git rebase main",
    "git reset --soft HEAD~1",
])
def test_recoverable_ref_moves_run_in_bypass_unattended(tmp_path, cmd):
    assert _gate(_perms(tmp_path), cmd) is None, cmd
