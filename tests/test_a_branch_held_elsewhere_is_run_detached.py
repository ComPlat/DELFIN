"""A branch checked out in another worktree is run detached, not refused.

Input: a checkout (here a ``git worktree add`` worktree) asked to move to a
branch that another worktree of the same repository has checked out.
git refuses that with "fatal: 'main' is already used by worktree at ...".
Output: the checkout runs the same commit with a detached HEAD, and the
message names the branch, the commit and the checkout that holds the name.
Uncommitted tracked changes refuse the move and are left as they are.
A later update moves a detached checkout forward to the branch's current
commit (the remote-tracking one when it is ahead). The label reads
``main @ <sha> (detached, held by another checkout)``, not ``HEAD``.
"""

from __future__ import annotations

import json
import os
import subprocess

import pytest

from delfin.agent import held_branch as HB

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args) -> str:
    out = subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=True, capture_output=True, text=True)
    return out.stdout.strip()


def _commit(cwd, name, text="x\n"):
    (cwd / name).write_text(text)
    _git(cwd, "add", name)
    _git(cwd, "commit", "-q", "-m", name)
    return _git(cwd, "rev-parse", "HEAD")


@pytest.fixture
def repos(tmp_path):
    """origin (bare), main checkout on ``main``, worktree on ``feat``."""
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    main = tmp_path / "main"
    _git(tmp_path, "clone", "-q", str(origin), str(main))
    _git(main, "checkout", "-q", "-b", "main")
    _commit(main, "a.txt")
    _git(main, "push", "-q", "-u", "origin", "main")
    wt = tmp_path / "dash"
    _git(main, "worktree", "add", "-q", "-b", "feat", str(wt))
    _commit(wt, "feat.txt")
    return origin, main, wt


def _head(cwd):
    return _git(cwd, "rev-parse", "HEAD")


def test_git_itself_refuses_the_held_branch(repos):
    """The control: plain git fails, so the behaviour below is ours."""
    _o, main, wt = repos
    r = subprocess.run(["git", "switch", "main"], cwd=str(wt), env=_ENV,
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert HB.parse_held_error(r.stderr) == ("main", str(main))


def test_switching_to_a_held_branch_ends_detached_at_its_commit(repos):
    _o, main, wt = repos
    res = HB.detach_at_branch(wt, "main")
    assert res.ok, res.message
    assert _head(wt) == _head(main)
    assert _git(wt, "branch", "--show-current") == ""
    short = _git(wt, "rev-parse", "--short", "HEAD")
    assert str(main) in res.message
    assert short in res.message
    assert "detached" in res.message
    assert "git checkout --detach" in res.message
    assert HB.held_elsewhere(wt, "main") == str(main)


def test_uncommitted_changes_refuse_and_are_kept(repos):
    _o, _main, wt = repos
    before = _head(wt)
    (wt / "feat.txt").write_text("edited\n")
    res = HB.detach_at_branch(wt, "main")
    assert not res.ok
    assert "not clean" in res.message
    assert _head(wt) == before
    assert _git(wt, "branch", "--show-current") == "feat"
    assert (wt / "feat.txt").read_text() == "edited\n"


def test_the_label_names_the_branch_and_the_holder(repos):
    _o, _main, wt = repos
    assert HB.checkout_label(wt) == "feat"
    HB.detach_at_branch(wt, "main")
    short = _git(wt, "rev-parse", "--short", "HEAD")
    assert HB.checkout_label(wt) == (
        f"main @ {short} (detached, held by another checkout)")


def test_update_follows_main_forward(repos):
    _o, main, wt = repos
    HB.detach_at_branch(wt, "main")
    new = _commit(main, "b.txt")          # local main moves on
    res = HB.update_detached(wt)
    assert res is not None and res.ok, res and res.message
    assert _head(wt) == new
    short = _git(wt, "rev-parse", "--short", "HEAD")
    assert HB.checkout_label(wt).startswith(f"main @ {short} (detached")


def test_update_takes_the_remote_state_when_it_is_ahead(repos, tmp_path):
    origin, _main, wt = repos
    HB.detach_at_branch(wt, "main")
    other = tmp_path / "other"
    _git(tmp_path, "clone", "-q", str(origin), str(other))
    remote_new = _commit(other, "c.txt")
    _git(other, "push", "-q", "origin", "main")
    res = HB.update_detached(wt)
    assert res is not None and res.ok, res and res.message
    assert _head(wt) == remote_new
    assert HB.checkout_label(wt).startswith("main @ ")


def test_update_refuses_a_checkout_that_is_not_behind(repos):
    """Commits made on the detached HEAD are not moved away from."""
    _o, _main, wt = repos
    HB.detach_at_branch(wt, "main")
    mine = _commit(wt, "local.txt")
    res = HB.update_detached(wt)
    assert res is not None and not res.ok
    assert _head(wt) == mine


def test_update_is_not_for_a_checkout_on_a_branch(repos):
    _o, _main, wt = repos
    assert HB.update_detached(wt) is None


# -- the dashboard ---------------------------------------------------------

def test_the_dashboard_status_label_shows_the_detached_branch(repos):
    from delfin.dashboard import _format_git_status_html
    _o, _main, wt = repos
    HB.detach_at_branch(wt, "main")
    short = _git(wt, "rev-parse", "--short", "HEAD")
    html = _format_git_status_html(wt)
    assert f"main @ {short} (detached, held by another checkout)" in html
    assert "(detached) @" not in html


def test_the_dashboard_write_lock_works_in_a_worktree(repos):
    """``.git`` is a file in a worktree; the lock lives in its git dir."""
    from delfin.dashboard import _dashboard_git_write_lock
    _o, _main, wt = repos
    with _dashboard_git_write_lock(wt, "probe"):
        pass


# -- the agent's shell -----------------------------------------------------

def test_the_agents_checkout_of_a_held_branch_ends_detached(repos):
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
    _o, main, wt = repos
    perms = KitToolPermissions(workspace=wt)
    out = json.loads(_DocToolExecutor().execute(
        "bash", {"command": "git checkout main", "description": "t"}, perms))
    assert out.get("exit_code") == 0, out
    assert _head(wt) == _head(main)
    assert str(main) in out.get("note", "")
    assert "detached" in out.get("note", "")


def test_the_agents_checkout_with_changes_still_fails(repos):
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor
    _o, _main, wt = repos
    before = _head(wt)
    (wt / "feat.txt").write_text("edited\n")
    perms = KitToolPermissions(workspace=wt)
    out = json.loads(_DocToolExecutor().execute(
        "bash", {"command": "git checkout main", "description": "t"}, perms))
    assert out.get("exit_code") != 0, out
    assert _head(wt) == before
    assert (wt / "feat.txt").read_text() == "edited\n"
    assert "not clean" in out.get("hint", "")


def test_only_the_plain_two_word_move_is_retried():
    assert HB.plain_branch_move("git checkout main") == "main"
    assert HB.plain_branch_move("git switch main") == "main"
    for cmd in ("git checkout main && rm -rf x", "git checkout -f main",
                "git checkout main -- f", "git -C /x checkout main",
                "git switch -c main", "echo git checkout main"):
        assert HB.plain_branch_move(cmd) is None, cmd
