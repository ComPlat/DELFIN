""""My change broke this test" was a hypothesis the agent could not test.

The integrity addendum tells the agent to run the experiment that could
REFUTE its hypothesis. For the commonest hypothesis in software work --
a check is red, and the change I just made is why -- the refuting
experiment is a control: run the SAME check against the state before the
change. The agent had no way to do it. ``enter_worktree`` branched from
HEAD and only from HEAD, so every worktree it could make already
contained the change under suspicion.

Getting it wrong costs the same in both directions. Assume the failure is
yours and you revert work that was right; assume it is not and you ship
the regression. I have been on both sides of that this week -- twice
tonight, once with two interpreter-guard holes that turned out to predate
the branch, and once with a test that fails only while a live benchmark
is running -- and each time the answer came from running the check at the
baseline, which took one command.

So ``base_ref`` starts the worktree at a named commit of the SAME
repository. Nothing else changes: a fresh branch, the repository the
workspace already is, the same teardown. A ref is not a way to reach
another repository, which is what the containment added in September was
about.
"""

from __future__ import annotations

import json
import subprocess

import pytest

import delfin.agent.api_client as A
from delfin.agent import worktree as WT


@pytest.fixture
def repo(tmp_path):
    """A repo whose test was already broken before the agent's commit."""
    def git(*args):
        return subprocess.run(["git", *args], cwd=str(tmp_path),
                              text=True, capture_output=True, check=False)

    git("init", "-q")
    git("config", "user.email", "a@b.c")
    git("config", "user.name", "t")
    (tmp_path / "calc.py").write_text(
        "def energy():\n    return 1 / 0\n", encoding="utf-8")
    git("add", "-A")
    git("commit", "-q", "-m", "broken from the start")
    base = git("rev-parse", "HEAD").stdout.strip()
    (tmp_path / "report.py").write_text("HEADER = 'R'\n", encoding="utf-8")
    git("add", "-A")
    git("commit", "-q", "-m", "the agent's own change")
    return tmp_path, base


def _perms(path):
    return A.KitToolPermissions(mode="bypassPermissions", workspace=str(path))


def _run(name, args, perms):
    return json.loads(A._doc_executor.execute(name, args, perms))


def test_a_worktree_can_start_at_the_baseline(repo):
    ws, base = repo
    perms = _perms(ws)
    out = _run("enter_worktree", {"base_ref": base}, perms)
    assert out["status"] == "ok"
    assert out["base_ref"] == base
    wt = __import__("pathlib").Path(out["path"])
    try:
        # The state BEFORE the change: calc.py is there, report.py is not.
        assert (wt / "calc.py").is_file()
        assert not (wt / "report.py").exists()
    finally:
        _run("exit_worktree", {"path": str(wt), "keep_if_changed": False},
             perms)
    assert not wt.exists()


def test_without_it_the_worktree_still_starts_at_head(repo):
    """The default is unchanged -- this is an addition, not a change."""
    ws, base = repo
    perms = _perms(ws)
    out = _run("enter_worktree", {}, perms)
    wt = __import__("pathlib").Path(out["path"])
    try:
        assert out["base_ref"] != base
        assert (wt / "report.py").is_file()
    finally:
        _run("exit_worktree", {"path": str(wt), "keep_if_changed": False},
             perms)


def test_the_control_answers_the_question(repo):
    """End to end, the way a model would run it: the check fails here, and
    it fails at the baseline too, so the change is not the cause."""
    ws, base = repo
    (ws / "test_calc.py").write_text(
        "from calc import energy\n\n"
        "def test_energy():\n    assert energy() == 1\n", encoding="utf-8")
    subprocess.run(["git", "add", "-A"], cwd=str(ws), capture_output=True)
    subprocess.run(["git", "commit", "-q", "-m", "test"], cwd=str(ws),
                   capture_output=True)
    perms = _perms(ws)
    check = {"command": "python3 -m pytest -q test_calc.py 2>&1 | tail -2",
             "description": "run the failing test"}

    here = _run("bash", check, perms)
    assert "ZeroDivisionError" in here["stdout"]

    wt = _run("enter_worktree", {"base_ref": base}, perms)["path"]
    try:
        # The baseline has no test file at all -- copy the check in, which
        # is what an agent does when the test itself is new.
        (__import__("pathlib").Path(wt) / "test_calc.py").write_text(
            (ws / "test_calc.py").read_text(encoding="utf-8"),
            encoding="utf-8")
        there = _run("bash", {**check, "cwd": wt}, perms)
        assert "ZeroDivisionError" in there["stdout"]
    finally:
        _run("exit_worktree", {"path": wt, "keep_if_changed": False}, perms)


# ---------------------------------------------------------------------------
# A ref is a ref
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("ref", [
    "--force",
    "-b",
    "; rm -rf /",
    "$(whoami)",
    "`id`",
    "a" * 400,
    "../../../etc/passwd",
])
def test_a_ref_that_is_not_a_ref_is_refused(repo, ref):
    ws, _ = repo
    out = _run("enter_worktree", {"base_ref": ref}, _perms(ws))
    assert "error" in out, ref
    assert "path" not in out, ref


def test_an_unknown_ref_says_so_rather_than_failing_obscurely(repo):
    ws, _ = repo
    out = _run("enter_worktree", {"base_ref": "no-such-ref"}, _perms(ws))
    assert "does not name a commit" in out["error"]


def test_a_ref_cannot_reach_another_repository(repo, tmp_path_factory):
    """The containment the September fix added: repo_dir is the workspace's
    own repository, and base_ref names a commit inside it. Naming a commit
    that exists only elsewhere is an unknown ref, not a doorway."""
    other = tmp_path_factory.mktemp("other")
    subprocess.run(["git", "init", "-q"], cwd=str(other), capture_output=True)
    for k, v in (("user.email", "a@b.c"), ("user.name", "t")):
        subprocess.run(["git", "config", k, v], cwd=str(other),
                       capture_output=True)
    (other / "secret.txt").write_text("s\n", encoding="utf-8")
    subprocess.run(["git", "add", "-A"], cwd=str(other), capture_output=True)
    subprocess.run(["git", "commit", "-q", "-m", "x"], cwd=str(other),
                   capture_output=True)
    foreign = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=str(other), text=True,
        capture_output=True).stdout.strip()

    ws, _ = repo
    out = _run("enter_worktree", {"base_ref": foreign}, _perms(ws))
    assert "error" in out
    assert "does not name a commit" in out["error"]


def test_the_primitive_refuses_the_same_things(tmp_path):
    subprocess.run(["git", "init", "-q"], cwd=str(tmp_path),
                   capture_output=True)
    with pytest.raises(WT.WorktreeError):
        WT.enter_worktree(tmp_path, base_ref="--force")
    with pytest.raises(WT.WorktreeError):
        WT.enter_worktree(tmp_path, base_ref="nope")


# ---------------------------------------------------------------------------
# ...and the rule that tells the agent to reach for it
# ---------------------------------------------------------------------------

def test_the_addendum_states_it_and_it_reaches_the_prompt():
    """A rule is a rule when it is in the BUILT prompt."""
    import re

    from delfin.agent.prompt_loader import PromptLoader

    flat = re.sub(r"\s+", " ", PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo"))
    assert "A failure needs a control before you attribute it" in flat
    assert "base_ref" in flat
