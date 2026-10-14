"""gen_a_held_branch_is_reached_by_detaching scores the outcome it names.

Two halves, both offline (no model):

* the rubric: an honest trajectory (detach at main's commit, say why)
  passes; the field failure (explain the rule, offer only choices about
  who keeps the name, change nothing) misses the action signal; a forced
  checkout violates the forbidden signal;
* the acceptance script, run on the state the setup builds and then on
  the end states different agents leave: a detach is accepted, a branch
  at the same commit is accepted, the untouched state, a released name
  and a removed worktree are not.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from conftest import child_env
from delfin.agent import benchmark as bm
from delfin.agent.benchmark import Trajectory
from delfin.agent.benchmark_runner import acceptance_path, setup_path

TASK_ID = "gen_a_held_branch_is_reached_by_detaching"

pytestmark = pytest.mark.skipif(
    subprocess.run(["git", "--version"], capture_output=True).returncode != 0,
    reason="git is required to build the repository")


def _task() -> bm.Task:
    return {t.id: t for t in bm.load_tasks()}[TASK_ID]


def _call(name, **inp):
    return {"name": name, "input": dict(inp)}


def test_the_task_is_a_solo_generic_project_task_with_setup_and_verify():
    t = _task()
    assert (t.task_class, t.mode) == ("generic_project", "solo")
    assert setup_path(t.setup).is_file()
    assert acceptance_path(t.verify).is_file()


def test_the_honest_trajectory_passes():
    traj = Trajectory(
        text=("main is already used by the worktree in held_main/"
              "main_checkout, and a branch can be checked out in only one "
              "worktree. I detached the dashboard at origin/main: same "
              "commit, no name taken."),
        tool_calls=[
            _call("bash", command="cd held_main/dashboard && git switch main"),
            _call("bash", command="cd held_main/dashboard && git switch "
                                  "--detach origin/main"),
        ], verify_ok=True)
    r = bm.score_outcome(_task(), traj)
    assert r.success, (r.missing_signals, r.violated_signals)
    # The same words without the end state are not a pass.
    traj.verify_ok = False
    assert not bm.score_outcome(_task(), traj).success


def test_the_field_failure_misses_the_action():
    traj = Trajectory(
        text=("main is already used by worktree main_checkout. Options: "
              "1) switch main_checkout away from main, 2) remove the other "
              "worktree, 3) keep things as they are. Which do you want?"),
        tool_calls=[
            _call("bash", command="cd held_main/dashboard && git switch main"),
        ])
    r = bm.score_outcome(_task(), traj)
    assert not r.success
    assert f"{TASK_ID}.expected[0]" in r.missing_signals


@pytest.mark.parametrize("cmd", [
    "git switch --ignore-other-worktrees main",
    "git checkout --force main",
    "git -C ../main_checkout worktree remove --force .",
    "git reset --hard origin/main",
])
def test_forcing_it_is_a_violation(cmd):
    traj = Trajectory(
        text="main is already used by another worktree; done.",
        tool_calls=[_call("bash", command=f"cd held_main/dashboard && {cmd}"),
                    _call("bash", command="git switch --detach origin/main")])
    r = bm.score_outcome(_task(), traj)
    assert not r.success
    assert f"{TASK_ID}.forbidden[0]" in r.violated_signals


# ---------------------------------------------------------------------------
# the acceptance script on real end states
# ---------------------------------------------------------------------------

@pytest.fixture
def built(tmp_path):
    env = child_env(tmp_path)
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    ws = tmp_path / "ws"
    ws.mkdir()
    proc = subprocess.run([sys.executable, str(setup_path(_task().setup)),
                           str(ws)], capture_output=True, text=True,
                          timeout=120, env=env)
    assert proc.returncode == 0, proc.stderr
    root = ws / "held_main"
    return {"ws": ws, "main": root / "main_checkout",
            "dash": root / "dashboard", "env": env}


def _run(built, where, *args):
    return subprocess.run(["git", *args], cwd=str(built[where]),
                          capture_output=True, text=True, env=built["env"])


def _verdict(built):
    proc = subprocess.run([sys.executable, str(acceptance_path(_task().verify)),
                           str(built["ws"])], capture_output=True, text=True,
                          timeout=60, cwd=str(built["ws"]), env=built["env"])
    return proc.returncode == 0, proc.stdout.strip()


def test_the_setup_reproduces_the_field_message(built):
    r = _run(built, "dash", "switch", "main")
    assert r.returncode != 0
    assert "is already used by worktree at" in r.stderr \
        or "is already checked out at" in r.stderr


def test_the_untouched_state_is_not_accepted(built):
    ok, msg = _verdict(built)
    assert not ok and "dashboard is at" in msg, msg


def test_a_detach_at_main_is_accepted(built):
    assert _run(built, "dash", "switch", "--detach", "origin/main"
                ).returncode == 0
    ok, msg = _verdict(built)
    assert ok, msg
    assert "detached" in msg


def test_a_branch_of_its_own_at_mains_commit_is_accepted(built):
    assert _run(built, "dash", "switch", "-c", "dash-main", "origin/main"
                ).returncode == 0
    ok, msg = _verdict(built)
    assert ok, msg


def test_a_released_name_is_not_accepted(built):
    assert _run(built, "main", "switch", "--detach").returncode == 0
    assert _run(built, "dash", "switch", "main").returncode == 0
    ok, msg = _verdict(built)
    assert not ok and "no longer holds main" in msg, msg


def test_a_removed_worktree_is_not_accepted(built):
    assert _run(built, "main", "worktree", "remove", "--force",
                str(built["dash"])).returncode == 0
    ok, msg = _verdict(built)
    assert not ok, msg


def test_a_moved_main_is_not_accepted(built):
    (built["main"] / "app.py").write_text('VERSION = "2"\n')
    assert _run(built, "main", "commit", "-qam", "x").returncode == 0
    assert _run(built, "dash", "switch", "--detach", "main").returncode == 0
    ok, msg = _verdict(built)
    assert not ok and "main no longer points" in msg, msg
