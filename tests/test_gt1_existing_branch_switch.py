"""GT1 phase 1: switching to an EXISTING local branch of the own worktree.

Two red controls for the wave-21 friction findings:

(A) routine shape missing: plain ``git checkout <branch>`` /
``git switch <branch>`` (no -b/-c) is judged by the model today, so the
QS sessions needed an approval for every branch switch. Red = verdict
is None instead of an allow.

(B) write gate misjudges the branch name as a path: api_client folds
patch_targets.targets(cmd) into the write targets (api_client.py
~24536-24540); patch_targets._is_path_shaped (:575-578) returns True
for any operand containing "/", so the branch name ``agent/s3-ds1b``
became the refused write target ``<worktree>/agent/s3-ds1b``.
Red = targets() names the branch operand as a path.
"""

from delfin.agent import action_classifier as ac
from delfin.agent import patch_targets as pt


def _verdict(action):
    return ac.screen(action, {"task": "GT1: existing-branch switch"})


# --- (A) routine shape --------------------------------------------------------

def test_checkout_existing_branch_is_routine(monkeypatch):
    monkeypatch.chdir(_tmp_branch_repo())
    v = _verdict({"tool": "bash", "command": "git checkout agent/s3-ds1b"})
    assert v is not None and v.decision == "allow"


def test_switch_existing_branch_is_routine(monkeypatch):
    monkeypatch.chdir(_tmp_branch_repo())
    v = _verdict({"tool": "bash", "command": "git switch agent/s3-ds1b"})
    assert v is not None and v.decision == "allow"


def test_switch_to_existing_local_branch_is_routine(monkeypatch):
    monkeypatch.chdir(_tmp_branch_repo())
    v = _verdict({"tool": "bash", "command": "git switch main"})
    assert v is not None and v.decision == "allow"


# --- (A) look-alikes stay judged ----------------------------------------------

def test_checkout_b_path_like_name_not_routine():
    # Creation with a path-shaped name is the EXISTING creation shape's
    # business, not the new switch-to-existing shape: any flag -> judged.
    v = _verdict({"tool": "bash",
                  "command": "git checkout -b agent/s3-ds1b"})
    assert v is not None and v.decision == "allow"   # existing creation shape


def test_checkout_orphan_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git checkout --orphan fresh"})
    assert v is None


def test_checkout_detach_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git checkout --detach origin/main"})
    assert v is None


def test_checkout_with_path_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git checkout main -- delfin/api.py"})
    assert v is None


def test_checkout_dashdash_not_routine():
    v = _verdict({"tool": "bash", "command": "git checkout -- file.py"})
    assert v is None


def test_checkout_f_not_routine():
    v = _verdict({"tool": "bash", "command": "git checkout -f main"})
    assert v is None


def test_checkout_b_upper_force_not_routine():
    v = _verdict({"tool": "bash", "command": "git checkout -B main"})
    assert v is None


def test_switch_c_still_creation_shape():
    v = _verdict({"tool": "bash", "command": "git switch -c new-branch"})
    assert v is not None and v.decision == "allow"


def test_switch_discard_changes_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git switch --discard-changes main"})
    assert v is None


def test_checkout_two_operands_not_routine():
    # `git checkout <tree-ish> <path>`: the second operand is a path.
    v = _verdict({"tool": "bash",
                  "command": "git checkout main delfin/api.py"})
    assert v is None


def test_checkout_shell_operator_not_routine():
    v = _verdict({"tool": "bash",
                  "command": "git checkout main && git push"})
    assert v is None


# --- (B) write-gate target parser --------------------------------------------

def _tmp_branch_repo():
    """A real tmp git repository that carries the local branches
    ``main`` and ``agent/s3-ds1b``. The fail-closed branch probe needs a
    ref store to ask; a bare directory would keep the operand a path
    forever. The routine-shape tests chdir here, because
    routine_actions asks git in os.getcwd() -- without the chdir they
    depend on whatever repo pytest starts in (the CI checkout of a PR
    head has neither branch, wave-22 CI shard 3)."""
    import os
    import subprocess
    import tempfile
    td = tempfile.mkdtemp(prefix="gt1-branch-")
    # The CI runner has no git identity, so a bare `git commit` exits
    # 128 there; sandbox the fixture from the host config either way
    # (HOME=tmp) and set identity via env, like the reviewer's git_ws
    # fixture in tests/test_gt1_adversarial_branch_switch.py.
    home = tempfile.mkdtemp(prefix="gt1-home-")
    env = dict(os.environ)
    env["HOME"] = home
    env["GIT_CONFIG_GLOBAL"] = os.path.join(home, ".gitconfig")
    env["GIT_CONFIG_SYSTEM"] = os.devnull
    env["GIT_AUTHOR_NAME"] = "gt1 fixture"
    env["GIT_AUTHOR_EMAIL"] = "gt1-fixture@example.invalid"
    env["GIT_COMMITTER_NAME"] = "gt1 fixture"
    env["GIT_COMMITTER_EMAIL"] = "gt1-fixture@example.invalid"
    for cmd in (["git", "init", "-q", "-b", "main", td],
                ["git", "-C", td, "commit", "-q", "--allow-empty",
                 "-m", "gt1 fixture"],
                ["git", "-C", td, "branch", "agent/s3-ds1b"]):
        subprocess.run(cmd, check=True, capture_output=True, env=env)
    return td


def test_checkout_branch_operand_is_not_named_as_a_write_path():
    wd = _tmp_branch_repo()
    got = pt.targets("git checkout agent/s3-ds1b", wd)
    assert got == [], got


def test_switch_branch_operand_is_not_named_as_a_write_path():
    got = pt.targets("git switch agent/s3-ds1b", "/some/worktree")
    assert got == [], got


def test_checkout_real_path_is_still_named():
    # No widening: a genuine pathspec checkout still produces its target.
    got = pt.targets("git checkout -- delfin/api.py", "/some/worktree")
    assert got == ["/some/worktree/delfin/api.py"], got


def test_checkout_branch_that_is_also_a_file_is_named():
    # A branch name that resolves to an existing file in the tree stays
    # a path target -- the file interpretation wins, fail-safe direction.
    import os
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        os.makedirs(os.path.join(td, "main"))
        got = pt.targets("git checkout main", td)
        assert got == [os.path.join(td, "main")], got


def test_checkout_paths_with_two_operands_still_named():
    got = pt.targets("git checkout main delfin/api.py", "/some/worktree")
    assert got == ["/some/worktree/delfin/api.py"], got
