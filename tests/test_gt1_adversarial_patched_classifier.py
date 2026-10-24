"""GT1 reviewer probe s19: verify the builder's patched patch_targets
copy against my OWN adversarial battery (independent of the builder's
probe). Loads .gate/px/s19_verify/new_patch_targets.py (my own fresh
origin/main copy + the builder's hunk, applied by me) via importlib --
no sys.modules pollution. Uncommitted probe (test_zz_)."""

import importlib.util
import subprocess
from pathlib import Path

import pytest

_COPY = Path(__file__).resolve().parents[1] / ".gate/px/s19_verify" / \
    "new_patch_targets.py"


def _patched():
    spec = importlib.util.spec_from_file_location("gt1_patched_pt", _COPY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def repo(tmp_path_factory):
    """A real git repo with a slashy local branch and an existing file
    of a branch-like name."""
    ws = tmp_path_factory.mktemp("gt1probe")
    env_git = ["git", "-C", str(ws)]
    subprocess.run(env_git + ["init", "-q"], check=True)
    (ws / "delfin").mkdir()
    (ws / "delfin" / "api.py").write_text("x = 1\n")
    subprocess.run(env_git + ["add", "."], check=True)
    subprocess.run(env_git + ["commit", "-q", "-m", "init"], check=True,
                   env={**__import__("os").environ,
                        "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
                        "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"})
    subprocess.run(env_git + ["branch", "agent/s3-ds1b"], check=True)
    subprocess.run(env_git + ["branch", "feat.dot.name"], check=True)
    (ws / "main").mkdir()          # a DIRECTORY named like a branch
    (ws / "main" / "f.txt").write_text("f\n")
    return ws


def _targets(cmd_path, base):
    return _patched().targets(cmd_path, str(base))


def test_branch_switch_slashy_names_no_target(repo):
    assert _targets("git checkout agent/s3-ds1b", repo) == []
    assert _targets("git switch agent/s3-ds1b", repo) == []


def test_branch_name_with_dot_is_confirmed(repo):
    assert _targets("git checkout feat.dot.name", repo) == []


def test_pathspec_after_dashes_still_named(repo):
    assert _targets("git checkout -- delfin/api.py", repo) == \
        [str(repo / "delfin" / "api.py")]


def test_branch_that_is_also_a_file_wins_the_file_reading(repo):
    # 'main' exists as a directory under base AND is a real branch:
    # the file reading must win -> target named.
    named = _targets("git checkout main", repo)
    assert [str(repo / "main")] == named


def test_two_operands_keep_the_path_reading(repo):
    # 'main' exists under base -> the conservative pre-existing reading
    # names ALL operands (main's behaviour at :506-509, unchanged by
    # GT1: with >=2 operands the module over-refuses by design).
    assert _targets("git checkout main delfin/api.py", repo) == \
        [str(repo / "main"), str(repo / "delfin" / "api.py")]


def test_no_such_branch_stays_a_path(repo):
    assert _targets("git checkout no/such/branch", repo) == \
        [str(repo / "no/such/branch")]


def test_detached_and_force_shapes_stay_untouched(repo):
    assert _targets("git checkout -f agent/s3-ds1b", repo) == []
    assert _targets("git checkout --detach agent/s3-ds1b", repo) == []
    assert _targets("git checkout -B newbr agent/s3-ds1b", repo) == []


def test_non_path_shaped_single_operand_was_never_a_target(repo):
    # 'deadbeef42': no '/', ':', glob chars -> _is_path_shaped False,
    # no file exists -> pre-existing main already returned [] for it
    # (git reads a bare operand as a branch). Not GT1's concern.
    assert _targets("git checkout deadbeef42", repo) == []


def test_quoted_spacing_shape_is_the_same_single_operand(repo):
    # shell quoting disappears before this parser; the operand itself
    # is unchanged -> branch confirmed
    assert _targets("git checkout 'agent/s3-ds1b'", repo) == []
