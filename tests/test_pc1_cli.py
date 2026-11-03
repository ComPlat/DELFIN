"""PC1 builder tests: the `delfin prcheck` CLI body (nacht-s38).

The real path: cli_prcheck.main() on a fixture repo, exercising exit
codes and output the way a scripted hand-over would see them.  The
argparse wiring in protected cli.py is the operator's .gate/px patch;
its red test lives in tests/test_zz_pc1_cli_wiring.py (uncommitted).
"""

import subprocess

import pytest

from delfin.agent import cli_prcheck
from delfin.agent import pre_pr_check


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True)


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "--allow-empty", "-qm", "base")
    _git(repo, "checkout", "-qb", "work")
    (repo / "tests").mkdir()
    (repo / "delfin").mkdir()
    (repo / "tests" / "test_r3_direct_scheduler_call_guard.py").write_text(
        "# guard test\n")
    _git(repo, "add", "-A")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "tree")
    monkeypatch.chdir(repo)
    return repo


def test_ready_exit_is_zero(repo, capsys):
    assert cli_prcheck.main([]) == 0
    out = capsys.readouterr().out
    assert out.startswith("READY")
    assert "tests/test_r3_direct_scheduler_call_guard.py" in out


def test_missing_partner_exit_is_one_with_blocker(repo, capsys):
    """A partner branch with its own unmerged commit is not contained."""
    _git(repo, "checkout", "-qb", "partner")
    (repo / "review_notes.md").write_text("adversarial\n")
    _git(repo, "add", "-A")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "adversarial tests")
    _git(repo, "checkout", "-q", "work")
    code = cli_prcheck.main(["--partner", "partner"])
    assert code == 1
    out = capsys.readouterr().out
    assert out.startswith("BLOCKED")
    assert "partner" in out


def test_base_flag_is_passed_through(repo, capsys):
    """--base names the branch point: with a fresh base the branch's own
    change selects nothing new, the verdict stays READY."""
    (repo / "delfin" / "new_thing.py").write_text("x = 1\n")
    _git(repo, "add", "-A")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "change")
    _git(repo, "branch", "fresh-base")
    assert cli_prcheck.main(["--base", "fresh-base"]) == 0


def test_bad_base_is_a_clean_error_not_a_traceback(repo, capsys):
    code = cli_prcheck.main(["--base", "no-such-branch"])
    out = capsys.readouterr().out
    assert code == 1
    assert out.startswith("BLOCKED")
    assert "unknown-base" in out
    assert "no-such-branch" in out
    assert "READY" not in out
