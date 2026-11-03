"""Adversarial tests for package PC1 pre_pr_check (reviewer nacht-s19).

These encode the package target from fixture repositories, not from the
builder's implementation. All are RED on main (the module does not exist
there) and are the red control for the builder's work. run() must:
- select tests from the TREE (not a fixed list),
- report blockers with file:line / test id,
- never write to the checked tree.

Fixtures are built at run time under tmp_path; the live tree is never
touched.
"""

import subprocess

import pytest

from delfin.agent import pre_pr_check


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True)


@pytest.fixture()
def base_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "--allow-empty", "-qm", "base")
    _git(repo, "checkout", "-qb", "work")
    (repo / "tests").mkdir()
    (repo / "tests" / "test_r3_direct_scheduler_call_guard.py").write_text(
        "# guard test\n"
    )
    (repo / "delfin").mkdir()
    _git(repo, "add", "-A")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "tree")
    return repo


def test_scheduler_call_change_must_run_the_scheduler_guard(base_repo):
    """Miss #6 (direct sacct call): the guard test must be selected."""
    (base_repo / "delfin" / "scheduler_thing.py").write_text("x = 1\n")
    _git(base_repo, "add", "-A")
    _git(base_repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "change")
    report = pre_pr_check.run(str(base_repo), base="main", partners=[])
    text = str(report)
    assert "test_r3_direct_scheduler_call_guard.py" in text


def test_new_skill_without_domains_is_a_blocker(base_repo):
    """Miss: skill without `domains` front matter must be reported."""
    skills = base_repo / "delfin" / "agent" / "pack" / "skills"
    skills.mkdir(parents=True)
    (skills / "new_skill.md").write_text(
        "---\ntitle: new skill\n---\nbody\n"
    )
    _git(base_repo, "add", "-A")
    _git(base_repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "skill")
    report = pre_pr_check.run(str(base_repo), base="main", partners=[])
    text = str(report)
    assert "new_skill.md" in text


def test_clean_branch_is_ready(base_repo):
    report = pre_pr_check.run(str(base_repo), base="main", partners=[])
    assert report.verdict == "READY"


def test_missing_partner_branch_is_a_blocker(base_repo):
    (base_repo / "note.txt").write_text("hi\n")
    _git(base_repo, "add", "-A")
    _git(base_repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", "note")
    report = pre_pr_check.run(str(base_repo), base="main",
                              partners=["missing-partner"])
    assert report.verdict != "READY"


def test_check_never_writes_to_the_tree(base_repo):
    before = subprocess.run(
        ["git", "status", "--porcelain"], cwd=str(base_repo),
        capture_output=True, text=True).stdout
    pre_pr_check.run(str(base_repo), base="main", partners=[])
    after = subprocess.run(
        ["git", "status", "--porcelain"], cwd=str(base_repo),
        capture_output=True, text=True).stdout
    assert before == after
