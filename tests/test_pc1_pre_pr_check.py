"""PC1 builder tests: pre_pr_check.run on fixture repos (nacht-s38).

Same-shape fixtures as the reviewer's adversarial file, plus the cases the
adversarial file does not cover: the scheduler-adjacent miss, the private
string, the setup-without-args miss, the count test, and the exact verdict
carrier the reviewer only greps for in text.
"""

import subprocess

import pytest

from delfin.agent import pre_pr_check


def _git(repo, *args):
    subprocess.run(["git", *args], cwd=repo, check=True, capture_output=True)


@pytest.fixture()
def repo(tmp_path):
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
    return repo


def _commit_all(repo, msg):
    _git(repo, "add", "-A")
    _git(repo, "-c", "user.email=t@e", "-c", "user.name=t",
         "commit", "-qm", msg)


def test_a_delfin_change_selects_the_tree_guard(repo):
    """The miss EF1 measured: a delfin/ change with a guard test in the tree
    that scans delfin/ selects that guard."""
    (repo / "delfin" / "scheduler_thing.py").write_text("x = 1\n")
    _commit_all(repo, "change")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert "tests/test_r3_direct_scheduler_call_guard.py" in (
        report.selected_tests)
    assert report.verdict == "READY"


def test_a_skill_without_domains_is_a_blocker(repo):
    skills = repo / "delfin" / "agent" / "pack" / "skills"
    skills.mkdir(parents=True)
    (skills / "new_skill.md").write_text(
        "---\ntitle: new skill\n---\nbody\n")
    _commit_all(repo, "skill")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert report.verdict == "BLOCKED"
    assert any("new_skill.md" in str(b) for b in report.blockers)


def test_a_skill_with_domains_is_ready(repo):
    skills = repo / "delfin" / "agent" / "pack" / "skills"
    skills.mkdir(parents=True)
    (skills / "new_skill.md").write_text(
        "---\nname: new\ndomains: [code]\n---\nbody\n")
    _commit_all(repo, "skill")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert report.verdict == "READY"


def test_a_private_string_in_the_diff_blocks(repo):
    (repo / "delfin" / "note.py").write_text(
        'api_key = "ZXkj8m2Pqrs9wvut0123456789abcdef"\n')
    _commit_all(repo, "token")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert report.verdict == "BLOCKED"
    assert any("private-string" in str(b) for b in report.blockers)


def test_a_merges_partner_branch_is_ready(repo):
    (repo / "note.txt").write_text("hi\n")
    _commit_all(repo, "note")
    _git(repo, "checkout", "-qb", "partner")
    (repo / "review.py").write_text("# review\n")
    _commit_all(repo, "review")
    _git(repo, "checkout", "work")
    _git(repo, "merge", "--no-edit", "partner")
    report = pre_pr_check.run(str(repo), base="main", partners=["partner"])
    assert report.verdict == "READY"


def test_check_writes_nothing(repo):
    before = subprocess.run(
        ["git", "status", "--porcelain"], cwd=str(repo),
        capture_output=True, text=True).stdout
    (repo / "delfin" / "x.py").write_text("y = 2\n")
    _commit_all(repo, "x")
    pre_pr_check.run(str(repo), base="main", partners=[])
    after = subprocess.run(
        ["git", "status", "--porcelain"], cwd=str(repo),
        capture_output=True, text=True).stdout
    assert before == after


def test_a_direct_scheduler_call_is_a_blocker(repo):
    (repo / "delfin" / "scheduler_thing.py").write_text(
        'subprocess.run(["sacct", "-X"])\n')
    _commit_all(repo, "sacct")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert report.verdict == "BLOCKED"
    assert any("direct-scheduler-call" in str(b) for b in report.blockers)


def test_a_clean_benchmark_setup_with_args_is_not_flagged(repo):
    setup = repo / "delfin" / "agent" / "pack" / "benchmark" / "setup"
    setup.mkdir(parents=True)
    (setup / "run.py").write_text(
        "import sys\nif len(sys.argv) < 2:\n"
        "    print('usage'); sys.exit(2)\nprint('ok')\n")
    _commit_all(repo, "setup")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert not any("setup" in b.kind for b in report.blockers)


def test_changed_files_are_listed(repo):
    (repo / "delfin" / "a.py").write_text("a = 1\n")
    (repo / "tests" / "test_a_thing.py").write_text("def t(): pass\n")
    _commit_all(repo, "files")
    report = pre_pr_check.run(str(repo), base="main", partners=[])
    assert "delfin/a.py" in report.changed_files
    assert "tests/test_a_thing.py" in report.changed_files


def test_guard_that_only_mentions_the_area_does_not_over_select(repo):
    """Reviewer finding (s19, on 73e6aaa6): a guard whose source only
    MENTIONS delfin/ in a comment must not be selected for every delfin/
    change — only a guard that scans the area does.  The guard file sits
    at the branch point (unchanged here), so only its scanned surface can
    select it, not its own being changed."""
    (repo / "tests" / "test_output_guard.py").write_text(
        '# a unit test file, only a comment mentions delfin/\n'
        'def t(): pass\n')
    (repo / "tests" / "test_r3_direct_scheduler_call_guard.py").write_text(
        # Fixture follows the REAL r3 guard's source shape
        # (tests/test_r3_direct_scheduler_call_guard.py:26 'import ast'): a
        # guard is selected because its SOURCE scans the area, not because
        # a comment says so — the forgeable '# guard test' marker branch
        # was removed in cdc2e668c (reviewer finding A04). The comment-only
        # file above stays the red case: it is NOT selected.
        'import pathlib\n'
        'for p in pathlib.Path("delfin").rglob("*.py"):\n'
        '    pass\n')
    _commit_all(repo, "guards")
    _git(repo, "branch", "guards")
    (repo / "delfin" / "agent").mkdir()
    (repo / "delfin" / "agent" / "new.py").write_text("x = 1\n")
    _commit_all(repo, "delfin change")
    report = pre_pr_check.run(str(repo), base="guards", partners=[])
    names = [t.rsplit("/", 1)[-1] for t in report.selected_tests]
    assert "test_output_guard.py" not in names
    assert "test_r3_direct_scheduler_call_guard.py" in names


def test_guard_scanning_the_area_still_selects(repo):
    """The converse of the finding: a guard that really scans delfin/ (its
    source walks the area) is selected for a delfin/ change."""
    (repo / "tests" / "test_output_guard.py").write_text(
        'import pathlib\n'
        'for p in pathlib.Path("delfin").rglob("*.py"):\n'
        '    pass\n')
    _commit_all(repo, "scanner")
    _git(repo, "branch", "scanner")
    (repo / "delfin" / "agent").mkdir()
    (repo / "delfin" / "agent" / "new.py").write_text("x = 1\n")
    _commit_all(repo, "delfin change")
    report = pre_pr_check.run(str(repo), base="scanner", partners=[])
    assert "tests/test_output_guard.py" in report.selected_tests
