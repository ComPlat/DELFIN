"""GA1: the readiness report asks whether the change's tests are honest.

Input:  a throwaway git repository with one commit beyond its base, driven
        through ``pre_pr_check.run`` — the function ``delfin prcheck`` and
        the swarm RPC both call.
Output: whether the honesty findings appear, whether they change the
        verdict, and whether the hunk ranges the scanner needs are right.
Semantics: ``test_honesty.check`` is a 320-line AST scanner built from
        wave 17's measurements. Its only caller in ``delfin/`` was
        ``verify_guard.scan_for_dishonest_tests``, and the only callers of
        THAT were tests — so the capability was unreachable from the
        product while the reachability baseline counted the module as
        referenced, because the wrapper named it. Nothing here touches the
        real repository: every path is inside ``tmp_path``.

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_ga1_the_readiness_report_asks_about_the_tests.py -q
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.agent import pre_pr_check

HONEST = '''"""A test that measures what it claims."""


def test_two_and_two():
    assert 2 + 2 == 4
'''

DISHONEST = '''"""A test of the kind wave 17 found committed."""


def test_a_feature_behind_a_switch():
    try:
        from delfin.agent import pre_pr_check
        assert pre_pr_check is not None
    except Exception:
        pass
'''


def _git(tree: Path, *args: str) -> str:
    out = subprocess.run(["git", *args], cwd=tree, capture_output=True,
                          text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    return out.stdout


@pytest.fixture
def repo(tmp_path):
    """A repository with a base commit and nothing else yet."""
    tree = tmp_path / "repo"
    (tree / "tests").mkdir(parents=True)
    _git(tmp_path, "init", "-q", "repo")
    _git(tree, "config", "user.email", "t@example.invalid")
    _git(tree, "config", "user.name", "t")
    (tree / "README.md").write_text("base\n")
    _git(tree, "add", "README.md")
    _git(tree, "commit", "-q", "-m", "base")
    base = _git(tree, "rev-parse", "HEAD").strip()
    return tree, base


def _commit(tree: Path, rel: str, body: str) -> None:
    (tree / rel).write_text(body)
    _git(tree, "add", rel)
    _git(tree, "commit", "-q", "-m", f"add {rel}")


def test_a_dishonest_test_is_named_in_the_report(repo):
    tree, base = repo
    _commit(tree, "tests/test_zz_probe.py", DISHONEST)
    report = pre_pr_check.run(str(tree), base)
    note = "\n".join(report.notes)
    assert "test_honesty" in note, report.notes
    assert "test_zz" in note, "the committed-probe finding"
    assert "except" in note, "the bare-swallow finding"


def test_the_finding_is_a_note_and_never_the_verdict(repo):
    """Its own report text says it blocks nothing. A finding among the
    blockers would turn READY into BLOCKED and stop work the author is
    entitled to do."""
    tree, base = repo
    _commit(tree, "tests/test_zz_probe.py", DISHONEST)
    report = pre_pr_check.run(str(tree), base)
    assert report.notes, "the finding was produced"
    kinds = {b.kind for b in report.blockers}
    assert "test-honesty" not in kinds, kinds
    assert not any("honesty" in b.detail.lower() for b in report.blockers)


def test_an_honest_test_produces_no_note(repo):
    """The other direction: a report that always says something says
    nothing. Without this the suite would pass for a scanner wired to
    emit unconditionally."""
    tree, base = repo
    _commit(tree, "tests/test_plain.py", HONEST)
    report = pre_pr_check.run(str(tree), base)
    assert report.notes == [], report.notes


def test_the_records_carry_one_hunk_range_each(repo):
    """The only scanner that reads the range asks whether NEW code sits
    inside a bare except, so a file-wide range would answer that for lines
    the change never touched."""
    tree, base = repo
    (tree / "a.py").write_text("x = 1\n" * 40)
    _git(tree, "add", "a.py")
    _git(tree, "commit", "-q", "-m", "forty lines")
    # Two separate edits, far apart: two hunks, not one span.
    lines = (tree / "a.py").read_text().splitlines()
    lines[2] = "x = 99"
    lines[35] = "x = 98"
    (tree / "a.py").write_text("\n".join(lines) + "\n")
    _git(tree, "add", "a.py")
    _git(tree, "commit", "-q", "-m", "two edits")

    records = pre_pr_check.change_records(str(tree), "HEAD~1")
    a_py = [r for r in records if r[0] == "a.py"]
    assert len(a_py) == 2, a_py
    starts = sorted(r[1] for r in a_py)
    assert starts == [3, 36], a_py
    for _path, start, end, source in a_py:
        assert start == end, "a one-line edit is a one-line range"
        assert source and "x = 99" in source, "the NEW text, not the old"


def test_a_deleted_file_contributes_a_quiet_record(repo):
    """A path the worktree no longer has must not raise and must not make
    the text scanners read a stale file."""
    tree, base = repo
    # Deleted relative to the BASE, so the diff names it. A file added and
    # removed inside the range is no net change and appears in no diff --
    # the first version of this test asserted on that and measured nothing.
    _git(tree, "rm", "-q", "README.md")
    _git(tree, "commit", "-q", "-m", "remove it")
    records = pre_pr_check.change_records(str(tree), base)
    gone = [r for r in records if r[0] == "README.md"]
    assert gone and gone[0][3] is None, gone


def test_the_product_reaches_the_scanner_not_only_the_suite():
    """The reachability check, which is what was actually wrong: the
    capability existed, was tested, and nothing in delfin/ called it."""
    root = Path(pre_pr_check.__file__).resolve().parent
    callers = [p.name for p in root.rglob("*.py")
               if p.name != "verify_guard.py"
               and "scan_for_dishonest_tests" in p.read_text()]
    assert "pre_pr_check.py" in callers, callers
    baseline = (root.parent.parent / "tests" / "fixtures"
                / "unmounted_capabilities.txt")
    if baseline.is_file():
        rows = baseline.read_text()
        assert "scan_for_dishonest_tests" not in rows
        assert "test_honesty_feedback" not in rows
