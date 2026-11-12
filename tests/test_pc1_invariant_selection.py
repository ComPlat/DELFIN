"""Invariant tests (test_every_*/test_delfins_own_*) are selectable.

RED CONTROL, written before the fix (operator finding on 4723e1b0):
the operator's build failed CI on tests/test_delfins_own_git_reads_
go_through_safe_git.py — an invariant test whose name carries neither
"guard" nor a changed file's stem, so pre_pr_check's selection never
picked it. The fix: _guard_tests_for also scans the source of
test_every_*.py / test_delfins_own_*.py files for what changed, the
same scan-vs-mention rule guards use.
"""

import ast
from pathlib import Path

import delfin.agent.pre_pr_check as ppc


def _select_for_the_safe_git_violation() -> list[str]:
    # The shape of the operator's finding: a delfin/ module calls git
    # directly. The invariant test whose source AST-scans delfin/ for
    # subprocess-run("git"...) sites must be selected for it.
    return ppc.select_tests(["delfin/agent/pre_pr_check.py"], worktree=".")


def test_safe_git_invariant_is_selected_for_a_delfin_change():
    tests = _select_for_the_safe_git_violation()
    assert ("tests/test_delfins_own_git_reads_go_through_safe_git.py"
            in tests), (
        "an invariant test that scans delfin/ must be selected for a "
        f"delfin change; selected: {tests}")


def test_mention_only_invariant_is_not_selected():
    # A test_every_* file that only MENTIONS the changed file's area in a
    # comment must not be selected (the scan-vs-mention rule).
    tests = _select_for_the_safe_git_violation()
    assert "tests/test_every_file_has_a_licence_by_reuse.py" not in tests, (
        "a mention-only invariant must not be selected")


def test_selection_rule_is_stated_in_the_docstring():
    doc = (Path(ppc.__file__).read_text(encoding="utf-8"))
    assert "test_every_" in doc and "test_delfins_own_" in doc, (
        "the invariant-prefixed families must be part of the documented "
        "selection rules")
