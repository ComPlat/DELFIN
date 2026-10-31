"""The repo-wide identity guards are selected for the changes they scan.

pre_pr_check.select_tests picks the tests a builder must run before
"PR ready". Two guards CI caught in wave 23 scanned tree areas their
NAMES do not suggest (test_no_account_name_is_committed.py scans the
whole tracked tree via git ls-files; test_orca_fixtures_carry_no_host_
identity.py scans tests/fixtures via rglob) and neither carries a
*guard* name nor an invariant prefix, so a builder changing delfin/ or
fixtures was told to run everything BUT the guards CI would fail on
(wave-23 real failure: CI red on guards nobody ran locally).

Control (red before the fix, measured on d48ccf6fe): both asserts
failed — the selection named safe_git/every_control_file/watchpost
etc. but neither identity guard.
"""
import sys

sys.path.insert(0, ".")

from delfin.agent import pre_pr_check  # noqa: E402

_WORKTREE = "."


def test_account_name_guard_selected_for_any_delfin_change():
    tests = pre_pr_check.select_tests(["delfin/agent/example.py"],
                                      _WORKTREE)
    assert "tests/test_no_account_name_is_committed.py" in tests, tests


def test_fixture_host_guard_selected_for_fixture_change():
    tests = pre_pr_check.select_tests(
        ["tests/fixtures/orca611/x/x.out"], _WORKTREE)
    assert "tests/test_orca_fixtures_carry_no_host_identity.py" in tests, \
        tests


def test_plain_delfin_change_also_still_names_the_invariant_families():
    """Nothing was removed by the name-filter change: the families the
    old filter selected for a delfin/ change are still there."""
    tests = pre_pr_check.select_tests(["delfin/agent/example.py"],
                                      _WORKTREE)
    assert "tests/test_every_control_file_delfin_shipped_still_runs.py" \
        in tests, tests
    assert "tests/test_delfins_own_git_reads_go_through_safe_git.py" \
        in tests, tests
