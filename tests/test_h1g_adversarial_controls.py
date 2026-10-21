"""counter_test.assess: a control counts on the fix's own changed test.

A completed control on an unrelated file or on another node of the same
file must not stand for the fix's own test that was only seen green or
never ran. A whole-file red plus green counts only on a changed test
file.
"""

from delfin.agent.counter_test import assess


class TestUnrelatedTargetDoesNotMask:
    """A controlled UNRELATED file cannot mask a green-only
    fix. changed test_thing.py (its ::test_a ran green only) + an
    unrelated test_other.py::test_b that is red+green -> the fix is NOT
    shown; its own test is named as never-seen-failing."""

    def test_uncontrolled_fix_test_not_masked_by_unrelated_control(self):
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [
            {"target": "tests/test_thing.py::test_a",
             "green": True},                      # fix's test: green only
            {"target": "tests/test_other.py::test_b",
             "red": True, "green": True},         # unrelated: controlled
        ]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py::test_a"
            " was never seen failing"
        )

    def test_unrelated_control_does_not_mask_changed_never_run(self):
        # Even a fully controlled unrelated file does not cover the fix's
        # changed-but-never-run test (rule-1 still fires first).
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_other.py::test_b",
                 "red": True, "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py"
            " was changed but never run"
        )


class TestSameFileOtherNodeDoesNotMask:
    """A red+green on an OLD node in the SAME file must not
    mask the fix's NEW counter-test node that is only green (never red)
    -- the "passes with and without the fix" anti-pattern."""

    def test_green_only_new_node_not_masked_by_old_node_control(self):
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [
            {"target": "tests/test_thing.py::test_old",
             "red": True, "green": True},         # old node: controlled
            {"target": "tests/test_thing.py::test_new",
             "green": True},                      # new counter-test: green only
        ]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py::test_new"
            " was never seen failing"
        )

    def test_whole_file_control_does_not_mask_new_node_green_only(self):
        # A whole-file red+green (single aggregate) is ambiguous and the
        # green-only NEW node still owns the verdict.
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [
            {"target": "tests/test_thing.py", "red": True, "green": True},
            {"target": "tests/test_thing.py::test_new", "green": True},
        ]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py::test_new"
            " was never seen failing"
        )


class TestWholeFileProtocol:
    """A whole-file target (no ::node) is never a completed
    control -- it cannot say which node carried the red."""

    def test_whole_file_single_aggregate_not_shown(self):
        changed = ["delfin/agent/thing.py"]
        runs = [{"target": "tests/test_thing.py", "red": True, "green": True}]
        result = assess(changed, runs)
        assert result != "shown"
        assert result.startswith("no control run:")

    def test_whole_file_red_green_in_scope_is_a_control(self):
        # In-scope (test file changed): a whole-file red+green on the fix's
        # own changed counter-test IS a completed control.
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_thing.py", "red": True, "green": True}]
        result = assess(changed, runs)
        assert result == "shown"

    def test_whole_file_green_only_still_never_seen_failing(self):
        # Regression guard: the whole-file green-only shape stays
        # sentenced.
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_thing.py", "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py"
            " was never seen failing"
        )
