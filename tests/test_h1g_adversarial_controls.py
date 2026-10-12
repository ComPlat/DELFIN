"""H1G phase-2 adversarial controls: the fix's own counter-test decides.

S28's three QS-confirmed findings (A, B, C) converge on one root: the
phase-2 `assess` returned "shown" for a red+green on ANY target, so a
fix whose own counter-test was only ever green (or never run) was masked
by a completed control on an UNRELATED target or an unrelated
same-file node. The mandate: a control counts only on the fix's own
changed counter-test (node-level); a whole-file target is never a
completed control.

These tests pin the fixed behaviour: each asserts the fix's green-only /
unrun counter-test is NOT masked, and a whole-file red+green does NOT
read as "shown".
"""

from delfin.agent.counter_test import assess


class TestUnrelatedTargetDoesNotMask:
    """Finding C: a controlled UNRELATED file cannot mask a green-only
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
    """Finding A: a red+green on an OLD node in the SAME file must not
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
        # green-only NEW node still owns the verdict (findings A + B).
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
    """Finding B: a whole-file target (no ::node) is never a completed
    control -- it cannot say which node carried the red."""

    def test_whole_file_single_aggregate_not_shown(self):
        changed = ["delfin/agent/thing.py"]
        runs = [{"target": "tests/test_thing.py", "red": True, "green": True}]
        result = assess(changed, runs)
        assert result != "shown"
        assert result.startswith("no control run:")

    def test_whole_file_red_green_not_shown_even_in_scope(self):
        # In-scope (test file changed) too: a whole-file red+green never
        # reads as a control for the fix.
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_thing.py", "red": True, "green": True}]
        result = assess(changed, runs)
        assert result != "shown"
        assert result.startswith("no control run:")

    def test_whole_file_green_only_still_never_seen_failing(self):
        # Regression guard: the whole-file green-only shape that phase-2
        # already sentenced stays sentenced.
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_thing.py", "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py"
            " was never seen failing"
        )
