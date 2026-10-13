"""counter_test.assess on hand-made run records.

Shapes: green only -> "was never seen failing"; red then green on one
target -> "shown"; no product file changed -> None; a changed test that
never ran -> "was changed but never run"; a red and a green on different
nodes of one file are no control; a red from an earlier turn plus a green
now is a control.
"""

from delfin.agent.counter_test import assess


class TestGreenOnlyFix:
    """A fix whose test was only ever seen green: no control on record."""

    def test_product_plus_test_changed_test_green_only(self):
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_thing.py", "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py was never seen failing"
        )

    def test_red_flag_absent_counts_as_never_red(self):
        changed = ["delfin/agent/thing.py"]
        runs = [{"target": "tests/test_thing.py::test_x", "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py::test_x"
            " was never seen failing"
        )


class TestRedThenGreen:
    """A control run is on record: the same target seen red, then green."""

    def test_two_runs_red_then_green_same_node(self):
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [
            {"target": "tests/test_thing.py::test_x", "red": True},
            {"target": "tests/test_thing.py::test_x", "green": True},
        ]
        assert assess(changed, runs) == "shown"

    def test_red_from_previous_turn_green_this_turn(self):
        # The red was recorded in a PREVIOUS turn.
        # The durable ledger (tool_trace) carries it; the caller puts it
        # into test_runs, so assess sees red-before + green-now.
        changed = ["delfin/agent/thing.py"]
        runs = [
            {"target": "tests/test_thing.py::test_x",
             "red": True},          # recorded last turn (durable)
            {"target": "tests/test_thing.py::test_x",
             "green": True},         # recorded this turn
        ]
        assert assess(changed, runs) == "shown"

    def test_control_wins_over_a_green_only_unrelated_target(self):
        changed = ["delfin/agent/thing.py"]
        runs = [
            {"target": "tests/test_thing.py::test_x", "red": True},
            {"target": "tests/test_thing.py::test_x", "green": True},
            {"target": "tests/test_other.py::test_y", "green": True},
        ]
        assert assess(changed, runs) == "shown"


class TestRedAndGreenOnDifferentNodesIsNotAControl:
    """Per node, not per file."""

    def test_red_on_one_node_green_on_another_same_file(self):
        # A red run of test_a and a green run of test_b in the SAME file
        # must NOT combine into a control.
        changed = ["delfin/agent/thing.py"]
        runs = [
            {"target": "tests/test_thing.py::test_a", "red": True},
            {"target": "tests/test_thing.py::test_b", "green": True},
        ]
        assert assess(changed, runs).startswith("no control run:")

    def test_whole_file_red_then_single_node_green_is_not_whole_control(self):
        # A whole-file red run does not pair with a green single-node run
        # into a control for the passing node (the file target and the
        # node target differ).
        changed = ["delfin/agent/thing.py"]
        runs = [
            {"target": "tests/test_thing.py", "red": True},
            {"target": "tests/test_thing.py::test_b", "green": True},
        ]
        assert assess(changed, runs).startswith("no control run:")


class TestNoProductFileChanged:
    """A turn that changed no product file adds nothing (None)."""

    def test_only_a_test_file_changed(self):
        changed = ["tests/test_thing.py"]
        runs = [{"target": "tests/test_thing.py", "green": True}]
        assert assess(changed, runs) is None

    def test_nothing_changed(self):
        assert assess([], []) is None

    def test_only_a_readme_changed(self):
        # A README is documentation, not the fix: no control is owed.
        changed = ["README.md"]
        runs = [{"target": "tests/test_thing.py", "red": True, "green": True}]
        assert assess(changed, runs) is None


class TestTestChangedButNeverRun:
    """A test file changed together with the fix but never run."""

    def test_changed_test_file_with_no_run(self):
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        assert assess(changed, []) == (
            "no control run: test tests/test_thing.py was changed but never run"
        )

    def test_changed_test_file_never_covered_but_other_test_run(self):
        # The fix's own counter-test sits unrun even though another ran.
        changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
        runs = [{"target": "tests/test_other.py", "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py was changed but never run"
        )


class TestRelevantCorners:
    """Red-but-never-green, and a product change with no run at all."""

    def test_seen_failing_but_never_passing_after(self):
        changed = ["delfin/agent/thing.py"]
        runs = [{"target": "tests/test_thing.py::test_x", "red": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py::test_x was seen "
            "failing but never passing after the change"
        )

    def test_product_change_with_no_test_run(self):
        changed = ["delfin/agent/thing.py"]
        assert assess(changed, []) == (
            "no control run: no test run was recorded for the change"
        )

    def test_garbage_runs_are_ignored_not_fatal(self):
        changed = ["delfin/agent/thing.py"]
        runs = [None, "junk", {"target": "tests/test_thing.py", "green": True}]
        assert assess(changed, runs) == (
            "no control run: test tests/test_thing.py was never seen failing"
        )
