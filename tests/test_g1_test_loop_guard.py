"""The fifth loop guard: a test that keeps failing identically.

The degenerate coding loop -- edit a line, rerun the same test, fail the same
way -- is invisible to the four existing guards: a failing ``run_tests`` is a
SUCCESSFUL tool call (not an ``{"error"...}``), and the model's edit differs
every round (so the identical ``(name, args)`` check never sees a repeat).
``test_loop_guard`` counts consecutive test results with the SAME failing-node
set AND the SAME first failure message, steering once at 3 and aborting at 5.

This file unit-tests the pure guard. The integration (the guard firing in a
scripted headless run) lives in ``test_g1_degenerate_loop_fires_in_a_run.py``.
"""

from __future__ import annotations

from delfin.agent.test_loop_guard import (
    SameFailureCounter,
    latest_signature,
    report_text_signature,
    result_signature,
    step,
)

FAILING = {
    "summary": {"passed": 0, "failed": 1},
    "failures": [
        {"node_id": "tests/test_x.py::test_a",
         "message": "assert 1 == 2"},
    ],
}


def test_signature_composes_nodes_and_first_message():
    sig = result_signature(FAILING)
    assert "tests/test_x.py::test_a" in sig
    assert "assert 1 == 2" in sig
    # A different message for the same node is a different failure.
    other = result_signature({"failures": [
        {"node_id": "tests/test_x.py::test_a",
         "message": "assert 1 == 3"}]})
    assert sig != other


def test_signature_is_order_independent_but_node_set_sensitive():
    ordered = result_signature({"failures": [
        {"node_id": "b.py::t2", "message": "m2"},
        {"node_id": "a.py::t1", "message": "m1"}]})
    reordered = result_signature({"failures": [
        {"node_id": "a.py::t1", "message": "m1"},
        {"node_id": "b.py::t2", "message": "m2"}]})
    assert ordered == reordered
    # Same nodes, one test additionally failing -> NOT the same failure.
    assert ordered != result_signature({"failures": [
        {"node_id": "a.py::t1", "message": "m1"},
        {"node_id": "b.py::t2", "message": "m2"},
        {"node_id": "c.py::t3", "message": "m3"}]})


def test_signature_empty_when_nothing_fails():
    assert result_signature({"summary": {"passed": 5, "failed": 0},
                             "failures": []}) == ""
    assert result_signature({"FAILURE": []}) == ""


def test_report_text_signature_reads_console_shape():
    text = (
        "________________ test_thing ________________\n"
        "tests/test_x.py:3: in test_thing\n"
        "    assert 1 == 2\n"
        "E   assert 1 == 2\n"
        "________________ test_other ________________\n"
        "____\n"
        "\n"
        "FAILED tests/test_x.py::test_thing - assert 1 == 2\n"
        "FAILED tests/test_y.py::test_other - assert 3 == 4\n"
    )
    sig = report_text_signature(text)
    assert "tests/test_x.py::test_thing" in sig
    assert "tests/test_y.py::test_other" in sig
    assert "assert 1 == 2" in sig


def test_report_text_ignores_ansi_and_unparsable():
    import re
    ansi = re.sub(r"E\s+", "E ", (
        "\x1b[31mFAILED\x1b[0m tests/test_x.py::test_a - x\n"))
    assert report_text_signature(ansi) != ""
    assert report_text_signature("no failure shape here") == ""


def test_latest_signature_is_the_newest_test_entry():
    ev = [
        {"tool": "read_file", "command": "", "exit_code": 0, "status": "ok",
         "passed": 0, "failed": 0, "ts": 0.0},
        {"tool": "run_tests", "test_sig": "S1", "ts": 0.0},
        {"tool": "bash", "ts": 0.0},
    ]
    assert latest_signature(ev) == "S1"
    assert latest_signature([]) == ""


def test_counter_steers_at_3_and_aborts_at_5_once():
    c = SameFailureCounter()
    assert c.observe("S") == ""
    assert c.observe("S") == ""
    assert c.observe("S") == "steer"     # the only steer
    assert c.observe("S") == ""
    assert c.observe("S") == "abort"     # and from here on it never steers
    assert c.observe("S") == "abort"


def test_counter_resets_on_a_different_failure():
    c = SameFailureCounter()
    c.observe("S")
    c.observe("S")
    c.observe("S")
    c.observe("T")                       # a new failure resets the run
    c.observe("T")
    assert c.observe("T") == "steer"     # three T's, third steers


def test_counter_resets_on_an_empty_signature():
    c = SameFailureCounter()
    c.observe("S")
    c.observe("S")
    c.observe("S")
    c.observe("")                        # a green / unknown run resets
    c.observe("S")
    c.observe("S")
    assert c.observe("S") == "steer"


def test_step_only_advances_when_a_new_test_entry_appeared():
    state: dict = {}
    ev: list = []
    assert step(ev, state) == ""
    ev.append({"tool": "run_tests", "test_sig": "S", "ts": 1.0})
    assert step(ev, state) == ""         # count 1
    ev.append({"tool": "run_tests", "test_sig": "S", "ts": 3.0})
    assert step(ev, state) == ""         # count 2
    ev.append({"tool": "run_tests", "test_sig": "S", "ts": 4.0})
    assert step(ev, state) == "steer"    # count 3
    ev.append({"tool": "run_tests", "test_sig": "S", "ts": 5.0})
    assert step(ev, state) == ""
    ev.append({"tool": "run_tests", "test_sig": "S", "ts": 6.0})
    assert step(ev, state) == "abort"    # count 5
