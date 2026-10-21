"""check_completion_claim with durable_test_runs: the control-run sentence.

Input: tool-trace-shaped entries ``{tool, input, output, ok, ts}`` fed
through ``counter_test.durable_test_runs``, exactly as the executor's
``_check_completion`` builds them. Output checked: the note of the tests
verdict.
"""
from delfin.agent import task_evidence as te
from delfin.agent import counter_test as ct


def _changes(*paths):
    return [{"path": p, "ts": 0.0, "created": True} for p in paths]


def _entry(cmd, output, ok, ts):
    """A tool_trace-shaped entry."""
    return {"tool": "bash", "input": {"command": cmd},
            "output": output, "ok": ok, "ts": ts}


def _durable(entries):
    """The durable ledger exactly as the wiring builds it from the trace."""
    return ct.durable_test_runs(entries)


def _run(target):
    """A live green per-turn ledger entry for the tests argument."""
    return {"tool": "run_tests", "command": "tests/test_x.py -q",
            "exit_code": 0, "status": "ok", "passed": 1, "failed": 0,
            "ts": 1000.0}


def _claim(*, entry_sets, live):
    return te.check_completion_claim(
        "run the tests for the fix",
        "",
        changes=_changes("delfin/x.py", "tests/test_x.py"),
        tests=live,
        durable_test_runs=_durable(entry_sets),
        window_start=0.0,
    )


# --------------------------------------------------------------------------
# Scenario (a): gate-tool whole-file red then green = shown
# --------------------------------------------------------------------------
def test_gate_whole_file_red_then_green_is_a_control():
    res = _claim(
        entry_sets=[
            _entry(".delfin/gate-tools/gate tests/test_x.py -q",
                   "1 failed, 0 passed in 0.1s\nFAILED tests/test_x.py::test_a",
                   False, 1.0),
            _entry(".delfin/gate-tools/gate tests/test_x.py -q",
                   "1 passed, 0 failed in 0.1s", True, 2.0),
        ],
        live=[_run("tests/test_x.py")],
    )
    assert "was changed but never run" not in (res.get("note") or ""), res


# --------------------------------------------------------------------------
# Scenario (b): node red then green = shown (definitive)
# --------------------------------------------------------------------------
def test_node_red_then_green_is_a_control():
    res = _claim(
        entry_sets=[
            _entry("pytest tests/test_x.py::test_a -q",
                   "1 failed, 0 passed in 0.1s\nFAILED tests/test_x.py::test_a",
                   False, 1.0),
            _entry("pytest tests/test_x.py::test_a -q",
                   "1 passed, 0 failed in 0.1s", True, 2.0),
        ],
        live=[_run("tests/test_x.py::test_a")],
    )
    assert res["verdict"] == "verified"
    assert "never seen failing" not in (res.get("note") or ""), res


# --------------------------------------------------------------------------
# Scenario (c): green only = "never seen failing"
# --------------------------------------------------------------------------
def test_green_only_surfaces_never_seen_failing():
    res = _claim(
        entry_sets=[
            _entry("pytest tests/test_x.py::test_a -q",
                   "1 passed, 0 failed in 0.1s", True, 2.0),
        ],
        live=[_run("tests/test_x.py::test_a")],
    )
    assert "never seen failing" in (res.get("note") or ""), res
    assert "changed but never run" not in (res.get("note") or ""), res


# --------------------------------------------------------------------------
# Scenario (d): green then red = NOT shown
# --------------------------------------------------------------------------
def test_green_then_red_is_not_shown():
    res = _claim(
        entry_sets=[
            _entry("pytest tests/test_x.py::test_a -q",
                   "1 passed, 0 failed in 0.1s", True, 1.0),
            _entry("pytest tests/test_x.py::test_a -q",
                   "1 failed, 0 passed in 0.1s\nFAILED tests/test_x.py::test_a",
                   False, 2.0),
        ],
        live=[_run("tests/test_x.py::test_a")],
    )
    # The change BROKE a previously passing test -> never a control.
    assert "never seen failing" not in (res.get("note") or ""), res


# --------------------------------------------------------------------------
# Scenario (e): bare pytest -q red then green, FAILED names the file
# = shown
# --------------------------------------------------------------------------
def test_bare_pytest_red_then_green_with_failed_lines_is_shown():
    res = _claim(
        entry_sets=[
            _entry("pytest -q", "1 failed, 0 passed in 0.1s\n"
                                "FAILED tests/test_x.py::test_a",
                   False, 1.0),
            _entry("pytest -q", "1 passed, 0 failed in 0.1s", True, 2.0),
        ],
        live=[_run("tests/test_x.py")],
    )
    assert "changed but never run" not in (res.get("note") or ""), res


# --------------------------------------------------------------------------
# A changed test file never run is named.
# --------------------------------------------------------------------------
def test_changed_never_run_test_names_the_test():
    res = _claim(entry_sets=[], live=[_run("tests/test_x.py")])
    assert "changed but never run" in (res.get("note") or ""), res
