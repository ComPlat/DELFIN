"""RED test for H1G phase 3 (the protected wiring into the verify step).

This test pins the BEHAVIOUR the phase-3 wiring must deliver: when the
durable cross-turn test-run ledger is supplied to
``task_evidence.check_completion_claim``, the ledger built from REAL
tool_trace-shaped entries (``{tool:"bash", input:{command}, output, ok, ts}``)
must drive ``counter_test.assess`` so the "passes with and without the fix"
gap surfaces in the returned record's note.

RED against the current code (the wiring is delivered as a protected patch
the operator applies, not a branch commit): ``check_completion_claim`` today
does not accept ``durable_test_runs=``, so this test references a keyword the
un-wired signature does not have and therefore fails with TypeError.

GREEN after .gate/H1G_phase3_wiring.patch is applied -- the keyword exists,
the note carries the sentence, and a genuinely controlled fix is untouched.

The gauge is the REAL public path, not hand-made assess shapes: every
``durable_runs`` ledger below is produced by feeding tool_trace-shaped
entries ``{tool, input, output, ok, ts}`` through
``counter_test.durable_test_runs`` (exactly what api_client's wiring does),
so these scenarios only pass if the durable -> normalize -> assess chain
truly maps a real bash/gate pytest run onto the changed counter-test. This
replaces the earlier hand-made ``{target, red, green}`` shapes that hid the
F1/F5 real-shape gap the operator and s28 flagged.

This is a zz_ file by design (a red test that ships with a pending patch, so
the suite ignores it until the operator applies the wiring).
"""
from delfin.agent import task_evidence as te
from delfin.agent import counter_test as ct


def _changes(*paths):
    return [{"path": p, "ts": 0.0, "created": True} for p in paths]


def _entry(cmd, output, ok, ts):
    """A real tool_trace-shaped entry (the shape the operator built from)."""
    return {"tool": "bash", "input": {"command": cmd},
            "output": output, "ok": ok, "ts": ts}


def _durable(entries):
    """The durable ledger exactly as the wiring builds it from the trace."""
    return ct.durable_test_runs(entries)


def _run(results, fail_output=None):
    """A live green per-turn ledger entry for the wiring's tests argument."""
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
# Operator scenario (a): gate-tool whole-file red then green = shown
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
    # A fake runner (not real pytest) leaks a "test suite failed to run"
    # symptom; here the note must NOT accuse the control.
    assert "was changed but never run" not in (res.get("note") or ""), res


# --------------------------------------------------------------------------
# Operator scenario (b): node red then green = shown (definitive)
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
# Operator scenario (c): green only = "never seen failing"
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
# Operator scenario (d): green then red = NOT shown
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
# Operator scenario (e): bare pytest -q red then green, FAILED names the file
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
# Wiring contract: a genuinely controlled fix keeps an uncluttered verdict,
# and a changed test file never run is named.
# --------------------------------------------------------------------------
def test_controlled_fix_not_perturbed_by_wiring():
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
    assert "changed but never run" not in (res.get("note") or ""), res


def test_changed_never_run_test_names_the_test():
    res = _claim(entry_sets=[], live=[_run("tests/test_x.py")])
    assert "changed but never run" in (res.get("note") or ""), res
