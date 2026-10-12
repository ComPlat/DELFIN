"""RED test for H1G phase 3 (the protected wiring into the verify step).

This test pins the BEHAVIOUR the phase-3 wiring must deliver: when the
durable cross-turn test-run ledger is supplied to
``task_evidence.check_completion_claim``, a green run on a fix whose own
counter-test was only ever green (never red) must NOT read as plain
"verified" -- the fix's "passes with and without the fix" gap must surface
in the returned record's note via ``counter_test.assess``.

RED against the current code (the wiring is delivered as a protected patch
the operator applies, not a branch commit): ``check_completion_claim`` today
does not accept ``durable_test_runs=``, so this test references a keyword the
un-wired signature does not have and therefore fails with TypeError.

GREEN after .gate/H1G_phase3_wiring.patch is applied -- the keyword exists,
the note carries the sentence, and a genuinely controlled fix is untouched.

The ledger shapes below mirror the real verify step (task_evidence.py): the
subject must match ``_TEST_TASK_RE`` ("tests"/"verify"), ``tests`` is a LIST
of per-run dicts ``{tool, command, exit_code, status, passed, failed, ts}``
(exactly `_observe_test_evidence` emits, api_client.py:11126), and
``durable_test_runs`` is the assess-shaped list
``{target, red, green}`` (normalize_test_runs output).

This is a zz_ file by design (a red test that ships with a pending patch, so
the suite ignores it until the operator applies the wiring; it is re-collected
the moment the signature gains the keyword).
"""
from delfin.agent import task_evidence as te


def _changes(*paths):
    return [{"path": p, "ts": 0.0, "created": True} for p in paths]


def _green_run(target):
    """A fresh green run record in the shape _observe_test_evidence emits."""
    return {"tool": "run_tests", "command": target, "exit_code": 0,
            "status": "ok", "passed": 1, "failed": 0, "ts": 100.0}


def _claim(*, changes, tests_runs, durable_runs):
    return te.check_completion_claim(
        "run the tests for the fix",
        "",
        changes=changes,
        tests=tests_runs,
        durable_test_runs=durable_runs,
        window_start=0.0,
    )


def test_green_even_with_durable_runs_surfaces_sentence():
    """A green run on a changed-fix test never seen RED is not a completed
    control: the note must name the gap."""
    # The fix changed the module AND its counter-test; the durable ledger
    # saw test_thing::test_x GREEN only (never red) -- "passes with and
    # without the fix".
    res = _claim(
        changes=_changes("delfin/agent/thing.py", "tests/test_thing.py"),
        tests_runs=[_green_run("tests/test_thing.py::test_x")],
        durable_runs=[{"target": "tests/test_thing.py::test_x",
                       "red": False, "green": True}],
    )
    assert res["verdict"] == "verified"
    assert "never seen failing" in (res.get("note") or ""), res


def test_controlled_fix_not_perturbed_by_wiring():
    """A genuinely controlled fix (counter-test node red AND green) keeps
    an uncluttered verified record."""
    res = _claim(
        changes=_changes("delfin/agent/thing.py", "tests/test_thing.py"),
        tests_runs=[_green_run("tests/test_thing.py::test_x")],
        durable_runs=[{"target": "tests/test_thing.py::test_x",
                       "red": True, "green": True}],
    )
    assert res["verdict"] == "verified"
    assert "never seen failing" not in (res.get("note") or ""), res


def test_changed_never_run_test_names_the_test():
    """A test file changed together with the fix but never run is named."""
    res = _claim(
        changes=_changes("delfin/agent/thing.py", "tests/test_thing.py"),
        tests_runs=[_green_run("tests/test_thing.py::test_x")],
        durable_runs=[],
    )
    assert "changed but never run" in (res.get("note") or ""), res
