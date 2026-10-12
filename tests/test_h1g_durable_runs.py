"""H1G phase 3: durable test_runs reconstruction (counter_test.durable_test_runs).

The verify step's live test ledger (``_observe_test_evidence``) is
PER-TURN and replaced at every turn start, so the "seen red before the
change" observation may live in a PREVIOUS turn. ``assess`` needs that
cross-turn red to decide a control. ``durable_test_runs`` rebuilds the
run records from the DURABLE tool_trace so a red from an earlier turn is
inside what ``assess`` sees.

Covered here: the pure reconstruction — per-run records, the red/green
split, the dropping of non-test and evidence-less entries, deterministic
timestamps, and the end-to-end durable control (red turn + green turn ->
``shown``).
"""

from delfin.agent import counter_test as ct


def _entry(tool="bash", input_="pytest tests/test_thing.py -q",
           output="2 passed, 0 failed in 0.2s", ok=True, ts=1000.0):
    return {"tool": tool, "input": input_, "output": output,
            "ok": ok, "ts": ts}


# --------------------------------------------------------------------------
# basic reconstruction
# --------------------------------------------------------------------------

def test_green_run_is_recorded_green():
    runs = ct.durable_test_runs([_entry(output="2 passed, 0 failed in 0.2s")])
    assert len(runs) == 1
    r = runs[0]
    assert r["command"] == "pytest tests/test_thing.py -q"
    assert r["failed"] == 0
    assert r["passed"] == 2
    assert r["exit_code"] == 0
    assert r["status"] == "ok"


def test_red_run_is_recorded_red():
    runs = ct.durable_test_runs([
        _entry(output="1 failed, 1 passed in 0.2s"
                     "\nFAILED tests/test_thing.py::test_y", ok=False)])
    assert len(runs) == 1
    r = runs[0]
    assert r["failed"] == 1
    assert r["status"] == "failed"
    assert r["exit_code"] == 1


def test_non_zero_exit_without_summary_still_red():
    # A failed run whose output lacks the "N failed" line (e.g. a crash or
    # an early exit) is still red via its non-zero exit (api_client :11179).
    runs = ct.durable_test_runs([
        _entry(output="INTERNALERROR: ...", ok=False)])
    assert len(runs) == 1
    assert runs[0]["failed"] == 1
    assert runs[0]["status"] == "failed"


# --------------------------------------------------------------------------
# dropping non-evidence
# --------------------------------------------------------------------------

def test_non_test_command_is_dropped():
    assert ct.durable_test_runs([_entry(input_="ls -la", output="")]) == []


def test_collect_only_without_summary_is_dropped():
    # A run that demonstrated no tests (no "N passed"/"N failed") is not
    # evidence and must not masquerade as a green run
    # (api_client :11189-11192).
    runs = ct.durable_test_runs([
        _entry(input_="echo hi", output="")])  # non-test command, dropped
    assert runs == []
    runs2 = ct.durable_test_runs([
        _entry(input_="python -m pytest --collect-only",
               output="collected 3 items", ok=True)])
    # "--collect-only" matches the pytest command regex but carries no
    # pass/fail summary, so it is dropped too.
    assert runs2 == []


def test_missing_output_is_dropped():
    assert ct.durable_test_runs([_entry(output="")]) == []


# --------------------------------------------------------------------------
# durable across turns + determinism
# --------------------------------------------------------------------------

def test_red_and_green_from_separate_entries_both_survive():
    # A red run from an earlier turn and a green run from this turn — both
    # must survive the reconstruction so assess can see the control.
    runs = ct.durable_test_runs([
        _entry(output="1 failed, 0 passed in 0.1s"
                      "\nFAILED tests/test_thing.py::test_new", ok=False,
               ts=1.0),
        _entry(output="1 passed, 0 failed in 0.1s",
                     ok=True, ts=2.0),
    ])
    assert len(runs) == 2
    assert runs[0]["ts"] == 1.0
    assert runs[1]["ts"] == 2.0
    assert runs[0]["failed"] == 1
    assert runs[1]["passed"] == 1


def test_now_pins_timestamps():
    runs = ct.durable_test_runs(
        [_entry(output="1 passed in 0.1s", ts=999.0)], now=500.0)
    assert runs[0]["ts"] == 500.0


# --------------------------------------------------------------------------
# end-to-end durable control -> assess says "shown"
# --------------------------------------------------------------------------

def test_durable_red_then_green_node_in_changed_test_is_shown():
    # The phase-3 payoff: a red run in an EARLIER turn (only visible via
    # the durable trace) plus a green run this turn, on the fix's own
    # changed counter-test node, must read as a completed control.
    changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
    runs = ct.durable_test_runs([
        _entry(output="1 failed, 0 passed in 0.1s"
                      "\nFAILED tests/test_thing.py::test_new", ok=False,
               ts=1.0),
        _entry(output="1 passed, 0 failed in 0.1s", ok=True, ts=2.0),
    ])
    # The per-run records only name the whole file in the command; a real
    # durable control run names the SAME node red (earlier turn) and green
    # (this turn). Refine both to the counter-test node so the control is
    # genuinely on one node.
    for r in runs:
        r["command"] = "tests/test_thing.py::test_new"
    test_runs = ct.normalize_test_runs(runs)
    assert ct.assess(changed, test_runs) == "shown"


def test_green_only_durable_node_is_not_shown():
    # A durable trace showing only green on the fix's own counter-test is
    # the classic "passes with and without the fix": not a control.
    changed = ["delfin/agent/thing.py", "tests/test_thing.py"]
    runs = ct.durable_test_runs([
        _entry(output="1 passed, 0 failed in 0.1s", ok=True, ts=1.0),
    ])
    for r in runs:
        r["command"] = "tests/test_thing.py::test_new"
    test_runs = ct.normalize_test_runs(runs)
    assert ct.assess(changed, test_runs).startswith(
        "no control run: test tests/test_thing.py::test_new was never")
