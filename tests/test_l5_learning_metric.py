"""L5 learning metric: does a second run of the same task do better?

Drives ``delfin.agent.pack.benchmark.learning_metric.learning_verdict``
through its public path with constructed run records (one JSONL row as
``write_run`` writes it, PREREG_W15.md). The rules under test are fixed in
PREREG_W15.md:

* correctness dominates -- a second run that is faster but wrong is a fail,
* identical runs are ``no_learning``, never ``learned``,
* a record with a missing field is an ``error``, never a silent zero,
* first-vs-second must be the SAME task (task_id must match).
"""

from delfin.agent.pack.benchmark import learning_metric as lm

VERD = lm.learning_verdict


def _record(*, task_id="l5_read_the_value_from_file", success=True,
            quality_0_100=90, tool_calls=6, input_tokens=300,
            output_tokens=150, duration_s=20.0, cost_usd=0.12, verify_ok=None,
            unmeasured=False):
    """A valid run record -- one asdict(BenchmarkResult) row."""
    return {
        "task_id": task_id,
        "success": success,
        "quality_0_100": quality_0_100,
        "tool_calls": tool_calls,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
        "duration_s": duration_s,
        "cost_usd": cost_usd,
        "verify_ok": verify_ok,
        "unmeasured": unmeasured,
    }


# ---------------------------------------------------------------------------
# Correctness dominates: faster-but-wrong is a fail whatever it saved
# ---------------------------------------------------------------------------


def test_faster_but_wrong_fails():
    """Second run saves on every efficiency proxy but is wrong -> wrong."""
    first = _record(success=True, tool_calls=8, input_tokens=500,
                    output_tokens=300, duration_s=40.0)
    second = _record(success=False, tool_calls=3, input_tokens=90,
                     output_tokens=40, duration_s=8.0)
    out = VERD(first, second)
    assert out["verdict"] == "wrong"
    assert out["error"] is None
    assert out["evidence"]["second_success"] is False


def test_correctness_regression_is_wrong():
    """First passed, second fails -> wrong, no matter what it saved."""
    first = _record(success=True)
    second = _record(success=False, tool_calls=1, input_tokens=10,
                     output_tokens=5, duration_s=1.0)
    out = VERD(first, second)
    assert out["verdict"] == "wrong"


def test_verify_flip_true_to_false_is_wrong():
    """Second run has verify_ok=False (artifact failed verification) even
    though success stayed True -> wrong; correct the correction."""
    first = _record(success=True, verify_ok=True)
    second = _record(success=True, verify_ok=False, tool_calls=2,
                     input_tokens=40, output_tokens=20, duration_s=5.0)
    out = VERD(first, second)
    assert out["verdict"] == "wrong"


def test_verify_not_flipped_does_not_fire_wrong():
    """verify_ok absent on both (no accept script) must NOT be wrong."""
    first = _record(success=True)
    second = _record(success=True, tool_calls=5, input_tokens=260,
                     output_tokens=130, duration_s=18.0)
    out = VERD(first, second)
    assert out["verdict"] != "wrong"


def test_second_verify_false_is_wrong_even_when_first_had_none():
    """A second run whose artifact failed verification is wrong even if the
    first run had no verify signal at all (verify_ok=None is not a pass)."""
    first = _record(success=True)                    # verify_ok=None
    second = _record(success=True, verify_ok=False, tool_calls=2,
                     input_tokens=30, output_tokens=15, duration_s=4.0)
    out = VERD(first, second)
    assert out["verdict"] == "wrong"


def test_verify_none_is_not_a_correctness_signal():
    """verify_ok=None (no verify script ran) is neither a pass nor a fail:
    an otherwise-identical run scores no_learning, not learned, and never
    trips wrong."""
    first = _record(success=True, verify_ok=None)
    second = _record(success=True, verify_ok=None)
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"
    assert out["verdict"] != "wrong"


# ---------------------------------------------------------------------------
# Identical runs are no_learning, never learned
# ---------------------------------------------------------------------------


def test_identical_runs_are_no_learning():
    """Same task, same efficiency, both correct -> no_learning, not learned."""
    first = _record()
    second = _record()
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"


def test_only_correctness_improvement_is_not_learned():
    """First wrong, second correct, BUT identical efficiency -> the second
    is not wrong, yet nothing is faster: no_learning, not learned."""
    first = _record(success=False)
    second = _record(success=True)
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"


# ---------------------------------------------------------------------------
# Missing/unknown field is an error, never a silent zero
# ---------------------------------------------------------------------------


def test_missing_tool_calls_is_error():
    """A record without tool_calls cannot be compared; not a silent zero."""
    first = _record()
    second = _record(tool_calls=99)
    second.pop("tool_calls")
    out = VERD(first, second)
    assert out["verdict"] == "error"
    assert "tool_calls" in out["error"]


def test_missing_cost_usd_is_error():
    """A record without cost_usd cannot be compared; not a silent zero."""
    first = _record()
    second = _record(cost_usd=0.5)
    second.pop("cost_usd")
    out = VERD(first, second)
    assert out["verdict"] == "error"
    assert "cost_usd" in out["error"]


def test_complete_record_does_not_error():
    """The valid fixed record must NOT trip an error guard."""
    out = VERD(_record(), _record())
    assert out["verdict"] != "error"
    assert out["error"] is None


def test_missing_correctness_success_is_error():
    """The correctness field must be present and non-null; missing -> error."""
    first = _record()
    second = _record(success=False)
    second.pop("success")
    out = VERD(first, second)
    assert out["verdict"] == "error"
    assert "success" in out["error"]


def test_non_bool_success_is_error():
    """success typed wrong (a string) is an error, not silently truthy."""
    first = _record()
    second = _record(success="yes")
    out = VERD(first, second)
    assert out["verdict"] == "error"


def test_non_dict_record_is_error():
    """A record that is not a dict (None, a list) is an error, never an
    AttributeError -- the function promises an error verdict for anything it
    cannot compare."""
    first = _record()
    assert VERD(None, first)["verdict"] == "error"
    assert VERD(first, [])["verdict"] == "error"
    assert VERD(first, "x")["verdict"] == "error"


def test_bool_quality_is_error():
    """quality_0_100=True (a bool) is an error like the other numeric
    fields, not accepted as an int."""
    first = _record()
    second = _record(quality_0_100=True)
    out = VERD(first, second)
    assert out["verdict"] == "error"


def test_unmeasured_run_is_error():
    """A run that never reached the model is not a slower run."""
    first = _record()
    second = _record(unmeasured=True)
    out = VERD(first, second)
    assert out["verdict"] == "error"
    assert "unmeasured" in out["error"]


# ---------------------------------------------------------------------------
# Same task, or error
# ---------------------------------------------------------------------------


def test_task_id_mismatch_is_error():
    """Learning compares two runs of the SAME task; a mismatch is an error."""
    first = _record(task_id="l5_read_the_value_from_file")
    second = _record(task_id="l5_report_absent_quantity_as_absent")
    out = VERD(first, second)
    assert out["verdict"] == "error"
    assert "task_id" in out["error"]


# ---------------------------------------------------------------------------
# The positive case: faster + correct is learned
# ---------------------------------------------------------------------------


def test_faster_correct_is_learned():
    """Second run strictly better on a proxy and still correct -> learned."""
    first = _record(success=True, tool_calls=6, input_tokens=300,
                    output_tokens=150, duration_s=20.0)
    second = _record(success=True, tool_calls=4, input_tokens=180,
                     output_tokens=90, duration_s=12.0)
    out = VERD(first, second)
    assert out["verdict"] == "learned"
    assert out["error"] is None
    assert "tool_calls" in out["evidence"]["reason"]


# ---------------------------------------------------------------------------
# Efficiency regression: slower + correct is regressed, not learned
# ---------------------------------------------------------------------------


def test_slower_correct_is_regressed():
    """Second run correct but slower on every proxy -> regressed."""
    first = _record(success=True, tool_calls=4, input_tokens=180,
                    output_tokens=90, duration_s=12.0)
    second = _record(success=True, tool_calls=9, input_tokens=600,
                     output_tokens=400, duration_s=45.0)
    out = VERD(first, second)
    assert out["verdict"] == "regressed"


def test_mixed_efficiency_is_regressed_not_learned():
    """One proxy better, one worse is not clean learning -> regressed."""
    first = _record(success=True, tool_calls=6, input_tokens=300,
                    output_tokens=150, duration_s=20.0)
    second = _record(success=True, tool_calls=3, input_tokens=600,
                     output_tokens=400, duration_s=25.0)
    out = VERD(first, second)
    assert out["verdict"] == "regressed"
