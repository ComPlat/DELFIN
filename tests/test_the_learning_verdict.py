"""The learning verdict on two constructed run records (PREREG_W15.md).

Drives ``learning_metric.learning_verdict``. Records are rows as
``write_run`` writes them for a task run with five replicates; unless a
test says otherwise every replicate has the same counts, so a difference in
the scalar counts is also a complete separation of the samples. The spread
rule itself is pinned in
tests/test_a_learning_verdict_needs_a_difference_beyond_the_spread.py.

Rules: correctness dominates (a wrong second run is ``wrong`` whatever it
saved); only ``tool_calls`` and ``tokens_total`` decide, never duration or
cost; equal counted work is ``no_learning``; any separated worsening is
``regressed``; a missing or mistyped field is ``error``, never a zero.
"""

from delfin.agent.pack.benchmark import learning_metric as lm

VERD = lm.learning_verdict
N = 5


def _rec(*, task_id="l5_read_the_value_from_file", success=True,
         quality_0_100=90, tool_calls=6, input_tokens=300,
         output_tokens=150, duration_s=20.0, cost_usd=0.12, verify_ok=None,
         unmeasured=False, **extra):
    """A valid five-replicate run record with identical replicates."""
    out = {
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
        "n_samples": N,
        "per_run_tool_calls": [tool_calls] * N,
        "per_run_tokens_total": (
            [input_tokens + output_tokens] * N
            if isinstance(input_tokens, int) and isinstance(output_tokens, int)
            else [0] * N),
    }
    out.update(extra)
    return out


# --- correctness dominates ---------------------------------------------------


def test_faster_but_wrong_is_wrong():
    first = _rec(tool_calls=8, input_tokens=500, output_tokens=300,
                 duration_s=40.0)
    second = _rec(success=False, tool_calls=3, input_tokens=90,
                  output_tokens=40, duration_s=8.0)
    out = VERD(first, second)
    assert out["verdict"] == "wrong"
    assert out["error"] is None
    assert out["evidence"]["second_success"] is False


def test_a_second_run_that_failed_verification_is_wrong():
    for first_verify in (True, None):
        first = _rec(verify_ok=first_verify)
        second = _rec(verify_ok=False, tool_calls=2, input_tokens=40,
                      output_tokens=20)
        assert VERD(first, second)["verdict"] == "wrong", first_verify


def test_verify_none_is_neither_a_pass_nor_a_fail():
    assert VERD(_rec(), _rec())["verdict"] == "no_learning"
    assert VERD(_rec(), _rec(tool_calls=5, input_tokens=260,
                             output_tokens=130))["verdict"] == "learned"


def test_explicit_verify_true_still_scores_on_counted_work():
    out = VERD(_rec(verify_ok=True),
               _rec(verify_ok=True, tool_calls=3, input_tokens=150,
                    output_tokens=70))
    assert out["verdict"] == "learned"


def test_only_the_second_runs_correctness_is_judged():
    first = _rec(verify_ok=False, tool_calls=6)
    second = _rec(verify_ok=None, tool_calls=4)
    assert VERD(first, second)["verdict"] == "learned"
    assert VERD(_rec(success=False), _rec())["verdict"] == "no_learning"


def test_a_lower_quality_under_success_is_evidence_not_a_verdict():
    first = _rec(quality_0_100=90)
    second = _rec(quality_0_100=40, tool_calls=3, input_tokens=150,
                  output_tokens=70)
    out = VERD(first, second)
    assert out["verdict"] == "learned"
    assert out["evidence"]["second_quality_0_100"] == 40


# --- only counted work decides -----------------------------------------------


def test_duration_and_cost_never_decide():
    for dur, cost in ((20.000037, 0.120004), (25.0, 0.2), (15.0, 0.05),
                      (19.99, 0.119), (20.0, 0.12)):
        out = VERD(_rec(), _rec(duration_s=dur, cost_usd=cost))
        assert out["verdict"] == "no_learning", (dur, cost)


def test_fewer_tool_calls_is_learned_whatever_the_clock_says():
    first = _rec(tool_calls=6, duration_s=20.0, cost_usd=0.12)
    second = _rec(tool_calls=4, input_tokens=260, output_tokens=130,
                  duration_s=25.0, cost_usd=0.15)
    out = VERD(first, second)
    assert out["verdict"] == "learned"
    assert "tool_calls" in out["evidence"]["reason"]
    ev = out["evidence"]
    assert ev["first"]["duration_s"] == 20.0
    assert ev["second"]["duration_s"] == 25.0
    assert ev["second"]["cost_usd"] == 0.15


def test_tokens_total_is_input_plus_output():
    # A split that nets equal is a tie...
    assert VERD(_rec(input_tokens=300, output_tokens=150),
                _rec(input_tokens=270, output_tokens=180)
                )["verdict"] == "no_learning"
    # ...a lower sum is learned.
    assert VERD(_rec(input_tokens=300, output_tokens=150),
                _rec(input_tokens=200, output_tokens=100)
                )["verdict"] == "learned"


def test_better_on_both_counted_fields_is_learned():
    out = VERD(_rec(tool_calls=6, input_tokens=300, output_tokens=150),
               _rec(tool_calls=4, input_tokens=200, output_tokens=100))
    assert out["verdict"] == "learned"


def test_worse_on_counted_work_is_regressed():
    assert VERD(_rec(tool_calls=4, input_tokens=180, output_tokens=90),
                _rec(tool_calls=9, input_tokens=600, output_tokens=400)
                )["verdict"] == "regressed"
    out = VERD(_rec(output_tokens=150), _rec(output_tokens=250))
    assert out["verdict"] == "regressed"
    assert out["error"] is None


def test_better_on_one_field_and_worse_on_the_other_is_regressed():
    out = VERD(_rec(tool_calls=6, output_tokens=150),
               _rec(tool_calls=4, output_tokens=400))
    assert out["verdict"] == "regressed"
    assert "worse" in out["evidence"]["reason"]


# --- an unreadable record is an error ----------------------------------------


def test_a_complete_record_is_not_an_error():
    out = VERD(_rec(), _rec())
    assert out["verdict"] != "error"
    assert out["error"] is None


def test_a_record_that_is_not_a_dict_is_an_error():
    for first, second, side in ((None, _rec(), "first"),
                                (["x"], _rec(), "first"),
                                (_rec(), [], "second"),
                                (_rec(), "x", "second")):
        out = VERD(first, second)
        assert out["verdict"] == "error"
        assert side in out["error"]


def test_a_missing_field_is_an_error_on_either_side():
    fields = ("success", "tool_calls", "input_tokens", "output_tokens",
              "duration_s", "cost_usd", "quality_0_100", "n_samples",
              "per_run_tool_calls", "per_run_tokens_total")
    for field in fields:
        for side in (0, 1):
            pair = [_rec(), _rec(tool_calls=4)]
            pair[side].pop(field)
            out = VERD(*pair)
            assert out["verdict"] == "error", (field, side)
            assert out["error"], (field, side)


def test_a_mistyped_field_is_an_error():
    for field, value in (("success", "yes"), ("success", 1),
                         ("quality_0_100", True), ("quality_0_100", "90"),
                         ("tool_calls", True), ("duration_s", True)):
        out = VERD(_rec(), _rec(**{field: value}))
        assert out["verdict"] == "error", (field, value)


def test_an_unmeasured_run_is_an_error_on_either_side():
    out = VERD(_rec(), _rec(unmeasured=True))
    assert out["verdict"] == "error"
    assert "unmeasured" in out["error"]
    assert VERD(_rec(unmeasured=True), _rec())["verdict"] == "error"


def test_two_different_tasks_are_an_error():
    out = VERD(_rec(task_id="l5_read_the_value_from_file"),
               _rec(task_id="l5_report_absent_quantity_as_absent"))
    assert out["verdict"] == "error"
    assert "task_id mismatch" in out["error"]


def test_an_unrelated_extra_field_changes_nothing():
    out = VERD(_rec(), _rec(tool_calls=4, walltime_h=3))
    assert out["verdict"] == "learned"
    assert out["error"] is None


def test_the_inputs_are_not_mutated():
    first = _rec()
    second = _rec(tool_calls=4, input_tokens=200, output_tokens=100)
    before = (repr(first), repr(second))
    VERD(first, second)
    assert (repr(first), repr(second)) == before
