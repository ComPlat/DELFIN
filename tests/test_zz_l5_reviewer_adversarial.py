"""Reviewer adversarial tests for package L5's learning metric.

Independent of tests/test_l5_learning_metric.py (the builder's suite). These
probe the contract fixed in PREREG_W15.md beyond the builder's 18 cases: the
mixed better+worse verdict, the tokens_total-derived-proxy semantics, the
verify_ok asymmetry (only the SECOND run's correctness is judged), the
pure-function guarantee (no input mutation), and the both-direction of every
missing/wrong-type guard.

Contract under test (pre-registered, PREREG_W15.md):
  verdict in {error, wrong, learned, no_learning, regressed};
  correctness dominates (a wrong second run fails whatever it saved);
  identical runs are no_learning, never learned;
  missing/unknown field is an error, never a silent zero;
  first-vs-second must be the SAME task (task_id must match);
  the metric is pure (reads only its two dicts, mutates neither).
"""

from delfin.agent.pack.benchmark import learning_metric as lm

VERD = lm.learning_verdict


def _rec(*, task_id="l5_read_the_value_from_file", success=True,
         quality_0_100=90, tool_calls=6, input_tokens=300,
         output_tokens=150, duration_s=20.0, cost_usd=0.12, verify_ok=None,
         unmeasured=False, **extra):
    """A valid run record -- one asdict(BenchmarkResult) row."""
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
    }
    out.update(extra)
    return out


# ---------------------------------------------------------------------------
# 1. Mixed better+worse -> regressed (the reviewer question the QS flagged)
# ---------------------------------------------------------------------------


def test_mixed_better_and_worse_is_regressed():
    """Second is faster (tool_calls down) but costlier (cost_usd up) -> the
    pre-registration names a partially-worse run 'regressed', not the all-
    equal 'no_learning' and not 'learned'."""
    first = _rec(tool_calls=6, cost_usd=0.12)
    second = _rec(tool_calls=4, cost_usd=0.20)  # better on calls, worse on cost
    out = VERD(first, second)
    assert out["verdict"] == "regressed"
    assert out["error"] is None
    assert "worse" in out["evidence"]["reason"]


# ---------------------------------------------------------------------------
# 2. tokens_total is the derived proxy input+output; a split that nets equal
# ---------------------------------------------------------------------------


def test_tokens_split_that_nets_equal_is_no_learning():
    """Input down, output up, sum equal -> the derived tokens_total proxy is a
    tie, and with all other proxies equal that is no_learning, not learned."""
    first = _rec(input_tokens=300, output_tokens=150)
    second = _rec(input_tokens=270, output_tokens=180)  # sum 450 both
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"
    assert out["error"] is None


def test_tokens_total_improves_to_learned():
    """The derived total strictly down with other proxies equal -> learned."""
    first = _rec(input_tokens=300, output_tokens=150)
    second = _rec(input_tokens=200, output_tokens=100)  # total 300 < 450
    out = VERD(first, second)
    assert out["verdict"] == "learned"


# ---------------------------------------------------------------------------
# 3. verify_ok asymmetry: only the SECOND run's artifact is judged
# ---------------------------------------------------------------------------


def test_first_verify_ok_false_not_wrong():
    """A first run that failed its own verify is still a valid baseline; the
    second run's correctness is what the verdict judges. second is correct and
    faster -> learned, not wrong on account of first's verify_ok False."""
    first = _rec(verify_ok=False, tool_calls=6)
    second = _rec(verify_ok=None, tool_calls=4)
    out = VERD(first, second)
    assert out["verdict"] == "learned"


def test_verify_ok_none_never_a_verdict():
    """verify_ok absent or None is never a correctness signal -- neither a
    pass nor a fail. A second run that is otherwise identical (no efficiency
    gain) must be no_learning, not skipped to learned."""
    first = _rec(verify_ok=None)
    second = _rec(verify_ok=None)
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"


# ---------------------------------------------------------------------------
# 4. The metric is a pure function: it must not mutate its inputs
# ---------------------------------------------------------------------------


def test_pure_no_input_mutation():
    """learning_verdict must leave both record dicts unchanged."""
    first = _rec()
    second = _rec(tool_calls=4, input_tokens=200, output_tokens=100,
                  duration_s=12.0, cost_usd=0.08)
    before_first = dict(first)
    before_second = dict(second)
    VERD(first, second)
    assert first == before_first
    assert second == before_second


# ---------------------------------------------------------------------------
# 5. Every guard in both directions (missing and wrong-type)
# ---------------------------------------------------------------------------


def test_success_missing_is_error():
    first = dict(_rec())
    first.pop("success")
    out = VERD(first, _rec(tool_calls=4))
    assert out["verdict"] == "error"
    assert out["error"]


def test_success_int_is_error():
    """success carried as int 1 (bool-like) is not a bool -> error."""
    first = _rec(success=1)
    out = VERD(first, _rec(tool_calls=4))
    assert out["verdict"] == "error"


def test_quality_wrong_type_is_error():
    second = _rec(quality_0_100="90")  # str, not int
    out = VERD(_rec(), second)
    assert out["verdict"] == "error"


def test_tool_calls_bool_is_error():
    """bool is an int subclass; _get_int must reject it as a mistyped value."""
    second = _rec(tool_calls=True)
    out = VERD(_rec(), second)
    assert out["verdict"] == "error"
    assert out["error"]


def test_missing_tokens_is_error_on_either_side():
    for missing in ("input_tokens", "output_tokens"):
        second = dict(_rec())
        second.pop(missing)
        assert VERD(_rec(), second)["verdict"] == "error", missing
        first = dict(_rec())
        first.pop(missing)
        assert VERD(first, _rec(tool_calls=4))["verdict"] == "error", missing


def test_missing_duration_and_cost_are_error():
    for field in ("duration_s", "cost_usd"):
        second = dict(_rec())
        second.pop(field)
        out = VERD(_rec(), second)
        assert out["verdict"] == "error", field


def test_duration_bool_is_error():
    first = _rec(duration_s=True)
    assert VERD(first, _rec(tool_calls=4))["verdict"] == "error"


def test_unmeasured_on_either_side_is_error():
    first = _rec(unmeasured=True)
    out = VERD(first, _rec(tool_calls=4))
    assert out["verdict"] == "error"


def test_task_id_mismatch_is_error():
    out = VERD(_rec(task_id="a"), _rec(task_id="b"))
    assert out["verdict"] == "error"
    assert "task_id mismatch" in out["error"]


# ---------------------------------------------------------------------------
# 6. Verdict boundaries on the economy strict-inequality rule
# ---------------------------------------------------------------------------


def test_second_better_on_all_proxies_is_learned():
    first = _rec(tool_calls=6, input_tokens=300, output_tokens=150,
                 duration_s=20.0, cost_usd=0.12)
    second = _rec(tool_calls=4, input_tokens=200, output_tokens=100,
                  duration_s=12.0, cost_usd=0.08)
    out = VERD(first, second)
    assert out["verdict"] == "learned"


def test_single_proxy_worse_is_regressed():
    """Only cost went up while everything else is equal -> regressed, because
    at least one proxy strictly worsened and no proxy improved."""
    first = _rec(cost_usd=0.12)
    second = _rec(cost_usd=0.18)
    out = VERD(first, second)
    assert out["verdict"] == "regressed"
    assert out["error"] is None


def test_qualitative_same_but_float_rounding_tie():
    """duration_s equal as floats (the runner records floats) is a tie -> the
    strict-inequality rule must not flip an equal float to 'learned'."""
    first = _rec(duration_s=20.5)
    second = _rec(duration_s=20.5)
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"


def test_unknown_extra_field_ignored_not_error():
    """An unrelated extra field is not a missing one; it must neither raise
    nor change the verdict (the metric reads the proxies it knows)."""
    second = _rec(tool_calls=4, walltime_h=3)
    out = VERD(_rec(), second)
    assert out["verdict"] == "learned"
    assert out["error"] is None
