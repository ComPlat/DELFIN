"""A learning verdict needs a difference beyond the run-to-run spread.

Measured on 65 disjoint 5-versus-5 windows of successful N=1 rows of the
same task and model from existing run files (no learning setup on any of
them): strict inequality on the medians of ``tool_calls`` and
``tokens_total`` gave 39 learned, 26 regressed, 0 no_learning. Complete
separation of the replicate samples gave 1 learned, 1 regressed, 63
no_learning.

The rule: each arm carries at least ``MIN_REPLICATES`` per-replicate
counts; a counted field is better only when every second-arm replicate is
below every first-arm replicate, worse only when every one is above. Under
no effect (exchangeable samples, no ties) the chance of a given field
separating in a given direction is 1 / C(10, 5) = 1/252 at N = 5 per arm.
"""
from __future__ import annotations

from itertools import combinations

from delfin.agent import benchmark as bm
from delfin.agent.pack.benchmark import learning_metric as lm


def _rec(tool_calls, tokens=None, **over):
    tokens = tokens if tokens is not None else [450] * len(tool_calls)
    rec = {"task_id": "l5_read_the_value_from_file", "success": True,
           "quality_0_100": 90, "tool_calls": sorted(tool_calls)[len(tool_calls) // 2],
           "input_tokens": 300, "output_tokens": 150, "duration_s": 20.0,
           "cost_usd": 0.12, "verify_ok": True,
           "n_samples": len(tool_calls),
           "per_run_tool_calls": list(tool_calls),
           "per_run_tokens_total": list(tokens)}
    rec.update(over)
    return rec


def test_aggregation_keeps_the_per_replicate_counts():
    reps = [bm.BenchmarkResult(task_id="t", task_class="c", model="m",
                               success=True, tool_calls=tc,
                               input_tokens=100 + tc, output_tokens=10)
            for tc in (3, 7, 5)]
    agg = bm.aggregate_replicates(reps)
    assert agg.per_run_tool_calls == [3, 7, 5]
    assert agg.per_run_tokens_total == [113, 117, 115]


def test_one_replicate_per_arm_is_not_evidence():
    out = lm.learning_verdict(_rec([6]), _rec([2]))
    assert out["verdict"] == "error"
    assert "N=1" in out["error"]


def test_fewer_replicates_than_the_minimum_is_not_evidence():
    out = lm.learning_verdict(_rec([6, 6, 6, 6]), _rec([2, 2, 2, 2]))
    assert out["verdict"] == "error"


def test_a_lower_median_inside_the_spread_is_no_learning():
    first = _rec([4, 6, 7, 9, 12])
    second = _rec([3, 5, 6, 8, 10])
    assert second["tool_calls"] < first["tool_calls"]
    assert lm.learning_verdict(first, second)["verdict"] == "no_learning"


def test_a_complete_separation_downwards_is_learned():
    out = lm.learning_verdict(_rec([7, 8, 9, 9, 12]), _rec([2, 3, 3, 4, 6]))
    assert out["verdict"] == "learned"


def test_a_complete_separation_upwards_is_regressed():
    out = lm.learning_verdict(_rec([2, 3, 3, 4, 6]), _rec([7, 8, 9, 9, 12]))
    assert out["verdict"] == "regressed"


def test_a_shared_value_at_the_boundary_is_not_a_separation():
    out = lm.learning_verdict(_rec([6, 7, 8, 9, 10]), _rec([2, 3, 4, 5, 6]))
    assert out["verdict"] == "no_learning"


def test_per_replicate_counts_that_do_not_match_n_are_an_error():
    bad = _rec([6, 6, 6, 6, 6], n_samples=6)
    assert lm.learning_verdict(bad, _rec([2] * 5))["verdict"] == "error"
    bad = _rec([6, 6, 6, 6, 6], per_run_tokens_total=[450] * 4)
    assert lm.learning_verdict(bad, _rec([2] * 5))["verdict"] == "error"
    bad = _rec([6, 6, 6, 6, True])
    assert lm.learning_verdict(bad, _rec([2] * 5))["verdict"] == "error"


def test_under_no_effect_one_split_in_252_reads_as_learned():
    values = list(range(1, 11))
    learned = regressed = 0
    for idx in combinations(range(10), 5):
        a = [values[i] for i in idx]
        b = [values[i] for i in range(10) if i not in idx]
        v = lm.learning_verdict(_rec(a), _rec(b))["verdict"]
        learned += v == "learned"
        regressed += v == "regressed"
    assert (learned, regressed) == (1, 1)


def test_the_evidence_carries_the_samples_and_n():
    out = lm.learning_verdict(_rec([7, 8, 9, 9, 12]), _rec([2, 3, 3, 4, 6]))
    ev = out["evidence"]
    assert ev["n_first"] == 5 and ev["n_second"] == 5
    assert ev["first"]["tool_calls_samples"] == [7, 8, 9, 9, 12]
    assert ev["second"]["tool_calls_samples"] == [2, 3, 3, 4, 6]
