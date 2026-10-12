"""L5 F1 (operator finding): wall-clock cost must not decide the verdict.

The metric decides learning by how many counted steps a second run cost. Wall
clock (``duration_s``) and its price proxy (``cost_usd``) are never equal
between two real runs, so if they decided, ``no_learning`` would be
unreachable and two runs that used the same counted ``tool_calls`` /
``tokens_total`` would come out ``learned`` or ``regressed`` on timing noise
alone.

Fixed rule (PREREG_W15.md dated amendment): the verdict is decided ONLY by the
counted proxies ``tool_calls`` and ``tokens_total``. ``duration_s`` and
``cost_usd`` are still validated (a record without them is an error) and still
reported in the evidence, but they never decide ``learned`` / ``regressed`` /
``no_learning``.

Red before the fix: duration.s on both lines still in the deciding set and
answers ``learned``/``regressed`` to these inputs.
"""

from delfin.agent.pack.benchmark import learning_metric as lm

VERD = lm.learning_verdict


def _record(*, task_id="l5_duration_cutoff", success=True,
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


def test_same_counted_fields_different_duration_is_no_learning():
    """Same tool_calls and same tokens_total, second 0.01 s faster -> NOT
    learned: the counted work is identical, so nothing measurable was saved.
    Wall-clock noise decides nothing (operator F1 case 1)."""
    first = _record(duration_s=20.00, cost_usd=0.120)
    second = _record(duration_s=19.99, cost_usd=0.119)   # 0.01 s / 0.001$ lower
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"


def test_same_counted_fields_different_duration_higher_is_no_learning():
    """Same counted fields, second 0.01 s SLOWER -> not regressed either: the
    counted work is identical, so a wall-clock wobble cannot count as a
    regression."""
    first = _record(duration_s=20.00, cost_usd=0.120)
    second = _record(duration_s=20.01, cost_usd=0.121)
    out = VERD(first, second)
    assert out["verdict"] == "no_learning"


def test_fewer_tool_calls_wins_over_higher_duration():
    """Second run used FEWER tool_calls (6 -> 4) but ran longer (20.0 -> 25.0):
    the counted work shrank, so it is learned. The longer wall clock must NOT
    flip it to regressed (operator F1 case 2)."""
    first = _record(tool_calls=6, input_tokens=300, output_tokens=150,
                    duration_s=20.0)
    second = _record(tool_calls=4, input_tokens=260, output_tokens=130,
                     duration_s=25.0)   # fewer tool_calls, higher duration
    out = VERD(first, second)
    assert out["verdict"] == "learned"


def test_identical_counted_fields_still_no_learning():
    """Exactly identical records (all fields equal) stay no_learning under the
    new rule -- the sanity invariant preserved."""
    out = VERD(_record(), _record())
    assert out["verdict"] == "no_learning"


def test_evidence_still_carries_duration_and_cost():
    """duration_s and cost_usd stay VALIDATED and REPORTED in the evidence,
    they just do not decide. A learned verdict's evidence still shows the
    second run's cost fields."""
    first = _record(tool_calls=6, duration_s=20.0, cost_usd=0.12)
    second = _record(tool_calls=4, duration_s=23.0, cost_usd=0.15)
    out = VERD(first, second)
    assert out["verdict"] == "learned"
    ev = out["evidence"]
    assert ev["first"]["duration_s"] == 20.0
    assert ev["second"]["duration_s"] == 23.0
    assert ev["second"]["cost_usd"] == 0.15
