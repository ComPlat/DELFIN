"""Adversarial red test (reviewer s14): the timing-noise defect (operator block).

The operator blocked L5: the verdict loop required ALL FOUR efficiency proxies
(tool_calls, tokens_total, duration_s, cost_usd) to be strictly equal for
``no_learning``, but ``duration_s`` is a float wall-clock that never round-trips
bit-identical across two real runs -- so ``no_learning`` is unreachable and the
verdict is dominated by duration/cost noise. Operator-required fix: the verdict
is decided by ``tool_calls`` + ``tokens_total`` ONLY.

This test is RED on the unfixed code (learning_metric.py:39-40, :154-172) and
must go GREEN once the fix lands: two runs equal on tool_calls and tokens_total,
differing only in duration_s (and cost_usd), must be ``no_learning`` -- NOT
``learned``/``regressed`` on a timing artefact.
"""
from __future__ import annotations

from delfin.agent.pack.benchmark import learning_metric as lm

VERD = lm.learning_verdict


def _row(*, task_id="l5_timing", success=True, quality_0_100=90, tool_calls=6,
         input_tokens=300, output_tokens=150, duration_s=20.0, cost_usd=0.12,
         verify_ok=True):
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
    }


def test_equal_effort_no_learning_not_driven_by_duration_noise():
    """tool_calls and tokens_total equal across the pair, duration_s differs
    by normde float noise: the verdict must be no_learning (effort did not
    change), never learned/regressed on a wall-clock artefact."""
    first = _row(duration_s=20.0, cost_usd=0.12)
    second = _row(duration_s=20.000037, cost_usd=0.120004)  # float jitter
    out = VERD(first, second)
    assert out["verdict"] == "no_learning", (
        f"got {out['verdict']!r}; equal effort (tool_calls+tokens) with "
        f"duration/cost float jitter must be no_learning, not a timing/cost "
        f"artefact")


def test_equal_effort_but_duration_regressed_is_not_regressed():
    """Duration alone must not flip the verdict. Same effort (tool_calls,
    tokens), second run slower on wall-clock only -> no_learning, not
    regressed."""
    first = _row(duration_s=20.0)
    second = _row(duration_s=25.0)  # slower wall-clock, same effort
    out = VERD(first, second)
    assert out["verdict"] == "no_learning", (
        f"got {out['verdict']!r}; duration-only difference must not count as "
        f"a regression when tool_calls+tokens_total are unchanged")


def test_equal_effort_but_duration_better_is_not_learned():
    """Duration alone must not flip the verdict. Same effort, second run faster
    on wall-clock only -> no_learning, not learned."""
    first = _row(duration_s=20.0)
    second = _row(duration_s=15.0)  # faster wall-clock, same effort
    out = VERD(first, second)
    assert out["verdict"] == "no_learning", (
        f"got {out['verdict']!r}; duration-only gain must not count as "
        f"learning when tool_calls+tokens_total are unchanged")
