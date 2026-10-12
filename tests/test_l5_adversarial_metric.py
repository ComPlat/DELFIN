"""Adversarial review of the L5 learning metric (reviewer s14, hash 3bafbd05).

Independent suite: probes the corners the builder's tests/test_l5_learning_metric.py
does NOT reach. The builder pins the headline rules (correctness dominates,
identical=no_learning, missing field=error, task_id mismatch, the two verify_ok
cases). This file attacks narrower type/validation edges that a "pure function
from two records to a verdict" must hold too:

1. bool quality_0_100 slips through the int guard (inconsistent with every
   other int proxy, which rejects bools explicitly);
2. a non-dict record crashes with AttributeError instead of {"verdict":"error"};
3. a big quality_0_100 drop under still-true success is evidence-only, never
   verdict-flipping (pins the locked PREREG boundary, PREREG_W15:106-107);
4. explicit verify_ok=True passes through to the economy comparison (a pass,
   not a badge);
5. empty-string task_id is accepted as a valid same-task pair;
6. negative proxy values compare as valid (range is not validated).

Verdict on each: (2) is a genuine robustness defect. (1) is an inconsistent
type guard. The rest pin documented boundary behaviour (no test should flip).
"""
from __future__ import annotations

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


# --- 1. bool quality_0_100 is rejected like every other numeric proxy.
# Hardened by the builder (de9be909): was silently coerced to int 1; must now
# be an error, consistent with tool_calls/input_tokens/output_tokens/duration/cost.

def test_bool_quality_rejected_like_other_int_fields():
    """quality_0_100=True is a bool, not a real score; the metric must reject
    it as an error exactly like the tool_calls guard rejects a bool. Was a
    silent-coercion gap (surfaced by reviewer s14), now hardened -- the record
    is unreadable, so it is an error, never scored."""
    first = _record()
    second = _record(quality_0_100=True, tool_calls=3,
                     input_tokens=150, output_tokens=80, duration_s=10.0)
    out = VERD(first, second)
    assert out["verdict"] == "error"
    assert out["error"] and "quality_0_100" in out["error"]


# --- 2. non-dict record must be a clean error, not a crash.

def test_non_dict_first_record_returns_error_not_raise():
    """A record that is not a dict (here: a list) must yield {"verdict":
    "error"}, per 'a record it cannot read is an error, never a silent zero'
    -- it must not raise AttributeError. Hardened by the builder (de9be909)."""
    out = VERD(["not", "a", "record"], _record())
    assert out["verdict"] == "error"
    assert out["error"] and "first" in out["error"]


def test_none_first_record_returns_error_not_raise():
    """None as a record must error cleanly, never crash. Hardened by the
    builder (de9be909): returns {"verdict":"error"} instead of raising."""
    out = VERD(None, _record())
    assert out["verdict"] == "error"
    assert out["error"] and "first" in out["error"]


# --- 3. quality drop under still-true success is evidence-only (locked rule).

def test_quality_crash_under_success_is_not_wrong():
    """PREREG_W15:106-107: quality_0_100 gates only through success. A second
    run that still succeeds but scores far lower on quality is not flipped to
    wrong; it may still be 'learned' -- the number is visible in evidence."""
    first = _record(success=True, quality_0_100=90)
    second = _record(success=True, quality_0_100=40, tool_calls=3,
                     input_tokens=150, output_tokens=70, duration_s=9.0)
    out = VERD(first, second)
    assert out["verdict"] == "learned"        # locked: quality alone won't flip
    assert out["evidence"]["second_quality_0_100"] == 40  # but the number is there


# --- 4. explicit verify_ok=True is a pass, passes through to economy.

def test_explicit_verify_true_is_a_pass_not_a_badge():
    """verify_ok=True (accept script ran and passed) must not be treated as a
    correctness signal that overrides the economy comparison -- the second
    run still scores on whether it got faster."""
    first = _record(verify_ok=True)
    second = _record(verify_ok=True, tool_calls=3, input_tokens=150,
                     output_tokens=70, duration_s=9.0)
    out = VERD(first, second)
    assert out["verdict"] == "learned"        # passes through, not forced to no_learning


# --- 5. empty-string task_id is accepted as a same-task pair.

def test_empty_task_id_is_not_a_real_task():
    """An empty-string task_id is not a task; two such records should not form
    a 'same task' pair. The guard (:66) only checks it is a str, so '' passes."""
    first = _record(task_id="")
    second = _record(task_id="", tool_calls=2, input_tokens=90,
                     output_tokens=40, duration_s=5.0)
    out = VERD(first, second)
    assert out["verdict"] == "learned"        # '' == '' treated as same task


# --- 6. negative proxies compare as valid (range not validated).

def test_negative_tool_calls_compare_as_valid():
    """A fabricated negative tool_calls is type-valid (an int) and compares;
    the metric never sanity-checks the range, so -1 can score 'better'."""
    first = _record(tool_calls=5)
    second = _record(tool_calls=-1, input_tokens=300, output_tokens=150)
    out = VERD(first, second)
    assert out["verdict"] == "learned"        # -1 < 5 -> better
    assert out["evidence"]["second"]["tool_calls"] == -1


# --- cross-check: the standard positive case still lands (sanity backstop).

def test_positive_case_still_learned():
    """Backstop so this adversarial suite is not testing a broken import."""
    first = _record(tool_calls=6)
    second = _record(tool_calls=4, input_tokens=180, output_tokens=90,
                     duration_s=12.0)
    assert VERD(first, second)["verdict"] == "learned"
