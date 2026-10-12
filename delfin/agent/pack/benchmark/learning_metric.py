"""learning_metric: does a second run of the same task do better than the first?

Wave 15a built a learning loop (skill use counted with its outcome, a
correction proposes a skill patch, a user model, cross-session search), but
nothing measured whether a *second* run of the same benchmark task is better
than the first. This module is that instrument: a single pure function from
two run records of the same task to a verdict.

A run record is one JSONL row as written by ``write_run``
(``delfin/agent/benchmark.py:2099``) -- ``asdict(BenchmarkResult)``. The
record has NO "rounds" and no "questions" counter; the honest per-run
efficiency proxies are ``tool_calls``, ``input_tokens``, ``output_tokens``
and ``duration_s``, and the correctness proxies are ``success``,
``quality_0_100`` and ``verify_ok``. The metric never invents a missing
field: a record it cannot read is an error, never a silent zero.

Rules fixed before any run (PREREG_W15.md):

* **Correctness dominates.** A second run that is faster but wrong is a fail
  whatever it saved (``wrong``).
* **Identical runs are not learning.** Equality on every efficiency proxy is
  ``no_learning``, never ``learned``.
* **Missing/unknown field is an error**, never a silent zero.

The CALLER decides which record is first and which is second: two rows from
one run file are replicates of one snapshot and are not a learning
comparison; first-vs-second requires two runs separated in time (two run
files). The metric is time-agnostic and pure -- it cannot see the caller's
intent, so that pairing rule lives in the caller (the runner hook).
"""

from __future__ import annotations

from typing import Any, Optional

# The efficiency proxies compared between the two runs; fewer is better.
# ``tokens_total`` is derived from the per-run ``input_tokens`` +
# ``output_tokens``, which is why it is not itself a record field.
_EFFICIENCY_FIELDS: tuple[str, ...] = ("tool_calls", "tokens_total",
                                       "duration_s", "cost_usd")

VERDICTS: tuple[str, ...] = ("error", "wrong", "learned", "no_learning",
                             "regressed")


def _get_int(record: dict, field: str) -> Optional[int]:
    """Read an integer-typed per-run field, or None if absent/mistyped.

    A field the record does not carry, or carries as a non-int, is treated
    as unreadable -- the caller asked for a comparison the record cannot
    answer, so the honest result is an error, not a fabricated zero.
    """
    value = record.get(field)
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return value


def _validate(record: dict, which: str) -> Optional[str]:
    """Validate one run record; return an error string or None if valid."""
    if not isinstance(record, dict):
        return f"{which} record is not a dict"
    unmeasured = record.get("unmeasured")
    if unmeasured is True:
        return (f"{which} run is unmeasured (the request never reached the "
                "model); a non-run is not a measurement")
    if "task_id" not in record or not isinstance(record.get("task_id"), str):
        return f"{which} record has no string task_id"
    if not isinstance(record.get("success"), bool):
        return f"{which} record has no bool success"
    q = record.get("quality_0_100")
    if isinstance(q, bool) or not isinstance(q, int):
        return f"{which} record has no int quality_0_100"
    for field in ("tool_calls", "input_tokens", "output_tokens"):
        if _get_int(record, field) is None:
            return (f"{which} record has no int {field} "
                    "(missing field is an error, never a silent zero)")
    dur = record.get("duration_s")
    if isinstance(dur, bool) or not isinstance(dur, (int, float)):
        return f"{which} record has no numeric duration_s"
    cost = record.get("cost_usd")
    if isinstance(cost, bool) or not isinstance(cost, (int, float)):
        return f"{which} record has no numeric cost_usd"
    return None


def learning_verdict(first: dict, second: dict) -> dict:
    """Compare two run records of the same task; how did the second do?

    Returns ``{"verdict": str, "evidence": dict, "error": str|None}``.
    ``verdict`` is one of ``("error", "wrong", "learned", "no_learning",
    "regressed")``. ``evidence`` carries first/second per efficiency proxy
    plus the second-run correctness, so a verdict can be quoted with its
    numbers. ``error`` is set (and ``verdict`` is ``"error"``) when the two
    records cannot be compared; never a silent zero.

    Correctness dominates: any second run that is not correct is ``wrong``
    whatever it saved. Equality on every efficiency proxy is
    ``no_learning``. See PREREG_W15.md for the fixed rules.
    """
    if not isinstance(first, dict):
        return {"verdict": "error", "evidence": {},
                "error": "first record is not a dict "
                         "(missing field is an error, never a silent zero)"}
    if not isinstance(second, dict):
        return {"verdict": "error", "evidence": {},
                "error": "second record is not a dict "
                         "(missing field is an error, never a silent zero)"}
    err_first = _validate(first, "first")
    err_second = _validate(second, "second")
    err = err_first or err_second
    if err is None and first["task_id"] != second["task_id"]:
        err = (f"task_id mismatch: {first['task_id']!r} != "
               f"{second['task_id']!r} -- a learning comparison needs two "
               "runs of the SAME task")
    if err is not None:
        return {"verdict": "error", "evidence": {}, "error": err}

    # Efficiency proxies, derived from the record.
    proxies_first = {
        "tool_calls": _get_int(first, "tool_calls"),
        "tokens_total": (_get_int(first, "input_tokens")
                         + _get_int(first, "output_tokens")),
        "duration_s": first["duration_s"],
        "cost_usd": first["cost_usd"],
    }
    proxies_second = {
        "tool_calls": _get_int(second, "tool_calls"),
        "tokens_total": (_get_int(second, "input_tokens")
                         + _get_int(second, "output_tokens")),
        "duration_s": second["duration_s"],
        "cost_usd": second["cost_usd"],
    }
    evidence = {
        "first": dict(proxies_first),
        "second": dict(proxies_second),
        "second_success": bool(second["success"]),
        "second_quality_0_100": int(second["quality_0_100"]),
    }

    # Correctness first. A wrong second answer fails whatever it saved.
    if not second["success"]:
        evidence["reason"] = "second run is not correct"
        return {"verdict": "wrong", "evidence": evidence, "error": None}
    # A second run whose artifact failed verification is wrong, whatever
    # ``success`` says. ``verify_ok is None`` (no verify script ran) is NOT a
    # correctness signal -- never scored as a pass, never as a fail: there
    # is no evidence a regression was proven, which is "no regression
    # proven", not "correct".
    if second.get("verify_ok") is False:
        evidence["reason"] = "second run did not verify (verify_ok False)"
        return {"verdict": "wrong", "evidence": evidence, "error": None}

    # Efficiency comparison across the efficiency proxies.
    better = []
    worse = []
    for field in _EFFICIENCY_FIELDS:
        a = proxies_first[field]
        b = proxies_second[field]
        if b < a:
            better.append(field)
        elif b > a:
            worse.append(field)

    if worse and not better:
        case = "regressed"
    elif better and worse:
        case = "regressed"
    elif better:
        case = "learned"
    else:
        case = "no_learning"
    evidence["reason"] = (
        f"better={better or 'none'}, worse={worse or 'none'}")
    return {"verdict": case, "evidence": evidence, "error": None}
