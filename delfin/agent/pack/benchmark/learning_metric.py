"""learning_metric: does a second run of the same task do better than the first?

Input: two run records of the same task, each one JSONL row as written by
``benchmark.write_run`` (``asdict(BenchmarkResult)``), the first made before
the second under one learning home (``bench run --learning-home``). Output:
a verdict with its evidence. Pure: reads its two dicts, mutates neither.

Rules (PREREG_W15.md):

* **Correctness dominates.** A second run that is not correct
  (``success`` False, or ``verify_ok`` False) is ``wrong`` whatever it
  saved.
* **N=1 is not evidence.** Each record must carry at least
  ``MIN_REPLICATES`` replicates with their per-replicate counts
  (``per_run_tool_calls``, ``per_run_tokens_total``); fewer is ``error``.
* **A difference must lie beyond the spread.** A counted field is
  ``better`` only when every second-arm replicate is below every
  first-arm replicate, ``worse`` only when every one is above (complete
  separation; a shared value is not a separation). ``learned``: at least
  one field better, none worse. ``regressed``: any field worse.
  ``no_learning``: neither field separates. Under no effect, with
  exchangeable samples and no ties, a given field separates in a given
  direction with probability 1 / C(2N, N): 1/252 at N = 5.
* **Only counted work decides.** ``tool_calls`` and ``tokens_total``
  (``input_tokens + output_tokens``) decide; ``duration_s`` and
  ``cost_usd`` are validated and reported but never decide.
* **A missing or mistyped field is an error**, never a silent zero.

Measured cost of the separation rule against strict inequality on the
medians, over 65 disjoint 5-versus-5 windows of successful rows of one
task and model from existing run files (no learning setup): medians gave
39 learned, 26 regressed, 0 no_learning; separation gave 1, 1, 63.
Rejected: a relative band on the medians (a free parameter with no
calibration behind it; across tasks the median coefficient of variation
is 0.33 for tokens and 0.37 for tool calls, with a wide spread, so any one
band is too wide for some tasks and too narrow for others).

Pairing (same model, profile, mode, learning home) is checked by the
caller, ``benchmark.compare_learning``; this function checks only that
both records are of the same task.
"""

from __future__ import annotations

from typing import Any, Optional

#: Replicates per arm below which no verdict but ``error`` is given.
MIN_REPLICATES = 5

#: The counted proxies the verdict is decided on (fewer is better), each
#: with the record field that holds its per-replicate samples.
_DECIDING_FIELDS: tuple[tuple[str, str], ...] = (
    ("tool_calls", "per_run_tool_calls"),
    ("tokens_total", "per_run_tokens_total"),
)

VERDICTS: tuple[str, ...] = ("error", "wrong", "learned", "no_learning",
                             "regressed")


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _get_int(record: dict, field: str) -> Optional[int]:
    """An integer field, or None if absent or mistyped (a bool is not an
    int here)."""
    value = record.get(field)
    return value if _is_int(value) else None


def _samples(record: dict, field: str) -> Optional[list[int]]:
    """The per-replicate samples of ``field``: a list of ints whose length
    equals ``n_samples``; None otherwise."""
    values = record.get(field)
    n = record.get("n_samples")
    if not isinstance(values, list) or not _is_int(n) or len(values) != n:
        return None
    if not all(_is_int(v) for v in values):
        return None
    return list(values)


def _validate(record: Any, which: str) -> Optional[str]:
    """Validate one run record; return an error string or None if valid."""
    if not isinstance(record, dict):
        return (f"{which} record is not a dict "
                "(missing field is an error, never a silent zero)")
    if record.get("unmeasured") is True:
        return (f"{which} run is unmeasured (the request never reached the "
                "model); a non-run is not a measurement")
    if not isinstance(record.get("task_id"), str):
        return f"{which} record has no string task_id"
    if not isinstance(record.get("success"), bool):
        return f"{which} record has no bool success"
    if _get_int(record, "quality_0_100") is None:
        return f"{which} record has no int quality_0_100"
    for field in ("tool_calls", "input_tokens", "output_tokens"):
        if _get_int(record, field) is None:
            return (f"{which} record has no int {field} "
                    "(missing field is an error, never a silent zero)")
    for field in ("duration_s", "cost_usd"):
        value = record.get(field)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return f"{which} record has no numeric {field}"
    n = record.get("n_samples")
    if not _is_int(n) or n < MIN_REPLICATES:
        return (f"{which} run has N={n} replicates; a learning verdict "
                f"needs at least {MIN_REPLICATES} per arm (N=1 is not "
                "evidence)")
    for _, samples_field in _DECIDING_FIELDS:
        if _samples(record, samples_field) is None:
            return (f"{which} record has no {samples_field} of {n} ints "
                    "(missing field is an error, never a silent zero)")
    return None


def _proxies(record: dict) -> dict:
    out: dict[str, Any] = {
        "tool_calls": record["tool_calls"],
        "tokens_total": record["input_tokens"] + record["output_tokens"],
        "duration_s": record["duration_s"],
        "cost_usd": record["cost_usd"],
    }
    for field, samples_field in _DECIDING_FIELDS:
        out[f"{field}_samples"] = _samples(record, samples_field)
    return out


def learning_verdict(first: dict, second: dict) -> dict:
    """Compare two run records of the same task; how did the second do?

    Returns ``{"verdict": str, "evidence": dict, "error": str|None}`` with
    ``verdict`` one of :data:`VERDICTS`. ``evidence`` carries, per arm, the
    median counts, the per-replicate samples, duration and cost, plus
    ``n_first`` / ``n_second`` and the second run's correctness. ``error``
    is set (and ``verdict`` is ``"error"``) when the records cannot be
    compared.
    """
    for which, record in (("first", first), ("second", second)):
        err = _validate(record, which)
        if err is not None:
            return {"verdict": "error", "evidence": {}, "error": err}
    if first["task_id"] != second["task_id"]:
        return {"verdict": "error", "evidence": {},
                "error": (f"task_id mismatch: {first['task_id']!r} != "
                          f"{second['task_id']!r} -- a learning comparison "
                          "needs two runs of the SAME task")}

    p_first, p_second = _proxies(first), _proxies(second)
    evidence: dict[str, Any] = {
        "first": p_first,
        "second": p_second,
        "n_first": first["n_samples"],
        "n_second": second["n_samples"],
        "second_success": bool(second["success"]),
        "second_quality_0_100": int(second["quality_0_100"]),
    }

    if not second["success"]:
        evidence["reason"] = "second run is not correct"
        return {"verdict": "wrong", "evidence": evidence, "error": None}
    # ``verify_ok is None`` (no acceptance script ran) is neither a pass
    # nor a fail.
    if second.get("verify_ok") is False:
        evidence["reason"] = "second run did not verify (verify_ok False)"
        return {"verdict": "wrong", "evidence": evidence, "error": None}

    better: list[str] = []
    worse: list[str] = []
    for field, _ in _DECIDING_FIELDS:
        a = p_first[f"{field}_samples"]
        b = p_second[f"{field}_samples"]
        if max(b) < min(a):
            better.append(field)
        elif min(b) > max(a):
            worse.append(field)

    if worse:
        case = "regressed"
    elif better:
        case = "learned"
    else:
        case = "no_learning"
    evidence["reason"] = (
        f"separated better={better or 'none'}, worse={worse or 'none'}")
    return {"verdict": case, "evidence": evidence, "error": None}
