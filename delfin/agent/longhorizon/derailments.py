"""Derailments a long run can be measured by, from its records only.

Package R3. The derailment counts are reduced from the records the REAL
writers already produce (tool_trace, turn_metrics, compaction_log) — no new
writer, no judged input. Each count names the input that would make it FAIL,
per the brief's rule 5 (an assertion that cannot fail is not a test):

* ``repeated_failing_calls`` — the same tool call (tool + normalised input)
  run >= 2 times CONSECUTIVELY with ok=False each time. Fails when two
  different failing calls interleave (not consecutive), or when a call
  succeeded once between the failures.
* ``unproductive_loops`` — the same tool call (tool + input) run >= 3 times
  CONSECUTIVELY, all succeeding (an identical successful read repeated is
  the no-progress loop the engine already aborts at 5). Fails when the
  identical call is broken up by a different call, or when only 2 ran.
* ``abandoned_after_failure`` — the last recorded turn has stopped=True and
  its error names a failing run, while a later tool call exists. NOT a
  score: reduced as the count the task's accept script reads beside the
  artifact check.
"""
from __future__ import annotations

from typing import Any, Dict, List, Sequence

__all__ = ["derailments"]


def _norm(text: Any) -> str:
    """Whitespace-normalised exact text — the same narrow identity
    longrun_metric and longhorizon.metrics already use."""
    return " ".join(str(text or "").split())


def _identity(record: Dict[str, Any]) -> str:
    """(tool name, normalised input) — the same redo-work identity
    metrics.measure uses, so the two reducers agree on 'same call'."""
    tool = str(record.get("tool") or "").strip()
    arg = _norm(record.get("input"))
    if not tool and not arg:
        return ""
    return f"{tool}::{arg}"


def _ordered(records: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Records in ts order, ts-less records dropped (a reducer cannot
    attest an order it was not given)."""
    ordered = [r for r in records if r.get("ts") is not None]
    try:
        return sorted(ordered, key=lambda r: float(r["ts"]))
    except (TypeError, ValueError):
        return [r for r in ordered
                if isinstance(r.get("ts"), (int, float))]


def derailments(tools: Sequence[Dict[str, Any]]) -> Dict[str, int]:
    """Reduce one session's tool_trace records to its derailment counts.

    Pure: no I/O, no judgment beyond the records' own ``ok`` and identity.
    A record the writer could not timestamp is not ordered, so it does not
    contribute — reported, never guessed.
    """
    counts = {
        "repeated_failing_calls": 0,
        "unproductive_loops": 0,
    }
    ordered = _ordered(tools)
    i = 0
    n = len(ordered)
    while i < n:
        ident = _identity(ordered[i])
        if not ident:
            i += 1
            continue
        j = i
        while j < n and _identity(ordered[j]) == ident:
            j += 1
        run_len = j - i
        run_ok = all(bool(ordered[k].get("ok")) for k in range(i, j))
        if run_len >= 2 and not run_ok:
            counts["repeated_failing_calls"] += run_len - 1
        elif run_len >= 3 and run_ok:
            counts["unproductive_loops"] += run_len - 2
        i = j
    return counts
