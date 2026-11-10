#!/usr/bin/env python3
"""lh_metrics_accept.py — grade one long-horizon run, accept or reject it.

Called by the DELFIN benchmark harness as::

    python lh_metrics_accept.py <workspace>

with exactly ONE argument (run_setup / run_acceptance contract,
benchmark_runner.py:1318/1350). Reads the workspace file
``.longhorizon_metrics.json`` — written by the F2 runner
(delfin/agent/longhorizon/runner.py) — grades the run's metrics against its
thresholds, prints a one-line verdict, and exits 0 on accept, 1 on reject.

This accept script is the single source of truth for "this long-horizon
run was good enough to adopt". A compaction change (F3, operator-B) may be
adopted only if it scores strictly better here (E2-holdout style): the
baseline's own acceptance run is the holdout that must not regress.

The JSON shape::

    {
      "task_id": str,
      "metrics": {
        "turns": int, "tool_calls": int, "compactions": int,
        "facts_recalled_late": int, "facts_total": int,
        "decisions_honoured_late": int, "decisions_total": int,
        "redone_work": int,
        "input_tokens": int, "output_tokens": int, "cached_tokens": int,
        "compactions_size": [int, ...], "total_ms": int
      },
      "thresholds": {
        "facts_total_min": int|None, "recall_late_min": float|None,
        "turns_min": int|None, "redone_work_max": int|None
      }
    }

``recall_late`` is computed = facts_recalled_late / facts_total (defensive:
0 recalled of 0 facts is a pass of that clause). Anything that is missing
or malformed is reported and counts as a reject — never silently passed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def _metrics_file(workspace: Path) -> Path:
    return Path(workspace) / ".longhorizon_metrics.json"


def _check(
    name: str,
    metric: int | float | None,
    op: str,
    bound: int | float | None,
    problems: list[str],
) -> None:
    """Append a violation to ``problems`` when ``metric`` breaks the bound."""
    if bound is None:
        return
    if metric is None:
        problems.append(f"{name} is missing (bound {op} {bound})")
        return
    if op == "min" and metric < bound:
        problems.append(f"{name}={metric} below min {bound}")
    elif op == "max" and metric > bound:
        problems.append(f"{name}={metric} above max {bound}")


def grade(workspace: Path) -> tuple[bool, str]:
    f = _metrics_file(workspace)
    if not f.is_file():
        return False, f"metrics file not found: {f}"
    try:
        data = json.loads(f.read_text("utf-8"))
    except (OSError, ValueError) as exc:
        return False, f"metrics file unreadable: {exc}"

    metrics = data.get("metrics")
    thresholds = data.get("thresholds") or {}
    if not isinstance(metrics, dict) or not isinstance(thresholds, dict):
        return False, "metrics or thresholds is not an object"

    recall_late = None
    facts_total = metrics.get("facts_total")
    facts_recalled = metrics.get("facts_recalled_late")
    if isinstance(facts_total, int) and isinstance(facts_recalled, int):
        if facts_total > 0:
            recall_late = facts_recalled / facts_total
        else:
            recall_late = 1.0

    problems: list[str] = []
    _check("recall_late", recall_late, "min", thresholds.get("recall_late_min"), problems)
    _check("turns", metrics.get("turns"), "min", thresholds.get("turns_min"), problems)
    _check("facts_total", facts_total, "min", thresholds.get("facts_total_min"), problems)
    _check("redone_work", metrics.get("redone_work"), "max", thresholds.get("redone_work_max"), problems)

    if problems:
        return False, "REJECT " + "; ".join(problems)
    return True, f"ACCEPT task={data.get('task_id')} recall_late={recall_late} turns={metrics.get('turns')}"


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: lh_metrics_accept.py <workspace>", file=sys.stderr)
        return 2
    ok, msg = grade(Path(sys.argv[1]))
    print(msg)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
