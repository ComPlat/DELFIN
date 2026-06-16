"""Run metrics — aggregate the run store into timing / success / cost signals.

Reads the persisted :class:`~delfin.tools._runtime.RunRecord` history and rolls
it up per status and per application, so a human or agent can see how the
platform is performing and route accordingly.
"""

from __future__ import annotations

from typing import Any, Dict, Optional


def aggregate_runs(store=None) -> Dict[str, Any]:
    """Aggregate all persisted runs into a metrics summary."""
    from delfin.tools._runtime import RunStore

    store = store or RunStore()
    records = store.list()

    by_status: Dict[str, int] = {}
    by_application: Dict[str, Dict[str, Any]] = {}
    total_elapsed = 0.0

    for r in records:
        by_status[r.status] = by_status.get(r.status, 0) + 1
        elapsed = float(r.metrics.get("elapsed_s", 0) or 0)
        total_elapsed += elapsed
        app = by_application.setdefault(
            r.name, {"count": 0, "success": 0, "elapsed_s": 0.0}
        )
        app["count"] += 1
        app["success"] += 1 if r.status == "success" else 0
        app["elapsed_s"] += elapsed

    for app in by_application.values():
        n = app["count"]
        app["success_rate"] = round(app["success"] / n, 3) if n else 0.0
        app["avg_elapsed_s"] = round(app["elapsed_s"] / n, 3) if n else 0.0
        app["elapsed_s"] = round(app["elapsed_s"], 3)

    return {
        "total_runs": len(records),
        "by_status": by_status,
        "total_elapsed_s": round(total_elapsed, 3),
        "by_application": by_application,
    }


__all__ = ["aggregate_runs"]
