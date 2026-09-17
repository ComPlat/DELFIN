"""Which tests are not telling the truth? A flake report.

A flake is a test that was green in one run and red in another with no
code change between them. The evidence already exists — nothing new is
recorded:

* ``delfin.agent.tool_trace`` (durable, on disk): every ``run_tests``
  call with its full JSON result, including the per-node ``failures``
  list, plus every file-mutating tool call (write/edit/patch) with its
  timestamp. Read through ``tool_trace.sessions()`` /
  ``tool_trace.read(..., root=...)`` — this module adds no store.

Two honest limits, both surfaced by the report itself:

* ``run_tests`` records only the FAILING node ids, not the full
  collected set. A node is therefore only a flake candidate once it has
  failed somewhere; a test that is silently green twice and red once is
  found, but one that is red once and never seen again is not paired.
* Only runs with a durable verdict count (``status`` in ok/failed).
  A ``stopped`` or ``timeout`` run says nothing about the suite (the
  runner says so itself in its ``note``), and a trace entry whose
  recorded output was truncated past readability (the trace caps at
  2000 chars; ~half of the real recorded run_tests outputs are cut like
  that) is counted as unreadable, never guessed from.

A code change is, conservatively, ANY recorded mutating tool call
(write_file / edit_file / multi_edit / apply_patch / notebook_edit /
create-* of files) between two runs of the same target — a pair with a
write anywhere between them is not evidence of a flake, even when the
write touched an unrelated file.

CLI: ``python -m delfin.agent.report_flakes``
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from . import tool_trace

# run_tests results that carry a verdict about the suite's state.
_VERDICT_STATUSES = frozenset({"ok", "failed"})

# Tool names (substring-matched, prefixes like mcp__kit-coding__ included)
# that mutate files. Any one of these between two runs disqualifies the
# pair from being flake evidence.
_MUTATING = ("write_file", "edit_file", "multi_edit", "apply_patch",
             "notebook_edit")

_MAX_EXAMPLES = 20


def _is_run_tests(tool: str) -> bool:
    return "run_tests" in str(tool or "")


def _is_mutating(entry: dict) -> bool:
    tool = str(entry.get("tool", "") or "")
    return any(m in tool for m in _MUTATING)


def _result_of(entry: dict) -> dict | None:
    """The parsed run_tests result of a trace entry, or None.

    None when the output was truncated past parseability (the trace
    caps output at 2000 chars) or is not a JSON object — counted as
    unreadable, never guessed from.
    """
    out = entry.get("output")
    if not isinstance(out, str):
        return None
    try:
        data = json.loads(out)
    except ValueError:
        return None
    return data if isinstance(data, dict) else None


def collect(*, runs: list[dict] | None = None,
            root: str | Path = "") -> dict[str, Any]:
    """Scan recorded test runs and return the flake analysis.

    ``runs`` overrides the trace scan with pre-collected run records
    (each: ``ts``, ``target``, ``status``, ``failures`` node ids,
    optional ``session``) plus their sibling ``edits`` — the shape used
    by the tests and by any caller that already holds the records.
    ``root`` reads a trace directory other than this machine's own.
    """
    if runs is not None:
        run_recs = [r for r in runs if isinstance(r, dict)]
        # Edits ride along per-run under "_edits" in the fabricated
        # form (each: ts of one mutating tool call); flatten them.
        all_edits: list[dict] = []
        for r in run_recs:
            all_edits.extend(e for e in (r.get("_edits") or [])
                             if isinstance(e, dict))
        n_unreadable = 0
        n_no_verdict = 0
        sessions = set()
    else:
        run_recs = []
        all_edits = []
        n_unreadable = 0
        n_no_verdict = 0
        sessions = set()
        for sid in tool_trace.sessions(root=root):
            sessions.add(sid)
            for entry in tool_trace.read(sid, root=root):
                if not isinstance(entry, dict):
                    continue
                if _is_mutating(entry):
                    all_edits.append(entry)
                    continue
                if not _is_run_tests(entry.get("tool", "")):
                    continue
                data = _result_of(entry)
                if data is None:
                    n_unreadable += 1
                    continue
                status = str(data.get("status", "") or "")
                if status not in _VERDICT_STATUSES:
                    n_no_verdict += 1
                    continue
                args = tool_trace.call_args(entry)
                run_recs.append({
                    "ts": float(entry.get("ts", 0.0) or 0.0),
                    "session": sid,
                    "target": str(args.get("target", "") or ""),
                    "status": status,
                    "exit_code": data.get("exit_code"),
                    "failures": [
                        str(f.get("node_id", "") or "")
                        for f in (data.get("failures") or [])
                        if isinstance(f, dict) and f.get("node_id")
                    ],
                })

    edit_ts = sorted(
        float(e.get("ts", 0.0) or 0.0) for e in all_edits if e.get("ts"))

    def _edits_between(t_lo: float, t_hi: float) -> int:
        lo, hi = 0, len(edit_ts)
        while lo < hi:
            mid = (lo + hi) // 2
            if edit_ts[mid] < t_lo:
                lo = mid + 1
            else:
                hi = mid
        first = lo
        lo, hi = first, len(edit_ts)
        while lo < hi:
            mid = (lo + hi) // 2
            if edit_ts[mid] < t_hi:
                lo = mid + 1
            else:
                hi = mid
        return lo - first

    # Group verdict runs per target, oldest first.
    by_target: dict[str, list[dict]] = {}
    for r in run_recs:
        by_target.setdefault(r.get("target", ""), []).append(r)
    for lst in by_target.values():
        lst.sort(key=lambda r: float(r.get("ts", 0.0) or 0.0))

    flakes: list[dict[str, Any]] = []
    for target, lst in sorted(by_target.items()):
        if len(lst) < 2:
            continue
        # Every node that failed anywhere in this target's runs.
        candidates: set[str] = set()
        for r in lst:
            candidates.update(r.get("failures") or [])
        for i in range(len(lst) - 1):
            a, b = lst[i], lst[i + 1]
            t_a, t_b = float(a["ts"]), float(b["ts"])
            if t_a == t_b:
                continue
            if _edits_between(t_a, t_b):
                continue  # a code change sits between the two runs
            red_a = set(a.get("failures") or [])
            red_b = set(b.get("failures") or [])
            for node in sorted(candidates):
                if node in red_a and node not in red_b:
                    flakes.append({
                        "node_id": node, "target": target,
                        "was": "red", "then": "green",
                        "ts_red": t_a, "ts_green": t_b,
                    })
                elif node not in red_a and node in red_b:
                    flakes.append({
                        "node_id": node, "target": target,
                        "was": "green", "then": "red",
                        "ts_red": t_b, "ts_green": t_a,
                    })

    return {
        "flakes": flakes,
        "n_runs": len(run_recs),
        "n_targets": len(by_target),
        "n_traces": len(sessions) if runs is None else None,
        "n_unreadable_runs": n_unreadable,
        "n_no_verdict_runs": n_no_verdict,
    }


def format_text(data: dict) -> str:
    """Render collect()'s output as a compact report."""
    flakes = data.get("flakes") or []
    lines: list[str] = []
    lines.append(
        f"Flaky tests over {data.get('n_runs', 0)} recorded verdict run(s), "
        f"{data.get('n_targets', 0)} target(s):"
    )
    lines.append(
        f"  unreadable run records: {data.get('n_unreadable_runs', 0)} "
        "(trace output truncated past parseability)"
    )
    lines.append(
        f"  runs without a verdict (stopped/timeout/error): "
        f"{data.get('n_no_verdict_runs', 0)}"
    )
    if not flakes:
        lines.append("")
        lines.append(
            "No test flipped colour between two runs of the same target "
            "with no code change recorded between them.")
        lines.append(
            "(Caveat: only tests that FAILED somewhere are candidates — "
            "run_tests records failing node ids only.)")
        return "\n".join(lines)
    lines.append("")
    lines.append(f"{len(flakes)} flip(s):")
    lines.append(f"  {'node_id':<64}{'was':>6}{'then':>7}")
    lines.append("  " + "-" * 77)
    for f in flakes[:_MAX_EXAMPLES]:
        lines.append(
            f"  {f['node_id'][:64]:<64}{f['was']:>6}{f['then']:>7}"
            f"  ({f['target'] or 'whole suite'})")
    if len(flakes) > _MAX_EXAMPLES:
        lines.append(f"  … and {len(flakes) - _MAX_EXAMPLES} more")
    lines.append("")
    lines.append(
        "A flip is red-then-green or green-then-red across two runs of "
        "the same target with no mutating tool call between them.")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    data = collect()
    print(format_text(data))
    return 0


if __name__ == "__main__":
    sys.exit(main())
