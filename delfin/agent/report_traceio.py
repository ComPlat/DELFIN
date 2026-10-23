"""How long does reading the tool traces take?

Since ``delfin.agent.tool_trace`` grew ``root=`` (``read(session, root=...)``,
``sessions(root=...)``) a reader no longer has to walk the directory and parse
the JSONL itself -- but before that, every reader did exactly that. This
report measures whether the difference matters, on the real store, both ways:

* the module's readers: ``sessions(root=...)`` + ``read(sid, root=...)``
* a plain ``open()`` + ``json.loads`` walk of the same ``*.jsonl`` files

The measurement is repeated (default 5 runs per way) and reported as median
and spread -- one run of each answers nothing. The report carries sizes,
counts and timings ONLY: no trace content ever enters it.

Usage: ``python -m delfin.agent.report_traceio [root]``
"""

from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path
from typing import Any

from . import tool_trace

REPEATS = 5


def _plain_walk(root: Path) -> list[int]:
    """Pre-``root=`` style reader: glob + open() + json.loads, no helpers.

    Returns the per-file entry counts (the content itself is never kept).
    """
    counts: list[int] = []
    for f in sorted(root.glob("*.jsonl")):
        n = 0
        try:
            with open(f, encoding="utf-8") as fh:
                for line in fh:
                    try:
                        json.loads(line)
                        n += 1
                    except ValueError:
                        continue
        except OSError:
            continue
        counts.append(n)
    return counts


def _timed(fn: Any, *, repeats: int) -> dict[str, float]:
    """Run ``fn`` ``repeats`` times; median / min / max wall time in ms."""
    times: list[float] = []
    for _ in range(int(repeats)):
        t0 = time.perf_counter()
        fn()
        times.append((time.perf_counter() - t0) * 1000.0)
    return {
        "median_ms": statistics.median(times),
        "min_ms": min(times),
        "max_ms": max(times),
        "spread_ms": max(times) - min(times),
        "runs": float(len(times)),
    }


def collect(*, root: "Path | str" = "", repeats: int = REPEATS) -> dict:
    """Measure both readers over the trace store at ``root``.

    Never raises on a missing store: an empty dict-shaped result with
    ``n_files=0`` is returned, and the timings are zero (nothing was read).
    """
    directory = Path(root or tool_trace._DIR)
    files = sorted(directory.glob("*.jsonl")) if directory.is_dir() else []
    n_entries = 0
    total_bytes = 0
    for f in files:
        try:
            total_bytes += f.stat().st_size
        except OSError:
            pass

    def module_read() -> None:
        nonlocal n_entries
        n = 0
        for sid in tool_trace.sessions(root=directory):
            n += len(tool_trace.read(sid, root=directory))
        n_entries = n

    plain_counts: list[int] = []

    def plain_read() -> None:
        nonlocal plain_counts
        plain_counts = _plain_walk(directory)

    module_t = _timed(module_read, repeats=repeats) if files else {}
    plain_t = _timed(plain_read, repeats=repeats) if files else {}
    n_entries = n_entries or sum(plain_counts)

    return {
        "root": str(directory),
        "n_files": len(files),
        "n_entries": n_entries,
        "total_mb": total_bytes / (1024.0 * 1024.0),
        "repeats": int(repeats),
        "module_read": module_t,
        "plain_walk": plain_t,
    }


def format_text(data: dict) -> str:
    """Render collect()'s output as a compact timing report."""
    if not data.get("n_files"):
        return f"No trace files under {data.get('root', '?')} -- nothing to measure."
    lines = [
        f"Tool-trace read timing over {data['root']}",
        (f"files: {data['n_files']}  entries: {data['n_entries']}  "
         f"size: {data['total_mb']:.2f} MiB  runs each: {data['repeats']}"),
        "",
        f"{'reader':<24}{'median':>12}{'min':>12}{'max':>12}{'spread':>12}  (ms)",
    ]
    lines.append("-" * len(lines[-1]))
    for label, key in (("module (root=)", "module_read"),
                       ("plain open()+loads", "plain_walk")):
        t = data.get(key) or {}
        if not t:
            lines.append(f"{label:<24}{'-':>12}")
            continue
        lines.append(
            f"{label:<24}{t['median_ms']:>12.1f}{t['min_ms']:>12.1f}"
            f"{t['max_ms']:>12.1f}{t['spread_ms']:>12.1f}"
        )
    mod = (data.get("module_read") or {}).get("median_ms")
    pln = (data.get("plain_walk") or {}).get("median_ms")
    if mod is not None and pln:
        lines.append("")
        lines.append(f"module/plain ratio: {mod / pln:.2f}x")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    root = args[0] if args else ""
    try:
        data = collect(root=root)
    except Exception as exc:  # never raise from the CLI
        print(f"report_traceio failed: {exc}", file=sys.stderr)
        return 1
    print(format_text(data))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
