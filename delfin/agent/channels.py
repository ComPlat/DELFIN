"""Counters that prove a mechanism under test actually fired.

Input: ``fire(name)`` at the site where a mechanism acts (a line
injected, a claim flagged, a command admitted). Output: ``snapshot()``, a
``{name: count}`` dict, read by the benchmark runner after each attempt and
stored with its result (``BenchmarkResult.channels``). ``reset()`` at the
start of an attempt.

Semantics: an attempt that finished says nothing about whether the path
under test was reached. A comparison of two arms is read together with
these counts, so "the mechanism did not help" and "the mechanism never
ran" are told apart.

Process-wide and thread-safe (tool calls run on worker threads). A
benchmark attempt runs one engine at a time in this process, which is
what makes a process-wide counter per-attempt. Every security event
recorded through ``api_client._record_security_event`` is counted too,
under ``security:<kind>``, so the gate's admissions and refusals need no
second call site.

Rejected: a field on the engine. The mechanisms live in api_client
module functions, the gate, the engine and the sub-agent runner; not all
of them can reach an engine object, and a counter that only some of them
could reach would report silence as zero.
"""

from __future__ import annotations

import threading
from collections import Counter

_LOCK = threading.Lock()
_COUNTS: Counter = Counter()


def fire(name: str, n: int = 1) -> None:
    """Count one firing of ``name``. Never raises."""
    try:
        with _LOCK:
            _COUNTS[str(name)] += int(n)
    except Exception:
        pass


def snapshot() -> dict[str, int]:
    """The counts since the last reset, as a plain dict."""
    with _LOCK:
        return {k: int(v) for k, v in sorted(_COUNTS.items()) if v}


def reset() -> None:
    with _LOCK:
        _COUNTS.clear()


__all__ = ["fire", "snapshot", "reset"]
