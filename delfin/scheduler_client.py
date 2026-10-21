"""One throttled, cached place to ask the SLURM scheduler.

Wave 11 finding 1: a scheduler query package would have called ``squeue``
every 5 s per idle terminal, and cluster operators had objected to
per-widget ``squeue`` load before. Scheduler queries are spread over
several modules, each with its own idea of how often to ask. This module
is the single address for those queries: a process-wide gate that lets a
*distinct query kind* reach the real scheduler at most once per
:data:`MIN_INTERVAL_S` seconds, and re-serves the last answer it got
instead of asking again while the gate is closed.

Two properties the cluster needs, both process-wide:

* **Throttle by kind.** ``squeue`` and ``sacct`` are different kinds; a
  watch loop that alternates them may ask each once per interval. The
  gate is per kind, so one busy component cannot starve another kind, and
  no kind is asked more often than the interval allows.
* **Shared cache.** The cache is process-global, keyed by
  ``(kind, argv)``. Two watchers asking the *same* question within the
  interval share one answer and one scheduler call. A fresh argv that the
  gate has closed returns :data:`THROTTLED` rather than inventing an
  answer — a caller must tell "throttled, not asked" apart from a real
  answer, exactly as ``job_monitor`` already distinguishes "not known"
  from "not asked".

The client never touches the scheduler itself: the lower-level runner is
injectable (``run_fn``), so tests drive it with fakes and the real
subprocess lives behind the default runner only.
"""
from __future__ import annotations

import os
import time
from typing import Callable, Optional, Sequence, Tuple

#: Minimum seconds between two *actual* queries of the same kind.
MIN_INTERVAL_S = 25.0

#: Sentinel for "the gate is closed; this query was NOT answered by the
#: scheduler". Distinct from ``None`` (the scheduler answered "no such
#: job" / "could not be run"), so a caller can count throttled asks.
THROTTLED = object()


def _default_runner(cmd: Sequence[str]) -> Optional[str]:
    """Run a scheduler query. ``None`` means it could not be run at all."""
    try:
        import subprocess
        out = subprocess.run(list(cmd), capture_output=True, text=True, timeout=20)
        return out.stdout if out.returncode == 0 else None
    except Exception:
        return None


class SchedulerClient:
    """Process-wide throttle + shared cache for scheduler queries.

    Not thread-locked on purpose: the scheduler is asked at most once per
    interval per kind, and a rare double-ask during a race is a harmless
    single extra query — the lock would buy little and complicate tests.
    """

    def __init__(
        self,
        min_interval_s: float = MIN_INTERVAL_S,
        clock: Optional[Callable[[], float]] = None,
        runner: Optional[Callable[[Sequence[str]], Optional[str]]] = None,
    ) -> None:
        self.min_interval_s = min_interval_s
        self._clock = clock or time.monotonic
        self._runner = runner or _default_runner
        #: key = (kind, argv) -> (last_run_ts, last_result)
        self._cache: dict[Tuple[str, Tuple[str, ...]], Tuple[float, object]] = {}
        #: kind -> last time an *actual* query of this kind ran
        self._kind_last_run: dict[str, float] = {}
        # stats for tests and diagnostics
        self.queries_run = 0
        self.served_from_cache = 0

    def _now(self) -> float:
        return self._clock()

    def query(
        self,
        kind: str,
        argv: Sequence[str],
        runner: Optional[Callable[[Sequence[str]], Optional[str]]] = None,
    ) -> object:
        """Return the answer for ``argv``, throttled per ``kind``.

        Returns the cached answer when the gate is open for this exact
        ``(kind, argv)``, :data:`THROTTLED` when the gate is closed and
        there is nothing cached for this argv, else runs the query (or
        the injected ``runner``) and caches the result.
        """
        key = (kind, tuple(argv))
        now = self._now()
        cached = self._cache.get(key)

        # Same question asked again within the interval: serve the answer
        # we already have, do not ask the scheduler again.
        if cached is not None and now - cached[0] < self.min_interval_s:
            self.served_from_cache += 1
            return cached[1]

        # A *different* argv of a kind whose gate is still closed must not
        # start a fresh scheduler call. Re-serve the last answer for this
        # exact argv if we have one; otherwise the gate is closed, not asked.
        last_kind = self._kind_last_run.get(kind)
        if last_kind is not None and now - last_kind < self.min_interval_s:
            if cached is not None:
                self.served_from_cache += 1
                return cached[1]
            return THROTTLED

        real_runner = runner or self._runner
        result = real_runner(list(argv))
        self._cache[key] = (now, result)
        self._kind_last_run[kind] = now
        self.queries_run += 1
        return result

    def last_result(self, kind: str, argv: Sequence[str]) -> object:
        """The cached answer for ``(kind, argv)`` or :data:`THROTTLED`."""
        return self._cache.get((kind, tuple(argv)), (0.0, None))[1]

    def clear(self) -> None:
        """Reset cache, gate and counters (tests, and a hedge against a
        long-running process whose kinds rotate)."""
        self._cache.clear()
        self._kind_last_run.clear()
        self.queries_run = 0
        self.served_from_cache = 0


def kind_of_argv(cmd: Sequence[str]) -> str:
    """The query kind of a scheduler command, from its basename.

    ``['/usr/bin/squeue', '-j', ...]`` -> ``"squeue"``. Unknown commands
    keep their basename, so a new scheduler tool is still throttled under
    its own kind rather than sharing ``squeue``'s gate.
    """
    base = os.path.basename(str(cmd[0]))
    return base or "unknown"


#: The process-wide client every module shares, so the throttle and the
#: cache are truly global and one component cannot reset another's by
#: constructing its own client.
default_client = SchedulerClient()


def query(
    kind: str,
    argv: Sequence[str],
    runner: Optional[Callable[[Sequence[str]], Optional[str]]] = None,
) -> object:
    """Module-level helper bound to :data:`default_client`."""
    return default_client.query(kind, argv, runner=runner)


def reset() -> None:
    """Clear the process-wide client (tests)."""
    default_client.clear()
