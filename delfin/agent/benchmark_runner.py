"""Live-engine runner for the canned-task benchmark suite.

The data plane (``benchmark.py``) defines Task / Trajectory / scoring; this
module bridges to the real ``AgentEngine`` and feeds it the prompts.  A
single function — ``run_suite`` — is the entire surface area:

    results = run_suite(load_tasks(), model="kit.qwen3.5-397b-A17b",
                        backend="api", provider="kit")
    path = write_run(results, model="kit.qwen3.5-397b-A17b")

An engine factory can be injected for unit-testing — the real factory
defers all heavy imports until invocation so a bench CLI parse fails
fast on a bad argument before pulling the engine stack.
"""

from __future__ import annotations

import re
import time
from typing import Any, Callable, Iterable

from .benchmark import (
    BenchmarkResult,
    Task,
    Trajectory,
    score_outcome,
)


# ---------------------------------------------------------------------------
# Trajectory builder — convert _run_once output into a Trajectory.
# ---------------------------------------------------------------------------


# Same five ACTION variants tab_agent.py accepts in the dashboard.
_ACTION_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^\s*ACTION\s*:\s*(/\S.*?)\s*$", re.MULTILINE),
    re.compile(r"^\s*ACTION\s+(/\S.*?)\s*$", re.MULTILINE),
    re.compile(r"^\s*Action\s*:\s*(/\S.*?)\s*$", re.MULTILINE),
    re.compile(r"`(/(?:tab|control|orca|effort|mode|provider|model)\s+\S[^`]*)`"),
)
_BARE_SLASH_PREFIXES = (
    "/tab ", "/control ", "/orca ", "/effort ", "/mode ", "/provider ",
    "/model ",
)


def extract_actions(text: str) -> list[str]:
    """Pull ACTION lines (canonical + tolerant variants) from text.

    Mirrors the dashboard's _action_cmd parser so the benchmark sees the
    same routing decisions an end user would.  Bare-slash lines (no
    ``ACTION:`` prefix) are accepted only when the slash command starts
    with a known dashboard prefix — that guards against false positives
    like ``/home/user/file.py``.
    """

    if not text:
        return []
    out: list[str] = []
    for rx in _ACTION_PATTERNS:
        for m in rx.finditer(text):
            cmd = (m.group(1) or "").strip()
            if cmd and cmd not in out:
                out.append(cmd)
    # Bare-slash lines (whitelisted prefixes only)
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("/"):
            continue
        if any(s.startswith(p) for p in _BARE_SLASH_PREFIXES):
            if s not in out:
                out.append(s)
    return out


def trajectory_from_run(raw: dict, *, duration_s: float, cost_usd: float = 0.0) -> Trajectory:
    """Convert ``_run_once``'s return dict into a Trajectory."""
    text = str(raw.get("text") or "")
    return Trajectory(
        text=text,
        actions=extract_actions(text),
        tool_calls=list(raw.get("tool_calls") or []),
        duration_s=float(duration_s),
        cost_usd=float(cost_usd),
        input_tokens=int(raw.get("input_tokens") or 0),
        output_tokens=int(raw.get("output_tokens") or 0),
        error=str(raw.get("error") or ""),
    )


# ---------------------------------------------------------------------------
# Engine factory + cost estimator
# ---------------------------------------------------------------------------


EngineFactory = Callable[[str, str, str, str], Any]
"""(model, backend, provider, mode) -> engine."""


def _default_engine_factory(model: str, backend: str, provider: str, mode: str) -> Any:
    """Build a real AgentEngine for the given config.

    AgentEngine creates its own client internally; we just hand it the
    resolved tuple and let it own the lifecycle.  Heavy imports happen
    here, not at module-import time.
    """
    import os as _os
    from pathlib import Path as _Path
    from .engine import AgentEngine

    return AgentEngine(
        repo_dir=_Path(_os.getcwd()).resolve(),
        backend=backend or "api",
        provider=provider,
        model=model,
        mode=mode or "solo",
    )


def _cost_delta(before: float, after: float) -> float:
    """Δ cost for a single turn — defends against engines that don't
    expose cost_usd."""
    try:
        return max(0.0, float(after) - float(before))
    except (TypeError, ValueError):
        return 0.0


# ---------------------------------------------------------------------------
# Public runner
# ---------------------------------------------------------------------------


def _run_task_once(
    task: Task,
    *,
    model: str,
    backend: str,
    provider: str,
    profile_name: str,
    engine_factory: EngineFactory,
    max_tokens: int,
    run_once: Callable[..., dict],
    clock: Callable[[], float],
) -> BenchmarkResult:
    """Single attempt — kept private so retry-aggregation logic lives
    in one place at the public ``run_task`` entry."""
    try:
        engine = engine_factory(model, backend, provider, task.mode)
    except Exception as exc:
        traj = Trajectory(error=f"engine init failed: {exc}")
        return score_outcome(
            task, traj, model=model, profile_name=profile_name, ts=clock(),
        )
    cost_before = float(getattr(engine, "cost_usd", 0.0) or 0.0)
    t0 = clock()
    try:
        raw = run_once(engine, task.prompt, max_tokens=max_tokens)
    except Exception as exc:
        raw = {"text": "", "tool_calls": [], "input_tokens": 0,
               "output_tokens": 0, "error": f"_run_once raised: {exc}"}
    t1 = clock()
    cost_after = float(getattr(engine, "cost_usd", 0.0) or 0.0)
    traj = trajectory_from_run(
        raw, duration_s=t1 - t0, cost_usd=_cost_delta(cost_before, cost_after),
    )
    return score_outcome(
        task, traj, model=model, profile_name=profile_name, ts=clock(),
    )


def run_task(
    task: Task,
    *,
    model: str,
    backend: str = "api",
    provider: str = "",
    profile_name: str = "",
    engine_factory: EngineFactory | None = None,
    max_tokens: int = 1024,
    repeats: int = 1,
    run_once: Callable[..., dict] | None = None,
    clock: Callable[[], float] | None = None,
    on_replicate: Callable[[int, BenchmarkResult], None] | None = None,
) -> BenchmarkResult:
    """Execute one task and return a scored result.

    ``repeats=1`` (default): single sample, returned as-is.
    ``repeats>1``: run the task N times against a FRESH engine each time
    (so prior-turn state doesn't bias the result) and return a single
    aggregated BenchmarkResult with median-quality + majority-success
    + quality_stdev as a noise indicator.  ``on_replicate`` (if given)
    is called after each individual run with ``(idx, result)`` so a
    caller can stream progress.

    The defaults wire up the real ``AgentEngine`` + ``_run_once`` from
    ``cli.py``; tests can inject stubs.
    """
    now = clock or time.time
    factory = engine_factory or _default_engine_factory
    if run_once is None:
        from .cli import _run_once as _real_run_once
        run_once = _real_run_once

    n = max(1, int(repeats))
    replicates: list[BenchmarkResult] = []
    for i in range(n):
        r = _run_task_once(
            task,
            model=model, backend=backend, provider=provider,
            profile_name=profile_name, engine_factory=factory,
            max_tokens=max_tokens, run_once=run_once, clock=now,
        )
        replicates.append(r)
        if on_replicate is not None:
            try:
                on_replicate(i, r)
            except Exception:
                pass

    if n == 1:
        return replicates[0]
    from .benchmark import aggregate_replicates
    return aggregate_replicates(replicates)


def run_suite(
    tasks: Iterable[Task],
    *,
    model: str,
    backend: str = "api",
    provider: str = "",
    profile_name: str = "",
    engine_factory: EngineFactory | None = None,
    max_tokens: int = 1024,
    repeats: int = 1,
    progress: Callable[[Task, BenchmarkResult], None] | None = None,
    on_replicate: Callable[[Task, int, BenchmarkResult], None] | None = None,
    run_once: Callable[..., dict] | None = None,
) -> list[BenchmarkResult]:
    """Run every task and return scored (optionally aggregated) results.

    ``repeats`` is forwarded to ``run_task`` so each task is sampled
    N times; the returned list still has one entry per task with the
    aggregate stats.  ``progress`` fires once per task with the final
    aggregated result; ``on_replicate`` fires for each individual
    sample as it lands.
    """

    out: list[BenchmarkResult] = []
    for task in tasks:
        per_task_replicate: Callable[..., None] | None = (
            (lambda t: (lambda i, r: on_replicate(t, i, r)))(task)
            if on_replicate is not None else None
        )
        result = run_task(
            task,
            model=model, backend=backend, provider=provider,
            profile_name=profile_name, engine_factory=engine_factory,
            max_tokens=max_tokens, repeats=repeats,
            run_once=run_once,
            on_replicate=per_task_replicate,
        )
        out.append(result)
        if progress is not None:
            try:
                progress(task, result)
            except Exception:
                pass  # best-effort progress reporting
    return out


# ---------------------------------------------------------------------------
# Profile-name resolution
# ---------------------------------------------------------------------------


def resolve_profile_name(model: str) -> str:
    """Return the profile's ``notes`` field (or 'default') for stamping
    onto results — useful for filtering historical runs by profile."""
    try:
        from .model_profiles import get_profile
        p = get_profile(model)
        return (p.notes or "default")[:80]
    except Exception:
        return "default"


__all__ = [
    "extract_actions",
    "trajectory_from_run",
    "run_task",
    "run_suite",
    "resolve_profile_name",
    "EngineFactory",
]
