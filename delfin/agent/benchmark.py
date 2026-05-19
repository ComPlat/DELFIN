"""Canned-task benchmark suite for iterative agent optimisation.

Workflow:

1. Define tasks in ``delfin/agent/pack/benchmark/tasks.yaml`` (prompt +
   expected/forbidden signals + budget caps).
2. Run the suite against a model+mode combo — the runner (separate
   module) sends each prompt and collects the trajectory.
3. ``score_outcome(task, trajectory)`` produces a deterministic
   ``BenchmarkResult`` (success bool + 0-100 quality + components).
4. ``write_run`` appends each result to a JSONL file under
   ``~/.delfin/benchmark_runs/<ts>_<model>.jsonl``.
5. ``compare_runs(baseline_path, candidate_path)`` produces a per-task
   delta table — the "did this profile knob change help?" verdict.

The data plane is deliberately model-free: it knows nothing about the
engine or providers, only about (prompt, trajectory, outcome).  The
runner that drives the engine and supplies the trajectory lives in a
separate module so we can unit-test the scoring rubric without a live
provider.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

try:
    import yaml as _yaml
except ImportError:                                            # pragma: no cover
    _yaml = None  # type: ignore[assignment]


_DEFAULT_TASKS_PATH = (
    Path(__file__).resolve().parent / "pack" / "benchmark" / "tasks.yaml"
)
_DEFAULT_RUNS_DIR = Path.home() / ".delfin" / "benchmark_runs"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Signal:
    """One regex signal against the trajectory."""

    pattern: str
    against: str = "any"        # text | action | tool_name | any
    optional: bool = False


@dataclass(frozen=True)
class Task:
    """One canned benchmark task."""

    id: str
    task_class: str
    mode: str                   # solo | dashboard | plan | quick | …
    prompt: str
    expected_signals: tuple[Signal, ...] = ()
    forbidden_signals: tuple[Signal, ...] = ()
    max_duration_s: float = 60.0
    max_cost_usd: float = 0.10
    max_tool_calls: int = 5


@dataclass
class Trajectory:
    """What the runner observed during a single task execution.

    The scorer reads only this shape — the runner can build it from the
    engine event stream however it likes."""

    text: str = ""                  # concatenated assistant text
    actions: list[str] = field(default_factory=list)
    tool_calls: list[dict] = field(default_factory=list)
    duration_s: float = 0.0
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    error: str = ""

    def as_string(self) -> str:
        parts = [self.text]
        for a in self.actions:
            parts.append(f"\nACTION: {a}")
        for c in self.tool_calls:
            parts.append(f"\nTOOL: {c.get('name', '')}({c.get('input', '')})")
        return "".join(parts)


@dataclass
class BenchmarkResult:
    """Score-card for one (task, model, profile) execution."""

    task_id: str
    task_class: str
    model: str
    profile_name: str = ""
    mode: str = ""
    ts: float = 0.0
    success: bool = False
    quality_0_100: int = 0
    components: dict[str, int] = field(default_factory=dict)
    duration_s: float = 0.0
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    tool_calls: int = 0
    matched_signals: list[str] = field(default_factory=list)
    violated_signals: list[str] = field(default_factory=list)
    missing_signals: list[str] = field(default_factory=list)
    budget_violations: list[str] = field(default_factory=list)
    error: str = ""


# ---------------------------------------------------------------------------
# Task loading
# ---------------------------------------------------------------------------


def _coerce_signal(raw: Any) -> Signal:
    if isinstance(raw, Signal):
        return raw
    if isinstance(raw, str):
        return Signal(pattern=raw)
    if isinstance(raw, dict):
        return Signal(
            pattern=str(raw.get("pattern", "")),
            against=str(raw.get("against", "any")),
            optional=bool(raw.get("optional", False)),
        )
    raise TypeError(f"Cannot coerce signal: {raw!r}")


def _coerce_task(raw: dict) -> Task:
    expected = tuple(_coerce_signal(s) for s in raw.get("expected_signals") or [])
    forbidden = tuple(_coerce_signal(s) for s in raw.get("forbidden_signals") or [])
    return Task(
        id=str(raw["id"]),
        task_class=str(raw.get("task_class", "")),
        mode=str(raw.get("mode", "solo")),
        prompt=str(raw.get("prompt", "")),
        expected_signals=expected,
        forbidden_signals=forbidden,
        max_duration_s=float(raw.get("max_duration_s", 60.0)),
        max_cost_usd=float(raw.get("max_cost_usd", 0.10)),
        max_tool_calls=int(raw.get("max_tool_calls", 5)),
    )


def load_tasks(path: Path | None = None) -> list[Task]:
    """Load tasks from a YAML file.  Defaults to the packaged suite."""
    p = path or _DEFAULT_TASKS_PATH
    if not p.exists():
        return []
    if _yaml is None:                                           # pragma: no cover
        raise RuntimeError("PyYAML is required to load benchmark tasks")
    raw = _yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    items = raw.get("tasks") if isinstance(raw, dict) else raw
    if not isinstance(items, list):
        return []
    return [_coerce_task(t) for t in items if isinstance(t, dict)]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


_SIGNAL_AGAINST_VALUES = {"text", "action", "tool_name", "any"}


def _signal_matches(signal: Signal, traj: Trajectory) -> bool:
    pat = signal.pattern
    if not pat:
        return False
    try:
        rx = re.compile(pat)
    except re.error:
        return False
    against = signal.against if signal.against in _SIGNAL_AGAINST_VALUES else "any"
    haystacks: list[str]
    if against == "text":
        haystacks = [traj.text]
    elif against == "action":
        haystacks = list(traj.actions)
    elif against == "tool_name":
        haystacks = [c.get("name", "") for c in traj.tool_calls]
    else:                                                       # any
        haystacks = [traj.as_string()]
    return any(rx.search(h or "") for h in haystacks)


def score_outcome(
    task: Task,
    traj: Trajectory,
    *,
    model: str = "",
    profile_name: str = "",
    ts: float | None = None,
) -> BenchmarkResult:
    """Deterministic rubric: success + quality 0-100 + component breakdown."""

    matched: list[str] = []
    missing: list[str] = []
    violated: list[str] = []
    budget_violations: list[str] = []

    # 1. Expected signals — each non-optional MUST match for success.
    success_required_ok = True
    optional_hits = 0
    optional_total = 0
    for idx, sig in enumerate(task.expected_signals):
        label = f"{task.id}.expected[{idx}]"
        hit = _signal_matches(sig, traj)
        if hit:
            matched.append(label)
            if sig.optional:
                optional_hits += 1
        else:
            if sig.optional:
                missing.append(label + ":optional")
            else:
                missing.append(label)
                success_required_ok = False
        if sig.optional:
            optional_total += 1

    # 2. Forbidden signals — any match flips success to False.
    for idx, sig in enumerate(task.forbidden_signals):
        label = f"{task.id}.forbidden[{idx}]"
        if _signal_matches(sig, traj):
            violated.append(label)

    success = bool(success_required_ok and not violated and not traj.error)

    # 3. Budget checks (don't flip success on their own — visible in score).
    if traj.duration_s > task.max_duration_s > 0:
        budget_violations.append(
            f"duration_s {traj.duration_s:.1f}>{task.max_duration_s:.1f}"
        )
    if traj.cost_usd > task.max_cost_usd > 0:
        budget_violations.append(
            f"cost_usd {traj.cost_usd:.4f}>{task.max_cost_usd:.4f}"
        )
    if len(traj.tool_calls) > task.max_tool_calls > 0:
        budget_violations.append(
            f"tool_calls {len(traj.tool_calls)}>{task.max_tool_calls}"
        )

    # 4. Component scoring (0-100):
    #    success_pts  (40)   binary did-it-work
    #    routing_pts  (20)   expected signals + optional bonus
    #    speed_pts    (15)   within duration budget, scaled
    #    cost_pts     (10)   within cost budget, scaled
    #    clean_pts    (15)   no forbidden, no error, within tool-call budget
    success_pts = 40 if success else 0
    n_required = sum(1 for s in task.expected_signals if not s.optional)
    n_required_hit = n_required - sum(
        1 for m in missing if not m.endswith(":optional")
    )
    if n_required > 0:
        routing_base = int(round(15 * n_required_hit / n_required))
    else:
        routing_base = 15
    routing_bonus = (
        int(round(5 * optional_hits / max(1, optional_total)))
        if optional_total > 0
        else 0
    )
    routing_pts = routing_base + routing_bonus

    if task.max_duration_s > 0 and traj.duration_s >= 0:
        speed_ratio = min(1.0, traj.duration_s / task.max_duration_s)
        speed_pts = int(round(15 * (1.0 - speed_ratio)))
    else:
        speed_pts = 15

    if task.max_cost_usd > 0 and traj.cost_usd >= 0:
        cost_ratio = min(1.0, traj.cost_usd / task.max_cost_usd)
        cost_pts = int(round(10 * (1.0 - cost_ratio)))
    else:
        cost_pts = 10

    clean_pts = 15
    if violated:
        clean_pts -= 8
    if traj.error:
        clean_pts -= 5
    if len(traj.tool_calls) > task.max_tool_calls > 0:
        clean_pts -= 2
    clean_pts = max(0, clean_pts)

    components = {
        "success_pts": success_pts,
        "routing_pts": routing_pts,
        "speed_pts": speed_pts,
        "cost_pts": cost_pts,
        "clean_pts": clean_pts,
    }
    quality = max(0, min(100, sum(components.values())))

    return BenchmarkResult(
        task_id=task.id,
        task_class=task.task_class,
        model=model,
        profile_name=profile_name,
        mode=task.mode,
        ts=ts if ts is not None else time.time(),
        success=success,
        quality_0_100=quality,
        components=components,
        duration_s=traj.duration_s,
        cost_usd=traj.cost_usd,
        input_tokens=traj.input_tokens,
        output_tokens=traj.output_tokens,
        tool_calls=len(traj.tool_calls),
        matched_signals=matched,
        violated_signals=violated,
        missing_signals=missing,
        budget_violations=budget_violations,
        error=traj.error,
    )


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------


def _slug(s: str) -> str:
    out = "".join(c if c.isalnum() or c in "-._" else "_" for c in s)
    return out[:64] or "unknown"


def write_run(
    results: Iterable[BenchmarkResult],
    *,
    model: str,
    runs_dir: Path | None = None,
    run_id: str | None = None,
) -> Path:
    """Persist all results from one benchmark run as JSONL."""

    d = runs_dir or _DEFAULT_RUNS_DIR
    d.mkdir(parents=True, exist_ok=True)
    ts = int(time.time())
    rid = run_id or f"{ts}_{_slug(model)}"
    path = d / f"{rid}.jsonl"
    with path.open("w", encoding="utf-8") as f:
        for r in results:
            f.write(json.dumps(asdict(r), ensure_ascii=False) + "\n")
    return path


def read_run(path: Path) -> list[dict]:
    if not path.exists():
        return []
    out: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out


# ---------------------------------------------------------------------------
# Aggregation + comparison
# ---------------------------------------------------------------------------


def summarise_run(results: list[dict] | list[BenchmarkResult]) -> dict[str, Any]:
    """Roll-up across a single run: pass-rate + avg quality + totals."""

    rows = [
        asdict(r) if isinstance(r, BenchmarkResult) else r for r in results
    ]
    if not rows:
        return {
            "n_tasks": 0, "n_pass": 0, "pass_rate": 0.0,
            "avg_quality": 0.0, "total_cost_usd": 0.0,
            "total_duration_s": 0.0, "total_tool_calls": 0,
        }
    n = len(rows)
    n_pass = sum(1 for r in rows if r.get("success"))
    return {
        "n_tasks": n,
        "n_pass": n_pass,
        "pass_rate": n_pass / n,
        "avg_quality": sum(int(r.get("quality_0_100") or 0) for r in rows) / n,
        "total_cost_usd": sum(float(r.get("cost_usd") or 0) for r in rows),
        "total_duration_s": sum(float(r.get("duration_s") or 0) for r in rows),
        "total_tool_calls": sum(int(r.get("tool_calls") or 0) for r in rows),
    }


def compare_runs(
    baseline: list[dict] | Path,
    candidate: list[dict] | Path,
) -> dict[str, Any]:
    """Per-task delta + roll-up.  Used to answer "did the knob help?"

    Direction conventions:
      - quality_0_100 / success → higher is better
      - cost_usd / duration_s / tool_calls → lower is better

    Returns ``{"per_task": [...], "summary": {...}, "verdict": "..."}``.
    Verdict is one of {"better", "worse", "mixed", "neutral", "thin"} where
    ``thin`` means fewer than 3 overlapping tasks (too little signal).
    """

    if isinstance(baseline, Path):
        baseline = read_run(baseline)
    if isinstance(candidate, Path):
        candidate = read_run(candidate)

    by_id_old = {r.get("task_id"): r for r in baseline}
    by_id_new = {r.get("task_id"): r for r in candidate}
    overlap = sorted(set(by_id_old) & set(by_id_new))

    per_task: list[dict] = []
    n_better = n_worse = n_neutral = 0
    for tid in overlap:
        o = by_id_old[tid]
        n = by_id_new[tid]
        d_quality = int(n.get("quality_0_100") or 0) - int(o.get("quality_0_100") or 0)
        d_cost = float(n.get("cost_usd") or 0) - float(o.get("cost_usd") or 0)
        d_dur = float(n.get("duration_s") or 0) - float(o.get("duration_s") or 0)
        d_tools = int(n.get("tool_calls") or 0) - int(o.get("tool_calls") or 0)
        # Classification per task:
        # better if quality strictly up OR (quality flat AND cost or duration
        # meaningfully down ≥10% or absolute 1s/1¢).
        improved = False
        regressed = False
        if d_quality > 0:
            improved = True
        elif d_quality < 0:
            regressed = True
        else:
            base_cost = float(o.get("cost_usd") or 0)
            base_dur = float(o.get("duration_s") or 0)
            cost_drop = d_cost < -0.01 or (
                base_cost > 0 and d_cost / base_cost <= -0.10
            )
            dur_drop = d_dur < -1.0 or (
                base_dur > 0 and d_dur / base_dur <= -0.10
            )
            cost_rise = d_cost > 0.01 or (
                base_cost > 0 and d_cost / base_cost >= 0.10
            )
            dur_rise = d_dur > 1.0 or (
                base_dur > 0 and d_dur / base_dur >= 0.10
            )
            if cost_drop or dur_drop:
                improved = True
            elif cost_rise or dur_rise:
                regressed = True
        if improved and not regressed:
            cls = "better"
            n_better += 1
        elif regressed and not improved:
            cls = "worse"
            n_worse += 1
        else:
            cls = "neutral"
            n_neutral += 1
        per_task.append({
            "task_id": tid,
            "class": cls,
            "old_quality": int(o.get("quality_0_100") or 0),
            "new_quality": int(n.get("quality_0_100") or 0),
            "d_quality": d_quality,
            "d_cost_usd": d_cost,
            "d_duration_s": d_dur,
            "d_tool_calls": d_tools,
            "old_success": bool(o.get("success")),
            "new_success": bool(n.get("success")),
        })

    summary_old = summarise_run(baseline)
    summary_new = summarise_run(candidate)
    summary = {
        "n_overlap": len(overlap),
        "old": summary_old,
        "new": summary_new,
        "n_better": n_better,
        "n_worse": n_worse,
        "n_neutral": n_neutral,
    }

    if len(overlap) < 3:
        verdict = "thin"
    elif n_better >= max(3, int(0.8 * len(overlap))) and n_worse == 0:
        verdict = "better"
    elif n_worse >= max(3, int(0.5 * len(overlap))):
        verdict = "worse"
    elif n_better > n_worse + 1:
        verdict = "better"
    elif n_worse > n_better + 1:
        verdict = "worse"
    elif n_better > 0 and n_worse > 0:
        verdict = "mixed"
    else:
        verdict = "neutral"

    return {"per_task": per_task, "summary": summary, "verdict": verdict}


__all__ = [
    "Signal",
    "Task",
    "Trajectory",
    "BenchmarkResult",
    "load_tasks",
    "score_outcome",
    "write_run",
    "read_run",
    "summarise_run",
    "compare_runs",
]
