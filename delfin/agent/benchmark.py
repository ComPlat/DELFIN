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
import os
import re
import subprocess
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
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
    """Score-card for one (task, model, profile) execution.

    With ``n_samples > 1`` this is an aggregated record over N
    repeated executions of the same task; ``quality_stdev`` and
    ``success_rate`` then expose run-to-run noise.  ``components``,
    ``duration_s`` etc are medians across the N replicates and
    ``cost_usd`` is the SUM of the N runs (the real spend), not a
    median — useful for budget tracking.
    """

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
    # --- replicate-aware fields (N=1 → trivially the only sample) ---
    n_samples: int = 1
    quality_stdev: float = 0.0
    success_rate: float = 0.0           # fraction of N samples with success=True
    per_run_quality: list[int] = field(default_factory=list)
    per_run_success: list[bool] = field(default_factory=list)
    # --- forensic fields for pattern-bug-vs-real-fail diagnosis ---
    text_excerpt: str = ""              # first ≤400 chars of model output
    tool_names: list[str] = field(default_factory=list)  # tools the model actually called


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

    excerpt = (traj.text or "")[:400].strip()
    tool_names = [str(c.get("name", "")) for c in traj.tool_calls
                  if c.get("name")]

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
        text_excerpt=excerpt,
        tool_names=tool_names,
    )


# ---------------------------------------------------------------------------
# Replicate aggregation
# ---------------------------------------------------------------------------


def _median(values: list[float]) -> float:
    """Numerically-stable median for small lists.  Empty → 0.0."""
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 0:
        return (s[mid - 1] + s[mid]) / 2.0
    return float(s[mid])


def _stdev(values: list[float]) -> float:
    """Sample standard deviation (Bessel's correction).  Empty/single → 0."""
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return var ** 0.5


def aggregate_replicates(
    results: list[BenchmarkResult],
) -> BenchmarkResult:
    """Collapse N per-run results for the SAME task into one aggregated
    BenchmarkResult.

    Rules:
      - ``success`` = majority vote (ceil(N/2) successes required to pass)
      - ``quality_0_100`` = median across N runs
      - ``duration_s`` / ``tool_calls`` / ``input_tokens`` / ``output_tokens``
        = median (typical observation, not skewed by one outlier)
      - ``cost_usd`` = SUM (real total spend across replicates)
      - ``components`` = median per component
      - ``matched_signals`` / ``missing_signals`` / ``violated_signals`` /
        ``budget_violations`` = union across replicates (so a flaky
        violation is still visible)
      - ``quality_stdev`` exposes run-to-run noise; ``success_rate``
        the fraction of passes; ``per_run_*`` keep the raw samples for
        deeper analysis.

    Raises ``ValueError`` on empty input or task_id mismatch.
    """
    if not results:
        raise ValueError("aggregate_replicates needs at least one result")
    first = results[0]
    if not all(r.task_id == first.task_id for r in results):
        raise ValueError("aggregate_replicates: all results must share task_id")

    n = len(results)
    qualities = [int(r.quality_0_100) for r in results]
    durations = [float(r.duration_s) for r in results]
    tool_counts = [int(r.tool_calls) for r in results]
    in_toks = [int(r.input_tokens) for r in results]
    out_toks = [int(r.output_tokens) for r in results]
    success_flags = [bool(r.success) for r in results]
    n_pass = sum(1 for f in success_flags if f)

    # Components: median per field across replicates
    comp_keys: set[str] = set()
    for r in results:
        comp_keys.update((r.components or {}).keys())
    components = {
        k: int(round(_median([float((r.components or {}).get(k, 0)) for r in results])))
        for k in comp_keys
    }

    # Union of signal labels (preserves order of first occurrence)
    def _union(field_name: str) -> list[str]:
        out: list[str] = []
        for r in results:
            for s in getattr(r, field_name, []) or []:
                if s not in out:
                    out.append(s)
        return out

    # Pick the first non-empty excerpt — gives forensic value without
    # bloating storage with N copies; tool_names unioned the same way
    # as signals (a flaky tool-call still surfaces).
    excerpt = ""
    for r in results:
        if r.text_excerpt:
            excerpt = r.text_excerpt
            break
    tool_names_union: list[str] = []
    for r in results:
        for n_name in (r.tool_names or []):
            if n_name and n_name not in tool_names_union:
                tool_names_union.append(n_name)

    return BenchmarkResult(
        task_id=first.task_id,
        task_class=first.task_class,
        model=first.model,
        profile_name=first.profile_name,
        mode=first.mode,
        ts=first.ts,
        success=(n_pass * 2 >= n),                   # majority (ties = pass)
        quality_0_100=int(round(_median([float(q) for q in qualities]))),
        components=components,
        duration_s=_median(durations),
        cost_usd=sum(float(r.cost_usd) for r in results),
        input_tokens=int(_median([float(x) for x in in_toks])),
        output_tokens=int(_median([float(x) for x in out_toks])),
        tool_calls=int(_median([float(x) for x in tool_counts])),
        matched_signals=_union("matched_signals"),
        violated_signals=_union("violated_signals"),
        missing_signals=_union("missing_signals"),
        budget_violations=_union("budget_violations"),
        error="; ".join(sorted({r.error for r in results if r.error}))[:500],
        n_samples=n,
        quality_stdev=_stdev([float(q) for q in qualities]),
        success_rate=n_pass / n,
        per_run_quality=list(qualities),
        per_run_success=list(success_flags),
        text_excerpt=excerpt,
        tool_names=tool_names_union,
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


# ---------------------------------------------------------------------------
# Markdown export + profile-commit linking
# ---------------------------------------------------------------------------


def runs_dir() -> Path:
    """Public accessor for the default runs directory."""
    return _DEFAULT_RUNS_DIR


def run_timestamp(path: Path) -> float:
    """Return the ts of the first record (≈ run start) or 0 if unreadable.

    Used to bracket ``git log`` so we can auto-find which profile commit
    sits between a baseline and a candidate run.
    """
    records = read_run(path)
    if records:
        try:
            return float(records[0].get("ts") or 0)
        except (TypeError, ValueError):
            return 0.0
    return 0.0


def find_profile_commits_between(
    old_ts: float,
    new_ts: float,
    *,
    profile_file: str = "delfin/agent/model_profiles.py",
    repo_root: str | Path | None = None,
    timeout_s: float = 5.0,
) -> list[str]:
    """Return ``<short-hash> <subject>`` lines for commits that touched
    the profile file in the (old_ts, new_ts] window.

    Used to auto-annotate compare reports so we know which knob change
    drove the observed Δ.  Returns ``[]`` on any git failure (missing
    repo, no commits, timeout) — never raises.
    """

    if old_ts <= 0 or new_ts <= 0 or new_ts <= old_ts:
        return []
    cwd = str(repo_root) if repo_root else os.getcwd()
    # Pad the window by 1 s on each side: git stores commit timestamps at
    # second precision, so a commit at exactly `new_ts` would otherwise
    # be excluded by `--until` (which is strict "before").
    try:
        old_iso = datetime.fromtimestamp(max(0.0, old_ts - 1.0)).isoformat()
        new_iso = datetime.fromtimestamp(new_ts + 1.0).isoformat()
    except (OSError, ValueError, OverflowError):
        return []
    try:
        out = subprocess.run(
            [
                "git", "log", "--oneline",
                f"--since={old_iso}", f"--until={new_iso}",
                "--", profile_file,
            ],
            capture_output=True, text=True, timeout=timeout_s, cwd=cwd,
        )
    except (subprocess.SubprocessError, OSError):
        return []
    if out.returncode != 0:
        return []
    return [line.strip() for line in out.stdout.splitlines() if line.strip()]


_VERDICT_EMOJI = {
    "better": "[+]", "worse": "[-]", "mixed": "[~]",
    "neutral": "[=]", "thin": "[?]",
}
_CLASS_MARK = {"better": "[+]", "worse": "[-]", "neutral": "[=]"}


def format_compare_markdown(
    cmp_result: dict,
    *,
    baseline_path: Path | None = None,
    candidate_path: Path | None = None,
    include_profile_commits: bool = True,
    repo_root: str | Path | None = None,
) -> str:
    """Render a ``compare_runs`` result as a markdown report.

    Suitable for pasting into a PR body, memory entry, or chat message.
    If both baseline_path and candidate_path are supplied AND
    ``include_profile_commits=True``, the report also lists any commits
    that touched ``delfin/agent/model_profiles.py`` between the two run
    timestamps — the "which knob change drove this delta?" annotation.
    """

    summary = cmp_result.get("summary") or {}
    verdict = str(cmp_result.get("verdict") or "neutral")
    old = summary.get("old") or {}
    new = summary.get("new") or {}

    lines: list[str] = ["# Benchmark Comparison Report", ""]
    lines.append(
        f"**Verdict: {verdict.upper()}** "
        f"{_VERDICT_EMOJI.get(verdict, '[?]')}"
    )
    lines.append(
        f"  {summary.get('n_better', 0)} better, "
        f"{summary.get('n_worse', 0)} worse, "
        f"{summary.get('n_neutral', 0)} neutral  "
        f"(n={summary.get('n_overlap', 0)})"
    )
    if baseline_path is not None:
        lines.append(f"  Baseline:  `{baseline_path}`")
    if candidate_path is not None:
        lines.append(f"  Candidate: `{candidate_path}`")
    lines.append("")

    # Summary table
    lines.append("## Summary")
    lines.append("")
    lines.append("|              | Baseline   | Candidate  | Δ          |")
    lines.append("|--------------|------------|------------|------------|")
    op = float(old.get("pass_rate") or 0)
    np_ = float(new.get("pass_rate") or 0)
    lines.append(
        f"| Pass rate    | {op:>9.0%}  | {np_:>9.0%}  | "
        f"{(np_ - op) * 100:+8.0f}pp |"
    )
    oq = float(old.get("avg_quality") or 0)
    nq = float(new.get("avg_quality") or 0)
    lines.append(
        f"| Avg quality  | {oq:>10.1f} | {nq:>10.1f} | {nq - oq:+10.1f} |"
    )
    oc = float(old.get("total_cost_usd") or 0)
    nc = float(new.get("total_cost_usd") or 0)
    lines.append(
        f"| Total cost   | ${oc:>9.4f} | ${nc:>9.4f} | "
        f"${nc - oc:+9.4f} |"
    )
    od = float(old.get("total_duration_s") or 0)
    nd = float(new.get("total_duration_s") or 0)
    lines.append(
        f"| Total time   | {od:>9.1f}s | {nd:>9.1f}s | {nd - od:+9.1f}s |"
    )
    lines.append("")

    # Per-task table
    per_task = cmp_result.get("per_task") or []
    if per_task:
        lines.append("## Per-task")
        lines.append("")
        lines.append(
            "| Task                          | Class   | Quality | Δq   | "
            "Δcost     | Δdur      |"
        )
        lines.append(
            "|-------------------------------|---------|---------|------|"
            "-----------|-----------|"
        )
        for row in per_task:
            mark = _CLASS_MARK.get(row.get("class") or "", "[?]")
            tid = str(row.get("task_id") or "")[:30]
            oq2 = int(row.get("old_quality") or 0)
            nq2 = int(row.get("new_quality") or 0)
            dq = int(row.get("d_quality") or 0)
            dc = float(row.get("d_cost_usd") or 0)
            dd = float(row.get("d_duration_s") or 0)
            lines.append(
                f"| `{tid:<29}` | {mark:<7} | {oq2:>3}→{nq2:<3} | "
                f"{dq:+5d} | {dc:+9.4f} | {dd:+8.2f}s |"
            )
        lines.append("")

    # Profile-commit annotation
    if (include_profile_commits and baseline_path is not None
            and candidate_path is not None):
        try:
            old_ts = run_timestamp(baseline_path)
            new_ts = run_timestamp(candidate_path)
            commits = find_profile_commits_between(
                old_ts, new_ts, repo_root=repo_root,
            )
        except Exception:
            commits = []
        if commits:
            lines.append("## Profile changes between runs")
            lines.append("")
            lines.append(
                "Commits that modified "
                "`delfin/agent/model_profiles.py`:"
            )
            lines.append("")
            for c in commits:
                lines.append(f"- `{c}`")
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# Audit — pattern-bug-vs-real-fail diagnosis
# ---------------------------------------------------------------------------


def audit_run(
    records: list[dict] | Path,
    *,
    tasks: list[Task] | None = None,
) -> list[dict]:
    """Return one audit entry per FAILED task in a run.

    A failure is anything with ``success == False``.  For each, the
    entry surfaces enough state to decide pattern-bug-vs-real-fail in
    one glance:

      - task_id / model
      - quality / success_rate / quality_stdev (low σ + low rate
        = deterministic; high σ = flaky)
      - missing / violated signal labels (which check failed)
      - text_excerpt (what the model actually wrote)
      - tool_names (which tools the model actually called)
      - hint_pattern_bug: True if the model output looks reasonable
        (excerpt is non-empty AND no error AND quality_stdev < 5),
        which means the signals are likely the problem, not the model.
    """
    if isinstance(records, Path):
        records = read_run(records)
    task_by_id = {t.id: t for t in (tasks or [])}
    if not task_by_id:
        try:
            task_by_id = {t.id: t for t in load_tasks()}
        except Exception:
            task_by_id = {}

    out: list[dict] = []
    for r in records:
        if r.get("success"):
            continue
        rate = float(r.get("success_rate") or 0)
        stdev = float(r.get("quality_stdev") or 0)
        excerpt = str(r.get("text_excerpt") or "")
        error = str(r.get("error") or "")
        violated = list(r.get("violated_signals") or [])
        # Pattern-bug heuristic: deterministic fail + non-empty reasonable-
        # looking output + no engine error + NO forbidden violation.  A
        # violated forbidden signal means the model did something explicitly
        # disallowed → that's real misbehaviour, never a pattern bug.
        hint_pattern_bug = (
            rate <= 0.34
            and stdev < 5.0
            and len(excerpt) >= 30
            and not error
            and not violated
        )
        task = task_by_id.get(r.get("task_id") or "")
        signal_defs: dict[str, str] = {}
        if task is not None:
            for idx, sig in enumerate(task.expected_signals):
                signal_defs[f"{task.id}.expected[{idx}]"] = (
                    f"{sig.pattern}   (against={sig.against})"
                )
            for idx, sig in enumerate(task.forbidden_signals):
                signal_defs[f"{task.id}.forbidden[{idx}]"] = (
                    f"{sig.pattern}   (against={sig.against})"
                )
        out.append({
            "task_id": r.get("task_id") or "",
            "model": r.get("model") or "",
            "quality": int(r.get("quality_0_100") or 0),
            "success_rate": rate,
            "quality_stdev": stdev,
            "missing_signals": list(r.get("missing_signals") or []),
            "violated_signals": violated,
            "tool_names": list(r.get("tool_names") or []),
            "text_excerpt": excerpt,
            "error": error,
            "signal_defs": signal_defs,
            "hint_pattern_bug": hint_pattern_bug,
        })
    return out


def format_audit_report(entries: list[dict]) -> str:
    """Render audit entries as a developer-friendly text report."""
    if not entries:
        return "No failed tasks — nothing to audit.\n"
    lines: list[str] = []
    bug_entries = [e for e in entries if e["hint_pattern_bug"]]
    real_entries = [e for e in entries if not e["hint_pattern_bug"]]
    lines.append(f"=== AUDIT: {len(entries)} failed task(s) "
                 f"({len(bug_entries)} likely pattern-bug, "
                 f"{len(real_entries)} likely real-fail) ===\n")
    for group_name, group in (
        ("SUSPECTED PATTERN BUG", bug_entries),
        ("LIKELY REAL FAIL OR FLAKY", real_entries),
    ):
        if not group:
            continue
        lines.append(f"--- {group_name} ---")
        for e in group:
            lines.append(f"\n  task   : {e['task_id']}    (model={e['model']})")
            lines.append(f"  quality: q={e['quality']}  "
                         f"rate={e['success_rate']:.0%}  "
                         f"σ={e['quality_stdev']:.1f}")
            if e["tool_names"]:
                lines.append(f"  tools  : {', '.join(e['tool_names'])}")
            if e["missing_signals"]:
                lines.append("  missing signals:")
                for label in e["missing_signals"]:
                    if label.endswith(":optional"):
                        continue
                    defn = e["signal_defs"].get(label, "(definition not found)")
                    lines.append(f"    {label}")
                    lines.append(f"      pattern: {defn}")
            if e["violated_signals"]:
                lines.append("  VIOLATED signals:")
                for label in e["violated_signals"]:
                    defn = e["signal_defs"].get(label, "(definition not found)")
                    lines.append(f"    {label}")
                    lines.append(f"      pattern: {defn}")
            if e["error"]:
                lines.append(f"  error  : {e['error'][:200]}")
            if e["text_excerpt"]:
                lines.append("  model output (≤400 chars):")
                for line in e["text_excerpt"].splitlines()[:8]:
                    lines.append(f"    │ {line[:120]}")
                if len(e["text_excerpt"].splitlines()) > 8:
                    lines.append("    │ …")
        lines.append("")
    return "\n".join(lines) + "\n"


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
    "aggregate_replicates",
    "runs_dir",
    "run_timestamp",
    "find_profile_commits_between",
    "format_compare_markdown",
    "audit_run",
    "format_audit_report",
]
