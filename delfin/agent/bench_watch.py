"""Unattended benchmark operations: rolling baseline, drift watch, recheck.

The benchmark suite (``benchmark.py`` scoring + ``benchmark_runner.py``
execution) produces one JSONL file per run in ``~/.delfin/benchmark_runs``.
Single samples are NOISY — measured across two identical same-day 48-task
runs of the same model: 6 tasks flipped pass/fail and the per-task quality
delta had a stdev of ~20 points.  A naive "compare to last run" alarm
would therefore fire on almost every run.  This module implements a
noise-aware three-stage pipeline instead:

1. ``load_history``  — per-task rolling baseline over the last k samples
   of a model (majority success + median quality + quality spread).
2. ``compare_run``   — classify each task of a run against the baseline
   as ``stable`` / ``improved`` / ``suspect_regression`` / ``new_task``.
   The suspect threshold ``max(8, 1.5 * baseline_spread)`` is derived
   from the measured noise (see ``suspect_threshold``).
3. ``recheck_suspects`` — re-run ONLY the suspects with repeats and let
   majority-success + median-quality separate ``confirmed_regression``
   from ``noise``, under hard task-count and cost caps.

``nightly`` chains the stages into one unattended cycle (run suite →
compare → recheck → markdown report → one attention event on confirmed
regressions) and never raises — partial failures are reported inside the
returned summary.  All live execution (suite run, suspect rerun,
attention emission) is injectable so tests run fully offline.

Scheduling is OPT-IN: nothing in this module registers a scheduler
entry.  ``scheduler add-bench`` (see ``cli.py``) creates one explicitly,
and ``scheduler_daemon`` executes it via its bench-entry hook.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

from . import benchmark as _bm


# ---------------------------------------------------------------------------
# Thresholds (grounded in measured run-to-run noise)
# ---------------------------------------------------------------------------

# Floor for the quality-drop threshold.  Between two identical runs the
# per-task deltas are bimodal: |delta| <= ~7 points of scoring jitter
# (speed/cost point components move with wall-clock and token noise) or
# |delta| >= ~20 points from a structural pass/fail flip.  8 sits in the
# gap: above the jitter band, well below any structural change.
QUALITY_DROP_FLOOR = 8.0

# Multiplier on the per-task baseline spread.  Flagging at 1.5 sigma
# (one-sided) accepts roughly 7% false suspects per noisy task per run —
# deliberate: stage 2 (recheck with repeats + majority/median) is sized
# to absorb a handful of false suspects cheaply, and a confirmed alarm
# then requires BOTH a flagged drop AND a failed multi-sample recheck.
SPREAD_FACTOR = 1.5

# Recheck economics: the 48-task KIT-Qwen suite costs ~$8 per run at
# repeats=1, i.e. ~$0.17 per task-sample.  Used only when a task has no
# cost history of its own.
FALLBACK_TASK_COST_USD = 0.17

DEFAULT_LAST_K = 5
DEFAULT_RECHECK_REPEATS = 3
DEFAULT_MAX_RECHECK_TASKS = 6
DEFAULT_MAX_RECHECK_COST_USD = 3.0


def suspect_threshold(spread: float) -> float:
    """Quality-drop threshold for one task: ``max(8, 1.5 * spread)``.

    ``spread`` is the task's baseline quality stdev.  See the module
    constants for why 8 and 1.5: the floor filters ubiquitous scoring
    jitter (<= ~7 points between identical runs), the multiplier scales
    the band for tasks whose own history is noisy, and the ~7% one-sided
    false-positive rate this implies is intentionally handed to the
    capped recheck stage instead of being suppressed here.
    """
    try:
        s = max(0.0, float(spread))
    except (TypeError, ValueError):
        s = 0.0
    return max(QUALITY_DROP_FLOOR, SPREAD_FACTOR * s)


# ---------------------------------------------------------------------------
# History: rolling per-task baseline
# ---------------------------------------------------------------------------


def _runs_dir(runs_dir: Path | str | None) -> Path:
    return Path(runs_dir).expanduser() if runs_dir else _bm.runs_dir()


def list_model_runs(
    model: str, *, runs_dir: Path | str | None = None,
) -> list[Path]:
    """Run files for ``model``, newest first.

    Only canonical ``<epoch>_<model-slug>.jsonl`` names are matched, the
    format ``benchmark.write_run`` produces — special-purpose files with
    a custom ``run_id`` prefix (e.g. behaviour baselines) are excluded
    so they cannot skew the rolling baseline.
    """
    d = _runs_dir(runs_dir)
    slug = _bm._slug(str(model or ""))
    out: list[tuple[int, Path]] = []
    try:
        candidates = list(d.glob(f"*_{slug}.jsonl"))
    except OSError:
        return []
    for p in candidates:
        stem = p.name[: -len(".jsonl")]
        prefix = stem[: -(len(slug) + 1)] if stem.endswith(f"_{slug}") else ""
        if prefix.isdigit():
            out.append((int(prefix), p))
    out.sort(key=lambda t: (t[0], t[1].name), reverse=True)
    return [p for _, p in out]


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    vs = sorted(values)
    n = len(vs)
    mid = n // 2
    return float(vs[mid]) if n % 2 else (vs[mid - 1] + vs[mid]) / 2.0


def _stdev(values: list[float]) -> float:
    """Sample stdev (n-1); 0.0 for fewer than two samples."""
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    return (sum((v - mean) ** 2 for v in values) / (n - 1)) ** 0.5


def _expand_samples(record: dict) -> list[tuple[float, bool, float]]:
    """(quality, success, cost_per_sample) samples inside one record.

    A record with ``n_samples > 1`` is already an aggregate; its
    ``per_run_quality`` / ``per_run_success`` lists are expanded so the
    baseline sees the true sample distribution.  ``cost_usd`` of an
    aggregate is the SUM over samples, so it is split evenly here.
    """
    try:
        n = int(record.get("n_samples") or 1)
    except (TypeError, ValueError):
        n = 1
    prq = record.get("per_run_quality") or []
    prs = record.get("per_run_success") or []
    cost = float(record.get("cost_usd") or 0.0)
    if n > 1 and len(prq) == n and len(prs) == n:
        per_cost = cost / n if n else 0.0
        return [
            (float(q), bool(s), per_cost) for q, s in zip(prq, prs)
        ]
    return [(
        float(record.get("quality_0_100") or 0),
        bool(record.get("success")),
        cost,
    )]


def load_history(
    model: str,
    *,
    runs_dir: Path | str | None = None,
    last_k: int = DEFAULT_LAST_K,
    exclude: Path | str | None = None,
) -> dict:
    """Rolling per-task baseline over the model's recent run files.

    Files are scanned newest-first; each task's baseline uses samples
    from the newest ``last_k`` files that actually CONTAIN the task.
    Real run dirs are dominated by partial, task-filtered runs — a
    window counted in files regardless of content would starve most
    baselines down to a single sample, so the window is per task.
    Torn lines and missing tasks are tolerated (``benchmark.read_run``
    skips unparseable lines).  ``exclude`` drops one file from the
    window — used when the run under comparison already sits in the
    runs dir and must not be its own baseline.

    Returns ``{"model", "runs_dir", "last_k", "files", "tasks":
    {task_id: {"n_files", "n_samples", "majority_success",
    "pass_fraction", "median_quality", "quality_spread",
    "median_cost_usd"}}}``.
    """
    k = max(1, int(last_k or 1))
    files = list_model_runs(model, runs_dir=runs_dir)
    if exclude is not None:
        skip = Path(exclude).expanduser().resolve()
        files = [p for p in files if p.resolve() != skip]
    per_task: dict[str, dict[str, Any]] = {}
    used_files: list[str] = []
    for path in files:
        records = _bm.read_run(path)
        touched = False
        for rec in records:
            tid = str(rec.get("task_id") or "")
            if not tid:
                continue
            bucket = per_task.setdefault(
                tid, {"n_files": 0, "samples": []})
            if bucket["n_files"] >= k:
                continue
            bucket["n_files"] += 1
            bucket["samples"].extend(_expand_samples(rec))
            touched = True
        if touched:
            used_files.append(str(path))
        if per_task and all(b["n_files"] >= k for b in per_task.values()):
            # Every known task has a full window; later (older) files
            # can only add tasks not seen in any newer run, which would
            # be stale baselines — stop here.
            break

    tasks: dict[str, dict[str, Any]] = {}
    for tid, bucket in per_task.items():
        samples = bucket["samples"]
        if not samples:
            continue
        qualities = [q for q, _, _ in samples]
        passes = sum(1 for _, s, _ in samples if s)
        n = len(samples)
        tasks[tid] = {
            "n_files": bucket["n_files"],
            "n_samples": n,
            # Strict majority: a 50/50 flaky task does NOT count as a
            # baseline pass, so its routine failures stay non-suspect.
            "majority_success": passes * 2 > n,
            "pass_fraction": passes / n,
            "median_quality": _median(qualities),
            "quality_spread": _stdev(qualities),
            "median_cost_usd": _median([c for _, _, c in samples]),
        }
    return {
        "model": str(model or ""),
        "runs_dir": str(_runs_dir(runs_dir)),
        "last_k": k,
        "files": used_files,
        "tasks": tasks,
    }


# ---------------------------------------------------------------------------
# Stage 2: classify a run against the baseline
# ---------------------------------------------------------------------------

STATUS_STABLE = "stable"
STATUS_IMPROVED = "improved"
STATUS_SUSPECT = "suspect_regression"
STATUS_NEW = "new_task"


def _classify(record: dict, baseline: dict | None) -> dict:
    """Classification entry for one task of the current run."""
    quality = float(record.get("quality_0_100") or 0)
    success = bool(record.get("success"))
    entry: dict[str, Any] = {
        "task_id": str(record.get("task_id") or ""),
        "task_class": str(record.get("task_class") or ""),
        "success": success,
        "quality": quality,
        "cost_usd": float(record.get("cost_usd") or 0.0),
        "error": str(record.get("error") or ""),
        "baseline": baseline,
        "status": STATUS_NEW,
        "reason": "",
        "threshold": QUALITY_DROP_FLOOR,
        "delta_quality": None,
    }
    if not baseline or not baseline.get("n_samples"):
        entry["reason"] = "no baseline for this task"
        return entry
    thr = suspect_threshold(baseline.get("quality_spread", 0.0))
    median = float(baseline.get("median_quality") or 0.0)
    delta = quality - median
    entry["threshold"] = thr
    entry["delta_quality"] = delta
    if baseline.get("majority_success") and not success:
        entry["status"] = STATUS_SUSPECT
        entry["reason"] = "fail_flip"
    elif delta <= -thr:
        entry["status"] = STATUS_SUSPECT
        entry["reason"] = "quality_drop"
    elif (success and not baseline.get("majority_success")) or delta >= thr:
        entry["status"] = STATUS_IMPROVED
        entry["reason"] = ("pass_flip"
                           if success and not baseline.get("majority_success")
                           else "quality_gain")
    else:
        entry["status"] = STATUS_STABLE
    return entry


def compare_run(
    run_path: Path | str | list[dict],
    *,
    history: dict,
) -> dict:
    """Classify every task of a run against the rolling baseline.

    ``run_path`` may be a JSONL run file (torn lines tolerated) or a
    pre-read record list.  Pure computation — no API spend.  Returns
    ``{"per_task": [...], "suspects": [...], "counts": {...},
    "run_summary": {...}, "n_tasks": int, "run_path": str}``.
    """
    if isinstance(run_path, (str, Path)):
        path = Path(run_path).expanduser()
        records = _bm.read_run(path)
        src = str(path)
    else:
        records = [r for r in (run_path or []) if isinstance(r, dict)]
        src = ""
    baselines = (history or {}).get("tasks") or {}
    per_task = [
        _classify(rec, baselines.get(str(rec.get("task_id") or "")))
        for rec in records
        if isinstance(rec, dict) and rec.get("task_id")
    ]
    counts = {s: 0 for s in
              (STATUS_STABLE, STATUS_IMPROVED, STATUS_SUSPECT, STATUS_NEW)}
    for e in per_task:
        counts[e["status"]] += 1
    return {
        "run_path": src,
        "n_tasks": len(per_task),
        "per_task": per_task,
        "suspects": [e for e in per_task if e["status"] == STATUS_SUSPECT],
        "counts": counts,
        "run_summary": _bm.summarise_run(records),
    }


# ---------------------------------------------------------------------------
# Stage 3: capped recheck of suspects
# ---------------------------------------------------------------------------


RerunFn = Callable[[str], dict]
"""task_id -> aggregated result record (dict shaped like a run line)."""


def _default_rerun(
    *, model: str, provider: str, backend: str, repeats: int,
    max_tokens: int = 1024,
) -> RerunFn:
    """LIVE rerun of one task via the existing runner machinery.

    Imports are deferred so offline tests that inject their own rerun
    never touch the engine stack.  ``run_task`` with ``repeats>1``
    already returns a majority-success / median-quality aggregate.
    """
    from dataclasses import asdict
    from .benchmark_runner import resolve_profile_name, run_task

    tasks = {t.id: t for t in _bm.load_tasks()}
    profile_name = resolve_profile_name(model)

    def _rerun(task_id: str) -> dict:
        task = tasks.get(task_id)
        if task is None:
            raise KeyError(f"unknown benchmark task: {task_id}")
        result = run_task(
            task,
            model=model, backend=backend or "api", provider=provider,
            profile_name=profile_name, max_tokens=max_tokens,
            repeats=max(1, int(repeats)),
        )
        return asdict(result)

    return _rerun


def _severity_key(entry: dict) -> tuple[int, float]:
    """Fail flips first, then by drop magnitude (largest first)."""
    flip = 0 if entry.get("reason") == "fail_flip" else 1
    delta = entry.get("delta_quality")
    try:
        drop = -float(delta) if delta is not None else 0.0
    except (TypeError, ValueError):
        drop = 0.0
    return (flip, -drop)


def _estimated_task_cost(entry: dict, repeats: int) -> float:
    baseline = entry.get("baseline") or {}
    per_sample = float(baseline.get("median_cost_usd") or 0.0)
    if per_sample <= 0:
        per_sample = float(entry.get("cost_usd") or 0.0)
    if per_sample <= 0:
        per_sample = FALLBACK_TASK_COST_USD
    return per_sample * max(1, int(repeats))


def recheck_suspects(
    suspects: Iterable[dict],
    *,
    model: str,
    provider: str = "",
    backend: str = "api",
    repeats: int = DEFAULT_RECHECK_REPEATS,
    max_tasks: int = DEFAULT_MAX_RECHECK_TASKS,
    max_cost_usd: float = DEFAULT_MAX_RECHECK_COST_USD,
    rerun: RerunFn | None = None,
) -> dict:
    """Re-run ONLY the suspect tasks and separate signal from noise.

    Each selected task is re-run with ``repeats`` fresh samples (via the
    existing ``run_task`` machinery unless ``rerun`` is injected); the
    aggregate's majority success + median quality decides
    ``confirmed_regression`` vs ``noise``:

      - a ``fail_flip`` suspect is confirmed when the recheck aggregate
        still fails;
      - a ``quality_drop`` suspect is confirmed when the recheck median
        quality is still below ``baseline_median - threshold`` (or the
        recheck fails a baseline-majority-pass task outright).

    Hard caps: at most ``max_tasks`` tasks (most severe first — fail
    flips, then largest drops) and a cost ceiling ``max_cost_usd``
    enforced BOTH on the pre-run estimate and on actual accumulated
    spend between tasks.  Everything not rechecked is listed in
    ``skipped`` with the cap that excluded it.  Never raises; a rerun
    exception becomes a ``recheck_error`` entry for that task.
    """
    entries = [dict(e) for e in (suspects or []) if isinstance(e, dict)]
    entries.sort(key=_severity_key)
    cap_n = max(0, int(max_tasks))
    budget = max(0.0, float(max_cost_usd))

    out: dict[str, Any] = {
        "model": str(model or ""),
        "repeats": max(1, int(repeats)),
        "caps": {"max_tasks": cap_n, "max_cost_usd": budget},
        "checked": [],
        "confirmed": [],
        "noise": [],
        "errors": [],
        "skipped": [],
        "spent_usd": 0.0,
    }
    if not entries:
        return out

    selected: list[dict] = []
    projected = 0.0
    for i, entry in enumerate(entries):
        tid = str(entry.get("task_id") or "")
        if len(selected) >= cap_n:
            out["skipped"].append({
                "task_id": tid,
                "reason": f"task cap reached (max_tasks={cap_n})",
            })
            continue
        est = _estimated_task_cost(entry, repeats)
        if projected + est > budget:
            out["skipped"].append({
                "task_id": tid,
                "reason": (f"cost cap reached (estimate "
                           f"${projected + est:.2f} > ${budget:.2f})"),
            })
            continue
        projected += est
        selected.append(entry)

    if not selected:
        return out

    runner = rerun
    if runner is None:
        try:
            runner = _default_rerun(
                model=model, provider=provider, backend=backend,
                repeats=repeats)
        except Exception as exc:
            out["errors"].append(f"recheck runner init failed: {exc}")
            for entry in selected:
                out["skipped"].append({
                    "task_id": str(entry.get("task_id") or ""),
                    "reason": "recheck runner unavailable",
                })
            return out

    spent = 0.0
    for entry in selected:
        tid = str(entry.get("task_id") or "")
        if spent >= budget:
            out["skipped"].append({
                "task_id": tid,
                "reason": (f"cost cap reached mid-recheck "
                           f"(spent ${spent:.2f} of ${budget:.2f})"),
            })
            continue
        try:
            rec = runner(tid) or {}
        except Exception as exc:
            out["errors"].append(f"{tid}: recheck failed: {exc}")
            out["checked"].append({
                "task_id": tid, "verdict": "recheck_error",
                "detail": str(exc)[:200],
            })
            continue
        cost = float(rec.get("cost_usd") or 0.0)
        spent += cost
        success = bool(rec.get("success"))
        quality = float(rec.get("quality_0_100") or 0)
        baseline = entry.get("baseline") or {}
        thr = float(entry.get("threshold") or QUALITY_DROP_FLOOR)
        median = float(baseline.get("median_quality") or 0.0)
        still_failing = bool(baseline.get("majority_success")) and not success
        still_low = quality <= median - thr
        confirmed = still_failing or still_low
        verdict = "confirmed_regression" if confirmed else "noise"
        checked = {
            "task_id": tid,
            "verdict": verdict,
            "suspect_reason": str(entry.get("reason") or ""),
            "recheck_success": success,
            "recheck_quality": quality,
            "recheck_n_samples": int(rec.get("n_samples") or 1),
            "baseline_median_quality": median,
            "threshold": thr,
            "cost_usd": cost,
        }
        out["checked"].append(checked)
        (out["confirmed"] if confirmed else out["noise"]).append(checked)
    out["spent_usd"] = spent
    return out


# ---------------------------------------------------------------------------
# Full unattended cycle
# ---------------------------------------------------------------------------


SuiteFn = Callable[..., Path]
"""(model, provider, backend, task_ids, repeats, runs_dir) -> run file."""


def _live_suite(
    model: str, provider: str, backend: str,
    task_ids: Optional[list[str]], repeats: int,
    runs_dir: Path | str | None,
    max_tokens: int = 1024,
) -> Path:
    """LIVE suite execution via the existing runner (deferred imports)."""
    from .benchmark_runner import resolve_profile_name, run_suite

    tasks = _bm.load_tasks()
    if task_ids:
        wanted = {str(t) for t in task_ids}
        tasks = [t for t in tasks if t.id in wanted]
    if not tasks:
        raise RuntimeError("no benchmark tasks matched")
    results = run_suite(
        tasks,
        model=model, backend=backend or "api", provider=provider,
        profile_name=resolve_profile_name(model),
        max_tokens=max_tokens, repeats=max(1, int(repeats)),
    )
    d = _runs_dir(runs_dir)
    return _bm.write_run(results, model=model, runs_dir=d)


def nightly(
    model: str,
    provider: str = "",
    backend: str = "api",
    *,
    task_ids: Optional[list[str]] = None,
    repeats: int = 1,
    recheck: bool = True,
    last_k: int = DEFAULT_LAST_K,
    runs_dir: Path | str | None = None,
    recheck_repeats: int = DEFAULT_RECHECK_REPEATS,
    max_recheck_tasks: int = DEFAULT_MAX_RECHECK_TASKS,
    max_recheck_cost_usd: float = DEFAULT_MAX_RECHECK_COST_USD,
    run_suite_fn: SuiteFn | None = None,
    rerun: RerunFn | None = None,
    emit: Callable[..., str] | None = None,
) -> dict:
    """One unattended benchmark cycle.  Never raises.

    run suite → compare vs rolling history → recheck suspects →
    markdown report under ``<runs_dir>/reports/<ts>_<model>.md`` →
    ONE ``run_failed`` attention event iff regressions were CONFIRMED
    by the recheck.  Injectables (``run_suite_fn``, ``rerun``,
    ``emit``) exist so tests and dry-runs never spend API money.
    Partial failures land in ``summary["errors"]``.
    """
    ts = time.time()
    summary: dict[str, Any] = {
        "ok": True,
        "model": str(model or ""),
        "provider": str(provider or ""),
        "backend": str(backend or "api"),
        "repeats": max(1, int(repeats)),
        "recheck_enabled": bool(recheck),
        "ts": ts,
        "run_path": "",
        "report_path": "",
        "history_files": [],
        "comparison": None,
        "recheck": None,
        "confirmed": [],
        "attention_id": "",
        "errors": [],
    }

    # 1. Baseline BEFORE the new run lands, so it never includes itself.
    history: dict = {"tasks": {}}
    try:
        history = load_history(model, runs_dir=runs_dir, last_k=last_k)
        summary["history_files"] = list(history.get("files") or [])
    except Exception as exc:
        summary["errors"].append(f"load_history failed: {exc}")

    # 2. Run the suite.
    run_path: Path | None = None
    try:
        suite = run_suite_fn or _live_suite
        run_path = Path(suite(
            summary["model"], summary["provider"], summary["backend"],
            list(task_ids) if task_ids else None,
            summary["repeats"], runs_dir,
        ))
        summary["run_path"] = str(run_path)
    except Exception as exc:
        summary["ok"] = False
        summary["errors"].append(f"suite run failed: {exc}")

    # 3. Compare against the baseline.
    suspects: list[dict] = []
    if run_path is not None:
        try:
            comparison = compare_run(run_path, history=history)
            summary["comparison"] = comparison
            suspects = comparison.get("suspects") or []
        except Exception as exc:
            summary["ok"] = False
            summary["errors"].append(f"compare failed: {exc}")

    # 4. Recheck suspects (skipped when disabled or nothing suspicious).
    if recheck and suspects:
        try:
            rc = recheck_suspects(
                suspects,
                model=summary["model"], provider=summary["provider"],
                backend=summary["backend"], repeats=recheck_repeats,
                max_tasks=max_recheck_tasks,
                max_cost_usd=max_recheck_cost_usd, rerun=rerun,
            )
            summary["recheck"] = rc
            summary["confirmed"] = list(rc.get("confirmed") or [])
            for err in rc.get("errors") or []:
                summary["errors"].append(f"recheck: {err}")
        except Exception as exc:
            summary["ok"] = False
            summary["errors"].append(f"recheck failed: {exc}")

    # 5. Markdown report (best-effort).
    try:
        reports = _runs_dir(runs_dir) / "reports"
        reports.mkdir(parents=True, exist_ok=True)
        report_path = reports / f"{int(ts)}_{_bm._slug(summary['model'])}.md"
        report_path.write_text(format_report(summary), encoding="utf-8")
        summary["report_path"] = str(report_path)
    except Exception as exc:
        summary["errors"].append(f"report write failed: {exc}")

    # 6. ONE attention event, only for CONFIRMED regressions.
    if summary["confirmed"]:
        try:
            emitter = emit
            if emitter is None:
                from . import attention as _att
                emitter = _att.emit_attention
            lines = [
                (f"{c['task_id']}: q={c['recheck_quality']:.0f} vs baseline "
                 f"median {c['baseline_median_quality']:.0f} "
                 f"({c['suspect_reason']})")
                for c in summary["confirmed"]
            ]
            detail = "; ".join(lines)
            if summary["report_path"]:
                detail += f" — report: {summary['report_path']}"
            summary["attention_id"] = emitter(
                "run_failed",
                title=(f"Benchmark regression: {summary['model']} — "
                       f"{len(summary['confirmed'])} confirmed"),
                detail=detail,
            )
        except Exception as exc:
            summary["errors"].append(f"attention emit failed: {exc}")
    return summary


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------


def format_report(summary: dict) -> str:
    """Compact markdown report for one ``nightly`` summary."""
    s = summary or {}
    when = time.strftime(
        "%Y-%m-%d %H:%M", time.localtime(float(s.get("ts") or time.time())))
    lines = [
        f"# Benchmark watch — {s.get('model', '?')}",
        "",
        f"- when: {when}",
        f"- run file: {s.get('run_path') or '(suite run failed)'}",
        f"- baseline: last-k window over "
        f"{len(s.get('history_files') or [])} run file(s)",
        f"- status: {'OK' if s.get('ok') else 'DEGRADED'}",
    ]
    comparison = s.get("comparison") or {}
    counts = comparison.get("counts") or {}
    if comparison:
        rs = comparison.get("run_summary") or {}
        lines += [
            "",
            "## Run vs rolling baseline",
            "",
            f"- tasks: {comparison.get('n_tasks', 0)} "
            f"(pass {rs.get('n_pass', 0)}/{rs.get('n_tasks', 0)}, "
            f"avg quality {float(rs.get('avg_quality') or 0):.1f}, "
            f"cost ${float(rs.get('total_cost_usd') or 0):.2f})",
            f"- stable: {counts.get(STATUS_STABLE, 0)}   "
            f"improved: {counts.get(STATUS_IMPROVED, 0)}   "
            f"suspect: {counts.get(STATUS_SUSPECT, 0)}   "
            f"new: {counts.get(STATUS_NEW, 0)}",
        ]
        suspects = comparison.get("suspects") or []
        if suspects:
            lines += ["", "### Suspects", ""]
            for e in suspects:
                base = e.get("baseline") or {}
                lines.append(
                    f"- `{e['task_id']}` ({e.get('reason', '?')}): "
                    f"q={e.get('quality', 0):.0f} vs median "
                    f"{float(base.get('median_quality') or 0):.0f} "
                    f"(threshold {float(e.get('threshold') or 0):.1f}, "
                    f"{int(base.get('n_samples') or 0)} baseline samples)")
    rc = s.get("recheck") or {}
    if rc:
        lines += [
            "",
            "## Recheck",
            "",
            f"- repeats: {rc.get('repeats', 0)}   "
            f"caps: {rc.get('caps', {}).get('max_tasks', 0)} tasks / "
            f"${float(rc.get('caps', {}).get('max_cost_usd') or 0):.2f}   "
            f"spent: ${float(rc.get('spent_usd') or 0):.2f}",
            f"- confirmed: {len(rc.get('confirmed') or [])}   "
            f"noise: {len(rc.get('noise') or [])}   "
            f"errors: {len(rc.get('errors') or [])}",
        ]
        for c in rc.get("checked") or []:
            lines.append(
                f"- `{c['task_id']}` → {c.get('verdict', '?')} "
                f"(recheck q={float(c.get('recheck_quality') or 0):.0f}, "
                f"success={c.get('recheck_success')})")
        skipped = rc.get("skipped") or []
        if skipped:
            lines += ["", "### Skipped by cap", ""]
            for sk in skipped:
                lines.append(f"- `{sk.get('task_id', '?')}`: "
                             f"{sk.get('reason', '?')}")
    elif s.get("recheck_enabled") and counts.get(STATUS_SUSPECT):
        lines += ["", "## Recheck", "", "- not executed (see errors)"]
    elif not s.get("recheck_enabled"):
        lines += ["", "## Recheck", "",
                  "- disabled — suspects above are UNCONFIRMED single-sample "
                  "observations"]
    confirmed = s.get("confirmed") or []
    lines += [
        "",
        "## Verdict",
        "",
        (f"- {len(confirmed)} CONFIRMED regression(s): "
         + ", ".join(f"`{c['task_id']}`" for c in confirmed))
        if confirmed else
        "- no confirmed regressions",
    ]
    errors = s.get("errors") or []
    if errors:
        lines += ["", "## Errors", ""]
        lines += [f"- {e}" for e in errors]
    return "\n".join(lines) + "\n"


__all__ = [
    "DEFAULT_LAST_K",
    "DEFAULT_MAX_RECHECK_COST_USD",
    "DEFAULT_MAX_RECHECK_TASKS",
    "DEFAULT_RECHECK_REPEATS",
    "FALLBACK_TASK_COST_USD",
    "QUALITY_DROP_FLOOR",
    "SPREAD_FACTOR",
    "compare_run",
    "format_report",
    "list_model_runs",
    "load_history",
    "nightly",
    "recheck_suspects",
    "suspect_threshold",
]
