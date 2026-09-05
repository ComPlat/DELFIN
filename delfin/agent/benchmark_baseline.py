"""A committed reference standard for the office benchmark, and the gate
that decides whether a new run has fallen below it.

Why this is neither ``compare_runs`` nor the rolling history
------------------------------------------------------------
Three things compare benchmark runs and they answer different questions.

``bench_watch.load_history`` builds a rolling baseline from the newest K
run files on this machine and asks "is today unusual compared to the last
few days?". It catches a sudden break well. It cannot catch slow decay --
the standard follows the code down, so a fifth of a percent lost per week
never registers -- and it lives in ``~/.delfin``, so a fresh clone has no
baseline at all and nobody reviews it when it moves.

This module asks the other question: "is this still as good as the
standard we agreed on?". The standard is committed, so it survives a
clone, and changing it shows up in a diff where somebody has to vouch for
the new numbers.

``compare_runs`` answers "did this knob help?" between two arbitrary runs,
and it classifies any drop in quality as *worse*. That is the right
question for an A/B experiment and the wrong one for a ratchet, because a
sampled benchmark is noisy: on the measured reference the ambiguity task
scored 45, 93, 93 across three runs of unchanged code. A gate that fires
on any drop would fire on that, every time, until it was ignored.

So the gate is built around the noise instead of pretending it is absent:

* **Suite pass rate** is the primary signal. Six tasks times N samples is
  far more evidence than any single task provides, and a real regression
  in the document layer shows up across tasks.
* **A task that dies outright** is caught individually. A task whose
  reference rate is at least ``_SOLID_RATE`` failing *every* sample is not
  noise -- at the reference rate of 1.0 observed over three samples, that
  outcome does not occur unless something changed.
* **Everything in between is a warning**, not a failure. A single sample
  flipping on a task already known to be unstable is exactly what the
  reference records as normal.

Quality scores ride along as a reported number and never gate on their
own. They are a graded judgement of one sampled answer; the pass/fail
signal is what the rubric actually decided.

What the reference costs to produce
-----------------------------------
Every entry comes from real model calls. Re-measuring is a spend decision,
never a side effect of running the suite, so nothing here writes a
reference implicitly -- unlike a first run that silently becomes the
standard, which records whatever state the code happened to be in.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

# A task the reference saw succeed at least this often is expected to work.
# Below it, the task is on record as unstable and only counts through the
# suite-level rate.
_SOLID_RATE = 0.9

# How far the suite may sit below the reference before the gate calls it a
# regression, in standard deviations of the reference's own spread. So a
# reference measured with more samples tightens the gate by itself, rather
# than a fixed percentage that is too loose at N=18 and absurd at N=3.
#
# Three, not two, because of the cadence this gate is meant for. A nightly
# run at two sigma turns red on chance alone roughly once a month, and a
# gate that cries wolf monthly stops being read. At three that becomes
# about once a year. The cost of the wider tolerance is low: a task that
# breaks outright is caught by its own rule below, regardless of sigma.
_SUITE_SIGMA = 3.0

# Floor for the suite tolerance: with a reference rate at or near 1.0 the
# standard error collapses to 0 and any single failure would trip the
# gate. One sample's worth of the suite is the smallest honest tolerance.
_MIN_SUITE_SLACK = 1.0

SCHEMA_VERSION = 1


# ---------------------------------------------------------------------------
# The reference standard
# ---------------------------------------------------------------------------

@dataclass
class TaskReference:
    """What one task did over N samples of unchanged code.

    The observation is stored as two integers -- how many samples ran and
    how many passed -- not as a rate. A rate has to be rounded to be
    written down, and a reference that does not survive its own round trip
    cannot be compared against exactly.
    """

    task_id: str
    n_samples: int
    n_passed: int
    quality_mean: float
    quality_stdev: float
    per_run_success: list[bool] = field(default_factory=list)

    @property
    def success_rate(self) -> float:
        return (self.n_passed / self.n_samples) if self.n_samples else 0.0

    @property
    def is_solid(self) -> bool:
        """True when the reference expects this task to work every time."""
        return self.success_rate >= _SOLID_RATE

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_samples": int(self.n_samples),
            "n_passed": int(self.n_passed),
            "quality_mean": round(float(self.quality_mean), 2),
            "quality_stdev": round(float(self.quality_stdev), 2),
            "per_run_success": [bool(x) for x in self.per_run_success],
        }


@dataclass
class Baseline:
    """The committed standard for one model."""

    model: str
    mode: str
    measured_at: str
    commit: str
    repeats: int
    tasks: dict[str, TaskReference]
    note: str = ""

    @property
    def suite_pass_rate(self) -> float:
        """Fraction of all (task, sample) pairs the reference passed."""
        total = sum(t.n_samples for t in self.tasks.values())
        if not total:
            return 0.0
        return sum(t.n_passed for t in self.tasks.values()) / total

    @property
    def total_samples(self) -> int:
        return sum(t.n_samples for t in self.tasks.values())

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_VERSION,
            "model": self.model,
            "mode": self.mode,
            "measured_at": self.measured_at,
            "commit": self.commit,
            "repeats": self.repeats,
            "note": self.note,
            "suite_pass_rate": round(self.suite_pass_rate, 4),
            "tasks": {tid: t.to_dict() for tid, t in sorted(self.tasks.items())},
        }



def observed(row: dict) -> tuple[int, int]:
    """How many samples a run row represents, and how many passed.

    Reading this wrong is silent and total: a single-sample run carries
    ``success_rate`` as the dataclass default 0.0 -- the field exists but
    is only ever filled when repeats > 1 -- so treating "present" as
    "meaningful" scores a perfect run as zero and reports every task as
    dead. The recorded per-sample outcomes are the observation; the rate
    is a summary that only exists above one sample; the bare flag is what
    a single run actually decided.
    """
    n = int(row.get("n_samples") or 1)
    per_run = list(row.get("per_run_success") or [])
    if per_run:
        return len(per_run), sum(1 for x in per_run if x)
    if n > 1 and row.get("success_rate") is not None:
        return n, int(round(float(row.get("success_rate") or 0.0) * n))
    return n, (n if row.get("success") else 0)


def baseline_from_results(
    rows: Iterable[dict],
    *,
    measured_at: str,
    commit: str = "",
    note: str = "",
) -> Baseline:
    """Build a reference from the rows of a completed benchmark run.

    ``rows`` are the dicts a run writes -- one per task, already
    aggregated over its replicates by ``run_task``.
    """
    rows = list(rows)
    if not rows:
        raise ValueError("cannot build a reference from an empty run")

    tasks: dict[str, TaskReference] = {}
    for r in rows:
        tid = str(r.get("task_id") or "")
        if not tid:
            continue
        n, n_passed = observed(r)
        per_run = list(r.get("per_run_success") or [])
        tasks[tid] = TaskReference(
            task_id=tid,
            n_samples=n,
            n_passed=n_passed,
            quality_mean=float(r.get("quality_0_100") or 0),
            quality_stdev=float(r.get("quality_stdev") or 0.0),
            per_run_success=[bool(x) for x in per_run],
        )
    if not tasks:
        raise ValueError("no task rows in the run")

    first = rows[0]
    return Baseline(
        model=str(first.get("model") or ""),
        mode=str(first.get("mode") or ""),
        measured_at=measured_at,
        commit=commit,
        repeats=max((t.n_samples for t in tasks.values()), default=1),
        tasks=tasks,
        note=note,
    )


def load_baseline(path: Path) -> Baseline | None:
    """Read a committed reference. Returns None when there is none."""
    if not path.exists():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    schema = int(data.get("schema") or 0)
    if schema != SCHEMA_VERSION:
        raise ValueError(
            f"reference at {path} has schema {schema}, this build speaks "
            f"{SCHEMA_VERSION} -- re-measure rather than reinterpreting "
            "numbers whose meaning may have changed")
    tasks = {
        tid: TaskReference(
            task_id=tid,
            n_samples=int(t.get("n_samples") or 1),
            n_passed=int(t.get("n_passed") or 0),
            quality_mean=float(t.get("quality_mean") or 0.0),
            quality_stdev=float(t.get("quality_stdev") or 0.0),
            per_run_success=[bool(x) for x in (t.get("per_run_success") or [])],
        )
        for tid, t in (data.get("tasks") or {}).items()
    }
    return Baseline(
        model=str(data.get("model") or ""),
        mode=str(data.get("mode") or ""),
        measured_at=str(data.get("measured_at") or ""),
        commit=str(data.get("commit") or ""),
        repeats=int(data.get("repeats") or 1),
        tasks=tasks,
        note=str(data.get("note") or ""),
    )


def save_baseline(baseline: Baseline, path: Path) -> Path:
    """Write a reference. Only ever called from a deliberate re-measure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(baseline.to_dict(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------

def _suite_tolerance(baseline: Baseline) -> float:
    """How far below the reference rate a candidate may land, in rate units.

    The spread is summed per task rather than pooled across the suite.
    Pooling assumes every sample shares one rate, which would take the one
    genuinely unstable task and smear its noise over five rock-solid ones
    -- inflating the tolerance exactly where the reference is most certain,
    and so hiding real breakage there. Each task contributes its own
    n*p*(1-p): a task the reference saw pass every time contributes none.

    Floored so a perfect reference still tolerates a single sample; at
    p=1 the spread is zero and any one failure anywhere would be fatal.
    """
    n = baseline.total_samples
    if n <= 0:
        return 1.0
    variance = sum(
        t.n_samples * t.success_rate * (1.0 - t.success_rate)
        for t in baseline.tasks.values()
    )
    spread_in_samples = math.sqrt(max(variance, 0.0))
    slack = max(_SUITE_SIGMA * spread_in_samples, _MIN_SUITE_SLACK)
    return slack / n


def compare_to_baseline(rows: Iterable[dict], baseline: Baseline) -> dict[str, Any]:
    """Judge a candidate run against the committed reference.

    Returns ``{"verdict", "regressions", "warnings", "notes", "summary"}``
    where verdict is one of:

    ``pass``      -- at or above the reference, within noise
    ``regressed`` -- a solid task died outright, or the suite rate fell
                     further below the reference than sampling explains
    ``thin``      -- too little overlap with the reference to judge
    """
    rows = list(rows)
    by_id = {str(r.get("task_id") or ""): r for r in rows}
    overlap = sorted(set(by_id) & set(baseline.tasks))

    regressions: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []

    missing = sorted(set(baseline.tasks) - set(by_id))
    for tid in missing:
        regressions.append(
            f"{tid}: in the reference but absent from this run -- a task that "
            "stops being measured is not a task that passes")

    new_tasks = sorted(set(by_id) - set(baseline.tasks))
    for tid in new_tasks:
        notes.append(f"{tid}: not in the reference, so it is reported but not gated")

    if len(overlap) < 3:
        return {
            "verdict": "thin",
            "regressions": regressions,
            "warnings": warnings,
            "notes": notes + [
                f"only {len(overlap)} tasks overlap the reference; "
                "not enough to judge the suite"],
            "summary": {},
        }

    cand_passed = 0.0
    cand_samples = 0
    for tid in overlap:
        row = by_id[tid]
        ref = baseline.tasks[tid]
        n, n_passed = observed(row)
        rate = (n_passed / n) if n else 0.0
        cand_passed += n_passed
        cand_samples += n

        if ref.is_solid and rate <= 0.0:
            regressions.append(
                f"{tid}: every sample failed ({n}/{n}), and the reference has "
                f"it passing {ref.success_rate:.0%} of {ref.n_samples} -- "
                "not a sampling accident")
        elif rate < ref.success_rate:
            warnings.append(
                f"{tid}: {rate:.0%} vs {ref.success_rate:.0%} in the reference "
                f"(n={n} against {ref.n_samples})")

    cand_rate = cand_passed / cand_samples if cand_samples else 0.0
    tolerance = _suite_tolerance(baseline)
    floor = baseline.suite_pass_rate - tolerance
    if cand_rate < floor:
        regressions.append(
            f"suite pass rate {cand_rate:.1%} is below the reference "
            f"{baseline.suite_pass_rate:.1%} by more than sampling explains "
            f"(tolerance {tolerance:.1%}, floor {floor:.1%})")

    summary = {
        "n_overlap": len(overlap),
        "candidate_pass_rate": round(cand_rate, 4),
        "baseline_pass_rate": round(baseline.suite_pass_rate, 4),
        "tolerance": round(tolerance, 4),
        "floor": round(max(floor, 0.0), 4),
        "candidate_samples": cand_samples,
        "baseline_samples": baseline.total_samples,
    }
    return {
        "verdict": "regressed" if regressions else "pass",
        "regressions": regressions,
        "warnings": warnings,
        "notes": notes,
        "summary": summary,
    }


def format_baseline_report(result: dict[str, Any], baseline: Baseline) -> str:
    """Human-readable rendering of a gate result."""
    verdict = str(result.get("verdict") or "?")
    s = result.get("summary") or {}
    lines = [
        f"Office benchmark vs reference ({baseline.model}, "
        f"measured {baseline.measured_at}"
        + (f", commit {baseline.commit[:8]}" if baseline.commit else "")
        + ")",
        f"  verdict: {verdict.upper()}",
    ]
    if s:
        lines.append(
            f"  pass rate {s['candidate_pass_rate']:.1%} vs reference "
            f"{s['baseline_pass_rate']:.1%} "
            f"(floor {s['floor']:.1%}, {s['candidate_samples']} samples)")
    for r in result.get("regressions") or []:
        lines.append(f"  REGRESSION  {r}")
    for w in result.get("warnings") or []:
        lines.append(f"  warn        {w}")
    for n in result.get("notes") or []:
        lines.append(f"  note        {n}")
    if verdict == "regressed":
        lines.append(
            "  If this is an intended change, re-measure deliberately rather "
            "than editing the numbers by hand.")
    return "\n".join(lines)


__all__ = [
    "SCHEMA_VERSION",
    "observed",
    "TaskReference",
    "Baseline",
    "baseline_from_results",
    "load_baseline",
    "save_baseline",
    "compare_to_baseline",
    "format_baseline_report",
]
