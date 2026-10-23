"""Nightly holdout runner (PL4): the schedule-side half of the
self-optimising loop.

Stage 1a builds the parts — CB0A run evidence, CB0B holdout sampling,
SC1 correctness gate, EF1 cluster accounting — but nothing runs them on
a schedule, keeps a score per tool family, or stops the loop. This
module is that runner and the report writer (:mod:`bench_report`).

Design contracts, from the PL4 brief:

* ``run(config) -> RunRecord`` draws instances through a sampler
  interface (CB0B's ``bench_holdout.sample(seed, n)`` once merged; any
  object with ``sample(seed, n)`` works here), runs each through the
  EXISTING benchmark runner (``benchmark_runner.run_task`` — protected,
  reused, never changed), grades with a grader interface, and records
  cost through an accounting interface.
* The holdout seed is read from a file OUTSIDE the repository whose
  path comes from a setting. It is never printed, logged, written into
  a report or a task workspace. Without the file the runner refuses to
  grade the holdout and says so; practice runs use the public seed.
* Budget is in code, not in prompts: before each instance the runner
  adds the instance's resource estimate to the day's spend and refuses
  to start it above the daily limit (default 50 core-hours).

Nothing here imports the engine stack at module level; the actual
runner function is injected in ``config["runner"]`` (defaults to
``benchmark_runner.run_task`` — imported lazily so tests parse fast
and the engine only loads when a real run starts).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional, Protocol

DEFAULT_DAILY_BUDGET_CORE_HOURS = 50.0


class Sampler(Protocol):
    """CB0B's bench_holdout.sample(seed, n) once merged."""

    def sample(self, seed: str, n: int): ...


class Grader(Protocol):
    """CB0A evidence + SC1 claims once merged: returns (success, quality)."""

    def grade(self, task, result) -> tuple[bool, int]: ...


class Accounting(Protocol):
    """EF1's cluster accounting: per-instance core-hour estimate + spend."""

    def estimate(self, task) -> float: ...
    def charge(self, core_hours: float) -> None: ...


@dataclass
class InstanceRecord:
    """One drawn task: whether it ran, why not, and its graded outcome."""

    task_id: str
    task_class: str
    ran: bool = False
    skip_reason: str = ""
    success: Optional[bool] = None
    quality: Optional[int] = None
    core_hours: float = 0.0
    duration_s: float = 0.0
    error: str = ""


@dataclass
class RunRecord:
    """Everything the report needs, and nothing the seed needs hiding in."""

    model: str = ""
    started_at: float = 0.0
    ended_at: float = 0.0
    mode: str = "holdout"
    instances: list[InstanceRecord] = field(default_factory=list)
    # Budget bookkeeping (EF1 interface, core-hours)
    daily_budget_core_hours: float = DEFAULT_DAILY_BUDGET_CORE_HOURS
    spent_core_hours: float = 0.0
    refused_over_budget: int = 0
    # Holdout gate
    holdout_ran: bool = False
    holdout_refusal_reason: str = ""
    # Flat-round bookkeeping for the stop verdict
    previous_rounds_gain: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def started(self) -> list[InstanceRecord]:
        return [i for i in self.instances if i.ran]


def _read_seed(seed_path: str) -> str:
    """Read the holdout seed from its file. The value stays in memory;
    no caller of this module prints or stores it."""
    with open(seed_path, "r", encoding="utf-8") as fh:
        return fh.read().strip()


def run(config: dict, *, clock: Callable[[], float] | None = None) -> RunRecord:
    """One nightly round: sample -> budget-gate -> run -> grade -> charge.

    ``config`` keys:
      ``model``                    the model the benchmark runs against
      ``runner``                   the task runner (default
                                   ``benchmark_runner.run_task``); a
                                   callable ``(task, **kwargs) -> result``
      ``sampler`` / ``grader`` /   the three stage-1a interfaces
      ``accounting``
      ``instances``                how many tasks to draw (default 10)
      ``daily_budget_core_hours``  the day's limit (default 50)
      ``holdout_seed_path``        PATH of the seed file outside the
                                   repository (a setting, not the seed)
      ``mode``                     "holdout" (default) or "practice"
      ``public_seed``              seed for practice runs
      ``previous_rounds_gain``     gains of earlier rounds, for the report

    The seed VALUE never enters the record: the sampler receives it as
    an argument, and ``InstanceRecord``/``RunRecord`` carry task ids and
    outcomes only.
    """
    now = clock or (lambda: __import__("time").time())
    runner = config.get("runner")
    if runner is None:
        from .benchmark_runner import run_task as runner
    sampler: Sampler = config["sampler"]
    grader: Grader = config["grader"]
    accounting: Accounting = config["accounting"]
    model: str = config.get("model", "")
    mode: str = config.get("mode", "holdout")
    n_instances = int(config.get("instances", 10))
    budget = float(config.get("daily_budget_core_hours",
                              DEFAULT_DAILY_BUDGET_CORE_HOURS))

    record = RunRecord(model=model, started_at=now(), mode=mode,
                       daily_budget_core_hours=budget,
                       previous_rounds_gain=list(
                           config.get("previous_rounds_gain", ())))

    # --- holdout seed gate -------------------------------------------------
    seed: Optional[str] = None
    seed_path = config.get("holdout_seed_path")
    if mode == "practice":
        seed = config.get("public_seed", "practice")
    elif seed_path and Path(seed_path).is_file():
        seed = _read_seed(seed_path)
        record.holdout_ran = True
    else:
        record.holdout_refusal_reason = (
            "holdout seed file not found at the configured setting path; "
            "refusing to grade the holdout (use mode='practice' with a "
            "public seed for practice runs)")

    if seed is None:
        record.ended_at = now()
        return record

    tasks = list(sampler.sample(seed, n_instances))
    spent = 0.0
    for task in tasks:
        inst = InstanceRecord(task_id=task.id, task_class=task.task_class)
        est = float(accounting.estimate(task))
        # Budget in code: the estimate joins the day's spend BEFORE the
        # instance starts; == the limit fits, above it is refused.
        if spent + est > budget:
            inst.skip_reason = (
                "over daily budget: estimate %.3f + spent %.3f exceeds "
                "limit %.3f core-hours" % (est, spent, budget))
            record.refused_over_budget += 1
            record.instances.append(inst)
            continue
        spent += est
        accounting.charge(est)
        inst.core_hours = est
        inst.ran = True
        try:
            t0 = now()
            result = runner(task, model=model)
            inst.duration_s = now() - t0
            inst.success, inst.quality = grader.grade(task, result)
        except Exception as exc:  # a failed run is recorded, not swallowed
            inst.error = "%s: %s" % (type(exc).__name__, exc)
        record.instances.append(inst)
    record.spent_core_hours = spent
    record.ended_at = now()
    return record
