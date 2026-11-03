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
    # HB2 D2: the run did not measure the model (e.g. the engine never
    # started: "engine init failed: No Anthropic API key found"). Such
    # an instance is not a failure of the model — it is no measurement
    # at all. Recorded with the error text, never graded, never scored.
    unmeasured: bool = False
    # Submit path (P3): how the instance ran. "local" (default,
    # in-process) or "slurm" when a submitter was configured.
    submit_path: str = "local"
    job_id: Optional[str] = None


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
    # What may be said about the seed: practice logs the public seed,
    # holdout only whether the seed file was present (QS ruling).
    seed_log: str = ""
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
    submitter = config.get("submitter")
    tool_python = config.get("tool_python")
    if submitter is not None:
        from ..cluster_bench.slurm_script import (
            _cbatch_staging_lines as staging_lines,
        )
    else:
        staging_lines = None
    sampler: Sampler = config["sampler"]
    grader: Grader = config["grader"]
    accounting: Accounting = config["accounting"]
    model: str = config.get("model", "")
    provider: str = config.get("provider", "")
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
        record.seed_log = seed
    elif seed_path and Path(seed_path).is_file():
        seed = _read_seed(seed_path)
        if not seed:
            # An existing-but-empty seed file is not a seed (s19 F3).
            seed = None
            record.seed_log = "holdout seed file empty"
            record.holdout_refusal_reason = (
                "holdout seed file is empty; refusing to grade the "
                "holdout on a meaningless seed")
        else:
            record.holdout_ran = True
            record.seed_log = "holdout seed file present"
    else:
        record.seed_log = "holdout seed file absent"
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
        try:
            est = float(accounting.estimate(task))
        except Exception as exc:
            # A broken stage-1a interface is recorded on the instance,
            # not allowed to kill the round (s19 F1).
            inst.skip_reason = ("accounting estimate failed: %s: %s"
                                % (type(exc).__name__, exc))
            record.instances.append(inst)
            continue
        if est < 0:
            # A negative estimate must never reduce the day's spend —
            # it would let later instances fit for free (s19 F4).
            inst.skip_reason = ("negative resource estimate %.3f refused"
                                % est)
            record.refused_over_budget += 1
            record.instances.append(inst)
            continue
        # Budget in code: the estimate joins the day's spend BEFORE the
        # instance starts; == the limit fits, above it is refused.
        if spent + est > budget:
            inst.skip_reason = (
                "over daily budget: estimate %.3f + spent %.3f exceeds "
                "limit %.3f core-hours" % (est, spent, budget))
            record.refused_over_budget += 1
            record.instances.append(inst)
            continue
        try:
            accounting.charge(est)
        except Exception as exc:
            # Same rule as the estimate: recorded, not fatal (s19 F2).
            inst.skip_reason = ("accounting charge failed: %s: %s"
                                % (type(exc).__name__, exc))
            record.refused_over_budget += 1
            record.instances.append(inst)
            continue
        spent += est
        # EF1's ledger was already charged above (once per started
        # instance, s19 F5) — the spent counter only mirrors it.
        inst.core_hours = est
        # P3: jobs go through DELFIN's own submit path when a submitter
        # is configured; the staging lines come from the REAL
        # _cbatch_staging_lines (reused, never copied) and show the
        # no-tool_python edge as [].
        if submitter is not None:
            inst.submit_path = "slurm"
            try:
                inst.job_id = submitter(
                    task,
                    staging_lines(tool_python),
                    "delfin-bench-%s-%s" % (record.mode, task.id))
            except Exception as exc:
                # A broken submitter is recorded, not fatal: the round
                # survives a failed stage interface (same rule as F1/F2).
                inst.submit_path = "slurm-failed"
                inst.error = ("submitter failed: %s: %s"
                              % (type(exc).__name__, exc))
                record.instances.append(inst)
                continue
            # F7: ran means "started", so only a SUBMITTED instance
            # counts as ran — a submission failure never ran.
            inst.ran = True
        else:
            inst.ran = True
        try:
            t0 = now()
            result = runner(task, model=model, provider=provider)
            inst.duration_s = now() - t0
            # HB2 D2: a result that says the run never happened is
            # recorded as such — error text kept, never graded (a dead
            # run is not a quality-0 failure of the model) and never
            # counted in the score (bench_report._score excludes it).
            unmeasured = bool(getattr(result, "unmeasured", False))
            if unmeasured:
                inst.unmeasured = True
                inst.error = str(getattr(result, "error", "") or
                                 "unmeasured (no error text)")
            else:
                inst.success, inst.quality = grader.grade(task, result)
                inst.error = str(getattr(result, "error", "") or "")
        except Exception as exc:  # a failed run is recorded, not swallowed
            inst.error = "%s: %s" % (type(exc).__name__, exc)
        record.instances.append(inst)
    record.spent_core_hours = spent
    record.ended_at = now()
    return record


# ---------------------------------------------------------------------------
# CLI entry point (HB1 phase 3)
# ---------------------------------------------------------------------------

def main(argv: Optional[list] = None) -> int:
    """``python -m delfin.agent.bench_nightly`` — the one entry point
    the operator runs for a nightly round (HB1).

    ::

        --model M            the model the benchmark runs against
        --instances N        how many tasks to draw (default 10)
        --seed-path PATH     the holdout seed file's PATH (a setting);
                             the VALUE is read into memory, never
                             printed, logged or written
        --practice           practice mode on the public seed; the
                             seed file is not consulted
        --budget H           daily core-hour budget (default 50)

    Without the seed file the holdout is REFUSED (run()'s seed gate)
    and the CLI exits non-zero, printing the refusal reason. Practice
    runs never touch the seed file.
    """
    import argparse
    import json
    import sys

    parser = argparse.ArgumentParser(
        prog="python -m delfin.agent.bench_nightly",
        description="One nightly benchmark round: sample, budget-gate, "
                    "run through benchmark_runner.run_task, grade with "
                    "DELFIN's own grading, record the cost.")
    parser.add_argument("--model", default="")
    parser.add_argument("--provider", default="",
                        help="provider backend for the model (kit, "
                             "anthropic, openai, ollama); default: "
                             "derived from the model id like the bench "
                             "CLI does")
    parser.add_argument("--instances", type=int, default=10)
    parser.add_argument("--seed-path", default="")
    parser.add_argument("--practice", action="store_true")
    parser.add_argument("--budget", type=float,
                        default=DEFAULT_DAILY_BUDGET_CORE_HOURS)
    args = parser.parse_args(argv)

    if not args.model:
        print("error: --model is required (the model the benchmark "
              "runs against)", file=sys.stderr)
        return 2

    from .bench_holdout_glue import build_default_config

    cfg = build_default_config(
        model=args.model,
        instances=args.instances,
        mode="practice" if args.practice else "holdout",
        holdout_seed_path=args.seed_path,
        provider=args.provider,
    )
    cfg["daily_budget_core_hours"] = args.budget
    record = run(cfg)
    print(json.dumps(record.to_dict(), indent=2, sort_keys=True,
                     default=str))
    if record.holdout_refusal_reason:
        print("refused: %s" % record.holdout_refusal_reason,
              file=sys.stderr)
        return 1
    # HB2 D2: a round where EVERY drawn instance is unmeasured did not
    # measure anything (the 31.10. defect: a dead baseline that looked
    # like a failed round). Exiting 0 would let a later round read it
    # as a real baseline; a healthy round still exits 0.
    started = record.started()
    if started and all(i.unmeasured for i in started):
        print("error: every drawn instance is UNMEASURED (the runs "
              "never happened; e.g. engine init failed) — this round "
              "measured nothing and is not a baseline", file=sys.stderr)
        for inst in started:
            print("  unmeasured: %s: %s" % (inst.task_id, inst.error),
                  file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":  # pragma: no cover — exercised via main()
    raise SystemExit(main())
