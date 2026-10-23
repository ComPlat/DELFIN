# Nightly benchmark round (`delfin.agent.bench_nightly`)

One scheduled round measures the agents' work on DELFIN's own
capabilities against a fixed yardstick: the same holdout seed draws
the same instances every round, so a score change means the work
changed, not the sample.

## Round flow (`bench_nightly.run(config) -> RunRecord`)

1. **Seed gate.** In holdout mode the seed is read from a file whose
   *path* comes from the `holdout_seed_path` setting — the seed value
   itself never appears in the record, the report, or a task
   workspace. A missing, empty or unreadable seed file refuses to
   grade the holdout with a stated reason in
   `RunRecord.holdout_refusal_reason`; practice runs use the
   `public_seed` instead and log it.
2. **Sampling.** `config["sampler"].sample(seed, n)` draws the round's
   tasks (stage-1a `bench_holdout` once merged; any object with that
   method works). The same seed draws the same instances in the same
   order.
3. **Budget gate, in code.** For each task,
   `config["accounting"].estimate(task)` joins the day's spend *before*
   the instance starts. An instance is refused when its estimate is
   negative or when `spent + est` exceeds
   `config.get("daily_budget_core_hours", 50)`; an estimate exactly
   equal to the remaining budget is allowed. Every refusal carries a
   stated reason on the instance and counts in
   `RunRecord.refused_over_budget`.
4. **Run + grade.** A started instance charges EF1's ledger once
   (`accounting.charge(est)`), then runs through
   `benchmark_runner.run_task` (or the configured `runner`) and is
   graded by `config["grader"].grade(task, result)`.
   `InstanceRecord.error` is reserved for a failed run;
   accounting/sampling failures land in `skip_reason` so one broken
   stage-1a interface cannot kill the round.
5. **Submit path.** By default the runner executes in-process and the
   record says `submit_path="local"`. With `config["submitter"]` set —
   a callable `(task, staging_lines, job_name) -> job_id` — the job
   goes through DELFIN's own submit path with the node-local staging
   lines reused from
   `delfin.cluster_bench.slurm_script._cbatch_staging_lines`
   (`tool_python` in config; the lines are `[]` without a tool
   interpreter). The record carries `submit_path="slurm"` and the job
   id; a submitter that raises is recorded (`"slurm-failed"`) and the
   round continues.

## Report (`bench_report.write(record, history)`)

A pure function of the record and earlier rounds' gains; JSON with
sorted keys, reproducible byte-identically:

- score per tool family and per correctness level (consistency,
  reproduction, reference, trap), cost beside every score;
- `weakest_queue`: the tasks under 60 quality, weakest first, as the
  next round's queue;
- `stop_verdict`: stop after two consecutive rounds without positive
  gain, or when the daily budget is exhausted — with the reason.

## Scheduling

The runner is a pure function; the nightly schedule drives it with a
fresh config each round (yesterday's report feeds
`previous_rounds_gain`). Jobs larger than the calc policy's small-
system limits must go through the submit path, never raw `sbatch`.
