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
   stage-1a interface cannot kill the round. A result whose `unmeasured`
   flag is set (HB2: the run never happened — engine-init failure) is
   NOT graded: the instance keeps the error text, is excluded from
   every score, and the report names it.
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
- `unmeasured` (HB2): the run never happened (e.g. the engine failed
  to start: "engine init failed: ..."). Such a result is recorded
  with its error text on the instance (`InstanceRecord.unmeasured`,
  `InstanceRecord.error`), is NEVER graded — a dead run is not a
  quality-0 failure of the model — and is excluded from every score.
  The report names it under `unmeasured`; the CLI exits **3** when
  every started instance is unmeasured (a silent all-0 baseline like
  the one of 2026-10-31 can no longer look like a completed round).
  A history round in that shape is marked `history_invalid` in the
  report and its gains are excluded from the stop verdict.
- `weakest_queue`: the tasks under 60 quality, weakest first, as the
  next round's queue;
- `stop`: stop after two consecutive rounds without positive
  gain, or when the daily budget is exhausted — with the reason.

## Scheduling

The runner is a pure function; the nightly schedule drives it with a
fresh config each round (yesterday's report feeds
`previous_rounds_gain`). Jobs larger than the calc policy's small-
system limits must go through the submit path, never raw `sbatch`.

## Command-line entry point (HB1)

The one command the person runs for a nightly round:

    python -m delfin.agent.bench_nightly --model M --instances N \
        --seed-path PATH [--provider P] [--practice] [--budget H]

- `--seed-path PATH` is the *path* of the holdout seed file outside
  the repository (a setting, not the seed). The value is read into
  memory only — never printed, logged or written; the record says
  only whether the file was present.
- Without the seed file the holdout is refused: the refusal reason is
  printed on stderr and the exit status is 1. `--practice` uses the
  public seed instead and never consults the seed file.
- `--provider P` (HB2) names the provider backend (kit, anthropic,
  openai, ollama). Default: derived from the model id the same way
  DELFIN's CLI derives it (single source: `pricing.infer_provider`),
  so `--model kit.<id>` runs on kit without a flag. Without a usable
  provider the engine is built as Anthropic and every run dies
  unmeasured with an engine-init error (defect D1, fixed 2026-11).
- exit status: 0 a completed round; 1 the holdout was refused; 2 bad
  usage; 3 the round measured nothing (every started instance
  unmeasured).
- stdout carries the full `RunRecord` as JSON — the input
  `bench_report.write` consumes.

The round's wiring comes from `bench_holdout_glue`
(`build_default_config`): the CB0B sampler draws through
`bench_holdout.sample`, each instance is adapted to a
`benchmark.Task` (`instance_to_task`, unique draw ids
`<instance_id>.<n>`), the runner goes through
`benchmark_runner.run_task`, the grader layers DELFIN's own verdicts
(the runner's rubric result, the CB0A run-evidence check
`grade_run`, the SC1 claim check `result_claims.check` — extra layers
only lower), and the accounting estimate comes from DELFIN's own
resource logic (`workflow_grader._CORE_HOUR_ESTIMATE` per family
step, capped by `cluster_efficiency.advise` PAL vs
`calc_policy.LOGIN_LOCAL_CORES`).

A real holdout baseline is run by the person with the real seed file;
practice runs (`--practice`) measure the wiring on fakes and never
touch the seed.
