# SA1 explore pilot — measurement harness and guard

## What this is

`delfin/agent/explore_pilot.py` measures whether ONE narrow use of the
existing read-only `explore` sub-agent pays off: the same small set of
read-heavy tasks ("find where DELFIN decides X, answer with file:line")
run twice — once plain, once allowed to call the explore sub-agent —
recording tokens, wall time, tool calls and citation correctness per
arm. The helper becomes a recommendation (the skill
`delfin/agent/pack/skills/explore-helper.md`) ONLY if the measurement
says it saves tokens or time at equal correctness; otherwise the PR
ships the harness and the negative result.

## The pilot setting (the guard)

`check_delegate` allows exactly:
- the `explore` sub-agent type (read-only preset, subagents.py:1260-1276),
- at most one running per session (the claimed slot, held until
  `reap_slot`, covering the background channel too),
- depth 0 only (no nesting — depth is a property of the slot; the
  delegate's own runner is built with depth+1 baked in via
  `child_runner`),
- no worktree isolation (the pilot runs read-only in the harness's own
  workspace).

Everything else is refused with the reason stated; a refused call never
reaches `run_subagent`. The wrapper never raises, matching the real
runner's contract (subagents.py:3103).

Enforcement point: for the pilot the harness installs this guard AS the
`subagent_runner` it hands the engine (the real call path,
engine.py:961-992 is where the engine binds it). No protected-file
change is part of the pilot.

## The measurement

- `run_pilot(tasks, runner, root)`: two arms per task, per arm records
  success, input/output tokens, tool calls, wall time, final text.
  In production `runner` is `benchmark_runner.run_task` (protected);
  the tests pin the shape with a double, never a live model.
- `verify_citations(text, root, expect)`: every `path:line` claim is
  checked against the real filesystem — file exists under root, line
  in range, optional content expectation. Violations are named.
- `pilot_verdict(rows)`: explore wins only with strictly fewer tool
  calls AND no more wall time at equal correctness across all tasks;
  mixed evidence leaves the verdict None (a win is never claimed from
  noisy data).

## Tests

- `tests/test_sa1_guard.py` (17 cases): the pilot setting, both
  directions; types enumerated from the real `SUBAGENT_PRESETS` dict.
- `tests/test_sa1_harness.py` (9 cases): harness + citation oracle.
- `tests/test_sa1_adversarial_guard.py` (15 cases, reviewer nacht-s17):
  real-path nesting, crash release, look-alike types, busy refusal.

All 41 green through the gate at bb3a8779.

## Status

The guard and the harness are built and tested. The live measurement
(benchmark runs against real models) is run by the operator, per wave
rules — sessions never benchmark against real models. Until those
numbers exist, the skill ships as a recommendation conditional on the
measurement, and the harness is the deliverable either way.

# SP1 — the measurement runs end to end

SP1 (wave 22, builder nacht-s13) closed the two gaps that stopped the
operator from running the live measurement: there was no runner
adapter and no task file.

## The runner adapter

`make_pilot_runner` (`explore_pilot.py:352`) wraps the protected
production runner `benchmark_runner.run_task` (`benchmark_runner.py:1450`)
into the harness contract `runner(task, arm) -> dict` with `success`,
`input_tokens`, `output_tokens`, `tool_calls`, `duration_s`,
`final_text`. Two correctness points:

1. **Full text, not the excerpt.** `BenchmarkResult` carries only
   `text_excerpt`, truncated to 400 characters on success
   (`benchmark.py:2012`). An answer graded through the excerpt is not
   a citation measurement, so the adapter captures the full model
   output through a `run_once` wrapper around `cli._run_once`
   (`cli.py:332`).
2. **The arms differ only in the delegation permission.** The explore
   arm's engine is built by `_engine_factory_builder` and gets
   `permissions.subagent_runner = run_pilot_subagent` with a fresh
   pilot slot (one at a time, depth 0, nesting refused on the real
   path — depth is `max(call depth, slot.depth)`,
   `explore_pilot.py:141-146`). The plain arm's engine stays unwired:
   no `subagent_runner` binding, no delegation (`engine.py:972-974`).
   Same model, same task, same budgets otherwise.

HB1 (bench-holdout glue, builder nacht-s20) shares the adapter shape
`runner(task, arm=...) -> dict` and reuses the full-text capture
(agreed 2026-10-25).

## The tasks

`delfin/agent/pack/benchmark/pilot_tasks_sp1.yaml`: six code-reading
questions about DELFIN itself (the pilot guard, CONTROL parsing, the
benchmark runner entry, the harness, the scoring, the wiring). Every
prompt REQUIRES `path:line` citations in the form `file.py:LINE` —
without that requirement in the prompt the arms are compared at an
unequal bar. The file is deliberately NOT named `tasks_*.yaml`, so
`benchmark.load_tasks`' default glob (`benchmark.py:688`) never pulls
the pilot suite into the production benchmark; it is loaded explicitly
(`benchmark.py:684-685`). A test re-checks every asked fact against
the current tree, so a drifted task fails before a measurement, not
during one.

## Running it

Dry run (fakes, no model call, no network — checks the plumbing and
prints the report shape):

```bash
python -m delfin.agent.explore_pilot --model MODEL --dry-run
```

Live measurement (operator, real model endpoint):

```bash
python -m delfin.agent.explore_pilot --model MODEL \
    --tasks delfin/agent/pack/benchmark/pilot_tasks_sp1.yaml
```

A live run without `--tasks` is refused (exit 2): the pilot never
silently falls back to the production benchmark suite. The report is
JSON with the per-task per-arm numbers and the overall verdict
(`explore_wins`, `equal_correctness`, `detail`).

## SP1 tests

- `tests/test_sp1_runner_adapter.py` (7): adapter contract on fakes.
- `tests/test_sp1_review_findings.py` (6): review round 1 — the
  citation bar (`claims > 0`), nesting refused on the real path,
  refusal paths return the real `SubagentResult`.
- `tests/test_sp1_pilot_tasks.py` (5): the task file loads through the
  real loader, prompts demand citations, every task answerable against
  the current tree.
- `tests/test_sp1_cli_entry.py` (5): dry run, `--model` required, live
  run refuses without `--tasks`, live path through the real functions.

All 23 green through the gate at 382af5d0. The live measurement is
the operator's run.
