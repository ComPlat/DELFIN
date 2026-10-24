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
