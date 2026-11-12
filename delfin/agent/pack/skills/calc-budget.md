---
domains: code
name: calc-budget
description: Estimate a job's core-hours, memory, wall time and queue wait from system size and method before submitting, run it against the session budget (go, or ask the person), report the wait honestly as "unknown" when the scheduler throttle is closed, and let the estimate learn from finished jobs.
---
# Calculate and budget a cluster job before submitting

> Prevent runaway calculations: estimate first, run it against the session
> budget, ask the person when over -- never silently shrink PAL/maxcore/wall
> time to make a job fit. DELFIN's own calculator is `delfin.agent.calc_budget`.

## Before submit, always estimate
1. Size from the system: atoms, electrons, basis functions; the method (a
   large DLPNO/CCSD(T) is far heavier than a DFT single point); the partition.
2. `estimate(job)` -> core-hours, memory MB, wall time, and the expected
   queue wait. The wait comes from the scheduler's OWN `--test-only` answer
   through the ONE throttled scheduler client -- never a per-loop or per-widget
   `squeue`/`sacct`. When the gate is closed (THROTTLED) or the answer states
   no start time, report the wait as `"unknown"`; never invent a number.
3. Compare against the session `Budget` remaining core-hours:
   within budget -> go; over budget -> ask the person with the estimate.
   The agent never shrinks PAL/maxcore/wall time silently to fit.

## After the job, let the estimate learn
`record_outcome(account_record, estimate)` -- where `account_record` is one
record from `cluster_efficiency.account()` (a real `sacct` row, never a
hand-typed value) -- appends actual/estimate per method family so the next
estimate for the same family is corrected. A job with no measured core-hours
is skipped, never treated as a zero.

## Honesty rules
- A wait, memory or core-hours figure has its unit and source (model, or the
  scheduler's own answer, or a learned factor over N finished jobs).
- An unmeasured or throttled answer is `"unknown"`, never `0`.
- Do not tune PAL/maxcore/time to hit a budget; ask the person instead.
