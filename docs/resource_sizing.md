# Resource sizing — measured, not guessed

DELFIN's pre-submit advice (`delfin.agent.cluster_efficiency.advise`)
sizes PAL/maxcore through DELFIN's own logic and, when the caller
hands in EF1 `account()` records, adds measured CPU-efficiency advice.
This page records the numbers the advice threshold rests on, each with
its source. A number without a recorded run or a labeled secondary
source does not belong in advice.

## Measured populations (SOLVATOR / xTB shell steps)

The primary artifacts live in job folders under `.gate/practice/`
(not committed). The exact summary lines are fixed as records in
`tests/test_rs1_resource_sizing.py`; the source is the operator's
transcription (wave 22, md5 of run 1 summary
`2f195b869915ee805aa9602ec31be629`).

Run 1 — `.gate/practice/sv1_tiny_solvator/delfin_7922181.out` (NSOLV 12):

```
CPUs: 2
Time Limit: 00:30:00
State: FAILED (exit code 138)
CPU Utilized: 00:01:31
CPU Efficiency: 1.53% of 01:39:08 core-walltime
Job Wall-clock time: 00:24:47
```

Run 2 — `.gate/practice/sv1_tiny_solvator_r2/delfin_7922985.out` (NSOLV 3):

```
CPUs: 1
Time Limit: 00:45:00
CPU Utilized: 00:02:33
CPU Efficiency: 35.75% of 00:07:08 core-walltime
Job Wall-clock time: 00:03:34
```

## What these numbers do and do not support

- **Core-walltime counts LOGICAL CPUs (hyperthreads included).**
  Run 1: 00:24:47 wall x 2 logical CPUs x 2 threads = 01:39:08
  core-walltime; run 2: 00:03:34 x 2 = 00:07:08. The percentages are
  of allocated *logical* CPU time.
- **Run 1 ended FAILED (exit code 138)** — not a clean finish. The
  1.53 % comes from the scheduler's own `CPU Efficiency` summary line,
  NOT from sacct `TotalCPU` (which is a known artefact on killed
  jobs; EF1 reports UNKNOWN for TIMEOUT/CANCELLED rows and excludes
  them from advice).
- **The two runs are different tasks** (NSOLV 12 vs NSOLV 3). The
  defensible claim is only: SOLVATOR at pal 2 used 1.53 % of its
  allocated CPU time on a 12-molecule shell (run 1, FAILED); at pal 1
  a smaller 3-molecule shell used 35.75 %. No same-task wall-time
  comparison between the two runs is made — 37.7 s was ORCA's compute
  time inside run 2, while the job ran 00:03:34.

## Threshold

`LOW_PARALLEL_EFFICIENCY = 0.10` (`cluster_efficiency.py`) sits between
the two measured populations — 1.53 % (pal 2, FAILED run) and 35.75 %
(pal 1, smaller shell) — far from either, so a step kind is classified
low-parallel only when its measured mean is closer to the SOLVATOR
population than to the well-parallelised one. The measured DFT
populations (opt/freq/sp, ~90 % in wave-21 runs) sit far above.

## Advice semantics (what advise() does and never does)

- Advice is derived only from `account()` records the caller passes in;
  no records means today's sizing, with `"efficiency": None`.
- Only COMPLETED rows with a measured `cpu_efficiency` count;
  TIMEOUT/CANCELLED rows carry the TotalCPU artefact and are excluded
  with a count, never read as 0 %; NaN/inf records are excluded the
  same way.
- Matching is by exact job name — 'crest' is not 'crest_preopt'.
- **Advice only.** The user's PAL, maxcore and time limit are never
  changed by advice; a low-parallel step gets a note naming the
  measured mean, the run count and the source, and leaves the
  decision to the user.
- Convergence settings are never touched by advice; a deviation from a
  workflow default needs its own stated reason.

## What the advice would have saved (wave-20/21 practice runs)

Run 1 above requested 2 cores x 00:30:00 for a step that used 1.53 %
of its allocated CPU time and was killed at its limit (FAILED, exit
138): 2 cores x 00:24:47 wall = 0.83 core-hours allocated, of which
0.013 core-hours actually computed (CPU Utilized 00:01:31 — the
hyperthread factor sits in the scheduler's denominator, not in the
allocation). At pal 1 with a time limit sized from the measured wall
time, the same submit would have held 1 core instead of 2 and its
limit would have matched reality. EF1's accounting rows that carry
this number did not exist before this package's phase 1; future
measured runs feed the same path through `account()`. No number on
this page comes from anywhere but a recorded run or a labeled
transcription.

