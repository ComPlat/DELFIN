---
domains: code
name: cluster-workflows
description: Run DELFIN's workflows efficiently on a SLURM cluster — submit through DELFIN, size from DELFIN's own logic, monitor throttled, recover via recalc, account for core-hours.
---
# Cluster workflows
> Submit through DELFIN, size from its own logic, monitor throttled, recover via recalc, report core-hours.

## 1. Pick the workflow, never invent one
DELFIN's workflows (CONTROL `method=`, `delfin/submit_templates/`)
encode what to compute; follow them, deviations need a reason.
Small systems first.

## 2. Submit through DELFIN — never raw sbatch
Use DELFIN's submit path: node-local I/O (`$TMPDIR`,
`submit_delfin.sh:217-222`), partitions (`slurm_submit.py:135`),
limits (`calc_policy.py:768`). Never call `sbatch` yourself.

## 3. Size from DELFIN's own logic
`cluster_efficiency.advise(step, system_size)` gives PAL/maxcore from
DELFIN's auto-config (`cluster_utils.py:181`) and configured partitions.
Rule: earliest start over ALL matching partitions as one `--partition`
list; never change PAL, maxcore or time to start earlier.

## 4. Monitor throttled
Job state through DELFIN's throttled listing (`list_jobs`); one query
per interval, shared process-wide; batch ids, never per id in a loop.

## 5. Recover through the recalc path
Recover via the recalc path (`agent/job_fix.py`): the old calculation
is never overwritten, deleted, or resubmitted by hand.

## 6. Report the core-hours
`cluster_efficiency.account(job_ids)` reports requested vs used
core/memory/wall, queue wait, core-hours used and wasted per job. Quote
numbers with the job id — a number without an id is not evidence.

## Safety
Never raw `sbatch`/`scancel`; never overwrite a prior calculation; never
loosen convergence to fit a queue.
