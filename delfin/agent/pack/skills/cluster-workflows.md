---
domains: code
name: cluster-workflows
description: Run DELFIN's workflows efficiently on a SLURM cluster — submit through DELFIN, size from DELFIN's own logic, monitor throttled, recover via recalc, account for core-hours.
---
# Cluster workflows
> Submit through DELFIN, size from its own logic, monitor throttled, recover via recalc, report core-hours.

## 0. Read in this order
1. `cluster_efficiency.advise(step, system_size)` — PAL/maxcore/partitions
   from DELFIN's own logic (`delfin/agent/cluster_efficiency.py:416`).
2. The job's state and wait reason through the throttled listing
   (`backend_slurm.get_pending_start_times()`, one squeue per interval).
3. On failure: the job's `.err`/`.out` tail, then `sacct` through
   `cluster_efficiency.account(job_ids)` — never a per-id loop.
4. The fix path (`agent/job_fix.py`): classify first, propose second.

## 1. Pick the workflow, never invent one
DELFIN's workflows (CONTROL `method=`, `delfin/submit_templates/`)
encode what to compute; follow them, deviations need a reason.
Small systems first.

## 2. Submit through DELFIN — never raw sbatch
Use DELFIN's submit path: node-local I/O (`$TMPDIR`,
`submit_delfin.sh:217-222`), partitions (`slurm_submit.py:135`),
limits (`calc_policy.py:768`). Never call `sbatch` yourself.
`--export=ALL` exactly — never with explicit variables (that form
caused env-retrieval failures on this site; backend_slurm.py:445-453).

## 3. Size from DELFIN's own logic
`cluster_efficiency.advise(step, system_size)` gives PAL/maxcore from
DELFIN's auto-config (`cluster_utils.py:181`) and configured partitions.
Rule: earliest start over ALL matching partitions as one `--partition`
list (bwUniCluster 3.0: cpu, cpu_il — measured 12 days apart with
`--test-only`); never change PAL, maxcore or time to start earlier.
On a host outside a scheduler job the auto-PAL is capped at
`calc_policy.LOGIN_LOCAL_CORES` (2) — an uncapped 192-core auto-PAL
was a measured defect (SL1).

## 4. Monitor throttled — the 25 s rule
Job state through DELFIN's throttled listing (`list_jobs`); one query
per interval, shared process-wide (25 s minimum between two squeue,
`backend_slurm.py:97`); batch ids, never per id in a loop. A pending
job's reason comes from the same listing (`get_pending_start_times`)
— `ReqNodeNotAvail, Reserved for maintenance` means WAIT, the nodes
are reserved, not the request wrong; never resize to start earlier.

## 5. Recover through the recalc path
Recover via the recalc path (`agent/job_fix.py`): the old calculation
is never overwritten, deleted, or resubmitted by hand. A walltime kill
(TIMEOUT, no "ORCA finished by error termination") proposes a
time_limit bump through recalc; an OOM kill without an ORCA error
proposes maxcore x1.5; an ORCA-internal error is DELFIN's in-run
recovery, not the agent's. Chemistry/convergence keywords are
off-limits (`job_fix._CHEMISTRY_GUARD`).

## 6. Report the core-hours
`cluster_efficiency.account(job_ids)` reports requested vs used
core/memory/wall, queue wait, core-hours used and wasted per job. Quote
numbers with the job id — a number without an id is not evidence.
sacct TotalCPU is ~0 % on TIMEOUT/CANCELLED jobs on this cluster — an
ARTEFACT, never a measurement; the honest cross-check is the job .out's
"CPU Efficiency: x % of ... core-walltime" line. Account records only
count COMPLETED runs' efficiency.

## 7. Fair share
Jobs on this account run in a fair-share group: a queue wait is often
your own group's usage, not a broken scheduler. Read the wait reason
(section 4) before touching anything.

## Safety
Never raw `sbatch`/`scancel`; never overwrite a prior calculation; never
loosen convergence to fit a queue; never lower pal/maxcore/time just to
start earlier; every change goes through the /control path with the
user's confirmation.
