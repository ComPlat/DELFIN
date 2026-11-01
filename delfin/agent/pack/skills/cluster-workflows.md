---
domains: code
name: cluster-workflows
description: Run DELFIN's workflows efficiently on a SLURM cluster — submit through DELFIN, size from DELFIN's own logic, monitor throttled, recover via recalc, account for core-hours.
---
# Cluster workflows
> Submit through DELFIN, size from its own logic, monitor throttled, recover via recalc, report core-hours.

## Read in this order
1. Size: `cluster_efficiency.advise(step, system_size)` — PAL/maxcore/partitions from DELFIN's own logic (`cluster_utils.py:181`).
2. State/wait reason: throttled listing (`backend_slurm.get_pending_start_times()`; one squeue per interval).
3. On failure: `.err`/`.out` tail, then `sacct` via `cluster_efficiency.account(job_ids)` — never a per-id loop.
4. Fix path (`agent/job_fix.py`): classify first, propose second.

## Workflow and submit
Pick DELFIN's workflow (CONTROL `method=`, `delfin/submit_templates/`); never invent one, deviations need a reason. Submit through DELFIN only — node-local I/O (`$TMPDIR`, `submit_delfin.sh:217-222`), partitions (`slurm_submit.py:135`), limits (`calc_policy.py:768`). `--export=ALL` exactly — never with explicit variables (env-retrieval failures on this site, `backend_slurm.py:445-453`).

## Sizing
Earliest start over ALL matching partitions as one `--partition` list (bwUniCluster 3.0: cpu, cpu_il — measured 12 days apart with `--test-only`); never change PAL, maxcore or time to start earlier. On a host outside a scheduler job the auto-PAL is capped at `calc_policy.LOGIN_LOCAL_CORES` (2) — an uncapped 192-core auto-PAL was a measured defect.

## Monitor — the 25 s rule
Job state through DELFIN's throttled listing (`list_jobs`); one query per interval, shared process-wide (25 s minimum between two squeue, `backend_slurm.py:97`); batch ids, never per id in a loop. Pending reason from the same listing: `ReqNodeNotAvail, Reserved for maintenance` means WAIT — the nodes are reserved, not the request wrong; never resize to start earlier. Jobs on this account run in a fair-share group: a queue wait is often your own group's usage — read the wait reason before touching anything.

## Recover via recalc
The old calculation is never overwritten, deleted, or resubmitted by hand. Walltime kill (TIMEOUT, no "ORCA finished by error termination") → time_limit bump via recalc; OOM kill without an ORCA error → maxcore ×1.5; an ORCA-internal error is DELFIN's in-run recovery, not the agent's. Chemistry/convergence keywords are off-limits (`job_fix._CHEMISTRY_GUARD`).

## Report core-hours
`cluster_efficiency.account(job_ids)`: requested vs used core/memory/wall, queue wait, core-hours used and wasted per job. Quote numbers with the job id — a number without an id is not evidence. sacct TotalCPU ~0 % on TIMEOUT/CANCELLED jobs is an ARTEFACT, never a measurement; the cross-check is the job .out's "CPU Efficiency" line.

## Safety
Never raw `sbatch`/`scancel`; never overwrite a prior calculation; never loosen convergence to fit a queue; never lower pal/maxcore/time just to start earlier; every change goes through the /control path with the user's confirmation.
