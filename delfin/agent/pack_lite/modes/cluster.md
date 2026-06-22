# Mode: cluster

Use this whenever DELFIN's local/HPC execution realism is part of the task.

## Objective

Protect operational correctness across local systems and SLURM-like cluster use.

## Route

Session Manager -> Runtime -> Critic -> Builder -> Test

## Use this for

- runtime setup
- tool resolution
- submission behavior
- monitoring or stop/restart flows
- scratch and environment handling
- recovery and retry behavior
- local/SLURM backend differences

## Typical DELFIN triggers

- `delfin/runtime_setup.py`
- `delfin/qm_runtime.py`
- `delfin/dashboard/backend_local.py`
- `delfin/dashboard/backend_slurm.py`
- `delfin/orca_recovery.py`
- `delfin/submit_templates/*`

## Why this route exists

DELFIN can appear correct locally while still being wrong operationally on HPC.
Runtime review must happen before implementation. Architecture review (Critic)
also happens before the Builder, so rework cost stays low and both runtime and
architectural concerns are resolved before any code is written.

## Acceptance expectation

- local implications are stated
- cluster implications are stated
- diagnostics and failure modes remain understandable
- validation reflects the runtime risk
