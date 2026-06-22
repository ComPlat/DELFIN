# Selection Guide

Choose the cheapest mode that still protects DELFIN.

## `quick`

Use when:

- the change is small or medium
- the task is local to one subsystem
- no scheduler, cluster, recovery, or public API risk is central
- the task can be closed in one cycle

Typical DELFIN examples:

- small bugfix in a reporting module
- contained CLI fix with existing tests
- helper extraction in a non-critical module
- docs or prompt pack maintenance

Route:

Session Manager -> Builder -> Test

## `reviewed`

Use when:

- architecture quality matters more than raw speed
- a refactor touches large orchestration files
- user-facing behavior may change
- there is meaningful regression risk but no cluster specialization is central

Typical DELFIN examples:

- refactor in `delfin/cli.py`
- changes in `delfin/pipeline.py`
- public API changes in `delfin/api.py`
- CONTROL parsing behavior changes in `delfin/config.py`

Route:

Session Manager -> Critic -> Builder -> Test

## `cluster`

Use when:

- the task affects runtime realism
- local and SLURM behavior may diverge
- scratch, submission, monitoring, PAL, tool detection, or recovery is involved

Typical DELFIN examples:

- `delfin/runtime_setup.py`
- `delfin/qm_runtime.py`
- `delfin/dashboard/backend_slurm.py`
- `delfin/dashboard/backend_local.py`
- submit templates
- restart, retry, or recovery behavior

Route:

Session Manager -> Runtime -> Critic -> Builder -> Test

## `full`

Use when:

- preparing a release or merge gate
- deciding whether a larger cycle is safe to close
- validating cross-cutting quality for user-facing delivery

Typical DELFIN examples:

- release hardening
- final review for large refactor batches
- pre-merge confidence cycle for runtime plus architecture changes

Route:

Chief -> Session Manager -> Runtime -> Critic -> Builder -> Test
