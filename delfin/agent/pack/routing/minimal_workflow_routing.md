# Minimal Workflow Routing

Use the smallest route that matches the task.
These routes align with the DELFIN_AGENT_LITE modes.

## quick (default)

Session Manager -> Builder -> Test

Use for: small bugfixes, docs, contained changes, agent-pack maintenance.

## reviewed

Session Manager -> Critic -> Builder -> Test

Use for: risky refactors, public API changes, config semantics, orchestration
file changes.

## cluster

Session Manager -> Runtime -> Critic -> Builder -> Test

Use for: runtime setup, tool resolution, submission, SLURM backends, recovery,
local/cluster divergence.

## full

Chief -> Session Manager -> Runtime -> Critic -> Builder -> Test

Use for: release hardening, milestone validation, cross-cutting changes.

## DELFIN-specific routing triggers

- Involve `Runtime` when touching:
  `runtime_setup.py`, `qm_runtime.py`, submit templates, dashboard backends,
  cluster detection, scratch handling, monitoring, restart, or recovery behavior
- Involve `Critic` when touching:
  `cli.py`, `pipeline.py`, `parallel_classic_manually.py`, `config.py`, `api.py`,
  or any orchestration-level file
- Involve `Chief` only for release gating or strategic decisions
