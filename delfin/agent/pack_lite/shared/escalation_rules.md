# Escalation Rules

Start in `quick` unless a real risk pushes the cycle higher.

## Escalate from `quick` to `reviewed` if any of these are true

- touching `delfin/cli.py`
- touching `delfin/pipeline.py`
- touching `delfin/parallel_classic_manually.py`
- changing `delfin/api.py` semantics
- changing CONTROL parsing or validation behavior
- changing result semantics visible to users

## Escalate from `quick` or `reviewed` to `cluster` if any of these are true

- touching `delfin/runtime_setup.py`
- touching `delfin/qm_runtime.py`
- touching local or SLURM dashboard backends
- touching submit scripts or cluster templates
- changing scratch, restart, monitoring, tool discovery, or recovery logic
- changing PAL, maxcore, scheduler/runtime coordination, or environment assumptions

## Escalate to `full` if any of these are true

- closing a larger milestone
- preparing a release candidate
- combining architecture and runtime changes
- change impact spans multiple user-facing entry points
- confidence must be strong enough for broad external use

## De-escalate when possible

- If runtime concerns are absent, do not include `Runtime`.
- If architecture concerns are minor, do not include `Critic`.
- If no strategic decision is needed, do not include `Chief`.
- If no external pattern choice matters, do not include `Research`.

## Cost rule

The best orchestration for DELFIN is not the largest one.
The best orchestration is the smallest route that still protects correctness,
reproducibility, runtime realism, and maintainability.
