# DELFIN_AGENT_LITE

Lightweight orchestration presets for DELFIN.

This pack is the recommended day-to-day operating mode when the full
`DELFIN_AGENT` pack would be too expensive or too heavy for the task.

It is not a separate philosophy.
It is a cheaper execution layer that reuses the same DELFIN priorities:

- platform quality over feature count
- bounded work cycles
- only one production code writer by default
- test as a hard gate unless explicitly waived
- extra agents only when risk justifies them

## Relationship to `DELFIN_AGENT`

Use `DELFIN_AGENT_LITE` as the default operating layer.
Escalate to the full `DELFIN_AGENT` pack only when:

- the task spans multiple cycles
- the architecture direction is disputed
- the work touches multiple high-risk subsystems at once
- a release-quality or strategy-quality decision is needed

This lite pack reuses the role prompts from `../pack/agents/`
and the shared context from `../pack/shared/`.

## Presets

- `quick`: normal implementation cycle
- `reviewed`: architecture-sensitive refactor or behavior change
- `cluster`: local/HPC, SLURM, tool resolution, recovery, submission
- `full`: release-hardening, merge gating, final confidence pass

## Recommended default

For most DELFIN tasks, use:

Session Manager -> Builder -> Test

Do not add Critic, Runtime, Research, or Chief unless there is a concrete reason.

## File map

- `manifest.yaml`: overview
- `modes/quick.md`
- `modes/reviewed.md`
- `modes/cluster.md`
- `modes/full.md`
- `shared/escalation_rules.md`
- `shared/selection_guide.md`
