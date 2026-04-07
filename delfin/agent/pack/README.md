# DELFIN_AGENT

Repo-tailored agent orchestration pack for DELFIN.

This pack is optimized for DELFIN as it exists today:

- large orchestration modules such as `delfin/cli.py`, `delfin/pipeline.py`,
  `delfin/parallel_classic_manually.py`, and `delfin/runtime_setup.py`
- runtime and tool resolution in `delfin/qm_runtime.py`
- validation in `delfin/common/control_validator.py`
- recovery in `delfin/orca_recovery.py`
- local and SLURM backends in `delfin/dashboard/backend_local.py` and
  `delfin/dashboard/backend_slurm.py`
- scheduler and resource coordination in `delfin/global_manager.py`,
  `delfin/global_scheduler.py`, and `delfin/dynamic_pool.py`

## Structure

- `agents/`: role-specific system prompts
- `shared/`: common rules, DELFIN context, input template, verdict schema
- `routing/`: task routing and handoff rules
- `examples/`: minimal kickoff template
- `manifest.yaml`: machine-readable overview of the pack

## Design principles

- bounded work cycles, not open-ended swarms
- only one production code writer unless explicitly overruled
- test is a hard gate unless explicitly waived by the requester
- DELFIN-first architecture decisions over generic agent behavior
- special care around scheduler, runtime, recovery, and cluster behavior

## Intended use

1. Load `shared/delfin_context.md`.
2. Load `shared/work_cycle_rules.md`.
3. Load the relevant role prompt from `agents/`.
4. Use `shared/universal_input_template.md` for every active agent.
5. Route work with `routing/minimal_workflow_routing.md`.
6. Require final verdicts in `shared/minimal_final_verdict.md`.

## Notes

- The prompts are intentionally short and strict.
- They are written for DELFIN, not for generic software projects.
- They assume DELFIN should evolve toward a robust workflow platform across
  local systems and HPC clusters.

## Operational recommendation

For daily work, prefer the lightweight preset pack in
`../pack_lite/`.

Use the full `DELFIN_AGENT` pack when:

- the cycle is strategically important
- the route spans several high-risk subsystems
- release-level judgment is needed
- you need the full governance layer rather than the cheaper default presets
