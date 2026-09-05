---
name: pipeline-build
domains: code
---
# Build a computational-chemistry pipeline
> Assemble, validate, persist and run a pipeline (Bausteine) through the delfin-tools MCP server.

The user wants a multi-step calculation assembled from existing building
blocks — or a missing block added and wired in, then saved as a reusable
Application and run. Work **only** through the MCP server `delfin-tools`;
never hand-roll pipeline logic.

If `delfin-tools` is not registered, say so and stop. Improvising a
pipeline without it produces a spec nothing can run.

## The loop

Start by calling **`get_guide`**, then:

1. **DISCOVER** — `describe_capability` / `catalog` / `describe_key` (use
   only the allowed enum/key values), `compatible_successors` to order
   the steps.
2. **BUILD** — write a spec per `schemas.pipeline_spec`; set only the
   params that are actually needed (defaults + autowiring fill the rest).
   Check what will really run with `resolve_spec`.
3. **VALIDATE** — `validate_spec`; implement every diagnostic fix
   concretely until it is clean.
4. **MISSING BLOCK** — nothing fits? `new_capability_template` → fill in
   the code → `register_module` (it checks + integrates), then back to
   BUILD.
5. **PERSIST** — `save_application`; it then appears in the Pipelines tab.
6. **RUN** — `submit_application` (local or SLURM); results land in
   `~/calc`.
7. **DIAGNOSE** — `run_diagnostics(run_id)`: read status / error / logs
   and iterate.

## Rules

- Always `validate_spec` before a run.
- ORCA and Turbomole are licensed — **never install them**; only state
  when they are required.
- Report the saved application name and the run id; a pipeline the user
  cannot find again was not delivered.
