---
name: workflow-integration
domains: code
description: Integrate or extend a DELFIN workflow the ONE-contract way — descriptor (kind/purpose/entry/CONTROL keys/steps/artifacts/extension points/resume/backend) extending Application, CONTROL keys validated against control_validator, steps from building blocks, acceptance test with fake_orca, registration, scaffold. Loaded for workflow tasks.
---

# Integrate a DELFIN workflow (ONE contract)

Every DELFIN workflow — classic CONTROL.txt, pipeline, user — is described by
**one** type: `WorkflowDescriptor`, which **extends** `delfin.tools._application.Application`
(no parallel type). Before integrating a workflow, use the existing `build-a-workflow`
skill to see the pieces; this skill is about wiring a NEW workflow in, or extending an
existing one, so the agent can see the whole repertoire, extend any workflow and
integrate new ones the same way.

## The contract

A descriptor carries exactly these fields (all filled for a real workflow):

| field             | meaning |
|-------------------|---------|
| `kind`            | `classic-control` \| `pipeline` \| `user` |
| `purpose`         | one sentence on what the workflow is for |
| `entry`           | how it is started (the dispatch gate / CLI) |
| `control_keys`    | the CONTROL.txt keys it reads (name,type,default,allowed,changes) |
| `steps`           | building blocks (pipelines) or the ORCA jobs a classic workflow writes |
| `artifacts`       | files/dirs the workflow produces |
| `extension_points`| where a step may be inserted (before/after/in place), and what it receives |
| `resume`          | resume/recovery behaviour |
| `backend`         | local vs Slurm |

## CONTROL keys must be validated, never invented

Derive every key from `delfin.common.control_validator.CONTROL_FIELD_SPECS`:

```python
from delfin.tools.workflows import control_key_spec
spec = control_key_spec("functional")   # None if control_validator does not know it
```

A descriptor that names a key control_validator does not know is refused
(raises). Type/default/help are copied from the validator, so the descriptor
can never drift from `validate_control_config`.

## Steps come from building blocks

`delfin.tools.list_steps()` is the whole block list (`orca_sp`, `orca_opt`,
`orca_freq`, `orca_tddft`, `smiles_to_xyz`, `occupier`, `esd`, `imag_fix`, …).
A step spec is `{"step": "<name>", "label": "<label>", "kwargs": {...}}`.
For a classic-control workflow the steps DESCRIBE the jobs it writes — the
descriptor describes, it does not reimplement; the classic engines are untouched.

## Extension

`add_step(descriptor, {"step": "orca_sp", "label": "x", "kwargs": {...}},
after="<step-label>")` returns a NEW spec: the workflow's steps unchanged,
then X. `before=` / `replace=` also work. An unknown anchor raises `KeyError`.
Always validate the result with the existing checks:

```python
from delfin.tools import user_workflows as uw
report = uw.check_spec(extended)
assert report["ok"], report["errors"]   # run nothing until this is green
```

## A NEW workflow: the documented path

1. **Scaffold.** `delfin-agent new-workflow NAME [--kind pipeline|classic-control|user]`
   (or `delfin.tools.scaffold.build_workflow_scaffold`) writes a module
   registering a `WorkflowDescriptor` plus an acceptance test. The scaffold
   validates the step spec through `check_spec` BEFORE writing, so it cannot
   emit a workflow that would not run.
2. **Descriptor, CONTROL keys, steps.** Fill the TODO `purpose`/`entry`; list the
   CONTROL keys it reads (from `CONTROL_FIELD_SPECS`); assemble `steps` from
   building blocks.
3. **Acceptance test with fake ORCA.** Point a run at `tests/scripts/fake_orca.sh`
   (invoked as `orca <input.inp>`); the test must get a real result, not a stub.
4. **Parity test** if the new workflow replaces an existing one: run the old and
   the new on the same input and record both numbers side by side.
5. **Register.** `register_application(wf)` for a pipeline/user workflow, or
   add it to `register_classic_workflows()` for a classic CONTROL workflow.
   Registration is descriptive only — it never changes how a classic workflow
   runs.

## Check the repertoire

`list_workflows()` and `build_manifest()["workflows"]` list classic-control,
pipeline and user workflows **together** — the same contract for the whole
repertoire. The built prompt for a non-chemistry task must be identical on
main and your branch (never add this knowledge to the always-on prompt).

## Acceptance criteria

- every classic workflow has a descriptor whose CONTROL keys match control_validator
- the catalog lists all three kinds together
- an extension `"W then X on the chosen state"` validates under `check_spec`
- a scaffold writes a workflow that passes its own acceptance test
- no classic workflow's behaviour changed (its existing tests stay green)
