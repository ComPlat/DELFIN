# Integrating a DELFIN workflow — ONE contract

Every DELFIN workflow — classic CONTROL.txt, pipeline, user — is described by
**one** type, [`WorkflowDescriptor`](../delfin/tools/workflows.py), which
**extends** the existing `Application` (`delfin/tools/_application.py`); there
is **no parallel workflow type**. This is the documented recipe for
integrating a NEW workflow, or extending an existing one, so the agent can see
the whole repertoire, extend any workflow and integrate new ones the same way.

For assembling a workflow from building blocks (rather than integrating a new
one) see [`build_your_workflow.md`](build_your_workflow.md).

## The descriptor

A `WorkflowDescriptor` extends `Application` with exactly these fields:

| field             | meaning |
|-------------------|---------|
| `kind`            | `classic-control` \| `pipeline` \| `user` |
| `purpose`         | one sentence on what the workflow is for |
| `entry`           | how the workflow is started (the dispatch gate / CLI) |
| `control_keys`    | the CONTROL.txt keys it reads, with type, default, allowed values and what each changes |
| `steps`           | the building blocks (pipelines) or the ORCA jobs a classic workflow writes |
| `artifacts`       | files / directories the workflow produces |
| `extension_points`| where a step may be inserted (before / after / in place) and what it receives |
| `resume`          | resume / recovery behaviour |
| `backend`         | local vs Slurm |

It is constructed exactly like an `Application` (`register_application` works),
and serializes round-trip through `to_dict` / `from_dict`.  Register nothing
twice: `register_classic_workflows()` is idempotent, and classic descriptors
are listed together with pipeline and user workflows by `list_workflows()` and
`build_manifest()["workflows"]`.

## CONTROL keys — validated, never invented

Every CONTROL key a descriptor declares comes from
`delfin.common.control_validator.CONTROL_FIELD_SPECS`.  Use the helper so the
descriptor can never drift from `validate_control_config`:

```python
from delfin.tools.workflows import control_key_spec   # None if unknown
from delfin.common.control_validator import CONTROL_FIELD_SPECS

names = ["method", "functional", "main_basisset", "metal_basisset", "PAL"]
keys  = tuple(control_key_spec(n) for n in names)      # raises if a name is unknown
```

A descriptor that names a key `control_validator` does not know raises at
registration time; type, default and help are copied from the validator.

## Steps come from building blocks

The whole block list is `delfin.tools.list_steps()` (`orca_sp`, `orca_opt`,
`orca_freq`, `orca_tddft`, `smiles_to_xyz`, `occupier`, `esd`, `imag_fix`, …).
A step is the PipelineSpec shape:

```python
{"step": "orca_sp", "label": "sp", "kwargs": {"charge": 0, "mult": 1,
  "method": "B3LYP"}}
```

For a classic-control workflow the `steps` DESCRIBE the jobs it writes; the
descriptor **describes**, it does not reimplement — the classic engines in
`delfin/workflows/pipeline.py` are completely untouched.

## Extension — "add step X after step Y of workflow W"

`add_step` returns a NEW spec: W's steps unchanged, then X (or X before W, or
in place):

```python
from delfin.tools.workflows import workflow_for, add_step
occupier = workflow_for("OCCUPIER")
extended = add_step(occupier,
    {"step": "orca_sp", "label": "dlpno_sp",
     "kwargs": {"charge": 0, "mult": 1, "method": "DLPNO-CCSD(T)"}},
    after="spin_states")           # or before= / replace=
```

An unknown anchor raises `KeyError`.  **Always** validate the result through
the existing checks before anything runs:

```python
from delfin.tools import user_workflows as uw
report = uw.check_spec(extended)
assert report["ok"], report["errors"]
```

## A NEW workflow — the documented path

1. **Scaffold.** `delfin-agent new-workflow NAME [--kind pipeline|classic-control|user]`
   (backed by `delfin.tools.scaffold.build_workflow_scaffold`) writes a module
   registering a `WorkflowDescriptor` and an acceptance test.  The scaffold
   validates the step spec through `check_spec` **before** writing anything, so
   it cannot emit a workflow that would not run.

2. **Descriptor, CONTROL keys, steps.** Fill the TODO `purpose` and `entry`;
   list the CONTROL keys it reads (all from `CONTROL_FIELD_SPECS`); assemble
   `steps` from building blocks.

3. **Acceptance test with fake ORCA.** Drive a run with
   [`tests/scripts/fake_orca.sh`](../tests/scripts/fake_orca.sh) on `PATH`.  The
   fake is invoked exactly as the real runner would be — `orca <input.inp>`
   with cwd = the run directory — and copies the output files of a chosen
   fixture run (`$FAKE_ORCA_FIXTURE`) next to the input.  The test must exercise
   the real DELFIN launch path (`delfin.orca.run_orca` → the adapter) and get a
   genuine result from the fixture, never a hand-returned stub.

4. **Parity test** if the new workflow replaces an existing one: run the old
   and the new on the same input and record both numbers side by side
   (`same-level-of-theory` rules apply — compare like with like).

5. **Register.** `register_application(wf)` for a pipeline / user workflow; add
   it to `register_classic_workflows()` for a classic CONTROL workflow.
   Registration is descriptive only — it never changes how a workflow runs.

## Acceptance criteria

- every classic workflow has a descriptor whose CONTROL keys match `control_validator`;
- `list_workflows()` / `build_manifest()["workflows"]` list classic-control,
  pipeline and user workflows together;
- an extension `"W then X on the chosen state"` validates under `check_spec`;
- a scaffold writes a workflow that passes its own acceptance test;
- no classic workflow's behaviour changed — its existing tests stay green
  (`grep tests/ for occupier, esd, imag, redox, co2, classic`).

The agent reaches this knowledge through this doc and the
`workflow-integration` skill on demand; it is never added to the always-on
prompt (the built prompt for a non-chemistry task is identical on main and any
branch).
