# Build your own workflow

DELFIN's agent can build a workflow for you — from scratch out of
DELFIN's catalog, or by extending and modifying an existing template or
one of your saved workflows. Every draft is checked before anything
runs, and the original template or workflow is never modified.

## The pieces

| Module | Role |
|---|---|
| `delfin/agent/workflow_builder.py` | request text → draft spec, per-step reasons, the mandatory checks, cost estimate |
| `delfin/tools/user_workflows.py` | save / list / load / re-validate / export your workflows |
| `delfin/tools/templates/` | the built-in parametric templates drafts start from |
| `delfin/workflows/contrib/` | contributed workflow wrappers; exported files PR here |
| `delfin/tools/mcp_server.py` | `list_capabilities`, `describe_capability`, `compatible_successors`, `resolve_spec`, `validate_spec`, `scientific_lint` over MCP |

## The five checks, in order

Before anything is saved or submitted, `workflow_builder.run_checks`
runs — in this fixed order:

1. **resolve** — `resolve_spec` fills each step's contract defaults and
   shows what will be auto-wired from upstream steps.
2. **validate** — `validate_spec` statically checks wiring and required
   parameters. Nothing executes.
3. **lint** — `scientific_lint` applies physics rules (a composite
   method with a redundant basis, relativity without a relativistic
   basis, a solvation model with no solvent, an impossible
   multiplicity). A lint *error* blocks; warnings are shown with the fix.
4. **dry run** — static validation with the geometry flag. This is a
   re-validation (`Pipeline.dry_run` is an alias of `Pipeline.validate`),
   not an execution — nothing runs at this stage either.
5. **cost** — `cluster_efficiency.advise` per step, DELFIN's own
   resource sizing (PAL/maxcore/partitions).

A failed check returns the finding verbatim. The fix is named; the
builder never applies it silently and never submits.

## Extending or modifying

```python
from delfin.tools import user_workflows as uw

# extend the classic template with a single point: NEW spec
spec = uw.extend_template(
    "classic_opt_freq",
    add=[{"step": "orca_sp", "kwargs": {"charge": 0}}],
    new_name="opt_freq_sp",
)
```

`add` steps go last; `remove`/`replace` match steps by label;
`overrides` extend the spec's defaults. The source template is read
only — compare `delfin.tools.templates.classic_opt_freq.to_dict()`
before and after.

## Saving and sharing

```python
uw.check_spec(spec, geometry=True)      # the check chain, returns ok
path = uw.save(spec, name="opt_freq_sp")  # ~/.delfin/workflows/<name>/vNNN.json
uw.load("opt_freq_sp")                  # identical spec back
uw.list_workflows()                     # everything you saved
uw.export_contrib("opt_freq_sp", out_dir="contrib/")  # PR-able file
```

Saving the same name again creates the next version; loading without a
version returns the newest. Re-running a saved workflow re-validates it
and submits through DELFIN's own submit path; a failed step recovers
through the recalc path — the old calculation is never overwritten.

## The skill

The agent-side playbook is the packaged skill `build-a-workflow`
(`delfin/agent/pack/skills/build-a-workflow.md`): workflows first, the
five checks as the gate, save/reload/export, and the safety rules.
