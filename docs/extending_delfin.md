# Extending DELFIN with new steps

DELFIN's pipelines are built from **steps**.  Each step is an *adapter*:
a small Python class that wraps a tool (an external program, a parser, an
analysis) and exposes it to the workflow engine through a declared
**contract**.  This page shows the three ways to add one, and the checks
that keep it honest.

## 1. Scaffold a new adapter (`delfin-agent new-adapter`)

The fastest supported path:

```bash
delfin-agent new-adapter my_tool --kind program   # or parser|analysis
```

This writes three files and prints their paths:

| File | Purpose |
|---|---|
| `delfin/tools/adapters/<name>.py` | the adapter skeleton |
| `tests/test_adapter_<name>.py` | a contract test the suite will run |
| `docs/extending_delfin.md` | a docs stub to fill in |

The skeleton is a complete `StepAdapter` — declared `category`, `params`,
`validate_params`, and an `execute` that returns a `StepResult` — ending
in `register(<Class>())`.  The **registry entry needs no hand edit**: the
name is in the registry the moment the module is imported, because
`delfin.tools._registry._discover_builtin_adapters` walks the adapters
directory on every discovery (delfin/tools/_registry.py).  Fill in the
`TODO` in `execute()` and you have a step.

The module name must be a plain Python identifier — the same rule
`delfin.tools.platform.module_name_error` enforces for agent-built
modules, because the file is imported immediately.

## 2. Add a step to the repository (by hand)

Scaffold with the command above and edit the skeleton, or write the
module directly under `delfin/tools/adapters/`.  The contract elements
(optional but expected — `tests/test_tool_contracts.py` requires every
adapter to declare a category):

```python
class MyToolAdapter(StepAdapter):
    name = "my_tool"                     # unique step name
    description = "One-line description"
    produces_geometry = False            # True if it writes an .xyz
    category = "program"                 # program | parser | analysis | ...
    params = (ParamSpec("input", "path", required=True),)
    produces = ("geometry",)             # capability tags emitted
    consumes = ()                        # capability tags required upstream
```

Registration is one line at the bottom of the module:
`register(MyToolAdapter())`.  Nothing else to edit — no list, no
manifest.

## 3. Ship an adapter from outside the package (user directory)

Drop a `.py` file into `~/.delfin/adapters/` (or the directory named by
`$DELFIN_ADAPTERS_DIR`).  `_load_user_adapters` imports it during
discovery; a broken file is skipped so it never blocks the rest.  This is
also how the agent integrates freshly built blocks at runtime
(`delfin.tools.platform.register_module`).

## Workflows and skills

- **Workflow**: a pipeline YAML/CONTROL that chains steps.  See the
  templates under `delfin/tools/templates/` for the vocabulary; the
  CONTROL keys are defined in `delfin/config.py`.
- **Skill**: a playbook under `.delfin/skills/<name>/SKILL.md`.  A skill
  may *describe* a workflow; it must never widen a permission, name a
  credential path, or tell the agent to skip a check.

## The checks that apply

1. **Contract tests** (`tests/test_tool_contracts.py`): every adapter
   builds a `StepContract`; required params are actually enforced;
   declared params accept realistic values.  Scaffolded adapters ship a
   test in this style.
2. **Science lint** (`mcp delfin-tools scientific_lint`): physics
   plausibility of a spec — composite-method/basis redundancy, solvent
   model with no solvent, impossible multiplicity.  Run it on any spec
   that composes your new step into a workflow.
3. **Self-check (PC1, `delfin-agent prcheck`)**: before opening a PR, the
   pre-PR self-check runs the guard/invariant tests that match the areas
   your change touches.

## Discovery guarantee

`tests/test_wb3_discovery.py` fails when any module under
`delfin/tools/adapters/` registers adapters that discovery never reaches,
and `tests/test_wb3_scaffold.py` fails when a newly added registering
module is invisible to `list_steps()`.  These are the tests that would
have caught the wave-20 OR1 defect (three adapters registered on import,
nothing imported their modules) — do not weaken them.
