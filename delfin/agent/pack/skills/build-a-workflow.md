---
domains: code
name: build-a-workflow
description: Build the user a new DELFIN workflow — draft from DELFIN's catalog or templates, run the five checks, save as a versioned user workflow, export as a contribution. Never submit un-checked.
---
# Build a workflow
> Start from DELFIN's own pieces; check before anything runs; the original stays untouched.

## 1. Use DELFIN's workflows first
A new workflow starts from a named template (`delfin/tools/templates/`,
e.g. `classic_opt_freq`, or the OPI chain `opi_opt_freq_dlpno`:
OPI-built `orca_opt_opi` → `orca_freq_opi` → DLPNO-CCSD(T) single point
with knowledge-store-validated inputs) or a saved user workflow — never a step list
invented in chat. Extending means `user_workflows.extend_template`:
add, remove, replace or override, and the result is a NEW spec; the
original template or saved workflow is never modified.

## 2. Draft with the builder
`workflow_builder.build_workflow(request)` produces a draft spec from
the catalog, each step with a one-line reason, and runs the mandatory
checks in this fixed order before anything is saved or submitted:
`resolve_spec` -> `validate_spec` -> `scientific_lint` -> dry run ->
cost estimate (`cluster_efficiency.advise`, DELFIN's own sizing).
A failed check returns the finding verbatim — fix the spec and rerun;
never edit around a finding or drop the failing step silently. The
"dry run" stage is a re-validation with the geometry flag, NOT an
execution — nothing ever runs during building; submission is a
separate explicit act.

## 3. Nothing runs before every check passes
`checks["ok"]` is the gate. A lint error (e.g. multiplicity < 1, a
composite method with a redundant basis) blocks saving and submitting.
Warnings are shown to the user with the named fix. `submitted` in the
builder's result is never set by the builder — submission is a separate
explicit act through DELFIN's own submit path.

## 4. Save, reload, share
`user_workflows.save` versions the checked spec in the user's DELFIN
data directory; `load` returns the identical spec. Re-run means
re-validating the saved spec and submitting through DELFIN — a failed
step recovers through DELFIN's recalc path, the old calculation is
never overwritten. `export_contrib` writes a contribution file the user
can PR into `delfin/workflows/contrib/`.

## Safety
Never submit a spec whose checks did not pass; never silently change a
spec to make a check green; never overwrite a prior calculation; never
invent a step name — vocabulary comes from the catalog
(`list_capabilities`, `compatible_successors`).
