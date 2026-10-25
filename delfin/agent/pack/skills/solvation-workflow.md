---
name: solvation-workflow
domains: [solvation, workflows]
description: >
  Build an explicit-solvation workflow (Packmol box + xTB relaxation +
  ORCA-xTB SOLVATOR shell) through DELFIN's checked workflow path, and
  know when implicit solvation is the right model instead.
---

# Explicit solvation the physically honest way

## When explicit, when implicit

- **Implicit** (CPCM/SMD in ORCA, ALPB/GBSA in xTB) is right when the
  solvent is a continuum correction to a free energy — cheap, no box
  size decisions, valid to typical ±1–2 kcal/mol accuracy against
  experiment for neutral solutes.
- **Explicit** (a real solvent box: `packmol_solvate`, then
  `xtb_solvator`) is right when the specific solute–solvent hydrogen
  bonding or coordination matters: pKa at the solvent interface,
  excited-state quenching, solute–solvent charge transfer. It costs
  far more per step and its errors are NOT the same as the implicit
  ones.
- Never mix both on one species: a CPCM/SMD model on a step AFTER the
  explicit shell counts the solvation free energy twice. DELFIN's
  lint flags this (`explicit-and-implicit-solvation-double-count`).

## The workflow DELFIN ships

`explicit_solvation` (contributed, `delfin/workflows/contrib/
w21_solvation_workflow.py`): packmol box → xTB relax → SOLVATOR shell.

```python
from delfin.workflows.registry import get
wf = get("explicit_solvation")
result = wf.run(config={}, geometry="solute.xyz",
                solvent_xyz="water.xyz", solvent="water",
                n_solvent=50, box_side=20.0, nsolv=12,
                solvator_method="GFN2-xTB")
```

## Declared limits (enforced, not advisory)

- `xtb_solvator.method` must be GFN0/1/2-xTB — SOLVATOR is an xTB
  feature (ORCA manual, %SOLVATOR block). A DFT name is an ERROR in
  the lint and a FAILED step in the adapter.
- The box must hold the solute + 2×tolerance (Packmol 2009,
  DOI 10.1002/jcc.21224); the adapter refuses with the minimum
  `box_side` in the message.
- `mult >= 1`, charge integer, `n_solvent`/`nsolv` >= 1 — the
  adapters refuse these before any binary runs.
- Charge and multiplicity belong to the SOLUTE; the workflow passes
  them into every step. A box of 50 waters is closed-shell — do not
  "spread" an open-shell solute over the solvent.
- The xTB box relaxation is structure preparation only. The result
  reports `levels` so nobody quotes an xTB energy as the shell energy.

## What DELFIN does NOT decide for you

- Whether 12 shell molecules (`nsolv`) is enough — it is a
  compromise (ORCA's own examples use ~10–20); a converged shell is
  checked by comparing nsolv 8/12/20, not assumed.
- Whether the box density is physical — `box_side` sets it; a 20 Å
  box with 50 waters is ~0.03 g/cm³, far below liquid water. Increase
  `n_solvent` until the density is right for the question asked.
