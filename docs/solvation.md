# Solvation in DELFIN

DELFIN offers implicit continuum models inside the standard ORCA/xTB
adapters and two explicit-solvation tools. This page states what each
is for, and which limits DELFIN enforces.

## Implicit models

| Engine | Model | Adapter | Key |
|---|---|---|---|
| ORCA | CPCM, SMD | `orca` | `solvent`, `solvent_model` |
| xTB | ALPB, GBSA | `xtb_native`, `crest` | `solvent`, `solvent_model` |

Implicit solvation is a continuum correction: cheap and adequate for
free energies of neutral solutes at typical (±1–2 kcal/mol) accuracy.
Choose it when you need a solvation *correction*, not a solvation
*structure*.

## Explicit solvation

Two tools, one contributed workflow:

- `packmol_solvate` — Packmol (Martínez/Andrade/Birgin/Martínez,
  J. Comput. Chem. 2009, 30, 2157–2164, DOI 10.1002/jcc.21224) packs a
  cubic box of solvent molecules around a solute. Params:
  `solvent_xyz` (required), `n_solvent` (default 50), `box_side`
  (default 20 Å), `tolerance` (default 2 Å).
- `xtb_solvator` — ORCA-xTB SOLVATOR adds an explicit shell of `nsolv`
  solvent molecules around a solute (ORCA manual, %SOLVATOR block).
  Params: `solvent` (required), `method` (default GFN2-xTB), `charge`,
  `mult`, `n_solvent` (default 12).
- `explicit_solvation` workflow (`delfin/workflows/contrib/
  w21_solvation_workflow.py`) — packmol box → xTB relaxation →
  SOLVATOR shell, with the level of each step reported.

### Enforced limits

Both adapters enforce their physics limits in `execute()` — the
product call path (`delfin/agent/campaign.py`) calls `execute()`
directly, so limits that lived only in `validate_params()` would never
be checked:

- `xtb_solvator`: `method` must be GFN0/1/2-xTB (SOLVATOR is an xTB
  feature; a DFT name would produce an invalid ORCA input only after
  queue time); `mult >= 1`; charge integer; `n_solvent >= 1`;
  non-empty solvent name.
- `packmol_solvate`: `n_solvent >= 1`, `box_side > 0`, `tolerance > 0`
  (the tolerance is a hard minimum inter-molecular distance); a box
  that cannot hold the solute extent plus 2×tolerance is refused up
  front with the minimum `box_side` in the message — Packmol would
  otherwise churn to its 300 s timeout.

### Lint rules (scientific lint)

- `xtb-solvator-method-not-xtb` (error): a non-xTB `method` on an
  `xtb_solvator` step.
- `explicit-and-implicit-solvation-double-count` (warning): a
  CPCM/SMD/… model on a step after the explicit shell — the solvation
  free energy would be counted twice. Justify the combination or drop
  the model.

## Which to use

| Question | Model |
|---|---|
| Free energy of a neutral molecule in water | implicit (SMD) |
| Solvent H-bonding to a specific site | explicit |
| Excited-state solvatochromism of a chromophore | explicit (ESD can combine both) |
| Quick screening of many candidates | implicit |

## Physical sanity of the box

`box_side` and `n_solvent` set the box density; check it against the
real solvent density (water: 1.0 g/cm³, 33.4 molecules/nm³). A 20 Å
box with 50 waters is ~0.03 g/cm³ — far below liquid water and wrong
for any question about solvation structure. The shell size (`nsolv`)
is a compromise; converge it by comparing 8/12/20, not by assuming.
