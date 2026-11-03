---
name: mlp-method
domains: code
description: >
  Choose and drive DELFIN's machine-learned potential adapters (ANI-2x,
  AIMNet2, MACE-OFF, CHGNet, ...) safely: which elements and charges each
  model covers, when an MLP is a pre-screen only, and how the adapters
  enforce these limits fail-closed. Read this before any mlp_* or ase_*
  step.
---

# Machine-learned potentials in DELFIN

DELFIN's MLP steps are pre-screening and pre-optimization tools. They are
**fast approximations** — never publication-grade energies — and the adapters
enforce each model's limits fail-closed: a request outside a model's declared
scope is REFUSED before any model runs, never computed silently wrong.

## Steps and backends

- `mlp_single_point` / `mlp_optimize` / `mlp_rank`
  (`delfin/tools/adapters/mlp.py`): backends `ani2x`, `aimnet2`, `mace_off`,
  `chgnet`, `m3gnet`, `schnetpack`, `nequip`, `alignn` (aliases accepted,
  case-insensitive: `mace` = `mace_off`, `ani-2x` = `ani2x`).
- `ase_optimize` / `ase_single_point` / `ase_md`
  (`delfin/tools/adapters/ase_adapter.py`): calculators `xtb`, `orca`,
  `emt`, `eam`, `lj`, `mace`, `ani` — same realism tables, same enforcement.
- `mlp_screen` workflow (`delfin/workflows/contrib/w21_mlp_screen.py`):
  ranks a conformer ensemble through `mlp_rank`, re-using the adapter
  tables (single source of truth).

## Agent tools (MCP, `mcp__delfin-ops__*`)

Registered like every other DELFIN tool; the agent calls them by name:

- `mlp_list_backends` — availability, element scope, version and install
  hint per backend. Read-only: passes the agent's gate untouched.
- `mlp_single_point`, `mlp_optimize`, `mlp_rank`, `mlp_preoptimize` — the
  run tools. They compute (a model runs, minutes on CPU), so they are NOT
  on the read-only allowlist: in default mode the side-effect gate asks
  before each call, like the calculation runners. Each takes an xyz
  string and answers JSON; an unknown backend or a missing package comes
  back as `ok=false` with `available`/`reason`/`install_hint` — a
  diagnosis, never a traceback. `mlp_preoptimize` returns the relaxed
  geometry as `xyz` plus `note`: the geometry is a starting guess for the
  DFT step, never a final structure.

## What each model covers (declared, with sources, enforced)

| Backend | Elements | Charge/mult | Valid for | NOT valid for |
|---|---|---|---|---|
| ani2x | H,C,N,O,S,F,Cl | no (refused) | neutral closed-shell molecular clusters | periodic, charged, open-shell, metals |
| aimnet2 | H,B,C,N,O,F,Si,P,S,Cl,As,Se,Br,I | **yes** | neutral AND charged organic molecules | extended periodic solids |
| mace_off | H,C,N,O,F,P,S,Cl,Br,I | no | organic/bio-organic/drug-like | periodic, metals |
| chgnet/m3gnet/alignn | model-file scope | no | materials/MP-trained scope | decided by the model file |
| schnetpack/nequip | user model scope | no | the trained model's training set | out-of-scope elements |
| ase xtb | all spd-block ≤ Z=86 | yes | molecular pre-screens | fully periodic bulk |
| ase emt | Al,Cu,Ag,Au,Ni,Pd,Pt + H,C,N,O | no | educational/test | serious production |
| ase lj/mace/eam | model-defined | no | generic/test potentials | element-specific physics |

Sources: ANI-2x Devereux et al., JCTC 16 (2020) 4192, DOI
10.1021/acs.jctc.0c00121; AIMNet2 Anstine, Zubatyuk, Isayev, Chem. Sci. 16
(2025) 10228, DOI 10.1039/D4SC08572H; MACE-OFF Kovacs et al., JACS 147
(2025) 9254, DOI 10.1021/jacs.4c07099; xTB GFN Grimme, Bannwarth, Shushkov,
JCTC 13 (2017) 1989, DOI 10.1021/acs.jctc.7b00118; emt/eam/lj from the
installed ASE calculator docstrings.

The tables live in `SUPPORTED_ELEMENTS` / `ELEMENT_TABLE_SOURCES` /
`VALIDITY` (mlp.py) and `ASE_SUPPORTED_ELEMENTS` / `ASE_SOURCES` /
`ASE_VALIDITY` (ase_adapter.py). A refusal names the element and cites the
source — that is the adapter telling you the physics, not a bug.

## How the enforcement behaves (so you can explain a refusal)

- Unknown backend/calculator, empty or non-string name: `ValueError` in
  `validate_params`, before anything runs.
- Element outside the declared set: checked on the ACTUAL structure in
  `execute()`, before any model is created; the result is a FAILED step
  with the unsupported element(s) named. Models with model-file scope are
  NOT element-refused (a refusal there would be wrong physics).
- Charge ≠ 0 or mult ≠ 1 on a backend that cannot model them (all except
  `aimnet2`; in ASE only `xtb`/`orca`): refused — the calculator would
  silently ignore it and produce a wrong number.
- Install errors (`torch`, `aimnet2calc`, `mace`, ... missing): surface as
  a FAILED step with the import error, not a traceback.

## When an MLP is enough, and when DFT must follow

- **Good**: pre-optimizing geometries before DFT, ranking conformer
  ensembles cheaply, fast scans over many candidates. State "MLP
  pre-screen" in the answer.
- **Not good**: reaction barriers < ~1 kcal/mol, binding energies for
  publication, spin-state splittings on metals, anything charged outside
  `aimnet2`'s scope — DFT (ORCA) must follow or replace it. Metals,
  open-shell species and charged species are out of scope for
  ani2x/mace_off outright (the adapter refuses them); only `aimnet2`
  accepts charge and multiplicity, and its spin physics is the model's,
  not a variational treatment — confirm with DFT whenever the spin state
  is part of the question. A minimum is only verified with zero imaginary
  frequencies from a Hessian — an MLP optimization alone never verifies
  one.
- **Always confirm final energies with DFT.** Every MLP number is a
  screening value in eV on the model's own scale — not comparable with
  DFT energies, not with another backend's energies, and not a
  publication number. A relaxed geometry is a *starting guess*: the DFT
  energy of that geometry is the only final number. Label every MLP
  result "MLP pre-screen" and follow it with the ORCA/DFT step that
  decides.
- CHGNet/M3GNet/ALIGNN are materials-scoped (periodic, MP/JARVIS-trained):
  for molecular chemistry prefer ani2x/aimnet2/mace_off.

## The `mlp_preopt` CONTROL key

`mlp_preopt = none | ani2x | aimnet2 | mace_off | ...` in CONTROL.txt
relaxes the start geometry with the chosen MLP BEFORE the ORCA steps
(off by default: nothing runs, the inputs are byte-identical). Use it
when the builder geometry is poor — large flexible molecules, many
rotatable bonds — and the system is element-clean and closed-shell
(elsewhere the adapter refuses). A failed or non-converged MLP step never
stops the run: the original geometry is kept, the run continues, and the
warning says so. FieldSpec `:2306` in
`delfin/common/control_validator.py`; notes in
`delfin/common/control_key_notes.json` (key `mlp_preopt`).

## ExtOpt providers: delfin-mlp-extopt (default) and orca-mlips

ORCA `! ExtOpt` runs can be driven by two providers:

- **`delfin-mlp-extopt`** (default, `delfin/mlp_tools/orca_extopt.py`): our
  own wrapper around the registry backends; one fresh process per ORCA
  step (torch import + model load each step — fine for preopt, slow for
  long optimizations).
- **`orca-mlips`** (`delfin/mlp_tools/orca_mlips.py`, optional package
  `pip install orca-mlips[...]`, MIT): resident model server, so repeated
  calls during an optimization are fast; supports xTB implicit solvent
  (`--solvent`, `--solvent-model alpb|cpcmx`) and the two-step analytical
  Hessian (`--dump-hessian` + `%geom inhess Read / InHessName`) for
  OptTS/Freq. UMA weights are gated (user's own HF login); MACE-OMOL-0
  weights are non-commercial (Academic Software License). Use the
  prefixed ProgExt names (`orca-mlips-uma`), not the short ones.

Plan jobs with `mlp_extopt_plan` (dry run — writes input + job script,
submits nothing); GPU partition configured -> one sbatch job with
`--nodes=1` so ORCA and the MLIP server share a node; no GPU -> local CPU
fallback (<= 2 cores, < 30 min, small systems only) or refusal.

## When an MLIP is enough, and when DFT must follow

### MLIP preopt / TS guess: sound when

- neutral, closed-shell organic molecules inside the model's element set
  (see the tables above);
- the goal is a STARTING GEOMETRY or a conformer ranking, and a DFT
  optimization/verification follows;
- for TS guesses: the MLIP Hessian is a *guess* for the `%geom inhess`
  first guess — the TS must be re-optimized and characterized with DFT.

### DFT check REQUIRED when

- the species is charged or open-shell outside `aimnet2`'s scope (adapter
  refuses), or contains metals outside the model's training set — an MLP
  number here is not just imprecise, it is silently wrong physics;
- reaction barriers or spin-state splittings are the question;
- an MLIP optimization reports convergence: still verify with a DFT
  Hessian — MLIP local minima can keep extra imaginary frequencies
  (orca-mlips README explicitly warns about FP32 precision there);
- any number is going into a publication: MLP energies are screening
  values on the model's own scale, never comparable with DFT.

### Licence discipline (LICENSE_GUARD.md family 4)

No weights or HF tokens ever enter DELFIN, the repo or the logs; gated
checkpoints are the user's own login. `refuse_noncommercial_weights=True`
(or the CONTROL/tool option) refuses non-commercially licensed weights.

## References

See the DOI table above; each adapter table also carries its source string
in `ELEMENT_TABLE_SOURCES` / `ASE_SOURCES`, so a refusal message cites the
paper it comes from.

Every quantitative statement in this skill is a method-characterisation
claim backed by the reference next to it, not a measured DELFIN number.
