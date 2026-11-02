# Machine-learned potentials (MLP) and ASE steps in DELFIN

Package ML1 (wave 21). DELFIN's MLP/ASE adapters are **physically realistic
workflow steps**: every limit a model has — supported elements, charge and
spin multiplicity, level-of-theory validity — is *declared in a table with a
verified source and enforced fail-closed*, not just documented.

## Steps

| Step | Module | Purpose |
|---|---|---|
| `mlp_single_point` | `delfin/tools/adapters/mlp.py` | SP energy/forces via an ML potential |
| `mlp_optimize` | `delfin/tools/adapters/mlp.py` | geometry optimization via an ML potential |
| `mlp_rank` | `delfin/tools/adapters/mlp.py` | rank a conformer ensemble by MLP energy |
| `ase_optimize` / `ase_single_point` / `ase_md` | `delfin/tools/adapters/ase_adapter.py` | any ASE calculator (xtb, orca, emt, eam, lj, mace, ani) |
| `mlp_screen` workflow | `delfin/workflows/contrib/w21_mlp_screen.py` | end-to-end ensemble ranking with realism enforced |

## Declared limits (the enforcement contract)

ML potential backends (`SUPPORTED_ELEMENTS`, `ELEMENT_TABLE_SOURCES`,
`VALIDITY` in `delfin/tools/adapters/mlp.py`):

- **ani2x** — {H,C,N,O,S,F,Cl}, neutral closed shell. Source: Devereux et
  al., *J. Chem. Theory Comput.* 16 (2020) 4192–4202,
  DOI 10.1021/acs.jctc.0c00121.
- **aimnet2** — {H,B,C,N,O,F,Si,P,S,Cl,As,Se,Br,I}; the only
  charge-aware, spin-aware backend. Source: Anstine, Zubatyuk, Isayev,
  *Chem. Sci.* 16 (2025) 10228–10244, DOI 10.1039/D4SC08572H.
- **mace_off** — {H,C,N,O,F,P,S,Cl,Br,I}, neutral organic. Source: Kovacs
  et al., *J. Am. Chem. Soc.* 147 (2025) 9254, DOI 10.1021/jacs.4c07099.
- **chgnet, m3gnet, schnetpack, nequip, alignn** — element scope is
  decided by the trained model file (`None` in the table). A structural
  element refusal would be wrong physics for these; charge/multiplicity
  is still refused, because the calculator wrapper ignores it.

ASE calculators (`ASE_SUPPORTED_ELEMENTS`, `ASE_SOURCES`, `ASE_VALIDITY`
in `delfin/tools/adapters/ase_adapter.py`):

- **xtb** — all spd-block elements up to Z=86, molecular systems; charge
  and multiplicity supported. Source: Grimme, Bannwarth, Shushkov,
  *J. Chem. Theory Comput.* 13 (2017) 1989–2009, DOI
  10.1021/acs.jctc.7b00118.
- **orca** — full periodic table, charge/multiplicity passed to ORCA.
- **emt/eam** — EMT standard metals Al,Cu,Ag,Au,Ni,Pd,Pt plus H,C,N,O
  ("not for any serious use", installed ASE docstring); EAM scope set by
  its potential file.
- **lj, mace** — generic/model-defined element scope; no element refusal
  (wrong physics), but charge/multiplicity is refused because the
  constructor silently ignores it.
- **ani (ASE)** — same ANI-2x scope as above.

## How the enforcement works

1. **Value check** (`validate_params`): the backend/calculator name is
   checked against the factory's alias set (case-insensitive, aliases
   accepted); unknown, empty or non-string values raise `ValueError`
   before anything runs.
2. **Charge/multiplicity check** (`validate_params`): nonzero charge or
   multiplicity ≠ 1 on a backend that cannot model it raises
   `ValueError` — the calculator would otherwise *silently ignore* it
   and return a physically wrong energy.
3. **Structure check** (`execute`): every structure is checked against
   the backend's declared element set *on the actual structure*, before
   any model is created; the refusal names the unsupported element(s)
   and cites the source. Model-file-scoped backends are exempt.
4. **Failure surfacing**: all refusals surface as a FAILED `StepResult`
   with the diagnosis (including missing ML-framework installs), never a
   silently wrong number and never a swallowed error.

## The `mlp_screen` workflow

```python
from delfin.workflows.registry import get
wf = get("mlp_screen")
out = wf.run(config={}, ensemble="conformers.xyz",
             backend="ani2x", charge=0, mult=1)
out["best_index"], out["best_energy_eV"], out["ranking"]
```

The workflow re-uses the adapter tables (single source of truth — it
carries no copied table) and refuses invalid requests before reading the
ensemble. Discoverable through the shared workflow registry
(`delfin/workflows/registry.py` auto-discovery of `contrib/`).

## Agent tools (MCP)

The same evaluation is exposed to the DELFIN agent as five MCP tools
(`delfin/api.py` handlers, registered in `delfin/ops_server/server.py`):

| Tool | Read-only? | Purpose |
|---|---|---|
| `mlp_list_backends` | yes (gate allowlist) | availability, element scope, install hint per backend |
| `mlp_single_point` | no (side-effect gate) | SP energy/forces for an xyz string |
| `mlp_optimize` | no | optimise an xyz; returns the geometry |
| `mlp_rank` | no | rank a multi-structure xyz ensemble |
| `mlp_preoptimize` | no | pre-optimise before ORCA; geometry is a starting guess |

All take xyz text and answer JSON. Failures are diagnoses, not
tracebacks: an unknown backend answers `ok=false` with the accepted
names; a missing backend package answers `ok=false`,
`available=false`, `reason`, `install_hint` (see
`delfin/mlp_tools/SETUP_GUIDE.txt`). The four run tools compute, so they
start denied by the agent's side-effect gate like every calculation
runner; `mlp_list_backends` passes the gate as read-only.

## The `mlp_preopt` CONTROL key (opt-in)

`mlp_preopt = none | ani2x | aimnet2 | mace_off | ...` relaxes the start
geometry with the chosen ML potential BEFORE the ORCA steps. The default
`none` is the shipped behaviour — nothing runs and the generated inputs
are byte-identical. Semantics: a screening step; the relaxed geometry is
a starting guess for the DFT optimisation, and every final energy comes
from DFT. A failed or non-converged MLP step never stops the run (the
original geometry is kept, the run continues, the log says why). Use it
for element-clean, closed-shell systems whose builder geometry is poor;
for metals, open-shell or charged systems keep it off unless the backend
models charge and multiplicity (only `aimnet2`). Declared in
`CONTROL_FIELD_SPECS` (`delfin/common/control_validator.py`), notes under
key `mlp_preopt` in `delfin/common/control_key_notes.json`.

## When to use an MLP, and when not to

- Use for: pre-optimization before DFT, conformer ranking, cheap scans
  over many candidates. Label results "MLP pre-screen".
- Do not use for: publication-grade energies, sub-kcal/mol barriers,
  spin-state splittings, charged species (outside aimnet2), periodic
  solids with molecular-scope models. A minimum is verified only by a
  Hessian with zero imaginary frequencies — never by an MLP optimization
  alone.
