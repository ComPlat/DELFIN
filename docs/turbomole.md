# Turbomole in DELFIN

Turbomole is commercial, licensed software. DELFIN never installs it
(`delfin/tools/_environment.py`); if it is present on the cluster or
workstation, DELFIN drives it through the step adapters and the
contributed workflow described here.

## Step adapters (`delfin/tools/adapters/turbomole.py`)

| Adapter | Program(s) | Purpose |
|---|---|---|
| `turbomole_sp` | `ridft` (RI) / `dscf` | Single-point energy |
| `turbomole_opt` | `jobex` | Geometry optimization (max `max_cycles`, default 200) |
| `turbomole_freq` | `ridft` + `aoforce` | Frequencies + Hessian (minimum check) |

Parameters (shared): `charge` (required, int), `mult` (2S+1, int ≥ 1),
`method` (Turbomole functional keyword, default `b3-lyp`), `basis`
(default `def2-SVP`), `ri` (bool, default true), `timeout` (s).
`turbomole_opt` adds `max_cycles`.

### Physical limits — declared and enforced

- `charge`/`mult` must be plain integers (bool rejected); `mult ≥ 1`.
  Violations return a FAILED result naming the offending value before
  anything runs.
- An open shell (`mult > 1`) is declared explicitly: `$unpaired N`
  (N = mult − 1) in the control file (manual path) and the define
  input. A singlet writes no `$unpaired` directive. The define input
  writes the `ri` block only when `ri` is true.
- Electron-count/multiplicity parity (WB1 rule 10,
  `delfin/tools/_science_rules/electronic.py`) is flagged by
  `scientific_lint` before compute is spent: nelec ≡ mult − 1 (mod 2).
- A missing binary returns an honest FAILED result naming the missing
  program (no exception, no install attempt).
- `turbomole_freq` surfaces the aoforce `hessian` file in
  `result.data["hessian"]`; if aoforce exits 0 without writing one the
  step FAILS with that message.

### Security posture

All Turbomole invocations run through an argv list without a shell
(`subprocess.run(argv, shell=False)`); `max_cycles` is validated as a
positive integer before it reaches the process, so no user-supplied
string can reach a shell (no injection surface).

## Contributed workflow (`delfin/workflows/contrib/w21_turbomole_workflow.py`)

Workflow name `turbomole`: takes an XYZ geometry, runs
`turbomole_opt` then `turbomole_freq`, and returns energy, optimized
geometry, hessian path, and the count of imaginary frequencies parsed
from the aoforce log — reported with a warning when > 0 (a
non-minimum is never hidden). There is no SMILES step: the workflow is
geometry-in; build a geometry with `smiles_to_xyz` first if needed.

## Discovery and availability

- Tool discovery: `TURBODIR` / `TURBOMOLE_HOME` / `TURBOMOLE_ROOT`
  environment variables, system dirs and module patterns
  (`delfin/qm_runtime.py`).
- License-restricted: never auto-installed; guidance only.

## Tests

- `tests/test_tm1_limits.py` — physical limits, open-shell directive,
  hessian surfacing, missing-binary honesty (10 cases).
- `tests/test_tm1_workflow.py` — registry discovery, WB2 check chain,
  workflow run/limits.
- `tests/test_tm1_adversarial_turbomole.py` — reviewer adversarial
  suite (nacht-s25).
