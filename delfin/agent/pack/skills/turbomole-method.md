---
name: turbomole-method
domains: code
description: >
  Run Turbomole through DELFIN's adapters and the contributed turbomole
  workflow: pass charge/multiplicity correctly ($unpaired), choose RI vs
  conventional DFT, know which binaries are needed (define/x2t optional,
  ridft/jobex/aoforce required), and know Turbomole's limits — it is
  license-restricted and never auto-installed. Read this before running
  any turbomole_* adapter or the turbomole workflow.
---

# Running Turbomole in DELFIN

DELFIN drives Turbomole three ways — never hand-write a raw control
file or shell command:

- **Adapters** (`delfin/tools/adapters/turbomole.py`): `turbomole_sp`
  (ridft/dscf), `turbomole_opt` (jobex), `turbomole_freq` (ridft +
  aoforce). Parameters: `charge` (required, int), `mult` (2S+1, int,
  >= 1), `method` (Turbomole functional keyword, default `b3-lyp`),
  `basis` (default `def2-SVP`), `ri` (bool, default true), `timeout`,
  `max_cycles` (opt only, positive int).
- **Contributed workflow** (`delfin/workflows/contrib/
  w21_turbomole_workflow.py`): `turbomole` — XYZ geometry in, opt +
  freq out, with energy, optimized geometry, hessian and an explicit
  imaginary-frequency count in the result.
- **Guidance only**: Turbomole is commercial, licensed software
  (`delfin/tools/_environment.py`) — DELFIN never installs it.

## Charge, multiplicity and open shells

- `charge` and `mult` must be plain integers; `mult >= 1`. A bool or
  string is rejected before anything runs.
- An open shell (`mult > 1`) writes `$unpaired N` (N = mult − 1) into
  the control file and the define input. A singlet writes none. A
  closed-shell run of what you asked as a triplet is a defect, not a
  fallback — if the run dies on `$unpaired`, fix the multiplicity, do
  not delete the directive.
- The electron-count/multiplicity parity rule (WB1, 
  `delfin/tools/_science_rules/electronic.py`) applies before the run:
  nelec ≡ mult − 1 (mod 2). Turbomole will stop a physically impossible
  combination; the lint flags it before compute is spent.

## RI or not

`ri=true` (default) runs `ridft` with the RI-J approximation — the
right default for DFT with a fitting basis. `ri=false` runs `dscf`
(conventional) — needed when no auxiliary basis is set or when RI is
not supported for the chosen functional. The define path writes the RI
block only when `ri` is true.

## Binaries

- Required at run time: `ridft`/`dscf`, `jobex` (opt), `aoforce` (freq).
- Optional: `define` (interactive control generator — without it the
  adapter writes a minimal control file directly), `x2t` (xyz ↔ coord
  conversion — without it a built-in converter is used).
- If a binary is missing the adapter returns an honest FAILED result
  naming the missing program; it does not raise and does not install.

## What a run must report (never hide)

- The freq step surfaces the aoforce `hessian` file path in
  `result.data["hessian"]`; if aoforce exits 0 without writing one the
  step FAILS with that message.
- The `turbomole` workflow counts imaginary frequencies from the
  aoforce log and reports them in `imaginary_frequencies` plus a
  warning when > 0 — a non-minimum is reported, never silently passed.
