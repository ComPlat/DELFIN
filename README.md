# DELFIN

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17208145.svg)](https://doi.org/10.5281/zenodo.17208145)

**Automated DFT-based prediction of preferred spin states and associated redox potentials**

This repository contains DELFIN, a comprehensive workflow tool for automated quantum chemistry calculations using ORCA, xTB, and CREST. DELFIN automates the identification of preferred electron configurations, tracks orbital occupation changes during redox processes, and calculates redox potentials.

## 🚀 Quick Install

```bash
pip install delfin-complat
```

**Requirements:**
- **Python 3.10+**
- **ORCA 6.1.0** in your `PATH` ([free for academic use](https://orcaforum.kofo.mpg.de/app.php/portal))
- **Optional:** CREST, xTB (for extended workflows)

> **Prereqs**
>
> * ORCA **6.1.0** in your `PATH` (`orca` and `orca_pltvib`)
> * Optional: `crest` (for CREST workflow), **xTB** if used (`xtb` and `crest` in `PATH`)
> * Python **3.10+** required

---

## Install

**PyPI Package**: https://pypi.org/project/delfin-complat/

From the `delfin` folder (the one containing `pyproject.toml`):

recommended (isolated) install
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```
regular install
```bash
pip install delfin-complat
```

All Python dependencies (for example `mendeleev` for covalent radii) are installed automatically. Using a virtual environment or tools such as `pipx` keeps the scientific software stack reproducible and avoids system-wide modifications.

This exposes the console command **`delfin`** and enables `python -m delfin`.

---

## Quick start

Create a working directory with at least these two files:

* `CONTROL.txt` — your control/config file
* `input.txt` — the starting geometry (XYZ body without the first two header lines)
* starting from a `XYZ` file is optional

Then run:

from the directory that contains CONTROL.txt and input.txt
```bash
delfin
```
alternatively
```bash
python -m delfin
```

**CLI shortcuts**

- `delfin --define[=input.xyz] [--overwrite]`
  creates/updates `CONTROL.txt` and optionally converts an XYZ into `input.txt`.
- `delfin --control /path/to/CONTROL.txt`
  runs the workflow from another directory while normalising all paths.
- `delfin --no-cleanup`
  keeps temporary files and scratch folders after the pipeline finishes.
- `delfin --cleanup`
  removes previously generated intermediates and exits immediately.
- `delfin cleanup --orca`
  stops running ORCA jobs in the current workspace, purges OCCUPIER scratch folders, and cleans leftover temporary files.
- `delfin --purge`
  clears the working directory (keeps CONTROL.txt and the configured input file only) after confirmation.
- `delfin --recalc`
  re-parses existing results and only restarts missing or incomplete jobs.
- `delfin --help`
  prints the full list of CLI flags, including the new pipeline/resource switches.

Results and reports are written to the current working directory,
e.g. `DELFIN.txt`, `OCCUPIER.txt`, and per-step folders.
---

## Excited-State Dynamics (ESD) Module

> **Status:** The ESD module is under construction and not yet production-ready.

The optional ESD pipeline optimises S0/S1/T1/T2 states and launches intersystem
crossing (ISC) / internal conversion (IC) calculations inside a dedicated `ESD/`
folder. Enable it via `ESD_modul=yes` in `CONTROL.txt`. It runs after the normal
redox workflows when `method` is set, or on its own when `method` is left empty.

Minimal CONTROL settings:

```ini
ESD_modul=yes
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1
ICs=S1>S0
```

Geometries, GBW and Hessian files are staged under `ESD/` and reused automatically.
All states inherit the system charge and use multiplicity 1 (triplets employ UKS internally).

> Note: the %ESD `DELE` entry is currently omitted; ORCA falls back to its internal defaults.
> Temperature is still taken from `temperature` (defaults to 298.15 K).

Directory snapshot:

```
ESD/
  S0.inp/.out/.xyz/.gbw/.hess
  S1.inp/.out/.xyz/.gbw/.hess
  T1.inp/.out/.xyz/.gbw/.hess
  T2.inp/.out/.xyz/.gbw/.hess
  S1_T1_ISC.inp/.out
  ...
```

---

## Project layout

```
delfin/
  __init__.py
  __main__.py       # enables `python -m delfin`
  cli.py            # main CLI entry point orchestrating the full workflow
  cli_helpers.py    # CLI argument parsing and helper functions
  cli_recalc.py     # recalc mode wrapper functions for computational tools
  cli_banner.py     # banner display and file validation utilities
  cli_calculations.py # redox potential calculation methods (M1, M2, M3)
  main.py           # optional small loader (may delegate to cli.main)
  pipeline.py       # high-level orchestration across workflow phases (classic/manually/OCCUPIER)
  esd_module.py     # optional excited-state dynamics workflow (states, ISC/IC scheduling)
  esd_input_generator.py # ORCA input builders for ESD states/ISC/IC jobs
  config_manager.py # shared configuration utilities used by new pipeline helpers
  safe.py           # lightweight sandbox helpers for robust filesystem ops
  define.py         # CONTROL template generator (+ .xyz → input.txt conversion, path normalisation + logging hooks)
  cleanup.py        # delete temporary files
  config.py         # CONTROL.txt parsing & helpers
  utils.py          # common helpers (transition metal scan, basis-set selection, electron counts)
  orca.py           # ORCA executable discovery & runs
  imag.py           # IMAG workflow (plotvib helpers, imaginary-mode loop, freq-first order for optional output blocks)
  xyz_io.py         # XYZ/ORCA-input read/write helpers (freq block comes before any optional %output sections)
  xtb_crest.py      # xTB / GOAT / CREST / ALPB solvation workflows
  energies.py       # extractors for energies (FSPE, Gibbs, ZPE, electronic energies)
  parser.py         # parser utilities for ORCA output files
  occupier.py       # OCCUPIER workflow (sequence execution + summary)
  copy_helpers.py   # file passing between OCCUPIER steps (prepare/copy/select)
  thread_safe_helpers.py  # thread-safe workflow execution with PAL management
  global_manager.py       # singleton global job manager for resource coordination
  dynamic_pool.py         # dynamic core pool for job scheduling
  parallel_classic_manually.py     # parallel execution for classic/manually modes
  parallel_occupier.py  # parallel OCCUPIER workflow integration
  verify_global_manager.py  # smoke tests for the global resource orchestration
  cluster_utils.py        # cluster resource detection (SLURM/PBS/LSF)
  api.py            # programmatic API (e.g. `delfin.api.run(...)` for notebooks/workflows)
  common/           # shared utilities
    __init__.py     # exposes common helpers
    banners.py      # CLI banner art + static strings
    logging.py      # logging configuration/get_logger helpers (cluster-friendly)
    orca_blocks.py  # reusable ORCA block assembly utilities
    paths.py        # central path & scratch-directory helpers (`DELFIN_SCRATCH` aware)
  reporting/        # modular report generation
    __init__.py     # reporting submodule exports
    occupier_reports.py  # OCCUPIER-specific report generation functions
    delfin_reports.py    # DELFIN-specific report generation functions
    occupier_selection.py # OCCUPIER selection helpers used by reports
```
---
## Typical workflow switches (in CONTROL.txt)

* `method = OCCUPIER | classic | manually` (leave empty for ESD-only runs)
* `calc_initial = yes | no`
* `oxidation_steps = 1,2,3` (string; steps to compute)
* `reduction_steps = 1,2,3` (string; steps to compute)
* `E_00 = yes | no`
* `absorption_spec = yes | no`
* `parallel_workflows = yes | no | auto` (parallelization)
* `pal_jobs = N` (number of parallel PAL processes; auto-detected from cluster if not set)
* `orca_parallel_strategy = auto | threads | serial` (ORCA parallelization mode)
  - `auto` (default): Use MPI + OpenMP; FoBs run in parallel
  - `threads`: Use only OpenMP (no MPI); FoBs still run in parallel, useful if MPI is unstable
  - `serial`: Force sequential execution; only 1 FoB at a time
* `IMAG_scope = initial | all`
  - `initial` (default): only the initial geometry will run IMAG elimination.
  - `all`: the IMAG workflow is executed for all configured redox steps.
