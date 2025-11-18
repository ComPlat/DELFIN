# DELFIN

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17208145.svg)](https://doi.org/10.5281/zenodo.17208145)

> 📄 **Preprint**: *Hartmann, M. et al. “DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials”*, ChemRxiv (2025). https://doi.org/10.26434/chemrxiv-2025-4c256

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
- `delfin cleanup [--dry-run] [--workspace PATH] [--scratch PATH]`
  offers finer control over workspace/scratch cleanup; combine with `--dry-run` to preview deletions.
- `delfin cleanup --orca`
  stops running ORCA jobs in the current workspace, purges OCCUPIER scratch folders, and cleans leftover temporary files.
- `delfin --purge`
  removes DELFIN-generated artifacts (OCCUPIER folders, ORCA inputs/outputs, logs) after confirmation while keeping CONTROL.txt, the configured input file, and any unrelated files.
- `delfin --recalc`
  re-parses existing results and only restarts missing or incomplete jobs.
- `delfin --report`
  re-calculates redox potentials from existing output files without launching new calculations.
- `delfin --imag`
  eliminates imaginary modes from existing ORCA results (`*.out`/`*.hess`) and regenerates the summary report.
- `delfin --version`
  prints the installed DELFIN version (`-V` shortcut).
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
  occupier_auto.py  # Auto OCCUPIER sequence management and tree navigation
  deep_auto_tree.py # Deep tree: adaptive BS evolution (reduction: BS(m-1,1) or BS(M±1,N); oxidation: pure only)
  deep2_auto_tree.py # Deep2 tree: only pure states (no BS), simple 3×3 branching
  generate_deep2_tree.py # Generator for deep2 (pure states only)
  generate_deep_tree.py # Generator for deep tree with adaptive BS evolution
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
* `OCCUPIER_method = auto | manually` (auto uses adaptive tree-based sequences)
* `OCCUPIER_tree = flat | deep2 | deep3 | deep | own` (auto tree mode; `own` builds a custom adaptive tree from the CONTROL `OCCUPIER_sequence_profiles` block)
  - `flat`: Legacy flat sequences with BS
  - `deep2`: Only pure states (no BS), simple testing
  - `deep`: Adaptive BS evolution (recommended)
    - Reduction: BS(m-1,1) when pure wins; BS(M±1,N), BS(M,N±1) when BS wins
    - Oxidation: Only pure states (no BS)
    - Depth: 3 levels (0 → ±1 → ±2 → ±3)
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
* `XTB_OPT = yes | no`
* `XTB_GOAT = yes | no`
* `CREST = yes | no`
* `XTB_SOLVATOR = yes | no`
* `ESD_modul = yes | no`
  * `states = S0,S1,T1,T2` (comma-separated subset)
  * `ISCs = S1>T1,...` (optional)
  * `ICs = S1>S0,...` (optional)
* `IMAG = yes | no`
  - enable imaginary-frequency elimination for pipeline runs and unlock `delfin --imag`
* `IMAG_scope = initial | all`
  - `initial` (default): only the starting geometry runs the IMAG loop; `all` processes every redox step
* `allow_imaginary_freq = N`
  - number of imaginary modes tolerated before IMAG restarts the step (default: 0)
* `IMAG_displacement_scale = float`
  - optional scaling factor for the `orca_pltvib` displacement amplitude (default: 1.0)
* `IMAG_sp_energy_window = float`
  - minimum single-point energy improvement (Hartree) required to accept a displaced geometry (default: 1e-6)
* `IMAG_optimize_candidates = yes | no`
  - when `yes`, each displaced structure is geometry-optimised before evaluating the single-point energy (default: `no`)

---

## Cluster & Workflow Integration

* **Scratch directory:** set `DELFIN_SCRATCH=/path/to/scratch` before launching jobs. Temporary files, markers, and runtime artefacts are written there (directories are created automatically).
* **Schema validation:** `delfin` validates `CONTROL.txt` on load (missing required keys, wrong types, inconsistent sequences) and aborts with a clear error message if something is off.
* **Logging:** the CLI now writes a timestamped `delfin_run.log` in every working directory and OCCUPIER subprocesses produce an additional `occupier.log` inside their folders. Custom drivers can still call `delfin.common.logging.configure_logging(level, fmt, stream)` to adjust handlers or formats.
* **Programmatic API:** use `delfin.api.run(control_file="CONTROL.txt")` for notebooks, workflow engines, or SLURM batch scripts. Add `cleanup=False` to preserve intermediates (`--no-cleanup`). Additional CLI flags can be provided through the `extra_args` parameter.
* **Alternate CONTROL locations:** supply `--control path/to/CONTROL.txt` (or the `control_file` argument in `delfin.api.run`) to stage input files outside the working directory.
* **XYZ geometry support:** if `input_file` in CONTROL (or the CLI/API) points to an `.xyz`, DELFIN converts it automatically to a matching `.txt` (header dropped) before the run.
* **Cluster templates:** see `examples/` for submit scripts:
  - `slurm_submit_example.sh` (SLURM)
  - `pbs_submit_example.sh` (PBS/Torque)
  - `lsf_submit_example.sh` (LSF)
* **Auto-resource detection:** DELFIN automatically detects available CPUs and memory on SLURM/PBS/LSF clusters and configures PAL/maxcore accordingly if not explicitly set in CONTROL.txt.

### Global Resource Management

DELFIN uses a **global job manager singleton** to coordinate all computational workflows and ensure that CPU resources (PAL) are never over-allocated:

* **Single source of truth:** PAL is read once from `CONTROL.txt` at startup and managed centrally throughout execution
* **Automatic PAL splitting:** When oxidation and reduction workflows run in parallel, DELFIN automatically splits available cores between them (e.g., PAL=12 → 6 cores per workflow)
* **Thread-safe execution:** All parallel workflows coordinate through a shared resource pool, preventing race conditions
* **Subprocess coordination:** OCCUPIER subprocesses receive their allocated PAL via environment variables and respect global limits
* **Sequential mode:** When workflows run sequentially (`parallel_workflows=no`), each workflow uses the full PAL

This architecture ensures:
- No double allocation of cores when ox/red workflows run simultaneously
- Consistent resource limits across all ORCA jobs spawned by DELFIN
- Proper coordination between main process and OCCUPIER subprocesses
- Efficient utilization of cluster resources without exceeding allocation

## Troubleshooting

* **`CONTROL.txt` not found**
  DELFIN exits gracefully and tells you what to do. Create it via `delfin --define` (or copy your own).

* **Input file not found**
  DELFIN exits gracefully and explains how to create/convert it.
  If you have a full `.xyz`, run: `delfin --define=your.xyz` → creates `input.txt` (drops the first two header lines) and sets `input_file=input.txt` in CONTROL.

* **ORCA not found**
  Ensure `orca` is callable in your shell: `which orca` (Linux/macOS) or `where orca` (Windows).
  Add the ORCA bin directory to your `PATH`.

* **`ModuleNotFoundError` for internal modules**
  Reinstall the package after copying files:


* **CREST/xTB tools missing**
  Disable the corresponding flags in `CONTROL.txt` or install the tools and put them in `PATH`.

---

## Dev notes

* Update CLI entry point via `pyproject.toml`
  `"[project.scripts] delfin = \"delfin.cli:main\""`
* Build a wheel: `pip wheel .` (inside `delfin/`).
* Run tests/workflow locally using a fresh virtual environment to catch missing deps.

---
## References

The generic references for ORCA, xTB and CREST are:

- Frank Neese. The ORCA program system. *Wiley Interdiscip. Rev. Comput. Mol. Sci.*, 2(1):73–78, 2012. doi:<https://doi.wiley.com/10.1002/wcms.81>.
- Frank Neese. Software update: the ORCA program system, version 4.0. *Wiley Interdiscip. Rev. Comput. Mol. Sci.*, 8(1):e1327, 2018. doi:<https://doi.wiley.com/10.1002/wcms.1327>.
- Frank Neese, Frank Wennmohs, Ute Becker, and Christoph Riplinger. The ORCA quantum chemistry program package. *J. Chem. Phys.*, 152(22):224108, 2020. doi:<https://aip.scitation.org/doi/10.1063/5.0004608>.
- Christoph Bannwarth, Erik Caldeweyher, Sebastian Ehlert, Andreas Hansen, Philipp Pracht, Jan Seibert, Sebastian Spicher, and Stefan Grimme. Extended tight-binding quantum chemistry methods. *WIREs Comput. Mol. Sci.*, 11:e1493, 2021. doi:<https://doi.org/10.1002/wcms.1493>. *(xTB & GFN methods)*
- Philipp Pracht, Stefan Grimme, Christoph Bannwarth, Florian Bohle, Sebastian Ehlert, Gunnar Feldmann, Jan Gorges, Max Müller, Timo Neudecker, Christoph Plett, Sebastian Spicher, Pascal Steinbach, Piotr A. Wesołowski, and Fabian Zeller. CREST — A program for the exploration of low-energy molecular chemical space. *J. Chem. Phys.*, 160:114110, 2024. doi:<https://doi.org/10.1063/5.0197592>. *(CREST)*

Please always check the output files—at the end, you will find a list of relevant papers for the calculations. Kindly cite them. Please do not only cite the above generic references, but also cite in addition the
[original papers](https://www.faccts.de/docs/orca/6.0/manual/contents/public.html) that report the development and ORCA implementation of the methods DELFIN has used! The publications that describe the functionality implemented in ORCA are
given in the manual.



# Dependencies and Legal Notice

**DISCLAIMER: DELFIN is a workflow tool that interfaces with external quantum chemistry software. Users are responsible for obtaining proper licenses for all required software.**

## ORCA Requirements
To use DELFIN, you must be authorized to use ORCA 6.1.0. You can download the latest version of ORCA here:
https://orcaforum.kofo.mpg.de/app.php/portal

***IMPORTANT: ORCA 6.1.0 requires a valid license and registration. Academic users can obtain free access, but commercial use requires a commercial license. Please carefully review and comply with ORCA's license terms before use.***
https://www.faccts.de/

**ORCA License Requirements:**
- Academic use: Free after registration and license agreement
- Commercial use: Requires commercial license
- Users must register and agree to license terms before downloading
- Redistribution of ORCA is prohibited
- Each user must obtain their own license
- DELFIN does not include or distribute ORCA
- ORCA is proprietary software owned by the Max Planck Institute for Coal Research
- End users must comply with ORCA's terms of service and usage restrictions
- DELFIN authors are not affiliated with or endorsed by the ORCA development team

## xTB Requirements
***xTB is free for academic use under the GNU General Public License (GPLv3).***
The code and license information are available here: https://github.com/grimme-lab/xtb
- Commercial use may require different licensing terms
- DELFIN does not include or distribute xTB

## CREST Requirements
***CREST is free for academic use under the GNU General Public License (GPLv3).***
The code and license information are available here: https://github.com/crest-lab/crest
- Commercial use may require different licensing terms
- DELFIN does not include or distribute CREST

**Legal Notice:** DELFIN itself is licensed under LGPL-3.0-or-later, but this does not grant any rights to use ORCA, xTB, or CREST. Users must comply with the individual license terms of each external software package.

## Warranty and Liability
DELFIN is provided "AS IS" without warranty of any kind. The authors disclaim all warranties, express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose. In no event shall the authors be liable for any damages arising from the use of this software.

---

## Please cite

If you use DELFIN in a scientific publication, please cite:

- Hartmann, M. (2025). *DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials* (v1.0.4). Zenodo. https://doi.org/10.5281/zenodo.17208145
- Hartmann, M. (2025). *DELFIN: Automated prediction of preferred spin states and redox potentials*. ChemRxiv. https://chemrxiv.org/engage/chemrxiv/article-details/68fa0e233e6156d3be78797a

### BibTeX
```bibtex
@software{hartmann2025delfin,
  author  = {Hartmann, Maximilian},
  title   = {DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials},
  version = {v1.0.4},
  year    = {2025},
  publisher = {Zenodo},
  doi     = {10.5281/zenodo.17208145},
  url     = {https://doi.org/10.5281/zenodo.17208145}
}

@article{hartmann2025chemrxiv,
  author  = {Hartmann, Maximilian},
  title   = {DELFIN: Automated prediction of preferred spin states and redox potentials},
  journal = {ChemRxiv},
  year    = {2025},
  url     = {https://chemrxiv.org/engage/chemrxiv/article-details/68fa0e233e6156d3be78797a}
}
```
## License

This project is licensed under the GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).

You should have received a copy of the GNU Lesser General Public License along with this repository in the files `COPYING` and `COPYING.LESSER`.  
If not, see <https://www.gnu.org/licenses/>.

Non-binding citation request:  
If you use this software in research, please cite the associated paper (see [CITATION.cff](./CITATION.cff)).


  
