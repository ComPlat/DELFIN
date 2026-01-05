# DELFIN

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17208145.svg)](https://doi.org/10.5281/zenodo.17208145)
[![PyPI version](https://img.shields.io/pypi/v/delfin-complat.svg)](https://pypi.org/project/delfin-complat/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/delfin-complat.svg)](https://pypistats.org/packages/delfin-complat)

> 📄 **Preprint**: *Hartmann, M. et al. “DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials”*, ChemRxiv (2025). https://doi.org/10.26434/chemrxiv-2025-4c256

This repository contains DELFIN, a comprehensive workflow tool for automated quantum chemistry calculations using ORCA, xTB, and CREST. DELFIN automates the identification of preferred electron configurations, tracks orbital occupation changes during redox processes, and calculates redox potentials.

## 🚀 Installation

**PyPI Package**: https://pypi.org/project/delfin-complat/

### Requirements
- **Python 3.10+**
- **ORCA 6.1.1** in your `PATH` (`orca` and `orca_pltvib`) — [free for academic use](https://orcaforum.kofo.mpg.de/app.php/portal)
- **Optional:** `crest` and `xtb` (for CREST/xTB workflows)

### Install Methods

**Standard install (recommended for most users):**
```bash
pip install delfin-complat
```

**Development install (from source):**
```bash
# Clone the repository
git clone https://github.com/ComPlat/DELFIN.git
cd DELFIN

# Create virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in editable mode
pip install -e .
```

All Python dependencies (e.g., `mendeleev` for covalent radii) are installed automatically. Using a virtual environment keeps the scientific software stack reproducible and avoids system-wide modifications.

This exposes the console command **`delfin`** and enables `python -m delfin`.

---

## Quick start

Create a working directory with at least these two files:

* `CONTROL.txt` — your control/config file
* `input.txt` — the starting geometry (XYZ body without the first two header lines)
* starting from a `XYZ` file is optional

Then run:

from the directory that contains `CONTROL.txt` and `input.txt`
```bash
delfin
```
alternatively
```bash
python -m delfin
```

You can also point DELFIN at a different workspace directory:
```bash
delfin /path/to/project
```

Results and reports are written to the current working directory,
e.g. `DELFIN.txt`, `OCCUPIER.txt`, and per-step folders.

---

## 📋 CLI Reference

### Basic Commands

- `delfin` — Run DELFIN workflow in current directory (requires `CONTROL.txt` and `input.txt`)
- `delfin /path/to/project` — Run DELFIN workflow in specified directory
- `python -m delfin` — Alternative way to run DELFIN

### Setup & Configuration

- `delfin --define[=input.xyz] [--overwrite]` — Create/update `CONTROL.txt` and optionally convert XYZ to `input.txt`
- `delfin /path/to/project --define[=input.xyz] [--overwrite]` — Same as above, but in a different workspace directory
- `delfin --control /path/to/CONTROL.txt` — Run workflow from another directory while normalising all paths

### Execution Control

- `delfin --recalc` — Re-parse existing results and only restart missing or incomplete jobs
- `delfin WORKSPACE --recalc --occupier-override STAGE=INDEX` — Force a specific OCCUPIER index for a stage during recalc (e.g., `--occupier-override red_step_2=1` uses index 1 for red_step_2_OCCUPIER instead of the automatically selected preferred index). Can be passed multiple times for different stages
- `delfin --report` — Re-calculate redox potentials from existing output files without launching new calculations
- `delfin --imag` — Eliminate imaginary modes from existing ORCA results (`*.out`/`*.hess`) and regenerate the summary report
- `delfin stop --workspace PATH` — Send a graceful stop signal to running DELFIN processes associated with a workspace

### Cleanup & Maintenance

- `delfin --no-cleanup` — Keep temporary files and scratch folders after the pipeline finishes
- `delfin --cleanup` — Remove previously generated intermediates and exit immediately
- `delfin cleanup [--dry-run] [--workspace PATH] [--scratch PATH]` — Fine control over workspace/scratch cleanup; combine with `--dry-run` to preview deletions
- `delfin cleanup --orca` — Stop running ORCA jobs in the current workspace, purge OCCUPIER scratch folders, and clean leftover temporary files
- `delfin --purge` — Remove DELFIN-generated artifacts (OCCUPIER folders, ORCA inputs/outputs, logs) after confirmation while keeping CONTROL.txt, the configured input file, and any unrelated files

### Export & Reporting

- `delfin --json` — Collect key results into `DELFIN_Data.json` (useful for scripts/notebooks)
- `delfin --report docx` — Generate a Word report in DOCX format
- `delfin --afp` — Generate an AFP spectrum plot (`AFP_spectrum.png`) from existing ESD/S0/S1/T1 results if present

### Information

- `delfin --version` (or `-V`) — Print the installed DELFIN version
- `delfin --help` — Print the full list of CLI flags, including the new pipeline/resource switches

### Specialized Workflows

- `delfin co2 ...` — CO2 Coordinator helper workflow for CO2 placement and scan generation/runs
- `delfin ESD` — Run excited-state dynamics workflow (requires `ESD_modul=yes` in CONTROL)

---

## ⚙️ Configuration (CONTROL.txt)

DELFIN is configured via `CONTROL.txt` in your working directory. Key settings include:

### Workflow Control
* `method = OCCUPIER | classic | manually` (leave empty for ESD-only runs)
* `OCCUPIER_method = auto | manually` (auto uses adaptive tree-based sequences)
* `OCCUPIER_tree = flat | deep2 | deep3 | deep | own` (auto tree mode)
* `calc_initial = yes | no`
* `oxidation_steps = 1,2,3` / `reduction_steps = 1,2,3`
* `parallel_workflows = yes | no | auto`

### Resource Management
* `pal_jobs = N` (number of parallel PAL processes; auto-detected from cluster if not set)
* `orca_parallel_strategy = auto | threads | serial`

### Optional Modules
* `XTB_OPT = yes | no` / `XTB_GOAT = yes | no` / `CREST = yes | no`
* `XTB_SOLVATOR = yes | no`
* `ESD_modul = yes | no` (excited-state dynamics)
  * `states = S0,S1,T1,T2`
  * `ISCs = S1>T1,...` / `ICs = S1>S0,...`
* `IMAG = yes | no` (imaginary frequency elimination)
  * `IMAG_scope = initial | all`
  * `allow_imaginary_freq = N`

### Error Recovery
* `enable_auto_recovery = yes | no` (intelligent ORCA error recovery with MOREAD)
* `max_recovery_attempts = N` (default: 1)
* `enable_job_timeouts = yes | no` (set to 'no' for unlimited runtime)

For the complete list of settings, see the "Typical workflow switches" section below.

---

## ✨ Features

### Excited-State Dynamics (ESD) Module

> **Status:** Under active development

The optional ESD pipeline optimizes S0/S1/T1/T2 states and calculates intersystem
crossing (ISC) / internal conversion (IC) inside a dedicated `ESD/` folder.

Enable via `ESD_modul=yes` in `CONTROL.txt`. Minimal settings:

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

---

## 🔧 Automatic Error Recovery & Retry System

**NEW**: DELFIN now includes an intelligent error recovery system that automatically detects and fixes common ORCA calculation failures.

### Features

- **Automatic Error Detection**: Identifies specific failure types (SCF convergence, TRAH crashes, geometry issues, etc.)
- **MOREAD-Based Continuation**: Uses last `.gbw` file to continue from last successful state (no restart from scratch!)
- **Progressive Escalation**: Applies increasingly aggressive fixes across retry attempts
- **Zero Manual Intervention**: Fully automatic recovery process
- **State Tracking**: Prevents infinite loops, tracks recovery history

### Quick Enable

**For 99% of users** - just enable this in your `CONTROL.txt`:

```ini
------------------------------------
Automatic Error Recovery & Retry:
enable_auto_recovery=yes            # ⭐ Intelligent recovery (RECOMMENDED!)
max_recovery_attempts=1             # Retry once (2 runs total)
------------------------------------
```

**Note**: Old parameter names (`orca_retry_enabled`, `orca_retry_max_attempts`) are still supported for backward compatibility.

**Advanced** - Optional configurations:
```ini
# Job timeout control:
enable_job_timeouts=yes             # Set to 'no' for unlimited runtime
job_timeout_hours=24
opt_timeout_hours=14
frequency_timeout_hours=36
sp_timeout_hours=3
```

💡 **What `auto_recovery` does**:
- Uses MOREAD to continue from last state (not from scratch!)
- Automatically updates geometry from latest .xyz file
- Creates GBW backup to prevent ORCA deletion on failure
- Modifies input with intelligent fixes based on error type
- Handles transient system errors with exponential backoff (2s, 4s, 8s, ...)
- Works for many cases (SCF, TRAH, geometry, MPI, memory, disk/network errors)

💡 **Timeout Management**:
- Set `enable_job_timeouts=no` for difficult systems requiring >24h runtime
- Jobs will run indefinitely until completion or error
- Useful for highly correlated calculations or large systems

### Supported Error Types

| Error | Automatic Fix |
|-------|--------------|
| **SCF not converged** | SlowConv → VerySlowConv + KDIIS|
| **TRAH segfault** | NoAutoTRAH |
| **DIIS errors** | Switch to KDIIS |
| **Geometry not converged** | Smaller trust radius → Loose criteria |
| **MPI crashes** | Reduce cores |
| **Memory errors** | Reduce maxcore and PAL |
| **Frequency failures** | Alternative methods or skip |
| **Transient system errors** | Exponential backoff retry (disk full, network timeout, I/O errors) |

### How It Works

```
ORCA fails → Detect error type → Modify input with fixes → Continue from last .gbw and xyz → Retry
```

**Example Recovery**:
```
Original:  complex.inp (SCF not converged)
    ↓
Recovery:  complex.retry1.inp (adds MOREAD + SlowConv)
    ↓
Success:   Converges using last wavefunction and xyz!
```

### Detailed Documentation

See **[docs/RETRY_LOGIC.md](docs/RETRY_LOGIC.md)** for:
- Complete configuration guide
- Recovery strategies for each error type
- Examples and troubleshooting
- Best practices

### When to Use

✅ **Enable auto_recovery when:**
- Running difficult systems (high spin, heavy metals, high redox states)
- Batch processing multiple jobs
- Using OCCUPIER or --recalc mode (automatically uses latest xyz/gbw)
- Geometry optimizations (automatically continues from last optimized structure)

✅ **Disable job_timeouts when:**
- Running very long calculations (>24h)
- Highly correlated methods (CCSD(T), NEVPT2, etc.)
- Large systems with slow SCF convergence

❌ **Disable auto_recovery when:**
- Testing new methods
- Need strict convergence criteria
- Publication-quality results requiring manual verification

---

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
To use DELFIN, you must be authorized to use ORCA 6.1.1. You can download the latest version of ORCA here:
https://orcaforum.kofo.mpg.de/app.php/portal

***IMPORTANT: ORCA 6.1.1 requires a valid license and registration. Academic users can obtain free access, but commercial use requires a commercial license. Please carefully review and comply with ORCA's license terms before use.***
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

---

## 📚 Appendix (For Developers)

### Project Layout

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

### Development Notes

* Update CLI entry point via `pyproject.toml`: `"[project.scripts] delfin = \"delfin.cli:main\""`
* Build a wheel: `pip wheel .` (inside `delfin/`)
* Run tests/workflow locally using a fresh virtual environment to catch missing deps
* Install development tools: `pip install -e '.[dev]'`
* Format code: `black .` / Lint code: `ruff check .`

---

## License

This project is licensed under the GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).

You should have received a copy of the GNU Lesser General Public License along with this repository in the files `COPYING` and `COPYING.LESSER`.  
If not, see <https://www.gnu.org/licenses/>.

Non-binding citation request:  
If you use this software in research, please cite the associated paper (see [CITATION.cff](./CITATION.cff)).


  
