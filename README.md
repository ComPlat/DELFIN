<h1>
  <img src="delfin/logo/DELFIN_logo.png" alt="DELFIN logo" width="44" align="center">
  DELFIN
</h1>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17208145.svg)](https://doi.org/10.5281/zenodo.17208145)
[![PyPI version](https://img.shields.io/pypi/v/delfin-complat.svg)](https://pypi.org/project/delfin-complat/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/delfin-complat.svg)](https://pypistats.org/packages/delfin-complat)

> 📄 **Preprint**: *Hartmann, M. et al. "DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials"*, ChemRxiv (2025). https://doi.org/10.26434/chemrxiv-2025-4c256

**DELFIN** is a modular computational chemistry platform that automates quantum chemical workflows — from SMILES input to realistic property predictions. It combines DFT, semi-empirical methods, ML potentials, and AI tools behind a unified interface with an interactive browser-based dashboard.

### What DELFIN can do

| Capability | Description |
|------------|-------------|
| **Redox Potentials** | Automated spin-state prediction and redox potential calculation via OCCUPIER/classic workflows |
| **Excited-State Dynamics** | ISC/RISC rates, fluorescence, phosphorescence, SOC coupling, E₀₀ energies, ΔE(S-T) gaps |
| **TADF Screening** | xTB-based singlet-triplet gap estimation for OLED material discovery |
| **Spectroscopy** | UV-Vis absorption, IR vibrational spectra, AFP (absorption/fluorescence/phosphorescence) plots |
| **Hyperpolarizability** | Static and frequency-dependent β tensors for NLO materials |
| **Structure Generation** | SMILES→3D for organic molecules and metal complexes (RDKit, Open Babel, architector, PSO builder) |
| **Conformer Sampling** | GUPPY multi-start sampling, CREST conformer search, XTB-GOAT global optimization |
| **ML Potentials** | 8 backends (ANI-2x, MACE, CHGNet, M3GNet, ...) for fast energy/force evaluation |
| **Crystal Structure Prediction** | Genarris integration for random crystal generation |
| **CO2 Coordination** | Automated CO2 placement, distance/rotation scans |
| **Reporting** | DOCX reports with embedded spectra, JSON export, text summaries |

---

## 🚀 Installation

**PyPI Package**: https://pypi.org/project/delfin-complat/

### Requirements
- **Python 3.10 or 3.11**
- **ORCA 6.1.1** in your `PATH` (`orca` and `orca_pltvib`) — [free for academic use](https://orcaforum.kofo.mpg.de/app.php/portal)
- **Optional:** `crest` and `xtb` (for CREST/xTB workflows)
- **Optional:** `xtb4stda`, `stda`, and `std2` plus the required `xtb4stda` runtime files (for xTB-based response/screening workflows)
- **Optional:** Any of the 88 supported computational tools — auto-detected via PATH, installable via Dashboard
- **Optional (Dashboard):** JupyterLab/Notebook or Voila for interactive UI usage

### Install Methods

**Standard install (recommended for most users):**
```bash
pip install delfin-complat
```

**Development install (from source):**
```bash
git clone https://github.com/ComPlat/DELFIN.git
cd DELFIN
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

All Python dependencies (e.g., RDKit/OpenBabel for SMILES workflows, ipywidgets/py3Dmol for dashboard visualisation) are installed automatically.

This exposes the console command **`delfin`** and enables `python -m delfin`.

External QM binaries are not installed automatically by `pip`. For the local binary-based setup of `xtb`, `crest`, `xtb4stda`, `stda`, `std2`, and the `xtb4stda` runtime bundle, see `delfin/qm_tools/README.txt`.

### External QM Tool Setup

After installing the Python package, initialise the bundled QM tool wrapper and validate the external binaries:

```bash
source delfin/qm_tools/env.sh
USE_SYSTEM_TOOLS=1 bash delfin/qm_tools/install_qm_tools.sh
bash delfin/qm_tools/check_qm_tools.sh
```

---

## Quick start

Create a working directory with at least these two files:

* `CONTROL.txt` — your control/config file
* `input.txt` — either
  * XYZ body (without the first two header lines), or
  * a single-line SMILES string
* starting from a full `.xyz` file is optional (`delfin --define=input.xyz` auto-converts)

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

## 🧪 Dashboard

Start in Jupyter/Voila and create the dashboard from Python:

```python
from delfin.dashboard import create_dashboard
ctx = create_dashboard(backend="auto")
```

`backend="auto"` selects SLURM if available, otherwise local execution.

### Dashboard Tabs

| Tab | Purpose |
|-----|---------|
| **Submit Job** | SMILES/XYZ input, 3D preview, `CONVERT SMILES`, `QUICK CONVERT`, `CONVERT + UFF`, `BUILD COMPLEX`, `ARCHITECTOR`, `SUBMIT GUPPY` |
| **Recalc** | Edit and resubmit existing CONTROL.txt |
| **ORCA Builder** | Interactive ORCA input generation with geometry preview |
| **TURBOMOLE Builder** | Turbomole define workflow (SLURM backends) |
| **Job Status** | Real-time queue monitoring (local/SLURM), resource usage, job cancellation |
| **Calculations** | File browser, search, recalculation trigger, energy statistics |
| **Archive** | Archive browser with statistics |
| **Settings** | Tool detection, per-tool Install/Update buttons, runtime configuration |

### GUI with Voila

```bash
delfin-voila                 # starts on port 8866
delfin-voila --port 9000     # custom port
delfin-voila --dark          # dark theme
```

Detailed documentation: [docs/SETTINGS_AND_SETUP.md](docs/SETTINGS_AND_SETUP.md)

---

## ✨ Features

### Redox Potential Prediction (OCCUPIER)

DELFIN's core workflow automates spin-state identification and redox potential calculation:

- **Multi-step redox**: Up to 3 sequential oxidation/reduction steps
- **Parallel workflows**: Oxidation and reduction run simultaneously with automatic PAL splitting
- **Smart recalc**: Fingerprint-based skip logic avoids unnecessary reruns

### Excited-State Dynamics (ESD)

The ESD module calculates photophysical properties:

- **Electronic states**: S0, S1, S2, T1, T2 geometry optimization
- **ISC / RISC rates**: Intersystem and reverse intersystem crossing - **Internal conversion (IC)**: Non-radiative decay rates
- **Fluorescence & phosphorescence rates**: Radiative lifetimes including per-sublevel phosphorescence
- **E₀₀ energies**: Adiabatic 0-0 transition energies
- **ΔE(S-T)**: Singlet-triplet energy gaps
```ini
ESD_modul=yes
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1
ICs=S1>S0
```

### TADF Screening (xTB-based)

Fast screening of TADF candidates using semi-empirical methods:

- S0/T1 optimization via xTB, S1 estimation via Stokes shift
- ΔE(S-T) gap calculation
- First allowed singlet state detection
- TADF efficiency metrics

### Spectroscopy & Properties

| Module | Output |
|--------|--------|
| **UV-Vis** | TD-DFT absorption spectra, oscillator strengths, transition analysis |
| **IR** | Vibrational frequencies, intensities, Lorentzian broadening, transmittance |
| **AFP** | Combined absorption/fluorescence/phosphorescence spectrum plot |
| **Hyperpolarizability** | Static & frequency-dependent β tensors (NLO), dipole moments |

### SMILES→3D Structure Generation

DELFIN provides multiple conversion methods for organic and metal-containing systems:

| Button | Method | Best for |
|--------|--------|----------|
| `CONVERT SMILES` | Full isomer/conformer search (RDKit + Open Babel) | Thorough exploration |
| `QUICK CONVERT` | Fast single-conformer generation | Quick previews |
| `CONVERT + UFF` | Full search + UFF force-field refinement | Refined geometries |
| `BUILD COMPLEX` | ORCA/XTB DOCKER stepwise complex assembly | Metal complexes (job submission) |
| `ARCHITECTOR` | [architector](https://github.com/lanl/Architector) automated 3D generation | Metal complexes (instant preview) |
| `SUBMIT GUPPY` | Multi-start XTB sampling with ranked trajectories | Robust start structures |

For coordination complexes, DELFIN combines Open Babel conformer pools, RDKit multi-seed embedding, topological isomer enumeration, and fragment sanity checks.

### ML Potentials & Unified Calculator Factory

34 computational backends accessible through a single interface returning standard ASE Calculator objects:

```python
from delfin.calculators import create_calculator

calc = create_calculator("ani2x", device="cuda")       # ML potential
calc = create_calculator("orca", method="B3LYP")       # DFT
calc = create_calculator("xtb")                        # Semi-empirical
calc = create_calculator("vasp", xc="PBE", kpts=[4,4,4])  # Periodic DFT

atoms.calc = calc
energy = atoms.get_potential_energy()
```

Automatic CUDA validation with CPU fallback. All backends are lazily loaded — only imported when actually used.

### AI/ML Tools Integration

21 AI/ML tools across molecular generation, retrosynthesis, screening, and metal complex design — all behind runtime validation with per-tool Install/Update buttons in the Dashboard Settings tab.

### Auto-Detection of External Programs

DELFIN automatically detects 88 computational chemistry programs via PATH search. This works seamlessly with cluster module systems (`module load gaussian/16` → DELFIN detects it). Programs that can't be pip-installed (Gaussian, VASP, Turbomole, ...) are detected and reported but not installed.

<details>
<summary><b>Full list of supported tools (88 total)</b></summary>

**ML Potentials (8):** ANI-2x, AIMNet2, MACE-OFF, CHGNet, M3GNet/MatGL, SchNetPack, NequIP/Allegro, ALIGNN

**QM Programs — Ab initio / DFT (11):** ORCA, Gaussian (g16/g09), Turbomole, NWChem, Q-Chem, GAMESS, Molpro, Dalton, Psi4, CFOUR, MRCC

**QM Programs — Periodic / Solid State (11):** VASP, Quantum ESPRESSO, CP2K, FHI-aims, CRYSTAL, SIESTA, GPAW, FLEUR, WIEN2k, Elk, ABINIT

**QM Programs — Multireference (3):** OpenMolcas, BAGEL, Columbus

**Semi-empirical (3):** xTB, MOPAC, Sparrow

**MD Engines (5):** GROMACS, LAMMPS, AMBER, NAMD, OpenMM

**AI/ML — Foundation Models (3):** MoLFormer, Uni-Mol, ChemBERTa

**AI/ML — Generative (2):** REINVENT4, SyntheMol

**AI/ML — Conformers (2):** GeoMol, torsional-diffusion

**AI/ML — Crystal Generation (2):** MatterGen, CDVAE

**AI/ML — Retrosynthesis (3):** AiZynthFinder, RXNMapper, LocalRetro

**AI/ML — Screening / ADMET (2):** DeepChem, ADMETlab

**AI/ML — Metal Complex ML (2):** molSimplify, architector

**Analysis / Post-Processing (12):** cclib, Multiwfn, CENSO, morfeus, nglview, Packmol, NBO, AIMAll, critic2, Chargemol, LOBSTER, phonopy

**Wrapper Libraries (4):** pymatgen, QCEngine, MDAnalysis, pymolpro

**Visualization (6):** plotly, VMD, Avogadro, Jmol, ChimeraX, IQmol

**Python-Only Backends (3):** PySCF, OpenMM, PLAMS

</details>

### Conformer Search & Sampling

- **CREST**: Conformer-rotamer ensemble generation and sorting
- **XTB-GOAT**: Gradient-based global optimization
- **GUPPY**: Multi-start SMILES sampling with parallel XTB optimization and energy-ranked trajectories

### Crystal Structure Prediction (CSP)

Genarris integration for random crystal structure generation with configurable space groups, Z values, and MPI-parallel execution.

### CO2 Coordinator

Automated CO2 placement around metal centers with relaxed distance scans (1.6–4.0 Å) and rotation scans (±180°).

### Reporting & Export

| Format | Command | Content |
|--------|---------|---------|
| **Text** | automatic | `DELFIN.txt` (redox potentials), `OCCUPIER.txt` (occupation tracking), `ESD_report.txt` |
| **DOCX** | `delfin --report docx` | Word document with embedded spectra plots and structured tables |
| **JSON** | `delfin --json` | Machine-readable `DELFIN_Data.json` with all energies, rates, spectra |
| **UV-Vis** | `delfin_ESD output.out` | Parsed UV-Vis spectrum with transitions and oscillator strengths |
| **IR** | `delfin_IR output.out` | Parsed IR spectrum with vibrational modes and intensities |
| **AFP** | `delfin --afp` | Absorption/fluorescence/phosphorescence combined plot |

---

## 📋 CLI Reference

### Basic Commands

- `delfin` — Run DELFIN workflow in current directory (requires `CONTROL.txt` and `input.txt`)
- `delfin /path/to/project` — Run DELFIN workflow in specified directory
- `python -m delfin` — Alternative way to run DELFIN

### Companion CLI Tools

- `delfin-build [input.txt]` — Build metal complexes stepwise from SMILES using ORCA/XTB Docker workflow
- `delfin-guppy [input.txt]` — Multi-start SMILES sampling workflow with repeated XTB optimization and ranked trajectories
- `delfin-voila` — Launch the DELFIN Dashboard as a standalone web app via Voila
- `delfin-json` — Collect DELFIN project outputs into JSON
- `delfin_ESD` — UV-Vis spectrum report from ORCA ESD output
- `delfin_IR` — IR spectrum report from ORCA frequency output

### Setup & Configuration

- `delfin --define[=input.xyz] [--overwrite]` — Create/update `CONTROL.txt` and optionally convert XYZ to `input.txt`
- `delfin /path/to/project --define[=input.xyz] [--overwrite]` — Same as above, but in a different workspace directory
- `delfin --control /path/to/CONTROL.txt` — Run workflow from another directory while normalising all paths

### Execution Control

- `delfin --recalc` — Re-parse existing results and only restart missing or incomplete jobs
- `delfin WORKSPACE --recalc --occupier-override STAGE=INDEX` — Force a specific OCCUPIER index for a stage during recalc
- `delfin --report` — Re-calculate redox potentials from existing output files without launching new calculations
- `delfin --imag` — Eliminate imaginary modes from existing ORCA results (`*.out`/`*.hess`) and regenerate the summary report
- `delfin stop --workspace PATH` — Send a graceful stop signal to running DELFIN processes associated with a workspace

### Cleanup & Maintenance

- `delfin --no-cleanup` — Keep temporary files and scratch folders after the pipeline finishes
- `delfin --cleanup` — Remove previously generated intermediates and exit immediately
- `delfin cleanup [--dry-run] [--workspace PATH] [--scratch PATH]` — Fine control over workspace/scratch cleanup
- `delfin cleanup --orca` — Stop running ORCA jobs, purge OCCUPIER scratch folders
- `delfin --purge` — Remove DELFIN-generated artifacts after confirmation

### Export & Reporting

- `delfin --json` — Collect key results into `DELFIN_Data.json`
- `delfin --report docx` — Generate a Word report in DOCX format
- `delfin --afp` — Generate an AFP spectrum plot from existing ESD results

### Specialized Workflows

- `delfin co2 ...` — CO2 Coordinator helper workflow
- `delfin ESD` — Run excited-state dynamics workflow (requires `ESD_modul=yes` in CONTROL)
- `delfin-guppy --runs N --parallel-jobs M --pal P` — Broadened start-structure sampling
- `delfin-build --goat [--no-ligand-goat]` — Ligand docking with optional GOAT optimization

---

## ⚙️ Configuration (CONTROL.txt)

DELFIN is configured via `CONTROL.txt` in your working directory. Key settings include:

### Workflow Control
* `method = OCCUPIER | classic | manually` (leave empty for ESD-only runs)
* `OCCUPIER_method = auto | manually` (auto uses adaptive tree-based sequences)
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

---

## 🏗️ Architecture

### Global Resource Management

DELFIN uses a **global job manager singleton** to coordinate all computational workflows:

* **Single source of truth:** PAL is read once from `CONTROL.txt` and managed centrally
* **Automatic PAL splitting:** Parallel ox/red workflows share cores (PAL=12 → 6+6)
* **Thread-safe execution:** Shared resource pool prevents race conditions
* **Subprocess coordination:** OCCUPIER subprocesses inherit limits via environment variables

### Modular Tool Architecture

```
delfin/
  calculators.py       ← Unified ASE calculator factory (34 backends)
  mlp_tools/           ← ML potential backends with lazy loading
  ai_tools/            ← AI/ML tools registry and per-tool installers
  analysis_tools/      ← Analysis wrappers (cclib, Packmol, Multiwfn, CENSO, morfeus)
  csp_tools/           ← Crystal structure prediction (Genarris)
  runtime_setup.py     ← Auto-detection of 88 external programs
  dashboard/           ← 10-tab interactive dashboard (Voila/JupyterLab)
```

All tool integrations follow the same pattern:
1. **Lazy loading** — `importlib.util.find_spec()` for availability checks, no imports until use
2. **Per-tool install** — Individual Install/Update buttons in Dashboard Settings
3. **Auto-detection** — `shutil.which()` for binaries, compatible with cluster module systems
4. **Unified interface** — All backends return standard ASE Calculator objects

### Cluster & HPC Integration

* **Backends:** Local execution, SLURM, SSH-based remote transfer
* **Auto-resource detection:** CPUs and memory on SLURM/PBS/LSF clusters
* **Scratch directory:** `DELFIN_SCRATCH=/path/to/scratch`
* **Logging:** `delfin_run.log` per workspace, `occupier.log` per subprocess
* **Programmatic API:** `delfin.api.run(control_file="CONTROL.txt")` for notebooks and workflow engines

**Cluster templates:** see `examples/` for SLURM, PBS, and LSF submit scripts.

---

## 🔧 Automatic Error Recovery & Retry System

DELFIN includes an intelligent error recovery system that automatically detects and fixes common ORCA calculation failures.

### Quick Enable

```ini
enable_auto_recovery=yes
max_recovery_attempts=1
```

### How It Works

```
ORCA fails → Detect error type → Modify input with fixes → Continue from last .gbw and xyz → Retry
```

| Error | Automatic Fix |
|-------|--------------|
| **SCF not converged** | SlowConv → VerySlowConv + KDIIS |
| **TRAH segfault** | NoAutoTRAH |
| **Geometry not converged** | Smaller trust radius → Loose criteria |
| **MPI crashes** | Reduce cores |
| **Memory errors** | Reduce maxcore and PAL |
| **Transient system errors** | Exponential backoff retry |

See **[docs/RETRY_LOGIC.md](docs/RETRY_LOGIC.md)** for the complete guide.

---

## Troubleshooting

* **`CONTROL.txt` not found** — Create it via `delfin --define` (or copy your own).
* **Input file not found** — Run `delfin --define=your.xyz` to auto-convert.
* **ORCA not found** — Ensure `orca` is in your PATH: `which orca`.
* **CREST/xTB tools missing** — Disable corresponding flags in `CONTROL.txt` or install and add to PATH.
* **Optional tool not detected** — Check Dashboard Settings tab or run `delfin --diagnostics`.

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

---

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

- Hartmann, M. (2026). *DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials* (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.17208145
- Hartmann, M. (2025). *DELFIN: Automated prediction of preferred spin states and redox potentials*. ChemRxiv. https://chemrxiv.org/engage/chemrxiv/article-details/68fa0e233e6156d3be78797a

### BibTeX
```bibtex
@software{hartmann2025delfin,
  author  = {Hartmann, Maximilian},
  title   = {DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials},
  version = {v1.1.0},
  year    = {2026},
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
  __main__.py              # enables `python -m delfin`
  cli.py                   # main CLI entry point
  pipeline.py              # high-level orchestration (classic/manually/OCCUPIER)
  config.py                # CONTROL.txt parsing & helpers
  calculators.py           # unified ASE calculator factory (34 backends)
  api.py                   # programmatic API for notebooks/workflows

  # ── Core Workflows ──
  occupier.py              # OCCUPIER workflow (sequence execution + summary)
  occupier_auto.py         # auto OCCUPIER sequence management and tree navigation
  deep_auto_tree.py        # adaptive BS evolution tree
  esd_module.py            # excited-state dynamics (ISC/IC/fluorescence/phosphorescence)
  esd_input_generator.py   # ORCA input builders for ESD states
  tadf_xtb.py              # TADF screening via xTB
  hyperpol.py              # hyperpolarizability (NLO)
  xtb_crest.py             # xTB / GOAT / CREST / ALPB solvation workflows
  imag.py                  # imaginary frequency elimination
  guppy_sampling.py        # multi-start SMILES sampling + ranking
  build_up_complex.py      # stepwise metal-complex assembly (ORCA/XTB DOCKER)
  build_up_complex2.py     # PSO-based metal complex builder

  # ── Resource Management ──
  global_manager.py        # singleton global job manager
  global_scheduler.py      # dynamic shared resource scheduling
  dynamic_pool.py          # dynamic core pool for job scheduling
  cluster_utils.py         # SLURM/PBS/LSF resource detection
  runtime_setup.py         # auto-detection of 88 external programs

  # ── Tool Integrations ──
  mlp_tools/               # ML potentials (ANI-2x, MACE, CHGNet, M3GNet, ...)
  ai_tools/                # AI/ML tools (21 tools across 9 categories)
  analysis_tools/          # analysis wrappers (cclib, Packmol, Multiwfn, CENSO, morfeus)
  csp_tools/               # crystal structure prediction (Genarris)
  qm_tools/                # external QM binary management

  # ── Dashboard ──
  dashboard/               # 10-tab interactive dashboard (Voila/JupyterLab)

  # ── Reporting ──
  reporting/               # DOCX, JSON, text report generation
  co2/                     # CO2 coordinator workflows

  # ── Shared ──
  common/                  # logging, paths, ORCA block assembly
```

### Development Notes

* CLI entry points are defined in `pyproject.toml`
* Build a wheel: `pip wheel .`
* Install development tools: `pip install -e '.[dev]'`
* Format code: `black .` / Lint code: `ruff check .`

---

## License

This project is licensed under the GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).

You should have received a copy of the GNU Lesser General Public License along with this repository in the files `COPYING` and `COPYING.LESSER`.
If not, see <https://www.gnu.org/licenses/>.

Non-binding citation request:
If you use this software in research, please cite the associated paper (see [CITATION.cff](./CITATION.cff)).
