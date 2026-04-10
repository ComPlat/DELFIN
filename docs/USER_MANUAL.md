# DELFIN User Manual

This manual is a practical guide for using DELFIN. It covers installation, workspace setup, all CLI tools, the full `CONTROL.txt` reference, dashboard usage, and troubleshooting.

For scientific methodology and validation details, see [methodology.md](methodology.md).
For implementation and architecture, see the [README](../README.md).

---

## Table of Contents

1. [Installation](#1-installation)
2. [Getting Started](#2-getting-started)
3. [Input Files](#3-input-files)
4. [CONTROL.txt Reference](#4-controltxt-reference)
5. [CLI Reference](#5-cli-reference)
6. [Workflow Modes](#6-workflow-modes)
7. [Optional Modules](#7-optional-modules)
8. [Structure Generation & Sampling](#8-structure-generation--sampling)
9. [Dashboard](#9-dashboard)
10. [Settings & Runtime Configuration](#10-settings--runtime-configuration)
11. [Error Recovery & Retry System](#11-error-recovery--retry-system)
12. [Reporting & Export](#12-reporting--export)
13. [Cluster & HPC Usage](#13-cluster--hpc-usage)
14. [Troubleshooting](#14-troubleshooting)
15. [Recipes & Examples](#15-recipes--examples)

---

## 1. Installation

### Requirements

- **Python 3.10 or 3.11**
- **ORCA 6.1.1** in your `PATH` — [free for academic use](https://orcaforum.kofo.mpg.de/app.php/portal)
- **Optional:** `xtb`, `crest` (for xTB/CREST workflows)
- **Optional:** `xtb4stda`, `stda`, `std2` (for xTB-based screening)
- **Optional:** `censo`, `anmr`, `c2anmr`, `nmrplot` (for ensemble NMR)
- **Optional:** JupyterLab/Notebook or Voila (for dashboard)

### Standard install

```bash
pip install delfin-complat
```

### Development install (from source)

```bash
git clone https://github.com/ComPlat/DELFIN.git
cd DELFIN
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Optional extras (source install)

```bash
pip install -e ".[agent,analysis,mlp]"
```

### External QM tool setup

After installing the Python package, set up the bundled QM tool wrappers:

```bash
source delfin/qm_tools/env.sh
USE_SYSTEM_TOOLS=1 bash delfin/qm_tools/install_qm_tools.sh
bash delfin/qm_tools/check_qm_tools.sh
```

This can also be done from the dashboard Settings tab (see [Section 10](#10-settings--runtime-configuration)).

---

## 2. Getting Started

DELFIN works with a **workspace directory** containing at least two files:

- `CONTROL.txt` — workflow and runtime configuration
- `input.txt` — geometry (XYZ body without header lines) or a SMILES string

### Create a workspace

```bash
mkdir my_project && cd my_project
delfin --define
```

This generates a template `CONTROL.txt` and an empty `input.txt`.

### Start from an XYZ file

```bash
delfin --define=structure.xyz
```

DELFIN removes the first two XYZ header lines, writes `input.txt`, and sets `input_file=input.txt` in `CONTROL.txt`.

### Create a workspace in another directory

```bash
delfin /path/to/project --define
```

### Edit CONTROL.txt

Set at minimum:

```ini
charge=0
calc_initial=yes
oxidation_steps=1
reduction_steps=1
method=OCCUPIER
```

### Run

```bash
delfin
```

or equivalently:

```bash
python -m delfin
```

or for a different workspace:

```bash
delfin /path/to/project
```

### Outputs

Typical results include:

| File | Content |
|------|---------|
| `DELFIN.txt` | Text summary with redox potentials |
| `DELFIN_Data.json` | Structured data export |
| `OCCUPIER.txt` | OCCUPIER workflow tracking |
| `delfin_run.log` | Run-level log |
| `initial_OCCUPIER/` | Initial state calculations |
| `ox_step_1_OCCUPIER/` | Oxidation step directories |
| `red_step_1_OCCUPIER/` | Reduction step directories |

---

## 3. Input Files

### input.txt

`input.txt` contains either:

**Atomic coordinates** (XYZ body without the first two header lines):

```
C     0.000000     0.000000     0.000000
H     0.000000     0.000000     1.089000
H     1.026719     0.000000    -0.363000
...
```

**A SMILES string** (single line):

```
[Fe+2]([N]1=CC=CC=1)([N]2=CC=CC=2)([N]3=CC=CC=3)([N]4=CC=CC=4)([N]5=CC=CC=5)[N]6=CC=CC=6
```

**QM/MM region splitting**: Insert a line containing only `$` to separate QM and MM regions. DELFIN preserves this splitting for all follow-up jobs automatically.

### CONTROL.txt

`CONTROL.txt` is the main configuration file. All keys use `key=value` format. Lines starting with `----` are section separators (ignored by the parser). Empty values are allowed and fall back to defaults.

See [Section 4](#4-controltxt-reference) for the complete reference.

---

## 4. CONTROL.txt Reference

### Input & Identity

| Key | Default | Description |
|-----|---------|-------------|
| `input_file` | `input.txt` | Geometry or SMILES input file |
| `NAME` | (empty) | Project name |
| `SMILES` | (empty) | Optional SMILES string (alternative to input file) |
| `charge` | (required) | System charge |

### Solvation

| Key | Default | Description |
|-----|---------|-------------|
| `implicit_solvation_model` | `CPCM` | Solvation model (`CPCM` or `SMD`) |
| `solvent` | (required) | Solvent name (e.g., `acetonitrile`, `dmf`, `dcm`, `thf`, `dmso`, `acetone`) |
| `XTB_SOLVATOR` | `no` | Enable xTB ALPB solvation |
| `number_explicit_solv_molecules` | `2` | Number of explicit solvent molecules |

### Global Geometry Optimisation

| Key | Default | Description |
|-----|---------|-------------|
| `xTB_method` | `XTB2` | xTB method for pre-optimisation |
| `smiles_converter` | (empty) | SMILES conversion method: `QUICK`, `NORMAL`, `GUPPY`, or `ARCHITECTOR` |
| `XTB_OPT` | `no` | Run xTB geometry optimisation before DFT |
| `XTB_GOAT` | `no` | Run XTB-GOAT global optimisation |
| `CREST` | `no` | Run CREST conformer search |
| `multiplicity_global_opt` | (empty) | Override multiplicity for pre-optimisation |

### Imaginary Frequency Elimination (IMAG)

| Key | Default | Description |
|-----|---------|-------------|
| `IMAG` | `yes` | Enable imaginary frequency elimination |
| `IMAG_scope` | `initial` | Scope: `initial` (initial state only) or `all` (all states) |
| `IMAG_option` | `2` | IMAG method option |
| `allow_imaginary_freq` | `0` | Number of allowed imaginary frequencies |
| `IMAG_sp_energy_window` | `1e-3` | Energy window (Hartree) for single-point comparison |
| `IMAG_optimize_candidates` | `no` | Optimize IMAG candidate structures |

### Properties of Interest

| Key | Default | Description |
|-----|---------|-------------|
| `calc_prop_of_interest` | `no` | Calculate additional properties |
| `properties_of_interest` | `IP,EA` | Properties to calculate |
| `reorganisation_energy` | `lambda_p,lambda_m` | Reorganisation energy types |

### Redox Workflow

| Key | Default | Description |
|-----|---------|-------------|
| `calc_initial` | `yes` | Calculate initial (neutral) state |
| `oxidation_steps` | (empty) | Number of oxidation steps (1–3) |
| `reduction_steps` | (empty) | Number of reduction steps (1–3) |
| `method` | (required) | Workflow method: `classic`, `manually`, or `OCCUPIER` |
| `calc_potential_method` | `2` | Potential calculation method |

### ESD Module (Excited-State Dynamics)

| Key | Default | Description |
|-----|---------|-------------|
| `ESD_modul` | `no` | Enable excited-state dynamics |
| `ESD_modus` | `TDDFT` | Method: `TDDFT`, `deltaSCF`, or `hybrid1` |
| `ESD_T1_opt` | (empty) | T1 optimisation method: `uks` or `tddft` |
| `ESD_frequency` | `yes` | Run frequency calculation for ESD states |
| `states` | `S1,T1,S2,T2` | Electronic states to compute |
| `ISCs` | `S1>T1,T1>S1` | Intersystem crossing rates |
| `ICs` | `S2>S1` | Internal conversion rates |
| `emission_rates` | `f,p` | Emission rates: `f` (fluorescence), `p` (phosphorescence) |
| `phosp_IROOT` | `1,2,3` | Phosphorescence IROOT sublevels |
| `phosp_keywords` | (empty) | Additional phosphorescence keywords |
| `fluor_keywords` | (empty) | Additional fluorescence keywords |
| `TROOTSSL` | `-1,0,1` | TROOT spin sublevels |
| `addition_S0` | (empty) | Additional S0 ORCA keywords |
| `DOHT` | `TRUE` | Duschinsky/Herzberg-Teller coupling |
| `ESD_LINES` | `LORENTZ` | Broadening function |
| `ESD_LINEW` | `50` | Line width |
| `ESD_INLINEW` | `250` | Input line width |
| `ESD_NPOINTS` | `131072` | Number of FFT points |
| `ESD_MAXTIME` | `12000` | Maximum time for ESD calculation |
| `hybrid1_geom_MaxIter` | `60` | Max geometry iterations for hybrid1 |

### xTB Hyperpolarizability (sTD-DFT-xTB)

| Key | Default | Description |
|-----|---------|-------------|
| `hyperpol_xTB` | `no` | Enable xTB hyperpolarizability |
| `hyperpol_xTB_xyz` | `start.txt` | Input geometry |
| `hyperpol_xTB_preopt` | `none` | Pre-optimisation method |
| `hyperpol_xTB_engine` | `std2` | Calculation engine |
| `hyperpol_xTB_bfw` | `no` | Bandwidth-filtered weights |
| `hyperpol_xTB_wavelengths` | (empty) | Wavelengths for frequency-dependent calculations |
| `hyperpol_xTB_energy_window` | `15.0` | Energy window (eV) |

### xTB TADF Screening

| Key | Default | Description |
|-----|---------|-------------|
| `tadf_xTB` | `no` | Enable TADF screening |
| `tadf_xTB_xyz` | `start.txt` | Input geometry |
| `tadf_xTB_preopt` | `none` | Pre-optimisation method |
| `tadf_xTB_excited_method` | `stda` | Excited-state method |
| `tadf_xTB_bfw` | `no` | Bandwidth-filtered weights |
| `tadf_xTB_energy_window` | `10.0` | Energy window (eV) |
| `tadf_xTB_run_t1_opt` | `yes` | Optimize T1 state |

### Thermodynamics

| Key | Default | Description |
|-----|---------|-------------|
| `thermodynamics` | `no` | Enable thermodynamics workflow |
| `thermodynamics_mode` | (empty) | Mode: `auto` or `reaction` |
| `thermodynamics_reaction` | (empty) | Reaction SMILES: `a*{SMILES}+b*{SMILES}...>>>c*{SMILES}+d*{SMILES}...` |
| `n_explicit_solvent` | `6` | Number of explicit solvent molecules |
| `logK_exp` | (empty) | Experimental log K for comparison |
| `thdy_smiles_converter` | (empty) | Converter: `QUICK`, `NORMAL`, `GUPPY`, or `ARCHITECTOR` |
| `thdy_preopt` | (empty) | Pre-optimisation: `none`, `xtb`, `crest`, or `goat` |

### Electrical Properties

| Key | Default | Description |
|-----|---------|-------------|
| `elprop_Dipole` | `no` | Calculate dipole moment |
| `elprop_Quadrupole` | `no` | Calculate quadrupole moment |
| `elprop_Hyperpol` | `no` | Calculate hyperpolarizability |
| `elprop_Polar` | `no` | Calculate polarizability |
| `elprop_PolarVelocity` | `no` | Velocity-gauge polarizability |
| `elprop_PolarDipQuad` | `no` | Dipole-quadrupole polarizability |
| `elprop_PolarQuadQuad` | `no` | Quadrupole-quadrupole polarizability |

### deltaSCF Settings

| Key | Default | Description |
|-----|---------|-------------|
| `deltaSCF_DOMOM` | `true` | Use MOM (Maximum Overlap Method) |
| `deltaSCF_PMOM` | `false` | Use PMOM |
| `deltaSCF_keepinitialref` | `true` | Keep initial reference |
| `deltaSCF_SOSCFHESSUP` | `LSR1` | SOSCF Hessian update method |
| `deltaSCF_keywords` | `FreezeAndRelease` | Additional deltaSCF keywords |
| `deltaSCF_maxiter` | `300` | Max SCF iterations |
| `deltaSCF_SOSCFConvFactor` | `500` | SOSCF convergence factor |
| `deltaSCF_SOSCFMaxStep` | `0.1` | SOSCF maximum step size |

### TD-DFT Settings

| Key | Default | Description |
|-----|---------|-------------|
| `TDDFT_TDDFT_maxiter` | `500` | Max TD-DFT iterations |
| `TDDFT_nroots` | `15` | Number of roots |
| `TDDFT_maxdim` | `30` | Max Davidson dimension |
| `TDDFT_TDA` | `FALSE` | Tamm-Dancoff approximation |
| `TDDFT_followiroot` | `true` | Follow IROOT during optimisation |
| `TDDFT_SOC` | `false` | Spin-orbit coupling |

### Manual Multiplicity Settings (for `method=manually`)

| Key | Description |
|-----|-------------|
| `multiplicity_0` | Multiplicity for initial state |
| `BrokenSym_0` | Broken symmetry specification for initial state |
| `multiplicity_ox1` / `BrokenSym_ox1` | First oxidation step |
| `multiplicity_ox2` / `BrokenSym_ox2` | Second oxidation step |
| `multiplicity_ox3` / `BrokenSym_ox3` | Third oxidation step |
| `multiplicity_red1` / `BrokenSym_red1` | First reduction step |
| `multiplicity_red2` / `BrokenSym_red2` | Second reduction step |
| `multiplicity_red3` / `BrokenSym_red3` | Third reduction step |

### Level of Theory

| Key | Default | Description |
|-----|---------|-------------|
| `functional` | `PBE0` | Exchange-correlation functional |
| `disp_corr` | `D4` | Dispersion correction |
| `ri_jkx` | `RIJCOSX` | RI approximation |
| `relativity` | `ZORA` | Relativistic method (`ZORA`, `X2C`, `DKH`) |
| `aux_jk` | `def2/J` | Auxiliary basis (non-relativistic) |
| `aux_jk_rel` | `SARC/J` | Auxiliary basis (relativistic) |
| `main_basisset` | `def2-SVP` | Main basis set |
| `main_basisset_rel` | `ZORA-def2-SVP` | Main basis set (relativistic) |
| `metal_basisset` | `def2-TZVP` | Metal basis set |
| `metal_basisset_rel` | `SARC-ZORA-TZVP` | Metal basis set (relativistic) |
| `first_coordination_sphere_metal_basisset` | `no` | Use metal basis set for first coordination sphere |
| `first_coordination_sphere_scale` | `1.3` | Bond detection scale factor |
| `geom_opt` | `OPT` | Geometry optimisation keyword |
| `freq_type` | `FREQ` | Frequency keyword |
| `initial_guess` | `PModel` | SCF initial guess |
| `temperature` | `298.15` | Temperature (K) |
| `maxiter` | `125` | Max SCF iterations |
| `qmmm_option` | `QM/PBEH-3c` | QM/MM level for QM/MM regions |

DELFIN automatically selects relativistic settings when 4d/5d transition metals are detected.

### Reference Values

| Key | Default | Description |
|-----|---------|-------------|
| `E_ref` | (auto) | Reference electrode potential override (V vs SHE) |
| `Literature_reference` | (empty) | Literature reference |
| `reference_CV` | `V Vs. Fc+/Fc` | Reference convention |
| `E_00_exp` | (empty) | Experimental E₀₀ for comparison |
| `E_red_exp` / `E_red_2_exp` / `E_red_3_exp` | (empty) | Experimental reduction potentials |
| `E_ox_exp` / `E_ox_2_exp` / `E_ox_3_exp` | (empty) | Experimental oxidation potentials |
| `*E_red_exp` / `*E_ox_exp` | (empty) | Excited-state experimental potentials |

### Prints

| Key | Default | Description |
|-----|---------|-------------|
| `print_MOs` | `no` | Print molecular orbital information |
| `print_Loewdin_population_analysis` | `no` | Print Loewdin population analysis |

### Resource Settings

| Key | Default | Description |
|-----|---------|-------------|
| `PAL` | `12` | Number of CPU cores for ORCA |
| `maxcore` | `6000` | Memory per core (MB) |
| `parallel_workflows` | `yes` | Run ox/red workflows in parallel (auto-splits PAL) |
| `pal_jobs` | `3` | Number of parallel job slots |
| `enable_job_timeouts` | `no` | Enable job timeouts |
| `job_timeout_hours` | `36` | Default job timeout |
| `opt_timeout_hours` | `14` | Geometry optimisation timeout |
| `frequency_timeout_hours` | `36` | Frequency calculation timeout |
| `sp_timeout_hours` | `3` | Single-point timeout |

### GUPPY Settings

| Key | Default | Description |
|-----|---------|-------------|
| `GUPPY_RUNS` | `20` | Number of GUPPY sampling runs |
| `GUPPY_GOAT` | `0` | Number of GOAT runs per GUPPY sample |
| `GUPPY_PARALLEL_JOBS` | `4` | Parallel GUPPY jobs |
| `GUPPY_SEED` | `31` | Random seed |

### Error Recovery

| Key | Default | Description |
|-----|---------|-------------|
| `enable_auto_recovery` | `yes` | Enable automatic ORCA error recovery |
| `max_recovery_attempts` | `3` | Maximum retry attempts per job |
| `enable_adaptive_parallelism` | `yes` | Adapt parallelism on resource failures |
| `enable_performance_metrics` | `yes` | Track performance metrics |

### OCCUPIER Settings

| Key | Default | Description |
|-----|---------|-------------|
| `OCCUPIER_method` | `auto` | OCCUPIER method: `auto` (adaptive tree) or `manually` |
| `OCCUPIER_tree` | `own` | Tree type |
| `OWN_TREE_PURE_WINDOW` | `3` | Pure spin state window |
| `OWN_progressive_from` | `no` | Progressive tree expansion |
| `fob_equal_weights` | `yes` | Equal weights for FoB calculations |
| `frequency_calculation_OCCUPIER` | `no` | Run frequencies during OCCUPIER |
| `occupier_selection` | `tolerance` | Selection method: `tolerance`, `truncation`, or `rounding` |
| `occupier_precision` | `3` | Precision for selection |
| `occupier_epsilon` | `5e-4` | Energy tolerance (Hartree) |
| `clean_override_window_h` | `0.002` | Clean state override window (Hartree) |
| `clean_quality_improvement` | `0.05` | Quality improvement threshold (Hartree) |
| `clean_quality_good` | `0.05` | Good quality threshold (Hartree) |
| `maxiter_occupier` | `125` | Max SCF iterations for OCCUPIER |
| `geom_opt_OCCUPIER` | `OPT` | Geometry optimisation keyword |
| `pass_wavefunction` | `no` | Pass wavefunction between OCCUPIER stages |
| `approximate_spin_projection_APMethod` | `2` | Spin projection method (Noodleman/Ruiz/Yamaguchi) |

### ORCA Base Overrides

| Key | Description |
|-----|-------------|
| `keyword:basename=[]` | Override ORCA keywords for a specific job basename |
| `additions:basename=[]` | Add ORCA input blocks for a specific job basename |

### CO2 Coordination

| Key | Default | Description |
|-----|---------|-------------|
| `co2_coordination` | `off` | CO2 coordination mode |
| `co2_species_delta` | `0` | CO2 species delta |

---

## 5. CLI Reference

### Main command: `delfin`

```bash
delfin [WORKSPACE] [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `WORKSPACE` | Workspace directory (default: current directory) |
| `-D`, `--define[=FILE]` | Generate CONTROL.txt + input file and exit. `.xyz` files are auto-converted |
| `--overwrite` | Overwrite existing CONTROL.txt and input file |
| `--control FILE` | Use a specific CONTROL.txt |
| `--recalc` | Only rerun missing/incomplete jobs |
| `--occupier-override STAGE=INDEX` | Force OCCUPIER index for a stage during recalc |
| `--report [text\|docx]` | Regenerate report from existing outputs |
| `--imag` | Run IMAG elimination on existing outputs, then generate report |
| `--json` | Generate DELFIN_Data.json and exit |
| `--json-output FILE` | Custom path for DELFIN_Data.json |
| `--afp` | Generate AFP spectrum plot |
| `--afp-fwhm NM` | FWHM for AFP Gaussian broadening (default: 50.0 nm) |
| `-C`, `--cleanup` | Remove intermediate files and exit |
| `--purge` | Remove all DELFIN-generated files (keeps CONTROL.txt + input) |
| `--no-cleanup` | Keep intermediate files after run |
| `-V`, `--version` | Show version and exit |

#### Subcommands

```bash
delfin cleanup [--dry-run] [--workspace PATH] [--scratch PATH] [--orca]
delfin stop --workspace PATH
delfin co2 --define
delfin qm_check [TOOL ...]
delfin qm_run TOOL -- [ARGS]
```

### Companion CLI tools

| Command | Description |
|---------|-------------|
| `delfin-voila` | Launch dashboard as a standalone web app via Voila |
| `delfin-build INPUT` | Build metal complex stepwise using ORCA XTB DOCKER |
| `delfin-guppy INPUT` | Multi-start SMILES sampling with ranked XTB trajectories |
| `delfin-step` | Run a single registered step |
| `delfin-pipeline YAML` | Execute a YAML-defined multi-step pipeline |
| `delfin-json` | Collect project outputs into JSON |
| `delfin_ESD OUTPUT` | UV-Vis spectrum report from ORCA ESD output |
| `delfin_IR OUTPUT` | IR spectrum report from ORCA frequency output |
| `delfin_NMR OUTPUT` | NMR spectrum report |

### delfin-voila

```bash
delfin-voila                     # starts on localhost:8866
delfin-voila --port 9000         # custom port
delfin-voila --ip 0.0.0.0       # bind to all interfaces
delfin-voila --dark              # dark theme
```

### delfin-build

Build metal complexes stepwise from SMILES using ORCA XTB DOCKER workflow.

```bash
delfin-build [INPUT] [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `INPUT` | `input.txt` | Input file with SMILES string |
| `-d`, `--directory` | `builder` | Output directory |
| `-m`, `--multiplicity` | `1` | Spin multiplicity |
| `-p`, `--pal` | `32` | CPU cores |
| `--maxcore` | `1000` | Memory per core (MB) |
| `--dry-run` | — | Only create input files |
| `--goat` | — | Run GOAT global optimisation |
| `--no-ligand-goat` | — | Skip GOAT for initial ligand structures |
| `--uphill-include-metals` | — | Include metals in UPHILLATOMS for GOAT |
| `--step-goat-preoptimized` | — | Use preoptimized trajectory for GOAT steps |
| `--step-goat-swarm` | — | Use swarm trajectory for GOAT steps |
| `-v`, `--verbose` | — | Verbose output |

### delfin-guppy

Multi-start SMILES sampling with parallel xTB optimisation and energy ranking.

```bash
delfin-guppy [INPUT] [OPTIONS]
```

Key options:

| Option | Description |
|--------|-------------|
| `--runs N` | Number of sampling runs |
| `--parallel-jobs M` | Parallel xTB jobs |
| `--pal P` | Cores per xTB job |

### delfin-step

Run a single registered computational step.

```bash
delfin-step --list                    # list available steps
delfin-step xtb_opt --geometry input.xyz --charge 0 --cores 8
delfin-step orca_sp --geometry input.xyz --charge 0 --method B3LYP --basis def2-SVP
```

### delfin-pipeline

Execute a YAML-defined multi-step pipeline.

```bash
delfin-pipeline workflow.yaml --cores 8
delfin-pipeline workflow.yaml --cores auto --scheduled
```

---

## 6. Workflow Modes

DELFIN supports three workflow methods, selected via `method=` in CONTROL.txt.

### OCCUPIER

Best for transition-metal complexes and systems with non-innocent ligands. OCCUPIER systematically explores spin states, handles broken-symmetry DFT, and manages wavefunction passing between stages.

**Key features:**
- Automated spin-state screening across multiplicities
- Broken-symmetry workflows for antiferromagnetic coupling
- Adaptive tree-based sequence navigation (`OCCUPIER_method=auto`)
- Per-atom basis set assignment for metals and first coordination sphere
- Approximate spin projection (Noodleman/Ruiz/Yamaguchi)
- Automatic PAL splitting for parallel ox/red branches

**How it works:**
1. DELFIN determines whether the system has even or odd electron count
2. Runs through a sequence of (multiplicity, BrokenSym) combinations
3. Ranks results by energy, checking spin contamination (⟨S²⟩)
4. Selects the preferred electronic state
5. Propagates wavefunctions to the next redox step

### classic

Standard DFT workflow for straightforward systems. Each redox state is calculated independently with automatic multiplicity detection from electron count.

**Best for:**
- Conventional organic redox systems
- Systems where spin state is not in question
- Simpler, easier-to-interpret per-state calculations

### manually

Expert-driven workflow where you specify multiplicities and broken-symmetry assignments explicitly for each state via `multiplicity_0`, `BrokenSym_0`, `multiplicity_ox1`, etc.

**Best for:**
- Reproducing known literature setups exactly
- Curated spin-state assignments
- Non-standard electronic configurations

---

## 7. Optional Modules

### Excited-State Dynamics (ESD)

Calculates photophysical properties: ISC/RISC rates, fluorescence and phosphorescence lifetimes, E₀₀ energies, and ΔE(S-T) gaps.

```ini
ESD_modul=yes
ESD_modus=TDDFT
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1
ICs=S2>S1
emission_rates=f,p
```

Creates a dedicated `ESD/` directory with ESD-specific ORCA jobs.

### Thermodynamics / Stability Constants

Predicts log K values from Born-Haber thermodynamic cycles for coordination complexes.

**Auto mode** (from complex SMILES):

```ini
thermodynamics=yes
thermodynamics_mode=auto
```

**Reaction mode** (explicit reaction):

```ini
thermodynamics=yes
thermodynamics_mode=reaction
thermodynamics_reaction=a*{SMILES}+b*{SMILES}...>>>c*{SMILES}+d*{SMILES}...
thdy_smiles_converter=ARCHITECTOR
thdy_preopt=xtb
```

### TADF Screening (xTB-based)

Fast screening of TADF candidates using semi-empirical methods:

```ini
tadf_xTB=yes
```

Estimates S0/T1 optimisation, S1 via Stokes shift, and ΔE(S-T) gaps.

### Hyperpolarizability (xTB-based)

Static and frequency-dependent β tensors for NLO materials:

```ini
hyperpol_xTB=yes
hyperpol_xTB_wavelengths=1064,532
```

### IMAG Elimination

Automatically detects and eliminates imaginary frequencies from converged structures:

```ini
IMAG=yes
IMAG_scope=initial   # or "all"
allow_imaginary_freq=0
```

Can also be run standalone: `delfin --imag`.

---

## 8. Structure Generation & Sampling

### SMILES Conversion (Dashboard)

The dashboard Submit tab offers multiple conversion methods:

| Button | Method | Best for |
|--------|--------|----------|
| `CONVERT SMILES` | Full isomer/conformer search (RDKit + Open Babel) | Thorough exploration |
| `QUICK CONVERT` | Fast single-conformer generation | Quick previews |
| `CONVERT + UFF` | Full search + UFF force-field refinement | Refined geometries |
| `BUILD COMPLEX` | ORCA/XTB DOCKER stepwise assembly | Metal complexes (job) |
| `ARCHITECTOR` | Architector automated 3D generation | Metal complexes (instant) |
| `SUBMIT GUPPY` | Multi-start xTB sampling with ranked trajectories | Robust start structures |

### delfin-build (ORCA/XTB DOCKER)

Stepwise metal-complex assembly using ORCA's XTB DOCKER workflow. Ligands are docked one-at-a-time onto the metal centre.

```bash
delfin-build input.txt --goat --pal 16
```

### delfin-guppy

Multi-start SMILES sampling: generates multiple chemically distinct starting geometries, runs parallel xTB optimisations, and ranks by energy.

```bash
delfin-guppy input.txt --runs 20 --parallel-jobs 4
```

### CREST / XTB-GOAT

Conformer search and global optimisation are triggered via CONTROL.txt:

```ini
XTB_OPT=yes    # xTB geometry optimisation
XTB_GOAT=yes   # XTB-GOAT global optimisation
CREST=yes      # CREST conformer search
```

---

## 9. Dashboard

### Starting the dashboard

**From Python / Jupyter:**

```python
from delfin.dashboard import create_dashboard
ctx = create_dashboard(backend="auto")
```

**As standalone app via Voila:**

```bash
delfin-voila
delfin-voila --port 9000 --ip 0.0.0.0 --dark
```

`backend="auto"` selects SLURM if `sbatch` is available, otherwise local execution.

### Dashboard tabs

| Tab | Purpose |
|-----|---------|
| **Submit Job** | SMILES/XYZ input, 3D preview, SMILES conversion buttons, job submission |
| **Recalc** | Edit and resubmit existing CONTROL.txt |
| **ORCA Builder** | Interactive ORCA input generation with geometry preview |
| **TURBOMOLE Builder** | Turbomole define workflow (SLURM backends) |
| **Job Status** | Real-time queue monitoring (local/SLURM), resource usage, job cancellation |
| **Calculations** | File browser, search, recalculation trigger, energy statistics, NMR/ANMR workflows |
| **Archive** | Archive browser with statistics |
| **Agent** | AI co-pilot with multi-agent pipelines, dashboard control, analysis, research |
| **Settings** | Tool detection, install/update buttons, runtime configuration, agent settings |

### Ensemble NMR via Calculations tab

The Calculations browser offers two NMR workflows:

- **Calc NMR**: Single-structure ORCA NMR
- **Calc ANMR**: End-to-end ensemble workflow: CREST → CENSO → c2anmr → ANMR → Boltzmann-weighted spectra

---

## 10. Settings & Runtime Configuration

User settings are stored in `~/.delfin_settings.json` (outside the git repo, not overwritten by updates).

### Settings tab sections

**Workspace Paths:**
- Calculations root (default: `~/calc`)
- Archive root (default: `~/archive`)

**Transfer Target:**
- SSH host, user, port, remote path for remote archive/transfer

**Runtime / Execution:**
- Backend: `Auto` / `Local` / `SLURM`
- ORCA path (global, local-only, SLURM-only overrides)
- qm_tools root (default: `~/.delfin/qm_tools`)
- Local CPU / RAM limits
- SLURM submit template directory
- Site profile (e.g., `bwunicluster3`)

### Setup buttons

| Button | Action |
|--------|--------|
| **Save Settings** | Persist current values to `~/.delfin_settings.json` |
| **Reload** | Reload settings from disk |
| **Validate Setup** | Run runtime diagnostics (read-only) |
| **Scan ORCA** | Search standard locations for ORCA installations |
| **Detect local resources** | Detect CPUs and RAM |
| **Prepare qm_tools** | Stage bundled QM tools to user area |
| **Install qm_tools** | Install QM tools from staged bundle |
| **Update qm_tools** | Refresh DELFIN-managed QM tools |
| **Setup bwUniCluster** | Configure DELFIN for bwUniCluster (lighter setup) |
| **Verify bwUniCluster** | Read-only cluster readiness check |
| **Full bwUni install** | Full bwUniCluster installation path |

### ORCA resolution order

If no explicit path is set, DELFIN resolves ORCA in this order:

1. `DELFIN_ORCA_BASE` environment variable
2. `ORCA_BINARY` environment variable
3. `ORCA_PATH` environment variable
4. `PATH`
5. Auto-detection in standard locations

### Tool detection

```bash
delfin qm_check                     # check all known tools
delfin qm_check xtb crest stda      # check specific tools
delfin qm_run xtb -- --version      # run a tool through DELFIN's resolver
```

For the full Settings documentation, see [SETTINGS_AND_SETUP.md](SETTINGS_AND_SETUP.md).

---

## 11. Error Recovery & Retry System

DELFIN can automatically detect and fix common ORCA failures when `enable_auto_recovery=yes`.

### Enable

```ini
enable_auto_recovery=yes
max_recovery_attempts=3
```

### How it works

```
ORCA fails → detect error type → modify input → continue from .gbw + latest .xyz → retry
```

### Error types and fixes

| Error | Automatic Fix |
|-------|--------------|
| **SCF not converged** | SlowConv → VerySlowConv + KDIIS + damping |
| **TRAH segfault** | NoAutoTRAH |
| **Geometry not converged** | Smaller trust radius → loose criteria |
| **MPI crash** | Reduce cores |
| **Memory error** | Reduce maxcore and PAL |
| **LEANSCF failure** | Fallback SCF strategy |
| **Frequency failure** | Skip frequency step |
| **Transient system error** | Exponential backoff retry |

Recovery state is tracked in `.delfin_recovery_state.json`.

Retry inputs are written as `input.retryN.inp` — preserving comments, geometry, inline basis sets, `%basis`, `%ecp`, and `$new_job` sections.

For complete details, see [RETRY_LOGIC.md](RETRY_LOGIC.md).

---

## 12. Reporting & Export

### Automatic outputs

| File | Content |
|------|---------|
| `DELFIN.txt` | Text summary: redox potentials, preferred spin states |
| `OCCUPIER.txt` | OCCUPIER stage tracking and state selection |
| `ESD_report.txt` | ESD results (when ESD module is enabled) |
| `delfin_run.log` | Run-level log with timing and job status |
| `occupier.log` | OCCUPIER subprocess log |

### On-demand reports

```bash
delfin --report text          # regenerate DELFIN.txt from existing outputs
delfin --report docx          # generate DELFIN.docx with embedded spectra/tables
delfin --json                 # generate DELFIN_Data.json
delfin --afp                  # generate AFP_spectrum.png
```

### Spectrum tools

```bash
delfin_ESD output.out         # UV-Vis spectrum with transitions and oscillator strengths
delfin_IR output.out          # IR spectrum with vibrational modes and intensities
delfin_NMR output.out         # NMR spectrum report
```

---

## 13. Cluster & HPC Usage

### SLURM backend

Set `backend=slurm` in the Settings tab or let `auto` detect `sbatch`.

DELFIN supports:
- SLURM submit templates (configurable path)
- Automatic resource detection on SLURM/PBS/LSF clusters
- Site profiles (e.g., `bwunicluster3`)
- Scratch directory via `DELFIN_SCRATCH` environment variable

### bwUniCluster setup

1. **Verify**: Settings → `Verify bwUniCluster` (read-only check)
2. **Setup**: Settings → `Setup bwUniCluster` (lightweight DELFIN-side config)
3. **Full install**: Settings → `Full bwUni install` (complete installation including OpenMPI, ORCA, venv)

### Environment variables

| Variable | Purpose |
|----------|---------|
| `DELFIN_SCRATCH` | Scratch directory for intermediate files |
| `DELFIN_ORCA_BASE` | ORCA base path override |
| `DELFIN_CHILD_GLOBAL_MANAGER` | Internal: subprocess PAL coordination |

### Resource management

When `parallel_workflows=yes`, DELFIN automatically splits PAL between parallel ox/red workflows (e.g., PAL=12 → 6+6). The global job manager ensures total CPU allocation is never exceeded.

### Example submission scripts

See `examples/example_Job_Submission_Scripts/` for SLURM, PBS, and LSF templates.

---

## 14. Troubleshooting

### CONTROL.txt not found

```bash
delfin --define
```

### Input file not found

```bash
delfin --define=your.xyz    # auto-converts XYZ to input.txt
```

### ORCA not found

```bash
which orca                  # check PATH
delfin qm_check             # DELFIN's own check
```

In the dashboard: Settings → `Scan ORCA` → select installation → `Save Settings`.

### xTB / CREST / STDA tools missing

```bash
delfin qm_check xtb crest xtb4stda stda std2
```

Or install via dashboard: Settings → `Install qm_tools`.

### A job failed — can it be recovered?

Check if auto-recovery is enabled:

```ini
enable_auto_recovery=yes
max_recovery_attempts=3
```

Recovery state: `.delfin_recovery_state.json`

### Recalculate incomplete jobs

```bash
delfin --recalc
```

Preserves existing results, only reruns jobs whose `.out` files are missing or incomplete.

### Override OCCUPIER state selection during recalc

```bash
delfin --recalc --occupier-override red_step_2_OCCUPIER=2
```

### Parallel jobs not behaving as expected

Check these CONTROL.txt keys:
- `PAL` — total cores
- `pal_jobs` — parallel job slots
- `parallel_workflows` — `yes`/`no`
- `enable_adaptive_parallelism`

See [job_prioritization.md](job_prioritization.md) for scheduler details.

### Cleanup after a run

```bash
delfin --cleanup              # remove intermediates
delfin cleanup --orca          # also terminate ORCA processes
delfin --purge                 # remove all generated files (keeps CONTROL.txt + input)
```

### Keep intermediates for debugging

```bash
delfin --no-cleanup
```

---

## 15. Recipes & Examples

### Basic redox workflow (organic molecule)

```ini
input_file=input.txt
charge=0
solvent=acetonitrile
calc_initial=yes
oxidation_steps=1
reduction_steps=1
method=classic
functional=PBE0
disp_corr=D4
main_basisset=def2-SVP
PAL=16
maxcore=4000
```

### Transition metal OCCUPIER workflow

```ini
input_file=input.txt
charge=2
solvent=acetonitrile
calc_initial=yes
oxidation_steps=1
reduction_steps=1
method=OCCUPIER
functional=PBE0
disp_corr=D4
main_basisset=def2-SVP
metal_basisset=def2-TZVP
first_coordination_sphere_metal_basisset=yes
PAL=16
maxcore=4000
parallel_workflows=yes
enable_auto_recovery=yes
max_recovery_attempts=3
```

### Excited-state dynamics (closed-shell)

```ini
input_file=input.txt
charge=0
solvent=acetonitrile
calc_initial=yes
method=classic
ESD_modul=yes
ESD_modus=TDDFT
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1
ICs=S2>S1
emission_rates=f,p
PAL=16
maxcore=4000
```

### Thermodynamics (stability constant from SMILES)

```ini
SMILES=[Cu+2]([N]1=CC=CC=1)([N]2=CC=CC=2)([N]3=CC=CC=3)([N]4=CC=CC=4)([OH2])([OH2])
charge=2
solvent=water
thermodynamics=yes
thermodynamics_mode=auto
thdy_smiles_converter=ARCHITECTOR
thdy_preopt=xtb
method=OCCUPIER
PAL=16
maxcore=4000
```

### SMILES → quick structure preview

```bash
# Write SMILES to input.txt
echo "c1ccccc1" > input.txt
delfin --define
# In dashboard: Submit → paste SMILES → QUICK CONVERT
```

### Build metal complex from SMILES

```bash
echo "[Fe+2]([N]1=CC=CC=1)([N]2=CC=CC=2)([N]3=CC=CC=3)([N]4=CC=CC=4)([N]5=CC=CC=5)[N]6=CC=CC=6" > input.txt
delfin-build input.txt --goat --pal 16 --maxcore 1000
```

### Multi-start sampling with GUPPY

```bash
echo "[Ru+2](N1=CC=CC=1)(N2=CC=CC=2)(N3=CC=CC=3)(N4=CC=CC=4)(Cl)(Cl)" > input.txt
delfin-guppy input.txt --runs 20 --parallel-jobs 4
```

### Generate reports from completed run

```bash
delfin --report text          # text summary
delfin --report docx          # Word document
delfin --json                 # JSON data export
delfin --afp                  # AFP spectrum plot
```

### Pipeline workflow (YAML)

```yaml
# workflow.yaml
steps:
  - name: preopt
    type: xtb_opt
    geometry: input.xyz
    charge: 0
  - name: dft_sp
    type: orca_sp
    geometry: preopt.xyz
    method: B3LYP
    basis: def2-SVP
    depends_on: preopt
```

```bash
delfin-pipeline workflow.yaml --cores 8
```

---

## Further Reading

- [README](../README.md) — Feature overview, architecture, tool list
- [Settings and Setup](SETTINGS_AND_SETUP.md) — Detailed settings documentation
- [Methodology](methodology.md) — Scientific methodology and validation
- [Retry Logic](RETRY_LOGIC.md) — Error recovery details
- [Job Prioritization](job_prioritization.md) — Scheduler behaviour
- [Example submission scripts](../examples/example_Job_Submission_Scripts/README.md) — HPC templates
