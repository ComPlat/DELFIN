# DELFIN User Manual

This manual is a practical guide for using DELFIN. It covers installation, workspace setup, all CLI tools, the full `CONTROL.txt` reference, dashboard usage, and troubleshooting.

For scientific methodology and validation details, see [methodology.md](methodology.md).
For implementation and architecture, see the [README](../README.md).

---

## Table of Contents

1. [Installation](#1-installation)
2. [Getting Started](#2-getting-started)
3. [Input Files](#3-input-files)
4. [How DELFIN works](#4-how-delfin-works)
5. [CONTROL.txt Reference](#5-controltxt-reference)
6. [CLI Reference](#6-cli-reference)
7. [Workflow Modes](#7-workflow-modes)
8. [Optional Modules](#8-optional-modules)
9. [Structure Generation & Sampling](#9-structure-generation--sampling)
10. [Dashboard](#10-dashboard)
11. [Settings & Runtime Configuration](#11-settings--runtime-configuration)
12. [The AI Agent](#12-the-ai-agent)
13. [Error Recovery & Retry System](#13-error-recovery--retry-system)
14. [Reporting & Export](#14-reporting--export)
15. [Cluster & HPC Usage](#15-cluster--hpc-usage)
16. [Troubleshooting](#16-troubleshooting)
17. [Recipes & Examples](#17-recipes--examples)

---

## 1. Installation

### Requirements

- **Python 3.10 or 3.11**
- **ORCA 6.1.1** in your `PATH` — [free for academic use](https://orcaforum.kofo.mpg.de/app.php/portal)
- **Optional:** `xtb`, `crest` (for xTB/CREST workflows)
- **Optional:** `xtb4stda`, `stda`, `std2` (for xTB-based screening)
- **Optional:** `censo`, `anmr`, `c2anmr`, `nmrplot` (for ensemble NMR)
- **Optional:** OpenMPI (for parallel ORCA), `g-xtb`, `dftb+`, `mopac`, `packmol`, `Multiwfn`
- **Optional:** JupyterLab/Notebook or Voila (for dashboard)

### The installer (recommended)

One script sets up DELFIN, wires up ORCA, builds OpenMPI as ORCA needs it, and
installs the QM and analysis tools. It needs no root and no module system.

```bash
git clone https://github.com/ComPlat/DELFIN.git ~/software/delfin
bash ~/software/delfin/install.sh              # DELFIN, ORCA wiring, QM/analysis tools
bash ~/software/delfin/install.sh --all        # everything, ML and AI stacks too (several GB)
bash ~/software/delfin/install.sh --only crest,gxtb
bash ~/software/delfin/install.sh --dry-run    # print the plan, change nothing
bash ~/software/delfin/install.sh --update     # update DELFIN and every installed tool
bash ~/software/delfin/install.sh --repair     # check everything, fix what is broken
```

ORCA is licensed and never downloaded; point the script at an unpacked copy with
`--orca DIR|TARBALL`. What was left out can be added later with `--only`, from the
dashboard's Settings tab, or automatically when a calculation first needs it.
`python -m delfin.installer --list` shows what can be installed, `--status` what is.

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
pip install -e ".[agent,docs,dev]"      # the agent, the MCP servers, the test tools
pip install -e ".[analysis,mlp]"        # analysis wrappers and ML potentials
```

### External QM tool setup

The installer above does this. By hand, **from a source checkout** (the paths are
relative to it, so this does not work after a plain `pip install`):

```bash
source delfin/qm_tools/env.sh
bash delfin/qm_tools/install_qm_tools.sh
bash delfin/qm_tools/check_qm_tools.sh
```

This can also be done from the dashboard Settings tab (see [Section 10](#11-settings--runtime-configuration)).

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

The template is written with placeholders in square brackets. Every one of them
has to be replaced before the run starts, or DELFIN stops with
`Missing required CONTROL values for: …`. At minimum:

```ini
charge=0
solvent=water
method=OCCUPIER
calc_initial=yes
oxidation_steps=1
reduction_steps=1
```

A run that builds its structure from a SMILES also needs `smiles_converter`.

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

See [Section 4](#5-controltxt-reference) for the complete reference.

---

## 4. How DELFIN works

What happens between typing `delfin` and reading `DELFIN.txt`. Read this once and the rest of the manual is a reference; skip it and the settings are a list of words.

### A run, in order

1. **CONTROL is read and validated.** Every key is checked against its allowed values, missing required keys are named, and everything you left out is filled from the template. A bad value stops the run here, before anything is computed.
2. **Resources are settled.** `PAL` and `maxcore` come from CONTROL, or from the node when CONTROL leaves them empty. They are read **once** and handed to a global job manager that hands out cores for the whole run.
3. **The structure is prepared.** A SMILES is built into 3D, an XYZ is taken as it is; then the optional xTB / GOAT / CREST steps refine it. The result is `start.txt`.
4. **The redox method runs.** `classic`, `manually` or `OCCUPIER` — this is the part that submits ORCA jobs.
5. **The extras run**, where enabled: excited states, stability constants, CO₂, hyperpolarizability, TADF.
6. **The reports are written**: `DELFIN.txt`, `ESD.txt`, `DELFIN_Data.json`, the DOCX.

### Two files, one geometry

| File | Who writes it | What it is |
|------|---------------|------------|
| `input.txt` | you | Your input, a geometry or a SMILES. **DELFIN never overwrites it.** |
| `start.txt` | DELFIN | The working geometry, in DELFIN's coordinate format (element and three numbers per line, no XYZ header). Every later step reads and rewrites this one. |

Confusing the two is the most common mistake: editing `input.txt` after a run changes nothing until the structure is built again.

### The structure stage

A SMILES goes through the converter named by `smiles_converter`:

| Value | What it does |
|-------|--------------|
| `QUICK` | One embedding. Fast, one structure. |
| `NORMAL` | Multi-seed embedding with force-field refinement. One structure. The fallback when nothing is set. |
| `MANTA` | Builds the coordination manifold of a metal complex, ranks it, and picks a winner. See *Structure generation*. |
| `ARCHITECTOR` | The external Architector builder, for metal complexes. |

Then, if enabled and in this order: `XTB_preOPT` (a quick xTB optimisation), `global_optimizer=GOAT|CREST` (a conformer search), `XTB_SOLVATOR` (explicit solvent shells). Each writes its result back over `start.txt`, so the next step starts from the last one.

Two exceptions worth knowing: MANTA can finish with a GOAT-refined winner, in which case the separate GOAT step is skipped; and in OCCUPIER the solvator runs *after* the first stage, on its winning geometry.

### OCCUPIER

The core idea: **the spin state of a metal complex is not known in advance, so DELFIN computes several and lets the energies decide** — for every charge state, and carrying what it learnt into the next one.

**Vocabulary.** A **stage** is one charge state, and it has a folder: `initial_OCCUPIER`, `ox_step_1_OCCUPIER`, `red_step_1_OCCUPIER`, and so on. Inside a stage, a **FoB** is one candidate electron configuration — a multiplicity, optionally a broken-symmetry label `M,N`, and the earlier FoB whose geometry it starts from. The list of FoBs is the stage's **sequence**.

**How a stage runs.** Each FoB is an ORCA optimisation; FoBs whose parents are independent run at the same time. When all of them are done, one selector job compares them: lowest energy first, ties broken on spin quality, with configurable windows that let a noticeably cleaner solution win against a marginally lower one. The result is written into the stage's `OCCUPIER.txt`, whose last two lines name the winner — that file is what DELFIN reads back.

**The hand-over.** The winner's geometry and orbitals are copied out of the stage folder to `input_<stage>_OCCUPIER.xyz` and `.gbw` beside CONTROL.txt. The next stage starts from them with its charge shifted by one.

**Where the sequence comes from.** With `OCCUPIER_method=auto` you give only the sequence for the neutral species (`even_seq` / `odd_seq`); every later stage is derived from the stage before it. If a pure state won, broken symmetry is tried next; if a broken-symmetry state `BS(M,N)` won, its neighbours `BS(M±1,N)` and `BS(M,N±1)` are tried. When several configurations are thermally populated rather than one clearly winning, all of them seed the next stage. The generated sequence is written into the stage's own CONTROL.txt, so you can read afterwards what was tried. With any other value you give every sequence yourself.

**The frequency job.** After the stages, one ORCA job per stage runs an optimisation *with* frequencies on the winning geometry, in the run root: `initial.inp`, `ox_step_1.inp`, and so on. **Its Gibbs energies are the ones the redox potentials come from** — not the FoB energies. FoBs get frequencies only when `OCCUPIER_compare=G`, and then only to rank configurations against each other.

### Running things at the same time

`PAL` is the total core budget; `pal_jobs` caps how many ORCA jobs run at once. One pool hands out cores: a job that can use more gets more when nothing else is waiting, and gives them back when something is. Jobs that many others depend on are started first.

Oxidation and reduction do not wait for each other — both branch off the initial stage, so they run side by side and share the cores. Excited-state jobs join the same pool rather than queueing behind the redox ladder.

Every ORCA job runs in its own directory, on the scratch disk when one is configured (`DELFIN_SCRATCH`), and its results are copied back.

### When a job fails

Failure detection is not just ORCA's exit: an optimisation that ran out of cycles, a collapsed excited-state root and an unphysical rate count as failures even when ORCA said it terminated normally.

With `enable_auto_recovery=yes`, DELFIN classifies the failure — SCF not converged, geometry not converged, memory, MPI, frequency, excited-state problems, transient system errors — and writes a **new** input `<name>.retry1.inp` with a targeted change: tighter or looser convergence, another SCF strategy, fewer cores, more memory. The original input is never edited, and the failed output is kept as `<name>.old1.out`.

A retry continues from the job's **own** last orbitals and geometry, never from another job's. The same error type recurring escalates the strategy rather than repeating it, up to `max_recovery_attempts`.

Two failures end a job immediately: an input ORCA refuses (no retry changes that) and an error DELFIN cannot classify. Without `enable_auto_recovery=yes` there is no retry at all.

### After the jobs

**Redox potentials** come from the Gibbs energies of the charge states in `initial.out`, `ox_step_*.out`, `red_step_*.out`. `calc_potential_method` selects how: from the neutral species, step by step, or the mean of both.

**Excited states**, with `ESD_modul=yes`, add their own jobs in `ESD/`: the states you listed, then the crossing, conversion and emission rates between them. Note that in classic mode ESD's own S0 job replaces the initial job.

**IMAG** runs after each frequency job when `IMAG=yes`: a structure with an imaginary mode is displaced along it in both directions, the lower side is re-optimised, and the frequencies are recomputed — up to `IMAG_max_rounds` times. The saddle's results are archived rather than deleted.

**The reports** are written last: `DELFIN.txt` (the potentials and the provenance), `ESD.txt` (the rates), `DELFIN_Data.json`, and the DOCX with the plots. A completed run also records the CONTROL it ran with, which is what makes the next point possible.

### Recalculation

`delfin --recalc` continues a run instead of starting over. A job is kept when its output is complete **and** its input has not changed. "Changed" is decided by a fingerprint over the input file — with the core and memory lines removed, because those do not change a result — and over the files it depends on.

What an edited CONTROL causes:

| You changed | What happens |
|-------------|--------------|
| Cores, memory, timeouts, the report name | Nothing is recomputed |
| The structure keys (SMILES, converter, charge, solvent, the pre-optimisation steps) | The structure is built again, and everything that follows from it |
| Anything else | The ORCA inputs it reaches are written anew and compared; a job runs again only if its input actually came out different |

Within OCCUPIER this is split once more: keys that only the stage frequency jobs are written from do not touch the configurations, and the excited-state, IMAG and side-module keys reach no OCCUPIER input at all.

`--occupier-override <stage>=<index>` forces a different winner for a stage and recomputes what depends on it.

### Names that are easy to confuse

| | |
|---|---|
| `input.txt` / `start.txt` | Your input / DELFIN's working geometry |
| `initial_OCCUPIER/input.xyz` / `input0.xyz` | The live geometry of the first configuration / the stage's untouched start |
| `input_initial_OCCUPIER.xyz` / `initial.xyz` | The stage's winner / the copy the frequency job uses |
| `OCCUPIER.txt` / `DELFIN.txt` | Which configuration won, per stage / the potentials, for the run |
| FoB frequencies / the stage frequency job | Ranking configurations / the energies the potentials use |

---

## 5. CONTROL.txt Reference

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
| `implicit_solvation_model` | `CPCM` | Solvation model: `CPCM` (also spelled `C-PCM`) or `SMD` |
| `solvent` | (required) | Solvent name (e.g., `acetonitrile`, `dmf`, `dcm`, `thf`, `dmso`, `acetone`) |
| `XTB_SOLVATOR` | `no` | Enable xTB ALPB solvation |
| `number_explicit_solv_molecules` | `2` | Number of explicit solvent molecules |

### Global Geometry Optimisation

| Key | Default | Description |
|-----|---------|-------------|
| `xTB_method` | `XTB2` | xTB method for pre-optimisation |
| `smiles_converter` | (required for a SMILES run) | `QUICK`, `NORMAL`, `MANTA` or `ARCHITECTOR` (own section). `GUPPY` is still read as a spelling of `MANTA` |
| `XTB_preOPT` | `no` | Run xTB geometry optimisation before DFT |
| `global_optimizer` | (empty) | Global optimisation: `GOAT`, `CREST`, or none |
| `multiplicity_global_opt` | (empty) | Override multiplicity for pre-optimisation |

The older spellings are still accepted, so existing CONTROL.txt files keep
working unchanged: `XTB_preOPT` is read as the older `XTB_OPT`, and `XTB_GOAT=yes` /
`CREST=yes` as `global_optimizer=GOAT` / `global_optimizer=CREST`.

### Imaginary Frequency Elimination (IMAG)

| Key | Default | Description |
|-----|---------|-------------|
| `IMAG` | `yes` | Move a structure whose frequencies show an imaginary mode to a minimum |
| `IMAG_scope` | `all` | `all`: the initial structure and every redox step; `initial`: the initial structure only. Excited states of the ESD module are always treated when `IMAG=yes`, since a rate needs both states at minima |
| `IMAG_option` | `2` | How OCCUPIER schedules IMAG |
| `allow_imaginary_freq` | `-50` | Imaginary frequencies smaller in magnitude than this (cm⁻¹, written ≤ 0) are numerical noise: IMAG leaves them and ESD rates are computed. The old template value `0` is read as `-50`; to treat every imaginary mode, write e.g. `-0.1` |
| `IMAG_sp_energy_window` | `1e-5` | A displaced single point must lie this far (Eh) below a single point at the saddle to be taken — a noise floor; which modes are worth removing is `allow_imaginary_freq`'s question. The old template value `1e-3` is read as `1e-5` |
| `IMAG_optimize_candidates` | `no` | Optimise the displaced structures instead of single points |
| `IMAG_max_rounds` | `2` | Most rounds per structure; each round costs one frequency calculation |

One round: the geometry in the structure's `.hess` is displaced along the
imaginary mode both ways (the atom that moves most by 0.3 Å, times
`IMAG_displacement_scale`), a single point of the same state is computed at
each, and the lower one — if it lies below the saddle — is re-optimised with
frequencies. If neither lies lower, the displacement is halved, at most twice.
Everything runs on the structure's own input, in place: the method, reference
and per-atom basis stay as written, jobs appended to the input run again at
the new geometry, and afterwards `.out`, `.xyz`, `.hess` and `.gbw` all belong
to the final geometry. Each round's saddle is kept in `<step>_IMAG/round<n>/`.

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
| `oxidation_steps` | (empty) | **Which** oxidation steps to run, as a list: `1`, `1,2` or `1,2,3`. Not a count — `oxidation_steps=2` runs only the second step, whose chain assumes the first one ran |
| `reduction_steps` | (empty) | Which reduction steps to run, same form |
| `method` | (required) | Workflow method: `classic`, `manually`, or `OCCUPIER` |
| `calc_potential_method` | `2` | Potential calculation method |

### ESD Module (Excited-State Dynamics)

| Key | Default | Description |
|-----|---------|-------------|
| `ESD_modul` | `no` | Enable excited-state dynamics |
| `ESD_modus` | `TDDFT` | Method: `TDDFT`, `deltaSCF`, or `hybrid1` |
| `ESD_T1_opt` | `uks` | T1 optimisation method: `uks` or `tddft` |
| `ESD_frequency` | `yes` | Run frequency calculation for ESD states |
| `states` | (empty) | Electronic states to compute, e.g. `S1,T1,S2,T2`. S1–S6 and T1–T6; a listed `S0` is dropped because S0 is always computed |
| `ISCs` | (empty) | Intersystem crossing rates, e.g. `S1>T1,T1>S1` |
| `ICs` | (empty) | Internal conversion rates, e.g. `S1>S0`. ORCA's ESD(IC) ends in the reference state, so singlets go `Sn>S0` and triplets `Tn>T1` |
| `emission_rates` | (empty) | Emission rates: `f` (fluorescence), `p` (phosphorescence) |

These four are **empty on purpose**: excited-state work is opt-in, and `ESD_modul=yes`
alone computes nothing. List the states and transitions you want. ESD also requires
`method=classic` — with `OCCUPIER` or `manually` the run is refused.
| `phosp_IROOT` | `1,2,3` | Phosphorescence IROOT sublevels |
| `phosp_keywords` | (empty) | Additional phosphorescence keywords |
| `fluor_keywords` | (empty) | Additional fluorescence keywords |
| `TROOTSSL` | `-1,0,1` | TROOT spin sublevels |
| `addition_S0` | (empty) | Additional S0 ORCA keywords |
| `DOHT` | `TRUE` | Duschinsky/Herzberg-Teller coupling |
| `ESD_LINES` | `LORENTZ` | Broadening function |
| `ESD_LINEW` | `50` | Line width |
| `ESD_INLINEW` | `250` | Input line width |
| `ESD_NPOINTS` | `auto` | Points of the ISC/IC correlation-function grid; `auto` lets ORCA choose |
| `ESD_MAXTIME` | `auto` | Time window of the ISC/IC correlation function (a.u.); `auto` lets ORCA choose it from the linewidth. The older template's `12000` (290 fs) is read as `auto`: it cut the correlation function off early (formaldehyde ISC 2.3x too fast, IC rate negative) |
| `hybrid1_geom_MaxIter` | `60` | Max geometry iterations for hybrid1 |

### xTB Hyperpolarizability (sTD-DFT-xTB)

| Key | Default | Description |
|-----|---------|-------------|
| `hyperpol_xTB` | `no` | Enable xTB hyperpolarizability |
| `hyperpol_xTB_xyz` | `start.txt` | Input geometry |
| `hyperpol_xTB_preopt` | `none` | Pre-optimisation method |
| `hyperpol_xTB_engine` | `std2` | Calculation engine |
| `hyperpol_xTB_bfw` | `no` | Bandwidth-filtered weights |
| `hyperpol_xTB_wavelengths` | `1064` | Wavelengths in nm, comma-separated for several (`1064,532`); `none` or `static` for the static tensor only |
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
| `thermodynamics_reaction` | the template's own pattern | Reaction SMILES: `a*{SMILES}+b*{SMILES}...>>>c*{SMILES}+d*{SMILES}...` |
| `n_explicit_solvent` | `6` | Number of explicit solvent molecules |
| `logK_exp` | (empty) | Experimental log K for comparison |
| `thdy_smiles_converter` | `NORMAL` | Converter: `QUICK`, `NORMAL`, `MANTA` or `ARCHITECTOR` |
| `thdy_preopt` | `xtb` | Pre-optimisation: `none`, `xtb`, `crest` or `goat` |

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

These keys reach every `%tddft` block DELFIN writes: the S0 absorption check,
every excited-state optimisation, the deltaSCF and hybrid1 check jobs, and the
ESD(ISC/IC/FLUOR/PHOSP) rate jobs. `KEY=?` in CONTROL.txt prints what a key
does and runs on the default.

| Key | Default | Description |
|-----|---------|-------------|
| `TDDFT_nroots` | `15` | Number of excited states (NRoots). Must reach the highest root a job asks for. |
| `TDDFT_maxdim` | `auto` | Davidson expansion space **in units of NRoots**: ORCA holds MaxDim × NRoots vectors. `auto` writes ORCA's own default, 10 (upper end of the 5–10 the ORCA manual recommends). |
| `TDDFT_maxiter` | `500` | Most Davidson iterations (MaxIter). `auto` leaves ORCA's own limit (100 in ORCA 6). |
| `TDDFT_TDA` | `TRUE` | Tamm-Dancoff approximation; `FALSE` is full TD-DFT (RPA). Applies to the rate jobs too, so states and rates share one level. |
| `TDDFT_followiroot` | `true` | Follow the optimised state by overlap when roots reorder. |
| `TDDFT_SOC` | `false` | DoSOC in state and check jobs. ISC and phosphorescence switch it on themselves. |
| `TDDFT_additions` | *(empty)* | Any other ORCA `%tddft` keyword, verbatim, `;`-separated, into every block — e.g. `DoNTO true; NTOThresh 1e-4; ETol 1e-7; EWin -5,10`. Keywords with their own key above, or set per job (iroot, irootmult, triplets, sroot, troot, trootssl, nacme, etf), are refused. |

Older spellings still work and mean the same key: `TDDFT_TDDFT_maxiter` and
`ESD_TDDFT_maxiter` → `TDDFT_maxiter`; `ESD_nroots`, `ESD_maxdim`, `ESD_TDA`,
`ESD_followiroot`, `ESD_SOC` → their `TDDFT_` names. A rate job can still get
its own root count with `ESD_ISC_NROOTS`, `ESD_IC_NROOTS`, `ESD_FLUOR_NROOTS`
or `ESD_PHOSP_NROOTS`.

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
| `relativity` | `ZORA` | Relativistic method: `ZORA`, `X2C`, `DKH`, `DKH2` or `none` |
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

`E_ref` sits in the Redox steps section. The experimental comparison keys below
are **optional** and are not written by `delfin --define` — paste in the ones
you want and the report gains an "Experimental properties" section.

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
| `PAL` | `48` | Number of CPU cores for ORCA |
| `maxcore` | `4500` | Memory per core (MB) |
| `parallel_workflows` | `yes` | Run ox/red workflows in parallel (auto-splits PAL) |
| `pal_jobs` | `4` | Number of parallel job slots |
| `enable_job_timeouts` | `no` | Enable job timeouts |
| `job_timeout_hours` | `36` | Default job timeout |
| `opt_timeout_hours` | `14` | Geometry optimisation timeout |
| `frequency_timeout_hours` | `36` | Frequency calculation timeout |
| `sp_timeout_hours` | `3` | Single-point timeout |

### MANTA Settings

The builder and its funnel. `delfin --define` writes the whole block with its
defaults; `KEY=?` in a CONTROL file prints a key's own explanation.

**Building the manifold**

| Key | Default | Description |
|-----|---------|-------------|
| `MANTA_QUALITY` | `extreme` | How thorough: `fast`, `normal`, `max`, `extreme` |
| `MANTA_SEEDS` / `MANTA_NUM_CONFS` | (empty) | Override the seed and conformer counts the quality implies |
| `MANTA_CONSTRUCTION` | `champion` | The flag set the builder runs with: `champion`, `builder`, `default` |
| `MANTA_MAX_ISOMERS` | `0` | `0` is the complete manifold; a number shrinks the search, not just the answer |
| `MANTA_BINDING_MODES` | `yes` | Enumerate binding modes |
| `MANTA_HAPTO` | `auto` | Hapticity: `auto`, `off` (a hapto SMILES then fails instead of being approximated) |
| `MANTA_UFF` | `yes` | Force-field refinement of the organic periphery (never at the metal) |
| `MANTA_DETERMINISTIC` | `yes` | Same SMILES in, same manifold out |
| `MANTA_ENV` | (empty) | `KEY=VALUE` pairs passed to the builder verbatim — the escape hatch |

**Gates** (each drops frames that are wrong, never frames that are merely unusual)

| Key | Default | Description |
|-----|---------|-------------|
| `MANTA_CLEAN_GATE` | `yes` | Drop collapsed bonds, inter-ligand clashes, decoordinated metals |
| `MANTA_TOPOLOGY_GATE` | `yes` | Drop frames with a bond the consensus says is not there |
| `MANTA_DEDUP` | `yes` | Remove permutation duplicates |

**The funnel: screen → optimise → refine**

| Key | Default | Description |
|-----|---------|-------------|
| `MANTA_SCREEN` | (empty) | Single-point method: `none`, `clash`, `gfnff`, `gfn0`, `gfn1`, `gfn2`, `gxtb` |
| `MANTA_SCREEN_ABOVE` | `30` | Screen first only above this many frames |
| `MANTA_SCREEN_KEEP` | `all` | How many survivors go on to the optimisation |
| `MANTA_OPT` | `xtb` | Optimise the survivors: `none` or `xtb` |
| `MANTA_OPT_METHOD` | (empty) | Override the xTB method used there |
| `MANTA_MULTIPLICITIES` | `auto` | Which multiplicities are ranked; the ranked unit is a (frame, multiplicity) pair |
| `MANTA_REFINE` | `goat` | Refine the best: `none`, `goat`, `crest` |
| `MANTA_REFINE_TOPK` | `0` | How many of the best go into the refinement (ceiling 10) |
| `MANTA_RMSD_CUTOFF` | `0.3` | Deduplication distance in Å |
| `MANTA_ENERGY_WINDOW` | `25.0` | Energy window in kcal/mol |
| `MANTA_PARALLEL_JOBS` | `auto` | Parallel frame jobs; `auto` is PAL/4 |
| `MANTA_TIME_BUDGET` | `3600` | Seconds the build may take before it is stopped |

Every `MANTA_*` key falls back to its `GUPPY_*` predecessor, so a CONTROL file
written before the rename keeps working and means the same thing. Still read:
`GUPPY_RUNS` (`20`), `GUPPY_SEED` (`31`), `GUPPY_GOAT` (the refinement top-k, not a
number of runs), `GUPPY_PARALLEL_JOBS`, and the bare `GUPPY=yes`, which means
`smiles_converter=MANTA`.

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
| `OWN_TREE_PURE_WINDOW` | `3` | Pure spin state window |
| `OWN_progressive_from` | `no` | Progressive tree expansion |
| `fob_equal_weights` | `yes` | Equal weights for FoB calculations |
| `OCCUPIER_compare` | `FSPE` | Energy OCCUPIER ranks configurations by: `FSPE` (electronic) or `G` (Gibbs; adds a frequency run per candidate). The former key `frequency_calculation_OCCUPIER` is still accepted. |
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
| `keyword:<job>=[...]` | ORCA `!` keywords for the jobs `<job>` names |
| `additions:<job>=[...]` | ORCA `%` blocks, one-line `%` settings or `!` lines for the jobs `<job>` names |

`<job>` is a job's name (its `%base`, else its input file name: `initial`,
`ox_step_1`, `S0`, `S1`, `T1_TDDFT`, `S1_T1_ISC_msp1`, ...), an OCCUPIER
folder (`initial_OCCUPIER`: every run in it), a pattern (`S*`, `*_ISC*`), or
`all`. Several entries reach a job in the order `all`, patterns, own name, so
the most specific wins.

Values are merged the way ORCA reads an input (ORCA 6.1.1 manual, section 2.1):

- A keyword replaces the member of its family the job has — `VeryTightSCF`
  for `TightSCF`, `DEFGRID3` for the grid, `B3LYP` for the functional,
  `TightOpt` for `Opt`. ORCA does not take the last of two such keywords.
- A block variable replaces the job's own value or is added to the block;
  `%scf maxiter 500 end` changes MaxIter, it does not add a second `%scf`.
- What DELFIN decides per job is left alone, with a hint: `%pal`/`PALn`,
  `%maxcore`, `%base`, `%moinp`/`MORead`, coordinates, the run type
  (`SP`, `OptTS`, ...), and the ESD module's `%tddft` roots
  (`iroot`, `irootmult`, `triplets`, ...).
- For `all` and patterns, a method keyword (functional, basis, dispersion,
  ...) only replaces what a job has — it is never added to an xTB job — and a
  block that starts a calculation of its own (`%tddft`, `%eprnmr`, a bare
  `%cpcm`, ...) only joins jobs that already run it. TD-DFT settings that
  have a `TDDFT_*` key are set there.

```ini
keyword:all=[VeryTightSCF DEFGRID3]
additions:initial=[
%scf
  maxiter 500
end
]
additions:S1_TDDFT=["%tddft maxiter 300 end"]
```

The Submit tab reports a value ORCA could not read as an error and a name no
job has, or a setting DELFIN keeps, as a hint.

### CO2 Coordination

| Key | Default | Description |
|-----|---------|-------------|
| `co2_coordination` | `off` | CO2 coordination mode |
| `co2_species_delta` | `0` | CO2 species delta |

---

## 6. CLI Reference

### Main command: `delfin`

```bash
delfin [WORKSPACE] [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `WORKSPACE` | Workspace directory (default: current directory) |
| `-D`, `--define[=FILE]` | Generate CONTROL.txt + input file and exit. `.xyz` files are auto-converted |
| `--overwrite` | Overwrite an existing CONTROL.txt (an existing input file is kept either way) |
| `--control FILE` | Use a specific CONTROL.txt |
| `--recalc` | Keep finished jobs; compute missing, incomplete and changed ones. Smart mode is **on by default**; `DELFIN_SMART_RECALC=0` falls back to skipping by output alone |
| `--occupier-override STAGE=INDEX` | Force OCCUPIER index for a stage during recalc |
| `--report [text\|docx]` | Regenerate report from existing outputs |
| `--imag` | Run IMAG elimination on existing outputs, then generate report |
| `--json` | Generate DELFIN_Data.json and exit |
| `--json-output FILE` | Custom path for DELFIN_Data.json |
| `--afp` | Generate AFP spectrum plot |
| `--afp-fwhm NM` | FWHM for AFP Gaussian broadening (default: 50.0 nm); also used by `--report docx` |
| `-C`, `--cleanup` | Remove intermediate files and exit |
| `--purge` | Remove all DELFIN-generated files (keeps CONTROL.txt + input). Asks before it does |
| `--no-cleanup` | Keep intermediate files after run |
| `-V`, `--version` | Show version and exit |

#### Subcommands

```bash
delfin cleanup [--dry-run] [--workspace PATH] [--scratch PATH] [--orca]
delfin stop [--workspace PATH] [--signal INT|TERM|KILL] [--dry-run] [--cleanup] [--wait-seconds S]
delfin co2 [--define] [--recalc] [--charge N] [--multiplicity M] [--solvent S] [--metal M] [--broken_sym B]
delfin run_orca [file.inp] [-i FILE] [-o FILE]
delfin tadf_xtb ...            # xTB/sTDA TADF screening
delfin hyperpol_xtb ...        # xTB hyperpolarizability
```

Checking an installation:

```bash
delfin doctor [--json] [--scratch PATH]   # ORCA, xTB, OpenMPI, scratch, SLURM, keys, doc index
delfin qm_check [TOOL ...]                # how each QM binary is resolved
delfin qm_run TOOL [--cwd DIR] [--capture] -- [ARGS]
delfin mlp_check                          # ML-potential backends, PyTorch, CUDA
delfin analysis_check                     # Multiwfn, CENSO, ANMR, morfeus
delfin csp_check                          # Genarris
```

Excited-state dynamics has no subcommand: set `ESD_modul=yes` in `CONTROL.txt`
and run `delfin`. A word that is not a subcommand is read as the workspace
directory, so a typo looks like a missing CONTROL file; an unknown `--flag`
stops the run with exit code 2.

#### `delfin doctor` — installation self-check

`delfin doctor` verifies that the local installation is usable before you
start a run: it probes ORCA, xTB and OpenMPI via `--version`, tests that the
scratch directory is writable, and reports whether SLURM, the KIT-Toolbox
API key and the doc-search index are configured. It never starts a
computation and never opens a network connection. Each check prints one
`OK` / `MISSING` / `BROKEN` line with a detail, and `MISSING`/`BROKEN`
results additionally print a fix hint. Only `BROKEN` (present but not
working) makes the command exit non-zero — a machine without SLURM or
without a KIT key is a valid setup.

```bash
delfin doctor               # human-readable report
delfin doctor --json        # raw results as JSON (name/status/detail/fix_hint)
delfin doctor --scratch /path/to/scratch
```

### Companion CLI tools

| Command | Description |
|---------|-------------|
| `delfin-voila` | Launch dashboard as a standalone web app via Voila |
| `delfin-build INPUT` | Build metal complex stepwise using ORCA XTB DOCKER |
| `delfin-guppy INPUT` | Multi-start sampling with ranked xTB structures (the MANTA builder under its old name) |
| `delfin-guppy-batch FILE` | The same over every SMILES in a batch file (`--row` for SLURM arrays) |
| `delfin-manta` | MANTA on its own: the coordination-isomer × conformer manifold from a metal SMILES |
| `delfin-fukui` | Atomic Fukui indices from three ORCA single points |
| `delfin-step` | Run a single registered step |
| `delfin-pipeline YAML` | Execute a YAML-defined multi-step pipeline |
| `delfin-app` | Application registry: `list`, `template`, `run <keyfile>`, `describe` |
| `delfin-agent` | The AI agent in the terminal (own chapter) |
| `delfin-json DIR` | Collect a project's outputs into JSON (the directory is required) |
| `delfin_ESD OUTPUT` | UV-Vis report from an ORCA output (writes a DOCX) |
| `delfin_IR OUTPUT` | IR report from an ORCA frequency output (DOCX + PNG) |
| `delfin_NMR OUTPUT` | ¹H NMR report from an ORCA NMR output (PNG) |
| `delfin-docs-index` | Build the documentation search index from `literature/` |
| `delfin-docs-server` · `delfin-ops-server` · `delfin-tools-server` | DELFIN's three MCP servers |

### delfin-voila

```bash
delfin-voila                     # starts on 127.0.0.1:8866
delfin-voila --port 9000         # custom port
delfin-voila --dark              # dark theme
delfin-voila --keep              # run inside a tmux session, survives a dropped terminal
delfin-voila --resume SID        # reopen a previous agent session

# On HPC/login nodes, keep 127.0.0.1 and use an SSH tunnel.
# Direct network binds require an explicit security override:
delfin-voila --ip 0.0.0.0 --allow-remote-bind
```

Voila prints a URL containing an access token; token authentication is mandatory
and cannot be switched off. If the port is taken, the next free one is used — read
the printed URL rather than assuming 8866.

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
| `--parallel-jobs M` | How many run at once; they share the core budget |
| `--pal P` | **Total** core budget for all runs together, not per job |
| `--screen` / `--optimise` / `--refine` | The funnel stages, as in the `MANTA_*` keys |
| `--max-isomers` / `--rmsd-cutoff` / `--energy-window-kcal` | How much of the manifold survives |

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

## 7. Workflow Modes

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

## 8. Optional Modules

### Excited-State Dynamics (ESD)

Calculates photophysical properties: ISC/RISC rates, fluorescence and phosphorescence lifetimes, E₀₀ energies, and ΔE(S-T) gaps.

```ini
ESD_modul=yes
ESD_modus=TDDFT
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1
ICs=S1>S0,T2>T1
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
IMAG_scope=all
allow_imaginary_freq=-50
```

Can also be run standalone: `delfin --imag`.

---

## 9. Structure Generation & Sampling

### Which builder a run uses

One CONTROL key decides, and a run that builds from a SMILES must set it:

```ini
smiles_converter=[QUICK|NORMAL|MANTA|ARCHITECTOR]
```

| Value | What it does |
|-------|--------------|
| `QUICK` | One embedding. Fast, one structure. |
| `NORMAL` | Multi-seed embedding with force-field refinement. One structure. The fallback when nothing decides otherwise. |
| `MANTA` | Builds the coordination manifold of a metal complex and ranks it down to one geometry — see below. |
| `ARCHITECTOR` | The external Architector builder; needs a metal-containing SMILES. |

`GUPPY` is still read as a spelling of `MANTA`, as is the older bare `GUPPY=yes`.

### MANTA

MANTA is DELFIN's own coordination builder: from a metal SMILES it enumerates the
coordination isomers by symmetry rather than searching for them, seats each one on
an ideal polyhedron with metal–donor distances from covalent radii, and expands it
into conformers. **No force field touches the metal** — UFF and MMFF have no
transition-metal parameters, which is what distorts metal–donor distances and
L–M–L angles in conventional builders.

The manifold then goes through a funnel that the `MANTA_*` keys steer:

```
build the manifold → screen (one single point per frame and multiplicity)
                   → optimise (xTB per survivor)
                   → refine (GOAT or CREST on the best)
                   → winner written to start.txt
```

The ranked unit is a (frame, multiplicity) pair, not a frame — the same geometry is
compared in several spin states. Working files land in a `GUPPY/` folder, the
winner in `best_coordniation.xyz` (the misspelling is the real file name).

**What it is for:** starting geometries, not production geometries. The result is
correct in topology and coordination, at roughly force-field quality — not
xTB- or DFT-accurate.

If MANTA's refinement already produced a GOAT-optimised winner, the separate
`global_optimizer=GOAT` step is skipped.

### In the dashboard

| Button | What it does | Best for |
|--------|--------------|----------|
| `CONVERT SMILES` | Isomer and conformer search (RDKit; deterministic by default) | Thorough exploration |
| `QUICK CONVERT SMILES` | One conformer | Quick previews |
| `CONVERT SMILES + UFF` | The same search with force-field refinement | Refined geometries |
| `MANTA` | The coordination manifold, ranked by GFN2 energy, with isomer navigation in the viewer | Metal complexes |
| `BUILD COMPLEX` | Stepwise assembly with ORCA's `%DOCKER` | Metal complexes (submitted as a job) |
| `ARCHITECTOR` | Architector 3D generation | Metal complexes (instant preview) |
| `SUBMIT GUPPY` | MANTA's sampling funnel as a submitted job | Robust start structures |

The first four live in the structure editor, which the Submit, Recalc and ORCA
Builder tabs all embed; the last three belong to the Submit tab.

### delfin-build (ORCA/XTB DOCKER)

Stepwise metal-complex assembly using ORCA's XTB DOCKER workflow. Ligands are docked one-at-a-time onto the metal centre.

```bash
delfin-build input.txt --goat --pal 16
```

### delfin-guppy

MANTA's sampling funnel on the command line, under its old name: it builds the
manifold, optimises the survivors in parallel with xTB and ranks them by energy.

```bash
delfin-guppy input.txt --runs 20 --parallel-jobs 4
delfin-guppy-batch batch.txt            # every SMILES in a file
delfin-manta "[Co](N)(N)(N)(N)(Cl)Cl"   # the builder alone, without the funnel
```

### CREST / XTB-GOAT

Conformer search and global optimisation are triggered via CONTROL.txt:

```ini
XTB_preOPT=yes           # xTB geometry optimisation
global_optimizer=GOAT    # XTB-GOAT global optimisation
global_optimizer=CREST   # ...or CREST conformer search
```

`XTB_OPT` / `XTB_GOAT` / `CREST` remain valid in older files.

---

## 10. Dashboard

### Starting the dashboard

**From Python / Jupyter:**

```python
from delfin.dashboard import create_dashboard
ctx = create_dashboard(backend="auto")
```

**As standalone app via Voila:**

```bash
delfin-voila
delfin-voila --port 9000 --dark
```

`backend="auto"` selects SLURM if `sbatch` is available, otherwise local execution.

### Dashboard tabs

| Tab | Purpose |
|-----|---------|
| **Submit Job** | SMILES/XYZ input, 3D preview, SMILES conversion buttons, job submission |
| **Recalc** | Edit and resubmit existing CONTROL.txt (smart recalc: what the edit changes is computed) |
| **ORCA Builder** | Interactive ORCA input generation with geometry preview |
| **TURBOMOLE Builder** | Turbomole define workflow (SLURM backends) |
| **Job Status** | Real-time queue monitoring (local/SLURM), resource usage, job cancellation |
| **Calculations** | File browser, search, recalculation trigger, energy statistics, NMR/ANMR workflows |
| **Archive** | Archive browser with statistics |
| **Remote Archive** | Browse and transfer results on another machine over SSH |
| **ChemDarwin** | Reaction-SMARTS structure enumeration and chemical-space map |
| **DELFIN Agent** | The AI agent — see its own chapter |
| **Agent Activity** | Running and recent agent and subagent work |
| **Office** | Documents and spreadsheets, with an agent that has no chemistry tools |
| **Literature** | The indexed literature corpus |
| **Tools** | The registered tool steps and their parameters |
| **Ketcher** | 2D structure drawing |
| **Reactions** | Reaction graph |
| **Pipelines** | Declarative step pipelines |
| **Settings** | Tool detection, install/update buttons, runtime configuration, agent settings |

Not every tab is visible: some are hidden until their module or key is present,
and the set can be changed under Settings → Dashboard Tabs. What you see is
therefore a subset of this list.

### Ensemble NMR via the Calculations tab

Select an `.xyz` file; the per-file dropdown then offers:

- **Calc NMR** — single-structure ORCA NMR
- **Calc CENSO/ANMR** — the ensemble chain: CREST → CENSO → c2anmr → ANMR

For an ORCA output file the dropdown offers **Print Mode**, **MO Plot** and
**Print NMR** instead.

---

## 11. Settings & Runtime Configuration

User settings are stored in `~/.delfin_settings.json` (outside the git repo, not overwritten by updates).

### Settings tab sections

**Workspace Paths:**
- Calculations root (default: `~/calc`, or the repository's own `calc/` when it exists;
  the archive is then its sibling)
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
| **Install qm_tools** | Download and build the QM tools (the staged bundle is what *Prepare* copies) |
| **Update qm_tools** | Refresh DELFIN-managed QM tools |
| **Setup cluster** | Configure DELFIN for the detected cluster (the light setup) |
| **Verify install** | Read-only readiness check |
| **Full install** | The complete installation, including OpenMPI, ORCA and the venv |

### ORCA resolution order

DELFIN takes the first of these that points at a usable ORCA:

1. A path given explicitly to the call
2. The **backend-specific** setting — `runtime.local.orca_base` for local runs,
   `runtime.slurm.orca_base` for SLURM. This beats the global setting.
3. The global `runtime.orca_base` setting
4. The environment: `DELFIN_ORCA_BASE`, then `ORCA_BINARY`, then `ORCA_PATH`
5. Installations the dashboard found beside the notebook directory (SLURM backend),
   preferring a 6.1.1
6. `orca` on `PATH`

### Tool detection

```bash
delfin qm_check                     # check all known tools
delfin qm_check xtb crest stda      # check specific tools
delfin qm_run xtb -- --version      # run a tool through DELFIN's resolver
```

For the full Settings documentation, see [SETTINGS_AND_SETUP.md](SETTINGS_AND_SETUP.md).

---

## 12. The AI Agent

A conversational agent that operates DELFIN and edits its code. It runs in the dashboard's **DELFIN Agent** tab and as the terminal program `delfin-agent`. Both surfaces share the same engine, modes and slash commands.

### Starting it

```bash
delfin-agent                     # chat in the current directory
delfin-agent run "<task>"        # one turn, headless
delfin-agent init                # scaffold AGENTS.md + .delfin/ in a project
delfin-agent doctor              # models, keys, sandbox, MCP servers
```

In the dashboard: the **DELFIN Agent** tab. Model, mode and permission level are chosen per session.

### Modes

A mode decides the role prompt and the tool surface. `plan` is **not** a mode but a permission level.

| Mode | Where | For |
|------|-------|-----|
| `dashboard` (default) | dashboard, terminal | Operating DELFIN: settings, submission, results |
| `solo`, shown as **Code** | dashboard, terminal | Reading and writing code |
| `office` | dashboard, terminal | Documents and spreadsheets, without the chemistry tools |
| `research` | terminal only | Literature and method research |

Switch with `/mode`. Older mode names (`quick`, `reviewed`, `tdd`, `cluster`, `full`, `pipeline`) are retired and resolve to `solo`; the multi-agent review pipeline they belonged to no longer exists.

### Models

| Provider | Authentication | Notes |
|----------|----------------|-------|
| `claude` (default) | `ANTHROPIC_API_KEY` | Uses the `claude` CLI when it is on PATH, otherwise the API |
| `openai` | `OPENAI_API_KEY` | Uses the `codex` CLI when it is on PATH, otherwise the API |
| `kit` | `KIT_TOOLBOX_API_KEY` | University-hosted, OpenAI-compatible |
| `ollama` | none | Any OpenAI-compatible endpoint: Ollama, vLLM, LM Studio, llama.cpp-server. `OLLAMA_HOST`, default `http://localhost:11434` |

Choose with `--provider` / `--model` or `/model`. DELFIN detects each model's real context window and its tool, vision and reasoning support — live from the endpoint where that is possible, otherwise from a built-in table — and sizes context management to it. For Ollama it sends the matching `num_ctx`, so a local model is not silently capped.

**A model that cannot call tools is refused.** Every non-trivial action is a tool call, so such a model could only talk. Small models automatically get a shorter prompt and a reduced tool surface.

### Permissions

| Level | Writes | Shell |
|-------|--------|-------|
| `plan` | refused | read-only |
| `default` | destructive actions ask first | asks, except for an allow-list |
| `diff_approval` | staged as a diff for `/approve` | asks |
| `acceptEdits` | allowed inside the workspace | asks |
| `bypassPermissions` | allowed, no prompts | no prompts |

Cycle with `/permissions` or Shift+Tab; the cycle deliberately never lands on bypass. Two things hold at **every** level: the sandbox and the deny lists stay in force, and a file that defines the agent's own permissions always requires explicit confirmation.

Settings, hook commands and MCP servers that come from a checked-out repository are ignored until you trust that directory (`/trust`).

### The shell sandbox

Every command runs through an allow-list and then, where available, bubblewrap or firejail. Credential directories (`~/.ssh`, `~/.aws`, `~/.gnupg`, …) are masked, the network is denied by default, and every command is appended to `~/.cache/delfin/agent-audit.jsonl`. `DELFIN_AGENT_SANDBOX={auto,bwrap,firejail,allowlist,off}` selects the mechanism; when isolation is unavailable the agent says so rather than running unprotected in silence.

### Subagents

Self-contained work is delegated to an isolated agent with its own tool loop and usually tighter permissions. The parent sees only the final answer.

| Preset | Permissions | For |
|--------|-------------|-----|
| `explore` | read-only | Investigation, reports findings |
| `plan` | read-only | A plan, no edits |
| `code-reviewer` | read-only | Independent review |
| `general-purpose` | inherited | A self-contained task |

Launch with `/explore`, `/review`, `/plan`, `/delegate <task>`, or let the agent delegate on its own. Several run in parallel — a writing preset gets its own git worktree, so concurrent edits cannot collide — or in the background, and a finished one can be continued with its context intact. Each run is bounded in wall-clock time, tool calls and output, and a subagent may not spawn subagents of its own unless that depth is raised; both in Settings.

Your own presets are markdown files with frontmatter in `~/.delfin/subagents/`.

### Memory and sessions

| What | Where | Lives |
|------|-------|-------|
| Facts and preferences | `~/.delfin/agent_memory.json` | Across sessions, recalled by relevance |
| Project notes | `~/.delfin/projects/<project>/memory/` | Per project |
| Project instructions | `DELFIN.md` or `AGENTS.md` in the project | Loaded automatically |
| Conversations | `~/.delfin/agent_sessions/` | Saved, restored, forked, exported |

`/remember`, `/memories`, `/forget` manage memory; a line starting with `#` is a note to memory. Memory is scoped, so a preference recorded while coding does not surface in an office session. Long conversations are compacted token-aware; `/compact` and `/context` show the state.

### Commands

| Group | Commands |
|-------|----------|
| Session | `/help` `/status` `/cost` `/usage` `/context` `/compact` `/export` `/undo` `/clear` |
| Setup | `/model` `/mode` `/permissions` `/effort` `/doctor` `/init` `/mcp` `/trust` `/tools` `/agents` `/skills` |
| History | `/session` `/rewind` `/tasks` `/memories` `/plans` |
| Changes | `/pending` `/approve` `/reject` `/undo-file` |
| Workspace | `/git` `/bash` `/jobs` `/attention` |
| DELFIN (dashboard) | `/control key <key> <value>` · `/orca set` · `/calc ls\|read\|info\|tree` · `/analyze energy\|convergence\|errors` · `/recalc check-all\|auto` · `/submit` `/cancel` |

`!` runs a shell command directly. A markdown file in `~/.delfin/commands/` (user-wide) or `<workspace>/.delfin/commands/` (project) becomes a slash command of its own; placeholders `$ARGUMENTS`, `$1`, `$2` are substituted. A **skill** (`~/.delfin/skills/<name>/SKILL.md`) can be invoked by you or chosen by the model; DELFIN ships skills for diagnosing a failed run, recalculating one, setting up TD-DFT and CASSCF, frequency thermochemistry, solvation, and tuning CONTROL.

### MCP servers

External Model Context Protocol servers are configured in `~/.delfin/mcp_servers.json` or `<workspace>/.delfin/mcp_servers.json` and reached over stdio or HTTP; their tools, resources and prompts join the agent's surface.

A tool reached through MCP runs in a process DELFIN did not start a command line for, so the shell sandbox is not around it. Give a stdio server `"roots": [...]` (read-write) or `"read_roots": [...]` (read-only) and it is started inside a namespace holding only those paths. A server that declares neither runs uncontained, and the startup banner, `/mcp` and `delfin-agent doctor` all name it as such.

DELFIN ships three servers of its own: `delfin-tools-server` (the tool platform), `delfin-docs-server` (literature and calculation search) and `delfin-ops-server` (typed runtime actions). They can take their roots from your settings with `agent.mcp_isolation: "builtin"`.

### Seeing what runs, and stopping it

```bash
delfin-agent stop-all --check    # what is open, kept, scheduled or running; changes nothing
delfin-agent sessions            # earlier sessions: id, age, model, task
delfin-agent where               # which node the dashboard runs on, and how to get back in
delfin-agent jobs                # watched cluster jobs
```

`stop-all --check` reads open dashboard sessions, kept sessions, schedules and PID files from the **shared home**, so they are visible from any login node; **agent processes** it can only see on the machine you run it on. Run it on every login node you have used.

Without `--check` the same command is the emergency stop: it ends every agent of yours, and it **disables every schedule**. That is permanent — a disabled schedule is not resumed, it has to be created again.

The dashboard's **Agent Activity** tab shows running and recent agent and subagent work.

---

## 13. Error Recovery & Retry System

DELFIN can automatically detect and fix common ORCA failures when `enable_auto_recovery=yes`.

### Enable

```ini
enable_auto_recovery=yes
max_recovery_attempts=3
```

### How it works

```
ORCA fails → detect error type → modify input → retry
```

Most strategies restart the job from its own last orbitals (`.gbw`) and, for an
optimisation, its own last geometry. Two do not: after an MPI crash the orbitals
are usually not written yet, and a memory error is not helped by reading more in.

### Error types and fixes

| Error | Automatic Fix |
|-------|--------------|
| **SCF not converged** | SlowConv → VerySlowConv + KDIIS + damping → SOSCF (`CNVSOSCF true`) |
| **TRAH segfault** | NoTRAH + SlowConv |
| **Geometry not converged** | Smaller trust radius → loose criteria |
| **MPI crash** | Fewer cores and OpenMPI without single-copy transport, finally serial (one core) |
| **Memory error** | Fewer cores, SOSCF removed, and MaxCore raised to what ORCA asked for when it named a figure |
| **LEANSCF failure** (every SCF in ORCA 6) | Outside a Hessian: the same escalation as SCF not converged. Inside one: its own escalation, and from the third attempt the frequency step is dropped |
| **Frequency failure** | NumFreq, then skip the frequency step |
| **CIS/TD-DFT failure** | TDA, then tighter SCF and a larger Davidson space, then half the cores |
| **TD-DFT root collapsed** (ORCA ends "normally") | Treated as a failure and rerun with the TD-DFT instability strategy |
| **DIIS error** | KDIIS with SlowConv, then SOSCF with damping, then stronger damping |
| **Optimisation ran out of cycles** (ORCA ends "normally") | Continue from the last geometry with more cycles |
| **ESD rate negative** (ORCA ends "normally") | More points, then more points again with a longer window |
| **ESD window truncated** (ORCA ends "normally") | The time window is handed back to ORCA; there is no second attempt |
| **Input ORCA refuses** (e.g. impossible multiplicity) | Named in the log, not retried |
| **An error DELFIN cannot classify** | Named in the log, not retried |
| **Transient system error** | Exponential backoff retry |

Recovery is **off** unless `enable_auto_recovery=yes` is set; without it a failed
ORCA job simply fails. The attempt counts are kept per job and error type in
`.delfin_recovery_state.json` in the job's own folder, capped by
`max_recovery_attempts`.

Each failed attempt leaves its output behind as `<job>.old<N>.out`. Together with
the rewritten `<job>.retryN.inp` those are the two files to read when you want to
see what recovery tried.

Retry inputs are written as `<job>.retryN.inp`. Only the job that failed is
changed (in a file with `$new_job`, the one ORCA was running); every other
line stays as it was — comments, inline basis sets (`NewGTO` on a metal's
line), `%basis`, a second `%scf` with `BrokenSym`. A retry writes its files
under the job's own name, restarts from the job's own `.gbw` (never another
job's), and an optimisation from its own last geometry. A retry ORCA could not
read is not written.

For complete details, see [RETRY_LOGIC.md](RETRY_LOGIC.md).

---

## 14. Reporting & Export

### Automatic outputs

Every run writes these, in this order, at the end:

| File | Where | Content |
|------|-------|---------|
| `DELFIN.txt` | run root | Redox potentials, preferred spin states, method provenance |
| `ESD.txt` | run root | Excited-state results, when the ESD module produced data |
| `DELFIN_Data.json` | run root | Every energy, rate and spectrum, machine-readable |
| `DELFIN.docx` | run root | The Word report, with the plots embedded |
| `delfin_run.log` | run root | Run-level log with timing and job status |
| `OCCUPIER.txt` | **each stage folder** | Which configuration won that stage, and why |
| `occupier.log` | **each stage folder** | That stage's log |

The DOCX also leaves its plots in the run root as PNG and DOCX files —
`AFP_spectrum.png`, `Energy_Level_Diagram.png`, `Vertical_Excitation_Energies.png`,
`Correlation_Diagram.png`, `Dipole_Moment_Visualization.png`, the UV/Vis and IR
documents, and the electrostatic-potential images.

### Producing a report again, without calculating

```bash
delfin --report text          # rewrite DELFIN.txt from the existing outputs
delfin --report docx          # rebuild DELFIN.docx
delfin --json                 # rebuild DELFIN_Data.json
delfin --afp                  # rebuild AFP_spectrum.png
```

`--report text` works in the **current** directory and leaves out the excited-state
section and the E₀₀ values; `--report docx` is the complete one.

### Spectrum tools

```bash
delfin_ESD output.out         # UV-Vis: writes Absorption_Spectrum_S0.docx
delfin_IR output.out          # IR: writes IR_<name>.docx and IR_<name>.png
delfin_NMR output.out         # 1H NMR: writes NMR_<name>.png
```

---

## 15. Cluster & HPC Usage

### SLURM backend

Set `backend=slurm` in the Settings tab or let `auto` detect `sbatch`.

DELFIN supports:
- SLURM submit templates (configurable path)
- Site profiles (e.g. `bwunicluster3`), detected from the login node
- SLURM, PBS and LSF are recognised; `PAL` and `maxcore` are taken from the
  allocation only when CONTROL leaves them empty — and the shipped template sets
  both, so in practice your CONTROL decides
- A scratch directory, see below

### Cluster setup from the dashboard

Settings offers three buttons, in increasing order of what they touch:

1. `Verify install` — read-only check
2. `Setup cluster` — the DELFIN-side configuration only
3. `Full install` — the complete installation, including OpenMPI, ORCA and the venv

### Scratch, and where your results are while the job runs

The submit script **sets `DELFIN_SCRATCH` itself** — to a BeeOND mount, else to
`$TMPDIR`, else to `/scratch/<user>/delfin_<jobid>` — so exporting it by hand
before `sbatch` has no effect. It is **deleted when the job ends**.

With `DELFIN_STAGE_IO=1` (the default) the workspace is copied to node-local disk
and synced back on an interval (`DELFIN_SYNC_INTERVAL`, 900 s). Until a sync fires,
the submit directory does not yet show the newest results — the usual surprise on
a first cluster run.

### Environment variables

| Variable | Purpose |
|----------|---------|
| `DELFIN_SCRATCH` | Scratch directory (set by the submit script on a cluster) |
| `DELFIN_ORCA_BASE` | ORCA base path override |
| `DELFIN_STAGE_IO` / `DELFIN_SYNC_INTERVAL` | Stage the workspace on node-local disk, and how often to sync back |
| `DELFIN_STAGE_ORCA` / `DELFIN_STAGE_VENV` | Stage ORCA and the venv on node-local disk |
| `DELFIN_PARTITION` / `DELFIN_HIGHMEM_PARTITION` | Force a partition |
| `DELFIN_NODE_CORES` / `DELFIN_NODE_MEM_MB` | Tell the script the node size |
| `DELFIN_VENV` / `DELFIN_REPO` / `DELFIN_OMPI_HOME` | Where the job finds DELFIN |
| `DELFIN_PAL` / `DELFIN_MAXCORE` / `DELFIN_CONTROL` | Override the run's resources and CONTROL file |
| `DELFIN_SMART_RECALC` | `0` turns smart recalc off |
| `DELFIN_SLURM_MEM_HEADROOM` | Fraction of the node's memory DELFIN may hand out (default 0.90) |
| `DELFIN_CHILD_GLOBAL_MANAGER` | Internal: subprocess PAL coordination |

The submit script's own header documents the rest.

### Resource management

`PAL` is the total core budget and `pal_jobs` caps how many ORCA jobs run at once.
Oxidation and reduction do not get a fixed half each: one pool hands out cores as
jobs become runnable, gives a lone job more, and takes them back when something
else can start. The total allocation is never exceeded.

### Example submission scripts

See `examples/example_Job_Submission_Scripts/` for SLURM, PBS, and LSF templates.

---

## 16. Troubleshooting

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

A recalc keeps finished jobs and computes what is missing or incomplete. What
else it computes depends on the mode:

- **Smart recalc** (default; `DELFIN_SMART_RECALC=1`; the Recalc tab and
  *Smart Recalc* in the Calculations browser): a finished job is also run
  again when its input changed. Every completed run records the CONTROL.txt
  and geometry input it was computed with in `.delfin_last_run.json`; the
  next recalc compares against it.
  - Only cores, memory, timeouts or the recovery budget changed: nothing is
    recomputed.
  - A key that reaches the calculations changed (functional, basis, OCCUPIER
    settings, `keyword:`/`additions:` overrides, …): the ORCA inputs,
    OCCUPIER's included, are written anew; a job whose input came out
    different runs again, and so does everything computed from it.
  - A key that builds the starting structure changed (SMILES, converter,
    `MANTA_*`, xTB pre-optimisation, GOAT/CREST, charge, solvent) or the
    geometry input itself: the structure is built again and everything
    follows. Otherwise the structure in `start.txt` is kept — MANTA does not
    build the same structure twice.
  - A job from before the record existed has none; editing its CONTROL.txt in
    the Recalc tab records the file being replaced, so the edit is known.
- **Classic recalc** (`DELFIN_SMART_RECALC=0`; *Recalc* in the Calculations
  browser): every finished job is kept, whatever changed.

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

## 17. Recipes & Examples

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
ICs=S1>S0,T2>T1
emission_rates=f,p
PAL=16
maxcore=4000
```

### Thermodynamics (stability constant from SMILES)

```ini
SMILES=[Cu+2]([N]1=CC=CC=C1)([N]2=CC=CC=C2)([N]3=CC=CC=C3)([N]4=CC=CC=C4)([OH2])([OH2])
charge=2
solvent=water
thermodynamics=yes
thermodynamics_mode=auto
smiles_converter=ARCHITECTOR
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
delfin --define          # writes a CONTROL.txt template; an existing input.txt is kept
```

`--define` writes the template with its placeholders still in it. Fill in at least
`charge`, `solvent`, `method` and `smiles_converter` before running `delfin`, or the
run stops with `Missing required CONTROL values for: …`.

In the dashboard: Submit → paste the SMILES → `QUICK CONVERT SMILES`.

### Build metal complex from SMILES

```bash
echo "[Fe+2]([N]1=CC=CC=C1)([N]2=CC=CC=C2)([N]3=CC=CC=C3)([N]4=CC=CC=C4)([N]5=CC=CC=C5)[N]6=CC=CC=C6" > input.txt
delfin-build input.txt --goat --pal 16 --maxcore 1000
```

### Multi-start sampling with GUPPY

```bash
echo "[Ru+2](N1=CC=CC=C1)(N2=CC=CC=C2)(N3=CC=CC=C3)(N4=CC=CC=C4)(Cl)(Cl)" > input.txt
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
name: preopt_then_sp
steps:
  - step: xtb_opt
    geometry: input.xyz
    charge: 0
  - step: orca_sp
    method: B3LYP
    basis: def2-SVP
```

A pipeline needs a top-level `name:`, and every entry names its step with `step:`.
(`type:` exists but is reserved for flow control — `if`, `loop`, `retry`, `map` and
so on.) Steps run in the order they are written; there is no `depends_on:`.
`delfin-step --list` prints the step names that can be used.

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
