<p align="center">
  <img src="delfin/logo/DELFIN_readme_demo.gif" alt="DELFIN logo animation" width="820">
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-LGPL--3.0--or--later-blue.svg?style=for-the-badge" alt="License"></a>
  <a href="https://github.com/ComPlat/DELFIN/graphs/commit-activity"><img src="https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=for-the-badge" alt="Maintenance"></a>
  <a href="https://github.com/ComPlat/DELFIN/issues"><img src="https://img.shields.io/github/issues/ComPlat/DELFIN.svg?style=for-the-badge" alt="GitHub issues"></a>
  <a href="https://github.com/ComPlat/DELFIN/graphs/contributors"><img src="https://img.shields.io/github/contributors/ComPlat/DELFIN.svg?style=for-the-badge" alt="GitHub contributors"></a>
  <a href="https://doi.org/10.5281/zenodo.21508364"><img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21508364-blue.svg?style=for-the-badge" alt="DOI"></a>
  <a href="https://pypi.org/project/delfin-complat/"><img src="https://img.shields.io/pypi/v/delfin-complat.svg?style=for-the-badge" alt="PyPI version"></a>
</p>

> 📄 **Preprint**: *Hartmann, M. et al. "DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials"*, ChemRxiv (2025). https://doi.org/10.26434/chemrxiv-2025-4c256 · https://www.cambridge.org/engage/chemrxiv/article-details/68fa0e233e6156d3be78797a

**DELFIN** is an open-source computational chemistry platform that automates molecular property prediction. From a SMILES string or an XYZ geometry it builds the structure, runs the quantum-chemistry workflow, and writes the results as a report — for organic molecules and for transition-metal complexes. It is driven from a browser dashboard, from the command line, or from Python, and it ships an AI agent that operates all three.

> 🧬 organic chemistry · 🧲 transition-metal complexes · 💡 photoactive materials · 🔋 redox systems · 🔬 spectroscopy · ⚛ excited-state dynamics

DELFIN is research infrastructure: it does not replace chemical judgement — it makes rigorous computational chemistry routine enough to inform it.

---

## Contents

1. [Overview](#overview) — [what it does](#what-delfin-does) · [ways to use it](#ways-to-use-it) · [design principles](#design-principles) · [use cases](#use-cases)
2. [Installation](#-installation)
3. [Quick start](#-quick-start) — [dashboard](#the-dashboard) · [command line](#the-command-line) · [Python](#python)
4. [Dashboard](#-dashboard)
5. [Workflows](#-workflows)
6. [Structure generation from SMILES](#-structure-generation-from-smiles)
7. [Configuration (CONTROL.txt)](#-configuration-controltxt)
8. [AI agent](#-ai-agent)
9. [CLI reference](#-cli-reference)
10. [Architecture](#-architecture)
11. [External programs](#-external-programs)
12. [Cluster and HPC](#-cluster-and-hpc)
13. [Troubleshooting](#-troubleshooting)
14. [Development](#-development)
15. [References](#references) · [Dependencies and legal notice](#dependencies-and-legal-notice) · [Please cite](#please-cite) · [License](#license)

---

## Overview

### What DELFIN does

Domain legend: 🧪 organic · 🧲 metal complex · 🔬 both / general · 🧱 solid state

| Capability | Domain | Description |
|------------|--------|-------------|
| **Redox potentials** | 🔬 both | Automated spin-state identification and redox-potential calculation via the OCCUPIER and classic workflows — for organic radicals and for transition-metal complexes |
| **Stability constants (log K)** | 🧲 metal | log K and ΔG for one coordination complex per run, from a Born-Haber-type cycle with OCCUPIER-aware metal, ligand and solvent sub-workflows |
| **Excited-state dynamics** | 🔬 both | ISC / RISC / IC rates, fluorescence and phosphorescence lifetimes (incl. per-sublevel), SOC coupling, E₀₀ and ΔE(S-T) — orchestrated from a list of states and transitions you give in `CONTROL.txt` |
| **TADF screening** | 🧪 organic | xTB/sTDA singlet-triplet gap estimation, oscillator strengths and estimated emission wavelength, for OLED candidate triage |
| **Spectroscopy** | 🔬 both | UV-Vis absorption, IR vibrational spectra, AFP (absorption/fluorescence/phosphorescence) overlay plots |
| **Ensemble NMR** | 🧪 organic | Boltzmann-weighted shieldings and couplings via CREST → CENSO → c2anmr → ANMR (all four are external programs) |
| **Hyperpolarizability** | 🔬 both | Static and frequency-dependent β tensors (β_HRS, β_VEC, DR) for NLO materials |
| **Fukui indices** | 🔬 both | Atomic Fukui indices from three ORCA single points (N, N±1) at fixed geometry, with ESP reporting |
| **Structure generation** | 🔬 both | SMILES→3D for organics (RDKit / Open Babel) and for metal complexes: **MANTA** (DELFIN's own force-field-free coordination builder), Architector, and stepwise ORCA/xTB `%DOCKER` assembly |
| **Conformer sampling** | 🔬 both | MANTA multi-start sampling with energy ranking, CREST conformer search, xTB-GOAT global optimization |
| **Imaginary-mode elimination** | 🔬 both | Iterative saddle escape: displaced single points on both sides of the imaginary mode, re-optimisation and frequency per round. ESD rates are refused on a structure still at a saddle |
| **ML potentials** | 🔬 both | 8 backends (ANI-2x, AIMNet2, MACE-OFF, CHGNet, M3GNet, SchNetPack, NequIP, ALIGNN) as ASE calculators, reachable as individual pipeline steps |
| **Crystal structure generation** | 🧱 solid state | Genarris integration for random organic crystal generation with configurable space groups, Z values and MPI-parallel execution |
| **CO₂ coordination** | 🧲 metal | Automated CO₂ placement around metal centres, relaxed distance (1.6–4.0 Å) and rotation (±180°) scans |
| **Reporting** | 🔬 both | DOCX report with embedded spectra and tables, generated at the end of every main run; JSON export and text summaries alongside |
| **Cluster execution** | 🔬 both | SLURM and SSH-remote backends, automatic resource detection, site profiles, node-local staging, queue monitoring and quota handling |
| **AI agent** | 🔬 both | Conversational coding agent on any **tool-capable** model (Claude / OpenAI / KIT Toolbox / local Ollama) — sandboxed bash, subagents, MCP, persistent memory, dashboard control |
| **Structure enumeration (ChemDarwin)** | 🔬 both | Dashboard tab that morphs a seed molecule with **user-supplied reaction SMARTS** over N iterations, with forbidden/protected substructure filters and ligand-wise mutation for metal complexes; draws the results and maps them in chemical space. Candidates are submitted to a DELFIN workflow by hand — there is no automated fitness loop |

### Ways to use it

- **Dashboard** — browser UI for job setup, submission, monitoring and result analysis. This is the usual way to run DELFIN.
- **Command line** — `delfin` runs a workflow from a `CONTROL.txt` in a working directory, plus 18 companion commands (see the [CLI reference](#-cli-reference)).
- **Python API** — `delfin.api` for notebooks and workflow engines.
- **Step and pipeline system** — `delfin-step` runs a single registered tool step, `delfin-pipeline` runs a declarative YAML pipeline with branching and parallel branches. This is how the non-ORCA engines and the ML potentials are reached.
- **MCP servers** — DELFIN ships three Model Context Protocol servers: `delfin-tools-server` (the tool platform), `delfin-docs-server` (literature and calculation search) and `delfin-ops-server` (typed runtime actions).
- **AI agent** — in the dashboard tab or as the standalone terminal program `delfin-agent`; it operates the dashboard, runs and interprets workflows, and edits code.

### Design principles

1. **One interface across the QM stack.** Each stage is routed to an established engine: RDKit / Open Babel / Architector / MANTA for structure, xTB / CREST / GOAT for conformer pre-screening, **ORCA for the DFT workflow**. Turbomole, the periodic codes and 30-odd further engines are reachable as individual pipeline steps, as are the ML potentials; the redox / OCCUPIER / ESD workflow itself is ORCA-based.
2. **Reproducible, documented runs.** Every run records its provenance — input SMILES, method choices, solver logs, spectra — and a DOCX report is generated at the end of each main run, ready for supplementary information.
3. **Geometry sanity checks.** Generated structures pass a default-on gate against severe covalent distortion and bad angles. A stricter topology gate (`DELFIN_TOPOLOGY_HARD_GATE=2`) additionally rejects detached metal–donor bonds, extra fragments, atom-count mismatches and atomic collisions; it is off by default. π-system planarity and hapticity are enforced constructively by the MANTA builder rather than checked after the fact.
4. **Structure variation on top of prediction.** The ChemDarwin tab enumerates structural variants from a seed molecule and a reaction SMARTS; their properties are then computed with the same DELFIN workflows.

### Use cases

**Catalysis & energy** 🧲
- Multi-step redox tuning of transition-metal complexes (up to 3 sequential oxidation and reduction steps per pipeline)
- Spin-state-resolved energetics with adaptive broken-symmetry configuration evolution
- Topology-aware structure generation across σ-donors, π-haptic modes (η²–η⁶) and mono- to multi-metallic complexes
- Born-Haber stability-constant (log K) cycles for individual complexes, for ligand-variant comparison

**Photophysics & emissive materials** 🧪 / 🧱
- TADF triage: ΔE(S-T) from an xTB/sTDA pre-screen, then ISC / RISC rates from DFT
- Phosphorescence design: per-sublevel lifetimes for Ir(III) / Pt(II) / Cu(I) and organic phosphors
- Excited-state-geometry-aware fluorescence and phosphorescence rates (IC, ISC, RISC; E₀₀ adiabatic energies)
- NLO chromophore design via β tensors; random organic crystal generation (Genarris)

**Pharmaceutical & medicinal chemistry** 🧪
- Conformer ensembles and free-energy ranking for drug-like scaffolds (CREST + xTB-GOAT + DFT)
- Redox-potential and spin-state prediction for metallodrug and metalloenzyme-mimic design
- Boltzmann-weighted ensemble NMR for stereo- and regioisomer assignment

**Spectroscopy & characterization** 🔬
- UV-Vis (TD-DFT) absorption, IR vibrational spectra, AFP overlay plots
- Boltzmann-weighted ensemble NMR (¹H / ¹³C shieldings and J-couplings)
- Hyperpolarizability β tensors and dipole moments; imaginary-frequency cleanup (`delfin --imag`)

---

## 🚀 Installation

**PyPI package**: https://pypi.org/project/delfin-complat/

### Requirements

- **Python 3.10 or 3.11**
- **ORCA 6.1.1** in your `PATH` (`orca` and `orca_pltvib`) — [free for academic use](https://orcaforum.kofo.mpg.de/app.php/portal)
- **Optional:** `crest` and `xtb` (CREST/xTB workflows)
- **Optional:** `censo`, `anmr`, `c2anmr`, `nmrplot` (ensemble NMR)
- **Optional:** `xtb4stda`, `stda`, `std2` plus the `xtb4stda` runtime files (TADF screening, xTB-based response)
- **Optional:** any of the 90+ supported programs — see [External programs](#-external-programs)
- **Optional (dashboard):** JupyterLab/Notebook or Voila

### Install methods

**Installer (any Linux workstation or cluster login node):**
```bash
git clone https://github.com/ComPlat/DELFIN.git ~/software/delfin
bash ~/software/delfin/install.sh              # DELFIN, ORCA wiring, QM/analysis tools, Ketcher
bash ~/software/delfin/install.sh --all        # everything, ML and AI stacks included (several GB)
bash ~/software/delfin/install.sh --only crest,gxtb
bash ~/software/delfin/install.sh --dry-run    # print the plan, change nothing
bash ~/software/delfin/install.sh --update     # pull DELFIN, update it and every installed tool
bash ~/software/delfin/install.sh --repair     # check everything, fix what is broken
```

It needs no root and no module system: a Python 3.10/3.11 is taken from the machine or fetched with micromamba, ORCA is found where it was unpacked (or pass `--orca DIR|TARBALL`; ORCA is licensed and never downloaded), and OpenMPI — the version ORCA's directory name gives, 4.1.8 when there is no ORCA yet, configured as ORCA needs it — is found or built. Nothing in it is specific to one site. The script itself lives in `delfin/installers/install_delfin.sh`.

Whatever is left out can be added later, three ways that share one list of tools (`delfin/installer.py`): `install.sh --only …`, the install buttons in the dashboard's Settings tab, or automatically the moment a calculation needs a tool (`DELFIN_AUTO_INSTALL_QM_TOOLS=0` switches that off). `python -m delfin.installer --list` shows what can be installed, `--status` what is.

**PyPI install:**
```bash
pip install delfin-complat
```

**Development install (from source):**
```bash
git clone https://github.com/ComPlat/DELFIN.git
cd DELFIN
python -m venv .venv
source .venv/bin/activate
pip install -e ".[agent,docs,dev]"
```

All Python dependencies (RDKit/Open Babel for SMILES workflows, ipywidgets/py3Dmol for dashboard visualisation, python-docx for reports) are installed automatically. This exposes the console command **`delfin`** and enables `python -m delfin`.

External QM binaries are not installed by `pip`. For the local binary-based setup of `xtb`, `crest`, `xtb4stda`, `stda`, `std2` and the `xtb4stda` runtime bundle, see `delfin/qm_tools/README.txt`.

### External QM tool setup

The installer above does this. By hand, after installing the Python package:

```bash
source delfin/qm_tools/env.sh
USE_SYSTEM_TOOLS=1 bash delfin/qm_tools/install_qm_tools.sh
bash delfin/qm_tools/check_qm_tools.sh
```

To check an installation at any time:

```bash
delfin doctor              # ORCA, xTB, OpenMPI, scratch, SLURM, doc index
delfin qm_check            # how each QM binary is resolved
delfin mlp_check           # ML-potential backends, PyTorch, CUDA
delfin analysis_check      # Multiwfn / CENSO / ANMR / morfeus
delfin csp_check           # Genarris
```

---

## 🏁 Quick start

### The dashboard

The usual way to run DELFIN. Start it as a standalone web app:

```bash
delfin-voila                 # starts on 127.0.0.1:8866
```

Voila prints a URL containing an access token — token authentication is mandatory and cannot be switched off. Open that URL, go to **Submit Job**, paste a SMILES or load an XYZ, choose the workflow settings, and submit. **Job Status** follows the queue; **Calculations** browses the results.

### The command line

A working directory needs two files:

* `CONTROL.txt` — the configuration (create one with `delfin --define`)
* `input.txt` — either an XYZ body (without the two header lines) or a single-line SMILES string

Starting from a full `.xyz` file works too: `delfin --define=input.xyz` converts it.

```bash
delfin                       # run in the current directory
python -m delfin             # same thing
delfin /path/to/project      # run in another workspace directory
```

Results are written into the workspace: `DELFIN.txt` (redox potentials), `OCCUPIER.txt` (occupation tracking), `ESD.txt` (excited-state results), `DELFIN_Data.json`, the DOCX report, and one folder per step.

### Python

```python
from delfin.api import run

run(control_file="CONTROL.txt")
```

```python
from delfin.dashboard import create_dashboard

ctx = create_dashboard(backend="auto")   # SLURM when available, else local
```

---

## 🧪 Dashboard

| Tab | Purpose |
|-----|---------|
| **Submit Job** | SMILES/XYZ input, 3D preview and editor, converter choice (`QUICK`, `NORMAL`, `MANTA`, `ARCHITECTOR`), `BUILD COMPLEX`, `MANTA`, `SUBMIT FUKUI`, `SUBMIT ONLY GOAT`, `SUBMIT DELFIN + CO2`, `VALIDATE CONTROL` |
| **Recalc** | Edit an existing `CONTROL.txt` and resubmit; a smart recalc recomputes only what the edit changed |
| **ORCA Builder** | Interactive ORCA input generation with geometry preview |
| **TURBOMOLE Builder** | Turbomole define workflow (SLURM backends) |
| **ChemDarwin** | Reaction-SMARTS structure enumeration and chemical-space map (off by default) |
| **DELFIN Agent** | The AI agent: subagents, MCP, sandboxed bash, dashboard control, result analysis, persistent memory |
| **Agent Activity** | Running and recent agent and subagent work |
| **Job Status** | Queue monitoring (local/SLURM), resource usage, fairshare and limits, job cancellation |
| **Calculations** | File browser, search, recalculation trigger, energy statistics, browser-launched workflows such as `Calc NMR` and `Calc CENSO/ANMR` |
| **Archive** | Archive browser with statistics |
| **Remote Archive** | Browse and transfer results on a remote machine over SSH |
| **Office** | Documents and spreadsheets, with an agent that has no chemistry tools |
| **Literature** | The indexed literature corpus |
| **Tools** | The registered tool steps and their parameters |
| **Ketcher** | 2D structure drawing |
| **Reactions** | Reaction graph |
| **Pipelines** | Declarative step pipelines |
| **Settings** | Tool detection, per-tool Install/Update buttons, runtime configuration, agent model settings |

Some tabs are hidden by default and are enabled in Settings.

### Running the dashboard

```bash
delfin-voila                 # 127.0.0.1:8866, prints a URL with an access token
delfin-voila --port 9000     # custom port
delfin-voila --dark          # dark theme
delfin-voila --keep          # run inside a tmux session, survives a dropped terminal
delfin-voila --resume SID    # reopen a previous agent session

# On HPC/login nodes keep the default 127.0.0.1 bind and use an SSH tunnel.
# A direct network bind needs an explicit override:
delfin-voila --ip 0.0.0.0 --allow-remote-bind
```

Detailed documentation: [docs/SETTINGS_AND_SETUP.md](docs/SETTINGS_AND_SETUP.md)

---

## ✨ Workflows

### Redox potentials (OCCUPIER)

DELFIN's core workflow automates spin-state identification and redox-potential calculation:

- **Multi-step redox**: up to 3 sequential oxidation and reduction steps
- **Parallel workflows**: oxidation and reduction run simultaneously with automatic PAL splitting (PAL=12 → 6+6)
- **Adaptive broken symmetry**: broken-symmetry configurations evolve along the OCCUPIER tree
- **Smart recalc**: a SHA-256 fingerprint of each input and its dependencies skips jobs an edit did not change

### Stability constants (log K)

A Born-Haber-type thermodynamic cycle for one coordination complex per run:

- **Automatic reaction analysis**: unique ligands, denticities, displaced solvent count and the matching metal-solvent reference complex are extracted from the complex SMILES
- **OCCUPIER-aware metal treatment**: the target complex and the solvated metal reference can use the same converter, pre-optimisation, multiplicity and broken-symmetry logic as the main workflow
- **Efficient sub-workflows**: duplicate ligands are computed once; ligand, solvent and metal-solvent jobs run in parallel
- **Output**: OPT+FREQ free energies for complex, solvated metal, ligand(s) and solvent, combined into ΔG and log K

The CONTROL keys are `stability_constant`, `stability_constant_mode` and `stability_reaction` (the older spellings `thermodynamics*` are still read).

### Excited-state dynamics (ESD)

Photophysical properties from a list of states and transitions that **you** specify — DELFIN orchestrates the jobs, it does not choose the states for you:

- **Electronic states**: S0, S1, S2, T1, T2 geometry optimisation
- **ISC / RISC rates**: intersystem and reverse intersystem crossing
- **Internal conversion (IC)**: non-radiative decay rates
- **Fluorescence and phosphorescence**: radiative lifetimes, including per-sublevel phosphorescence
- **E₀₀**: adiabatic 0-0 transition energies
- **ΔE(S-T)**: singlet-triplet gaps

```ini
ESD_modul=yes
states=[S1,T1,S2,T2]
ISCs=[S1>T1,T1>S1]
ICs=[S1>S0]
emission_rates=[f,p]
```

There is no `delfin ESD` subcommand: a normal `delfin` run performs the ESD work when `ESD_modul=yes`. Rates are refused on a structure that still sits at a saddle point — see *imaginary-mode elimination* below.

### TADF screening (xTB-based)

Fast triage of TADF candidates with semi-empirical methods (needs `xtb`, `xtb4stda` with its runtime files, and `stda`/`std2`):

- S0/T1 optimisation via xTB, S1 estimated from the Stokes shift
- ΔE(S-T), vertical and as a relaxed estimate
- First allowed and brightest singlet, with oscillator strengths
- Estimated absorption and emission wavelengths, T1 adiabatic and vertical

### Spectroscopy and properties

| Module | Output |
|--------|--------|
| **UV-Vis** | TD-DFT absorption spectra, oscillator strengths, transition analysis |
| **IR** | Vibrational frequencies, intensities, Lorentzian broadening, transmittance |
| **AFP** | Combined absorption/fluorescence/phosphorescence plot (`--afp-fwhm`, default 50 nm) |
| **Hyperpolarizability** | Static and frequency-dependent β tensors (β_HRS, β_VEC, DR), dipole moments |
| **Fukui** | Atomic Fukui indices from three ORCA single points (N, N±1) at fixed geometry, plus ESP reporting |

### Ensemble NMR via CREST + CENSO + ANMR

A browser-launched ensemble NMR workflow for `.xyz` inputs, started from the **Calculations** tab:

- **`Calc NMR`** — single-structure ORCA NMR workflow
- **`Calc CENSO/ANMR`** — end-to-end ensemble workflow: CREST conformer sampling with `-nmr` → CENSO re-ranking and refinement → ORCA shieldings and couplings for the survivors → Boltzmann-weighted ANMR spectra → `PNG`, `PDF`, `SVG` and JSON/text summaries

The four helper programs (`censo`, `anmr`, `c2anmr`, `nmrplot`) are external; they are detected in DELFIN's runtime layer and can be installed from the dashboard's Settings tab, or automatically on first use.

### Imaginary-mode elimination (IMAG)

A structure with an imaginary frequency is not a minimum, and rate constants computed on it are meaningless. `delfin --imag` (or `IMAG=yes` in CONTROL) walks out of the saddle iteratively: it displaces the geometry along the imaginary mode in both directions, takes single points, re-optimises and recomputes the frequencies, round after round, until the mode is gone or the round cap is reached. Modes below `allow_imaginary_freq` count as numerical noise. The ESD module refuses to compute rates on a structure that is still at a saddle.

### ML potentials and the calculator factory

34 computational backends behind one interface, each returning a standard ASE `Calculator`:

```python
from delfin.calculators import create_calculator

calc = create_calculator("ani2x", device="cuda")          # ML potential
calc = create_calculator("orca", method="B3LYP")          # DFT
calc = create_calculator("xtb")                           # Semi-empirical
calc = create_calculator("vasp", xc="PBE", kpts=[4,4,4])  # Periodic DFT

atoms.calc = calc
energy = atoms.get_potential_energy()
```

All backends are lazily loaded — nothing is imported until it is used — and the ML backends validate CUDA with a CPU fallback. Of the 8 ML potentials, six ship usable pre-trained weights; **SchNetPack and NequIP need a trained model of your own** (`model_path`). These calculators and the non-ORCA engines are used through the step and pipeline system (`delfin-step`, `delfin-pipeline`); the `CONTROL.txt` workflow itself runs on ORCA.

### Crystal structure generation (Genarris)

Random organic crystal generation with configurable space groups, Z values and MPI-parallel execution. This is generation, not polymorph ranking: no lattice energies are computed and no structures are scored. It needs the external `gnrs` binary and `mpirun`, and is driven from a Genarris configuration file through the step adapter.

### CO₂ coordinator

Automated CO₂ placement around metal centres with relaxed distance scans (1.6–4.0 Å) and rotation scans (±180°). Run with `delfin co2 …` or from the Submit tab.

### Reporting and export

| Format | Command | Content |
|--------|---------|---------|
| **Text** | automatic | `DELFIN.txt` (redox potentials), `OCCUPIER.txt` (occupation tracking), `ESD.txt` (excited-state results) |
| **DOCX** | automatic, or `delfin --report docx` | Word document with embedded spectra plots and structured tables |
| **JSON** | automatic, or `delfin --json` | Machine-readable `DELFIN_Data.json` (`--json-output` for another path) |
| **UV-Vis** | `delfin_ESD output.out` | Parsed UV-Vis spectrum with transitions and oscillator strengths |
| **IR** | `delfin_IR output.out` | Parsed IR spectrum with vibrational modes and intensities |
| **NMR** | `delfin_NMR output.out` | ¹H NMR spectrum plot and report from an ORCA NMR output |
| **AFP** | `delfin --afp` | Absorption/fluorescence/phosphorescence combined plot |

The DOCX and JSON reports are written at the end of every main run; a failure there is logged as a warning and does not fail the run.

---

## 🧬 Structure generation from SMILES

Which builder runs is decided by one CONTROL key, and a SMILES run must set it:

```ini
smiles_converter=[QUICK|NORMAL|MANTA|ARCHITECTOR]
```

| Value | What it does | Best for |
|-------|--------------|----------|
| `QUICK` | One ETKDGv3 embedding; for metals it tries stk → RDKit → unsanitised → Open Babel and returns on the first success. One structure. | Quick previews |
| `NORMAL` | RDKit ETKDG with multi-seed embedding and optional UFF refinement (RDKit / Open Babel). One structure. **This is the default when nothing is set.** | Organic molecules |
| `MANTA` | DELFIN's own coordination builder: the enumerated coordination-isomer × conformer manifold, then screening, optimisation and refinement down to one geometry. | Metal complexes |
| `ARCHITECTOR` | [architector](https://github.com/lanl/Architector) automated 3D generation; needs a metal-containing SMILES. | Metal complexes, instant preview |

In the Submit tab the same choice appears as the converter dropdown, next to `BUILD COMPLEX` (stepwise assembly using ORCA's `%DOCKER`, submitted as a job) and the `CONVERT SMILES`, `QUICK CONVERT SMILES` and `CONVERT SMILES + UFF` buttons of the structure editor.

### MANTA

MANTA is DELFIN's force-field-free coordination-structure engine — *construct, don't search*. From a transition-metal SMILES it enumerates the coordination isomers by Burnside–Pólya counting over the donor set and the coordination polyhedron, seats each isomer on an ideal polyhedron with metal–donor distances from covalent radii, and expands each into conformers. **No force field touches the metal**: UFF and MMFF have no transition-metal parameters, which is what systematically distorts M–D lengths and L–M–L angles in conventional builders.

The manifold then goes through a funnel that `CONTROL.txt` steers:

```
build manifold → screen (one single point per frame and multiplicity)
               → optimise (xTB per survivor)
               → refine (GOAT or CREST on the best)
               → winner written to start.txt
```

The ranked unit is a (frame, multiplicity) pair, not a frame. Roughly 26 `MANTA_*` keys control quality, gates, screening method, optimisation, refinement and time budget — `delfin --define` writes them all with their defaults. `delfin-manta` runs the builder on its own.

**What MANTA is for:** starting geometries, not production geometries. The output is a constructive geometry of roughly force-field quality — correct in topology and coordination, not xTB- or DFT-accurate.

**GUPPY and MANTA.** GUPPY was the earlier name for MANTA plus the energy-rank-and-pick driver around it. `smiles_converter=MANTA` is the current spelling; `GUPPY` and the older bare `GUPPY=yes` are still read and mean the same thing, and every `MANTA_*` key falls back to its `GUPPY_*` predecessor, so CONTROL files written before the rename keep working. The driver module, its working directory and the `delfin-guppy` / `delfin-guppy-batch` commands still carry the GUPPY name.

### Conformer search and sampling

- **CREST** — conformer-rotamer ensemble generation and sorting
- **xTB-GOAT** — gradient-based global optimisation
- **MANTA sampling** — multi-start sampling with parallel xTB optimisation and energy-ranked trajectories (`delfin-guppy`, `delfin-guppy-batch` for a batch of SMILES)

---

## ⚙️ Configuration (CONTROL.txt)

DELFIN is configured by `CONTROL.txt` in the working directory. `delfin --define` writes a fully commented template with every key and its default; the groups below are the ones you touch most.

### Workflow control
* `method = OCCUPIER | classic | manually` (leave empty for ESD-only runs)
* `OCCUPIER_method = auto | manually` (auto uses adaptive tree-based sequences)
* `calc_initial = yes | no`
* `oxidation_steps = 1,2,3` / `reduction_steps = 1,2,3`
* `parallel_workflows = yes | no | auto`

### Structure
* `smiles_converter = QUICK | NORMAL | MANTA | ARCHITECTOR` (required for a SMILES run)
* `MANTA_*` — 26 keys steering the MANTA builder and its funnel
* `XTB_preOPT = yes | no` / `global_optimizer = GOAT | CREST`
  * the older spellings `XTB_OPT` / `XTB_GOAT` / `CREST` are still accepted
* `XTB_SOLVATOR = yes | no`

### Resource management
* `PAL = N` / `maxcore = MB`
* `pal_jobs = N` (parallel PAL processes; auto-detected on a cluster when unset)
* `orca_parallel_strategy = auto | threads | serial`

### Optional modules
* `ESD_modul = yes | no` with `states`, `ISCs`, `ICs`, `emission_rates`
* `IMAG = yes | no`
  * `IMAG_scope = all | initial` (default `all`: every redox step too)
  * `allow_imaginary_freq = -50` (cm⁻¹, ≤ 0: smaller imaginary modes count as numerical noise)
  * `IMAG_max_rounds = 2` (each round is one frequency calculation)
* `stability_constant`, `co2_coordination`, `calc_prop_of_interest`, `reorganisation_energy`
* `TDDFT_*`, `deltaSCF_*`, `elprop_*`, `tadf_xTB*`, `hyperpol_xTB*`

### Error recovery
* `enable_auto_recovery = yes | no` (ORCA error recovery, continuing from the last `.gbw` and geometry)
* `max_recovery_attempts = N` (default: 3)
* `enable_job_timeouts = yes | no` (`no` for unlimited runtime)

---

## 🤖 AI agent

DELFIN ships a conversational AI agent: it reads and edits code, runs sandboxed shell commands, drives the dashboard, researches methods, analyses results and delegates self-contained work to subagents. It runs in the **DELFIN Agent** dashboard tab and as the standalone terminal program `delfin-agent`. It is a single, direct agent — there is no fixed review pipeline; when work benefits from extra hands it spins up subagents on demand.

### Modes

A mode decides which role prompt and which surface the agent gets:

| Mode | Where | For |
|------|-------|-----|
| `dashboard` (default) | dashboard, terminal | Operating DELFIN: settings, submissions, results |
| `solo`, shown as **Code** | dashboard, terminal | Reading and writing code |
| `office` | dashboard, terminal | Documents and spreadsheets, without the chemistry tools |
| `research` | terminal only | Literature and method research |

`plan` is **not** a mode but a permission profile — see below.

### Works with any tool-capable model

The agent is provider-agnostic. A model that cannot call tools is refused with an explanation, because every non-trivial action is a tool call.

| Provider | Primary backend | Notes |
|----------|------------------|-------|
| **Claude** (default) | API via `ANTHROPIC_API_KEY` | [Claude Code CLI](https://claude.ai/code) when `claude` is on PATH |
| **OpenAI / Codex** | API via `OPENAI_API_KEY` | [Codex CLI](https://github.com/openai/codex) when `codex` is on PATH |
| **KIT Toolbox** | API via `KIT_TOOLBOX_API_KEY` | university-hosted, OpenAI-compatible (vLLM) |
| **Ollama / vLLM / LM Studio** (local & open source) | OpenAI-compatible endpoint via `OLLAMA_HOST` (default `http://localhost:11434`) | **no API key, runs offline.** Any pulled model, e.g. `ollama pull qwen2.5-coder`. Tool-capable models (`qwen2.5-coder`, `qwen3-coder`, `llama3.3`) are the ones to pick. |

A model-capability layer detects each model's real context window and its tool, vision and reasoning support — live from Ollama's `/api/show` and from OpenAI-compatible `/v1/models`, with a curated fallback table — and sizes context management to it. For Ollama it sends the correct `num_ctx`, so local models are not silently capped at 2–4k; it skips parameters local servers reject and strips `<think>` leaks. Weak or small models automatically get a slimmer prompt and a reduced tool surface.

Available providers are auto-detected from environment variables and PATH. DELFIN treats each LLM the way it treats ORCA, xTB or CREST: a dependency that you install and authenticate yourself. DELFIN bundles no LLM and redistributes none; it calls each provider's official API, OpenAI-compatible endpoint or CLI binary. Hosted usage runs under your own account and billing; local models run on your own hardware at no per-token cost.

### Subagents

For parallel research, read-only audits or planning that must not edit, the agent delegates to isolated subagents — each with its own tool loop and usually tighter permissions:

- `explore` — read-only investigation, reports findings
- `plan` — step-by-step plan, makes no edits
- `code-reviewer` — independent read-only review
- `general-purpose` — self-contained task, inherits permissions

Your own presets can be added as markdown files in `~/.delfin/subagents/`. Subagents run in parallel (a writer preset gets its own git worktree, so concurrent edits cannot clobber) or in the background; a finished one can be continued with its context intact. Per-run limits — 40 tool calls, 900 s wall clock, 16000 output tokens by default — are configurable in Settings, and a subagent may not spawn subagents of its own unless that depth is raised. Launch one with `/explore`, `/review`, `/plan` or `/delegate <task>`.

### Tools, MCP and safety

- **Built-in tools** (73): read/edit/write files, grep, sandboxed bash with long-running background jobs, code navigation, test runner, notebooks, documents and spreadsheets, web search and fetch, task tracking, scheduling, git worktrees, delegation, plus DELFIN-specific calculation and manual search.
- **MCP (Model Context Protocol)**: connect external MCP servers over stdio or Streamable HTTP — their tools, resources and prompts join the agent's surface (configured in `~/.delfin/mcp_servers.json`). A tool reached through MCP runs in a process DELFIN did not launch a command line for, so the shell sandbox is not around it: give a stdio server `"roots": [...]` (read-write) or `"read_roots": [...]` (read-only) and it starts inside a namespace holding only those paths. Servers that declare neither run uncontained, and the startup banner, `/mcp` and `delfin-agent doctor` all say so. DELFIN's own servers can take their roots from your settings instead — `agent.mcp_isolation: "builtin"` binds the calculations, office, workspace and state folders read-write and the archive and runtime trees read-only; it is off by default because the roots are inferred.
- **Sandboxed execution**: every shell command goes through a layered defence — allow-list, then bubblewrap or firejail — with credential directories (`~/.ssh`, `~/.aws`, `~/.gnupg`, …) masked, network denied by default, and every command appended to `~/.cache/delfin/agent-audit.jsonl`. Configurable with `DELFIN_AGENT_SANDBOX={auto,bwrap,firejail,allowlist,off}`.
- **Permission modes**: `plan` (read-only) · `default` (destructive actions ask) · `diff_approval` (writes stage a diff for `/approve`) · `acceptEdits` (writes allowed, bash still asks) · `bypassPermissions` (no prompts; the sandbox and the deny lists still apply). Per-pattern allow-list rules can be remembered across sessions. Files that define the agent's own permissions always require confirmation, in every mode.
- **Workspace trust**: settings, hook commands and MCP servers that come from a checked-out repository are ignored until you trust that directory (`/trust`).

### Slash commands

Built-in commands cover the session (`/help`, `/status`, `/cost`, `/context`, `/compact`, `/export`, `/undo`), the setup (`/model`, `/mode`, `/permissions`, `/effort`, `/doctor`, `/mcp`, `/trust`, `/tools`, `/agents`, `/skills`), the history (`/session`, `/rewind`, `/tasks`, `/memories`, `/plans`), pending changes (`/pending`, `/approve`, `/reject`) and the workspace (`/git`, `/bash`, `/jobs`). The dashboard adds direct control of DELFIN:

```
/control key <key> <value>   — change a single CONTROL key
/orca set <param> <value>    — configure the ORCA Builder
/calc ls, /calc read, /calc info, /calc tree — browse calculations
/analyze energy, /analyze convergence, /analyze errors — analyse ORCA outputs
/recalc check-all, /recalc auto — smart recalculation
/submit, /cancel — job management (with confirmation)
```

Beyond those, a markdown file in `~/.delfin/commands/` or `<workspace>/.delfin/commands/` becomes a slash command of its own, and a skill (`~/.delfin/skills/<name>/SKILL.md`) can be invoked either by you or by the model. Fourteen skills ship with DELFIN, among them `diagnose-failed-run`, `recalc-failed`, `tddft-excited-states`, `freq-thermochemistry`, `solvation-setup` and `tune-control`.

### Memory and sessions

- **Persistent memory**: facts and preferences survive across sessions and are recalled by relevance (`/remember`, `/memories`, `/forget`). Memory is scoped by domain, so a preference recorded while coding does not surface in an office session.
- **Project instructions**: a `DELFIN.md` or `AGENTS.md` in the project is loaded automatically.
- **Long sessions**: token-aware compaction keeps long conversations coherent; `/compact` and `/context` expose the state.
- **Agent workspace**: `~/agent_workspace/` for uploaded reference files.
- **Cost tracking**: per-session token usage and cost in real time (local models are free).
- **Session persistence**: conversations are saved, restored, forked and exported as Markdown.

---

## 📋 CLI reference

### Running a workflow

- `delfin` — run the workflow in the current directory (needs `CONTROL.txt` and `input.txt`)
- `delfin /path/to/project` — run it in another directory
- `python -m delfin` — the same entry point

### Setup and configuration

- `delfin --define[=input.xyz] [--overwrite]` (short: `-D`) — create or update `CONTROL.txt`, optionally converting an XYZ to `input.txt`
- `delfin /path/to/project --define[=input.xyz]` — the same, in another workspace
- `delfin --control /path/to/CONTROL.txt` — run with a CONTROL file from elsewhere; it sets the workspace root
- `delfin --version` (short: `-V`)

### Diagnostics

- `delfin doctor [--json] [--scratch DIR]` — installation self-check: ORCA, xTB, OpenMPI, scratch directory, SLURM, API keys, documentation index. Runs no calculation and needs no network.
- `delfin qm_check [TOOLS…]` — how DELFIN resolves each QM binary
- `delfin qm_run TOOL [--cwd DIR] [--capture] -- ARGS` — run one QM tool through DELFIN's resolver
- `delfin mlp_check` — ML-potential backends, PyTorch and CUDA
- `delfin analysis_check` — Multiwfn, CENSO, ANMR, morfeus
- `delfin csp_check` — Genarris availability

### Execution control

- `delfin --recalc` — re-parse existing results and restart only what is missing or incomplete
- `delfin WORKSPACE --recalc --occupier-override STAGE=INDEX` — force an OCCUPIER index for a stage during a recalc
- `delfin --report` — recompute redox potentials from existing outputs (bare: text; `--report docx` for Word)
- `delfin --imag` — eliminate imaginary modes from existing ORCA results and regenerate the summary
- `delfin run_orca [file.inp] [-i FILE] [-o FILE]` — run a single ORCA input through DELFIN's runner
- `delfin stop --workspace PATH [--signal INT|TERM|KILL] [--dry-run] [--cleanup] [--wait-seconds S]` — stop running DELFIN processes for a workspace

### Cleanup and maintenance

- `delfin --no-cleanup` — keep temporary files and scratch folders after the run
- `delfin --cleanup` (short: `-C`) — remove intermediates and exit
- `delfin cleanup [--dry-run] [--workspace PATH] [--scratch PATH]` — fine control over workspace and scratch cleanup
- `delfin cleanup --orca` — stop running ORCA jobs and purge OCCUPIER scratch folders
- `delfin --purge` — remove DELFIN-generated artifacts after confirmation

### Export and reporting

- `delfin --json [--json-output FILE]` — collect results into `DELFIN_Data.json`
- `delfin --report docx` — Word report with embedded plots
- `delfin --afp [--afp-fwhm NM]` — AFP spectrum plot from existing ESD results (default FWHM 50 nm)

### Specialised workflows

- `delfin co2 [--define] [--recalc] [--charge N] [--multiplicity M] [--solvent S] [--metal M] [--broken_sym B]` — CO₂ coordinator
- `delfin tadf_xtb …` — xTB/sTDA TADF screening batch runner (`--smiles`, `--xyz-file`, `--crest`, `--goat`, `--t1-opt`, `--pal`, `--parallel-jobs`)
- `delfin hyperpol` / `delfin hyperpol_xtb …` — xTB-based hyperpolarizability batch runner

Excited-state dynamics has no subcommand of its own: set `ESD_modul=yes` in `CONTROL.txt` and run `delfin`.

### Companion commands

| Command | What it does |
|---------|--------------|
| `delfin-voila` | Launch the dashboard as a standalone web app |
| `delfin-agent` | The AI agent in the terminal (`chat`, `run`, `init`, `doctor`, `session`, `approvals`, `scheduler`, …) |
| `delfin-manta` | MANTA: the coordination-isomer × conformer manifold from a metal SMILES |
| `delfin-guppy [input.txt]` | Multi-start sampling with repeated xTB optimisation and ranked trajectories (`--runs`, `--parallel-jobs`, `--pal`, `--refine`, `--screen`, `--max-isomers`, …) |
| `delfin-guppy-batch` | The same over every SMILES in a batch file (`--row` for SLURM arrays) |
| `delfin-build [input.txt]` | Build metal complexes stepwise from SMILES using ORCA's `%DOCKER` (`--goat`, `--no-ligand-goat`, `--dry-run`, …) |
| `delfin-fukui` | Atomic Fukui indices from three ORCA single points |
| `delfin-step` | Run one registered tool step (`--list`, `--describe`, `--slurm`) |
| `delfin-pipeline` | Run a declarative YAML pipeline (`--param`, `--scheduled`, `--dry-run`) |
| `delfin-app` | Application registry: `list`, `template`, `run <keyfile>`, `describe` |
| `delfin-json <project_dir>` | Collect project outputs into JSON (`-o FILE`) |
| `delfin_ESD output.out` | UV-Vis spectrum report from an ORCA output |
| `delfin_IR output.out` | IR spectrum report from an ORCA frequency output |
| `delfin_NMR output.out` | ¹H NMR spectrum report from an ORCA NMR output |
| `delfin-docs-index` | Build the documentation search index from `literature/` |
| `delfin-docs-server` · `delfin-ops-server` · `delfin-tools-server` | The three MCP servers |

---

## 🏗️ Architecture

### Package layout

```
delfin/
  cli.py                   # command-line entry point and subcommand dispatch
  cli_helpers.py           # the top-level argument parser
  api.py                   # programmatic API for notebooks and workflow engines
  config.py                # CONTROL.txt parsing, aliases and defaults
  define.py                # the CONTROL.txt template (every key and its default)
  calculators.py           # unified ASE calculator factory (34 backends)

  # ── Orchestration ──
  workflows/               # where the workflow logic lives
    pipeline.py            # high-level orchestration (classic / manually / OCCUPIER)
    engine/                # classic.py, occupier.py, scheduler.py
    scheduling/            # manager.py (global job manager), pool.py, priority.py
    contrib/               # co2, esd, hyperpol, imag, occupier, tadf_xtb workflows
    registry.py, types.py

  # ── Core workflows ──
  occupier.py              # OCCUPIER sequence execution and summary
  occupier_auto.py         # adaptive sequence rules and populated-state branching
  bs_evolution.py          # broken-symmetry evolution along the OCCUPIER tree
  esd_module.py            # excited-state dynamics (ISC / IC / fluorescence / phosphorescence)
  esd_input_generator.py   # ORCA input builders for ESD states
  stability_constant.py    # Born-Haber cycle, log K
  tadf_xtb.py              # TADF screening via xTB/sTDA
  hyperpol.py              # hyperpolarizability (NLO)
  fukui.py, cli_fukui.py   # atomic Fukui indices
  xtb_crest.py             # xTB / GOAT / CREST / ALPB solvation workflows
  imag.py, cli_imag.py     # imaginary-frequency elimination
  co2/                     # CO2 coordinator

  # ── Structure generation ──
  manta/                   # the coordination-structure engine (~90 modules)
  cli_manta.py             # the delfin-manta CLI
  smiles_converter.py      # SMILES→3D entry points, orchestrates delfin.manta
  guppy_sampling.py        # multi-start sampling, ranking and refinement funnel
  build_up_complex.py      # stepwise metal-complex assembly (ORCA %DOCKER)
  class_modules/           # chemistry-domain SMILES→XYZ specialists
  topology_hard_gate.py    # optional SMILES-vs-XYZ topology gate

  # ── ORCA and recovery ──
  orca.py                  # running ORCA, retries
  orca_recovery.py         # error detection and input repair
  common/orca_input.py     # reading and writing ORCA inputs

  # ── Recalculation ──
  smart_recalc.py          # fingerprints: what may be skipped
  recalc_control.py        # what a CONTROL edit invalidates

  # ── Resource management and clusters ──
  cluster_utils.py         # SLURM / PBS / LSF resource detection
  slurm_submit.py, scheduler_profiles.py, submit_templates/
  ssh_transfer_jobs.py, remote_archive.py, quota.py
  runtime_setup.py         # auto-detection of external programs

  # ── Tool integrations ──
  tools/                   # step and pipeline system, 25 adapters, mcp_server.py
  mlp_tools/               # ML potentials (ANI-2x, MACE, CHGNet, M3GNet, …)
  ai_tools/                # AI/ML tool registry (21 tools across 9 categories)
  analysis_tools/          # cclib, Packmol, Multiwfn, CENSO, morfeus wrappers
  csp_tools/               # crystal structure generation (Genarris)
  qm_tools/                # external QM binary management

  # ── Interfaces ──
  dashboard/               # the browser dashboard (Voila / JupyterLab)
  agent/                   # the AI agent (engine, clients, tools, sandbox, memory)
  doc_server/              # MCP documentation server
  ops_server/              # MCP operations server
  installers/              # install_delfin.sh
  reporting/               # DOCX, JSON and text reports
  common/                  # logging, paths, ORCA block assembly
```

`delfin/pipeline.py`, `global_manager.py`, `global_scheduler.py`, `dynamic_pool.py`, `job_priority.py` and `parallel_*.py` still exist as backward-compatible shims that re-export from `delfin/workflows/`.

### Global resource management

DELFIN coordinates all computational work through a global job-manager singleton:

* **Single source of truth** — PAL is read once from `CONTROL.txt` and managed centrally
* **Automatic PAL splitting** — parallel oxidation and reduction workflows share cores (PAL=12 → 6+6)
* **Thread-safe execution** — a shared resource pool prevents race conditions
* **Subprocess coordination** — OCCUPIER subprocesses inherit their limits through the environment

### Tool integration pattern

1. **Lazy loading** — availability is checked with `importlib.util.find_spec()`; nothing is imported until used
2. **Per-tool install** — individual Install/Update buttons in the dashboard's Settings tab
3. **Auto-detection** — `shutil.which()` for binaries, compatible with cluster module systems
4. **Unified interface** — the calculator backends all return standard ASE `Calculator` objects

### Automatic error recovery

ORCA failures are detected and repaired automatically, continuing from the last `.gbw` and geometry:

```
ORCA fails → detect the error type → modify the input → continue from the last .gbw and xyz → retry
```

| Error | Automatic fix |
|-------|---------------|
| **SCF not converged** | SlowConv → VerySlowConv + KDIIS |
| **TRAH segfault** | NoAutoTRAH |
| **Geometry not converged** | Smaller trust radius → looser criteria |
| **MPI crashes** | Fewer cores |
| **Memory errors** | Lower maxcore and PAL |
| **Transient system errors** | Exponential backoff |

Enable it with `enable_auto_recovery=yes`; `max_recovery_attempts` (default 3) bounds the attempts per error type. See **[docs/RETRY_LOGIC.md](docs/RETRY_LOGIC.md)** for the complete guide.

---

## 🔌 External programs

DELFIN detects 90+ supported programs: about 55 external binaries through `PATH`, `$MODULEPATH` and HPC fallback directories (so `module load gaussian/16` is enough), and about 35 Python packages by import probe (ML potentials, AI/ML models, wrapper libraries). Programs that cannot be pip-installed — ORCA, Gaussian, VASP, TURBOMOLE, Multiwfn and other licensed or externally managed software — are detected and reported, never installed.

Install and update buttons for the pip-installable integrations are in the dashboard under `Settings → Tool Installation`.

Of the AI/ML tools listed below, DELFIN detects and installs all of them, and **architector** is additionally wired into structure generation; the others are made available for you to use, not called by a DELFIN workflow.

<details>
<summary><b>Linked overview of supported tools</b></summary>

**ML Potentials (8):**
[ANI-2x](https://github.com/aiqm/torchani),
[AIMNet2](https://github.com/isayevlab/AIMNet2),
[MACE-OFF](https://github.com/ACEsuit/mace),
[CHGNet](https://github.com/CederGroupHub/chgnet),
[M3GNet/MatGL](https://github.com/materialyzeai/matgl),
[SchNetPack](https://github.com/atomistic-machine-learning/schnetpack),
[NequIP/Allegro](https://github.com/mir-group/nequip),
[ALIGNN](https://github.com/usnistgov/alignn)

**QM Programs — Ab initio / DFT (11):**
[ORCA](https://orcaforum.kofo.mpg.de/app.php/portal),
[Gaussian (g16/g09)](https://gaussian.com/),
[TURBOMOLE](https://www.turbomole.org/),
[NWChem](https://www.nwchem-sw.org/),
[Q-Chem](https://www.q-chem.com/),
[GAMESS](https://www.msg.chem.iastate.edu/GAMESS/),
[Molpro](https://www.molpro.net/),
[Dalton](https://daltonprogram.org/),
[Psi4](https://psicode.org/),
[CFOUR](https://cfour.uni-mainz.de/),
[MRCC](https://www.mrcc.hu/)

**QM Programs — Periodic / Solid State (11):**
[VASP](https://www.vasp.at/),
[Quantum ESPRESSO](https://www.quantum-espresso.org/),
[CP2K](https://www.cp2k.org/),
[FHI-aims](https://fhi-aims.org/),
[CRYSTAL](https://www.crystal.unito.it/),
[SIESTA](https://siesta-project.org/),
[GPAW](https://gpaw.readthedocs.io/),
[FLEUR](https://www.flapw.de/),
[WIEN2k](https://www.tuwien.at/en/tch/tc/home-of-wien2k),
[Elk](https://elk.sourceforge.io/),
[ABINIT](https://www.abinit.org/)

**QM Programs — Multireference (3):**
[OpenMolcas](https://openmolcas.org/),
[BAGEL](https://nubakery.org/),
[Columbus](https://www.univie.ac.at/columbus/)

**Semi-empirical & Workflow Helpers (8):**
[xTB](https://github.com/grimme-lab/xtb),
[CREST](https://github.com/crest-lab/crest),
[MOPAC](https://openmopac.net/),
[Sparrow](https://scine.ethz.ch/download/sparrow),
[DFTB+](https://dftbplus.org/),
[xTB4STDA](https://github.com/grimme-lab/xtb4stda),
[sTDA](https://github.com/grimme-lab/std2),
[sTD2](https://github.com/grimme-lab/std2)

**MD Engines (5):**
[GROMACS](https://www.gromacs.org/),
[LAMMPS](https://www.lammps.org/),
[AMBER](https://ambermd.org/),
[NAMD](https://www.ks.uiuc.edu/Research/namd/),
[OpenMM](https://openmm.org/)

**AI/ML — Foundation Models (3):**
[MoLFormer](https://github.com/IBM/molformer),
[Uni-Mol](https://github.com/deepmodeling/Uni-Mol),
[ChemBERTa](https://huggingface.co/seyonec/ChemBERTa-zinc-base-v1)

**AI/ML — Generative (2):**
[REINVENT4](https://github.com/MolecularAI/REINVENT4),
[SyntheMol](https://github.com/swansonk14/SyntheMol)

**AI/ML — Conformers (2):**
[GeoMol](https://github.com/PattanaikL/GeoMol),
[torsional-diffusion](https://github.com/gcorso/torsional-diffusion)

**AI/ML — Crystal Generation (2):**
[MatterGen](https://github.com/microsoft/mattergen),
[CDVAE](https://github.com/txie-93/cdvae)

**AI/ML — Retrosynthesis (3):**
[AiZynthFinder](https://github.com/MolecularAI/aizynthfinder),
[RXNMapper](https://github.com/rxn4chemistry/rxnmapper),
[LocalRetro](https://github.com/kaist-amsg/LocalRetro)

**AI/ML — Screening / ADMET (2):**
[DeepChem](https://deepchem.io/),
[ADMETlab](https://admetlab3.scbdd.com/)

**AI/ML — Metal Complex ML (2):**
[molSimplify](https://molsimplify.mit.edu/),
[architector](https://github.com/lanl/Architector)

**Analysis / Post-Processing (15):**
[cclib](https://cclib.github.io/),
[Multiwfn](http://sobereva.com/multiwfn/),
[CENSO](https://github.com/grimme-lab/CENSO),
[ANMR](https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo_nmr.html),
[c2anmr](https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo_nmr.html),
[nmrplot](https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo_nmr.html),
[morfeus](https://github.com/digital-chemistry-laboratory/morfeus),
[nglview](https://github.com/nglviewer/nglview),
[Packmol](https://m3g.github.io/packmol/userguide.shtml),
[NBO](https://nbo6.chem.wisc.edu/),
[AIMAll](https://aim.tkgristmill.com/),
[critic2](https://aoterodelaroza.github.io/),
[Chargemol](https://sourceforge.net/projects/ddec/),
[LOBSTER](http://www.cohp.de/),
[phonopy](https://phonopy.github.io/phonopy/)

**Wrapper Libraries (5):**
[ASE](https://ase-lib.org/),
[pymatgen](https://pymatgen.org/),
[QCEngine](https://github.com/MolSSI/QCEngine),
[MDAnalysis](https://www.mdanalysis.org/),
[pymolpro](https://github.com/molpro/pymolpro)

**Visualization (6):**
[plotly](https://plotly.com/python/),
[VMD](https://www.ks.uiuc.edu/Research/vmd/),
[Avogadro](https://avogadro.cc/),
[Jmol](https://jmol.sourceforge.net/),
[ChimeraX](https://www.cgl.ucsf.edu/chimerax/),
[IQmol](https://www.iqmol.org/index.html)

**Python-Only Backends (2):**
[PySCF](https://pyscf.org/),
[PLAMS](https://www.scm.com/doc/plams/)

**Crystal Structure Prediction (1):**
[Genarris](https://github.com/Yi5817/Genarris)

</details>

---

## 🖥️ Cluster and HPC

* **Backends:** local execution, SLURM, SSH-based remote transfer
* **Automatic resource detection:** CPUs and memory on SLURM/PBS/LSF
* **Site profiles:** module names, node sizes and partitions per site, with automatic detection of the login node
* **Node-local staging:** the virtual environment and ORCA are staged onto node-local disk so Python does not start off a network home
* **Scratch directory:** `DELFIN_SCRATCH=/path/to/scratch`
* **Queue insight:** partitions, pending reasons, fairshare and account limits in the Job Status tab
* **Logging:** `delfin_run.log` per workspace, `occupier.log` per subprocess

**Cluster templates:** see `examples/` for SLURM, PBS and LSF submit scripts.

---

## 🩺 Troubleshooting

* **`CONTROL.txt` not found** — create it with `delfin --define` (or copy your own).
* **Input file not found** — run `delfin --define=your.xyz` to convert an XYZ.
* **ORCA not found** — check that `orca` is on your PATH (`which orca`), then `delfin doctor`.
* **CREST/xTB tools missing** — install them and add them to PATH, or disable the corresponding flags in `CONTROL.txt`.
* **Optional tool not detected** — check the dashboard's Settings tab, or run `delfin doctor` (and `delfin qm_check`, `delfin mlp_check`, `delfin analysis_check` for the specific families).
* **A SMILES run stops immediately** — `smiles_converter` must be set to `QUICK`, `NORMAL`, `MANTA` or `ARCHITECTOR`.

---

## 🛠️ Development

* CLI entry points are defined in `pyproject.toml`
* Development install: `pip install -e ".[agent,docs,dev]"`
* Run the fast test suite: `pytest -m "not slow"`
* Format: `black .` · Lint: `ruff check .`
* Build a wheel: `pip wheel .`

Contributions are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

---

## References

The generic references for ORCA, xTB and CREST are:

- Frank Neese. The ORCA program system. *Wiley Interdiscip. Rev. Comput. Mol. Sci.*, 2(1):73–78, 2012. doi:<https://doi.wiley.com/10.1002/wcms.81>.
- Frank Neese. Software update: the ORCA program system, version 4.0. *Wiley Interdiscip. Rev. Comput. Mol. Sci.*, 8(1):e1327, 2018. doi:<https://doi.wiley.com/10.1002/wcms.1327>.
- Frank Neese, Frank Wennmohs, Ute Becker, and Christoph Riplinger. The ORCA quantum chemistry program package. *J. Chem. Phys.*, 152(22):224108, 2020. doi:<https://aip.scitation.org/doi/10.1063/5.0004608>.
- Christoph Bannwarth, Erik Caldeweyher, Sebastian Ehlert, Andreas Hansen, Philipp Pracht, Jan Seibert, Sebastian Spicher, and Stefan Grimme. Extended tight-binding quantum chemistry methods. *WIREs Comput. Mol. Sci.*, 11:e1493, 2021. doi:<https://doi.org/10.1002/wcms.1493>. *(xTB & GFN methods)*
- Philipp Pracht, Stefan Grimme, Christoph Bannwarth, Florian Bohle, Sebastian Ehlert, Gunnar Feldmann, Jan Gorges, Max Müller, Timo Neudecker, Christoph Plett, Sebastian Spicher, Pascal Steinbach, Piotr A. Wesołowski, and Fabian Zeller. CREST — A program for the exploration of low-energy molecular chemical space. *J. Chem. Phys.*, 160:114110, 2024. doi:<https://doi.org/10.1063/5.0197592>. *(CREST)*

Please always check the output files — at the end you will find a list of relevant papers for the calculations. Kindly cite them. Please do not only cite the generic references above, but also the
[original papers](https://www.faccts.de/docs/orca/6.0/manual/contents/public.html) that report the development and ORCA implementation of the methods DELFIN has used. The publications that describe the functionality implemented in ORCA are given in the manual.

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

- Hartmann, M. (2026). *DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials* (v1.2.0). Zenodo. https://doi.org/10.5281/zenodo.21508364
- Hartmann, M. (2025). *DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials*. ChemRxiv. https://doi.org/10.26434/chemrxiv-2025-4c256

### BibTeX
```bibtex
@software{hartmann2025delfin,
  author  = {Hartmann, Maximilian},
  title   = {DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials},
  version = {v1.2.0},
  year    = {2026},
  publisher = {Zenodo},
  doi     = {10.5281/zenodo.21508364},
  url     = {https://doi.org/10.5281/zenodo.21508364}
}

@article{hartmann2025chemrxiv,
  author  = {Hartmann, Maximilian},
  title   = {DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials},
  journal = {ChemRxiv},
  year    = {2025},
  doi     = {10.26434/chemrxiv-2025-4c256},
  url     = {https://www.cambridge.org/engage/chemrxiv/article-details/68fa0e233e6156d3be78797a}
}
```

---

## License

This project is licensed under the GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later).

You should have received a copy of the GNU Lesser General Public License along with this repository in the file [LICENSE](LICENSE).
If not, see <https://www.gnu.org/licenses/>.

Non-binding citation request:
If you use this software in research, please cite the associated paper (see [CITATION.cff](./CITATION.cff)).
