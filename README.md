<h1>
  <img src="delfin/logo/DELFIN_logo.png" alt="DELFIN logo" width="60" align="center">
  <span style="color:#1976d2;">DELFIN</span>
</h1>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17208145.svg)](https://doi.org/10.5281/zenodo.17208145)
[![PyPI version](https://img.shields.io/pypi/v/delfin-complat.svg)](https://pypi.org/project/delfin-complat/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/delfin-complat.svg)](https://pypistats.org/packages/delfin-complat)

> Preprint: Hartmann, M. et al. "DELFIN: Automated DFT-based prediction of preferred spin states and corresponding redox potentials", ChemRxiv (2025). https://doi.org/10.26434/chemrxiv-2025-4c256

DELFIN is a computational chemistry platform for automated molecular workflows. It combines file-based chemistry workflows, a browser dashboard, a workflow/pipeline layer, and an optional multi-agent assistant in one codebase.

The repository has grown beyond the original redox workflow and now covers:

- automated redox and spin-state workflows (`OCCUPIER`, `classic`, `manually`)
- thermodynamics / stability-constant workflows
- excited-state dynamics, spectra, and report generation
- SMILES-to-3D conversion, GUPPY sampling, xTB / CREST / GOAT pre-optimisation
- local and SLURM execution backends
- a Jupyter / Voila dashboard with job management and setup tooling
- a step-based tool layer and YAML-defined pipelines
- optional AI-agent support inside the dashboard

## Main Entry Points

| Interface | Purpose |
| --- | --- |
| `delfin` | Main CLI for workspace-based DELFIN runs from `CONTROL.txt` |
| `delfin-voila` | Launch the dashboard as a browser app |
| `python -m delfin` | Module entry point equivalent to `delfin` |
| `delfin-step` | Run one registered step such as `xtb_opt` or `orca_sp` |
| `delfin-pipeline` | Execute a YAML-defined multi-step pipeline |
| `from delfin.dashboard import create_dashboard` | Start the dashboard from Python / Jupyter |
| `from delfin.api import run` | Programmatic execution with CLI semantics |

## Capabilities

### Core chemistry workflows

- redox potentials with automated spin-state handling
- OCCUPIER-based broken-symmetry workflows for transition-metal chemistry
- classic and manually configured redox sequences
- thermodynamics / reaction-mode stability-constant workflows
- excited-state dynamics (ESD), TD-DFT, delta-SCF, and AFP / IR / NMR reporting
- xTB-based TADF screening and hyperpolarizability workflows

### Structure generation and sampling

- SMILES to XYZ conversion for organic molecules and metal complexes
- quick conversion, normal conversion, GUPPY multi-start generation, and ARCHITECTOR-based build options
- xTB optimisation, XTB_GOAT, CREST, and xTB solvator integrations
- build-up workflow for stepwise metal-complex assembly

### Execution and orchestration

- local execution backend with configurable resource limits
- SLURM backend with submit-template support
- global resource coordination across parallel workflows
- automatic retry / recovery logic for ORCA failures
- job prioritisation for workflow bottlenecks

### Interfaces and automation

- browser dashboard for submission, monitoring, recalc, archive browsing, and settings
- agent tab for dashboard control, analysis, and coding / research workflows
- Python step registry and declarative pipeline builder
- YAML pipelines and single-step CLI execution

## Installation

PyPI package: https://pypi.org/project/delfin-complat/

### Requirements

- Python `3.10` or `3.11`
- ORCA `6.1.1` available via `PATH` or configured in the DELFIN runtime settings
- optional for some workflows: `xtb`, `crest`, `xtb4stda`, `stda`, `std2`
- optional for some analysis/report workflows: `censo`, `anmr`, `c2anmr`, `nmrplot`
- optional for dashboard usage: JupyterLab / Notebook or Voila

### Standard install

```bash
pip install delfin-complat
```

### Editable install from source

```bash
git clone https://github.com/ComPlat/DELFIN.git
cd DELFIN
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Optional extras

Source installs expose additional extras, for example:

```bash
pip install -e ".[agent,analysis,mlp]"
```

DELFIN does not bundle licensed or external QM binaries. Those tools must be installed separately and made visible to the runtime.

### QM tool bundle helpers

DELFIN ships a managed `qm_tools` helper layer. After installing the Python package, you can stage and validate the shared tool environment with:

```bash
source delfin/qm_tools/env.sh
USE_SYSTEM_TOOLS=1 bash delfin/qm_tools/install_qm_tools.sh
bash delfin/qm_tools/check_qm_tools.sh
```

The dashboard Settings tab exposes the same runtime setup flows graphically.

## Quick Start

Create a workspace with `CONTROL.txt` and `input.txt`:

```bash
mkdir my_project
cd my_project
delfin --define
```

If you already have an XYZ file, DELFIN can build the matching input template:

```bash
delfin --define=structure.xyz
```

Then edit `CONTROL.txt` and set at minimum:

- `charge`
- redox steps such as `oxidation_steps` / `reduction_steps`
- workflow method (`classic`, `manually`, or `OCCUPIER`)
- level-of-theory and resource settings as needed

Run DELFIN inside the workspace:

```bash
delfin
```

Equivalent module entry point:

```bash
python -m delfin
```

You can also point DELFIN at another project directory:

```bash
delfin /path/to/project
```

Typical outputs include:

- `DELFIN.txt`
- `DELFIN_Data.json`
- step-specific `.inp`, `.out`, `.xyz`, and report files
- workflow directories such as `initial_OCCUPIER`, `ox_step_1_OCCUPIER`, or `thermodynamics`

## Common CLI Tasks

```bash
delfin --help
delfin --recalc
delfin --imag
delfin --report text
delfin --report docx
delfin --json
delfin cleanup --help
delfin qm_check
delfin-step --list
delfin-pipeline workflow.yaml --cores 8
```

## Dashboard

Start the dashboard in Jupyter:

```python
from delfin.dashboard import create_dashboard

ctx = create_dashboard(backend="auto")
```

`backend="auto"` uses SLURM when `sbatch` is available and otherwise falls back to local execution.

Start the Voila app:

```bash
delfin-voila
delfin-voila --port 9000
delfin-voila --ip 0.0.0.0
```

Main dashboard areas:

- Submit Job
- Recalc
- Job Status
- ORCA Builder
- TURBOMOLE Builder (backend-dependent)
- Calculations Browser
- Archive
- Agent
- Settings

## Agent System

The dashboard includes an optional agent tab for:

- dashboard control through natural language
- result inspection and workflow troubleshooting
- literature / method research
- code-editing and implementation workflows

Supported provider families in the current codebase include Claude, OpenAI / Codex, and KIT Toolbox. DELFIN treats these like external runtime dependencies: the user installs the CLI or provides the API key, and DELFIN detects what is available.

## Documentation

- [User Manual](docs/USER_MANUAL.md)
- [Settings and Setup](docs/SETTINGS_AND_SETUP.md)
- [Methodology](docs/methodology.md)
- [Retry Logic](docs/RETRY_LOGIC.md)
- [Job Prioritization](docs/job_prioritization.md)
- [Example submission scripts](examples/example_Job_Submission_Scripts/README.md)

## Repository Layout

| Path | Purpose |
| --- | --- |
| `delfin/` | main Python package |
| `delfin/dashboard/` | dashboard UI and backends |
| `delfin/workflows/` | workflow orchestration and scheduling |
| `delfin/tools/` | step registry, adapters, declarative pipelines |
| `delfin/agent/` | agent orchestration, prompts, session state |
| `docs/` | user and technical documentation |
| `examples/` | sample projects and cluster submission scripts |
| `tests/` | automated test suite |

## Citation

If DELFIN is useful in academic work, please cite the ChemRxiv preprint above and the repository / release metadata in [CITATION.cff](CITATION.cff).

## License

LGPL-3.0-or-later. See [LICENSE](LICENSE).
