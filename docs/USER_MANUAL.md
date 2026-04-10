# DELFIN User Manual

This manual is a practical guide for running DELFIN as a user. It focuses on the current user-facing workflows in the repository: workspace-based CLI runs, dashboard usage, runtime setup, and the newer step / pipeline interfaces.

For implementation details and scientific background, see:

- [README](../README.md)
- [Settings and Setup](SETTINGS_AND_SETUP.md)
- [Methodology](methodology.md)
- [Retry Logic](RETRY_LOGIC.md)
- [Job Prioritization](job_prioritization.md)

## 1. What DELFIN Expects

DELFIN is organised around a workspace directory. A typical workspace contains:

- `CONTROL.txt`: workflow and runtime configuration
- `input.txt`: geometry or SMILES input in DELFIN's expected format
- generated `.inp`, `.out`, `.xyz`, and report files

Depending on the workflow, DELFIN may also create subdirectories such as:

- `initial_OCCUPIER`
- `ox_step_1_OCCUPIER`
- `red_step_1_OCCUPIER`
- `thermodynamics`
- `ESD`

The main mental model is simple:

1. prepare a workspace
2. define the chemistry and resources in `CONTROL.txt`
3. run DELFIN
4. inspect results, rerun missing pieces, or generate reports

## 2. Installation and Prerequisites

### Python package

Standard install:

```bash
pip install delfin-complat
```

Editable source install:

```bash
git clone https://github.com/ComPlat/DELFIN.git
cd DELFIN
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Required external software

For the main DFT workflows you need:

- Python `3.10` or `3.11`
- ORCA `6.1.1`

Optional tools become relevant depending on the workflow:

- `xtb`, `crest` for pre-optimisation and conformer search
- `xtb4stda`, `stda`, `std2` for xTB excited-state / screening workflows
- `censo`, `anmr`, `c2anmr`, `nmrplot` for ensemble NMR workflows
- Jupyter or Voila for the dashboard

### Runtime helpers

DELFIN includes a managed `qm_tools` helper bundle:

```bash
source delfin/qm_tools/env.sh
USE_SYSTEM_TOOLS=1 bash delfin/qm_tools/install_qm_tools.sh
bash delfin/qm_tools/check_qm_tools.sh
```

The same setup actions are available in the dashboard Settings tab.

## 3. Creating Your First Workspace

### Empty template

```bash
mkdir my_project
cd my_project
delfin --define
```

This creates:

- `CONTROL.txt`
- `input.txt`

### Start from an XYZ file

```bash
delfin --define=structure.xyz
```

If the input file ends in `.xyz`, DELFIN converts it into `input.txt` by dropping the first two XYZ header lines and updates `CONTROL.txt` accordingly.

### Start in a different directory

```bash
delfin /path/to/project --define
```

## 4. Input Files

### `input.txt`

`input.txt` can contain either:

- atomic coordinates in DELFIN's text format
- a one-line SMILES string

If you use a full XYZ file, DELFIN can convert it during `--define` or when a `CONTROL.txt` file points to an `.xyz` geometry.

### `CONTROL.txt`

`CONTROL.txt` is the main user interface for DELFIN workflows. The default template is extensive because it covers multiple modules. In practice, only a smaller subset is needed for most jobs.

The most important sections are:

- input and system identity
- redox workflow selection
- level of theory
- solvation
- resource settings
- optional modules such as thermodynamics, ESD, or xTB screening

## 5. Minimum `CONTROL.txt` Settings

For a normal redox run, review at least these keys:

```ini
input_file=input.txt
charge=0
calc_initial=yes
oxidation_steps=1
reduction_steps=1
method=OCCUPIER
functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
solvent=acetonitrile
PAL=16
maxcore=4000
parallel_workflows=yes
```

Important notes:

- `method` selects the workflow family: `OCCUPIER`, `classic`, or `manually`
- `oxidation_steps` and `reduction_steps` define how many redox steps DELFIN should run
- `PAL` and `maxcore` control ORCA resource usage
- DELFIN's scheduler can split resources across parallel branches when `parallel_workflows=yes`

## 6. Main Workflow Modes

### OCCUPIER

Use `OCCUPIER` for automated spin-state exploration and broken-symmetry handling, especially for transition-metal systems and cases where electronic-state assignment is part of the problem.

Typical strengths:

- spin-state screening
- broken-symmetry workflows
- automated oxidation / reduction branch management
- basis-set special handling for metals and the first coordination sphere

### classic

Use `classic` when you want straightforward DFT calculations for each state without the full OCCUPIER machinery.

Typical strengths:

- conventional redox workflows
- simpler systems
- easier interpretation of per-state jobs

### manually

Use `manually` when you want to control multiplicities and broken-symmetry choices explicitly.

Typical strengths:

- expert-driven setups
- curated spin-state assignments
- reproducing a known workflow exactly

## 7. Running DELFIN

### Standard run

```bash
delfin
```

### Run a different workspace

```bash
delfin /path/to/project
```

### Use a specific `CONTROL.txt`

```bash
delfin --control /path/to/CONTROL.txt
```

### Keep intermediates

```bash
delfin --no-cleanup
```

This is useful when you want to inspect scratch-stage inputs and intermediate files after the run.

## 8. Common CLI Operations

### Recalculate only missing or incomplete jobs

```bash
delfin --recalc
```

This keeps existing results and reruns only jobs whose outputs are missing or clearly incomplete.

### Run only the IMAG cleanup workflow

```bash
delfin --imag
```

Use this when you already have converged jobs and want DELFIN to eliminate imaginary frequencies and regenerate the summary.

### Rebuild reports from existing outputs

```bash
delfin --report text
delfin --report docx
delfin --json
delfin --afp
```

### Cleanup or purge

```bash
delfin cleanup --help
delfin --cleanup
delfin --purge
```

`--purge` is intentionally stricter and removes recognised DELFIN-generated artefacts while preserving `CONTROL.txt` and the referenced input file.

### Inspect tool resolution

```bash
delfin qm_check
delfin qm_check xtb crest stda
delfin qm_run xtb -- --version
```

## 9. Understanding Outputs

Common top-level outputs include:

- `DELFIN.txt`: text summary
- `DELFIN.docx`: Word report when requested
- `DELFIN_Data.json`: structured data export
- `delfin_run.log`: run-level log file
- ORCA input and output files for individual jobs

For OCCUPIER-based workflows, you will often also see:

- `OCCUPIER.txt`
- `occupier.log`
- OCCUPIER stage directories

For thermodynamics:

- `thermodynamics/`

For excited-state work:

- `ESD/`

## 10. Optional Workflow Modules

### Thermodynamics / stability constants

The control template still exposes the user-facing `thermodynamics` keys:

```ini
thermodynamics=yes
thermodynamics_mode=auto
```

Reaction mode is also available:

```ini
thermodynamics=yes
thermodynamics_mode=reaction
thermodynamics_reaction=a*{SMILES}+b*{SMILES}...>>>c*{SMILES}+d*{SMILES}...
thdy_smiles_converter=ARCHITECTOR
thdy_preopt=xtb
```

Internally, the code normalises these settings onto the `stability_constant` workflow family. From a user perspective, you can continue using the `thermodynamics*` keys shown in the default template.

### ESD module

The excited-state dynamics module is controlled through keys such as:

```ini
ESD_modul=yes
ESD_modus=TDDFT
states=S1,T1,S2,T2
ISCs=S1>T1,T1>S1
ICs=S2>S1
```

This creates a dedicated ESD workflow on top of the ground-state DELFIN run.

### Hyperpolarizability and TADF xTB workflows

These modules are enabled through:

```ini
hyperpol_xTB=yes
tadf_xTB=yes
```

They are intended for specialised xTB-based property and screening workflows and rely on additional external binaries.

## 11. Structure Generation and Sampling

DELFIN can start from SMILES and build candidate geometries through several routes:

- quick conversion
- normal conversion
- GUPPY multi-start sampling
- ARCHITECTOR-based complex generation

Relevant user-facing tools:

```bash
delfin-guppy --help
delfin-build --help
```

`delfin-guppy` is useful when one SMILES string may correspond to several chemically relevant starting geometries and you want DELFIN to rank them after xTB optimisation.

## 12. Dashboard Usage

### Start from Python / Jupyter

```python
from delfin.dashboard import create_dashboard

ctx = create_dashboard(backend="auto")
```

### Start with Voila

```bash
delfin-voila
delfin-voila --port 9000
delfin-voila --ip 0.0.0.0
```

### Main dashboard tabs

- `Submit Job`: create or edit inputs, convert SMILES, preview structures, trigger submission
- `Recalc`: modify and resubmit an existing project
- `Job Status`: inspect local or SLURM queue state
- `ORCA Builder`: prepare ORCA inputs interactively
- `TURBOMOLE Builder`: available on supported backends
- `Calculations Browser`: browse files, inspect outputs, trigger browser-side workflows
- `Archive`: inspect archived jobs and statistics
- `Agent`: use the built-in assistant
- `Settings`: configure paths, runtime, and tool setup

### Backends

The dashboard supports:

- `auto`
- `local`
- `slurm`

`auto` chooses SLURM when `sbatch` is available and otherwise falls back to the local backend.

## 13. Settings and Runtime Configuration

User-specific runtime settings are stored in:

```text
~/.delfin_settings.json
```

The Settings tab controls:

- calculations and archive root directories
- transfer target settings
- backend choice
- ORCA path overrides
- `qm_tools` root
- local CPU and RAM limits
- SLURM profile and submit-template locations

The most important setup actions are:

- `Validate Setup`
- `Scan ORCA`
- `Prepare qm_tools`
- `Install qm_tools`
- `Update qm_tools`
- `Setup bwUniCluster`
- `Verify bwUniCluster`

For details, see [SETTINGS_AND_SETUP.md](SETTINGS_AND_SETUP.md).

## 14. Single Steps and YAML Pipelines

The repository now includes a newer step / pipeline layer in addition to the classic workspace-based CLI.

### Run one step

```bash
delfin-step --list
delfin-step xtb_opt --geometry input.xyz --charge 0 --cores 8
delfin-step orca_sp --geometry input.xyz --charge 0 --method B3LYP --basis def2-SVP
```

### Run a YAML pipeline

```bash
delfin-pipeline workflow.yaml --cores 8
delfin-pipeline workflow.yaml --cores auto --scheduled
```

This interface is useful for reusable programmatic workflows, screening jobs, and experiments that are easier to express as a step graph than as one large `CONTROL.txt`.

## 15. Python API

If you want CLI-like execution from Python:

```python
from delfin.api import run

run(control_file="CONTROL.txt", recalc=False, cleanup=True)
```

For dashboard use:

```python
from delfin.dashboard import create_dashboard

ctx = create_dashboard(backend="auto")
```

For custom pipelines, use the `delfin.tools` and `delfin.workflows` layers.

## 16. Troubleshooting

### DELFIN cannot find ORCA

Check:

- `which orca`
- dashboard `Validate Setup`
- dashboard `Scan ORCA`
- runtime overrides in `~/.delfin_settings.json`

### xTB / CREST / STDA tools are missing

Check:

```bash
delfin qm_check
delfin qm_check xtb crest xtb4stda stda std2
```

### A job failed but should be recoverable

DELFIN includes retry logic controlled by:

```ini
enable_auto_recovery=yes
max_recovery_attempts=3
```

State is tracked in:

```text
.delfin_recovery_state.json
```

See [RETRY_LOGIC.md](RETRY_LOGIC.md).

### Parallel jobs are not behaving as expected

Review:

- `PAL`
- `pal_jobs`
- `parallel_workflows`
- `orca_parallel_strategy`

See [job_prioritization.md](job_prioritization.md) for scheduler-specific behaviour.

### The dashboard behaves differently on local and cluster systems

That is expected. DELFIN resolves the backend dynamically and supports both local and SLURM-specific runtime configuration. The Settings tab is the canonical place to inspect what backend and paths are currently active.

## 17. Further Reading

- [README](../README.md)
- [Settings and Setup](SETTINGS_AND_SETUP.md)
- [Methodology](methodology.md)
- [Retry Logic](RETRY_LOGIC.md)
- [Job Prioritization](job_prioritization.md)
- [Example job submission scripts](../examples/example_Job_Submission_Scripts/README.md)
