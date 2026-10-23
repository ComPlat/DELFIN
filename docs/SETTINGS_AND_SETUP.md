# Settings, Runtime, and Setup

This document explains the DELFIN Settings tab, the runtime model, and the setup helpers for local use and SLURM clusters.

## Overview

DELFIN stores user-specific settings in:

`~/.delfin_settings.json`

This file is outside the git repository and is not overwritten by `git pull`.

The Settings tab is designed around a simple split:

- `Workspace Paths`: where DELFIN reads and writes calculations and archive data.
- `Transfer Target`: SSH target for remote archive / transfer features.
- `Runtime / Execution`: how DELFIN resolves ORCA, qm tools, local resources, and SLURM / cluster behavior.

## Workspace Paths

### Calculations

The root directory used for active DELFIN jobs and project folders.

Default:

- `~/calc`

### Archive

The root directory used for archived projects.

Default:

- `~/archive`

If these fields are left empty, DELFIN uses the defaults. If they are set, the values are stored in `~/.delfin_settings.json`.

## Transfer Target

This section configures the SSH target used by Remote Archive and transfer-related actions.

Fields:

- `Host`
- `User`
- `Port`
- `Remote Path`

Notes:

- Passwords are not stored in DELFIN.
- SSH keys or `ssh-agent` should be used.
- The SSH helper text in the Settings tab is only guidance; it does not configure SSH automatically.

## Runtime / Execution

This section controls how DELFIN finds external tools and how it behaves in local or SLURM environments.

### Backend

Available values:

- `Auto`
- `Local`
- `SLURM`

Behavior:

- `Auto` chooses `SLURM` if `sbatch` is available.
- Otherwise DELFIN falls back to `Local`.

### ORCA path

Global ORCA override used by DELFIN if set.

If empty, DELFIN prefers:

- `DELFIN_ORCA_BASE`
- `ORCA_BINARY`
- `ORCA_PATH`
- `PATH`
- DELFIN auto-detection

### Detected ORCA

This dropdown is filled by `Scan ORCA`.

It searches:

- current DELFIN runtime candidates
- `PATH`
- `~/software`
- `~/apps`
- `~/local`
- `/opt`

Selecting an entry fills the global ORCA path field automatically.

### qm_tools root

Optional override for the DELFIN qm tool bundle.

Typical user-managed location:

- `~/.delfin/qm_tools`

This is where DELFIN-managed `xtb`, `crest`, `dftb+`, `stda`, and `xtb4stda` assets are expected when you use the setup buttons.

### Local ORCA

Optional ORCA override used only for the local backend.

### Local max cores / Local max RAM

Local resource limits used by the local backend queue.

These are now system-dependent by default and can be refreshed with `Detect local resources`.

### SLURM ORCA

Optional ORCA override used only for the SLURM backend.

### Submit templates

Optional override for the SLURM submit template directory.

If empty, DELFIN uses:

- the configured repo templates when available
- otherwise the packaged default submit templates

### Site profile

Optional label for the current SLURM site profile, for example:

- `bwunicluster3`

### Runtime validation

This section shows the effective runtime resolution and checks what DELFIN currently sees, for example:

- backend
- ORCA
- local runner or `sbatch`
- `xtb`
- `crest`
- `std2`
- `stda`
- `xtb4stda`
- `dftb+`
- `qm_tools_root`

`Validate Setup` updates this check view without saving settings.

## Buttons

### Save Settings

Writes the current Settings tab values to `~/.delfin_settings.json`.

This is a writing action.

### Reload

Reloads the current settings from disk into the UI.

This is a read-only action.

### Validate Setup

Runs the runtime diagnostics shown in the `Runtime validation` box.

This is a read-only action.

### Scan ORCA

Searches for ORCA installations in standard locations and updates the `Detected ORCA` dropdown.

This is a read-only action.

### Detect local resources

Detects CPU count and total RAM from the current system and updates the local resource fields in the UI.

This does not save automatically. Use `Save Settings` if you want to persist the detected values.

### Prepare qm_tools

Copies the bundled DELFIN qm tool bundle into the user area, usually:

- `~/.delfin/qm_tools`

This is a writing action.

It does not try to download or update everything. It only stages the packaged bundle.

### Install qm_tools

Runs the DELFIN qm tool installer against the user qm tools directory.

This is a writing action.

Use this to make a local DELFIN installation more complete without touching the Python package directory itself.

### Update qm_tools

Refreshes DELFIN-managed qm tool assets and updates DELFIN-managed local environments where applicable.

Important:

- tools already provided by `PATH` or an activated `.venv` are still accepted
- DELFIN does not aggressively rewrite arbitrary external environments

This is a writing action.

### Setup cluster

Prepares an existing DELFIN installation for the SLURM cluster it runs on.

It configures DELFIN-side behavior, for example:

- runtime backend `slurm`
- the site profile, when DELFIN knows the site (`bwunicluster3` on bwUniCluster 3.0; none elsewhere)
- submit template wiring
- `calc` / `archive`
- `~/.delfin_env.sh`, including the OpenMPI that ORCA uses
- shell sourcing

There is no venv tarball to prepare any more: a SLURM job packs one from the venv as it is when the job starts, keyed by the venv's contents, and caches it in `~/.cache/delfin/venv` (`DELFIN_VENV_CACHE_DIR`). A changed venv gets a new tar; an unchanged one is reused.

This is a writing action.

It does **not** perform the full system installation.

### Partitions a job may start in

`runtime.slurm.partitions` in `~/.delfin_settings.json` is a comma-separated list of partitions a CPU job may start in, for example `"cpu,cpu_il"`. SLURM starts the job in whichever can run it first. Before submitting, DELFIN asks SLURM with `sbatch --test-only` which of them can hold the request and lists only those. Empty uses the site profile's list (`cpu,cpu_il` on bwUniCluster 3.0), or the submit template's own partition elsewhere. `DELFIN_SLURM_PARTITIONS` sets the same from the environment. `runtime.slurm.gpu_partitions` (or `DELFIN_SLURM_GPU_PARTITIONS`) does the same for jobs that ask for a GPU. Empty uses the site profile's list, or every GPU partition that `scontrol show partition` reports, apart from `dev_*` test partitions. Each is checked with `sbatch --test-only --gres=gpu:1` against the job's time, cores and memory; when none fits, the job runs on CPUs.

### Verify install

Performs a read-only readiness check.

It checks items such as:

- installer availability
- repo / `.venv`
- OpenMPI
- ORCA / `orca_plot`
- `calc` / `archive`
- submit templates
- `sbatch`

This does **not** modify files or repair anything.

### Full install

Runs `delfin/installers/install_delfin.sh --profile core`, the same installer as `install.sh`.

This is the heavy setup action.

It can handle things such as:

- OpenMPI build/install
- ORCA reuse if already installed
- ORCA extraction from a matching tarball if present
- repo / `.venv` setup
- `pip install -e`
- runtime defaults
- `~/.delfin_env.sh`

This is a writing action and may take a while.

Important:

- ORCA is still treated as external software
- DELFIN does not fetch ORCA from the internet
- it only reuses an existing ORCA installation or a user-provided ORCA tarball

### Update all tools / Repair all tools

**Update all tools** fetches every installed tool again: QM programs, analysis, ML and AI packages, Ketcher. **Repair all tools** checks each installed tool and puts right what does not work -- a link into a deleted environment repointed, a program that does not start or a package that does not import reinstalled -- and proves each repair afterwards. A tool that is not installed is left alone.

Both go through `delfin.installer`, like `install.sh --update` and `install.sh --repair`. Updating DELFIN itself (git pull and pip) is `install.sh --update`.

## Local Usage

Typical local workflow:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
delfin-voila
```

Then in Settings:

1. `Scan ORCA`
2. choose ORCA if needed
3. `Install qm_tools`
4. `Detect local resources`
5. `Save Settings`

This usually makes a local DELFIN installation fully usable except for ORCA itself, which must already exist.

## PyPI Usage

Typical PyPI workflow:

```bash
pip install delfin-complat
delfin-voila
```

The same Settings workflow applies:

1. resolve ORCA
2. install or update `qm_tools`
3. validate

The packaged installer and submit templates are included so DELFIN can still perform setup actions without requiring a git checkout.

## Cluster Usage

From a shell on the login node, `install.sh` does everything (see the README). In the dashboard there are two intended modes:

### Existing installation already present

Use:

- `Setup cluster`

This is the lighter DELFIN-side preparation path.

### Full cluster-oriented installation

Use:

- `Full install`

This is the heavier path that runs the installer itself.

Recommended prerequisites:

- run it on the cluster's login node
- provide ORCA either as an existing installation or as its tarball in `$HOME` or `~/software`

### Safe check before changing anything

Use:

- `Verify install`

This is the recommended first step when you want to understand what is present or missing without changing the system.

## Compatibility Notes

The current design is intended to remain compatible with:

- editable installs: `pip install -e .`
- PyPI / wheel installs
- existing repo-based layouts
- existing `software/` layouts
- existing local setups
- existing bwUniCluster setups

The main goal is not to remove the old workflows, but to make runtime resolution and setup behavior explicit and accessible from the UI.
