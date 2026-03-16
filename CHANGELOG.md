# Changelog

All notable changes to DELFIN will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- _No unreleased entries yet._

## [1.1.9] - 2026-03-16

### Added
- New runtime and setup controls in the Settings tab, including configurable `Calculations`/`Archive` paths, backend selection, ORCA path overrides, ORCA scan/selection, and system-aware local resource detection.
- New packaged `qm_tools` workflow controls in the Settings tab: `Prepare qm_tools`, `Install qm_tools`, and `Update qm_tools`.
- New bwUniCluster actions in the Settings tab: `Setup bwUniCluster`, `Verify bwUniCluster`, and `Full bwUni install`.
- New packaged runtime resources for wheel/PyPI installs, including generic submit templates and the packaged bwUniCluster installer.
- New English setup documentation in `docs/SETTINGS_AND_SETUP.md`.

### Changed
- DELFIN runtime configuration is now stored in `~/.delfin_settings.json`, keeping user-specific paths and runtime choices outside the git repository.
- Local execution now uses a bundled Python local runner fallback when shell-based local runner scripts are not present.
- `delfin-voila` now stages packaged notebooks safely, avoids fragile server-side browser assumptions, and behaves better in remote/server sessions.
- Local and cluster setup flows are now exposed directly in the GUI while preserving compatibility with editable installs, existing repo-based installations, and established `software/` layouts.
- Packaged installs now include the resources needed for runtime setup instead of assuming a neighboring source checkout.

### Fixed
- Fixed packaged `delfin-voila` launch failures caused by notebook root-directory issues, hidden staging paths, and server/browser launch behavior.
- Fixed missing packaged runtime resources for PyPI/wheel installs, including local runner and setup assets that were previously only available from a source checkout.
- Fixed Settings-tab issues around saving configurable workspace/runtime paths and improved masking/handling of local sensitive transfer settings.

### Documentation
- Added detailed documentation for Settings, runtime controls, setup buttons, local installs, PyPI installs, and bwUniCluster workflows.

## [1.1.0] - 2026-03-03

### Added
- Dashboard/GUI expansion across Submit/Recalc/ORCA Builder/TURBOMOLE Builder/Calculations/Archive tabs, including improved archive-to-calculations workflows and MO plotting controls.
- Explicit SMILES conversion modes in the Submit tab (`CONVERT SMILES`, `QUICK CONVERT SMILES`, `CONVERT SMILES + UFF`) with improved handling of metal complexes.
- `delfin-guppy` workflow improvements: quick conversion as additional start geometry and post-XTB topology validation.
- Extended coordination-chemistry support in SMILES conversion (including additional coordination numbers and improved topology/isomer enumeration robustness).
- `ESD_T1_opt` toggle support for controlling UKS vs TDDFT T1 optimization behavior.
- Smart recalc controls and fingerprint-based skip logic for avoiding unnecessary ORCA reruns.

### Changed
- Isomer deduplication logic was tightened and then made workflow-aware: dashboard and GUPPY now preserve broader labeled variant diversity for SMILES isomer sets.
- Dashboard molecule viewer and trajectory UX refined (layout stability, playback reliability, print mode controls, cube/MO interaction refinements).
- Cleanup/copy-back behavior around `.orca_iso*` and GOAT/XTB side products hardened for better restart safety.
- HPC runtime I/O overhead reduced via caching and scratch/runtime optimizations.

### Fixed
- Multiple SMILES conversion regressions affecting metal complexes (fragment topology, ring-count checks, hydrogen handling, and fallback strategy behavior).
- Archive/statistics browser issues (path extraction, folder navigation toggles, and clipboard/export actions).
- Several ORCA input/output propagation issues in submit/recalc flows (including missing copied `.inp` and timeout/abort recovery paths).
- Additional CO2 coordinator input-generation edge cases (basis assignment and convergence keyword handling).

## [1.0.4] - 2025-02-XX

### IMAG improvements
- CLI: add `delfin --imag` entry point to rerun the IMAG workflow on existing results and regenerate DELFIN.txt.
- IMAG inputs now reuse the original ORCA template (incl. `%QMMM`, per-atom `NewGTO`) while stripping `MORead` / `%moinp` hints to force fresh guesses.
- `%pal` and `%maxcore` are inherited from the source frequency job so IMAG iterations use the same resource allocation.
- Classic, manually, and OCCUPIER pipelines call the refactored IMAG routine for initial and redox steps, copying original inputs/outputs for traceability.
- README updated with `--imag` documentation and CONTROL option `IMAG_scope` (initial/all).

### Added
- Optional excited-state dynamics (ESD) module with standalone or post-redox execution, including S0/S1/T1/T2 optimisations plus ISC/IC scheduling in a dedicated `ESD/` directory.
- ORCA input builders for ESD workflows (`esd_module.py`, `esd_input_generator.py`) and CLI/pipeline switches (`ESD_modul`, `states`, `ISCs`, `ICs`) with documentation.
- Automatic per-run logging: the CLI attaches a global `delfin_run.log` and OCCUPIER subprocesses emit an `occupier.log` alongside ORCA outputs.
- CLI help now documents QM/XTB splitting via `$` markers and clarifies how `parallel_workflows` toggles between parallel and sequential scheduling modes.
- `delfin --purge` command to wipe all intermediates (keeps CONTROL.txt + main input) with interactive confirmation.

### Changed
- Persist QMMM split detection (`$` separator) via shared cache so OCCUPIER/Classic/Manually runs keep the `QM/XTB` flag even if later geometries omit the marker.
- All geometry writers now pass the geometry path to the splitter so cache lookups work across workflow steps.
- OCCUPIER FoB scheduling always uses the global job manager; sequential runs simply cap `max_jobs` to one, ensuring consistent PAL enforcement.
- `delfin --purge` now deletes only recognized DELFIN artifacts (OCCUPIER folders, ORCA inputs/outputs, logs) and refuses to touch unknown files.
- Global job manager initialization is now idempotent and refuses to tear down an active pool when jobs are still running; configuration differences are detected via resource signatures.
- Improved banner and logging output for the global scheduler, including explicit parallel-mode reporting.
- Added `orca_parallel_strategy` CONTROL option allowing OCCUPIER runs to force serial ORCA execution (`threads`/`serial`) when MPI backends are unstable.

### Fixed
- Ensured newly spawned OCCUPIER steps reuse the cached QM range, restoring the `%QMMM` block for oxidation/reduction inputs.
- Sanitized PAL/maxcore/pal_jobs parsing so subprocesses inherit consistent limits from CONTROL.txt and `parallel_workflows`.

## [1.0.3] - 2025-01-XX

### Added
- **Global job manager singleton** for centralized resource coordination across all workflows
- **Automatic PAL splitting** for parallel oxidation/reduction workflows (e.g., PAL=12 → 6 cores per workflow)
- **Thread-safe workflow execution** with proper locking and coordination mechanisms
- **Subprocess PAL coordination** via `DELFIN_CHILD_GLOBAL_MANAGER` environment variable
- **Bootstrap mechanism** for OCCUPIER subprocesses to respect global resource limits
- New module: `global_manager.py` - Singleton GlobalJobManager class
- New module: `thread_safe_helpers.py` - Thread-safe workflow preparation functions
- New module: `cluster_utils.py` - Automatic cluster resource detection (SLURM/PBS/LSF)
- **Verification script** `verify_global_manager.py` for testing singleton behavior

### Changed
- **Refactored parallel execution** to use shared `DynamicCorePool` across all workflows
- **OCCUPIER subprocesses** now receive allocated PAL via environment variables instead of reading original CONTROL.txt
- **PAL management** centralized - read once at startup, propagated consistently to all jobs
- **CONTROL.txt updates** in subprocess folders now include reduced PAL values
- Improved logging with visual banners for parallel workflow execution
- Updated `parallel_classic_manually.py` to use global pool with intelligent fallback
- Updated `parallel_occupier.py` for global pool integration

### Fixed
- **Fixed critical issue**: Double core allocation when ox/red workflows run in parallel
- **Resolved race conditions** in parallel workflow execution through global coordination
- **Prevented PAL over-allocation** - total CPU usage never exceeds configured PAL
- Improved thread safety in file operations (CONTROL.txt, geometry files)
- Fixed subprocess resource management to respect global limits

### Documentation
- Added "Global Resource Management" section to README.md
- Expanded "Parallel Processing" section in docs/methodology.md with implementation details
- Updated project layout in README.md with new modules
- Added architectural details for global job manager pattern

## [1.0.2] - 2025-01-26

### Added
- Initial public release
- OCCUPIER, classic, and manually calculation modes
- Automated spin-state prediction and redox potential calculations
- Integration with ORCA 6.1.1, xTB, and CREST
- Support for up to 3 sequential oxidation/reduction steps
- Parallel workflow execution for oxidation/reduction steps
- Cluster environment detection (SLURM/PBS/LSF)
- CPCM/SMD implicit solvation models
- TD-DFT excited state calculations
- Broken-symmetry DFT for transition metal complexes

### Documentation
- Comprehensive README.md with installation and usage instructions
- Detailed methodology.md documentation
- Example job submission scripts for SLURM/PBS/LSF
- CITATION.cff for academic citation

---

## Release Notes

### Version 1.0.3 - Global Resource Management

This release introduces a major architectural improvement: **global job management**. Previously, when oxidation and reduction workflows ran in parallel, each could potentially allocate the full PAL (number of CPU cores), leading to double allocation (e.g., 2×12=24 cores when only 12 were available).

**Key Improvements:**
- Single source of truth for PAL allocation
- Automatic core splitting when workflows run in parallel
- Guaranteed compliance with cluster resource limits
- Thread-safe coordination between all workflows

**Impact:**
- More efficient cluster resource utilization
- No more job failures due to over-subscription
- Predictable performance in parallel execution mode

**Migration:**
No changes required to existing CONTROL.txt files or workflows. The global job manager is automatically initialized and works transparently.
