# Changelog

All notable changes to DELFIN will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.4] - 2025-02-XX

### Added
- Automatic per-run logging: the CLI attaches a global `delfin_run.log` and OCCUPIER subprocesses emit an `occupier.log` alongside ORCA outputs.
- CLI help now documents QM/XTB splitting via `$` markers and clarifies how `parallel_workflows` toggles between parallel and sequential scheduling modes.

### Changed
- Persist QMMM split detection (`$` separator) via shared cache so OCCUPIER/Classic/Manually runs keep the `QM/XTB` flag even if later geometries omit the marker.
- All geometry writers now pass the geometry path to the splitter so cache lookups work across workflow steps.
- OCCUPIER FoB scheduling always uses the global job manager; sequential runs simply cap `max_jobs` to one, ensuring consistent PAL enforcement.
- Global job manager initialization is now idempotent and refuses to tear down an active pool when jobs are still running; configuration differences are detected via resource signatures.
- Improved banner and logging output for the global scheduler, including explicit parallel-mode reporting.

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
- Integration with ORCA 6.1.0, xTB, and CREST
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
