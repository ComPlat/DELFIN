# DELFIN Context

Use this context in every DELFIN agent cycle.

## Product direction

DELFIN is moving toward becoming the best platform for automation,
orchestration, and interconnection of quantum chemistry workflows across local
systems and HPC clusters.

The platform goal is more important than feature count.

## What DELFIN already is

- Python package with CLI entry points
- workflow engine with job dependency scheduling
- shared global resource coordination
- runtime/tool discovery and installation helpers
- local and SLURM dashboard backends
- CONTROL-based configuration and validation
- ORCA-centric recovery and retry logic
- reporting and JSON collection
- existing test suite and changelog

## High-value DELFIN modules

- `delfin/cli.py`: CLI orchestration and user-facing execution path
- `delfin/pipeline.py`: high-level workflow orchestration
- `delfin/parallel_classic_manually.py`: workflow scheduler and job execution
- `delfin/global_manager.py`: shared resource coordination
- `delfin/global_scheduler.py`: aggregated ORCA scheduling
- `delfin/dynamic_pool.py`: core allocation and pool execution
- `delfin/config.py`: CONTROL parsing and validation entry path
- `delfin/common/control_validator.py`: config schema and normalization
- `delfin/qm_runtime.py`: tool resolution contracts
- `delfin/runtime_setup.py`: runtime setup and diagnostics
- `delfin/orca_recovery.py`: retry and recovery logic
- `delfin/dashboard/backend_local.py`: local execution backend
- `delfin/dashboard/backend_slurm.py`: SLURM execution backend
- `delfin/api.py`: public API surface, currently thin

## Current architectural pressure points

- Monolithic orchestration files
- CLI semantics leaking into the API layer
- mutable config dictionaries carrying runtime state
- runtime and cluster behavior that must remain diagnosable
- scheduler correctness and PAL/resource safety
- recovery logic complexity

## Preferred direction

- improve architecture while solving the task
- favor explicit state and contracts over hidden config mutation
- keep public surfaces stable or make stability boundaries clearer
- protect local and HPC behavior equally
- make failures more diagnosable, not more silent
- keep changes testable and bounded

## Key subsystems

### Config Layer

#### CONTROL file format (`delfin/config.py`)
The `CONTROL` file is the central configuration. Key=value pairs, one per line.
Common keys: `functional`, `basis`, `charge`, `mult`, `PAL`, `maxcore`,
`dispersion`, `solvent`, `grid`, `tightscf`, `freq`, `opt`, `neb`, `scan`.
Key functions: `read_control_file()`, `validate_control_text()`,
`parse_control_text()`, `OCCUPIER_parser()`, `get_E_ref()`.
Multi-phase parsing: literal eval → sequence blocks → template defaults merge.
Handles OCCUPIER sequence blocks, GUPPY legacy syntax, ONIOM detection.

#### CONTROL Validator (`delfin/common/control_validator.py`)
Comprehensive field validation with type coercion and constraint checking.
`FieldSpec` class defines each field's type, range, default, required status,
and custom validator. ~100+ field specs covering functionals, basis sets,
solvents, dispersion, SCF, and all workflow keys.
`validate_control_config()` is the main entry: iterates all specs, applies
context-aware rules (e.g. ESD fields only validated when `ESD_modul=yes`,
OCCUPIER_tree skipped when method is `manually`).

### Orchestration Layer

#### Pipeline (`delfin/pipeline.py`)
Orchestrates phase-based pipeline: OCCUPIER → classic → manual → ESD → optional.
Key classes: `FileBundle`, `PipelineContext`, `SummaryResults`.
Key functions: `run_classic_phase()`, `run_manual_phase()`, `run_occuper_phase()`,
`run_esd_phase()`, `compute_summary()`.
Context object carries state through phases; handles E_ref computation,
Gibbs energy collection, redox analysis.

#### Job Scheduler (`delfin/parallel_classic_manually.py`)
Schedules dependent ORCA jobs with intelligent parallelism.
Key classes: `_WorkflowManager`, `WorkflowJob`, `WorkflowRunResult`.
Key functions: `execute_classic_workflows()`, `execute_manually_workflows()`,
`_populate_classic_jobs()`, `_populate_manual_jobs()`.
Job graph with dependencies; determines effective slot count; manages
PAL/maxcore resource blocks; tracks completion/failures/skipped.
Integrates with global `DynamicCorePool`; handles step-based workflows
(oxidation/reduction steps); updates MOINP blocks for continuation.

#### Global Scheduler (`delfin/global_scheduler.py`)
Thin wrapper around `_WorkflowManager` to aggregate ORCA jobs from multiple
workflow phases into a single run. `GlobalOrcaScheduler.add_jobs()` →
`.run()` → returns `WorkflowRunResult`. Logs job graph with per-job
cores (min/opt/max).

### Resource Management

#### Global Manager (`delfin/global_manager.py`)
Singleton `GlobalJobManager` ensures all workflows share a PAL ceiling and
prevent core over-allocation. Thread-safe with signal handling (SIGINT/SIGTERM).
Key functions: `initialize()`, `get_pool()`, `resolve_job_resources()`,
`shutdown()`, `ensure_signal_handlers()`.
Manages subprocess tracking for cleanup; adaptive job count based on PAL;
stdin interrupt monitor for tmux/screen compatibility.

#### Dynamic Core Pool (`delfin/dynamic_pool.py`)
ThreadPoolExecutor-based pool with intelligent allocation and timeout management.
Key classes: `DynamicCorePool`, `PoolJob`, `JobPriority`.
Thread-local job context via `get_current_job_id()`, `get_current_job_cores()`.
Priority queues; parent-child job tracking; starvation prevention; adaptive
timeouts (opt, frequency, general). Supports nested/spawned jobs.

### Error Handling

#### ORCA Recovery (`delfin/orca_recovery.py`)
Automatic error detection and recovery for ORCA failures.
Key classes: `OrcaErrorType` (enum), `OrcaErrorDetector`, `OrcaErrorRecoverer`.
Error pattern matching with priorities: TRAH crashes > LEANSCF > MPI.
Progressive fix escalation: MOREAD blocks, geometry adjustments.
Detects SCF non-convergence, segfaults, MPI crashes; applies continuation
with old.gbw + MOREAD; prevents infinite retry loops.

### Tool Integration

#### QM Runtime (`delfin/qm_runtime.py`)
Discovers and resolves quantum chemistry tool binaries (ORCA, Gaussian, xTB,
CREST, CENSO, ANMR, etc.).
Key classes: `ToolSpec`, `ResolvedTool`.
Key functions: `resolve_tool()`, `check_tools()`, `canonical_tool_name()`.
Multi-strategy resolution: env vars → `which` → system dirs → module patterns.
Caches results; supports tool aliases.

#### Runtime Setup (`delfin/runtime_setup.py`)
Configures runtime environment: symlinks, module discovery, environment
variable overrides. Handles xtb4stda Fortran 80-char path limit via
`_shorten_xtb4stda_path()`. Integrates with module systems (lmod, EasyBuild).
Key functions: `apply_runtime_environment()`, `collect_runtime_diagnostics()`,
`discover_orca_installations()`.

### Interface Layer

#### Dashboard (`delfin/dashboard/`)
Voilà/ipywidgets-based. Tabs: Home, Settings, Agent, Results.
`tab_agent.py` is the AI agent interface with multi-agent pipeline,
message queuing, per-role model routing, and cost tracking.

#### Local Backend (`delfin/dashboard/backend_local.py`)
JSON-based job queue with background worker thread. `LocalJobBackend` class.
Thread-based daemon worker; persistent JSON queue; resource budgeting with
oversubscribe factor; live CPU/RAM load monitoring.

#### SLURM Backend (`delfin/dashboard/backend_slurm.py`)
sbatch/squeue/scancel integration for cluster submission. `SlurmJobBackend`.
Site-specific profiles (bwUniCluster3, JUSTUS2); auto-detects tools; injects
environment via sbatch. Profile-based config (DELFIN_MODULES, DELFIN_STAGE_ORCA,
DELFIN_NODE_CORES).

#### Public API (`delfin/api.py`)
Programmatic entry points: `run()`, `prepare()`. Wraps CLI with Python kwargs
(cleanup, recalc, overwrite, define). Enables Python scripts to invoke DELFIN.

#### CLI (`delfin/cli.py`)
Command-line interface with lazy loading of heavy dependencies.
`_load_full_cli_dependencies()` imports ~30+ submodules only when needed.
Integrates all components: config parsing, pipeline, job management, reporting.

### Domain-Specific Modules

#### CO2 Coordinator (`delfin/co2/`)
Places CO2 molecules near metal complexes. Uses Fibonacci sphere sampling for
direction optimization, orientation scanning, and xTB/ORCA single-point ranking.
Key file: `CO2_Coordinator6.py`. CONTROL keys: `co2_coordination`, `place_axis`,
`co2_species_delta`.

#### Metal Complex Builder (`delfin/build_up_complex2.py`)
Parses SMILES → metals + ligands → enumerates coordination templates (octahedral,
tetrahedral, square_planar, etc.) → PSO swarm placement (6 DOF per ligand) →
optional BFGS refinement → optional xTB GFN2 re-ranking. CLI: `delfin-build2`.
Uses VdW radii for clash detection, Procrustes init for bidentate+ ligands.

#### ORCA Interface (`delfin/orca/`)
Generates ORCA input files, parses output (energies, geometries, frequencies).
Recovery logic in `orca_recovery.py` handles crashes and restarts.

### Testing
Tests in `tests/`. Run with `python -m pytest tests/ -x -q`.
NEVER run real ORCA/xTB computations in tests. Use mocks or small test data.

## Avoid

- generic advice detached from DELFIN modules
- widening already large monoliths without need
- changing runtime behavior without validation
- assuming local success implies SLURM success
- turning a small cycle into a broad rewrite
