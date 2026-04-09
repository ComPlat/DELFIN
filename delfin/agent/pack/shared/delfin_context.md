# DELFIN Context

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

## High-value modules

- `delfin/cli.py`: CLI orchestration
- `delfin/pipeline.py`: workflow orchestration
- `delfin/parallel_classic_manually.py`: job scheduler
- `delfin/config.py`: CONTROL parsing and validation
- `delfin/qm_runtime.py`: tool resolution
- `delfin/runtime_setup.py`: runtime diagnostics
- `delfin/orca_recovery.py`: retry and recovery
- `delfin/dashboard/backend_local.py`: local execution
- `delfin/dashboard/backend_slurm.py`: SLURM execution
- `delfin/dashboard/tab_agent.py`: AI agent interface
- `delfin/build_up_complex2.py`: metal complex builder (PSO swarm)
- `delfin/smiles_converter.py`: SMILES → 3D coordinates

## Architectural pressure points

- Monolithic orchestration files (pipeline.py, parallel_classic_manually.py)
- CLI semantics leaking into the API layer
- Mutable config dicts carrying runtime state
- Scheduler correctness and PAL/resource safety
- Recovery logic complexity

## Preferred direction

- Explicit state and contracts over hidden config mutation
- Keep public surfaces stable
- Protect local and HPC behavior equally
- Make failures diagnosable, not silent
- Keep changes testable and bounded

## Testing

Tests in `tests/`. Run with `python -m pytest tests/ -x -q`.
NEVER run real ORCA/xTB computations. Use mocks or test data.

## Avoid

- Generic advice detached from DELFIN modules
- Widening monoliths without need
- Changing runtime behavior without validation
- Assuming local success implies SLURM success
- Turning a small cycle into a broad rewrite
