"""Uniform adapter layer for chaining chemistry tools.

Single-step usage::

    from delfin.tools import run_step

    result = run_step("smiles_to_xyz", smiles="C")
    result = run_step("xtb_opt", geometry=result.geometry, charge=0, cores=4)
    result = run_step("orca_sp", geometry=result.geometry, charge=0, cores=8)
    print(result.data["electronic_energy_Eh"])

Pipeline usage::

    from delfin.tools import Pipeline

    pipe = Pipeline("my_workflow")
    pipe.add("smiles_to_xyz", smiles="CCO")
    pipe.add("xtb_opt", charge=0)
    pipe.add("orca_opt", charge=0, bang="! B3LYP def2-SVP")
    pipe.add("orca_freq", charge=0, bang="! B3LYP def2-SVP")
    results = pipe.run(cores=8)
    print(results.summary())
"""

from delfin.tools._types import ErrorKind, StepError, StepResult, StepStatus
from delfin.tools._runner import run_step, step_as_workflow_job
from delfin.tools._registry import list_steps, register
from delfin.tools._spec import DataKeySpec, ParamSpec, StepContract
from delfin.tools._catalog import catalog, compatible_successors, describe
from delfin.tools._serialize import PipelineSerializationError, register_callable
from delfin.tools._application import (
    Application,
    OutputSpec,
    list_applications,
    register_application,
)
from delfin.tools.pipeline import Pipeline, PipelineInserter, PipelineResult, PipelineTemplate

__all__ = [
    "run_step",
    "step_as_workflow_job",
    "Pipeline",
    "PipelineInserter",
    "PipelineTemplate",
    "PipelineResult",
    "StepResult",
    "StepStatus",
    "StepError",
    "ErrorKind",
    "list_steps",
    "register",
    # Self-describing contract + discovery
    "ParamSpec",
    "DataKeySpec",
    "StepContract",
    "describe",
    "catalog",
    "compatible_successors",
    # Serialization (blocks as data)
    "register_callable",
    "PipelineSerializationError",
    # Applications (workflows with a contract) — facade in delfin.tools.platform
    "Application",
    "OutputSpec",
    "register_application",
    "list_applications",
]
