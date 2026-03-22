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

from delfin.tools._types import StepError, StepResult, StepStatus
from delfin.tools._runner import run_step, step_as_workflow_job
from delfin.tools._registry import list_steps, register
from delfin.tools.pipeline import Pipeline, PipelineResult, PipelineTemplate

__all__ = [
    "run_step",
    "step_as_workflow_job",
    "Pipeline",
    "PipelineTemplate",
    "PipelineResult",
    "StepResult",
    "StepStatus",
    "StepError",
    "list_steps",
    "register",
]
