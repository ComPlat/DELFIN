"""Workflow execution engine: job scheduling, dependency resolution, ORCA dispatch.

Note: ``occupier`` submodule is NOT imported eagerly to avoid circular imports
(occupier → copy_helpers → delfin.occupier → parallel_classic_manually).
Use ``from delfin.workflows.engine.occupier import ...`` when needed.
"""

from delfin.workflows.engine.classic import (
    WorkflowJob,
    WorkflowRunResult,
    _WorkflowManager,
    execute_classic_workflows,
    execute_manually_workflows,
    normalize_parallel_token,
    estimate_parallel_width,
    determine_effective_slots,
)
from delfin.workflows.engine.scheduler import GlobalOrcaScheduler

__all__ = [
    "WorkflowJob",
    "WorkflowRunResult",
    "_WorkflowManager",
    "execute_classic_workflows",
    "execute_manually_workflows",
    "normalize_parallel_token",
    "estimate_parallel_width",
    "determine_effective_slots",
    "GlobalOrcaScheduler",
]
