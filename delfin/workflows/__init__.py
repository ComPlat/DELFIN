"""DELFIN Workflows -- modular computational chemistry workflow framework.

Subpackages
-----------
scheduling/  -- job scheduling, core pool, resource management
engine/      -- workflow execution, dependency resolution, ORCA dispatch
contrib/     -- built-in workflow wrappers (TADF-xTB, hyperpol, ...)
"""

from delfin.workflows.types import Workflow
from delfin.workflows.registry import (
    get as get_workflow,
    list_workflows,
    register as register_workflow,
)

__all__ = [
    "Workflow",
    "get_workflow",
    "list_workflows",
    "register_workflow",
]
