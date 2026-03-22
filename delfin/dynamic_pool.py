"""Backward-compat shim. Canonical location: delfin.workflows.scheduling.pool"""
from delfin.workflows.scheduling.pool import *  # noqa: F401,F403
from delfin.workflows.scheduling.pool import (  # noqa: F401
    _job_context,
    _set_current_job_id,
)
