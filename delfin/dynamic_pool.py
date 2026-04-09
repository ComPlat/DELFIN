"""Backward-compat shim. Canonical location: delfin.workflows.scheduling.pool."""

from delfin.workflows.scheduling.pool import *  # noqa: F401,F403
from delfin.workflows.scheduling.pool import (  # noqa: F401
    _set_current_job_cores,
    _set_current_job_id,
)
