"""Backward-compat shim. Canonical location: delfin.workflows.engine.classic."""

from delfin.workflows.engine.classic import *  # noqa: F401,F403
from delfin.workflows.engine.classic import (  # noqa: F401
    _WorkflowManager,
    _LocalWorkflowPool,
    _update_pal_block,
    _add_moinp_block,
    _verify_orca_output,
    _populate_classic_jobs,
    _populate_manual_jobs,
    _extract_manual_broken_sym,
    _parse_int,
    _normalize_tokens,
    _step_enabled,
    _ACTIVE_MANAGERS,
    _ACTIVE_MANAGERS_LOCK,
)
