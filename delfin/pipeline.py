"""Backward-compat shim. Canonical location: delfin.workflows.pipeline"""

from delfin.workflows.pipeline import *  # noqa: F401,F403
from delfin.workflows.pipeline import (  # noqa: F401
    _is_truthy_token,
    _run_guppy_for_smiles,
    _skip_xtb_goat_after_guppy,
)
