"""OPI-backed ORCA geometry-optimisation / frequency adapter for delfin.tools.

Thin registry adapter: the implementation (input assembly with OPI's typed
blocks, launch through ``delfin.orca.run_orca``, convergence-gated read-back)
lives in ``delfin.tools.orca_opt_opi``; this module only re-exports it so the
step registry discovers it under the name ``orca_opt_opi`` automatically.
"""

from delfin.tools.orca_opt_opi import (  # noqa: F401
    OrcaOptOpiAdapter,
    run_irc,
    run_neb,
    run_opt_freq,
    run_opt_only,
    run_scan,
    validate_opt_method,
)

__all__ = ["OrcaOptOpiAdapter"]
