"""Registered step adapter for the OPI-backed CO2 coordinator input builder.

This adapter makes ``delfin.tools.co2_coordinator_opi`` reachable through
the step registry (``run_step`` / ``list_adapters``): it writes the same
ORCA ``.inp`` files the classic CO2 Coordinator writes, but with OPI's
input classes.  The classic coordinator remains the workflow's driver —
this step is the parity-proof reference for its inputs.

The step builds inputs only; it does not launch ORCA.  It requires OPI
(``requires_python=("opi",)``), imported lazily so the module imports
without OPI present.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import StepResult, StepStatus

logger = get_logger(__name__)


class Co2CoordinatorOpiAdapter(StepAdapter):
    """Build the CO2-coordinator ORCA inputs with OPI's input classes.

    Mirrors ``Co2Coordinator6``'s two input writers (orientation
    single-points and the relaxed distance scan) into a registered step.
    The inputs written are equivalent to the classic writer's — the parity
    test ``tests/test_opi_co2_coordinator_parity.py`` proves it.
    """

    name = "co2_coordinator_opi"
    description = "OPI-backed ORCA input builder for the CO2 coordinator"
    category = "dft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("orca_input",)
    requires_python = ("opi",)

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        **kwargs: Any,
    ) -> StepResult:
        del geometry  # classic writers take atoms + paths as keyword args
        start = time.monotonic()
        try:
            from delfin.tools.co2_coordinator_opi import (
                write_co2_sp_opi_input,
            )

            outdir = kwargs.get("outdir", work_dir / "calc")
            tag = kwargs.get("tag", "calc")
            inp_path = write_co2_sp_opi_input(
                atoms=kwargs["atoms"],
                xyz_path=str(kwargs.get("xyz_path", outdir / "input.xyz")),
                outdir=str(outdir),
                orca_keywords=kwargs.get("orca_keywords", "PBE0 def2-SVP SP"),
                broken_sym=kwargs.get("broken_sym", ""),
                charge=kwargs.get("charge", 0),
                multiplicity=kwargs.get("mult", 1),
                PAL=kwargs.get("PAL", 4),
                maxcore=kwargs.get("maxcore", 1000),
                tag=tag,
                metal_symbol=kwargs.get("metal_symbol"),
                metal_basis=kwargs.get("metal_basis"),
                control_args=kwargs.get("control_args"),
                qmmm_range=kwargs.get("qmmm_range"),
                co2_indices=kwargs.get("co2_indices"),
            )
        except KeyError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"missing required argument {exc}",
            )
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            output_file=inp_path,
            data={"input_file": str(inp_path)},
        )


register(Co2CoordinatorOpiAdapter())
