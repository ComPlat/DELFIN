"""Registered step adapter for the OPI-backed ESD input writer.

This adapter makes ``delfin.tools.esd_opi`` reachable through the step
registry (``run_step`` / ``list_adapters``): it writes the same ORCA
``.inp`` files the classic ESD writer (``delfin.esd_input_generator``)
writes, but with OPI's input classes.  The classic writer remains the
workflow's driver — this step is the parity-proof reference for its
inputs, exactly like ``delfin.tools.adapters.occupier_opi`` for OCCUPIER.

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


class EsdOpiAdapter(StepAdapter):
    """Build ESD state/rate ORCA inputs with OPI's input classes.

    Mirrors the classic ESD writer's inputs
    (``delfin.esd_input_generator.create_state_input`` and its rate
    counterparts).  The written inputs are equivalent to the classic
    writer's — the parity test ``tests/test_opi_esd_parity.py`` proves it
    over the repo's CONTROL matrix.
    """

    name = "esd_opi"
    description = "OPI-backed ORCA input builder for ESD state and rate jobs"
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
        start = time.monotonic()
        try:
            from delfin.tools.esd_opi import (
                write_esd_rate_opi_input,
                write_esd_state_opi_input,
            )

            esd_dir = Path(kwargs.get("esd_dir") or work_dir)
            out_dir = Path(kwargs.get("out_dir") or work_dir)
            kind = kwargs.get("kind")
            if kind:
                inp_path = write_esd_rate_opi_input(
                    esd_dir=esd_dir,
                    out_dir=out_dir,
                    kind=str(kind),
                    pair=kwargs.get("pair"),
                    charge=kwargs.get("charge", 0),
                    solvent=kwargs.get("solvent"),
                    metals=kwargs.get("metals") or [],
                    main_basisset=kwargs.get("main_basisset", ""),
                    metal_basisset=kwargs.get("metal_basisset"),
                    config=kwargs.get("config") or {},
                    trootssl=kwargs.get("trootssl"),
                    initial_state=kwargs.get("initial_state", "S1"),
                )
            else:
                inp_path = write_esd_state_opi_input(
                    esd_dir=esd_dir,
                    out_dir=out_dir,
                    state=kwargs.get("state", "S0"),
                    charge=kwargs.get("charge", 0),
                    solvent=kwargs.get("solvent"),
                    metals=kwargs.get("metals") or [],
                    main_basisset=kwargs.get("main_basisset", ""),
                    metal_basisset=kwargs.get("metal_basisset"),
                    config=kwargs.get("config") or {},
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


register(EsdOpiAdapter())
