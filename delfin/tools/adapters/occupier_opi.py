"""Registered step adapter for the OPI-backed OCCUPIER input builder.

This adapter makes ``delfin.tools.occupier_opi`` reachable through the
step registry (``run_step`` / ``list_adapters``): it writes the same
ORCA ``.inp`` file the classic OCCUPIER writer
(``delfin.occupier.read_and_modify_file_OCCUPIER``) writes, but with
OPI's input classes.  The classic writer remains the workflow's driver —
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


class OccupierOpiAdapter(StepAdapter):
    """Build the OCCUPIER FoB ORCA input with OPI's input classes.

    Mirrors the classic OCCUPIER writer's per-fragment-of-basis input
    (``delfin.occupier.read_and_modify_file_OCCUPIER``).  The input
    written is equivalent to the classic writer's — the parity test
    ``tests/test_opi_occupier_parity.py`` proves it.
    """

    name = "occupier_opi"
    description = "OPI-backed ORCA input builder for OCCUPIER fragments-of-basis"
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
            from delfin.tools.occupier_opi import write_occupier_fob_opi_input

            inp_path = write_occupier_fob_opi_input(
                work_dir,
                config=kwargs.get("config") or {},
                charge=kwargs.get("charge", 0),
                multiplicity=kwargs.get("multiplicity", 1),
                PAL=kwargs.get("PAL", 4),
                maxcore=kwargs.get("maxcore", 1000),
                solvent=kwargs.get("solvent"),
                found_metals=kwargs.get("found_metals") or [],
                broken_sym=kwargs.get("broken_sym", ""),
                from_index=kwargs.get("from_index", 1),
                base_name=kwargs.get("base_name", "input"),
                atoms=kwargs.get("atoms"),
                xyz_path=kwargs.get("xyz_path", geometry),
                qmmm_range=kwargs.get("qmmm_range"),
                pass_wavefunction=kwargs.get("pass_wavefunction", ""),
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


register(OccupierOpiAdapter())
