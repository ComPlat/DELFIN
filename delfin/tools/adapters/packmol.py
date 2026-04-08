"""Packmol solvation adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class PackmolSolvateAdapter(StepAdapter):
    name = "packmol_solvate"
    description = "Build solvation box via Packmol"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "solvent_xyz" not in kwargs:
            raise ValueError("'solvent_xyz' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry (solute) is required")

        from delfin.analysis_tools.packmol_wrapper import solvate

        output_xyz = work_dir / "solvated.xyz"
        result_path = solvate(
            solute_xyz=str(geometry),
            solvent_xyz=str(kwargs["solvent_xyz"]),
            n_solvent=kwargs.get("n_solvent", 50),
            box_side=kwargs.get("box_side", 20.0),
            tolerance=kwargs.get("tolerance", 2.0),
            output_xyz=str(output_xyz),
            working_dir=str(work_dir),
        )

        if result_path is not None and Path(result_path).is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=Path(result_path),
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error="Packmol solvation failed",
        )


register(PackmolSolvateAdapter())
