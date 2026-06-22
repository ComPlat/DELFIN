"""Imaginary frequency elimination adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class ImagFixAdapter(StepAdapter):
    name = "imag_fix"
    description = "Eliminate imaginary frequencies via iterative IMAG optimization"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        for key in ("charge", "mult", "solvent", "metals", "main_basisset", "metal_basisset"):
            if key not in kwargs:
                raise ValueError(f"'{key}' parameter is required")
        if "hess_file" not in kwargs:
            raise ValueError("'hess_file' parameter is required (path to .hess)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        # Copy input files into work_dir
        local_geom = self._copy_geometry_to_workdir(geometry, work_dir, "input.xyz")
        hess_src = Path(kwargs["hess_file"])
        local_hess = work_dir / hess_src.name
        if not local_hess.exists():
            import shutil
            shutil.copy2(hess_src, local_hess)

        config = kwargs.get("config", {})
        config.setdefault("PAL", cores)

        from delfin.imag import run_IMAG

        run_IMAG(
            input_file=str(local_geom),
            hess_file=str(local_hess),
            charge=kwargs["charge"],
            multiplicity=kwargs["mult"],
            solvent=kwargs["solvent"],
            metals=kwargs["metals"],
            config=config,
            main_basisset=kwargs["main_basisset"],
            metal_basisset=kwargs["metal_basisset"],
            broken_sym=kwargs.get("broken_sym", False),
            step_name=kwargs.get("step_name", "imag_fix"),
            pal_override=cores,
        )

        # IMAG writes the corrected geometry back to input_file
        if local_geom.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=local_geom,
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error="IMAG did not produce output geometry",
        )


register(ImagFixAdapter())
