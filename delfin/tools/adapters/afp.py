"""AFP (Absorption-Fluorescence-Phosphorescence) plot adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class AfpPlotAdapter(StepAdapter):
    name = "afp_plot"
    description = "Generate combined absorption/fluorescence/phosphorescence spectrum plot"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "workspace_dir" not in kwargs:
            raise ValueError("'workspace_dir' parameter is required (path to ESD workspace)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.afp_plot import create_afp_plot

        workspace_dir = Path(kwargs["workspace_dir"])
        output_png = work_dir / "afp_spectrum.png"
        fwhm = kwargs.get("fwhm", 50.0)

        result_path = create_afp_plot(
            workspace_dir=workspace_dir,
            output_png=output_png,
            fwhm=fwhm,
        )

        if result_path is not None and Path(result_path).is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                data={"plot": str(result_path)},
                artifacts={"afp_plot": Path(result_path)},
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error="AFP plot generation failed",
        )


register(AfpPlotAdapter())
