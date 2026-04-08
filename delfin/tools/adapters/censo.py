"""CENSO conformer sorting adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class CensoSortAdapter(StepAdapter):
    name = "censo_sort"
    description = "CENSO conformer ensemble sorting and ranking"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "ensemble" not in kwargs:
            raise ValueError("'ensemble' parameter is required (path to conformer ensemble XYZ)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.analysis_tools.censo_wrapper import run_censo

        ensemble = kwargs["ensemble"]
        result_dict = run_censo(
            ensemble,
            charge=kwargs.get("charge", 0),
            uhf=kwargs.get("uhf", 0),
            solvent=kwargs.get("solvent", ""),
            functional=kwargs.get("functional", "r2scan-3c"),
            prescreening=kwargs.get("prescreening", True),
            screening=kwargs.get("screening", True),
            optimization=kwargs.get("optimization", False),
            nprocs=cores,
            omp=kwargs.get("omp", min(2, cores)),
            threshold=kwargs.get("threshold", 4.0),
            timeout=kwargs.get("timeout", 7200),
            working_dir=str(work_dir),
        )

        success = result_dict.get("returncode", 1) == 0
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir,
            start,
            data=result_dict,
            error=None if success else f"CENSO failed (rc={result_dict.get('returncode')})",
        )


register(CensoSortAdapter())
