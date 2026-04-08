"""Genarris crystal structure prediction adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class GenarrisAdapter(StepAdapter):
    name = "genarris_csp"
    description = "Random crystal structure generation via Genarris"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "config" not in kwargs:
            raise ValueError("'config' parameter is required (path to Genarris INI config)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.csp_tools.genarris_wrapper import run_genarris_cli

        config_path = Path(kwargs["config"])
        np = kwargs.get("np", cores)
        mpirun = kwargs.get("mpirun", "mpirun")
        timeout = kwargs.get("timeout", None)

        try:
            proc = run_genarris_cli(
                config_path,
                np=np,
                cwd=str(work_dir),
                mpirun=mpirun,
                capture_output=True,
                timeout=timeout,
            )
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"Genarris failed: {exc}",
            )

        success = proc.returncode == 0
        data = {"returncode": proc.returncode}

        # Collect generated structure files
        cif_files = sorted(work_dir.glob("*.cif"))
        if cif_files:
            data["n_structures"] = len(cif_files)
            data["structure_files"] = [str(f.name) for f in cif_files]

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            data=data,
            error=None if success else f"Genarris failed (rc={proc.returncode})",
        )


register(GenarrisAdapter())
