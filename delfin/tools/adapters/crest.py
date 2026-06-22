"""CREST conformer search adapter."""

from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class CrestConformersAdapter(StepAdapter):
    name = "crest_conformers"
    description = "CREST conformer/isomer search"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        solvent = kwargs.get("solvent", "")
        ewin = kwargs.get("ewin", 6.0)

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        # Copy geometry into work_dir
        local_geom = self._copy_geometry_to_workdir(geometry, work_dir, "input.xyz")

        args = [
            str(local_geom),
            "--chrg", str(charge),
            "--uhf", str(mult - 1),
            "-T", str(cores),
            "--ewin", str(ewin),
        ]
        if solvent:
            args.extend(["--alpb", solvent])

        # Extra args passthrough
        extra_args = kwargs.get("extra_args", [])
        if extra_args:
            args.extend(extra_args)

        from delfin.qm_runtime import run_tool
        env = {"OMP_NUM_THREADS": str(cores), "OMP_STACKSIZE": "1G"}

        result = run_tool(
            "crest", args,
            cwd=str(work_dir),
            env=env,
            capture_output=True,
        )

        # Write log
        out_path = work_dir / "crest.out"
        out_path.write_text(result.stdout or "")

        best_xyz = work_dir / "crest_best.xyz"
        ensemble = work_dir / "crest_conformers.xyz"

        artifacts: dict[str, Path] = {}
        if ensemble.is_file():
            artifacts["ensemble"] = ensemble

        if result.returncode == 0 and best_xyz.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=best_xyz,
                output_file=out_path,
                artifacts=artifacts,
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            output_file=out_path,
            error=f"CREST failed (rc={result.returncode})",
        )


register(CrestConformersAdapter())
