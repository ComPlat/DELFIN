"""Guppy SMILES sampling adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class GuppySamplingAdapter(StepAdapter):
    name = "guppy_sampling"
    description = "SMILES-based conformer/isomer sampling with xTB ranking"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "smiles" not in kwargs:
            raise ValueError("'smiles' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        smiles = kwargs["smiles"]
        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        max_isomers = kwargs.get("max_isomers", 50)

        # Build argv for guppy_sampling.main()
        argv = [
            smiles,
            "--charge", str(charge),
            "--mult", str(mult),
            "--pal", str(cores),
            "--output-dir", str(work_dir),
            "--max-isomers", str(max_isomers),
        ]
        extra_args = kwargs.get("extra_args", [])
        if extra_args:
            argv.extend(extra_args)

        from delfin.guppy_sampling import main as guppy_main
        rc = guppy_main(argv)

        best_xyz = work_dir / "best.xyz"
        # Guppy writes various output files; look for the best one
        if not best_xyz.is_file():
            # Try common alternative names
            for candidate in ("guppy_best.xyz", "ranked_1.xyz"):
                alt = work_dir / candidate
                if alt.is_file():
                    best_xyz = alt
                    break

        if rc == 0 and best_xyz.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=best_xyz,
                data={"smiles": smiles},
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error=f"Guppy sampling failed (rc={rc})",
        )


register(GuppySamplingAdapter())
