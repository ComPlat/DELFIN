"""xTB SOLVATOR adapter (explicit solvation shell via ORCA-xTB)."""

from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class XtbSolvatorAdapter(StepAdapter):
    name = "xtb_solvator"
    description = "Add explicit solvation shell via xTB SOLVATOR"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "solvent" not in kwargs:
            raise ValueError("'solvent' parameter is required (e.g. 'water', 'acetonitrile')")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        solvent = kwargs["solvent"]
        n_solv = kwargs.get("n_solvent", 12)
        method = kwargs.get("method", "GFN2-xTB")

        # Copy geometry, strip XYZ header to get coordinate block
        local_geom = self._copy_geometry_to_workdir(geometry, work_dir, "input.xyz")
        with open(local_geom, "r", encoding="utf-8") as fh:
            lines = fh.readlines()
        # Strip 2-line XYZ header if present
        if len(lines) > 2:
            try:
                int(lines[0].strip())
                coord_lines = lines[2:]
            except ValueError:
                coord_lines = lines
        else:
            coord_lines = lines

        # Build ORCA input
        inp_file = work_dir / "XTB_SOLVATOR.inp"
        out_file = work_dir / "output_XTB_SOLVATOR.out"
        inp_text = (
            f"!{method} ALPB({solvent})\n"
            f"%SOLVATOR NSOLV {n_solv} END\n"
            f"%pal nprocs {cores} end\n"
            f"*xyz {charge} {mult}\n"
        )
        inp_text += "".join(coord_lines)
        if not inp_text.endswith("\n"):
            inp_text += "\n"
        inp_text += "*\n"
        inp_file.write_text(inp_text, encoding="utf-8")

        from delfin.orca import run_orca
        run_orca(str(inp_file), str(out_file), working_dir=str(work_dir))

        solvated_xyz = work_dir / "XTB_SOLVATOR.solvator.xyz"
        if solvated_xyz.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=solvated_xyz,
                data={"solvent": solvent, "n_solvent": n_solv},
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error="XTB_SOLVATOR did not produce solvated geometry",
        )


register(XtbSolvatorAdapter())
