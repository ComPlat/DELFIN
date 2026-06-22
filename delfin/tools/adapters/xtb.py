"""xTB optimization adapters (via ORCA's built-in xTB methods)."""

from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


def _read_xyz_coords(xyz_path: Path) -> str:
    """Read an XYZ file and return only the coordinate lines (no header)."""
    lines = xyz_path.read_text().strip().splitlines()
    if len(lines) < 3:
        return "\n".join(lines)
    # Standard XYZ: line 0 = atom count, line 1 = comment, rest = coords
    try:
        int(lines[0].strip())
        return "\n".join(lines[2:])
    except ValueError:
        # Already raw coords (no header)
        return "\n".join(lines)


def _build_orca_xtb_input(
    method: str,
    keyword: str,
    charge: int,
    mult: int,
    cores: int,
    maxcore: int,
    coords: str,
) -> str:
    """Build an ORCA input file for xTB calculations."""
    return (
        f"!{method} {keyword}\n"
        f"%maxcore {maxcore}\n"
        f"%pal nprocs {cores} end\n"
        f"*xyz {charge} {mult}\n"
        f"{coords}\n"
        f"*\n"
    )


class XtbOptAdapter(StepAdapter):
    name = "xtb_opt"
    description = "GFN2-xTB geometry optimization via ORCA"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "XTB2")
        maxcore = kwargs.get("maxcore", 1000)

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        coords = _read_xyz_coords(geometry)
        inp_content = _build_orca_xtb_input(method, "OPT", charge, mult, cores, maxcore, coords)

        inp_path = work_dir / "xtb_opt.inp"
        out_path = work_dir / "xtb_opt.out"
        xyz_path = work_dir / "xtb_opt.xyz"
        inp_path.write_text(inp_content)

        from delfin.orca import run_orca
        success = run_orca(str(inp_path), str(out_path), working_dir=work_dir)

        if success and xyz_path.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=xyz_path,
                output_file=out_path,
                data={"method": method},
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            output_file=out_path if out_path.is_file() else None,
            error="xTB optimization failed or output XYZ not produced",
        )


class XtbGoatAdapter(StepAdapter):
    name = "xtb_goat"
    description = "GOAT conformer search via ORCA's built-in xTB"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "XTB2")
        maxcore = kwargs.get("maxcore", 1000)

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        coords = _read_xyz_coords(geometry)
        inp_content = _build_orca_xtb_input(method, "GOAT", charge, mult, cores, maxcore, coords)

        inp_path = work_dir / "xtb_goat.inp"
        out_path = work_dir / "xtb_goat.out"
        xyz_path = work_dir / "xtb_goat.globalminimum.xyz"
        inp_path.write_text(inp_content)

        from delfin.orca import run_orca
        success = run_orca(str(inp_path), str(out_path), working_dir=work_dir)

        if success and xyz_path.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=xyz_path,
                output_file=out_path,
                data={"method": method},
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            output_file=out_path if out_path.is_file() else None,
            error="GOAT search failed or global minimum XYZ not produced",
        )


register(XtbOptAdapter())
register(XtbGoatAdapter())
