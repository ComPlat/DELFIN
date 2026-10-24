"""xTB SOLVATOR adapter (explicit solvation shell via ORCA-xTB)."""

from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._keys import key


class XtbSolvatorAdapter(StepAdapter):
    name = "xtb_solvator"
    description = "Add explicit solvation shell via xTB SOLVATOR"
    produces_geometry = True
    category = "solvation"
    consumes = ("geometry",)
    requires_binaries = ("orca",)
    params = (
        key("solvent", required=True),
        key("charge", default=0), key("mult"),
        ParamSpec("n_solvent", "int", default=12, description="Number of solvent molecules"),
        ParamSpec("method", "str", default="GFN2-xTB", description="xTB method"),
    )
    data_keys = (DataKeySpec("solvent", "str"), DataKeySpec("n_solvent", "int"))

    def validate_params(self, **kwargs: Any) -> None:
        if "solvent" not in kwargs:
            raise ValueError("'solvent' parameter is required (e.g. 'water', 'acetonitrile')")

        # ── Declared physical limits (SV1) ─────────────────────────────
        # SOLVATOR is an xTB feature (ORCA manual, %SOLVATOR block): the
        # method must be an xTB Hamiltonian. A DFT name here would only
        # surface as an ORCA error after queue time.
        method = str(kwargs.get("method", "GFN2-xTB")).strip()
        if method.upper() not in ("GFN0-XTB", "GFN1-XTB", "GFN2-XTB"):
            raise ValueError(
                f"method {method!r} is not an xTB Hamiltonian; %SOLVATOR "
                "(ORCA manual, %SOLVATOR block) supports GFN0-xTB, "
                "GFN1-xTB or GFN2-xTB")

        solvent = kwargs.get("solvent")
        if not isinstance(solvent, str) or not solvent.strip():
            raise ValueError(
                "solvent must be a non-empty string naming an ALPB/ALPHA "
                "solvent (ORCA manual §4.5), e.g. 'water' or 'acetonitrile'")

        mult = kwargs.get("mult", 1)
        if isinstance(mult, bool) or not isinstance(mult, int) or mult < 1:
            raise ValueError(
                f"mult {mult!r} is invalid; spin multiplicity (2S+1) must "
                "be an integer >= 1")

        charge = kwargs.get("charge", 0)
        if isinstance(charge, bool) or not isinstance(charge, int):
            raise ValueError(
                f"charge {charge!r} is invalid; the molecular charge must "
                "be an integer")

        n_solv = kwargs.get("n_solvent", 12)
        if isinstance(n_solv, bool) or not isinstance(n_solv, int) or n_solv < 1:
            raise ValueError(
                f"n_solvent {n_solv!r} is invalid; NSOLV counts explicit "
                "solvent molecules and must be an integer >= 1")

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
