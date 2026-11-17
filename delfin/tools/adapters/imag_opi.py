"""Registered step adapter for the OPI-backed IMAG input writer.

This adapter makes ``delfin.tools.imag_opi`` reachable through the step
registry (``run_step`` / ``list_adapters``): it writes the same ORCA ``.inp``
files the classic IMAG writer (``delfin.imag._candidate_job`` and the rerun
``_with_base(_with_geometry(...))``) writes, but re-emits them with OPI's
input classes.  The classic writer remains the workflow's driver — this step
is the parity-proof reference for its inputs (the OCCUPIER selection shape:
``tests/test_opi_occupier_parity.py`` proves the parity, an ``imag_opi`` step
is selected the same way ``occupier_opi`` is — registered adapter, lazily
imported OPI).

The step builds inputs only; it does not launch ORCA.  It requires OPI
(``requires_python=("opi",)``), imported lazily so the module imports without
OPI present.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import StepResult, StepStatus

logger = get_logger(__name__)


class ImagOpiAdapter(StepAdapter):
    """Re-emit IMAG's transformed ORCA input with OPI's input classes.

    Mirrors the classic IMAG writer's two outputs for one step: the
    single-point candidate (``delfin.imag._candidate_job``) and the
    re-optimisation rerun (``_with_base(_with_geometry(job, atoms)) + rest``).
    Parity is proven by ``tests/test_opi_imag_parity.py``.
    """

    name = "imag_opi"
    description = "OPI-backed ORCA input re-emission for IMAG candidate/rerun inputs"
    category = "dft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("orca_input",)
    requires_python = ("opi",)

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()
        try:
            from delfin.tools.imag_opi import (
                write_imag_candidate_opi_input,
                write_imag_rerun_opi_input,
            )

            mode = str(kwargs.get("mode", "candidate"))
            source_text = kwargs["source_text"]
            atoms = kwargs.get("atoms")
            if atoms is None:
                from delfin.imag import displaced_geometries

                hess_path = kwargs.get("hess_path")
                mode_i = int(kwargs.get("mode_index", 0))
                disps = displaced_geometries(
                    Path(hess_path), mode_i, 0.3)
                atoms = next(iter(disps.values()))
            base = str(kwargs.get("base") or "input")
            if mode == "rerun":
                inp_path = write_imag_rerun_opi_input(
                    work_dir,
                    source_text=source_text,
                    atoms=atoms,
                    base=base,
                )
            else:
                inp_path = write_imag_candidate_opi_input(
                    work_dir,
                    source_text=source_text,
                    atoms=atoms,
                    base=base,
                    optimise=bool(kwargs.get("optimise", False)),
                )
        except KeyError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"missing required argument {exc}",
            )
        except ValueError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=str(exc),
            )
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            output_file=inp_path,
            data={"input_file": str(inp_path)},
        )


register(ImagOpiAdapter())
