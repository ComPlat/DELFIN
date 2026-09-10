"""Imaginary frequency elimination adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec


class ImagFixAdapter(StepAdapter):
    name = "imag_fix"
    description = "Eliminate imaginary frequencies via iterative IMAG optimization"
    produces_geometry = True
    category = "dft"
    params = (
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("mult", "int", required=True, description="Spin multiplicity (2S+1)"),
        ParamSpec("solvent", "str", required=True, description="Implicit solvent name"),
        ParamSpec("metals", "list", required=True, description="Metal elements in the system"),
        ParamSpec("main_basisset", "str", required=True, description="Main basis set"),
        ParamSpec("metal_basisset", "str", required=True, description="Basis set for metals"),
        ParamSpec("hess_file", "path", required=True,
                  description="Path to the upstream .hess (hessian capability)"),
        ParamSpec("broken_sym", "bool", default=False, description="Use broken symmetry"),
        ParamSpec("step_name", "str", default="imag_fix", description="Label for the IMAG run"),
    )
    consumes = ("geometry",)   # the hessian is supplied via the explicit hess_file param
    requires_binaries = ("orca",)

    def validate_params(self, **kwargs: Any) -> None:
        for key in ("charge", "mult", "solvent", "metals", "main_basisset", "metal_basisset"):
            if key not in kwargs:
                raise ValueError(f"'{key}' parameter is required")
        if "hess_file" not in kwargs:
            raise ValueError("'hess_file' parameter is required (path to .hess)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        """IMAG on the calculation that wrote ``hess_file``.

        The Hessian's calculation is ``<name>.inp`` / ``<name>.out`` beside
        ``<name>.hess`` -- what the ORCA adapters write.  IMAG re-runs that
        input, so it needs it; this used to hand IMAG the geometry file in the
        output's place, IMAG found no imaginary mode in an xyz, and the step
        reported success without having done anything.
        """
        import shutil

        from delfin.imag import _pipeline_run_orca, _referenced_files, eliminate_imaginary_modes

        start = time.monotonic()
        hess_src = Path(kwargs["hess_file"])
        inp_src, out_src = hess_src.with_suffix(".inp"), hess_src.with_suffix(".out")
        missing = [p.name for p in (hess_src, inp_src, out_src) if not p.is_file()]
        if missing:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"IMAG needs the Hessian's calculation beside it; missing: {', '.join(missing)}",
            )

        # Work on copies: the upstream step's files stay what that step produced.
        work_dir.mkdir(parents=True, exist_ok=True)
        for src in (inp_src, out_src, hess_src):
            shutil.copy2(src, work_dir / src.name)
        for name in _referenced_files(inp_src.read_text(encoding="utf-8")):
            if (hess_src.parent / name).is_file():
                shutil.copy2(hess_src.parent / name, work_dir / Path(name).name)

        result = eliminate_imaginary_modes(
            label=kwargs.get("step_name", "imag_fix"),
            input_path=work_dir / inp_src.name,
            output_path=work_dir / out_src.name,
            config=kwargs.get("config") or {},
            run_orca=_pipeline_run_orca,
            pal=cores,
        )
        refined = work_dir / f"{hess_src.stem}.xyz"
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=refined if refined.is_file() else geometry,
            output_file=work_dir / out_src.name,
            data={"n_imaginary": len(result.remaining), "imag_rounds": result.rounds,
                  "imag_reason": result.reason or ""},
            artifacts={"hess": work_dir / hess_src.name},
        )


register(ImagFixAdapter())
