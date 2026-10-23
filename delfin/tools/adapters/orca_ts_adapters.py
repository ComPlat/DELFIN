"""OR1 phase 1: ORCA transition-state, IRC and frequency-verification adapters.

Built on the vocabulary of ``delfin/tools/adapters/orca.py`` (the ORCA
input builder, ``run_orca``) and ``delfin/imag.py`` (imaginary-mode
parsing).  No new ORCA keywords are invented: OptTS, AnFreq/NumFreq and
IrcForward/IrcBackward are ORCA's own keywords, and the imaginary-mode
check reuses DELFIN's IMAG parser.

Usage::

    from delfin.tools import run_step

    result = run_step("orca_ts_opt", geometry="ts_guess.xyz",
                      charge=0, method="B3LYP", basis="def2-SVP",
                      hessian_mode="AnFreq")
"""
from __future__ import annotations

import shutil
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec

# ORCA's own Hessian strategies for an OptTS job.
_HAISSIAN_MODES = ("AnFreq", "NumFreq")


def _read_xyz_coords(xyz_path: Path) -> str:
    """Coordinate lines of an XYZ file, no header (same contract as xtb.py)."""
    lines = xyz_path.read_text().strip().splitlines()
    if len(lines) < 3:
        return "\n".join(lines)
    try:
        int(lines[0].strip())
        return "\n".join(lines[2:])
    except ValueError:
        return "\n".join(lines)


def _build_orca_ts_input(
    bang: str,
    charge: int,
    mult: int,
    cores: int,
    maxcore: int,
    coords: str,
    extra_blocks: str = "",
) -> str:
    """A full ORCA input in the shape _build_full_orca_input writes."""
    parts = [f"! {bang}\n", f"%maxcore {maxcore}\n", f"%pal nprocs {cores} end\n"]
    if extra_blocks:
        parts.append(extra_blocks.rstrip() + "\n")
    parts.append(f"*xyz {charge} {mult}\n{coords}\n*\n")
    return "".join(parts)


def _run_orca_job(inp_path: Path, out_path: Path, work_dir: Path, timeout=None) -> bool:
    """One ORCA run through the sanctioned runner."""
    from delfin.orca import run_orca
    return run_orca(str(inp_path), str(out_path), timeout=timeout, working_dir=work_dir)


def _collect_artifacts(work_dir: Path, base_name: str) -> dict:
    """gbw/hess artifacts of a finished ORCA job, as orca.py does."""
    artifacts: dict[str, Path] = {}
    for key, ext in (("gbw", ".gbw"), ("hess", ".hess")):
        f = work_dir / f"{base_name}{ext}"
        if f.is_file():
            artifacts[key] = f
    return artifacts


class OrcaTsOptAdapter(StepAdapter):
    name = "orca_ts_opt"
    description = "ORCA OptTS transition-state optimization with an explicit Hessian strategy"
    produces_geometry = True
    category = "dft"
    params = (
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("mult", "int", default=1, description="Spin multiplicity (2S+1)"),
        ParamSpec("method", "str", default="B3LYP", description="Functional"),
        ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
        ParamSpec("hessian_mode", "str", default="AnFreq", enum=_HAISSIAN_MODES,
                  description="Hessian strategy for the TS search (AnFreq | NumFreq)"),
        ParamSpec("maxcore", "int", default=1000, unit="MB", description="Memory per core"),
        ParamSpec("base_name", "str", default="calc", description="Base name for input/output files"),
    )
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "hessian")
    requires_binaries = ("orca",)

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")
        if kwargs.get("hessian_mode", "AnFreq") not in _HAISSIAN_MODES:
            raise ValueError(
                f"hessian_mode must be one of {', '.join(_HAISSIAN_MODES)}; "
                f"got {kwargs.get('hessian_mode')!r}"
            )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        self.validate_params(**kwargs)
        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "B3LYP")
        basis = kwargs.get("basis", "def2-SVP")
        hessian_mode = kwargs.get("hessian_mode", "AnFreq")
        maxcore = kwargs.get("maxcore", 1000)
        base_name = kwargs.get("base_name", "calc")
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        coords = _read_xyz_coords(geometry)
        bang = f"{method} {basis} OptTS {hessian_mode}"
        inp_content = _build_orca_ts_input(bang, charge, mult, cores, maxcore, coords)
        inp_path = work_dir / f"{base_name}.inp"
        out_path = work_dir / f"{base_name}.out"
        inp_path.write_text(inp_content)

        success = _run_orca_job(inp_path, out_path, work_dir, timeout=kwargs.get("timeout"))

        geom_out = work_dir / f"{base_name}.xyz"
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            geometry=geom_out if (success and geom_out.is_file()) else None,
            output_file=out_path if out_path.is_file() else None,
            data={},
            artifacts=_collect_artifacts(work_dir, base_name),
            error=None if success else "ORCA OptTS calculation failed",
        )


class OrcaTsFreqCheckAdapter(StepAdapter):
    """Frequency calculation with an exactly-one-imaginary-mode verdict.

    A TS is a first-order saddle point.  The verdict is the step's
    result: 1 imaginary mode -> SUCCESS with ``n_imag``/``imag_mode``;
    0 or >= 2 -> FAILED, naming the count, so a "TS" with two imaginary
    modes can never pass as a transition state.
    """
    name = "orca_ts_freq_check"
    description = "ORCA frequency run with exactly-one-imaginary-mode TS verification"
    produces_geometry = False
    category = "dft"
    params = (
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("mult", "int", default=1, description="Spin multiplicity (2S+1)"),
        ParamSpec("method", "str", default="B3LYP", description="Functional"),
        ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
        ParamSpec("maxcore", "int", default=1000, unit="MB", description="Memory per core"),
        ParamSpec("base_name", "str", default="calc", description="Base name for input/output files"),
    )
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "hessian")
    requires_binaries = ("orca",)

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        self.validate_params(**kwargs)
        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "B3LYP")
        basis = kwargs.get("basis", "def2-SVP")
        maxcore = kwargs.get("maxcore", 1000)
        base_name = kwargs.get("base_name", "calc")
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        coords = _read_xyz_coords(geometry)
        bang = f"{method} {basis} FREQ"
        inp_content = _build_orca_ts_input(bang, charge, mult, cores, maxcore, coords)
        inp_path = work_dir / f"{base_name}.inp"
        out_path = work_dir / f"{base_name}.out"
        inp_path.write_text(inp_content)

        success = _run_orca_job(inp_path, out_path, work_dir, timeout=kwargs.get("timeout"))

        # --- the verdict ---
        modes: list = []
        if out_path.is_file():
            from delfin.imag import collect_imaginary_modes
            modes = collect_imaginary_modes(str(out_path))
        n_imag = len(modes)
        data: dict[str, Any] = {"n_imag": n_imag}
        if modes:
            data["imag_mode"] = modes[0][1]

        verdict_ok = success and n_imag == 1
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if verdict_ok else StepStatus.FAILED,
            work_dir, start,
            output_file=out_path if out_path.is_file() else None,
            data=data,
            artifacts=_collect_artifacts(work_dir, base_name),
            error=None if verdict_ok else (
                f"not a TS: {n_imag} imaginary mode(s) found "
                f"(expected exactly 1); ORCA success={success}"
            ),
        )


class OrcaTsIrcAdapter(StepAdapter):
    """IRC from a TS, forward and backward, continuing from the TS orbitals.

    Both directions are two separate ORCA runs (IrcForward / IrcBackward),
    each restarting from the TS's .gbw via MOREAD.  Without the gbw the
    step refuses to run — an IRC without the TS orbitals would silently
    restart from scratch and trace a different path.
    """
    name = "orca_ts_irc"
    description = "ORCA intrinsic reaction coordinate in both directions from a TS"
    produces_geometry = True
    category = "dft"
    params = (
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("mult", "int", default=1, description="Spin multiplicity (2S+1)"),
        ParamSpec("method", "str", default="B3LYP", description="Functional"),
        ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
        ParamSpec("moread", "path", required=True, description="TS orbital (.gbw) file for MOREAD"),
        ParamSpec("maxcore", "int", default=1000, unit="MB", description="Memory per core"),
    )
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "geometry")
    requires_binaries = ("orca",)

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")
        if not kwargs.get("moread"):
            raise ValueError(
                "moread is required: IRC continues from the TS orbitals (.gbw); "
                "without it ORCA would silently restart from scratch"
            )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "B3LYP")
        basis = kwargs.get("basis", "def2-SVP")
        maxcore = kwargs.get("maxcore", 1000)

        moread_src = Path(kwargs["moread"])
        if not moread_src.suffix:
            moread_src = moread_src.with_suffix(".gbw")
        if not moread_src.is_file():
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"TS gbw not found: {moread_src}")
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry is required")

        coords = _read_xyz_coords(geometry)
        endpoints: list[str] = []
        for direction in ("forward", "backward"):
            keyword = {"forward": "IrcForward", "backward": "IrcBackward"}[direction]
            bang = f"{method} {basis} {keyword} MOREAD"
            inp_content = _build_orca_ts_input(bang, charge, mult, cores, maxcore, coords)
            inp_path = work_dir / f"irc_{direction}.inp"
            out_path = work_dir / f"irc_{direction}.out"
            inp_path.write_text(inp_content)
            shutil.copy2(moread_src, work_dir / moread_src.name)
            success = _run_orca_job(inp_path, out_path, work_dir, timeout=kwargs.get("timeout"))
            if not success:
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    output_file=out_path if out_path.is_file() else None,
                    data={"irc_endpoints": endpoints or None},
                    error=f"IRC {direction} run failed",
                )
            end_xyz = work_dir / f"irc_{direction}.xyz"
            if end_xyz.is_file():
                endpoints.append(str(end_xyz))
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir, start,
            geometry=Path(endpoints[-1]) if endpoints else None,
            data={"irc_endpoints": endpoints or None},
            error=None,
        )


register(OrcaTsOptAdapter())
register(OrcaTsFreqCheckAdapter())
register(OrcaTsIrcAdapter())
