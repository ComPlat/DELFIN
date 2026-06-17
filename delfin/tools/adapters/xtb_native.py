"""Native xTB building blocks — the standalone ``xtb`` binary, full feature set.

These wrap the free, open-source Grimme ``xtb`` program directly (not via ORCA),
exposing all of its capabilities as composable building blocks:

* methods — GFN0/1/2-xTB and the GFN-FF force field
* run types — single point, optimization, Hessian/thermo (ohess), MD, Fukui
* implicit solvation — ALPB or GBSA with any xTB solvent name
* the usual knobs — charge, spin (via multiplicity → ``--uhf``), accuracy,
  electronic temperature, optimization level, plus a raw ``extra_args`` escape
  hatch for anything not surfaced as a named parameter.

They run the binary in a subprocess with ``cwd=work_dir`` (thread-safe, no
``os.chdir``), mirroring the other engine-wrapping adapters. The ORCA-driven
``xtb_opt`` / ``xtb_goat`` adapters are left untouched; these are additive and
named distinctly (``xtb_sp`` / ``xtb_optimize`` / ``xtb_hessian`` / ``xtb_md`` /
``xtb_fukui``). ``xtb`` is free software, so it may be installed automatically.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, List, Optional, Sequence

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import StepResult, StepStatus

# GFN-xTB level → xtb CLI flags.
_METHOD_FLAGS = {
    "gfn0": ["--gfn", "0"],
    "gfn1": ["--gfn", "1"],
    "gfn2": ["--gfn", "2"],
    "gfnff": ["--gfnff"],
}
_OPT_LEVELS = ("crude", "sloppy", "loose", "normal", "tight", "vtight", "extreme")

# Shared, declarative parameters across the native xtb blocks.
_COMMON_PARAMS = (
    ParamSpec("charge", "int", required=True, description="Molecular charge"),
    ParamSpec("mult", "int", default=1, description="Spin multiplicity (2S+1) → --uhf (mult-1)"),
    ParamSpec("method", "str", default="gfn2", enum=tuple(_METHOD_FLAGS),
              description="xTB level: GFN0/1/2-xTB or the GFN-FF force field"),
    ParamSpec("solvent", "str", description="Implicit solvent name (xTB/ALPB list), e.g. water"),
    ParamSpec("solvent_model", "str", default="alpb", enum=("alpb", "gbsa"),
              description="Implicit solvation model"),
    ParamSpec("accuracy", "float", description="xtb --acc accuracy (smaller = tighter)"),
    ParamSpec("etemp", "float", unit="K", description="Electronic temperature (--etemp)"),
    ParamSpec("extra_args", "list", description="Raw extra xtb CLI args (escape hatch)"),
)


def build_xtb_cmd(
    input_name: str,
    *,
    method: str = "gfn2",
    run: Optional[str] = None,
    charge: int = 0,
    uhf: int = 0,
    solvent: Optional[str] = None,
    solvent_model: str = "alpb",
    accuracy: Optional[float] = None,
    etemp: Optional[float] = None,
    opt_level: Optional[str] = None,
    cores: int = 1,
    extra_args: Sequence[str] = (),
) -> List[str]:
    """Assemble the ``xtb`` command line (pure, so it is unit-testable).

    ``run`` selects the job: ``None`` single point, ``"opt"``, ``"ohess"``
    (opt + Hessian + thermo), ``"hess"``, ``"md"`` (``--omd``), ``"fukui"``.
    """
    cmd: List[str] = ["xtb", input_name]
    cmd += _METHOD_FLAGS.get(str(method).lower(), ["--gfn", "2"])
    cmd += ["--chrg", str(int(charge)), "--uhf", str(int(uhf))]
    if run == "opt":
        cmd += ["--opt", opt_level or "normal"]
    elif run == "ohess":
        cmd += ["--ohess", opt_level or "normal"]
    elif run == "hess":
        cmd += ["--hess"]
    elif run == "md":
        cmd += ["--omd"]
    elif run == "fukui":
        cmd += ["--vfukui"]
    if solvent:
        cmd += [f"--{solvent_model}", solvent]
    if accuracy is not None:
        cmd += ["--acc", str(accuracy)]
    if etemp is not None:
        cmd += ["--etemp", str(etemp)]
    cmd += ["-P", str(int(cores))]
    cmd += [str(a) for a in extra_args]
    return cmd


def _parse_float(pattern: str, text: str) -> Optional[float]:
    m = re.search(pattern, text, re.IGNORECASE)
    return float(m.group(1)) if m else None


def _parse_energy(text: str) -> Optional[float]:
    return _parse_float(r"total energy\s+(-?\d+\.\d+)\s*Eh", text)


def _parse_gap(text: str) -> Optional[float]:
    return _parse_float(r"HOMO-LUMO gap\s+(-?\d+\.\d+)\s*eV", text)


def _parse_thermo(text: str) -> dict:
    out: dict = {}
    g = _parse_float(r"total free energy\s+(-?\d+\.\d+)\s*Eh", text)
    if g is not None:
        out["free_energy_Eh"] = g
    z = _parse_float(r"zero point energy\s+(-?\d+\.\d+)\s*Eh", text)
    if z is not None:
        out["zpve_Eh"] = z
    m = re.search(r"#\s*imaginary freq\.\s+(\d+)", text, re.IGNORECASE)
    if m:
        out["n_imaginary"] = int(m.group(1))
    return out


def _run_xtb(
    adapter: StepAdapter,
    work_dir: Path,
    geometry: Optional[Path],
    start: float,
    *,
    run: Optional[str],
    cores: int,
    kwargs: dict,
    md_control: Optional[str] = None,
):
    """Copy the geometry in, run xtb in *work_dir*, return (ok, out_path, text, fail).

    ``fail`` is a ready :class:`StepResult` on a setup/exec error (the adapter
    returns it directly), else ``None``.
    """
    def _fail(msg: str) -> StepResult:
        return adapter._make_result(adapter.name, StepStatus.FAILED, work_dir, start, error=msg)

    if geometry is None:
        return False, None, "", _fail("geometry is required")
    if shutil.which("xtb") is None:
        return False, None, "", _fail("xtb binary not found on PATH")
    inp = work_dir / "xtb_input.xyz"
    shutil.copy2(geometry, inp)
    args = list(kwargs.get("extra_args") or [])
    if md_control is not None:
        (work_dir / "md.inp").write_text(md_control)
        args = ["--input", "md.inp", *args]
    cmd = build_xtb_cmd(
        inp.name, method=kwargs.get("method", "gfn2"), run=run,
        charge=int(kwargs["charge"]), uhf=int(kwargs.get("mult", 1)) - 1,
        solvent=kwargs.get("solvent"), solvent_model=kwargs.get("solvent_model", "alpb"),
        accuracy=kwargs.get("accuracy"), etemp=kwargs.get("etemp"),
        opt_level=kwargs.get("opt_level"), cores=cores, extra_args=args,
    )
    env = {**os.environ, "OMP_NUM_THREADS": str(cores), "MKL_NUM_THREADS": str(cores)}
    out_path = work_dir / f"{adapter.name}.out"
    try:
        proc = subprocess.run(cmd, cwd=str(work_dir), env=env, text=True,
                              capture_output=True, timeout=kwargs.get("timeout", 3600))
    except (subprocess.TimeoutExpired, OSError) as exc:
        return False, out_path, "", _fail(f"xtb run failed: {exc}")
    out_path.write_text((proc.stdout or "") + (proc.stderr or ""))
    return proc.returncode == 0, out_path, out_path.read_text(), None


class _NativeXtb(StepAdapter):
    """Common validation for the native xtb blocks."""

    category = "semiempirical"
    consumes = ("geometry",)
    requires_binaries = ("xtb",)

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")


class XtbSpAdapter(_NativeXtb):
    name = "xtb_sp"
    description = "Native xTB single-point energy (GFN0/1/2/FF, ALPB/GBSA solvation)"
    produces_geometry = False
    params = _COMMON_PARAMS
    produces = ("qc_output",)
    data_keys = (DataKeySpec("energy_Eh", "float", "Eh", "Total xTB energy"),
                 DataKeySpec("homo_lumo_gap_eV", "float", "eV", "HOMO-LUMO gap"))

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        ok, out_path, text, fail = _run_xtb(self, work_dir, geometry, start,
                                            run=None, cores=cores, kwargs=kwargs)
        if fail is not None:
            return fail
        data: dict = {}
        e = _parse_energy(text)
        if e is not None:
            data["energy_Eh"] = e
        gap = _parse_gap(text)
        if gap is not None:
            data["homo_lumo_gap_eV"] = gap
        return self._make_result(
            self.name, StepStatus.SUCCESS if ok else StepStatus.FAILED, work_dir, start,
            output_file=out_path, data=data,
            error=None if ok else "xtb single point failed")


class XtbOptimizeAdapter(_NativeXtb):
    name = "xtb_optimize"
    description = "Native xTB geometry optimization (--opt, any GFN level / solvent)"
    produces_geometry = True
    params = _COMMON_PARAMS + (
        ParamSpec("opt_level", "str", default="normal", enum=_OPT_LEVELS,
                  description="Optimization convergence level"),
    )
    produces = ("qc_output",)
    data_keys = (DataKeySpec("energy_Eh", "float", "Eh", "Optimized xTB energy"),)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        ok, out_path, text, fail = _run_xtb(self, work_dir, geometry, start,
                                            run="opt", cores=cores, kwargs=kwargs)
        if fail is not None:
            return fail
        xyz = work_dir / "xtbopt.xyz"
        data: dict = {}
        e = _parse_energy(text)
        if e is not None:
            data["energy_Eh"] = e
        ok = ok and xyz.is_file()
        return self._make_result(
            self.name, StepStatus.SUCCESS if ok else StepStatus.FAILED, work_dir, start,
            geometry=xyz if xyz.is_file() else None, output_file=out_path, data=data,
            error=None if ok else "xtb optimization failed or no xtbopt.xyz")


class XtbHessianAdapter(_NativeXtb):
    name = "xtb_hessian"
    description = "Native xTB opt + Hessian + thermochemistry (--ohess)"
    produces_geometry = True
    params = _COMMON_PARAMS + (
        ParamSpec("opt_level", "str", default="normal", enum=_OPT_LEVELS,
                  description="Optimization convergence level"),
    )
    produces = ("qc_output", "hessian")
    data_keys = (
        DataKeySpec("energy_Eh", "float", "Eh", "Total xTB energy"),
        DataKeySpec("free_energy_Eh", "float", "Eh", "Total free energy (RRHO)"),
        DataKeySpec("zpve_Eh", "float", "Eh", "Zero-point vibrational energy"),
        DataKeySpec("n_imaginary", "int", "", "Number of imaginary frequencies"),
    )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        ok, out_path, text, fail = _run_xtb(self, work_dir, geometry, start,
                                            run="ohess", cores=cores, kwargs=kwargs)
        if fail is not None:
            return fail
        xyz = work_dir / "xtbopt.xyz"
        hess = work_dir / "hessian"
        data: dict = {}
        e = _parse_energy(text)
        if e is not None:
            data["energy_Eh"] = e
        data.update(_parse_thermo(text))
        artifacts = {"hess": hess} if hess.is_file() else {}
        return self._make_result(
            self.name, StepStatus.SUCCESS if ok else StepStatus.FAILED, work_dir, start,
            geometry=xyz if xyz.is_file() else None, output_file=out_path,
            data=data, artifacts=artifacts,
            error=None if ok else "xtb ohess failed")


class XtbMdAdapter(_NativeXtb):
    name = "xtb_md"
    description = "Native xTB molecular dynamics (--omd), trajectory as an ensemble"
    produces_geometry = True
    params = _COMMON_PARAMS + (
        ParamSpec("temperature_K", "float", default=300.0, unit="K", description="MD temperature"),
        ParamSpec("time_ps", "float", default=10.0, unit="ps", description="MD simulation time"),
        ParamSpec("step_fs", "float", default=2.0, unit="fs", description="MD time step"),
        ParamSpec("dump_fs", "float", default=50.0, unit="fs", description="Trajectory dump interval"),
    )
    produces = ("qc_output", "ensemble")

    def _md_control(self, kwargs: dict) -> str:
        return (
            "$md\n"
            f"   temp={kwargs.get('temperature_K', 300.0)}\n"
            f"   time={kwargs.get('time_ps', 10.0)}\n"
            f"   step={kwargs.get('step_fs', 2.0)}\n"
            f"   dump={kwargs.get('dump_fs', 50.0)}\n"
            "$end\n"
        )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        ok, out_path, text, fail = _run_xtb(self, work_dir, geometry, start,
                                            run="md", cores=cores, kwargs=kwargs,
                                            md_control=self._md_control(kwargs))
        if fail is not None:
            return fail
        traj = work_dir / "xtb.trj"
        xyz = work_dir / "xtbopt.xyz"
        artifacts = {"ensemble": traj} if traj.is_file() else {}
        return self._make_result(
            self.name, StepStatus.SUCCESS if ok else StepStatus.FAILED, work_dir, start,
            geometry=xyz if xyz.is_file() else None, output_file=out_path,
            artifacts=artifacts,
            error=None if ok else "xtb MD failed")


class XtbFukuiAdapter(_NativeXtb):
    name = "xtb_fukui"
    description = "Native xTB Fukui indices (--vfukui) for reactivity analysis"
    produces_geometry = False
    category = "analysis"
    params = _COMMON_PARAMS
    produces = ("qc_output",)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        ok, out_path, text, fail = _run_xtb(self, work_dir, geometry, start,
                                            run="fukui", cores=cores, kwargs=kwargs)
        if fail is not None:
            return fail
        return self._make_result(
            self.name, StepStatus.SUCCESS if ok else StepStatus.FAILED, work_dir, start,
            output_file=out_path,
            error=None if ok else "xtb Fukui calculation failed")


register(XtbSpAdapter())
register(XtbOptimizeAdapter())
register(XtbHessianAdapter())
register(XtbMdAdapter())
register(XtbFukuiAdapter())
