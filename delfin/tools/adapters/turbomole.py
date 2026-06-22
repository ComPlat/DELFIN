"""Turbomole adapter — geometry optimization, single-point, frequencies.

Requires Turbomole binaries (ridft, jobex, aoforce, etc.) in PATH or
configured via environment variables.

Usage::

    run_step("turbomole_sp", geometry="mol.xyz", charge=0, method="b3-lyp",
             basis="def2-SVP", cores=4)

    run_step("turbomole_opt", geometry="mol.xyz", charge=0, method="ri-bp86",
             basis="def2-SVP")
"""

from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


def _run_tm_command(cmd: str, work_dir: Path, cores: int,
                    timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    """Run a Turbomole command in work_dir."""
    import os
    env = os.environ.copy()
    env["PARNODES"] = str(cores)
    env["OMP_NUM_THREADS"] = str(cores)

    return subprocess.run(
        cmd, shell=True, cwd=work_dir, env=env,
        capture_output=True, text=True, timeout=timeout,
    )


def _setup_coord(geometry: Path, work_dir: Path) -> Path:
    """Convert XYZ to Turbomole coord format."""
    coord_path = work_dir / "coord"
    # Use x2t if available, else manual conversion
    x2t = shutil.which("x2t")
    if x2t:
        result = subprocess.run(
            [x2t, str(geometry)], cwd=work_dir,
            capture_output=True, text=True,
        )
        coord_path.write_text(result.stdout)
    else:
        # Fallback: simple xyz → coord conversion
        lines = geometry.read_text().strip().splitlines()
        n_atoms = int(lines[0])
        coord_lines = ["$coord"]
        bohr = 1.8897259886  # Å → Bohr
        for line in lines[2:2 + n_atoms]:
            parts = line.split()
            x, y, z = float(parts[1]) * bohr, float(parts[2]) * bohr, float(parts[3]) * bohr
            coord_lines.append(f"  {x:20.14f}  {y:20.14f}  {z:20.14f}  {parts[0].lower()}")
        coord_lines.append("$end")
        coord_path.write_text("\n".join(coord_lines) + "\n")
    return coord_path


def _setup_control(work_dir: Path, method: str, basis: str, charge: int,
                   mult: int, cores: int, **kwargs: Any) -> None:
    """Generate Turbomole control file via define or manually."""
    define = shutil.which("define")
    if define:
        # Use define with input script
        define_input = f"""
a coord
*
no
b all {basis}
*
eht
y
{charge}
y
*
dft
on
func {method}
*
ri
on
*
*
"""
        proc = subprocess.run(
            [define], input=define_input, cwd=work_dir,
            capture_output=True, text=True, timeout=60,
        )
    else:
        # Minimal manual control file
        control = f"""$title
$operating system unix
$symmetry c1
$coord file=coord
$atoms
basis={basis}
$charge {charge}
$dft
  functional {method}
$ricore 500
$maxcor 500
$scfconv 7
$end
"""
        (work_dir / "control").write_text(control)


def _extract_energy_tm(work_dir: Path) -> Optional[float]:
    """Extract energy from Turbomole energy file."""
    energy_file = work_dir / "energy"
    if not energy_file.is_file():
        return None
    try:
        lines = energy_file.read_text().strip().splitlines()
        for line in reversed(lines):
            line = line.strip()
            if line and not line.startswith("$"):
                parts = line.split()
                if len(parts) >= 2:
                    return float(parts[1])
    except Exception:
        pass
    return None


def _coord_to_xyz(work_dir: Path) -> Optional[Path]:
    """Convert Turbomole coord back to XYZ."""
    coord_path = work_dir / "coord"
    if not coord_path.is_file():
        return None
    try:
        bohr = 1.8897259886
        lines = coord_path.read_text().strip().splitlines()
        atoms = []
        for line in lines:
            line = line.strip()
            if line.startswith("$") or not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                x = float(parts[0]) / bohr
                y = float(parts[1]) / bohr
                z = float(parts[2]) / bohr
                elem = parts[3].capitalize()
                atoms.append(f"{elem}  {x:.8f}  {y:.8f}  {z:.8f}")
        if atoms:
            xyz_path = work_dir / "output.xyz"
            xyz_path.write_text(f"{len(atoms)}\nTurbomole output\n" + "\n".join(atoms) + "\n")
            return xyz_path
    except Exception:
        pass
    return None


class TurbomoleSpAdapter(StepAdapter):
    name = "turbomole_sp"
    description = "Turbomole single-point energy (DFT/RI)"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "b3-lyp")
        basis = kwargs.get("basis", "def2-SVP")
        timeout = kwargs.get("timeout")

        _setup_coord(geometry, work_dir)
        _setup_control(work_dir, method, basis, charge, mult, cores, **kwargs)

        # Run ridft (RI-DFT) or dscf (conventional)
        cmd = "ridft" if kwargs.get("ri", True) else "dscf"
        proc = _run_tm_command(cmd, work_dir, cores, timeout=timeout)

        success = proc.returncode == 0
        energy = _extract_energy_tm(work_dir)
        data = {}
        if energy is not None:
            data["energy_Eh"] = energy

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            data=data,
            error=None if success else f"Turbomole {cmd} failed: {proc.stderr[:500]}",
        )


class TurbomoleOptAdapter(StepAdapter):
    name = "turbomole_opt"
    description = "Turbomole geometry optimization"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "b3-lyp")
        basis = kwargs.get("basis", "def2-SVP")
        max_cycles = kwargs.get("max_cycles", 200)
        timeout = kwargs.get("timeout")

        _setup_coord(geometry, work_dir)
        _setup_control(work_dir, method, basis, charge, mult, cores, **kwargs)

        cmd = f"jobex -ri -c {max_cycles}"
        proc = _run_tm_command(cmd, work_dir, cores, timeout=timeout)

        success = proc.returncode == 0
        energy = _extract_energy_tm(work_dir)
        geom_out = _coord_to_xyz(work_dir)

        data = {}
        if energy is not None:
            data["energy_Eh"] = energy

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            geometry=geom_out,
            data=data,
            error=None if success else f"Turbomole jobex failed: {proc.stderr[:500]}",
        )


class TurbomoleFreqAdapter(StepAdapter):
    name = "turbomole_freq"
    description = "Turbomole frequency calculation (aoforce)"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        charge = kwargs["charge"]
        mult = kwargs.get("mult", 1)
        method = kwargs.get("method", "b3-lyp")
        basis = kwargs.get("basis", "def2-SVP")
        timeout = kwargs.get("timeout")

        _setup_coord(geometry, work_dir)
        _setup_control(work_dir, method, basis, charge, mult, cores, **kwargs)

        # First SP, then aoforce
        _run_tm_command("ridft", work_dir, cores, timeout=timeout)
        proc = _run_tm_command("aoforce", work_dir, cores, timeout=timeout)

        success = proc.returncode == 0
        energy = _extract_energy_tm(work_dir)

        data = {}
        if energy is not None:
            data["energy_Eh"] = energy

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            data=data,
            error=None if success else f"Turbomole aoforce failed: {proc.stderr[:500]}",
        )


register(TurbomoleSpAdapter())
register(TurbomoleOptAdapter())
register(TurbomoleFreqAdapter())
