"""ASE (Atomic Simulation Environment) adapters.

Provides access to any ASE calculator (ORCA, xTB, VASP, GPAW, …) and
ASE's optimization/dynamics infrastructure.

Usage::

    run_step("ase_optimize", geometry="mol.xyz", calculator="xtb",
             fmax=0.05, optimizer="BFGS")

    run_step("ase_single_point", geometry="mol.xyz",
             calculator="orca", orca_simple_input="B3LYP def2-SVP")

    run_step("ase_md", geometry="mol.xyz", calculator="xtb",
             temperature_K=300, timestep_fs=1.0, steps=1000)
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


def _get_calculator(name: str, work_dir: Path, cores: int, **kwargs: Any):
    """Instantiate an ASE calculator by name."""
    name_lower = name.lower()

    if name_lower == "xtb":
        from ase.calculators.xtb import XTB
        return XTB(method=kwargs.get("method", "GFN2-xTB"))

    elif name_lower == "orca":
        from ase.calculators.orca import ORCA
        simple_input = kwargs.get("orca_simple_input", "B3LYP def2-SVP")
        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        return ORCA(
            label="calc",
            orcasimpleinput=simple_input,
            orcablocks=f"%pal nprocs {cores} end\n%maxcore {kwargs.get('maxcore', 1000)}",
            charge=charge, mult=mult,
            directory=str(work_dir),
        )

    elif name_lower in ("emt", "eam"):
        from ase.calculators.emt import EMT
        return EMT()

    elif name_lower == "lj":
        from ase.calculators.lj import LennardJones
        return LennardJones()

    elif name_lower == "mace":
        from mace.calculators import MACECalculator
        model = kwargs.get("model", "medium")
        return MACECalculator(model_paths=model, device="cpu")

    elif name_lower == "ani":
        import torchani
        return torchani.models.ANI2x().ase()

    else:
        raise ValueError(
            f"Unknown ASE calculator '{name}'. "
            f"Supported: xtb, orca, emt, lj, mace, ani"
        )


class AseOptimizeAdapter(StepAdapter):
    name = "ase_optimize"
    description = "Geometry optimization via ASE (any calculator: xTB, ORCA, MACE, …)"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "calculator" not in kwargs:
            raise ValueError("'calculator' parameter is required (e.g. 'xtb', 'orca')")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        from ase.io import read, write
        from ase.optimize import BFGS, LBFGS, FIRE

        atoms = read(str(geometry))
        calc = _get_calculator(kwargs["calculator"], work_dir, cores, **kwargs)
        atoms.calc = calc

        optimizer_name = kwargs.get("optimizer", "BFGS").upper()
        optimizers = {"BFGS": BFGS, "LBFGS": LBFGS, "FIRE": FIRE}
        opt_class = optimizers.get(optimizer_name, BFGS)

        fmax = kwargs.get("fmax", 0.05)
        max_steps = kwargs.get("max_steps", 500)

        log_path = work_dir / "opt.log"
        traj_path = work_dir / "opt.traj"

        try:
            opt = opt_class(atoms, logfile=str(log_path), trajectory=str(traj_path))
            converged = opt.run(fmax=fmax, steps=max_steps)
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"ASE optimization failed: {e}")

        out_xyz = work_dir / "optimized.xyz"
        write(str(out_xyz), atoms)

        energy = atoms.get_potential_energy()
        data = {"energy_eV": energy, "converged": bool(converged),
                "n_steps": opt.nsteps, "calculator": kwargs["calculator"]}

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz, output_file=log_path,
            data=data,
        )


class AseSinglePointAdapter(StepAdapter):
    name = "ase_single_point"
    description = "Single-point energy/forces via ASE (any calculator)"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "calculator" not in kwargs:
            raise ValueError("'calculator' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        from ase.io import read

        atoms = read(str(geometry))
        calc = _get_calculator(kwargs["calculator"], work_dir, cores, **kwargs)
        atoms.calc = calc

        try:
            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"ASE single-point failed: {e}")

        data = {
            "energy_eV": float(energy),
            "max_force_eV_A": float(forces.max()),
            "calculator": kwargs["calculator"],
        }

        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=data)


class AseMdAdapter(StepAdapter):
    name = "ase_md"
    description = "Molecular dynamics via ASE (any calculator)"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "calculator" not in kwargs:
            raise ValueError("'calculator' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        from ase.io import read, write
        from ase.md.langevin import Langevin
        from ase import units

        atoms = read(str(geometry))
        calc = _get_calculator(kwargs["calculator"], work_dir, cores, **kwargs)
        atoms.calc = calc

        temperature_K = kwargs.get("temperature_K", 300)
        timestep_fs = kwargs.get("timestep_fs", 1.0)
        steps = kwargs.get("steps", 1000)
        traj_path = work_dir / "md.traj"

        try:
            dyn = Langevin(atoms, timestep_fs * units.fs, temperature_K=temperature_K,
                           friction=0.01, logfile=str(work_dir / "md.log"),
                           trajectory=str(traj_path))
            dyn.run(steps)
        except Exception as e:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"ASE MD failed: {e}")

        out_xyz = work_dir / "final.xyz"
        write(str(out_xyz), atoms)

        data = {"steps": steps, "temperature_K": temperature_K,
                "final_energy_eV": float(atoms.get_potential_energy()),
                "calculator": kwargs["calculator"]}

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz, data=data,
        )


register(AseOptimizeAdapter())
register(AseSinglePointAdapter())
register(AseMdAdapter())
