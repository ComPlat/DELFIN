"""ML potential adapters (single-point, optimize, rank)."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class MlpSinglePointAdapter(StepAdapter):
    name = "mlp_single_point"
    description = "Single-point energy/forces via ML potential (ANI-2x, MACE, CHGNet, ...)"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "backend" not in kwargs:
            raise ValueError("'backend' parameter is required (e.g. 'ani2x', 'mace_off', 'chgnet')")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry (XYZ) is required")

        from ase.io import read as ase_read
        from delfin.mlp_tools.calculators import single_point

        atoms = ase_read(str(geometry))
        backend = kwargs["backend"]
        result_dict = single_point(
            atoms,
            backend=backend,
            device=kwargs.get("device", "cpu"),
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
        )

        # forces is ndarray, convert for JSON
        data = {
            "energy_eV": result_dict["energy"],
            "backend": result_dict["backend"],
            "n_atoms": result_dict["n_atoms"],
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=data,
        )


class MlpOptimizeAdapter(StepAdapter):
    name = "mlp_optimize"
    description = "Geometry optimization via ML potential"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "backend" not in kwargs:
            raise ValueError("'backend' parameter is required (e.g. 'ani2x', 'mace_off', 'chgnet')")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry (XYZ) is required")

        from ase.io import read as ase_read, write as ase_write
        from delfin.mlp_tools.calculators import optimize_geometry

        atoms = ase_read(str(geometry))
        backend = kwargs["backend"]
        result_dict = optimize_geometry(
            atoms,
            backend=backend,
            device=kwargs.get("device", "cpu"),
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            fmax=kwargs.get("fmax", 0.05),
            steps=kwargs.get("steps", 200),
            optimizer=kwargs.get("optimizer", "LBFGS"),
        )

        out_xyz = work_dir / "optimized.xyz"
        ase_write(str(out_xyz), atoms)

        data = {
            "converged": result_dict["converged"],
            "energy_eV": result_dict["energy"],
            "n_steps": result_dict["n_steps"],
            "backend": result_dict["backend"],
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz,
            data=data,
        )


class MlpRankAdapter(StepAdapter):
    name = "mlp_rank"
    description = "Rank conformer ensemble by ML potential energy"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "backend" not in kwargs:
            raise ValueError("'backend' parameter is required")
        if "ensemble" not in kwargs:
            raise ValueError("'ensemble' parameter is required (path to multi-structure XYZ)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from ase.io import read as ase_read
        from delfin.mlp_tools.calculators import rank_structures

        ensemble_path = Path(kwargs["ensemble"])
        structures = ase_read(str(ensemble_path), index=":")
        backend = kwargs["backend"]

        ranking = rank_structures(
            structures,
            backend=backend,
            device=kwargs.get("device", "cpu"),
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
        )

        data = {
            "ranking": [{"index": idx, "energy_eV": e} for idx, e in ranking],
            "n_structures": len(structures),
            "backend": backend,
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=data,
        )


register(MlpSinglePointAdapter())
register(MlpOptimizeAdapter())
register(MlpRankAdapter())
