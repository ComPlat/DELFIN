"""Packmol solvation adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec


class PackmolSolvateAdapter(StepAdapter):
    name = "packmol_solvate"
    description = "Build solvation box via Packmol"
    produces_geometry = True
    category = "solvation"
    consumes = ("geometry",)
    requires_binaries = ("packmol",)
    params = (
        ParamSpec("solvent_xyz", "path", required=True, description="Solvent molecule XYZ"),
        ParamSpec("n_solvent", "int", default=50),
        ParamSpec("box_side", "float", default=20.0, unit="A"),
        ParamSpec("tolerance", "float", default=2.0, unit="A"),
    )

    def validate_params(self, **kwargs: Any) -> None:
        if "solvent_xyz" not in kwargs:
            raise ValueError("'solvent_xyz' parameter is required")

        # ── Declared physical limits (SV1) ─────────────────────────────
        # Packmol, Martínez/Andrade/Birgin/Martínez, J. Comput. Chem.
        # 2009, 30, 2157-2164 (DOI 10.1002/jcc.21224): `tolerance` is a
        # hard lower bound on the inter-molecular distance and the box
        # must be able to hold the solute plus at least one solvent
        # molecule.
        n_solvent = kwargs.get("n_solvent", 50)
        if isinstance(n_solvent, bool) or not isinstance(n_solvent, int) or n_solvent < 1:
            raise ValueError(
                f"n_solvent {n_solvent!r} is invalid; it counts solvent "
                "molecules and must be an integer >= 1")

        box_side = kwargs.get("box_side", 20.0)
        if isinstance(box_side, bool) or not isinstance(box_side, (int, float)) or box_side <= 0:
            raise ValueError(
                f"box_side {box_side!r} is invalid; the cubic box side "
                "must be a positive number of Angstrom")

        tolerance = kwargs.get("tolerance", 2.0)
        if isinstance(tolerance, bool) or not isinstance(tolerance, (int, float)) or tolerance <= 0:
            raise ValueError(
                f"tolerance {tolerance!r} is invalid; the packing "
                "tolerance is a minimum inter-molecular distance "
                "(Packmol 2009, DOI 10.1002/jcc.21224) and must be a "
                "positive number of Angstrom")

    def _solute_extent(self, geometry: Path) -> Optional[float]:
        """Largest extent (A) of the solute from its XYZ file, or None
        if the file is unreadable — the box-fit check then stays silent
        rather than blocking a solvable packing problem."""
        try:
            lines = Path(geometry).read_text(encoding="utf-8").splitlines()
            natoms = int(lines[0].strip())
            coords = []
            for line in lines[2:2 + natoms]:
                parts = line.split()
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
        except (ValueError, IndexError, OSError):
            return None
        if not coords:
            return None
        xs, ys, zs = zip(*coords)
        return max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="geometry (solute) is required")

        # ── Box must hold the solute (SV1) ─────────────────────────────
        # The solute is centred in the box (packmol_wrapper.solvate,
        # `fixed 0. 0. 0.` + `center`), so its largest extent must fit
        # inside box_side minus twice the tolerance. A smaller box is a
        # physically meaningless packing problem: Packmol would churn to
        # its 300 s timeout and return "Packmol solvation failed".
        box_side = float(kwargs.get("box_side", 20.0))
        tolerance = float(kwargs.get("tolerance", 2.0))
        extent = self._solute_extent(geometry)
        if extent is not None and extent + 2.0 * tolerance > box_side:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=(
                    f"box_side {box_side} A is too small for the solute: "
                    f"its largest extent is {extent:.2f} A and the "
                    f"tolerance {tolerance} A must be kept clear on both "
                    "sides (Packmol 2009, DOI 10.1002/jcc.21224). "
                    "Increase box_side to at least "
                    f"{extent + 2.0 * tolerance:.2f} A"))

        from delfin.analysis_tools.packmol_wrapper import solvate

        output_xyz = work_dir / "solvated.xyz"
        result_path = solvate(
            solute_xyz=str(geometry),
            solvent_xyz=str(kwargs["solvent_xyz"]),
            n_solvent=kwargs.get("n_solvent", 50),
            box_side=kwargs.get("box_side", 20.0),
            tolerance=kwargs.get("tolerance", 2.0),
            output_xyz=str(output_xyz),
            working_dir=str(work_dir),
        )

        if result_path is not None and Path(result_path).is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                geometry=Path(result_path),
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error="Packmol solvation failed",
        )


register(PackmolSolvateAdapter())
