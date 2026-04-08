"""SMILES-to-XYZ adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class SmilesToXyzAdapter(StepAdapter):
    name = "smiles_to_xyz"
    description = "Convert SMILES string to 3D XYZ geometry via RDKit"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if "smiles" not in kwargs:
            raise ValueError("'smiles' parameter is required")

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()
        smiles = kwargs["smiles"]
        apply_uff = kwargs.get("apply_uff", True)
        hapto_approx = kwargs.get("hapto_approx", None)

        from delfin.smiles_converter import smiles_to_xyz

        xyz_content, error_msg = smiles_to_xyz(
            smiles,
            apply_uff=apply_uff,
            hapto_approx=hapto_approx,
        )

        if error_msg or xyz_content is None:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=error_msg or "smiles_to_xyz returned None",
            )

        # Write XYZ to work_dir
        out_path = work_dir / "structure.xyz"

        # smiles_to_xyz returns raw coordinates without header —
        # count atoms and add standard XYZ header
        lines = [l for l in xyz_content.strip().splitlines() if l.strip()]
        n_atoms = len(lines)
        full_xyz = f"{n_atoms}\nGenerated from {smiles}\n" + "\n".join(lines) + "\n"
        out_path.write_text(full_xyz)

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_path,
            data={"smiles": smiles, "n_atoms": n_atoms},
        )


register(SmilesToXyzAdapter())
