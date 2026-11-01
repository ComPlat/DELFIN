"""SMILES-to-XYZ adapter."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec


class SmilesToXyzAdapter(StepAdapter):
    name = "smiles_to_xyz"
    description = "Convert SMILES string to 3D XYZ geometry via RDKit"
    produces_geometry = True
    category = "structure"
    params = (
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        ParamSpec("apply_uff", "bool", default=True,
                  description="Apply a UFF pre-optimisation"),
        ParamSpec("hapto_approx", "str",
                  description="Hapticity approximation hint (optional)"),
    )
    # Generates a geometry from scratch — consumes no upstream capability.
    data_keys = (
        DataKeySpec("smiles", "str", "", "Echo of the input SMILES"),
        DataKeySpec("n_atoms", "int", "", "Number of atoms in the generated structure"),
        DataKeySpec(
            "findings", "list", "",
            "Post-conversion checker findings (dicts with kind/detail/severity); "
            "empty list when the geometry is clean.  On failure carries the "
            "conversion_failed finding.",
        ),
    )
    requires_python = ("rdkit",)

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

        from delfin.smiles_converter import smiles_to_xyz, check_converted_complex

        xyz_content, error_msg = smiles_to_xyz(
            smiles,
            apply_uff=apply_uff,
            hapto_approx=hapto_approx,
        )

        if error_msg or xyz_content is None:
            # The caller owns the conversion result tuple: a conversion crash
            # is not a bare error string but a structured finding, so pipeline
            # consumers see WHY the structure is absent (QS arbitration,
            # nacht-s11; the caller drops the error today).
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=error_msg or "smiles_to_xyz returned None",
                data={
                    "smiles": smiles,
                    "findings": [{
                        "kind": "conversion_failed",
                        "detail": error_msg or "smiles_to_xyz returned None",
                        "severity": "error",
                    }],
                },
            )

        # Write XYZ to work_dir
        out_path = work_dir / "structure.xyz"

        # smiles_to_xyz returns raw coordinates without header —
        # count atoms and add standard XYZ header
        lines = [l for l in xyz_content.strip().splitlines() if l.strip()]
        n_atoms = len(lines)
        full_xyz = f"{n_atoms}\nGenerated from {smiles}\n" + "\n".join(lines) + "\n"
        out_path.write_text(full_xyz)

        # Post-conversion check (findings only, never a silent fix): the
        # phase-1 checker runs on the written geometry; the caller surfaces
        # its findings in data["findings"].
        try:
            findings = check_converted_complex(smiles, full_xyz)
        except Exception as exc:  # findings must never fail the conversion
            findings = [{
                "kind": "check_error",
                "detail": f"check_converted_complex raised: {exc}",
                "severity": "info",
            }]

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_path,
            data={"smiles": smiles, "n_atoms": n_atoms, "findings": findings},
        )


register(SmilesToXyzAdapter())
