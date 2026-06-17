"""morfeus steric descriptor adapters."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec


class _MorfeusBase(StepAdapter):
    """Shared base for all morfeus adapters."""

    produces_geometry = False
    category = "analysis"
    consumes = ("geometry",)
    requires_python = ("morfeus",)

    def validate_params(self, **kwargs: Any) -> None:
        if "metal_index" not in kwargs:
            raise ValueError("'metal_index' parameter is required")


class MorfeusBuriedVolumeAdapter(_MorfeusBase):
    name = "morfeus_buried_volume"
    description = "Buried volume (%Vbur) via morfeus-ml"
    params = (
        ParamSpec("metal_index", "int", required=True, description="0-based index of the metal centre"),
        ParamSpec("radius", "float", default=3.5, unit="A", description="Sphere radius"),
        ParamSpec("include_hydrogens", "bool", default=True),
        ParamSpec("excluded_atoms", "list", description="Atom indices to exclude"),
    )
    data_keys = (
        DataKeySpec("percent_buried_volume", "float", "%", "%Vbur"),
        DataKeySpec("buried_volume", "float", "A^3"),
        DataKeySpec("free_volume", "float", "A^3"),
        DataKeySpec("total_volume", "float", "A^3"),
    )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.morfeus_wrapper import buried_volume

        result_dict = buried_volume(
            str(geometry),
            metal_index=kwargs["metal_index"],
            radius=kwargs.get("radius", 3.5),
            include_hydrogens=kwargs.get("include_hydrogens", True),
            excluded_atoms=kwargs.get("excluded_atoms"),
        )
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


class MorfeusConeAngleAdapter(_MorfeusBase):
    name = "morfeus_cone_angle"
    description = "Cone angle via morfeus-ml"
    params = (ParamSpec("metal_index", "int", required=True,
                        description="0-based index of the metal centre"),)
    data_keys = (DataKeySpec("cone_angle", "float", "deg", "Tolman cone angle"),)

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.morfeus_wrapper import cone_angle

        result_dict = cone_angle(str(geometry), metal_index=kwargs["metal_index"])
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


class MorfeusSterimolAdapter(StepAdapter):
    name = "morfeus_sterimol"
    description = "Sterimol parameters (L, B1, B5) via morfeus-ml"
    produces_geometry = False
    category = "analysis"
    consumes = ("geometry",)
    requires_python = ("morfeus",)
    params = (
        ParamSpec("atom1_index", "int", required=True, description="Dummy/attachment atom index"),
        ParamSpec("atom2_index", "int", required=True, description="Substituent atom index"),
    )
    data_keys = (
        DataKeySpec("L", "float", "A", "Sterimol L"),
        DataKeySpec("B1", "float", "A", "Sterimol B1"),
        DataKeySpec("B5", "float", "A", "Sterimol B5"),
    )

    def validate_params(self, **kwargs: Any) -> None:
        if "atom1_index" not in kwargs or "atom2_index" not in kwargs:
            raise ValueError("'atom1_index' and 'atom2_index' are required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.morfeus_wrapper import sterimol

        result_dict = sterimol(
            str(geometry),
            atom1_index=kwargs["atom1_index"],
            atom2_index=kwargs["atom2_index"],
        )
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


class MorfeusFullReportAdapter(_MorfeusBase):
    name = "morfeus_full_report"
    description = "Full steric report (Vbur, cone angle, Sterimol) via morfeus-ml"
    params = (
        ParamSpec("metal_index", "int", required=True, description="0-based index of the metal centre"),
        ParamSpec("donor_indices", "list", description="Donor atom indices (enables bite angle)"),
        ParamSpec("radius", "float", default=3.5, unit="A", description="Sphere radius"),
    )
    data_keys = (
        DataKeySpec("buried_volume", "dict", "", "Buried-volume sub-report"),
        DataKeySpec("cone_angle", "dict", "", "Cone-angle sub-report"),
        DataKeySpec("bite_angle", "dict", "", "Bite-angle sub-report (if donor_indices given)"),
    )

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.morfeus_wrapper import full_steric_report

        result_dict = full_steric_report(
            str(geometry),
            metal_index=kwargs["metal_index"],
            donor_indices=kwargs.get("donor_indices"),
            radius=kwargs.get("radius", 3.5),
        )
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


register(MorfeusBuriedVolumeAdapter())
register(MorfeusConeAngleAdapter())
register(MorfeusSterimolAdapter())
register(MorfeusFullReportAdapter())
