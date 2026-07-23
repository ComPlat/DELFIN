"""Solution-entropy adapters for the DELFIN tools layer."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import StepResult, StepStatus


_SOLVENT_PARAMS = (
    ParamSpec("solvent", "str", default="benzene", description="Solvent name or custom label"),
    ParamSpec("temperature", "float", default=298.15, unit="K", description="Temperature"),
    ParamSpec("concentration_standard", "str", default="1M",
              enum=("1M", "liquid"), description="Condensed-phase standard state"),
    ParamSpec("symmetry_number", "int", description="Manual rotational symmetry number"),
    ParamSpec("density_g_ml", "float", unit="g/mL", description="Custom solvent density"),
    ParamSpec("molar_mass_g_mol", "float", unit="g/mol", description="Custom solvent molar mass"),
    ParamSpec("solvent_vdw_volume_A3", "float", unit="A^3", description="Custom solvent vdW volume"),
    ParamSpec("permittivity", "float", description="Custom solvent relative permittivity"),
    ParamSpec("volume_samples", "int", default=60000,
              description="Deterministic Monte-Carlo samples for vdW volume"),
)

_ENTROPY_KEYS = (
    DataKeySpec("S_solv", "float", "cal/mol/K", "Solvation entropy correction"),
    DataKeySpec("S_soln", "float", "cal/mol/K", "Final solution-phase entropy"),
    DataKeySpec("S_trans", "float", "cal/mol/K", "Solution translational entropy"),
    DataKeySpec("S_rot", "float", "cal/mol/K", "Solution rotational entropy"),
    DataKeySpec("S_vib_qrrho", "float", "cal/mol/K", "qRRHO vibrational entropy"),
    DataKeySpec("S_cav", "float", "cal/mol/K", "Cavitation entropy"),
    DataKeySpec("S_conc", "float", "cal/mol/K", "Standard-state entropy correction"),
    DataKeySpec("vdw_volume_A3", "float", "A^3", "Solute vdW volume"),
    DataKeySpec("radius_gyration_A", "float", "A", "Radius of gyration"),
    DataKeySpec("report_path", "str", "", "JSON report path"),
)


class SolutionEntropyAdapter(StepAdapter):
    name = "solution_entropy"
    description = "Solution-phase entropy correction from geometry and optional QC thermochemistry"
    produces_geometry = False
    category = "solvation"
    consumes = ("geometry", "qc_output")
    wires = {"qc_output": "filepath"}
    params = (
        ParamSpec("filepath", "path", description="ORCA/Gaussian output for thermochemistry"),
    ) + _SOLVENT_PARAMS
    produces = ("solution_entropy_report",)
    data_keys = _ENTROPY_KEYS

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="geometry is required",
            )

        from delfin.solvation_entropy import (
            calculate_solution_entropy,
            parse_thermochemistry,
            write_report,
        )

        thermo = None
        filepath = kwargs.get("filepath")
        if filepath:
            try:
                thermo = parse_thermochemistry(filepath)
            except Exception as exc:
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    error=f"thermochemistry parse failed: {exc}",
                )

        try:
            components = calculate_solution_entropy(
                geometry,
                solvent=kwargs.get("solvent", "benzene"),
                temperature_K=float(kwargs.get("temperature", 298.15)),
                thermochemistry=thermo,
                concentration_standard=kwargs.get("concentration_standard", "1M"),
                symmetry_number=kwargs.get("symmetry_number"),
                density_g_ml=kwargs.get("density_g_ml"),
                molar_mass_g_mol=kwargs.get("molar_mass_g_mol"),
                solvent_vdw_volume_A3=kwargs.get("solvent_vdw_volume_A3"),
                permittivity=kwargs.get("permittivity"),
                volume_samples=int(kwargs.get("volume_samples", 60000)),
            )
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"solution entropy calculation failed: {exc}",
            )

        report = write_report(work_dir / "solution_entropy.json", components)
        data = components.to_dict()
        data["report_path"] = str(report)
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            data=data,
            artifacts={"solution_entropy_report": report},
        )


class ReactionSolutionEntropyAdapter(StepAdapter):
    name = "reaction_solution_entropy"
    description = "Stoichiometric reaction or barrier correction from solution-entropy reports"
    produces_geometry = False
    category = "solvation"
    params = (
        ParamSpec("species", "dict", required=True,
                  description="Map species name to report path, S_solv value, or report dict"),
        ParamSpec("stoichiometry", "dict", required=True,
                  description="Map species name to coefficient; products positive"),
        ParamSpec("temperature", "float", default=298.15, unit="K"),
        ParamSpec("uncorrected_delta_g_kcal_mol", "float", unit="kcal/mol"),
    )
    produces = ("reaction_entropy_report",)
    data_keys = (
        DataKeySpec("delta_S_solv_cal_mol_K", "float", "cal/mol/K"),
        DataKeySpec("delta_G_entropy_corr_kcal_mol", "float", "kcal/mol"),
        DataKeySpec("delta_G_corrected_kcal_mol", "float", "kcal/mol"),
    )

    def validate_params(self, **kwargs: Any) -> None:
        if "species" not in kwargs:
            raise ValueError("'species' parameter is required")
        if "stoichiometry" not in kwargs:
            raise ValueError("'stoichiometry' parameter is required")

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()
        import json
        from delfin.solvation_entropy import calculate_reaction_entropy_correction, read_report

        species_arg = _coerce_mapping(kwargs["species"], "species")
        stoich_arg = _coerce_mapping(kwargs["stoichiometry"], "stoichiometry")
        species_in = dict(species_arg)
        species: dict[str, Any] = {}
        for name, value in species_in.items():
            if isinstance(value, (str, Path)) and Path(value).is_file():
                species[name] = read_report(value)
            else:
                species[name] = value
        try:
            correction = calculate_reaction_entropy_correction(
                species,
                dict(stoich_arg),
                temperature_K=float(kwargs.get("temperature", 298.15)),
                uncorrected_delta_g_kcal_mol=kwargs.get("uncorrected_delta_g_kcal_mol"),
            )
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"reaction entropy correction failed: {exc}",
            )
        data = correction.to_dict()
        report = work_dir / "reaction_solution_entropy.json"
        report.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            data=data,
            artifacts={"reaction_entropy_report": report},
        )


def _coerce_mapping(value: Any, name: str) -> dict[str, Any]:
    """Accept dicts or JSON object strings from CLIs/keyfiles."""

    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        import json

        try:
            parsed = json.loads(value)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{name} must be a JSON object string or dict") from exc
        if isinstance(parsed, dict):
            return parsed
    raise ValueError(f"{name} must be a dict")


register(SolutionEntropyAdapter())
register(ReactionSolutionEntropyAdapter())
