"""UV-Vis and IR spectrum parsing adapters."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class UvVisParseAdapter(StepAdapter):
    name = "uv_vis_parse"
    description = "Parse UV-Vis absorption spectrum from ORCA output"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "output_file" not in kwargs:
            raise ValueError("'output_file' parameter is required (path to ORCA .out)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.uv_vis_spectrum import parse_absorption_spectrum

        output_file = Path(kwargs["output_file"])
        transitions = parse_absorption_spectrum(output_file)

        if not transitions:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="No UV-Vis transitions found in output file",
            )

        data = {
            "n_transitions": len(transitions),
            "transitions": [
                {
                    "from_state": t.from_state,
                    "to_state": t.to_state,
                    "energy_eV": t.energy_ev,
                    "wavelength_nm": t.wavelength_nm,
                    "fosc": t.fosc,
                    "readable": t.readable_transition,
                }
                for t in transitions
            ],
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=data,
        )


class IrParseAdapter(StepAdapter):
    name = "ir_parse"
    description = "Parse IR vibrational spectrum from ORCA frequency output"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "output_file" not in kwargs:
            raise ValueError("'output_file' parameter is required (path to ORCA freq .out)")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.ir_spectrum import parse_ir_spectrum

        output_file = Path(kwargs["output_file"])
        modes = parse_ir_spectrum(output_file)

        if not modes:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="No IR modes found in output file",
            )

        data = {
            "n_modes": len(modes),
            "modes": [
                {
                    "mode_number": m.mode_number,
                    "frequency_cm1": m.frequency_cm1,
                    "intensity_km_mol": m.intensity_km_mol,
                    "transmittance_percent": m.transmittance_percent,
                }
                for m in modes
            ],
        }

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data=data,
        )


register(UvVisParseAdapter())
register(IrParseAdapter())
