"""cclib output parsing adapters."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class CclibParseAdapter(StepAdapter):
    name = "cclib_parse"
    description = "Parse QC output file (ORCA, Gaussian, ...) via cclib"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "filepath" not in kwargs:
            raise ValueError("'filepath' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.cclib_wrapper import parse_output

        filepath = kwargs["filepath"]
        result_dict = parse_output(filepath)

        if result_dict is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="cclib could not parse file")

        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


class CclibEnergiesAdapter(StepAdapter):
    name = "cclib_energies"
    description = "Extract energies from QC output via cclib"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "filepath" not in kwargs:
            raise ValueError("'filepath' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.cclib_wrapper import extract_energies

        result_dict = extract_energies(kwargs["filepath"])

        if result_dict is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="No energies found")

        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


class CclibVibrationsAdapter(StepAdapter):
    name = "cclib_vibrations"
    description = "Extract vibrational frequencies from QC output via cclib"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "filepath" not in kwargs:
            raise ValueError("'filepath' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        from delfin.analysis_tools.cclib_wrapper import extract_vibrations

        result_dict = extract_vibrations(kwargs["filepath"])

        if result_dict is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="No vibrations found")

        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=result_dict)


register(CclibParseAdapter())
register(CclibEnergiesAdapter())
register(CclibVibrationsAdapter())
