"""Reporting adapters (JSON data collection, DOCX report generation)."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class ReportJsonAdapter(StepAdapter):
    name = "report_json"
    description = "Collect DELFIN project data into JSON"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "project_dir" not in kwargs:
            raise ValueError("'project_dir' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.reporting.delfin_collector import save_esd_data_json

        project_dir = Path(kwargs["project_dir"])
        output_json = work_dir / "DELFIN_Data.json"

        try:
            json_path = save_esd_data_json(project_dir, output_json=output_json)
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"JSON collection failed: {exc}",
            )

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data={"json_path": str(json_path)},
            artifacts={"json_report": json_path},
        )


class ReportDocxAdapter(StepAdapter):
    name = "report_docx"
    description = "Generate DELFIN Word report (DOCX)"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "project_dir" not in kwargs:
            raise ValueError("'project_dir' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None, cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()

        from delfin.reporting.delfin_collector import save_esd_data_json
        from delfin.reporting.delfin_docx_report import generate_combined_docx_report

        project_dir = Path(kwargs["project_dir"])

        # First collect JSON data
        json_path = work_dir / "DELFIN_Data.json"
        if not json_path.exists():
            try:
                json_path = save_esd_data_json(project_dir, output_json=json_path)
            except Exception as exc:
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    error=f"JSON collection failed: {exc}",
                )

        # Generate DOCX
        output_docx = work_dir / "DELFIN.docx"
        try:
            docx_path = generate_combined_docx_report(
                project_dir=project_dir,
                json_path=json_path,
                output_docx=output_docx,
            )
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"DOCX generation failed: {exc}",
            )

        if docx_path is not None and Path(docx_path).is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                data={"docx_path": str(docx_path)},
                artifacts={"docx_report": Path(docx_path)},
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            error="DOCX generation produced no output",
        )


register(ReportJsonAdapter())
register(ReportDocxAdapter())
