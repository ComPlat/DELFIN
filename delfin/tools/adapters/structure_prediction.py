"""The structure-prediction engine, as a building block.

Input: a FASTA file (``sequence_file``, the ``sequence`` port).  Output: a
:class:`~delfin.tools._types.StepResult` whose ``data`` holds the sequence
the file named and one entry per prediction engine — detected or not, and
the command line that would run it on this sequence.

Semantics: the sequence is read through :mod:`delfin.tools.bio` (a DNA or
RNA record is refused by name — a prediction cannot be run on it), and the
engines are only *detected*, through ``shutil.which`` under the ``PATH`` the
caller hands in (``_prediction_env``, default the process environment).  No
engine is executed, installed or downloaded by this block; the command is
described so a caller who has read the engine's documentation can run it
where they choose.  Detection without running is pinned by a test whose
stand-in would fail the test if it were executed.

Licences of the programs this block names (code licences; weights carry
their own terms): ColabFold MIT, ESMFold MIT, Boltz MIT, AlphaFold 3
CC BY-NC-SA 4.0 (non-commercial code; model terms must be accepted).

Failures are kinds, not prose: a missing sequence file is ``MISSING_INPUT``,
an unreadable FASTA is ``PARSE``.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.bio import parse_fasta, prediction_engines
from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._reports import report_file_name


class StructurePredictionEngineAdapter(StepAdapter):
    """Detect the prediction engines and describe the command for a sequence."""

    name = "structure_prediction_engine"
    description = (
        "Detect protein-structure prediction engines (ColabFold, ESMFold, "
        "Boltz, AlphaFold 3) on PATH and describe the command each would run "
        "on a sequence file.  Detection only: no engine is run, installed "
        "or downloaded"
    )
    category = "structure"
    produces_geometry = False
    consumes = ("sequence",)
    produces = ("json_report",)
    params = (
        ParamSpec("sequence_file", "path", required=True,
                  description="FASTA file with the protein sequence"),
        ParamSpec("report_name", "str", default="prediction_engines.json",
                  description="File name of the JSON report written into "
                              "work_dir"),
    )
    data_keys = (
        DataKeySpec("sequence_id", "str", description="First record's id"),
        DataKeySpec("sequence_length", "int", description="First record's length"),
        DataKeySpec("n_engines", "int", description="Engines the table holds"),
        DataKeySpec("engines", "dict",
                    description="Per engine: detected, path, licence, "
                                "described command"),
        DataKeySpec("detected_engines", "list",
                    description="Names of the engines found on PATH"),
    )
    # A later block that hands on a sequence file fills this step's input.
    wires = (("sequence", "sequence_file"),)

    def validate_params(self, **kwargs: Any) -> None:
        if not kwargs.get("sequence_file"):
            raise ValueError(
                "'sequence_file' parameter is required (a FASTA file)"
            )
        report_file_name(kwargs.get("report_name"), default="prediction_engines.json")

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            report_name = report_file_name(
                kwargs.get("report_name"), default="prediction_engines.json",
            )
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )

        raw = kwargs.get("sequence_file")
        if not raw:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="'sequence_file' is required (a FASTA file)",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        path = Path(raw)
        if not path.is_file():
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"sequence file not found: {path.name}",
                error_kind=ErrorKind.MISSING_INPUT,
            )
        try:
            records = parse_fasta(
                path.read_text(encoding="utf-8"), source_name=path.name)
        except StructureFormatError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.PARSE,
            )
        first = records[0]

        env = kwargs.get("_prediction_env")
        entries = prediction_engines(env=env, fasta=str(path))
        detected = [e["name"] for e in entries if e["detected"]]

        data: Dict[str, Any] = {
            "sequence_id": first.id,
            "sequence_length": first.length,
            "n_engines": len(entries),
            "engines": entries,
            "detected_engines": detected,
            "source_name": path.name,
        }

        report_path = Path(work_dir) / report_name
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(data, indent=2, sort_keys=False),
                               encoding="utf-8")

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            output_file=report_path,
            data=data,
            artifacts={"json_report": report_path},
        )


from delfin.tools._registry import register  # noqa: E402

register(StructurePredictionEngineAdapter())
