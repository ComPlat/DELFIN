"""Pose RMSD, as a building block: how far each pose sits from a reference.

Input: a pose set (a multi-record SDF, one placement per record) and a
reference ligand file — the co-crystallised ligand for a redocking check,
or any chosen pose.  Output: the heavy-atom RMSD of every pose against the
reference, in Angstrom, without re-aligning — the docking frame is the
frame — plus how many poses are within 2 A.

Semantics: the pairing of the atoms is symmetry-aware (equal atoms of the
molecule's own symmetry — ring atoms, methyl hydrogens' heavy partners —
are matched to whichever pairing gives the smallest RMSD), so a pose
written in another atom order is not far away because of it.  Hydrogens
are excluded: the pose quality asked about is the heavy-atom geometry.
A pose set holding a different molecule than the reference is refused:
an RMSD across molecules is a number nobody asked for.  Poses the pairing
cannot use are reported with their error, not dropped — a skipped pose in
a table reads as a good pose.

RDKit is required and declared (``requires_python``); it is
BSD-3-licensed.  No external docking program is needed or used: this
block measures what a docking run left on disk.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, List, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.bio import (
    StructureFormatError,
    read_pose_set,
    read_small_molecule,
    pose_rmsd,
)

_REPORT_NAME = "pose_rmsd.json"
UNDER_2A = 2.0


def _output_file_name(value: Any) -> str:
    """The report's file name, or raise ``ValueError`` (see the ligand
    adapter: a separator would let a parameter name any file)."""
    name = str(value or "").strip() or _REPORT_NAME
    if (name in (".", "..") or "/" in name or "\\" in name
            or Path(name).name != name):
        raise ValueError(f"'report_name' must be a plain file name: {name!r}")
    return name


class PoseRmsdAdapter(StepAdapter):
    """RMSD of every pose of a pose set against one reference ligand."""

    name = "pose_rmsd"
    description = (
        "Symmetry-aware heavy-atom RMSD of every pose in a pose set "
        "against a reference ligand (redocking check), without "
        "re-alignment; reports how many poses are within 2 A"
    )
    category = "structure"
    produces_geometry = False
    consumes = ("pose_set", "ligand")
    produces = ("json_report",)
    params = (
        ParamSpec("pose_file", "path", required=True,
                  description="Multi-record SDF, one record per pose"),
        ParamSpec("reference_file", "path", required=True,
                  description="SDF/MOL2 holding the reference ligand"),
        ParamSpec("score_tag", "str",
                  description="SD tag of the pose score, when the file "
                              "carries one (repeated into the report)"),
        ParamSpec("cutoff", "float", default=UNDER_2A,
                  description="RMSD under which a pose counts as redocked "
                              "(Angstrom; default 2.0, the common "
                              "redocking criterion)"),
        ParamSpec("report_name", "str", default=_REPORT_NAME,
                  description="File name of the JSON report in work_dir"),
    )
    data_keys = (
        DataKeySpec("n_poses", "int"),
        DataKeySpec("n_reported", "int",
                    description="Poses whose RMSD was computed"),
        DataKeySpec("n_failed", "int",
                    description="Poses the pairing refused, with their "
                                "errors kept in the report"),
        DataKeySpec("n_under_cutoff", "int",
                    description="Poses with RMSD under the cutoff"),
        DataKeySpec("cutoff", "float", unit="Angstrom"),
        DataKeySpec("best_rmsd", "float", unit="Angstrom",
                    description="The smallest RMSD among the reported "
                                "poses, or None when none reported"),
        DataKeySpec("alignment", "str",
                    description="Always 'none': the docking frame is "
                                "the frame"),
    )
    wires = (("pose_set", "pose_file"), ("ligand", "reference_file"),)

    def validate_params(self, **kwargs: Any) -> None:
        for key in ("pose_file", "reference_file"):
            if not kwargs.get(key):
                raise ValueError(f"'{key}' parameter is required")
        _output_file_name(kwargs.get("report_name"))
        cutoff = kwargs.get("cutoff")
        if cutoff is not None:
            try:
                value = float(cutoff)
            except (TypeError, ValueError):
                raise ValueError(f"'cutoff' must be a number: {cutoff!r}")
            if value <= 0:
                raise ValueError(f"'cutoff' must be positive: {cutoff!r}")

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
            report_name = _output_file_name(kwargs.get("report_name"))
            cutoff = (float(kwargs["cutoff"])
                      if kwargs.get("cutoff") is not None else UNDER_2A)
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )

        paths = {}
        for key in ("pose_file", "reference_file"):
            path = Path(kwargs.get(key))
            if not path.is_file():
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    error=f"{key} not found: {path.name}",
                    error_kind=ErrorKind.MISSING_INPUT,
                )
            paths[key] = path

        try:
            pose_set = read_pose_set(
                paths["pose_file"],
                score_tag=(str(kwargs["score_tag"])
                           if kwargs.get("score_tag") else None),
            )
            reference_records = read_small_molecule(paths["reference_file"])
            if len(reference_records) != 1:
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    error=f"the reference file {paths['reference_file'].name} "
                          f"holds {len(reference_records)} records; one "
                          "reference ligand is required",
                    error_kind=ErrorKind.UNSUITABLE_INPUT,
                )
            reference = reference_records[0]
        except StructureFormatError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.PARSE,
            )

        rows: List[dict] = []
        failures: List[dict] = []
        for pose in pose_set.poses:
            row = {"index": pose.index, "name": pose.name,
                   "score": pose.score}
            try:
                row["rmsd_A"] = pose_rmsd(pose.molecule, reference)
                rows.append(row)
            except StructureFormatError as exc:
                row["error"] = str(exc)
                failures.append(row)

        under = sum(1 for r in rows if r["rmsd_A"] < cutoff)
        best = min((r["rmsd_A"] for r in rows), default=None)
        report = {
            "pose_file": paths["pose_file"].name,
            "reference_file": paths["reference_file"].name,
            "alignment": "none",
            "cutoff_A": cutoff,
            "n_poses": len(pose_set.poses),
            "n_reported": len(rows),
            "n_failed": len(failures),
            "n_under_cutoff": under,
            "best_rmsd_A": best,
            "poses": rows,
            "failed_poses": failures,
        }
        report_path = Path(work_dir) / report_name
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            output_file=report_path,
            data={
                "n_poses": len(pose_set.poses),
                "n_reported": len(rows),
                "n_failed": len(failures),
                "n_under_cutoff": under,
                "cutoff": cutoff,
                "best_rmsd": best,
                "alignment": "none",
            },
            artifacts={"json_report": report_path},
        )


register(PoseRmsdAdapter())
