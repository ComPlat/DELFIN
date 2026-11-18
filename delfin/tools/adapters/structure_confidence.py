"""Confidence of a predicted model, as a building block.

Input: a PDB or mmCIF file (``structure_file``, the ``macromolecule`` port)
and optionally the confidence JSON the AlphaFold Database publishes beside a
model (``pae_file``).  Output: a :class:`~delfin.tools._types.StepResult`
whose ``data`` holds what the file says about its own certainty and whose
``json_report`` artifact holds the same as a file.

Reported, and why each piece:

* ``mean_plddt`` and ``mean_plddt_per_chain`` — the per-residue pLDDT the
  file carries, averaged overall and per chain.  A chain average is what a
  builder looks at first: one floppy loop does not make a model unusable, a
  whole chain at 40 does.
* ``ranges_below_70`` and ``ranges_below_50`` — the runs of consecutive
  residues under each band, as chain/start/end/length.  A range is
  actionable (delete it, model a ligand elsewhere); 200 single-residue flags
  are not.  The bands are AlphaFold's own: 70 confident, 50 low, below that
  very low.
* ``has_confidence`` / ``plddt_source`` / ``statement`` — an experimental
  structure carries no confidence, and this block says so plainly, in a
  statement a report can quote, instead of zeros that read as "the model is
  0% sure".
* ``has_pae`` / ``pae_max`` / ``pae_n_residues`` — the predicted-aligned-
  error matrix when a ``pae_file`` is given, read with its documented keys
  (``predicted_aligned_error``, ``max_predicted_aligned_error``; Váradi et
  al., Nucleic Acids Research 52:D1, 2024, doi:10.1093/nar/gkad1011).

Semantics: the pLDDT is read from the file itself — the mmCIF
``ma_qa_metric_local`` category where the file names the metric, else the
B-factor column where a prediction wrote the same per-residue value on every
atom.  B-factors that differ between the atoms of a residue, or sit outside
0-100, or are identical across the whole file are *no confidence*, never
numbers: a crystallographic B-factor read as a pLDDT would be a plausible
wrong answer about how sure a model is.

The block runs no external program (``requires_binaries`` and
``requires_python`` are empty) and reads through :mod:`delfin.tools.bio`.

Failures are kinds, not prose: a missing file is ``MISSING_INPUT``, an
unreadable structure is ``PARSE``, an unreadable PAE file is ``PARSE``.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from delfin.tools._base import StepAdapter
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.bio import PLDDT_LOW, PLDDT_VERY_LOW, read_pae, read_structure
from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._reports import report_file_name


def _ranges(values, threshold: float) -> List[Dict[str, Any]]:
    """Consecutive runs of residues whose pLDDT is below *threshold*.

    Consecutive means consecutive sequence ids within one chain; a chain
    break, an insertion code or a gap in the numbering ends a range.
    """
    ranges: List[Dict[str, Any]] = []
    current: Optional[Dict[str, Any]] = None
    for residue, value in values:
        if value >= threshold:
            current = None
            continue
        same_run = (
            current is not None
            and current["chain_id"] == residue.chain_id
            and residue.seq_id == current["end_seq_id"] + 1
            and not residue.insertion_code
        )
        if same_run:
            current["end_seq_id"] = residue.seq_id
            current["n_residues"] += 1
        else:
            current = {"chain_id": residue.chain_id,
                       "start_seq_id": residue.seq_id,
                       "end_seq_id": residue.seq_id,
                       "n_residues": 1}
            ranges.append(current)
    return ranges


def _statement(values, source: str) -> str:
    """A sentence a report can quote about what the file says."""
    if not values:
        return ("The file carries no per-residue confidence: this reads as an "
                "experimental structure, not a predicted model.")
    scored = [v for _, v in values]
    mean = sum(scored) / len(scored)
    below50 = sum(1 for v in scored if v < PLDDT_VERY_LOW)
    below70 = sum(1 for v in scored if v < PLDDT_LOW)
    text = (f"Per-residue pLDDT read from the file ({source}): mean "
            f"{mean:.1f} over {len(scored)} residues, "
            f"{below70} below {PLDDT_LOW:.0f} and {below50} below "
            f"{PLDDT_VERY_LOW:.0f}.")
    return text


class StructureConfidenceAdapter(StepAdapter):
    """Read what a structure file says about its own certainty."""

    name = "structure_confidence"
    description = (
        "Confidence of a predicted model, read from the file itself: per-"
        "residue pLDDT (mmCIF ma_qa_metric_local or PDB B-factor), mean per "
        "chain, residue ranges below 70 and 50, and the AlphaFold PAE JSON "
        "when given.  An experimental structure is reported as carrying no "
        "confidence"
    )
    category = "structure"
    produces_geometry = False
    consumes = ("macromolecule",)
    produces = ("json_report",)
    params = (
        ParamSpec("structure_file", "path", required=True,
                  description="Macromolecular structure file to read "
                              "(PDB or mmCIF)"),
        ParamSpec("pae_file", "path",
                  description="Optional PAE JSON (predicted_aligned_error / "
                              "max_predicted_aligned_error) beside the model"),
        ParamSpec("report_name", "str", default="structure_confidence.json",
                  description="File name of the JSON report written into "
                              "work_dir"),
    )
    data_keys = (
        DataKeySpec("has_confidence", "bool",
                    description="Whether the file carries per-residue "
                                "confidence"),
        DataKeySpec("plddt_source", "str",
                    description="mmcif_ma_qa_metric_local, pdb_b_factor or "
                                "none"),
        DataKeySpec("mean_plddt", "float", "pLDDT",
                    description="Mean pLDDT over the scored residues"),
        DataKeySpec("mean_plddt_per_chain", "dict", "pLDDT",
                    description="Chain id -> mean pLDDT"),
        DataKeySpec("n_residues_scored", "int"),
        DataKeySpec("ranges_below_70", "list",
                    description="Consecutive residue runs below 70: chain, "
                                "start, end, length"),
        DataKeySpec("ranges_below_50", "list",
                    description="Consecutive residue runs below 50: chain, "
                                "start, end, length"),
        DataKeySpec("has_pae", "bool"),
        DataKeySpec("pae_max", "float", "A",
                    description="max_predicted_aligned_error"),
        DataKeySpec("pae_n_residues", "int"),
        DataKeySpec("statement", "str",
                    description="A sentence a report can quote"),
    )
    # A later block that hands on a predicted model fills this step's input.
    wires = (("macromolecule", "structure_file"),)

    def validate_params(self, **kwargs: Any) -> None:
        if not kwargs.get("structure_file"):
            raise ValueError(
                "'structure_file' parameter is required (a PDB or mmCIF file)"
            )
        report_file_name(kwargs.get("report_name"), default="structure_confidence.json")

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
                kwargs.get("report_name"), default="structure_confidence.json",
            )
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )

        raw = kwargs.get("structure_file")
        if not raw:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="'structure_file' is required (a PDB or mmCIF file)",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        path = Path(raw)
        if not path.is_file():
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"structure file not found: {path.name}",
                error_kind=ErrorKind.MISSING_INPUT,
            )

        try:
            structure = read_structure(path)
        except StructureFormatError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.PARSE,
            )

        values = structure.plddt_per_residue()
        if values:
            source = ("mmcif_ma_qa_metric_local"
                      if structure.source_format == "mmcif"
                      and structure.qa_metric_local else "pdb_b_factor")
        else:
            source = "none"

        scored = [v for _, v in values]
        per_chain: Dict[str, float] = {}
        chain_values: Dict[str, List[float]] = {}
        for residue, value in values:
            chain_values.setdefault(residue.chain_id, []).append(value)
        for chain_id, vals in chain_values.items():
            per_chain[chain_id] = sum(vals) / len(vals)
        mean = sum(scored) / len(scored) if scored else None

        data: Dict[str, Any] = {
            "has_confidence": bool(values),
            "plddt_source": source,
            "mean_plddt": mean,
            "mean_plddt_per_chain": per_chain,
            "n_residues_scored": len(scored),
            "ranges_below_70": _ranges(values, PLDDT_LOW),
            "ranges_below_50": _ranges(values, PLDDT_VERY_LOW),
            "has_pae": False,
            "pae_max": None,
            "pae_n_residues": None,
            "source_name": structure.source_name,
            "source_format": structure.source_format,
        }

        pae_raw = kwargs.get("pae_file")
        if pae_raw:
            pae_path = Path(pae_raw)
            if not pae_path.is_file():
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    error=f"PAE file not found: {pae_path.name}",
                    error_kind=ErrorKind.MISSING_INPUT,
                )
            try:
                pae = read_pae(pae_path)
            except StructureFormatError as exc:
                return self._make_result(
                    self.name, StepStatus.FAILED, work_dir, start,
                    error=str(exc), error_kind=ErrorKind.PARSE,
                )
            data["has_pae"] = True
            data["pae_max"] = pae.max_predicted_aligned_error
            data["pae_n_residues"] = pae.n_residues

        data["statement"] = _statement(values, source)

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

register(StructureConfidenceAdapter())
