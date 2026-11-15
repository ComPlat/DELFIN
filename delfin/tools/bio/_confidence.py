"""Confidence of a predicted model, read from the file itself.

Input: the parsed structure (``plddt_per_residue`` on
:class:`~delfin.tools.bio.MacromolecularStructure`) or the JSON file the
AlphaFold Database publishes beside a model (``read_pae``).  Output: the
per-residue pLDDT values the file carries, and the predicted-aligned-error
matrix as :class:`PredictedAlignedError`.

Semantics, and the decision that keeps the number honest:

* **pLDDT is a per-residue value; the B-factor column is only a prediction's
  way of storing it.**  An AlphaFold-style PDB writes the same B-factor for
  every atom of one residue, in 0-100.  An experimental structure does not:
  its B-factors differ between the atoms of a residue, or are zero across an
  NMR ensemble, or sit outside the 0-100 range.  All of those are read as
  *no confidence*, never as numbers — a crystallographic B-factor read as a
  pLDDT would be a plausible wrong answer about how sure a model is.
* **mmCIF carries the value where the format says confidence lives.**  The
  ``ma_qa_metric_local`` category, rows of one ``pLDDT`` metric, is read
  directly; it outranks the B-factor column when both are present, because
  the file names the metric rather than borrowing a column.
* **The PAE JSON is read with its documented keys** (``predicted_aligned_error``,
  ``max_predicted_aligned_error``; Váradi et al., Nucleic Acids Research
  52:D1, 2024, doi:10.1093/nar/gkad1011).  A matrix that is not square, a
  file without the keys, or values that are not numbers are refused — half
  of a confidence statement is worse than none.

Rejected alternative: guessing a prediction from a REMARK or from the file
name.  A depositor's label is not the file's own arithmetic, and the labels
differ between producers.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import List, Tuple, Union

from delfin.tools.bio._errors import StructureFormatError

_PATHLIKE = Union[str, Path]

#: pLDDT below this is "low confidence" (AlphaFold's own band).
PLDDT_LOW = 70.0
#: pLDDT below this is "very low confidence" (AlphaFold's own band).
PLDDT_VERY_LOW = 50.0

_DOCUMENTED_KEYS = ("predicted_aligned_error", "max_predicted_aligned_error")


# ---------------------------------------------------------------------------
# pLDDT from the structure
# ---------------------------------------------------------------------------


def _plddt_from_atoms(residues) -> List[Tuple[object, float]]:
    """Per-residue pLDDT values, or nothing when the file carries none.

    A residue yields a value when every one of its primary-atom B-factors is
    the same number within 0.01 (a prediction writes the per-residue value on
    each atom) and inside 0-100.  The first residue that breaks the pattern
    ends the reading for the whole structure: one residue with
    crystallographic B-factors means the file's B-factors are not per-residue
    confidence, and a half-read would report the unmodified part of a
    structure as if it were a prediction.
    """
    out: List[Tuple[object, float]] = []
    for residue in residues:
        atoms = residue.primary_atoms
        if not atoms:
            continue
        values = [atom.b_factor for atom in atoms]
        if any(v < -0.01 or v > 100.01 for v in values):
            return []
        value = values[0]
        if any(abs(v - value) > 0.01 for v in values):
            return []
        if value <= 0.0 and len(out) == 0 and all(v == 0.0 for v in values):
            # An all-zero first residue is an NMR structure until another
            # residue proves otherwise; a prediction writes real values.
            continue
        out.append((residue, value))
    return out


def _plddt_from_mmcif_metric(rows_by_chain, residues):
    """``ma_qa_metric_local`` rows keyed by (chain, seq_id) -> value."""
    by_key = {}
    for chain_id, seq_id, value in rows_by_chain:
        by_key[(chain_id, seq_id)] = value
    out = []
    for residue in residues:
        value = by_key.get((residue.chain_id, residue.seq_id))
        if value is None:
            return []          # an un-scored residue: no half-read
        out.append((residue, value))
    return out


def confidence_residues(structure, model_id=None):
    """The residues the structure scores, with their pLDDT, or ``()``.

    The mmCIF ``ma_qa_metric_local`` category is preferred where the parser
    kept it (``structure.qa_metric_local``); the B-factor column is the PDB
    route.  An empty tuple means the file carries no per-residue confidence
    that this reader can vouch for.
    """
    model = structure.first_model if model_id is None else structure.model(model_id)
    residues = [r for c in model.chains for r in c.polymer_residues]
    metric = getattr(structure, "qa_metric_local", None)
    if metric:
        values = _plddt_from_mmcif_metric(metric, residues)
        if values:
            return tuple(values)
    return tuple(_plddt_from_atoms(residues))


# ---------------------------------------------------------------------------
# the PAE JSON
# ---------------------------------------------------------------------------


class PredictedAlignedError:
    """The PAE matrix and its maximum, in Angstrom, as the file stated them."""

    __slots__ = ("matrix", "max_predicted_aligned_error", "source_name")

    def __init__(self, matrix: List[List[float]],
                 max_predicted_aligned_error: float,
                 source_name: str = "") -> None:
        self.matrix = matrix
        self.max_predicted_aligned_error = max_predicted_aligned_error
        self.source_name = source_name

    @property
    def n_residues(self) -> int:
        return len(self.matrix)

    def mean_off_diagonal(self) -> float:
        """Mean PAE over all off-diagonal entries (the diagonal is 0 by definition)."""
        n = self.n_residues
        total = 0.0
        count = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    total += self.matrix[i][j]
                    count += 1
        return total / count if count else 0.0


def parse_pae_json(data, source_name: str = "<json>") -> PredictedAlignedError:
    """Read the PAE object with its documented keys; refuse anything else."""
    if not isinstance(data, dict):
        raise StructureFormatError(
            f"{source_name}: the confidence JSON is not an object with the "
            "documented AlphaFold keys"
        )
    for key in _DOCUMENTED_KEYS:
        if key not in data:
            raise StructureFormatError(
                f"{source_name}: missing '{key}' — the documented keys of an "
                "AlphaFold confidence JSON are predicted_aligned_error and "
                "max_predicted_aligned_error"
            )
    matrix = data["predicted_aligned_error"]
    if (not isinstance(matrix, list) or not matrix
            or not all(isinstance(row, list) for row in matrix)):
        raise StructureFormatError(
            f"{source_name}: 'predicted_aligned_error' is not a matrix of numbers"
        )
    n = len(matrix)
    parsed: List[List[float]] = []
    for row in matrix:
        if len(row) != n:
            raise StructureFormatError(
                f"{source_name}: the PAE matrix is not square "
                f"({n} rows, a row of {len(row)}) — a predicted-aligned-error "
                "matrix is one error per residue pair"
            )
        values = []
        for entry in row:
            if isinstance(entry, bool) or not isinstance(entry, (int, float)):
                raise StructureFormatError(
                    f"{source_name}: the PAE matrix holds a value that is not "
                    f"a number ({entry!r})"
                )
            values.append(float(entry))
        parsed.append(values)
    maximum = data["max_predicted_aligned_error"]
    if isinstance(maximum, bool) or not isinstance(maximum, (int, float)) \
            or not math.isfinite(maximum):
        raise StructureFormatError(
            f"{source_name}: 'max_predicted_aligned_error' is not a number"
        )
    return PredictedAlignedError(parsed, float(maximum), source_name)


def read_pae(path: _PATHLIKE) -> PredictedAlignedError:
    """Read the PAE JSON the AlphaFold Database publishes beside a model."""
    p = Path(path)
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except OSError as exc:
        raise StructureFormatError(f"{p.name}: cannot be read ({exc})") from None
    except json.JSONDecodeError as exc:
        raise StructureFormatError(
            f"{p.name}: not a JSON file ({exc})"
        ) from None
    return parse_pae_json(data, source_name=p.name)


__all__ = [
    "PLDDT_LOW",
    "PLDDT_VERY_LOW",
    "PredictedAlignedError",
    "confidence_residues",
    "parse_pae_json",
    "read_pae",
]
