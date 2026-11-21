"""Per-residue confidence with an explicit source interpretation.

Named mmCIF pLDDT metrics are mapped to a model, chain, residue and insertion
code. B-factors are interpreted as pLDDT only when the caller explicitly opts
in. An unavailable metric establishes neither experimental nor predicted
origin. PAE inputs are finite, nonnegative square matrices with a stated upper
bound; the adapter separately checks their size against the chosen model.
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
    ends the reading for the whole structure to avoid reporting partial
    confidence. This check validates an explicitly requested interpretation;
    it does not establish whether coordinates were predicted or experimental.
    """
    out: List[Tuple[object, float]] = []
    for residue in residues:
        atoms = residue.primary_atoms
        if not atoms:
            continue
        values = [atom.b_factor for atom in atoms]
        if any(not math.isfinite(value) or not 0 <= value <= 100 for value in values):
            raise StructureFormatError("Explicit pLDDT B-factors must be finite values in 0–100")
        value = values[0]
        if any(abs(v - value) > 0.01 for v in values):
            return []
        out.append((residue, value))
    return out


def _plddt_from_mmcif_metric(rows_by_chain, residues, model_id, model_count):
    """``ma_qa_metric_local`` rows keyed by (chain, seq_id) -> value."""
    by_key = {}
    for row in rows_by_chain:
        chain_id, seq_id, value = row[:3]
        metric_model, insertion_code = row[3:] if len(row) == 5 else (None, "")
        if metric_model is None and model_count != 1:
            raise StructureFormatError("pLDDT metrics need a model id for a multi-model structure")
        if metric_model is not None and metric_model != model_id:
            continue
        if not math.isfinite(value) or not 0 <= value <= 100:
            raise StructureFormatError("pLDDT metric values must be finite numbers in 0–100")
        key = (chain_id, seq_id, insertion_code)
        if key in by_key:
            raise StructureFormatError(f"Duplicate pLDDT metric for residue {key}")
        by_key[key] = value
    out = []
    for residue in residues:
        value = by_key.get((residue.chain_id, residue.seq_id, residue.insertion_code))
        if value is None:
            raise StructureFormatError(f"No pLDDT metric for residue {residue.label} in model {model_id}")
        out.append((residue, value))
    return out


def confidence_residues(structure, model_id=None, *, b_factor_is_plddt=False):
    """The residues the structure scores, with their pLDDT, or ``()``.

    The mmCIF ``ma_qa_metric_local`` category is preferred where the parser
    kept it (``structure.qa_metric_local``); the B-factor column is the PDB
    route.  An empty tuple means the file carries no per-residue confidence
    that this reader can vouch for.
    """
    if not isinstance(b_factor_is_plddt, bool):
        raise ValueError("b_factor_is_plddt must be a boolean")
    model = structure.first_model if model_id is None else structure.model(model_id)
    residues = [r for c in model.chains for r in c.polymer_residues]
    metric = getattr(structure, "qa_metric_local", None)
    if metric:
        return tuple(_plddt_from_mmcif_metric(metric, residues, model.model_id, len(structure.models)))
    if b_factor_is_plddt:
        return tuple(_plddt_from_atoms(residues))
    return ()


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
    if isinstance(data, list) and len(data) == 1:
        data = data[0]
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
            if not math.isfinite(entry) or entry < 0:
                raise StructureFormatError(f"{source_name}: PAE entries must be finite and nonnegative")
            values.append(float(entry))
        parsed.append(values)
    maximum = data["max_predicted_aligned_error"]
    if isinstance(maximum, bool) or not isinstance(maximum, (int, float)) \
            or not math.isfinite(maximum) or maximum < 0:
        raise StructureFormatError(
            f"{source_name}: 'max_predicted_aligned_error' is not a number"
        )
    if any(value > maximum + 1e-6 for row in parsed for value in row):
        raise StructureFormatError(f"{source_name}: PAE entry exceeds max_predicted_aligned_error")
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
