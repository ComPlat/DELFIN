"""Ligand preparation, as a building block.

Input: a SMILES string or one small-molecule file (SDF/MOL or MOL2).
Output: an SDF file in ``work_dir`` holding the prepared ligand — 3D
coordinates, explicit hydrogens, Gasteiger charges — emitted as the
``ligand`` capability port, so a later block (docking, a receptor
connection, a QM calculation) consumes it without knowing how it was made.

Semantics: the embedding method and seed are named in the result data
(``ETDGv3``, fixed seed) so the geometry is reproducible; the file
ligand's coordinates are kept as they are and only completed with
hydrogens.  The protonation state is taken as given — no pKa model is
applied — and the result data says so.  An input holding more than one
fragment is refused unless ``keep_largest`` asks for the largest one, in
which case the removed fragments are named in the data.  A file with more
than one record is refused: which record is the ligand is a question for
the caller, not a default.

RDKit is required and declared in the contract (``requires_python``), so
the environment probe reports this block honestly on a machine without
it.  RDKit is BSD-3-licensed.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.bio import StructureFormatError, prepare_ligand, write_sdf

_OUTPUT_NAME = "prepared_ligand.sdf"


def _output_file_name(value: Any) -> str:
    """The output's file name, or raise ``ValueError``.

    The ligand is written into ``work_dir`` and nowhere else; a separator
    or a ``..`` in the name would let a parameter choose any file to
    overwrite, so they are refused rather than stripped.
    """
    name = str(value or "").strip() or _OUTPUT_NAME
    if (name in (".", "..") or "/" in name or "\\" in name
            or Path(name).name != name):
        raise ValueError(f"'output_name' must be a plain file name: {name!r}")
    return name


class LigandPrepareAdapter(StepAdapter):
    """Prepare one 3D ligand with explicit hydrogens and Gasteiger charges."""

    name = "ligand_prepare"
    description = (
        "Prepare a ligand from SMILES or an SDF/MOL2 file: 3D coordinates "
        "(ETDGv3, fixed seed; file coordinates kept), explicit hydrogens, "
        "Gasteiger charges as an SD property; the protonation state is "
        "taken as given (no pKa model)"
    )
    category = "structure"
    produces_geometry = False
    consumes = ()
    produces = ("ligand",)
    params = (
        ParamSpec("smiles", "str",
                  description="SMILES of the ligand (exactly one of "
                              "smiles / ligand_file)"),
        ParamSpec("ligand_file", "path",
                  description="SDF/MOL or MOL2 file holding one record "
                              "(exactly one of smiles / ligand_file)"),
        ParamSpec("keep_largest", "bool", default=False,
                  description="When the input is several fragments, keep "
                              "the largest and name the removed ones in "
                              "the data"),
        ParamSpec("seed", "int",
                  description="Embedding seed (a fixed default, so the "
                              "geometry is reproducible)"),
        ParamSpec("output_name", "str", default=_OUTPUT_NAME,
                  description="File name of the prepared SDF in work_dir"),
    )
    data_keys = (
        DataKeySpec("n_atoms", "int", description="Atoms of the prepared "
                    "ligand, hydrogens included"),
        DataKeySpec("n_heavy_atoms", "int"),
        DataKeySpec("n_hydrogens", "int"),
        DataKeySpec("method", "str", description="ETDGv3 (embedded from a "
                    "SMILES) or coordinates-kept (file coordinates)"),
        DataKeySpec("seed", "int"),
        DataKeySpec("protonation", "str",
                    description="Always 'as-given': no pKa model is applied"),
        DataKeySpec("removed_fragments", "list",
                    description="Fragments removed by keep_largest, named "
                                "by formula and size; empty otherwise"),
        DataKeySpec("output_file", "str",
                    description="Name of the prepared SDF in work_dir"),
    )
    requires_python = ("rdkit",)
    wires = (("ligand", "ligand_file"),)

    def validate_params(self, **kwargs: Any) -> None:
        if not kwargs.get("smiles") and not kwargs.get("ligand_file"):
            raise ValueError(
                "exactly one of 'smiles' or 'ligand_file' is required"
            )
        _output_file_name(kwargs.get("output_name"))
        seed = kwargs.get("seed")
        if seed is not None and not int(seed) == float(seed):
            raise ValueError(f"'seed' must be an integer: {seed!r}")

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
            output_name = _output_file_name(kwargs.get("output_name"))
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )
        try:
            seed = kwargs.get("seed")
            seed = int(seed) if seed is not None else None
        except (TypeError, ValueError):
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"'seed' must be an integer: {kwargs.get('seed')!r}",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        options = dict(
            smiles=kwargs.get("smiles") or None,
            ligand_file=kwargs.get("ligand_file") or None,
            keep_largest=bool(kwargs.get("keep_largest", False)),
            seed=seed if seed is not None else None,
        )
        options = {k: v for k, v in options.items() if v is not None}

        try:
            prepared = prepare_ligand(**options)
        except StructureFormatError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.UNSUITABLE_INPUT,
            )

        out_path = Path(work_dir) / output_name
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_sdf(prepared, out_path)

        molecule = prepared.molecule
        data = {
            "n_atoms": molecule.n_atoms,
            "n_heavy_atoms": sum(1 for a in molecule.atoms
                                 if not a.is_hydrogen),
            "n_hydrogens": sum(1 for a in molecule.atoms if a.is_hydrogen),
            "method": prepared.method,
            "seed": prepared.seed,
            "protonation": "as-given",
            "removed_fragments": list(prepared.removed_fragments),
            "output_file": output_name,
            "n_records_refused": 0,
        }
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            output_file=out_path,
            data=data,
            artifacts={"ligand": out_path},
        )


register(LigandPrepareAdapter())
