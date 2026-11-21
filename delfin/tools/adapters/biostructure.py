"""Composition of a macromolecular structure file, as a building block.

Input: a PDB or mmCIF file (``structure_file``), optionally a model number.
Output: a :class:`~delfin.tools._types.StepResult` whose ``data`` holds the
composition as flat numbers and whose ``json_report`` artifact holds the full
parsed summary — per chain, the hetero components, the residues the file
reports as unobserved, and its connection records.

Semantics: this block answers the question every biomolecular workflow asks
before it starts, and that no file answers by being opened: how many chains,
how many residues, what else is in there, is it protonated, what does the
depositor say is missing, is anything covalently attached.  It is also the
first consumer of the ``macromolecule`` capability port, so a later block that
prepares or fetches a receptor connects to it without either side being
changed.

It runs no external program and reads through :mod:`delfin.tools.bio`, which
has no dependencies, so the step is ready on every machine DELFIN runs on —
``requires_binaries`` and ``requires_python`` are empty and the environment
probe reports it ready.

Failures are kinds, not prose: a file that is not there is ``MISSING_INPUT``,
a file (or a suffix) that cannot be read as a structure is ``PARSE``, and a
model the file does not contain is ``MISSING_INPUT`` naming the model numbers
that exist — a composition for the wrong model would be a plausible answer to
a question nobody asked.

Truncation is refused where the file shows it: a record cut inside its
coordinates, a ``MODEL`` never closed, an mmCIF loop ending mid-row. A file
cut exactly at a line boundary carries no such sign and reads as a smaller
structure; neither format marks its own end reliably enough to demand it.

Rejected alternatives: reporting a partial composition when the file shows
it is truncated (the numbers would be indistinguishable from a complete
file's);
defaulting the model to 1 (a file whose models start elsewhere would fail for
no reason, so the default is the file's first model whatever its number);
writing the summary into ``data`` only (a report is the thing a later block
and a reader want to keep, which is what an artifact is for).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.bio import (
    MacromolecularStructure,
    StructureFormatError,
    StructureTooLargeError,
    read_structure,
)


def _chain_rows(structure: MacromolecularStructure,
                model_id: Optional[int]) -> List[Dict[str, Any]]:
    """Per-chain counts, in the order the chains appear in the file."""
    model = (structure.first_model if model_id is None
             else structure.model(model_id))
    rows = []
    for chain in model.chains:
        rows.append({
            "chain_id": chain.chain_id,
            "n_polymer_residues": len(chain.polymer_residues),
            "n_hetero_residues": len(chain.hetero_residues),
            "n_waters": len(chain.waters),
            "n_atoms": len(chain.atoms),
            "n_hydrogens": sum(1 for a in chain.atoms if a.is_hydrogen),
        })
    return rows


def _report(structure: MacromolecularStructure,
            model_id: Optional[int]) -> Dict[str, Any]:
    """The full parsed summary written as the ``json_report`` artifact."""
    composition = structure.composition(model_id)
    return {
        "source_name": structure.source_name,
        "source_format": structure.source_format,
        "entry_id": structure.entry_id,
        "composition": composition.to_dict(),
        "chains": _chain_rows(structure, model_id),
        "missing_residues": [
            {"label": m.label, "chain_id": m.chain_id, "seq_id": m.seq_id,
             "insertion_code": m.insertion_code, "name": m.name,
             "model_id": m.model_id}
            for m in structure.missing_residues
        ],
        "bonded_ligands": [
            {"kind": b.kind, "distance": b.distance,
             "atom1": {"label": b.atom1.label, "chain_id": b.atom1.chain_id,
                       "residue_name": b.atom1.residue_name,
                       "seq_id": b.atom1.seq_id,
                       "insertion_code": b.atom1.insertion_code,
                       "atom_name": b.atom1.atom_name},
             "atom2": {"label": b.atom2.label, "chain_id": b.atom2.chain_id,
                       "residue_name": b.atom2.residue_name,
                       "seq_id": b.atom2.seq_id,
                       "insertion_code": b.atom2.insertion_code,
                       "atom_name": b.atom2.atom_name}}
            for b in structure.bonded_ligands
        ],
    }


def _report_file_name(value: Any) -> str:
    """The report's file name, or raise ``ValueError``.

    Input: the ``report_name`` parameter. Output: a plain file name, the
    default one when the value is empty. Semantics: the report is written
    into ``work_dir`` and nowhere else. ``work_dir / "/abs/x"`` is
    ``/abs/x`` and ``work_dir / "../x"`` leaves it, so a separator, ``.``
    or ``..`` would let a parameter choose any file to overwrite; they are
    refused rather than stripped, because a silently renamed report is
    one the caller then cannot find.
    """
    name = str(value or "").strip() or "structure_summary.json"
    if (name in (".", "..") or "/" in name or "\\" in name
            or Path(name).name != name):
        raise ValueError(
            f"'report_name' must be a file name inside work_dir, not a path: "
            f"{name!r}")
    return name


def _model_number(value: Any) -> Optional[int]:
    """The requested model number, ``None`` for the file's first model, or
    raise ``ValueError``. Only an integer or its decimal text is a model
    number: ``int(1.5)`` would quietly report model 1."""
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        raise ValueError(f"'model' must be an integer model number: {value!r}")
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip().lstrip("-").isdigit():
        return int(value.strip())
    raise ValueError(f"'model' must be an integer model number: {value!r}")


class BiostructureSummaryAdapter(StepAdapter):
    """Read a macromolecular structure file and report what it holds."""

    name = "biostructure_summary"
    description = (
        "Composition of a macromolecular structure file (PDB/mmCIF): chains, "
        "residue counts, hetero components, waters, missing residues, bonded "
        "ligand records and whether hydrogens are present"
    )
    category = "structure"
    produces_geometry = False
    consumes = ("macromolecule",)
    produces = ("json_report",)
    params = (
        ParamSpec("structure_file", "path", required=True,
                  description="Macromolecular structure file to read (PDB or mmCIF)"),
        ParamSpec("model", "int",
                  description="Model number to report; the file's first model "
                              "when not given"),
        ParamSpec("report_name", "str", default="structure_summary.json",
                  description="File name of the JSON report written into work_dir"),
    )
    data_keys = (
        DataKeySpec("n_models", "int", description="Models in the file"),
        DataKeySpec("model_id", "int", description="Model the counts describe"),
        DataKeySpec("n_chains", "int"),
        DataKeySpec("n_polymer_residues", "int",
                    description="Residues written as ATOM records, waters excluded"),
        DataKeySpec("n_modified_residues", "int",
                    description="Polymer residues declared modified "
                                "(MODRES / _pdbx_struct_mod_residue)"),
        DataKeySpec("modified_residues", "list",
                    description="Modified residues: chain, sequence id, "
                                "insertion code, own name and parent name"),
        DataKeySpec("n_hetero_residues", "int",
                    description="Hetero components, waters excluded"),
        DataKeySpec("hetero_components", "dict",
                    description="Hetero component name -> count"),
        DataKeySpec("n_waters", "int"),
        DataKeySpec("n_atoms", "int", description="Atom records, alternate "
                                                  "locations included"),
        DataKeySpec("n_atoms_primary", "int",
                    description="Atoms of the primary conformation"),
        DataKeySpec("n_hydrogens", "int"),
        DataKeySpec("has_hydrogens", "bool"),
        DataKeySpec("n_missing_residues", "int",
                    description="Residues the file reports as unobserved"),
        DataKeySpec("n_bonded_ligand_links", "int",
                    description="Connection records attaching a residue or ligand"),
        DataKeySpec("source_format", "str", description="pdb or mmcif"),
    )
    # A later block that hands on a receptor file fills this step's input.
    wires = (("macromolecule", "structure_file"),)

    def validate_params(self, **kwargs: Any) -> None:
        if not kwargs.get("structure_file"):
            raise ValueError(
                "'structure_file' parameter is required (a PDB or mmCIF file)"
            )
        _report_file_name(kwargs.get("report_name"))
        _model_number(kwargs.get("model"))

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

        raw = kwargs.get("structure_file")
        if not raw:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="'structure_file' is required (a PDB or mmCIF file)",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        try:
            report_name = _report_file_name(kwargs.get("report_name"))
            model_id = _model_number(kwargs.get("model"))
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
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
        except StructureTooLargeError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.UNSUITABLE_INPUT,
            )
        except StructureFormatError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.PARSE,
            )

        try:
            composition = structure.composition(model_id)
        except KeyError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc).strip("'"), error_kind=ErrorKind.MISSING_INPUT,
            )

        report = _report(structure, model_id)
        report_path = Path(work_dir) / report_name
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=False),
                               encoding="utf-8")

        data = {
            "n_models": composition.n_models,
            "model_id": composition.model_id,
            "n_chains": len(composition.chain_ids),
            "chain_ids": list(composition.chain_ids),
            "n_polymer_residues": composition.n_polymer_residues,
            "n_modified_residues": composition.n_modified_residues,
            "modified_residues": [
                {"chain_id": c, "seq_id": s, "insertion_code": i,
                 "name": n, "parent_name": p}
                for c, s, i, n, p in composition.modified_residues
            ],
            "n_hetero_residues": composition.n_hetero_residues,
            "hetero_components": dict(composition.hetero_components),
            "n_waters": composition.n_waters,
            "n_atoms": composition.n_atoms,
            "n_atoms_primary": composition.n_atoms_primary,
            "n_hydrogens": composition.n_hydrogens,
            "has_hydrogens": composition.has_hydrogens,
            "n_missing_residues": composition.n_missing_residues,
            "n_bonded_ligand_links": composition.n_bonded_ligand_links,
            "altloc_residues": list(composition.altloc_residues),
            "source_format": structure.source_format,
            "source_name": structure.source_name,
        }
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            output_file=report_path,
            data=data,
            artifacts={"json_report": report_path},
        )


register(BiostructureSummaryAdapter())
