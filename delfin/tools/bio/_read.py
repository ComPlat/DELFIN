"""File-level entry points: open a path, pick the reader, name the failure.

Input: a path to a structure (``.pdb``, ``.ent``, ``.cif``, ``.mmcif``) or to
a small molecule (``.sdf``, ``.mol``, ``.mol2``).  Output: the parsed object
from the matching reader in this package.

Semantics: the format is taken from the file-name suffix, not sniffed from
the content.  A suffix is what the producer of the file declared; sniffing
guesses, and a wrong guess on these formats yields a plausible structure
rather than an error.  Every failure — a missing file, an unreadable byte, an
unknown suffix — is raised as :class:`StructureFormatError`, so a caller
needs one ``except`` clause for "this file cannot be used".

Rejected alternative: returning ``None`` for an unreadable file, which makes
the caller's next attribute access the error message.
"""

from __future__ import annotations

from pathlib import Path
from typing import Tuple, Union

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._ligand import parse_mol2, parse_sdf
from delfin.tools.bio._mmcif import parse_mmcif
from delfin.tools.bio._model import MacromolecularStructure, SmallMolecule
from delfin.tools.bio._pdb import parse_pdb

PDB_SUFFIXES = (".pdb", ".ent", ".pdb1")
MMCIF_SUFFIXES = (".cif", ".mmcif", ".pdbx")
SDF_SUFFIXES = (".sdf", ".sd", ".mol")
MOL2_SUFFIXES = (".mol2",)

_PathLike = Union[str, Path]


def _text(path: _PathLike) -> Tuple[str, str]:
    """(file text, file name).  Raises :class:`StructureFormatError` on any
    read failure, so a missing or unreadable file fails like a malformed one."""
    p = Path(path)
    try:
        return p.read_text(encoding="utf-8", errors="replace"), p.name
    except OSError as exc:
        raise StructureFormatError(f"{p.name}: cannot be read ({exc.strerror})") from None


def read_pdb(path: _PathLike) -> MacromolecularStructure:
    """Read a PDB file."""
    text, name = _text(path)
    return parse_pdb(text, source_name=name)


def read_mmcif(path: _PathLike) -> MacromolecularStructure:
    """Read an mmCIF (PDBx) file."""
    text, name = _text(path)
    return parse_mmcif(text, source_name=name)


def read_structure(path: _PathLike) -> MacromolecularStructure:
    """Read a macromolecular structure, choosing the reader by suffix."""
    suffix = Path(path).suffix.lower()
    if suffix in PDB_SUFFIXES:
        return read_pdb(path)
    if suffix in MMCIF_SUFFIXES:
        return read_mmcif(path)
    raise StructureFormatError(
        f"{Path(path).name}: {suffix or 'no suffix'} is not a macromolecular "
        f"structure format; expected one of "
        f"{', '.join(PDB_SUFFIXES + MMCIF_SUFFIXES)}"
    )


def read_sdf(path: _PathLike) -> Tuple[SmallMolecule, ...]:
    """Read every record of an SDF/MOL file."""
    text, name = _text(path)
    return parse_sdf(text, source_name=name)


def read_mol2(path: _PathLike) -> Tuple[SmallMolecule, ...]:
    """Read every record of a MOL2 file."""
    text, name = _text(path)
    return parse_mol2(text, source_name=name)


def read_small_molecule(path: _PathLike) -> Tuple[SmallMolecule, ...]:
    """Read a small-molecule file, choosing the reader by suffix."""
    suffix = Path(path).suffix.lower()
    if suffix in SDF_SUFFIXES:
        return read_sdf(path)
    if suffix in MOL2_SUFFIXES:
        return read_mol2(path)
    raise StructureFormatError(
        f"{Path(path).name}: {suffix or 'no suffix'} is not a small-molecule "
        f"format; expected one of {', '.join(SDF_SUFFIXES + MOL2_SUFFIXES)}"
    )


__all__ = [
    "PDB_SUFFIXES",
    "MMCIF_SUFFIXES",
    "SDF_SUFFIXES",
    "MOL2_SUFFIXES",
    "read_pdb",
    "read_mmcif",
    "read_structure",
    "read_sdf",
    "read_mol2",
    "read_small_molecule",
]
