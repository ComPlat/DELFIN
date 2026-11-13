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

import gzip
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

#: Most bytes a ``*.gz`` structure file may inflate to before it is refused.
#:
#: The PDB's largest entry is a fraction of a gigabyte uncompressed; a
#: structure file that decompresses past 1 GiB is either not a structure or a
#: decompression bomb (a 10 kB ``.gz`` can expand by a factor of ~1000), and
#: reading it whole would hang the Pi this runs on.  The cap is a hard
#: constant, not a parameter, because a caller that wants more than 1 GiB of
#: structure text is the caller to ask first.
MAX_DECOMPRESSED_BYTES = 1 << 30

_PathLike = Union[str, Path]


def _text(path: _PathLike) -> Tuple[str, str]:
    """(file text, file name).  Raises :class:`StructureFormatError` on any
    read failure, so a missing or unreadable file fails like a malformed one.

    A name ending in ``.gz`` is decompressed with :mod:`gzip` from the
    standard library; the suffix before ``.gz`` is the format's and the
    readers choose on it.  A stream that inflates beyond
    :data:`MAX_DECOMPRESSED_BYTES` is refused with the file's name, before
    the whole text is held in memory."""
    p = Path(path)
    name = p.name
    if name.lower().endswith(".gz"):
        return _text_gz(p, name)
    try:
        return p.read_text(encoding="utf-8", errors="replace"), name
    except OSError as exc:
        raise StructureFormatError(f"{name}: cannot be read ({exc.strerror})") from None


def _text_gz(p: Path, name: str) -> Tuple[str, str]:
    """Decompress *p* under the byte cap, naming the file on every failure."""
    try:
        with gzip.open(p, "rb") as fh:
            data = fh.read(MAX_DECOMPRESSED_BYTES + 1)
    except (OSError, gzip.BadGzipFile, EOFError) as exc:
        raise StructureFormatError(
            f"{name}: cannot be read as gzip ({exc})"
        ) from None
    if len(data) > MAX_DECOMPRESSED_BYTES:
        raise StructureFormatError(
            f"{name}: decompresses to more than "
            f"{MAX_DECOMPRESSED_BYTES} bytes, the cap for a gzipped "
            "structure file; the file is refused before it is read whole"
        )
    return data.decode("utf-8", errors="replace"), name


def read_pdb(path: _PathLike) -> MacromolecularStructure:
    """Read a PDB file."""
    text, name = _text(path)
    return parse_pdb(text, source_name=name)


def read_mmcif(path: _PathLike) -> MacromolecularStructure:
    """Read an mmCIF (PDBx) file."""
    text, name = _text(path)
    return parse_mmcif(text, source_name=name)


def _format_suffix(path: _PathLike) -> str:
    """The format suffix of a structure path: the suffix before ``.gz``.

    ``x.pdb.gz`` → ``.pdb``; a plain ``x.pdb`` → ``.pdb``.  A bare ``.gz``
    (``x.gz``) stays ``.gz``, which no reader list holds and so it is
    refused with the same message as any other unknown suffix.
    """
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".gz":
        stem = p.stem
        return Path(stem).suffix.lower()
    return suffix


def read_structure(path: _PathLike) -> MacromolecularStructure:
    """Read a macromolecular structure, choosing the reader by suffix."""
    suffix = _format_suffix(path)
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
    suffix = _format_suffix(path)
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
