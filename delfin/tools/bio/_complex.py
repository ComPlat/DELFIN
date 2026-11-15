"""A complex descriptor: one JSON file naming a macromolecule and its ligands.

Input: the path of a JSON descriptor whose keys are ``macromolecule`` (the
file name of a PDB or mmCIF structure) and ``ligands`` (a list of small-molecule
file names).  Output: a :class:`Complex` whose ``macromolecule`` and
``ligands`` are resolved :class:`~pathlib.Path` objects, in descriptor order.

Semantics: the file names are relative and resolved against the directory the
descriptor itself lives in, so a complex can be moved as a folder and the
names still point at the same files.  An absolute name, or a relative name
whose resolution leaves that directory, is refused as
:class:`~delfin.tools.bio.StructureFormatError`: a complex that reaches
outside its own folder is not portable, and a receptor resolved from
somewhere else is a plausible-looking answer to a question nobody asked.
Every named file must exist.  The reader opens and type-checks the descriptor
only — the named structure files are resolved, not parsed, so a caller can
decide which reader (``read_structure``, ``read_small_molecule``) to apply.

Rejected alternatives: resolving against the process working directory
(breaks the moment the descriptor is read from anywhere else); accepting
absolute names (makes the descriptor machine-bound).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Union

from delfin.tools.bio._errors import StructureFormatError

_PathLike = Union[str, Path]


@dataclass(frozen=True)
class Complex:
    """A macromolecule plus the small molecules bound (or to be bound) to it.

    ``macromolecule`` and every entry of ``ligands`` are absolute
    :class:`~pathlib.Path` objects inside the descriptor's directory.
    """

    macromolecule: Path
    ligands: Tuple[Path, ...] = ()
    source_name: str = ""

    @property
    def descriptor_directory(self) -> Path:
        """The directory the descriptor lives in; the root of the complex."""
        return self.macromolecule.parent


def _relative_inside(name: str, source: str, directory: Path) -> Path:
    """Resolve *name* against *directory*, refusing anything outside it."""
    if Path(name).is_absolute():
        raise StructureFormatError(
            f"{source}: {name!r} is absolute; a complex descriptor names its "
            "files relative to its own directory"
        )
    resolved = (directory / name).resolve()
    root = directory.resolve()
    if root != resolved and root not in resolved.parents:
        raise StructureFormatError(
            f"{source}: {name!r} resolves outside the descriptor's directory"
        )
    return resolved


def read_complex(path: _PathLike) -> Complex:
    """Read a complex descriptor.  Raises :class:`StructureFormatError`."""
    p = Path(path)
    try:
        text = p.read_text(encoding="utf-8")
    except OSError as exc:
        raise StructureFormatError(
            f"{p.name}: cannot be read ({exc.strerror})"
        ) from None
    try:
        descriptor = json.loads(text)
    except json.JSONDecodeError as exc:
        raise StructureFormatError(
            f"{p.name}: not a JSON complex descriptor ({exc.msg} at line "
            f"{exc.lineno})"
        ) from None
    if not isinstance(descriptor, dict):
        raise StructureFormatError(
            f"{p.name}: a complex descriptor is a JSON object with the keys "
            "'macromolecule' and 'ligands'"
        )
    source = p.name
    macromolecule = descriptor.get("macromolecule")
    if not isinstance(macromolecule, str) or not macromolecule.strip():
        raise StructureFormatError(
            f"{source}: the key 'macromolecule' is missing or empty; a "
            "complex names the macromolecular structure file it is about"
        )
    ligands = descriptor.get("ligands", [])
    if ligands is None:
        ligands = []
    if isinstance(ligands, str) or not isinstance(ligands, list):
        raise StructureFormatError(
            f"{source}: 'ligands' must be a list of file names, "
            f"got {type(ligands).__name__}"
        )

    directory = p.resolve().parent
    macro_path = _relative_inside(macromolecule, source, directory)
    if not macro_path.is_file():
        raise StructureFormatError(
            f"{source}: the named macromolecule file is not there: "
            f"{macromolecule}"
        )
    resolved_ligands = []
    for name in ligands:
        if not isinstance(name, str) or not name.strip():
            raise StructureFormatError(
                f"{source}: every entry of 'ligands' is a file name, "
                f"got {name!r}"
            )
        ligand_path = _relative_inside(name, source, directory)
        if not ligand_path.is_file():
            raise StructureFormatError(
                f"{source}: the named ligand file is not there: {name}"
            )
        resolved_ligands.append(ligand_path)

    return Complex(
        macromolecule=macro_path,
        ligands=tuple(resolved_ligands),
        source_name=source,
    )


__all__ = ["Complex", "read_complex"]
