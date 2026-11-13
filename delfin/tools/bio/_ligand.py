"""Readers for small-molecule records: MDL SDF/MOL (V2000) and SYBYL MOL2.

Input: the text of an SDF, MOL or MOL2 file, which may hold several records.
Output: a tuple of :class:`~delfin.tools.bio._model.SmallMolecule`.

Semantics: both formats state up front how many atoms and bonds a record
holds, and both are routinely cut short by a failed download or a partial
write.  The declared count is therefore treated as a contract: a record that
promises three atoms and supplies two is refused, naming the count it
promised.  SDF atom records are read by column (the V2000 block is
column-oriented); MOL2 records are read by token, which is what that format
specifies.  The element of a MOL2 atom comes from its SYBYL type before the
dot ("C.ar" → "C"), never from the atom name, which is arbitrary.

Only V2000 is read.  V3000 writes its atoms in a tagged block with a
different layout; a V3000 counts line is refused by name rather than parsed
as if it were V2000, which would silently read zero atoms.

Rejected alternatives: trusting the record to end where the atom block stops
making sense (turns a truncated file into a smaller molecule); taking the
element from the atom name (names like "CL1" and "C12" are indistinguishable
from the chlorine they are not).
"""

from __future__ import annotations

import re
from typing import Dict, List, Tuple

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._model import LigandAtom, SmallMolecule, normalise_element

_SDF_TAG = re.compile(r"<([^>]*)>")
_MOL2_RECORD = "@<TRIPOS>MOLECULE"


def _counts(line: str, source: str, line_no: int) -> Tuple[int, int]:
    """Atom and bond count of a V2000 counts line (columns 1-3 and 4-6)."""
    if "V3000" in line.upper():
        raise StructureFormatError(
            f"{source}: line {line_no}: this is a V3000 record; this reader "
            "handles V2000 molfiles only"
        )
    try:
        return int(line[0:3]), int(line[3:6])
    except ValueError:
        raise StructureFormatError(
            f"{source}: line {line_no}: not a V2000 counts line: {line[:6]!r}"
        ) from None


def _read_sdf_record(lines: List[str], start: int, source: str
                     ) -> Tuple[SmallMolecule, int]:
    """One SDF record beginning at *start* → (molecule, index after the record)."""
    if start + 3 >= len(lines):
        raise StructureFormatError(
            f"{source}: the record starting on line {start + 1} has no counts "
            "line — the file ends inside its header"
        )
    name = lines[start].strip()
    n_atoms, n_bonds = _counts(lines[start + 3], source, start + 4)

    atoms: List[LigandAtom] = []
    index = start + 4
    for k in range(n_atoms):
        if index >= len(lines):
            raise StructureFormatError(
                f"{source}: the record starting on line {start + 1} promises "
                f"{n_atoms} atoms, the file ends after {k}"
            )
        line = lines[index]
        try:
            x, y, z = float(line[0:10]), float(line[10:20]), float(line[20:30])
        except ValueError:
            raise StructureFormatError(
                f"{source}: line {index + 1}: atom record {k + 1} of {n_atoms} "
                f"is not an atom line: {line[:34]!r}"
            ) from None
        symbol = line[31:34].strip()
        atoms.append(LigandAtom(name=f"{symbol}{k + 1}",
                                element=normalise_element(symbol),
                                x=x, y=y, z=z))
        index += 1

    bonds: List[Tuple[int, int, int]] = []
    for k in range(n_bonds):
        if index >= len(lines):
            raise StructureFormatError(
                f"{source}: the record starting on line {start + 1} promises "
                f"{n_bonds} bonds, the file ends after {k}"
            )
        line = lines[index]
        try:
            first, second = int(line[0:3]), int(line[3:6])
            order = int(line[6:9]) if line[6:9].strip() else 1
        except ValueError:
            raise StructureFormatError(
                f"{source}: line {index + 1}: bond record {k + 1} of {n_bonds} "
                f"is not a bond line: {line[:9]!r}"
            ) from None
        bonds.append((first - 1, second - 1, order))
        index += 1

    properties: Dict[str, str] = {}
    while index < len(lines):
        line = lines[index]
        if line.strip() == "$$$$":
            index += 1
            break
        if line.startswith(">"):
            match = _SDF_TAG.search(line)
            index += 1
            values: List[str] = []
            while index < len(lines) and lines[index].strip():
                values.append(lines[index].strip())
                index += 1
            if match:
                properties[match.group(1)] = "\n".join(values)
            continue
        index += 1

    return SmallMolecule(name=name, atoms=tuple(atoms), bonds=tuple(bonds),
                         properties=properties, source_format="sdf",
                         source_name=source), index


def parse_sdf(text: str, *, source_name: str = "<sdf>") -> Tuple[SmallMolecule, ...]:
    """Read every record of an SDF/MOL file.  Raises :class:`StructureFormatError`."""
    lines = text.splitlines()
    molecules: List[SmallMolecule] = []
    index = 0
    while index < len(lines):
        if not any(line.strip() for line in lines[index:]):
            break
        molecule, index = _read_sdf_record(lines, index, source_name)
        molecules.append(molecule)
    if not molecules:
        raise StructureFormatError(f"{source_name}: no molecule record in the file")
    return tuple(molecules)


def _read_mol2_record(lines: List[str], start: int, end: int, source: str
                      ) -> SmallMolecule:
    """One MOL2 record (lines[start:end], starting at its MOLECULE tag)."""
    body = lines[start:end]
    if len(body) < 3:
        raise StructureFormatError(
            f"{source}: line {start + 1}: the MOLECULE record has no name and "
            "counts line"
        )
    name = body[1].strip()
    fields = body[2].split()
    try:
        n_atoms = int(fields[0])
        n_bonds = int(fields[1]) if len(fields) > 1 else 0
    except (IndexError, ValueError):
        raise StructureFormatError(
            f"{source}: line {start + 3}: not a MOL2 counts line: {body[2]!r}"
        ) from None

    def section(tag: str) -> int:
        for offset, line in enumerate(body):
            if line.strip().upper() == f"@<TRIPOS>{tag}":
                return offset + 1
        return -1

    atom_start = section("ATOM")
    if atom_start < 0:
        raise StructureFormatError(
            f"{source}: the record starting on line {start + 1} has no "
            "@<TRIPOS>ATOM section"
        )
    atoms: List[LigandAtom] = []
    for k in range(n_atoms):
        offset = atom_start + k
        if offset >= len(body) or body[offset].startswith("@"):
            raise StructureFormatError(
                f"{source}: the record starting on line {start + 1} promises "
                f"{n_atoms} atoms, its @<TRIPOS>ATOM section holds {k}"
            )
        fields = body[offset].split()
        if len(fields) < 6:
            raise StructureFormatError(
                f"{source}: line {start + offset + 1}: atom record {k + 1} of "
                f"{n_atoms} has {len(fields)} fields, at least 6 are required"
            )
        try:
            x, y, z = float(fields[2]), float(fields[3]), float(fields[4])
        except ValueError:
            raise StructureFormatError(
                f"{source}: line {start + offset + 1}: atom record {k + 1} of "
                f"{n_atoms} has no three coordinates: {body[offset]!r}"
            ) from None
        charge = 0.0
        if len(fields) >= 9:
            try:
                charge = float(fields[8])
            except ValueError:
                charge = 0.0
        atoms.append(LigandAtom(name=fields[1],
                                element=normalise_element(fields[5].split(".")[0]),
                                x=x, y=y, z=z, charge=charge))

    bonds: List[Tuple[int, int, int]] = []
    bond_start = section("BOND")
    if bond_start >= 0:
        for k in range(n_bonds):
            offset = bond_start + k
            if offset >= len(body) or body[offset].startswith("@"):
                raise StructureFormatError(
                    f"{source}: the record starting on line {start + 1} promises "
                    f"{n_bonds} bonds, its @<TRIPOS>BOND section holds {k}"
                )
            fields = body[offset].split()
            if len(fields) < 4:
                raise StructureFormatError(
                    f"{source}: line {start + offset + 1}: bond record {k + 1} of "
                    f"{n_bonds} has {len(fields)} fields, at least 4 are required"
                )
            try:
                first, second = int(fields[1]), int(fields[2])
            except ValueError:
                raise StructureFormatError(
                    f"{source}: line {start + offset + 1}: bond record {k + 1} of "
                    f"{n_bonds} does not name two atoms: {body[offset]!r}"
                ) from None
            # Aromatic ("ar"), amide ("am") and unspecified ("un") orders have
            # no integer in this format; they are reported as order 0 rather
            # than rounded to a single or double bond.
            order = int(fields[3]) if fields[3].isdigit() else 0
            bonds.append((first - 1, second - 1, order))

    return SmallMolecule(name=name, atoms=tuple(atoms), bonds=tuple(bonds),
                         source_format="mol2", source_name=source)


def parse_mol2(text: str, *, source_name: str = "<mol2>") -> Tuple[SmallMolecule, ...]:
    """Read every record of a MOL2 file.  Raises :class:`StructureFormatError`."""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines)
              if line.strip().upper() == _MOL2_RECORD]
    if not starts:
        raise StructureFormatError(
            f"{source_name}: no {_MOL2_RECORD} record in the file"
        )
    bounds = starts + [len(lines)]
    return tuple(
        _read_mol2_record(lines, bounds[k], bounds[k + 1], source_name)
        for k in range(len(starts))
    )


__all__ = ["parse_sdf", "parse_mol2"]
