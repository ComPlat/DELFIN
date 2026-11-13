"""Reader for the PDB coordinate format (column-oriented, format v3.3).

Input: the text of a PDB file (``ATOM`` / ``HETATM`` / ``MODEL`` / ``ENDMDL``
/ ``LINK`` / ``REMARK 465`` records).  Output: a
:class:`~delfin.tools.bio._model.MacromolecularStructure`.

Semantics: fields are taken from their **columns**, never by splitting on
whitespace.  The format allows every field to touch its neighbour — a
four-character atom name fills columns 13-16 with no space before the
residue name, a five-digit residue number touches the chain identifier — so a
whitespace split silently mis-assigns exactly the files that need care.  A
record shorter than 54 columns (the end of the z coordinate) is refused with
its line number; so is a ``MODEL`` block that the file never closes.

Rejected alternatives: splitting on whitespace (above); accepting a short
record by padding it with blanks (invents a coordinate of 0.0); reading the
element from the atom name even when columns 77-78 are filled (the name is
ambiguous, the element column is not).
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._model import (
    AtomRef,
    BioAtom,
    BondedLigandLink,
    Chain,
    MacromolecularStructure,
    MissingResidue,
    Residue,
    StructureModel,
    normalise_element,
)

# Element symbols that are two letters long, upper-cased for lookup.  Used
# only when the file carries no element column: a name starting in column 13
# is then a two-letter element unless it is a hydrogen name.
_TWO_LETTER_ELEMENTS = frozenset("""
HE LI BE NE NA MG AL SI CL AR CA SC TI CR MN FE CO NI CU ZN GA GE AS SE BR
KR RB SR ZR NB MO TC RU RH PD AG CD IN SN SB TE XE CS BA LA CE PR ND PM SM
EU GD TB DY HO ER TM YB LU HF TA RE OS IR PT AU HG TL PB BI PO AT RN FR RA
AC TH PA NP PU AM CM BK CF ES FM MD NO LR RF DB SG BH HS MT DS RG CN NH FL
MC LV TS OG
""".split())

# "<number><optional insertion code>", the way REMARK 465 writes a residue.
_SEQ_ICODE = re.compile(r"^(-?\d+)([A-Za-z]?)$")

_MIN_COORD_COLUMNS = 54          # through the end of the z field


def _infer_element(name_field: str) -> str:
    """Element symbol for an atom whose record carries no element column.

    Input: the raw four-column atom-name field (columns 13-16), spaces kept.
    Output: the canonical element symbol, or "" when the field is blank.

    The format places a one-letter element in column 14 and a two-letter one
    in column 13, so the indentation carries the information.  Two exceptions
    are common enough to matter: a hydrogen name may be led by a digit
    ("1HB") and a four-character hydrogen name starts in column 13 ("HG11"),
    which would otherwise read as hafnium.
    """
    field = (name_field + "    ")[:4]
    stripped = field.strip()
    if not stripped:
        return ""
    if stripped[0].isdigit():
        for ch in stripped[1:]:
            if ch.isalpha():
                return normalise_element(ch)
        return ""
    if field[0] == " ":
        return normalise_element(stripped[0])
    if len(stripped) == 4 and stripped[0].upper() in ("H", "D"):
        return normalise_element(stripped[0])
    candidate = field[:2].strip().upper()
    if candidate in _TWO_LETTER_ELEMENTS:
        return normalise_element(candidate)
    return normalise_element(stripped[0])


def _float(text: str, default: float) -> float:
    try:
        return float(text)
    except ValueError:
        return default


class _ModelBuilder:
    """Groups ATOM/HETATM records of one model into chains and residues.

    Residues are keyed by (chain, sequence number, insertion code, name,
    record type) and kept in the order they first appear, so the structure
    reads back in file order even when records of one chain are interleaved.
    """

    def __init__(self, model_id: int) -> None:
        self.model_id = model_id
        self._residues: Dict[Tuple, List[BioAtom]] = {}
        self._chain_order: List[str] = []

    def __bool__(self) -> bool:
        return bool(self._residues)

    def add(self, key: Tuple, atom: BioAtom) -> None:
        if key not in self._residues:
            self._residues[key] = []
            if key[0] not in self._chain_order:
                self._chain_order.append(key[0])
        self._residues[key].append(atom)

    def build(self) -> StructureModel:
        per_chain: Dict[str, List[Residue]] = {c: [] for c in self._chain_order}
        for (chain_id, seq_id, icode, name, is_hetero), atoms in self._residues.items():
            per_chain[chain_id].append(
                Residue(chain_id=chain_id, seq_id=seq_id, insertion_code=icode,
                        name=name, is_hetero=is_hetero, atoms=tuple(atoms))
            )
        return StructureModel(
            model_id=self.model_id,
            chains=tuple(Chain(chain_id=c, residues=tuple(per_chain[c]))
                         for c in self._chain_order),
        )


def _parse_atom(line: str, line_no: int, source: str) -> Tuple[Tuple, BioAtom]:
    """One ATOM/HETATM record → (residue key, atom).  Raises on a short line."""
    if len(line) < _MIN_COORD_COLUMNS:
        raise StructureFormatError(
            f"{source}: line {line_no}: {line[:6].strip()} record is truncated — "
            f"{len(line)} columns, the z coordinate ends at column "
            f"{_MIN_COORD_COLUMNS}"
        )
    padded = line.ljust(80)
    name_field = padded[12:16]
    try:
        x = float(padded[30:38])
        y = float(padded[38:46])
        z = float(padded[46:54])
    except ValueError:
        raise StructureFormatError(
            f"{source}: line {line_no}: coordinate columns 31-54 are not three "
            f"numbers: {padded[30:54]!r}"
        ) from None
    try:
        serial = int(padded[6:11])
    except ValueError:
        serial = 0
    element = normalise_element(padded[76:78]) or _infer_element(name_field)
    is_hetero = padded.startswith("HETATM")
    try:
        seq_id = int(padded[22:26])
    except ValueError:
        raise StructureFormatError(
            f"{source}: line {line_no}: residue sequence number columns 23-26 "
            f"are not a number: {padded[22:26]!r}"
        ) from None
    atom = BioAtom(
        serial=serial,
        name=name_field.strip(),
        element=element,
        x=x, y=y, z=z,
        altloc=padded[16].strip(),
        occupancy=_float(padded[54:60], 1.0),
        b_factor=_float(padded[60:66], 0.0),
        is_hetero=is_hetero,
    )
    key = (padded[21].strip(), seq_id, padded[26].strip(),
           padded[17:20].strip(), is_hetero)
    return key, atom


def _parse_link(line: str) -> Optional[BondedLigandLink]:
    """One LINK record → a bonded-ligand record, or None when unreadable."""
    padded = line.ljust(80)

    def seq(cols: slice) -> Optional[int]:
        try:
            return int(padded[cols])
        except ValueError:
            return None

    seq1, seq2 = seq(slice(22, 26)), seq(slice(52, 56))
    if seq1 is None or seq2 is None:
        return None
    distance: Optional[float]
    try:
        distance = float(padded[73:78])
    except ValueError:
        distance = None
    return BondedLigandLink(
        atom1=AtomRef(chain_id=padded[21].strip(), residue_name=padded[17:20].strip(),
                      seq_id=seq1, insertion_code=padded[26].strip(),
                      atom_name=padded[12:16].strip(), altloc=padded[16].strip()),
        atom2=AtomRef(chain_id=padded[51].strip(), residue_name=padded[47:50].strip(),
                      seq_id=seq2, insertion_code=padded[56].strip(),
                      atom_name=padded[42:46].strip(), altloc=padded[46].strip()),
        distance=distance,
        kind="link",
    )


def _parse_missing(payload: str) -> Optional[MissingResidue]:
    """One REMARK 465 table row → a missing residue, or None when it is prose.

    The table is written as ``[model] RES C SSSEQI``; the sequence number and
    the insertion code are one token ("2A").  A row that does not have that
    shape is not a residue and is skipped rather than guessed at.
    """
    tokens = payload.split()
    model_id = 1
    if len(tokens) == 4 and tokens[0].isdigit():
        model_id, tokens = int(tokens[0]), tokens[1:]
    if len(tokens) != 3:
        return None
    name, chain, seq_icode = tokens
    match = _SEQ_ICODE.match(seq_icode)
    if match is None or not name.isalnum():
        return None
    return MissingResidue(name=name, chain_id=chain, seq_id=int(match.group(1)),
                          insertion_code=match.group(2), model_id=model_id)


def parse_pdb(text: str, *, source_name: str = "<pdb>") -> MacromolecularStructure:
    """Read PDB *text* into a structure.  Raises :class:`StructureFormatError`."""
    models: List[StructureModel] = []
    builder = _ModelBuilder(1)
    open_model_line = 0
    next_implicit_id = 1
    missing: List[MissingResidue] = []
    links: List[BondedLigandLink] = []
    in_missing_table = False

    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r\n")
        record = line[:6]
        if record.startswith(("ATOM  ", "HETATM")):
            key, atom = _parse_atom(line, line_no, source_name)
            builder.add(key, atom)
        elif record.startswith("MODEL "):
            if open_model_line:
                raise StructureFormatError(
                    f"{source_name}: line {line_no}: MODEL starts while the model "
                    f"opened on line {open_model_line} has no ENDMDL"
                )
            if builder:
                models.append(builder.build())
            try:
                model_id = int(line[10:14])
            except ValueError:
                model_id = next_implicit_id
            builder = _ModelBuilder(model_id)
            open_model_line = line_no
            next_implicit_id = model_id + 1
        elif record.startswith("ENDMDL"):
            if not open_model_line:
                raise StructureFormatError(
                    f"{source_name}: line {line_no}: ENDMDL without a MODEL"
                )
            models.append(builder.build())
            builder = _ModelBuilder(next_implicit_id)
            open_model_line = 0
        elif record.startswith("LINK"):
            link = _parse_link(line)
            if link is not None:
                links.append(link)
        elif line.startswith("REMARK 465"):
            payload = line[10:].strip()
            if "SSSEQI" in payload:
                in_missing_table = True
            elif in_missing_table and payload:
                entry = _parse_missing(payload)
                if entry is not None:
                    missing.append(entry)

    if open_model_line:
        raise StructureFormatError(
            f"{source_name}: the MODEL opened on line {open_model_line} is never "
            "closed by ENDMDL — the file stops inside it"
        )
    if builder:
        models.append(builder.build())
    if not models:
        raise StructureFormatError(
            f"{source_name}: no atom records (ATOM/HETATM) in the file"
        )
    return MacromolecularStructure(
        models=tuple(models),
        missing_residues=tuple(missing),
        bonded_ligands=tuple(links),
        source_format="pdb",
        source_name=source_name,
    )


__all__ = ["parse_pdb"]
