"""Reader for the mmCIF (PDBx) coordinate format.

Input: the text of an mmCIF file.  Output: a
:class:`~delfin.tools.bio._model.MacromolecularStructure` built from the same
rules as the PDB reader, so the two formats of one system read as one system.

Semantics: a small, complete-enough CIF 1.1 tokenizer (quoted values,
semicolon text fields, ``loop_`` tables, comments) feeds a category table;
the structure is then assembled from ``_atom_site``, with
``_pdbx_unobs_or_zero_occ_residues`` for the residues the depositor reports
as unobserved and ``_struct_conn`` for connection records.  Chain, residue
number, residue name and atom name are read from the ``auth_*`` items with
the ``label_*`` items as fallback: the ``auth_*`` numbering is the one the
PDB file of the same entry uses and the one literature and users quote, so
reading ``label_seq_id`` would silently renumber every residue.  ``.`` and
``?`` are CIF nulls, not an altloc "." or an insertion code "?".

A ``loop_`` whose last row is incomplete is refused — a file cut off there
yields a structure missing its last atoms, which nothing downstream can see.

Rejected alternatives: a regular expression over ``_atom_site`` lines (breaks
on quoted values and on any item order but the usual one); reading
``label_*`` numbering (renumbers the system); treating ``.``/``?`` as data
(invents altlocs and insertion codes that then split residues).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._model import (
    AtomRef,
    _apply_modifications,
    BioAtom,
    BondedLigandLink,
    Chain,
    MacromolecularStructure,
    MissingResidue,
    Residue,
    StructureModel,
    normalise_element,
)

# Connection types that attach a ligand or residue to another one.  A
# hydrogen bond or a salt bridge is a contact, not an attachment, and is not
# reported as a bonded ligand.
_BONDED_CONN_TYPES = frozenset({"covale", "metalc", "disulf", "covale_base",
                                "covale_phosphate", "covale_sugar"})

_CIF_KEYWORDS = frozenset({"loop_", "stop_", "global_", "save_"})
_NULLS = ("", ".", "?")

_Token = Tuple[str, bool, int]          # (value, is_literal, line number)


def _split_line(text: str) -> List[Tuple[str, bool]]:
    """Split one CIF line into (value, is_literal) tokens, honouring quotes.

    A quoted value is literal: it never starts a tag or a keyword even when
    it begins with an underscore.  A ``#`` outside quotes starts a comment.
    """
    out: List[Tuple[str, bool]] = []
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch.isspace():
            i += 1
            continue
        if ch == "#":
            break
        if ch in "'\"":
            quote, j, buf = ch, i + 1, []
            while j < n:
                if text[j] == quote and (j + 1 >= n or text[j + 1].isspace()):
                    break
                buf.append(text[j])
                j += 1
            out.append(("".join(buf), True))
            i = j + 1
            continue
        j = i
        while j < n and not text[j].isspace():
            j += 1
        out.append((text[i:j], False))
        i = j
    return out


def _tokenize(text: str, source: str) -> List[_Token]:
    """The whole file as tokens.  Raises on an unterminated text field."""
    lines = text.splitlines()
    tokens: List[_Token] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            start = i + 1
            buf, closed = [line[1:]], False
            i += 1
            while i < len(lines):
                if lines[i].startswith(";"):
                    closed = True
                    i += 1
                    break
                buf.append(lines[i])
                i += 1
            if not closed:
                raise StructureFormatError(
                    f"{source}: the text field (';') opened on line {start} is "
                    "never closed — the file stops inside it"
                )
            tokens.append(("\n".join(buf).strip(), True, start))
            continue
        for value, literal in _split_line(line):
            tokens.append((value, literal, i + 1))
        i += 1
    return tokens


def _categories(text: str, source: str) -> Dict[str, List[Dict[str, str]]]:
    """Parse CIF *text* into ``{category: [row, …]}``.

    Key-value items of one category collapse into a single row, which is how
    a non-looped category differs from a looped one in this format.
    """
    tokens = _tokenize(text, source)
    out: Dict[str, List[Dict[str, str]]] = {}
    i, n = 0, len(tokens)
    while i < n:
        value, literal, line_no = tokens[i]
        if not literal and value.lower() == "loop_":
            i += 1
            names: List[str] = []
            while i < n and not tokens[i][1] and tokens[i][0].startswith("_"):
                names.append(tokens[i][0])
                i += 1
            if not names:
                raise StructureFormatError(
                    f"{source}: line {line_no}: loop_ without column names"
                )
            values: List[str] = []
            while i < n:
                tok, lit, _ = tokens[i]
                if not lit and (tok.startswith("_")
                                or tok.lower() in _CIF_KEYWORDS
                                or tok.lower().startswith("data_")):
                    break
                values.append(tok)
                i += 1
            category = names[0].split(".")[0]
            if len(values) % len(names):
                raise StructureFormatError(
                    f"{source}: the {category} loop ends in the middle of a row: "
                    f"{len(values)} values for {len(names)} columns "
                    f"({len(values) % len(names)} left over)"
                )
            rows = out.setdefault(category, [])
            for start in range(0, len(values), len(names)):
                rows.append(dict(zip(names, values[start:start + len(names)])))
            continue
        if not literal and value.startswith("_"):
            if i + 1 >= n:
                raise StructureFormatError(
                    f"{source}: line {line_no}: item {value} has no value"
                )
            category = value.split(".")[0]
            rows = out.setdefault(category, [])
            if not rows:
                rows.append({})
            rows[0][value] = tokens[i + 1][0]
            i += 2
            continue
        i += 1
    return out


def _item(row: Dict[str, str], *names: str, default: str = "") -> str:
    """First non-null value among *names*; CIF nulls count as absent."""
    for name in names:
        value = row.get(name, "")
        if value not in _NULLS:
            return value
    return default


def _int(row, names, source, what, default=None) -> int:
    text = _item(row, *names)
    if text in _NULLS and default is not None:
        return default
    try:
        return int(text)
    except ValueError:
        raise StructureFormatError(
            f"{source}: {what} is not a number: {text!r}"
        ) from None


def _float_item(row, names, source, what, default=None) -> float:
    text = _item(row, *names)
    if text in _NULLS and default is not None:
        return default
    try:
        return float(text)
    except ValueError:
        raise StructureFormatError(
            f"{source}: {what} is not a number: {text!r}"
        ) from None


def _build_models(rows, source) -> Tuple[StructureModel, ...]:
    """Group ``_atom_site`` rows into models, chains and residues (file order)."""
    per_model: Dict[int, Dict[Tuple, List[BioAtom]]] = {}
    model_order: List[int] = []
    chain_order: Dict[int, List[str]] = {}
    for index, row in enumerate(rows, start=1):
        model_id = _int(row, ("_atom_site.pdbx_PDB_model_num",), source,
                        f"atom {index}: model number", default=1)
        is_hetero = _item(row, "_atom_site.group_PDB",
                          default="ATOM").upper() == "HETATM"
        seq_id = _int(row, ("_atom_site.auth_seq_id", "_atom_site.label_seq_id"),
                      source, f"atom {index}: residue sequence number")
        atom = BioAtom(
            serial=_int(row, ("_atom_site.id",), source,
                        f"atom {index}: id", default=index),
            name=_item(row, "_atom_site.auth_atom_id", "_atom_site.label_atom_id"),
            element=normalise_element(_item(row, "_atom_site.type_symbol")),
            x=_float_item(row, ("_atom_site.Cartn_x",), source, f"atom {index}: x"),
            y=_float_item(row, ("_atom_site.Cartn_y",), source, f"atom {index}: y"),
            z=_float_item(row, ("_atom_site.Cartn_z",), source, f"atom {index}: z"),
            altloc=_item(row, "_atom_site.label_alt_id", "_atom_site.pdbx_auth_alt_id"),
            occupancy=_float_item(row, ("_atom_site.occupancy",), source,
                                  f"atom {index}: occupancy", default=1.0),
            b_factor=_float_item(row, ("_atom_site.B_iso_or_equiv",), source,
                                 f"atom {index}: B factor", default=0.0),
            is_hetero=is_hetero,
        )
        key = (
            _item(row, "_atom_site.auth_asym_id", "_atom_site.label_asym_id"),
            seq_id,
            _item(row, "_atom_site.pdbx_PDB_ins_code"),
            _item(row, "_atom_site.auth_comp_id", "_atom_site.label_comp_id"),
            is_hetero,
        )
        if model_id not in per_model:
            per_model[model_id] = {}
            chain_order[model_id] = []
            model_order.append(model_id)
        residues = per_model[model_id]
        if key not in residues:
            residues[key] = []
            if key[0] not in chain_order[model_id]:
                chain_order[model_id].append(key[0])
        residues[key].append(atom)

    models: List[StructureModel] = []
    for model_id in model_order:
        per_chain: Dict[str, List[Residue]] = {c: [] for c in chain_order[model_id]}
        for (chain_id, seq_id, icode, name, is_hetero), atoms in per_model[model_id].items():
            per_chain[chain_id].append(
                Residue(chain_id=chain_id, seq_id=seq_id, insertion_code=icode,
                        name=name, is_hetero=is_hetero, atoms=tuple(atoms))
            )
        models.append(StructureModel(
            model_id=model_id,
            chains=tuple(Chain(chain_id=c, residues=tuple(per_chain[c]))
                         for c in chain_order[model_id]),
        ))
    return tuple(models)


def _missing_residues(categories, source) -> Tuple[MissingResidue, ...]:
    out: List[MissingResidue] = []
    for row in categories.get("_pdbx_unobs_or_zero_occ_residues", []):
        name = _item(row, "_pdbx_unobs_or_zero_occ_residues.auth_comp_id",
                     "_pdbx_unobs_or_zero_occ_residues.label_comp_id")
        if not name:
            continue
        out.append(MissingResidue(
            name=name,
            chain_id=_item(row, "_pdbx_unobs_or_zero_occ_residues.auth_asym_id",
                           "_pdbx_unobs_or_zero_occ_residues.label_asym_id"),
            seq_id=_int(row, ("_pdbx_unobs_or_zero_occ_residues.auth_seq_id",
                              "_pdbx_unobs_or_zero_occ_residues.label_seq_id"),
                        source, f"unobserved residue {name}: sequence number"),
            insertion_code=_item(row, "_pdbx_unobs_or_zero_occ_residues.PDB_ins_code"),
            model_id=_int(row, ("_pdbx_unobs_or_zero_occ_residues.PDB_model_num",),
                          source, "unobserved residue: model number", default=1),
        ))
    return tuple(out)


def _bonded_links(categories, source) -> Tuple[BondedLigandLink, ...]:
    out: List[BondedLigandLink] = []
    for row in categories.get("_struct_conn", []):
        kind = _item(row, "_struct_conn.conn_type_id").lower()
        if kind not in _BONDED_CONN_TYPES:
            continue

        def ref(side: str) -> Optional[AtomRef]:
            name = _item(row, f"_struct_conn.ptnr{side}_auth_comp_id",
                         f"_struct_conn.ptnr{side}_label_comp_id")
            seq_text = _item(row, f"_struct_conn.ptnr{side}_auth_seq_id",
                             f"_struct_conn.ptnr{side}_label_seq_id")
            try:
                seq_id = int(seq_text)
            except ValueError:
                return None
            return AtomRef(
                chain_id=_item(row, f"_struct_conn.ptnr{side}_auth_asym_id",
                               f"_struct_conn.ptnr{side}_label_asym_id"),
                residue_name=name,
                seq_id=seq_id,
                insertion_code=_item(row, f"_struct_conn.ptnr{side}_PDB_ins_code"),
                atom_name=_item(row, f"_struct_conn.ptnr{side}_label_atom_id",
                                f"_struct_conn.ptnr{side}_auth_atom_id"),
                altloc=_item(row, f"_struct_conn.pdbx_ptnr{side}_label_alt_id"),
            )

        first, second = ref("1"), ref("2")
        if first is None or second is None:
            continue
        distance: Optional[float]
        text = _item(row, "_struct_conn.pdbx_dist_value")
        try:
            distance = float(text)
        except ValueError:
            distance = None
        out.append(BondedLigandLink(atom1=first, atom2=second,
                                    distance=distance, kind=kind))
    return tuple(out)


def _modified_residues(categories) -> Dict[Tuple[str, int, str], Tuple[str, str]]:
    """``_pdbx_struct_mod_residue`` rows → {(asym, seq, icode): (comp, parent)}.

    The category keys a modification by label numbering; auth numbering is
    used when present because the coordinates are keyed by auth numbering
    (see :func:`_build_models`).  A row without a sequence number is skipped:
    nothing to map to.
    """
    out: Dict[Tuple[str, int, str], Tuple[str, str]] = {}
    for row in categories.get("_pdbx_struct_mod_residue", []):
        comp = _item(row, "_pdbx_struct_mod_residue.label_comp_id",
                     "_pdbx_struct_mod_residue.SIF_name")
        parent = _item(row, "_pdbx_struct_mod_residue.parent_comp_id")
        chain = _item(row, "_pdbx_struct_mod_residue.auth_asym_id",
                      "_pdbx_struct_mod_residue.label_asym_id")
        seq_text = _item(row, "_pdbx_struct_mod_residue.auth_seq_id",
                         "_pdbx_struct_mod_residue.label_seq_id")
        icode = _item(row, "_pdbx_struct_mod_residue.PDB_ins_code")
        if not (comp and parent) or seq_text in _NULLS:
            continue
        try:
            seq = int(seq_text)
        except ValueError:
            continue
        out[(chain, seq, icode)] = (comp, parent)
    return out


def parse_mmcif(text: str, *, source_name: str = "<mmcif>") -> MacromolecularStructure:
    """Read mmCIF *text* into a structure.  Raises :class:`StructureFormatError`."""
    categories = _categories(text, source_name)
    rows = categories.get("_atom_site")
    if not rows:
        raise StructureFormatError(
            f"{source_name}: no _atom_site category — the file carries no coordinates"
        )
    models = _build_models(rows, source_name)
    modifications = _modified_residues(categories)
    models = [_apply_modifications(m, modifications) for m in models]
    entry = categories.get("_entry", [{}])
    return MacromolecularStructure(
        models=models,
        missing_residues=_missing_residues(categories, source_name),
        bonded_ligands=_bonded_links(categories, source_name),
        entry_id=entry[0].get("_entry.id", "") if entry else "",
        source_format="mmcif",
        source_name=source_name,
    )


__all__ = ["parse_mmcif"]
