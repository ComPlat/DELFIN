"""A protein sequence as a FASTA document, read in full or refused by name.

Input: FASTA text (``parse_fasta``) or a path to a FASTA file
(``read_fasta``).  Output: :class:`SequenceRecord` — one identifier and one
sequence per record, in file order.

Semantics, and the three decisions that keep this reader honest:

* **Every record is read or none is.**  A FASTA file with one wrong letter in
  its fifth record is not four good sequences plus a footnote; a structure
  predictor run on it would report a model nobody asked about.  The reader
  raises :class:`~delfin.tools.bio.StructureFormatError` naming the record and
  the position, the same refusal the structure layer makes for a file it
  cannot read in full.
* **The alphabet is the 20 standard amino acids plus ``X``.**  The
  ambiguity/unknown code ``X`` is what every predictor accepts in place of a
  residue it cannot place; letters like ``B``, ``Z``, ``U`` and ``O`` are
  valid FASTA codes but not in this alphabet, and are refused like any other
  letter, with the position named.  Lowercase letters are read uppercase,
  which is what the format means by them.
* **DNA and RNA records are recognised and refused by name.**  ``A`` ``C``
  ``G`` ``T`` (and ``U``) are valid letters in both alphabets; a record made
  only of them is a nucleic-acid record, and feeding one to a protein
  predictor produces a plausible nonsense.  The refusal says DNA or RNA, not
  "bad letter", so the caller can point at the right file.  The first letter
  outside the DNA/RNA alphabet decides: a protein with one ``X`` is a protein
  with an unknown residue, never DNA.

Whitespace, blank lines and ``;``-comment lines are ignored, which is what
the format allows.  A record without a ``>`` header is refused: an
identifier per record is part of what the port promises.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Union

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._read import _text

#: The 20 standard amino acids plus the unknown-residue code ``X``.
STANDARD_AMINO_ACIDS = frozenset(
    "ACDEFGHIKLMNPQRSTVWY"
)  # noqa: F401  (the frozenset below is the alphabet this module enforces)

_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYX")

#: The DNA alphabet (plus ``U`` for RNA).  A record made only of these
#: letters is recognised as a nucleic-acid record and refused for protein
#: prediction.
DNA_ALPHABET = frozenset("ACGTURYKMSWBDHVN")

_PATHLIKE = Union[str, "Path"]


class SequenceRecord:
    """One FASTA record: its identifier line and its sequence, uppercased."""

    __slots__ = ("id", "sequence", "source_name")

    def __init__(self, record_id: str, sequence: str,
                 source_name: str = "") -> None:
        self.id = record_id
        self.sequence = sequence
        self.source_name = source_name

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return (f"SequenceRecord(id={self.id!r}, "
                f"length={len(self.sequence)})")


def _judge_record(record_id: str, letters: str, source: str) -> None:
    """Refuse a record whose sequence is not a protein sequence.

    A bad letter is refused with the position named.  DNA and RNA are decided
    by composition, not by single letters: ``A`` ``C`` ``G`` ``T`` are also
    the amino acids alanine, cysteine, glycine and threonine, so no single
    letter proves a nucleic-acid record.  A record whose every letter (the
    unknown code ``X`` aside) is drawn from the four DNA bases plus ``U`` is
    recognised as a nucleic-acid record and refused by name — a protein long
    enough to predict a structure from is never spelled in alanine, cysteine,
    glycine, threonine and unknowns alone.
    """
    for position, letter in enumerate(letters, start=1):
        upper = letter.upper()
        # ``U`` is judged by composition, not per letter: it is the one-letter
        # code of selenocysteine in a protein and of uracil in RNA.
        if upper not in _ALPHABET and upper != "U":
            raise StructureFormatError(
                f"{source}: record '{record_id}' position {position}: "
                f"'{letter}' is not a standard amino acid or X — the 20 "
                "standard amino acids plus X are the alphabet a protein "
                "prediction reads"
            )
    dna = DNA_ALPHABET
    if all(letter.upper() in dna or letter.upper() == "X" for letter in letters):
        raise StructureFormatError(
            f"{source}: record '{record_id}' is a DNA or RNA sequence, not a "
            "protein — a protein structure prediction cannot be run on it"
        )


def parse_fasta(text: str, *, source_name: str = "") -> Tuple[SequenceRecord, ...]:
    """Read FASTA *text* into records.  Raises :class:`StructureFormatError`."""
    source = source_name or "<fasta>"
    records: List[Tuple[str, List[str]]] = []
    current_id: str = ""
    current_letters: List[str] = []
    have_record = False
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if have_record:
                records.append((current_id, current_letters))
            header = line[1:].strip()
            if not header:
                raise StructureFormatError(
                    f"{source}: line {line_no}: a '>' record header without "
                    "an identifier — every record carries its id"
                )
            current_id = header
            current_letters = []
            have_record = True
        else:
            if not have_record:
                raise StructureFormatError(
                    f"{source}: line {line_no}: sequence letters before the "
                    "first '>' record header — every record carries its "
                    "identifier"
                )
            current_letters.append(line.upper())
    if have_record:
        records.append((current_id, current_letters))
    if not records:
        raise StructureFormatError(
            f"{source}: no FASTA records — the file holds no '>' header"
        )
    out = []
    for record_id, letters in records:
        joined = "".join(letters)
        _judge_record(record_id, joined, source)
        out.append(SequenceRecord(record_id, joined, source_name))
    return tuple(out)


def read_fasta(path: _PATHLIKE) -> Tuple[SequenceRecord, ...]:
    """Read every record of a FASTA file.  Raises :class:`StructureFormatError`."""
    text, name = _text(path)
    return parse_fasta(text, source_name=name)


__all__ = [
    "DNA_ALPHABET",
    "STANDARD_AMINO_ACIDS",
    "SequenceRecord",
    "parse_fasta",
    "read_fasta",
]
