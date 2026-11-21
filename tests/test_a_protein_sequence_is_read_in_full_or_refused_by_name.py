"""A protein sequence enters DELFIN as its own file port, read in full.

Structure prediction is driven by a sequence, not a structure file, and a
sequence is its own document: a FASTA file.  The quantum-chemistry blocks
compose because one emits ``gbw`` and the next declares it consumes it; a
prediction workflow needs the same for ``sequence`` — the port is declared in
``delfin.tools._spec`` and read through :mod:`delfin.tools.bio`, so a block
that predicts connects to an upstream block that reads without either knowing
the other.

The reader is honest about what it accepts: the 20 standard amino acids plus
``X``, one identifier per record, every record read or none.  A letter outside
the alphabet is refused with the record and the position named — a plausible
sequence with one wrong letter is worse than a refusal, because a predictor
would run on it and report a structure nobody asked about.  DNA and RNA
records are recognised and refused for protein prediction, by name, for the
same reason: the letters are valid FASTA, they are just not a protein.

Licences of the programs this port is the input to are named in
:mod:`delfin.tools.adapters.structure_prediction`.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.tools import _spec, describe
from delfin.tools.bio import (SequenceRecord, StructureFormatError,
                              parse_fasta, read_fasta)


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------


def _write(tmp_path: Path, text: str, name: str = "query.fasta") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_one_record_with_an_id_is_read():
    records = parse_fasta(">sp|P0DTC2|SPIKE spike protein\nMFVFLVLLPLV\n")
    assert len(records) == 1
    record = records[0]
    assert record.id == "sp|P0DTC2|SPIKE spike protein"
    assert record.sequence == "MFVFLVLLPLV"
    assert record.source_name == ""


def test_multiple_records_keep_their_order_and_ids():
    records = parse_fasta(
        ">chain A lysozyme\nKVFGRCELAA\n\n>chain B trypsin\nIVGGYTCGAN\n"
    )
    assert [r.id for r in records] == ["chain A lysozyme", "chain B trypsin"]
    assert [r.sequence for r in records] == ["KVFGRCELAA", "IVGGYTCGAN"]


def test_blank_lines_and_whitespace_inside_a_sequence_are_ignored():
    records = parse_fasta(">rec\nKVF\n\n  GRC\nELAA\n")
    assert records[0].sequence == "KVFGRCELAA"


def test_a_record_without_an_id_is_refused():
    with pytest.raises(StructureFormatError) as excinfo:
        parse_fasta("KVFGRCELAA\n")
    message = str(excinfo.value)
    assert "record header" in message
    assert "identifier" in message


def test_a_letter_outside_the_alphabet_names_the_record_and_the_position():
    with pytest.raises(StructureFormatError) as excinfo:
        parse_fasta(">good\nKVFGRCELAA\n>bad\nKVFGR1ELAA\n")
    message = str(excinfo.value)
    assert "bad" in message
    assert "position 6" in message   # position of '1' inside the sequence
    assert "'1'" in message          # the letter itself


def test_lowercase_letters_are_read_uppercase():
    records = parse_fasta(">rec\nkvfgrc\n")
    assert records[0].sequence == "KVFGRCELAA"[:6]


def test_the_file_is_refused_when_it_holds_nothing():
    with pytest.raises(StructureFormatError):
        parse_fasta("\n\n  \n")


# ---------------------------------------------------------------------------
# DNA / RNA records are recognised and refused for protein prediction
# ---------------------------------------------------------------------------


def test_a_dna_record_is_refused_by_name():
    with pytest.raises(StructureFormatError) as excinfo:
        parse_fasta(">gene fragment\nATGCTGATCGA\n")
    assert "DNA" in str(excinfo.value)
    assert "gene fragment" in str(excinfo.value)


def test_an_rna_record_is_refused_by_name():
    with pytest.raises(StructureFormatError) as excinfo:
        parse_fasta(">mrna\nAUGCAUGCAU\n")
    assert "RNA" in str(excinfo.value)


def test_a_record_that_is_mostly_dna_with_one_gap_letter_is_dna_not_protein():
    """A DNA alphabet hit on the first letters decides; X does not rescue it."""
    with pytest.raises(StructureFormatError) as excinfo:
        parse_fasta(">orf\nATGCATGCX\n")
    assert "DNA" in str(excinfo.value)


def test_dna_is_refused_only_for_whole_records_of_dna_letters():
    """'DNASTAR' etc. inside an identifier does not make a protein DNA."""
    records = parse_fasta(">PDB_DNA_binding_protein\nKVFGRCELAA\n")
    assert records[0].sequence == "KVFGRCELAA"


# ---------------------------------------------------------------------------
# the file reader
# ---------------------------------------------------------------------------


def test_read_fasta_takes_the_text_from_a_path(tmp_path):
    path = _write(tmp_path, ">rec\nKVFGRCELAA\n", name="query.fasta")
    records = read_fasta(path)
    assert records[0].id == "rec"
    assert records[0].source_name == "query.fasta"


def test_read_fasta_refuses_a_missing_file(tmp_path):
    with pytest.raises(StructureFormatError):
        read_fasta(tmp_path / "nope.fasta")


# ---------------------------------------------------------------------------
# the port
# ---------------------------------------------------------------------------


def test_the_sequence_port_is_declared_in_the_capability_vocabulary():
    """A ``sequence`` block consumes a sequence file the ``gbw`` way."""
    text = Path(_spec.__file__).read_text(encoding="utf-8")
    assert "sequence" in text
    assert "FASTA" in text


def test_a_prediction_block_declares_it_consumes_sequence():
    contract = describe("structure_prediction_engine")
    assert "sequence" in contract.consumes
