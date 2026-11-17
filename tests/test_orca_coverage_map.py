"""The ORCA coverage map stays complete and grounded.

Three guarantees:

1. Schema and grounding: every record carries all required fields, a manual
   section and a probe verdict from :mod:`delfin.orca_knowledge.coverage`.
2. Nothing unexplained: every ``%block`` a DELFIN writer emits is mapped -- a
   new block in a workflow cannot appear without a coverage record.
3. The map is the source the builders trust (loadable, keys unique).
"""
from __future__ import annotations

import pathlib

import pytest

from delfin.orca_knowledge import (  # type: ignore
    CoverageError, covers_block, emitted_blocks, is_covered,
    load_coverage, record,
)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]

#: The writers that emit ORCA input as literal %blocks.  A block written here
#: but missing a coverage record fails test_every_emitted_block_is_mapped.
WRITER_FILES = [
    "delfin/common/orca_blocks.py",
    "delfin/common/solvation.py",
    "delfin/common/tddft_settings.py",
    "delfin/xyz_io.py",
    "delfin/esd_input_generator.py",
    "delfin/workflows/engine/classic.py",
    "delfin/workflows/engine/occupier.py",
    "delfin/build_up_complex.py",
    "delfin/calculators.py",
]

#: A record must name its purpose, a CONTROL key, at least one workflow, an OPI
#: answer and the file:line it came from.
_REQUIRED_FIELDS = frozenset({"key", "kind", "purpose", "manual", "probe",
                              "control", "opi", "workflows", "source"})
_KNOWN_PROBES = frozenset({"accepted", "rejected", "parsed", "finished",
                           "timeout", "refused", "unavailable"})


def test_coverage_file_loads_and_records_are_unique():
    data = load_coverage()
    records = data["records"]
    assert records, "coverage.json carries no records"
    # keys are unique by construction of a dict, but be explicit about the contract
    assert len(set(records)) == len(records)


def test_every_record_is_complete_and_grounded():
    """Every record names all fields, a manual section and a probe verdict."""
    for key, rec in load_coverage()["records"].items():
        assert rec["key"] == key, f"{key!r}: key field differs from dict key"
        missing = _REQUIRED_FIELDS - set(rec)
        assert not missing, f"{key!r} lacks {sorted(missing)}"
        assert str(rec["purpose"]).strip(), f"{key!r}: empty purpose"
        assert str(rec["manual"]).strip(), f"{key!r}: no manual section"
        assert rec["probe"] in _KNOWN_PROBES, f"{key!r}: probe {rec['probe']!r}"
        assert rec["kind"] in {"simple_keyword", "block_option"}, f"{key!r}: kind"
        assert str(rec["source"]).strip(), f"{key!r}: no file:line source"
        assert rec["workflows"], f"{key!r}: no workflows listed"


def test_every_emitted_block_is_mapped():
    """A %block written by a DELFIN writer must have a coverage record.

    This is the 'a new keyword in a workflow cannot appear unexplained' guard:
    any writer emitting a new ``%block ... end`` that no one documented fails
    the suite until a coverage record explains it.
    """
    coverage = load_coverage()["records"]
    uncovered: list[tuple[str, str]] = []
    for rel in WRITER_FILES:
        path = (REPO_ROOT / rel)
        if not path.exists():
            continue
        for block in emitted_blocks(path.read_text(encoding="utf-8")):
            # record key starts '%<block>' (option records name their block)
            if not any(k.strip().lower().startswith(block.lower() + " ")
                       or k.strip().lower() == block.lower()
                       for k in coverage):
                uncovered.append((rel, block))
    assert not uncovered, (
        "writers emit %blocks with no coverage record:\n" +
        "\n".join(f"  {rel}: {b}" for rel, b in uncovered) +
        "\nadd a record to delfin/orca_knowledge/coverage.json explaining it")


def test_a_writer_block_is_actually_found_by_the_grep():
    """The grep is not vacuous: it sees the %pal/%maxcore/%scf blocks orca_blocks writes."""
    path = REPO_ROOT / "delfin/common/orca_blocks.py"
    blocks = emitted_blocks(path.read_text(encoding="utf-8"))
    assert "%pal" in blocks
    assert "%maxcore" in blocks
    assert covers_block("%pal")
    assert covers_block("%scf")
    assert covers_block("%esd")
    assert covers_block("%tddft")


def test_unknown_block_is_not_covered():
    assert not covers_block("%notablock")
    assert is_covered("RIJCOSX")
    assert record("RIJCOSX")["kind"] == "simple_keyword"


def test_loader_rejects_a_corrupt_coverage_file(tmp_path):
    import delfin.orca_knowledge.coverage as cov

    original = cov.COVERAGE_FILE
    bad = tmp_path / "coverage.json"
    bad.write_text('{"recognized": "nothing really", "records": {}}', encoding="utf-8")
    cov.COVERAGE_FILE = bad  # type: ignore[assignment]
    load_coverage.cache_clear()
    try:
        with pytest.raises(CoverageError):
            cov.load_coverage()
    finally:
        cov.COVERAGE_FILE = original
        load_coverage.cache_clear()
