"""Tests for the ORCA keyword dictionary loader.

The loader (``delfin/orca_knowledge/dictionary``) merges every share's
``<share>.json`` into one record set and exposes ``explain`` (case-blind,
alias-resolved) and ``records`` (category filter). This is a shared artifact:
any share may add records, but the record contract stays the same.

These tests run on the CI interpreter with no OPI dependency and never start
ORCA (the parser verdicts are recorded data, not live probes).
"""
from __future__ import annotations

import pytest

from delfin.orca_knowledge.dictionary import explain, records

REQUIRED_FIELDS = (
    "id", "keyword", "category", "purpose", "manual", "parser", "opi", "status",
)

TEXT_FIELDS = ("purpose", "typical_use", "defaults", "requires", "excludes",
               "pairs_with")


def pytest_generate_tests(metafunc):
    # Every record is parameterised: a schema addition is tested on all records.
    if "record" in metafunc.fixturenames:
        metafunc.parametrize("record", records(), ids=lambda r: r["id"])


def test_records_present():
    """The dictionary is non-empty and every record is addressable."""
    all_recs = records()
    assert len(all_recs) >= 10
    # a record known to belong to the properties share must be present
    assert any(r["id"] == "property.nmr" for r in all_recs)


def test_required_fields_present(record):
    for field in REQUIRED_FIELDS:
        assert field in record, f"missing field {field!r} in {record['id']}"
        assert record[field] not in (None, ""), \
            f"empty required field {field!r} in {record['id']}"


def test_manual_has_page_and_section(record):
    manual = record["manual"]
    assert isinstance(manual, dict)
    assert manual.get("page"), f"no manual.page in {record['id']}"
    assert manual.get("section"), f"no manual.section in {record['id']}"


def test_parser_verdict_is_known(record):
    assert record["parser"] in ("accepted", "rejected", "unavailable"), \
        f"bad parser verdict in {record['id']}"


def test_status_is_known(record):
    assert record["status"] in ("grounded", "draft"), \
        f"bad status in {record['id']}"


def test_no_record_text_field_exceeds_400_chars(record):
    for field in TEXT_FIELDS:
        value = record.get(field)
        if isinstance(value, str):
            assert len(value) <= 400, \
                f"text field {field!r} too long ({len(value)} ch) in {record['id']}"


def test_explain_is_case_insensitive():
    rec = explain("NMR")
    assert rec is not None and rec["keyword"].lower() == "nmr"


def test_explain_unknown_word_returns_none():
    assert explain("this-is-not-an-orca-keyword-xyzzy") is None


def test_explain_resolves_aliases():
    # any record with an alias list resolves the alias to the master keyword
    aliased = [r for r in records() if len(r.get("aliases", [])) > 1]
    if not aliased:
        pytest.skip("no share has aliased records yet")
    rec = aliased[0]
    alias = rec["aliases"][1]
    resolved = explain(alias)
    assert resolved is not None and resolved["id"] == rec["id"]


def test_records_category_filter():
    categories = {r["category"] for r in records()}
    assert categories, "no categories present"
    one = next(iter(categories))
    subset = records(category=one)
    assert subset, "category filter returned nothing for a live category"
    assert all(r["category"] == one for r in subset)
