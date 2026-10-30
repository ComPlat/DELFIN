"""Tests for the shared ORCA keyword dictionary loader and its records.

Uses only the standard library plus the dependency-free loader
(delfin.orca_knowledge.dictionary) -- no OPI import. Verifies the schema,
field presence, text-length budget and the explain() contract.
"""
import json
import os

import delfin.orca_knowledge.dictionary as dictionary

DATADIR = os.path.dirname(os.path.abspath(dictionary.__file__))


def _all_records():
    """Flatten every record across a share's JSON files."""
    out = []
    for name in sorted(os.listdir(DATADIR)):
        if not name.endswith(".json"):
            continue
        with open(os.path.join(DATADIR, name), "r", encoding="utf-8") as fh:
            data = json.load(fh)
        out.extend(r for r in data.values() if isinstance(r, dict))
    return out


REQUIRED = ["keyword", "category", "purpose", "requires", "excludes",
            "pairs_with", "defaults", "typical_use", "literature", "manual",
            "parser", "opi", "status"]


def test_dictionary_files_exist():
    assert os.path.exists(os.path.join(DATADIR, "__init__.py")), \
        "shared loader __init__.py missing"
    jsons = [n for n in os.listdir(DATADIR) if n.endswith(".json")]
    assert jsons, "at least one share JSON must exist"


def test_every_record_has_required_fields():
    for rec in _all_records():
        assert isinstance(rec.get("manual"), dict), rec.get("keyword")
        for field in REQUIRED:
            assert field in rec, f"{rec.get('keyword')}: missing field {field}"


def test_record_text_within_budget():
    """No record text field longer than 400 chars (no pasted manual text)."""
    text_fields = ["purpose", "requires", "excludes", "pairs_with", "defaults",
                   "typical_use", "literature"]
    for rec in _all_records():
        for field in text_fields:
            val = rec.get(field)
            if val:
                assert len(str(val)) <= 400, \
                    f"{rec.get('keyword')}.{field} is {len(str(val))} chars"


def test_every_record_has_manual_section_and_parser():
    for rec in _all_records():
        # manual field must exist with a page string; every record has a parser verdict
        assert isinstance(rec["manual"], dict), rec.get("keyword")
        assert "page" in rec["manual"], rec.get("keyword")
        assert "section" in rec["manual"], rec.get("keyword")
        assert rec.get("parser") in ("accepted", "rejected", "unavailable"), \
            f"{rec.get('keyword')}: parser verdict {rec.get('parser')!r}"
        # a record marked grounded MUST point at a real manual section; draft
        # records are allowed an empty page (keyword not documented verbatim).
        if rec["status"] == "grounded":
            assert rec["manual"].get("page"), \
                f"grounded record {rec.get('keyword')} has no manual page"


def test_status_is_valid():
    for rec in _all_records():
        assert rec["status"] in ("grounded", "draft"), rec.get("keyword")


def test_explain_rijcosx_caseblind():
    rec = dictionary.explain("rijcosx")
    assert rec is not None
    assert rec["keyword"] == "rijcosx"
    assert dictionary.explain("RIJCOSX")["keyword"] == "rijcosx"


def test_explain_unknown_returns_none():
    assert dictionary.explain("zzz_no_such_keyword") is None


def test_records_filter_by_category():
    allr = dictionary.records()
    assert len(allr) > 0
    approx = dictionary.records("approximation")
    assert approx
    assert all(r["category"] == "approximation" for r in approx)


def test_no_duplicate_keywords():
    keys = [r["keyword"] for r in dictionary.records()]
    assert len(keys) == len(set(keys)), "duplicate keyword in dictionary"
