"""Tests for the basis/ECP/relativistic share (basis_ecp_relativistic.json).

Orbital/auxiliary basis-set families, ECP keywords, basis-option keywords,
the %basis and %rel blocks, and the relativistic simple keywords: every
record must carry the shared record contract, a manual section, a parser
verdict, and (for grounded records) a verified quote. Runs on the CI
interpreter with no OPI; never starts ORCA.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pytest

from delfin.orca_knowledge.dictionary import explain, records

SHARE_FILE = Path(__file__).resolve().parent.parent / (
    "delfin/orca_knowledge/dictionary/basis_ecp_relativistic.json"
)

# The licensed manual text lives outside the repository; tests that need it
# skip where it is absent (CI).
MANUAL_ROOT = Path(os.environ.get(
    "DELFIN_ORCA_MANUAL_DIR",
    str(Path.home() / "delfin-knowledge" / "orca" / "manual-6.1-text"),
))

QUOTE_FILE = Path(__file__).resolve().parent.parent / (
    "scratch/quotes_basis_ecp_relativistic.json"
)

SHARE_RECORDS = json.loads(SHARE_FILE.read_text())
SHARE_BY_ID = {r["id"]: r for r in SHARE_RECORDS}

REQUIRED = ("id", "keyword", "category", "purpose", "manual", "parser",
            "opi", "status")
TEXT_FIELDS = ("purpose", "typical_use", "defaults", "requires", "excludes",
               "pairs_with")


def test_share_has_records():
    assert len(SHARE_RECORDS) >= 50, "the basis share is unexpectedly small"


@pytest.mark.parametrize("record", SHARE_RECORDS, ids=lambda r: r["id"])
def test_record_contract(record):
    for field in REQUIRED:
        assert record.get(field) not in (None, ""), \
            f"missing/empty {field!r} in {record['id']}"
    for field in TEXT_FIELDS:
        value = record.get(field)
        if isinstance(value, str):
            assert len(value) <= 400, \
                f"text field {field!r} too long in {record['id']}"
    assert record["parser"] in ("accepted", "rejected", "unavailable")
    assert record["status"] in ("grounded", "draft")


@pytest.mark.parametrize("record", SHARE_RECORDS, ids=lambda r: r["id"])
def test_no_pasted_manual_text(record):
    """Purpose fields must be our own words, not manual prose (>=10 word
    verbatim runs of a manual page are rejected)."""
    page = MANUAL_ROOT / (record["manual"]["page"] + ".txt")
    if not page.exists():  # bibliography-only quotes are legal
        return
    page_norm = re.sub(r"\s+", " ", page.read_text(errors="replace")).lower()
    purpose = re.sub(r"\s+", " ", record["purpose"]).lower()
    words = purpose.split()
    for start in range(0, max(1, len(words) - 10)):
        run = " ".join(words[start:start + 10])
        if len(run) > 40 and run in page_norm:
            pytest.fail(f"record {record['id']}: purpose contains a verbatim "
                        f"manual run: {run!r}")


def test_explain_works_for_share_items():
    # case-blind resolution of representative simple keywords
    for word in ("def2-svp", "def2/J", "cc-pVTZ", "SARC", "AutoAux",
                 "zora", "X2C", "sarc2-zora-qzvp", "x2c-tzvpall-s",
                 "somf(1x)", "uncontractaux"):
        rec = explain(word)
        assert rec is not None, f"explain({word!r}) failed"
    # block options are addressable
    for word in ("NewGTO", "FragBasis", "PCDTrimAuxJ", "PictureChange",
                 "ScaleZORA"):
        assert explain(word) is not None, f"explain({word!r}) failed"


@pytest.mark.parametrize("rid", sorted(SHARE_BY_ID), ids=str)
def test_grounded_records_have_verified_quotes(rid):
    record = SHARE_BY_ID[rid]
    if record["status"] != "grounded":
        return
    if not (QUOTE_FILE.exists() and MANUAL_ROOT.is_dir()):
        pytest.skip("manual text or share quotes file not on this machine")
    quotes = json.loads(QUOTE_FILE.read_text())
    norm = lambda s: re.sub(r"\s+", " ", s).strip().lower()  # noqa: E731
    pages: dict[str, str] = {}

    def page_text(page: str) -> str:
        rel = page + ".txt"
        if rel not in pages:
            pages[rel] = norm((MANUAL_ROOT / rel).read_text(errors="replace"))
        return pages[rel]

    for entry in quotes.get(rid, []):
        if norm(entry["quote"]) in page_text(entry["page"]):
            return
    pytest.fail(f"grounded record {rid} has no verified quote")


def test_rejected_parser_verdicts_are_explained():
    for record in SHARE_RECORDS:
        if record["parser"] == "rejected":
            assert record.get("purpose"), record["id"]


def test_rel_alias_group_is_one_record():
    rec = explain("SARC-ZORA-TZVP")
    assert rec is not None and "sarc-dkh-tzvp" in \
        [a.lower() for a in rec["aliases"]]


def test_records_category_accessor():
    cats = {r["category"] for r in SHARE_RECORDS}
    for cat in cats:
        got = records(cat)
        assert got, f"records({cat!r}) returned nothing"
        assert all(r["category"] == cat for r in got)
