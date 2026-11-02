"""Tests for the correlated-wavefunction share (wft_correlated.json).

MP2/RI-MP2/DLPNO/CC/MDCI/CASSCF/NEVPT2/CASPT2/AVAS/AUTOCI/ICE keywords and
block options: every record must carry the shared record contract, a manual
section, a parser verdict, and (for grounded records) a verified quote.
Runs on the CI interpreter with no OPI; never starts ORCA.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from delfin.orca_knowledge.dictionary import explain

SHARE_FILE = Path(__file__).resolve().parent.parent / (
    "delfin/orca_knowledge/dictionary/wft_correlated.json"
)

MANUAL_ROOT = Path(
    "/home/qmchem_max/delfin-knowledge/orca/manual-6.1-text"
)

QUOTE_FILE = Path(__file__).resolve().parent.parent / (
    "scratch/quotes_wft_correlated.json"
)

SHARE_RECORDS = json.loads(SHARE_FILE.read_text())
SHARE_BY_ID = {r["id"]: r for r in SHARE_RECORDS}

REQUIRED = ("id", "keyword", "category", "purpose", "manual", "parser",
            "opi", "status")
TEXT_FIELDS = ("purpose", "typical_use", "defaults", "requires", "excludes",
               "pairs_with")


def test_share_has_records():
    assert len(SHARE_RECORDS) >= 50, "the wft share is unexpectedly small"


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
    for word in ("mp2", "RI-MP2", "dlpno-ccsd(t)", "casscf", "NEVPT2",
                 "avas(valence-d)", "TightPNO", "LED", "HFLD"):
        rec = explain(word)
        assert rec is not None, f"explain({word!r}) failed"
    # aliases resolve
    assert explain("RIMP2") is not None
    assert explain("OO-RIMP2") is not None
    # block options are addressable
    assert explain("TCutPNO") is not None


@pytest.mark.parametrize("rid", sorted(SHARE_BY_ID), ids=str)
def test_grounded_records_have_verified_quotes(rid):
    record = SHARE_BY_ID[rid]
    if record["status"] != "grounded":
        return
    assert QUOTE_FILE.exists(), "share quotes file missing"
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


def test_avas_alias_group_is_one_record():
    rec = explain("avas(valence-f)")
    assert rec is not None and "avas(double-f)" in rec["aliases"]
