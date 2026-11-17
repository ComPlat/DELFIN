"""Tests for the solvation/workflows share (solvation_workflows.json).

Solvation, semiempirical/xTB, force fields, QM/MM, GOAT/NEB/IRC/DOCKER/
SOLVATOR/fragments/EDA keywords and block options: every record must carry
the shared record contract, a manual section, a parser verdict, and (for
grounded records) a verified quote. Runs on the CI interpreter with no OPI;
never starts ORCA.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pytest

from delfin.orca_knowledge.dictionary import explain

SHARE_FILE = Path(__file__).resolve().parent.parent / (
    "delfin/orca_knowledge/dictionary/solvation_workflows.json"
)

MANUAL_ROOT = Path(os.environ.get(
    "DELFIN_ORCA_MANUAL_DIR",
    str(Path.home() / "delfin-knowledge" / "orca" / "manual-6.1-text"
        / "contents"),
))

QUOTE_FILE = Path(__file__).resolve().parent.parent / (
    "scratch/quotes_solvation_workflows.json"
)

SHARE_RECORDS = json.loads(SHARE_FILE.read_text())
SHARE_BY_ID = {r["id"]: r for r in SHARE_RECORDS}

REQUIRED = ("id", "keyword", "category", "purpose", "manual", "parser",
            "opi", "status")
TEXT_FIELDS = ("purpose", "typical_use", "defaults", "requires", "excludes",
               "pairs_with")


def test_share_has_records():
    assert len(SHARE_RECORDS) >= 200, "the share is unexpectedly small"


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
        if isinstance(value, list):
            assert all(isinstance(x, str) and len(x) <= 400 for x in value), \
                f"text field {field!r} has an overlong item in {record['id']}"
    assert record["parser"] in ("accepted", "rejected", "unavailable")
    assert record["status"] in ("grounded", "draft")
    assert record["manual"].get("page") and record["manual"].get("section"), \
        f"manual section missing in {record['id']}"


@pytest.mark.parametrize("record", SHARE_RECORDS, ids=lambda r: r["id"])
def test_no_pasted_manual_text(record):
    """Purpose fields must be our own words, not manual prose (>=10 word
    verbatim runs of a manual page are rejected)."""
    page = MANUAL_ROOT / (record["manual"]["page"] + ".txt")
    if not page.exists():
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
    for word in ("cpcm", "SMD18", "DRACO", "AM1", "XTB", "GOAT", "NEB-TS",
                 "IRC", "QMMM", "QUICKDOCK", "EXTOPT", "GFN-FF", "NATIVE-GFN-XTB"):
        rec = explain(word)
        assert rec is not None, f"explain({word!r}) failed"
    # aliases resolve case-blind
    assert explain("xtb2") is not None
    assert explain("c-pcm") is not None
    # block options are addressable
    assert explain("epsilon") is not None
    assert explain("nimages") is not None


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
    entries = quotes.get(rid)
    assert entries, f"grounded record {rid} has no quotes"
    for entry in entries:
        page = entry["page"]
        if page not in pages:
            p = MANUAL_ROOT / (page + ".txt")
            assert p.exists(), f"quote page missing for {rid}: {page}"
            pages[page] = norm(p.read_text(errors="replace"))
        assert norm(entry["quote"]) in pages[page], \
            f"quote not on {page} for {rid}"
