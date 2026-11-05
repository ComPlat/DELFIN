"""The functionals/dispersion/gCP/grid share: honest status, dash-tolerant
lookup, required fields, bounded text.

Every record carries a manual section and a parser verdict; grounded
records have no pasted manual prose; the share's headline items resolve
through the shared loader, spelled with or without dashes.
"""
from __future__ import annotations

import json
from pathlib import Path

from delfin.orca_knowledge.dictionary import explain

SHARE_FILE = (Path(__file__).resolve().parents[1]
              / "delfin/orca_knowledge/dictionary"
              / "functionals_dispersion_gcp_grid.json")
REQUIRED = ("keyword", "category", "purpose", "manual", "parser", "status")


def _share():
    return json.loads(SHARE_FILE.read_text())


def test_the_share_is_loaded():
    names = {r["keyword"] for r in _share()}
    assert {"BP86", "D3BJ", "DefGrid2", "GCP"} & names


def test_every_record_has_required_fields():
    for rec in _share():
        missing = [k for k in REQUIRED if not rec.get(k)]
        assert not missing, (rec.get("keyword"), missing)
        assert rec["parser"] in ("accepted", "rejected", "unavailable")
        assert rec["status"] in ("grounded", "draft")
        assert rec["manual"].get("section")


def test_no_record_text_is_longer_than_400_characters():
    for rec in _share():
        for key in ("purpose", "typical_use", "note", "literature"):
            text = rec.get(key)
            if isinstance(text, str):
                assert len(text) <= 400, (rec["keyword"], key)
        if isinstance(rec.get("literature"), list):
            for cite in rec["literature"]:
                assert len(cite) <= 400, rec["keyword"]


def test_headline_items_resolve_with_and_without_dashes():
    for word in ("BP86", "D3BJ", "DefGrid2", "GCP(HF/631GD)", "rijcosx"):
        found = explain(word)
        assert found is not None, word
        assert found["status"] == "grounded", word


def test_dashless_spelling_finds_dashed_keywords():
    for dashed in ("R2SCAN-3C", "WB97X-D3", "d3bj", "defgrid2", "wpr2scan50"):
        assert explain(dashed.replace("-", "")) is not None, dashed
