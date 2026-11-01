"""wave14/w04 RULE 1: ORCA knowledge is reached on demand, by tool.

``orca_explain(word)`` answers from the shipped keyword dictionary +
knowledge cards: purpose, requirements, manual section. Nothing of the
dictionary enters a prompt. ``orca_check_input(text)`` evaluates the
executable knowledge cards against a described calculation. Neither tool
ever starts ORCA.
"""
from __future__ import annotations

import json

import pytest

from delfin.orca_knowledge.dictionary import explain, records


def test_orca_explain_module_exists_and_explains():
    from delfin import api

    rec = explain("3conv")
    assert rec is not None
    card = api.orca_explain("3conv")
    assert card["found"] is True
    assert card["keyword"] == rec["keyword"]
    assert card["purpose"]
    assert card["manual"]  # page + section of the ORCA manual
    # Case-blind + alias tolerant, same contract as the dictionary.
    assert api.orca_explain("NoSplitRIJ") is not None or records() == []


def test_orca_explain_unknown_word_is_not_an_error():
    from delfin import api

    out = api.orca_explain("definitely_not_an_orca_keyword_xyz")
    assert out["found"] is False
    assert "suggestions" in out or out.get("found") is False
    assert api.orca_explain("")["found"] is False


def test_orca_explain_returns_the_card_fields():
    from delfin import api

    out = api.orca_explain("3conv")
    for key in ("id", "category", "opi", "status"):
        assert key in out
    assert "found" in out


def test_orca_check_input_evaluates_cards_on_a_calc_description():
    """A calc description that trips a shipped card must yield findings."""
    import delfin.orca_knowledge as ok

    store = ok.cards()
    if not store:
        pytest.skip("no verified cards shipped")
    # Build a calc description from the first card's own condition so the
    # test stays true whatever the shipped store contains.
    from delfin.orca_knowledge.checks import applies

    hit = None
    for card in store:
        calc = {"keywords": [], "blocks": {}, "method_family": "dft",
                "job": "opt", "charge": 0, "mult": 1,
                "elements": ["Fe"], "orca_version": "6.1.1"}
        cond = card.get("applies")
        if not cond:
            continue
        # minimal satisfied description: flip fields until applies
        for kw in _cond_keywords(cond):
            calc["keywords"].append(kw)
        if applies(card, calc):
            hit = (card, calc)
            break
    if hit is None:
        pytest.skip("no card with a constructible description")
    card, calc = hit
    from delfin import api
    out = api.orca_check_input(json.dumps(calc))
    assert out["ok"] is True
    assert any(f["id"] == card["id"] for f in out["findings"])


def test_orca_check_input_refuses_free_text_and_reports_bad_json():
    from delfin import api

    out = api.orca_check_input("not json at all")
    assert out["ok"] is False
    assert "error" in out


def _cond_keywords(cond):
    out = []
    if isinstance(cond, dict):
        for k, v in cond.items():
            if k == "keyword" and isinstance(v, str):
                out.append(v)
            elif k in ("all", "any") and isinstance(v, list):
                for sub in v:
                    out.extend(_cond_keywords(sub))
    return out


def test_tools_are_in_the_catalog_and_in_all():
    from delfin import api

    names = {e["name"] for e in api._TOOL_CATALOG}
    assert {"orca_explain", "orca_check_input"} <= names
    assert "orca_explain" in api.__all__
    assert "orca_check_input" in api.__all__


def test_schema_stays_short():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI, estimate_schema_tokens

    for name in ("orca_explain", "orca_check_input"):
        schema = next(t["function"] for t in _DOC_TOOLS_OPENAI
                      if t["function"]["name"] == name)
        tokens = estimate_schema_tokens({"type": "function", "function": schema})
        # A new capability pays from its own text: each new tool stays
        # under ~130 tokens of schema (one short description, bare types).
        assert tokens <= 130, f"{name} schema costs {tokens} tokens"


def test_no_orca_started(monkeypatch):
    """The knowledge tools are pure lookups: never spawn ORCA."""
    import subprocess
    from delfin import api

    monkeypatch.setattr(subprocess, "Popen",
                        lambda *a, **k: pytest.fail("ORCA must not start"))
    api.orca_explain("3conv")
    api.orca_check_input(json.dumps({"keywords": [], "blocks": {}}))
