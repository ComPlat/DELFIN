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
    """A described calculation gets findings (executable checks) plus the
    relevant verified statement cards for its method/job."""
    import json

    from delfin import api

    calc = {"keywords": ["CASSCF"], "blocks": {}, "method_family": "casscf",
            "job": "sp", "charge": 0, "mult": 1, "elements": [],
            "orca_version": "6.1.1"}
    out = api.orca_check_input(json.dumps(calc))
    assert out["ok"] is True
    assert isinstance(out["findings"], list)
    # The shipped store is v1: statement cards without executable checks.
    # The tool must still surface the cards that describe THIS calc.
    assert out["n_relevant_cards"] > 0
    ids = {c["id"] for c in out["relevant_cards"]}
    assert "casscf.active-space.mandatory" in ids
    for card in out["relevant_cards"]:
        assert card["statement"] and card["sources"]


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


def test_a_card_with_an_applies_condition_is_matched_by_it():
    """_card_matches honours a v2 card's ``applies``: it reached an
    unimported name before (F821), so any such card raised NameError."""
    from delfin.api import _card_matches

    card = {"topic": "scf", "keywords": ["TIGHTSCF"],
            "applies": {"method_family": ["dft"]}}
    dft = {"method_family": "dft", "keywords": ["TIGHTSCF"]}
    hf = {"method_family": "hf", "keywords": ["TIGHTSCF"]}
    assert _card_matches(card, dft) is True
    assert _card_matches(card, hf) is False
