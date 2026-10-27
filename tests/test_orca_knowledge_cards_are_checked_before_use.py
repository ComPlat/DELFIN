"""The ORCA knowledge store: cards that say how ORCA is used correctly.

A card is one statement about ORCA -- a requirement, a recommendation, an
incompatibility, an error and its cause -- with the ORCA versions it holds
for and at least one source (a manual section, a paper by DOI, an ORCA
message). Cards are data in ``delfin/orca_knowledge/cards/*.json``; the
checks of a calculation read them.

* A card that breaks the schema is refused, with every reason: unknown
  topic or kind, too short a statement, no source, a manual source without
  its page, a paper without a DOI, a version range that is no range.
* Every ORCA keyword a card names is one ORCA knows: checked against the
  keyword catalog of the pinned OPI (``opi_keywords.json``, generated from
  OPI, so the check needs no OPI installed). A misspelt or invented keyword
  cannot enter.
* A card is "verified" only with a verifier who is not its author; only
  verified cards are used. Ids are unique.
"""

from __future__ import annotations

import copy
import json

import pytest

from delfin import orca_knowledge as ok


def _card(**over):
    card = {
        "id": "scf.slowconv.open-shell-tm",
        "topic": "scf",
        "kind": "recommendation",
        "statement": ("For open-shell transition-metal complexes SlowConv damps the "
                      "early SCF iterations and often makes a diverging SCF converge."),
        "orca": ">=6.0",
        "keywords": ["SLOWCONV"],
        "sources": [{"type": "manual", "version": "6.1",
                     "page": "contents/essentialelements/scf.html",
                     "section": "SCF convergence"}],
        "status": "proposed",
        "author": "agent-07",
    }
    card.update(over)
    return card


def test_a_well_formed_card_passes():
    assert ok.validate(_card()) == []


@pytest.mark.parametrize("change, words", [
    ({"topic": "astrology"}, "topic"),
    ({"kind": "rumour"}, "kind"),
    ({"statement": "short"}, "statement"),
    ({"sources": []}, "source"),
    ({"sources": [{"type": "manual", "version": "6.1"}]}, "page"),
    ({"sources": [{"type": "paper", "title": "x"}]}, "doi"),
    ({"sources": [{"type": "paper", "doi": "not-a-doi"}]}, "doi"),
    ({"orca": "newest"}, "orca"),
    ({"id": "Has Spaces"}, "id"),
    ({"status": "verified"}, "verified_by"),
    ({"status": "verified", "verified_by": "agent-07"}, "verified_by"),
])
def test_what_breaks_the_schema_is_named(change, words):
    problems = ok.validate(_card(**change))
    assert problems and any(words in p for p in problems), problems


def test_a_keyword_orca_does_not_know_is_refused():
    problems = ok.validate(_card(keywords=["SLOWCONV", "SUPERFASTCONV"]))
    assert any("SUPERFASTCONV" in p for p in problems)


def test_keywords_are_matched_as_orca_reads_them():
    assert ok.known_keyword("SlowConv") and ok.known_keyword("def2-TZVP")
    assert ok.known_keyword("RIJCOSX") and ok.known_keyword("DLPNO-CCSD(T)")
    assert not ok.known_keyword("B3LYPP")


def test_the_keyword_catalog_comes_from_the_pinned_opi():
    from delfin import doctor
    catalog = json.loads(ok.KEYWORDS_FILE.read_text())
    assert catalog["opi"] == doctor.OPI_VERSION
    assert len(catalog["keywords"]) > 1000


def test_the_catalog_matches_opi_where_it_is_installed():
    pytest.importorskip("opi")
    from importlib import metadata
    from delfin import doctor
    if metadata.version("orca-pi") != doctor.OPI_VERSION:
        pytest.skip("another OPI version is installed")
    assert ok.keywords_from_opi() == json.loads(ok.KEYWORDS_FILE.read_text())["keywords"]


def test_the_shipped_cards_are_valid_unique_and_only_verified_ones_are_used(tmp_path, monkeypatch):
    problems = ok.check_store()
    assert problems == [], problems
    good = _card(id="a.one", status="verified", verified_by="agent-11")
    draft = _card(id="a.two")
    (tmp_path / "x.json").write_text(json.dumps([good, draft]))
    monkeypatch.setattr(ok, "CARDS_DIR", tmp_path)
    assert [c["id"] for c in ok.cards()] == ["a.one"]
    assert {c["id"] for c in ok.cards(include_unverified=True)} == {"a.one", "a.two"}
    (tmp_path / "y.json").write_text(json.dumps([copy.deepcopy(good)]))
    assert any("duplicate" in p for p in ok.check_store())
