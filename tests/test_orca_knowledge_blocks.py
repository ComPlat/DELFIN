"""Structured block options in cards: checked names, and an input renderer.

The expert audits found block options written as prose ("%mdci TCutPNO
3.33e-7" with no "end", IntAcc put into %grid) -- advice ORCA cannot read.
A card carries its block advice structured::

    "recommend": {"blocks": {"mdci": {"TCutPNO": "3.33e-7"}}}
    "avoid": {"blocks": {"grid": {"IntAcc": 7}}}          # wrong block, refused

Every block name and option name is checked against the ``blocks`` section
of ``opi_keywords.json`` (generated from the pinned OPI's typed block
models), so a misspelt option or a block OPI does not model is refused with
its name. Options OPI does not model but ORCA's parser accepts live in
``blocks_orca_parser`` and are kept when probes measure them.
"""

from __future__ import annotations

import json

from delfin import orca_knowledge as ok


def _card(**over):
    card = {
        "id": "mdci.tcutpno.tighter",
        "topic": "dlpno",
        "kind": "recommendation",
        "statement": ("For tighter DLPNO-CCSD thresholds set TCutPNO in the %mdci "
                      "block to 3.33e-7; the default 3.33e-6 loses weak correlation."),
        "orca": ">=6.0",
        "keywords": [],
        "sources": [{"type": "manual", "version": "6.1",
                     "page": "contents/methods/mdci.html"}],
        "status": "proposed",
        "author": "agent-sb1",
    }
    card.update(over)
    return card


def test_a_well_formed_blocks_recommendation_passes():
    assert ok.validate(_card(
        recommend={"blocks": {"mdci": {"TCutPNO": "3.33e-7"}}})) == []


def test_a_prose_block_option_is_not_the_structured_form():
    # The old failure mode: prose that ORCA cannot read stays text, but the
    # structured field is what the checks read -- an option that appears only
    # as "%mdci TCutPNO ..." inside the statement is caught by the keyword
    # scan, not silently trusted. The renderer is the way blocks go into an
    # input, and it writes the "end".
    text = ok.to_input({"mdci": {"TCutPNO": "3.33e-7"}})
    assert "%mdci" in text and "TCutPNO 3.33e-7" in text
    assert text.strip().splitlines()[-1].strip().lower() == "end"


def test_to_input_renders_complete_blocks():
    text = ok.to_input({"scf": {"MaxIter": 200}, "mdci": {"TCutPNO": "3.33e-7"}})
    blocks = [line.strip().lower() for line in text.strip().splitlines()]
    assert blocks[0] == "%scf" and "maxiter 200" in blocks
    assert "%mdci" in blocks and "tcutpno 3.33e-7" in blocks
    assert blocks.count("end") == 2


def test_an_unknown_block_is_refused_with_its_name():
    problems = ok.validate(_card(recommend={"blocks": {"mdcit": {"TCutPNO": "1"}}}))
    assert any("mdcit" in p for p in problems)


def test_an_unknown_option_is_refused_with_its_name():
    problems = ok.validate(_card(
        recommend={"blocks": {"mdci": {"TCutPNNO": "3.33e-7"}}}))
    assert any("TCutPNNO" in p for p in problems)


def test_an_option_in_the_wrong_block_is_refused():
    # IntAcc is a %method option; putting it in %grid is what the audit found.
    grid = _card(recommend={"blocks": {"grid": {"IntAcc": 7}}})
    assert any("IntAcc" in p or "grid" in p for p in ok.validate(grid))


def test_avoid_blocks_are_checked_too():
    assert ok.validate(_card(avoid={"blocks": {"scf": {"MaxIter": 100}}})) == []
    problems = ok.validate(_card(avoid={"blocks": {"nosuch": {"X": 1}}}))
    assert any("nosuch" in p for p in problems)


def test_blocks_orca_parser_options_are_accepted(tmp_path, monkeypatch):
    # OPI does not model it, but ORCA's parser measured it: accepted.
    assert ok.known_block_option("mdci", "SomeProbedOption") is False
    catalog = json.loads(ok.KEYWORDS_FILE.read_text())
    catalog.setdefault("blocks_orca_parser", []).append("mdci.SomeProbedOption")
    monkeypatch.setattr(ok, "KEYWORDS_FILE", tmp_path / "catalog.json")
    (tmp_path / "catalog.json").write_text(json.dumps(catalog))
    ok._CATALOG = None  # force a re-read of the patched file
    assert ok.known_block_option("mdci", "SomeProbedOption") is True
    assert ok.known_block_option("mdci", "OtherProbedOption") is False
    assert ok.validate(_card(
        recommend={"blocks": {"mdci": {"SomeProbedOption": 1}}})) == []
    ok._CATALOG = None


def test_block_option_conditions_are_checked():
    # Conditions describe a calculation: lowercase names, as the calc dict
    # stores them.
    assert ok.validate(_card(applies={"block_option": ["mdci", "tcutpno"]})) == []
    problems = ok.validate(_card(applies={"block_option": ["mdci", "nosuchopt"]}))
    assert any("nosuchopt" in p for p in problems)
    problems = ok.validate(_card(applies={"block_option": ["nosuchblock", "x"]}))
    assert any("nosuchblock" in p for p in problems)


def test_requires_block_option_in_check_is_checked():
    card = _card(check={"requires_block_option": ["mdci", "tcutpno"]})
    assert ok.validate(card) == []
    card = _card(check={"requires_block_option": ["mdci", "tcutpnno"]})
    assert any("tcutpnno" in p for p in ok.validate(card))
    card = _card(check={"forbids_combination": ["mdci.tcutpno", "SlowConv"]})
    assert ok.validate(card) == []
    card = _card(check={"forbids_combination": ["mdci.nosuchopt", "SlowConv"]})
    assert any("nosuchopt" in p for p in ok.validate(card))
