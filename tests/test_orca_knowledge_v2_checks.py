"""Executable knowledge cards (v2): applies, checks, relations.

Covers every condition operator, nesting, evaluation to findings, the graph
problems, and that every shipped v1 card stays valid under v2 validation.
"""

from __future__ import annotations

import json


from delfin.orca_knowledge import CARDS_DIR, known_keyword, validate
from delfin.orca_knowledge.checks import (applies, check_graph, check_problems,
                                          condition_problems,
                                          evaluate_check, findings,
                                          relations_problems)


def base_card(**extra) -> dict:
    """A minimal valid v1 card, with v2 fields added via ``extra``."""
    card = {
        "id": "test.v2.card",
        "topic": "scf",
        "kind": "check",
        "statement": "A test card used by the v2 checks tests.",
        "orca": ">=6.0",
        "keywords": ["TightSCF"],
        "severity": "warning",
        "sources": [{"type": "manual", "version": "6.1",
                     "page": "contents/essentialelements/scf.html"}],
        "status": "proposed",
        "author": "i04",
    }
    card.update(extra)
    return card


def calc(**kw) -> dict:
    base = {"keywords": ["RKS", "TIGHTSCF"],
            "blocks": {"scf": {"maxiter": 200, "convergence": "tight"}},
            "method_family": "dft", "job": "opt",
            "charge": 0, "mult": 1,
            "elements": ["Fe", "O"],
            "orca_version": "6.1.1"}
    base.update(kw)
    return base


# ------------------------------------------------------------- v1 stays valid

def test_all_shipped_cards_validate():
    problems = []
    for path in sorted(CARDS_DIR.glob("*.json")):
        rows = json.loads(path.read_text())
        for card in rows:
            problems += [f"{path.name} {card.get('id')}: {p}"
                         for p in validate(card)]
    assert problems == []


# ------------------------------------------------------------ condition forms

def test_condition_problems_unknown_key():
    assert condition_problems({"keywordz": "TightSCF"}, known_keyword)
    assert condition_problems({"mult": ">>1"}, known_keyword)
    assert condition_problems({"mult": "1"}, known_keyword)
    assert condition_problems({"orca": ">six"}, known_keyword)
    assert condition_problems({"keyword": "NOTAKEYWORD"}, known_keyword)
    assert condition_problems({"all": []}, known_keyword)
    assert condition_problems({"any": "x"}, known_keyword)
    assert condition_problems({"not": {}}, known_keyword)
    assert condition_problems({"block_option": ["scf", "MAXITER"]}, known_keyword)
    assert condition_problems({"method_family": []}, known_keyword)
    assert condition_problems("dft", known_keyword)


def test_condition_problems_accept():
    assert condition_problems({"keyword": "TightSCF"}, known_keyword) == []
    assert condition_problems({"all": [{"keyword": "TightSCF"},
                                       {"not": {"mult": "==1"}}]},
                              known_keyword) == []
    assert condition_problems({"orca": ">=6.0,<7"}, known_keyword) == []
    assert condition_problems({"block_option": ["scf", "maxiter"]}) == []
    assert condition_problems({"mult": ">1"}) == []
    assert condition_problems({"charge": "<=0"}) == []
    assert condition_problems({"method_family": ["dft", "casscf"]}) == []
    assert condition_problems({"job": ["opt", "ts"]}) == []
    assert condition_problems({"elements_in": ["3d", "4f"]}) == []
    assert condition_problems({"elements_in": ["Fe"]}) == []


def test_applies_each_operator():
    assert applies({"applies": {"keyword": "TightSCF"}}, calc())
    assert not applies({"applies": {"keyword": "SLOWCONV"}}, calc())
    # case-blind
    assert applies({"applies": {"keyword": "tightscf"}}, calc())
    assert applies({"applies": {"block_option": ["scf", "maxiter"]}}, calc())
    assert not applies({"applies": {"block_option": ["scf", "nel"]}}, calc())
    assert applies({"applies": {"method_family": ["casscf", "dft"]}}, calc())
    assert not applies({"applies": {"method_family": ["casscf"]}}, calc())
    assert applies({"applies": {"job": ["ts", "opt"]}}, calc())
    assert not applies({"applies": {"job": ["freq"]}}, calc())
    assert not applies({"applies": {"mult": ">1"}}, calc())
    assert applies({"applies": {"mult": "<=1"}}, calc())
    assert applies({"applies": {"mult": "==1"}}, calc())
    assert applies({"applies": {"charge": ">=0"}}, calc())
    assert not applies({"applies": {"charge": ">0"}}, calc())
    assert applies({"applies": {"elements_in": ["3d"]}}, calc())
    assert applies({"applies": {"elements_in": ["Fe"]}}, calc())
    assert applies({"applies": {"elements_in": ["O"]}}, calc())
    assert not applies({"applies": {"elements_in": ["4f"]}}, calc())
    assert applies({"applies": {"orca": ">=6.0"}}, calc())
    assert applies({"applies": {"orca": ">=6.1.1"}}, calc())
    assert not applies({"applies": {"orca": ">=6.2"}}, calc())
    assert not applies({"applies": {"orca": ">=6.0,<6.1"}}, calc())


def test_applies_nesting():
    assert applies({"applies": {"all": [{"keyword": "TightSCF"},
                                        {"method_family": ["dft"]}]}}, calc())
    assert not applies({"applies": {"all": [{"keyword": "TightSCF"},
                                            {"method_family": ["casscf"]}]}}, calc())
    assert applies({"applies": {"any": [{"method_family": ["casscf"]},
                                        {"keyword": "TightSCF"}]}}, calc())
    assert applies({"applies": {"not": {"keyword": "SLOWCONV"}}}, calc())
    assert not applies({"applies": {"not": {"method_family": ["dft"]}}}, calc())
    deep = {"all": [{"any": [{"keyword": "NOKEYWORD"}, {"mult": "<=2"}]},
                    {"not": {"any": [{"charge": ">4"}, {"elements_in": ["5f"]}]}}]}
    assert applies({"applies": deep}, calc())


def test_applies_without_condition_is_v1():
    assert applies({"id": "x"}, calc()) is True
    assert applies({}, calc()) is True


def test_applies_versionless_calc():
    # unknown version: open lower bounds cannot exclude it, upper ones can
    assert applies({"applies": {"orca": ">=6.0"}}, calc(orca_version=None))
    assert not applies({"applies": {"orca": "<7"}}, calc(orca_version=None))


# ------------------------------------------------------------------- checks

def test_check_problems():
    assert check_problems({"requires_keyword": ["TightSCF"]}, known_keyword) == []
    assert check_problems({"forbids_keyword": ["SLOWCONV"]}, known_keyword) == []
    assert check_problems({"requires_block_option": [["casscf", "nel"]]},
                          known_keyword) == []
    assert check_problems({"requires_block_option": ["casscf.nel"]},
                          known_keyword) == []
    assert check_problems({"forbids_combination": ["SLOWCONV", "scf.maxiter"]},
                          known_keyword) == []
    assert check_problems({"requires_keyword": ["NOTAKEYWORD"]}, known_keyword)
    assert check_problems({"forbids_keyword": ["NOTAKEYWORD"]}, known_keyword)
    assert check_problems({"requires_keyword": []}, known_keyword)
    assert check_problems({"bogus_key": ["TightSCF"]}, known_keyword)
    assert check_problems({"requires_block_option": [["scf"]]}, known_keyword)
    assert check_problems({"forbids_combination": ["TightSCF"]}, known_keyword)
    assert check_problems({"requires_keyword": "TightSCF"}, known_keyword)


def test_evaluate_check_findings():
    card = base_card(check={"requires_keyword": ["SlowConv"]},
                     severity="warning")
    assert applies(card, calc())  # no applies condition
    f = evaluate_check(card, calc())
    assert len(f) == 1
    assert f[0]["id"] == card["id"]
    assert f[0]["severity"] == "warning"
    assert f[0]["statement"] == card["statement"]
    assert f[0]["sources"] == card["sources"]
    assert "SlowConv" in f[0]["detail"]

    card_ok = base_card(check={"requires_keyword": ["TightSCF"]})
    assert evaluate_check(card_ok, calc()) == []

    forbid = base_card(check={"forbids_keyword": ["TightSCF"]})
    f = evaluate_check(forbid, calc())
    assert len(f) == 1 and "TightSCF" in f[0]["detail"]

    req_block = base_card(check={"requires_block_option": ["casscf", "nel"]})
    assert len(evaluate_check(req_block, calc())) == 1
    ok_block = base_card(check={"requires_block_option": ["scf", "maxiter"]})
    assert evaluate_check(ok_block, calc()) == []

    combo = base_card(check={"forbids_combination": ["TightSCF", "scf.maxiter"]})
    assert len(evaluate_check(combo, calc())) == 1
    combo_ok = base_card(check={"forbids_combination": ["SlowConv", "scf.nel"]})
    assert evaluate_check(combo_ok, calc()) == []


def test_evaluate_check_respects_applies():
    card = base_card(applies={"method_family": ["casscf"]},
                     check={"requires_keyword": ["SlowConv"]})
    assert evaluate_check(card, calc()) == []  # dft calc: card does not apply
    casscf = calc(method_family="casscf", keywords=["CASSCF"])
    assert len(evaluate_check(card, casscf)) == 1


def test_findings_default_store_and_malformed_card():
    good = base_card(check={"requires_keyword": ["NOTAKEYWORD"]})
    # a malformed card must surface, not vanish silently
    malformed = base_card(id="test.v2.malformed",
                          applies={"what": "ever"},
                          check={"requires_keyword": ["TightSCF"]})
    out = findings(calc(), cards=[good, malformed])
    assert all(f["id"] in ("test.v2.card", "test.v2.malformed") for f in out)
    assert any(f["id"] == "test.v2.malformed" and "malformed" in f["detail"]
               for f in out)
    # the well-formed card is still evaluated normally (its required keyword
    # is genuinely absent, so it reports -- with the keyword as written)
    assert any(f["id"] == "test.v2.card" and "NOTAKEYWORD" in f["detail"]
               for f in out)


def test_findings_full_shipped_store_runs():
    # every shipped card must be evaluatable against a calc without raising
    out = findings(calc())
    assert isinstance(out, list)


# ----------------------------------------------------------------- relations

def test_relations_problems():
    assert relations_problems({"requires": ["scf.slowconv-tm-default"]}) == []
    assert relations_problems({"see_also": ["a.b-c_d", "x.y"]}) == []
    assert relations_problems({"bogus": ["x"]})
    assert relations_problems({"requires": []})
    assert relations_problems({"requires": ["BAD ID"]})
    assert relations_problems({"requires": "scf.x"})


# ------------------------------------------------------------------- graph

def _graph_cards(*cards: dict) -> list[dict]:
    return cards


def test_check_graph_unknown_id_and_mutuality():
    a = base_card(id="test.v2.a", relations={"requires": ["test.v2.missing"]})
    problems = check_graph(_graph_cards(a))
    assert any("test.v2.a" in p and "missing" in p for p in problems)

    b = base_card(id="test.v2.b", relations={"conflicts_with": ["test.v2.c"]})
    c = base_card(id="test.v2.c")
    problems = check_graph(_graph_cards(b, c))
    assert any("not mutual" in p for p in problems)

    c2 = base_card(id="test.v2.c", relations={"conflicts_with": ["test.v2.b"]})
    assert check_graph(_graph_cards(b, c2)) == []


def test_check_graph_requires_cycle():
    a = base_card(id="test.v2.a", relations={"requires": ["test.v2.b"]})
    b = base_card(id="test.v2.b", relations={"requires": ["test.v2.a"]})
    problems = check_graph(_graph_cards(a, b))
    assert any("cycle" in p for p in problems)

    chain = [base_card(id="test.v2.a", relations={"requires": ["test.v2.b"]}),
             base_card(id="test.v2.b", relations={"requires": ["test.v2.c"]}),
             base_card(id="test.v2.c")]
    assert check_graph(chain) == []


def test_check_graph_require_forbid_overlap():
    req = base_card(id="test.v2.req",
                    applies={"method_family": ["dft"]},
                    check={"requires_keyword": ["TightSCF"]},
                    status="verified", verified_by="reviewer")
    forb = base_card(id="test.v2.forb",
                     applies={"job": ["opt"]},
                     check={"forbids_keyword": ["TightSCF"]},
                     status="verified", verified_by="reviewer")
    problems = check_graph(_graph_cards(req, forb))
    assert any("tightscf" in p.lower() and "test.v2.req" in p and "test.v2.forb" in p
               for p in problems)


def test_check_graph_require_forbid_disjoint_conditions_ok():
    # method families provably disjoint: no problem is reported
    req = base_card(id="test.v2.req", applies={"method_family": ["dft"]},
                    check={"requires_keyword": ["TightSCF"]},
                    status="verified", verified_by="reviewer")
    forb = base_card(id="test.v2.forb", applies={"method_family": ["casscf"]},
                     check={"forbids_keyword": ["TightSCF"]},
                     status="verified", verified_by="reviewer")
    assert check_graph(_graph_cards(req, forb)) == []


def test_check_graph_require_forbid_mult_disjoint_ok():
    req = base_card(id="test.v2.req", applies={"mult": "==1"},
                    check={"requires_keyword": ["TightSCF"]},
                    status="verified", verified_by="reviewer")
    forb = base_card(id="test.v2.forb", applies={"mult": ">1"},
                     check={"forbids_keyword": ["TightSCF"]},
                     status="verified", verified_by="reviewer")
    assert check_graph(_graph_cards(req, forb)) == []


def test_check_graph_not_keyword_exclusion_detected():
    # one card requires TightSCF; another applies only when TightSCF is NOT
    # set: provably disjoint, so no problem
    req = base_card(id="test.v2.req", applies={"keyword": "TightSCF"},
                    check={"requires_keyword": ["SlowConv"]},
                    status="verified", verified_by="reviewer")
    forb = base_card(id="test.v2.forb",
                     applies={"not": {"keyword": "TightSCF"}},
                     check={"forbids_keyword": ["SlowConv"]},
                     status="verified", verified_by="reviewer")
    assert check_graph(_graph_cards(req, forb)) == []


def test_check_graph_unverified_pair_not_flagged():
    req = base_card(id="test.v2.req", applies={"method_family": ["dft"]},
                    check={"requires_keyword": ["TightSCF"]},
                    status="proposed")
    forb = base_card(id="test.v2.forb", applies={"job": ["opt"]},
                     check={"forbids_keyword": ["TightSCF"]},
                     status="proposed")
    assert check_graph(_graph_cards(req, forb)) == []


def test_check_graph_full_shipped_store():
    assert check_graph() == []


# ------------------------------------------------------------------ validate

def test_validate_refuses_bad_v2_fields():
    card = base_card(applies={"keywordz": "TightSCF"})
    assert any("unknown condition key" in p for p in validate(card))
    card = base_card(applies={"keyword": "NOTAKEYWORD"})
    assert any("NOTAKEYWORD" in p for p in validate(card))
    card = base_card(check={"forbids_keyword": ["NOTAKEYWORD"]})
    assert any("NOTAKEYWORD" in p for p in validate(card))
    card = base_card(check={"nonsense": []})
    assert any("unknown check key" in p for p in validate(card))
    card = base_card(relations={"requires": ["nonexistent card"]})
    assert any("card ids" in p for p in validate(card))


def test_validate_accepts_full_v2_card():
    card = base_card(
        applies={"all": [{"method_family": ["casscf", "nevpt2"]},
                         {"mult": ">1"},
                         {"elements_in": ["3d", "4d"]}]},
        check={"requires_keyword": ["SlowConv"],
               "forbids_combination": ["TightSCF", "scf.maxiter"],
               "requires_block_option": ["casscf", "nel"]},
        relations={"requires": ["scf.slowconv-tm-default"],
                   "see_also": ["scf.veryslowconv"]},
    )
    assert validate(card) == []


def test_v2_fields_are_optional_in_validate():
    assert validate(base_card()) == []
