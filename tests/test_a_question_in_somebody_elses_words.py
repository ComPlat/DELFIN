"""A retrieval fixture whose questions do not quote the passages.

Input: an index and a callable that answers a prompt. Output: a fixture of
(question, source document), with every question that reuses its passage's
language refused.

Why it is needed. The corpus-derived set in retrieval_eval measures
whether a document can be reached by language taken from it, which is what
caught a chunker indexing a quarter of each paper. On that set adding
character n-grams or LSA moved recall@1 from 0.975 to 0.959 and 0.943 --
the set answering the question it was built for, not evidence about
meaning-based retrieval. Judging that needs questions nobody copied from
the text.

Why the filter is mechanical. A model told not to reuse the passage's
words reuses them anyway often enough that a trusted fixture would quietly
become a second copy of the verbatim set, at the cost of a model call
each. The overlap is measured with the same token function the memory
store dedups by.

Model-free by construction: `ask` is the caller's. The passages are the
user's literature, and sending them anywhere is a decision for whoever
runs this, with the backend in front of them -- not a default buried here.
"""

from __future__ import annotations

import json

import pytest

from delfin.doc_server import paraphrase_eval as P

_PASSAGE = (
    "Broken-symmetry DFT was used to describe the antiferromagnetic coupling "
    "between the two iron centres, with the B3LYP functional and a def2-TZVP "
    "basis set on all atoms. The coupling constant was extracted from the "
    "energy difference between the high-spin and broken-symmetry solutions, "
    "following the standard spin-projection treatment. Geometries were "
    "optimised without symmetry constraints and confirmed as minima by a "
    "frequency calculation; solvent effects were included through an "
    "implicit continuum model with the dielectric constant of acetonitrile."
)
assert len(_PASSAGE) >= 400, "the builder skips sections shorter than that"


def _index(**docs):
    return {"built_at": "2026-09-25T00:00:00Z", "documents": {
        did: {"title": did, "sections": {
            f"s{i}": {"title": f"Section {i}", "text": t}
            for i, t in enumerate(texts)}}
        for did, texts in docs.items()}}


# -- what may enter the fixture -------------------------------------------

def test_a_question_that_quotes_the_passage_is_refused():
    ok, why = P.accept(
        "Broken-symmetry DFT was used to describe the antiferromagnetic "
        "coupling between the two iron centres", _PASSAGE)
    assert not ok and "quotes" in why


def test_a_question_in_the_researchers_words_is_accepted():
    ok, why = P.accept(
        "how do you model two metal centres whose spins cancel each other",
        _PASSAGE)
    assert ok, why


@pytest.mark.parametrize("bad, expect", [
    ("", "empty"),
    ("what method", "too short"),
    (" ".join(["word"] * 40), "too long"),
    ("a reasonable question about magnetism\nand a second line", "one line"),
])
def test_a_malformed_question_is_refused_with_its_reason(bad, expect):
    ok, why = P.accept(bad, _PASSAGE)
    assert not ok and expect in why


def test_the_threshold_is_strict_on_purpose():
    """A fixture a lexical scorer can win by matching words it was shown
    measures nothing the verbatim set does not already measure."""
    assert P._MAX_OVERLAP <= 0.70


def test_containment_and_not_jaccard(monkeypatch):
    """Measured, and the reason this was changed: Jaccard is symmetric, so
    a verbatim fragment quoted out of a LONG passage scores 0.19 -- under
    any threshold strict enough to be useful -- while the same fragment
    out of a short one scores 0.57. Containment reads 1.00 for both, and
    0.18 / 0.27 for a real question. The filter has to hold where passages
    are long, which is every real one.
    """
    quote = ("Broken-symmetry DFT was used to describe the antiferromagnetic "
             "coupling between the two iron centres")
    short = quote + " with B3LYP and def2-TZVP."
    ok_short, _ = P.accept(quote, short)
    ok_long, _ = P.accept(quote, _PASSAGE)
    assert not ok_short and not ok_long, "a quotation must fail either way"


# -- building it -----------------------------------------------------------

def test_every_section_contributes_at_most_one_question():
    idx = _index(paper=[_PASSAGE, _PASSAGE.replace("iron", "cobalt")])
    asked = []

    def ask(prompt):
        asked.append(prompt)
        return "how do you model two metal centres whose spins cancel"

    out = P.build(idx, ask, per_document=6)
    assert len(out["cases"]) == 2
    assert len(asked) == 2
    assert all(c["doc_id"] == "paper" for c in out["cases"])


def test_a_refusal_is_recorded_and_skipped_not_retried():
    idx = _index(paper=[_PASSAGE])

    def ask(prompt):
        return "Broken-symmetry DFT was used to describe the antiferromagnetic"

    out = P.build(idx, ask)
    assert out["cases"] == []
    assert out["rejected"] and "quotes" in out["rejected"][0]["reason"]
    # The rejected text is kept: it is the only way to tell "the model
    # quotes" from "this passage holds no question".
    assert out["rejected"][0]["question"]


def test_a_model_that_raises_does_not_stop_the_build():
    idx = _index(a=[_PASSAGE], b=[_PASSAGE])
    calls = {"n": 0}

    def ask(prompt):
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("backend down")
        return "how do you model two metal centres whose spins cancel"

    out = P.build(idx, ask)
    assert len(out["cases"]) == 1
    assert any("ask failed" in r["reason"] for r in out["rejected"])


def test_short_sections_are_not_asked_about():
    """A 40-character section has no question in it, and asking costs a
    model call to find that out."""
    idx = _index(paper=["too short to ask about"])
    out = P.build(idx, lambda p: pytest.fail("should not have asked"))
    assert out["cases"] == []


def test_the_build_is_reproducible():
    idx = _index(paper=[_PASSAGE + f" run {i}" for i in range(8)])

    def ask(prompt):
        return "how do you model two metal centres whose spins cancel"

    a = P.build(idx, ask, per_document=3, seed=7)
    b = P.build(idx, ask, per_document=3, seed=7)
    assert [c["section_id"] for c in a["cases"]] == \
           [c["section_id"] for c in b["cases"]]


# -- using it --------------------------------------------------------------

def test_a_fixture_for_other_documents_is_refused(tmp_path):
    idx = _index(paper=[_PASSAGE])
    fixture = {"document_ids": ["something_else"],
               "cases": [{"question": "a question here now", "doc_id": "x"}]}
    assert P.cases(fixture, idx) == [], (
        "a set built for other documents produces a number about nothing")
    assert len(P.cases(fixture)) == 1, "without an index it is the caller's"


def test_it_survives_a_round_trip(tmp_path):
    path = tmp_path / "fixture.json"
    P.save({"version": 1, "cases": [{"question": "q", "doc_id": "d"}]}, path)
    assert json.loads(path.read_text())["cases"][0]["doc_id"] == "d"
    assert P.load(path)["cases"][0]["question"] == "q"


def test_a_missing_fixture_is_empty_not_an_error(tmp_path):
    assert P.load(tmp_path / "nothing.json") == {}


def test_the_summary_says_why_questions_were_refused():
    fixture = {"cases": [{"question": "q", "doc_id": "d"}],
               "rejected": [{"reason": "quotes the passage (overlap 0.40)"},
                            {"reason": "too short (2 words)"},
                            {"reason": "quotes the passage (overlap 0.55)"}]}
    s = P.summary(fixture)
    assert s == {"accepted": 1, "rejected": 3,
                 "reasons": {"quotes the passage": 2, "too short": 1}}


# -- the boundary ----------------------------------------------------------

def test_the_module_chooses_no_backend_and_reaches_no_network():
    """The passages are the user's literature. Which service sees them is
    a decision for whoever runs this, not a default in a library."""
    import inspect

    src = inspect.getsource(P)
    for forbidden in ("requests", "urllib", "http://", "https://",
                      "openai", "anthropic", "api_key", "load_settings"):
        assert forbidden not in src, forbidden
