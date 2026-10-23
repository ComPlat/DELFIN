"""A document is indexed whole, and retrieval is measured before it changes.

Input: PDF pages. Output: sections covering all of the extracted text.

The defect: the chunker kept the first 12000 characters of a document with
no numbered headings and dropped the rest without saying so, truncated any
heading section over 12000 with a marker, and never emitted the text
BEFORE the first heading at all -- in a paper, the title, abstract and
introduction. Measured over the corpus, as the share of extracted
characters that reached the index:

    1-s2.0-S0021979724009044     40.0%
    orca_manual_6_1_1            83.4%    (its numbered headings are found)
    partial-to-total             24.4%
    s13321-025-01008-1           48.6%

No retrieval method finds text that was never indexed, which is why this
comes before any change to the ranking.
"""

from __future__ import annotations

from delfin.doc_server import retrieval_eval as E
from delfin.doc_server.indexer import _chunk_pdf_into_sections, _windows


def _pages(*texts):
    return [{"page": i + 1, "text": t} for i, t in enumerate(texts)]


def _all_text(sections):
    return " ".join(s["text"] for s in sections)


# -- nothing is dropped ----------------------------------------------------

def test_a_document_without_headings_is_kept_whole():
    body = "\n\n".join(f"Paragraph {i} with enough words to matter." * 8
                       for i in range(120))
    sections = _chunk_pdf_into_sections(_pages(body))
    assert len(sections) > 1, "one section means it was truncated again"
    assert "Paragraph 119" in _all_text(sections), "the end was dropped"
    assert "Paragraph 0" in _all_text(sections)


def test_the_text_before_the_first_heading_is_kept():
    """In a paper this is the title, the abstract and the introduction."""
    front = "A distinctive opening sentence nobody else writes.\n\n"
    sections = _chunk_pdf_into_sections(
        _pages(front + "1 Methods\n\nWhat was done here."))
    assert "distinctive opening sentence" in _all_text(sections)
    assert any(s["section_id"].startswith("front_") for s in sections)


def test_a_long_heading_section_is_split_not_cut():
    long_body = "\n\n".join(f"Sentence {i} about a measured quantity." * 12
                            for i in range(200))
    sections = _chunk_pdf_into_sections(_pages("1 Results\n\n" + long_body))
    assert "[... truncated]" not in _all_text(sections)
    assert "Sentence 199" in _all_text(sections)
    assert len([s for s in sections if s["title"].startswith("1 Results")]) > 1


def test_no_pages_is_still_a_failure_not_an_empty_document():
    """Unchanged, and the reason is in the indexer: an empty section made
    the manual read as indexed and answer nothing, forever."""
    assert _chunk_pdf_into_sections([]) == []


# -- the windows themselves -----------------------------------------------

def test_windows_keep_every_paragraph():
    text = "\n\n".join(f"para {i}" for i in range(50))
    parts = _windows(text, target=60)
    assert len(parts) > 1
    for i in range(50):
        assert f"para {i}" in " ".join(parts)


def test_a_paragraph_longer_than_the_target_is_not_cut():
    """An over-long part costs ranking; a cut sentence costs meaning."""
    long_para = "word " * 2000
    parts = _windows(long_para.strip(), target=100)
    assert len(parts) == 1
    assert parts[0].count("word") == 2000


def test_empty_text_makes_no_parts():
    assert _windows("") == []
    assert _windows("   \n\n  ") == []


# -- the measurement -------------------------------------------------------

def _index(**docs):
    return {"documents": {
        did: {"title": did, "sections": {
            f"s{i}": {"title": t, "text": x}
            for i, (t, x) in enumerate(secs)}}
        for did, secs in docs.items()}}


def test_a_query_is_made_from_a_title_and_from_the_body():
    idx = _index(paper=[("Electrocatalytic reduction of carbon dioxide",
                         " ".join(f"word{i}" for i in range(80)))])
    qs = E.queries(idx)
    assert any(q.startswith("Electrocatalytic") for q, _ in qs)
    assert any("word4" in q for q, _ in qs), "nothing came from the body"


def test_the_body_query_comes_from_the_middle():
    """The start is what a truncating indexer keeps, so asking about the
    start cannot see what such an indexer lost."""
    text = " ".join(["START"] * 40 + ["MIDDLE"] * 40 + ["END"] * 40)
    phrase = E._body_phrase(text)
    assert "MIDDLE" in phrase and "START" not in phrase


def test_a_heading_that_names_nothing_makes_no_query():
    for junk in ("Introduction", "Part 3", "4.2", "Conclusions", "References"):
        assert not E._usable(junk), junk
    assert E._usable("Broken-symmetry DFT for antiferromagnetic coupling")


def test_one_large_document_cannot_decide_the_score():
    """The ORCA manual has 2373 sections against a paper's dozen."""
    idx = _index(
        big=[(f"A distinctive heading number {i} here", "body " * 60)
             for i in range(200)],
        small=[("Another distinctive heading entirely", "body " * 60)])
    counts = {}
    for _, did in E.queries(idx, per_document=12):
        counts[did] = counts.get(did, 0) + 1
    assert counts["big"] <= 24, counts


def test_a_document_that_cannot_be_found_scores_zero_and_is_named():
    class _Engine:
        def search(self, q, max_results=10):
            return {"results": [{"doc_id": "other"}]}

    out = E.score(_Engine(), [("a query", "wanted")], depth=10)
    assert out["recall@1"] == 0.0 and out["mrr"] == 0.0
    assert out["unreachable"] == ["a query"], "a miss has to be readable"


def test_an_engine_that_raises_is_a_miss_not_a_crash():
    class _Engine:
        def search(self, q, max_results=10):
            raise RuntimeError("no index")

    assert E.score(_Engine(), [("q", "d")])["recall@1"] == 0.0


def test_no_queries_is_reported_rather_than_divided_by():
    assert E.score(object(), [])["queries"] == 0
