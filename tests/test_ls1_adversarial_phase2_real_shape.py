"""LS1 phase-2 adversarial tests, phase-2 red findings.

F4: on a REAL index a correct page citation is falsely flagged as wrong.

`check_answer_citations` resolves a ``«doc p.N»`` marker by comparing
``_locator_value(locator_for(doc, s)) == ("page", N)`` (cited_retrieval.py:277-288).
But the real indexer STRIPS the ``--- PAGE N ---`` marker from stored PDF
section text (doc_server/indexer.py:184) and stores no page field, so
``locator_for`` on a real-shaped PDF section returns a SECTION locator,
never a page. Every correct page citation therefore fails to resolve and
is flagged ``claim-wrong-locator``.

The builder's own phase-2 tests build a corpus that re-injects the marker
into section text (tests/test_ls1_answer_claims.py:_corpus_with) — a shape
the real producer never writes. These two tests use the real shape (marker
already stripped, source_type "pdf") so a genuine index behaves like the
indexer, not like the test fixture.

Both are RED on the current build (marker-less PDF section -> page citation
unresolvable). That is the finding: a real page citation cannot be
verified on a real index.
"""

from delfin.agent.cited_retrieval import check_answer_citations


def _real_pdf_corpus():
    """A pdf section EXACTLY as the real indexer stores it: no page marker
    in the text (indexer.py:184 strips it), source_type present, section_id
    + title present. This is what cite()/the audit reads on a real index.
    """
    return {
        "documents": {
            "paper-a": {
                "source_type": "pdf",
                "title": "Paper A",
                "sections": {
                    "sec_3": {
                        "section_id": "sec_3",
                        "title": "3 Results",
                        "text": "B3LYP reproduces the geometry well.",
                    }
                },
            }
        }
    }


def test_f4_correct_page_citation_is_clean_on_real_shape_id():
    # A real index has no page marker in the stored text; a correct page
    # citation must still be accepted, not flagged as pointing nowhere.
    corpus = _real_pdf_corpus()
    # Right page marker that a reader could go to in the source PDF.
    answer = "We use B3LYP as in \u00abpaper-a p.3\u00bb."
    claims = check_answer_citations(answer, corpus=corpus)
    assert not any(c.code == "claim-wrong-locator" for c in claims), (
        "a correct page citation must not be flagged wrong on a real-shaped index; "
        f"got {[c.code for c in claims]}"
    )
