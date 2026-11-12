"""LS1 phase 1: ``cite(query)`` returns passages with a locator.

Every passage an answer would quote from must carry somewhere the reader
can check: the page (for a PDF section, from the ``--- PAGE N ---`` marker
the indexer writes) or the section id + title. A passage that has neither
is never silently returned.

These tests drive the public path (``delfin.agent.cited_retrieval.cite``)
with a fully synthetic index dict, so they never touch a real index.
"""

from __future__ import annotations

from delfin.agent.cited_retrieval import cite


def _pdf_doc(sections=(("1 Introduction", "introduction text", 3),
                       ("2 Methods", "methods text", 7))):
    """A synthetic PDF-like document: page markers inside section text."""
    return {
        "doc-abc.pdf": {
            "title": "A Synthetic Paper",
            "source_type": "pdf",
            "reference": "S. Author, J. Synthetic 2024, 1, 2.",
            "sections": {
                f"sec_{i}": {
                    "section_id": f"sec_{i}",
                    "title": title,
                    "level": 1,
                    "text": f"--- PAGE {page} ---\n{text}",
                }
                for i, (title, text, page) in enumerate(sections)
            },
        }
    }


def _md_doc():
    """A synthetic markdown-like document: no page markers, sections only."""
    return {
        "own-report.md": {
            "title": "Our Own Report",
            "source_type": "markdown",
            "reference": "ComPlat, DELFIN internal report.",
            "sections": {
                "sec_0": {
                    "section_id": "sec_0",
                    "title": "Background",
                    "level": 1,
                    "text": "background text about embedding",
                }
            },
        }
    }


def _index(pdf=True, md=True, no_locator=False):
    docs = {}
    docs.update(_pdf_doc())
    docs.update(_md_doc())
    if no_locator:
        # A section with neither a page marker nor a section id/title:
        # the indexer could never write this, but cite must not hand it
        # back anyway.
        docs["broken-doc"] = {
            "title": "Broken",
            "source_type": "text",
            "reference": "",
            "sections": {
                "sec_0": {
                    "section_id": "",
                    "title": "",
                    "level": 0,
                    "text": "orphan text with no locator",
                }
            },
        }
    return {"version": 1, "built_at": "", "document_count": len(docs),
            "documents": docs}


def test_a_pdf_passage_is_returned_with_a_page_locator():
    # The PDF section carries a page marker; the passage must name that page.
    out = cite("introduction", index=_index())
    sources = out["sources"]
    assert sources, "expected at least one source for 'introduction'"
    pc = [s for s in sources if s["doc_id"] == "doc-abc.pdf"]
    assert pc, "pdf doc should be among the sources"
    assert pc[0]["locator"]["kind"] == "page"
    assert str(pc[0]["locator"]["value"]).isdigit()
    assert pc[0]["title"] == "A Synthetic Paper"


def test_a_section_passage_is_returned_with_a_section_locator():
    # The markdown doc has no pages; its locator is the section id + title.
    out = cite("embedding background", index=_index())
    md = [s for s in out["sources"] if s["doc_id"] == "own-report.md"]
    assert md, "markdown doc should be among the sources"
    assert md[0]["locator"]["kind"] == "section"
    assert md[0]["locator"]["value"]
    assert md[0]["reference"]


def test_no_result_without_a_locator():
    # The broken section (no marker, no id/title) must not be returned bare.
    out = cite("orphan text with no locator", index=_index(no_locator=True))
    assert all(s["locator"]["kind"] in ("page", "section")
               for s in out["sources"])
    broken = [s for s in out["sources"] if s["doc_id"] == "broken-doc"]
    assert not broken, "a passage with no locator must be dropped, not quoted"


def test_every_source_carries_a_locator_even_when_the_snippet_matches():
    # Run a query that definitely matches; each source must have a locator.
    out = cite("introduction", index=_index())
    assert out["sources"]
    for s in out["sources"]:
        assert "locator" in s, f"source without locator: {s['doc_id']!r}"
        assert s["locator"]["kind"] in ("page", "section")


def test_empty_index_says_no_source_found():
    empty = {"version": 1, "built_at": "", "document_count": 0, "documents": {}}
    out = cite("anything", index=empty)
    assert out["sources"] == []
