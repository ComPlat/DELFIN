"""LS1 phase-1 fix: the locator's kind comes from the document's type.

Reviewer finding F2/F3 (nacht-s29): ``locator_for`` derived the page
locator from the section's stored *text*, so a pageless markdown document
whose body quotes ``--- PAGE 9 ---`` was mislabelled as page 9, and
``--- PAGE 0 ---`` produced page 0. The kind must come from the
document's ``source_type`` (never from text content), and a page number
must be a real page (>= 1).
"""

from __future__ import annotations

from delfin.agent.cited_retrieval import locator_for


def _pdf_sec(text):
    return {
        "section_id": "sec_0",
        "title": "Methods",
        "level": 1,
        "text": text,
    }


def _doc(source_type):
    return {"title": "X", "source_type": source_type, "reference": ""}


def test_a_markdown_page_marker_is_not_a_page_locator():
    # A markdown document's body quoting "--- PAGE 9 ---" is content, not a
    # real page: the locator must be a section, never a page.
    loc = locator_for(_doc("markdown"),
                      _pdf_sec("--- PAGE 9 ---\nquoted text"))
    assert loc is not None
    assert loc["kind"] != "page", "markdown body must not become a page locator"
    assert loc["kind"] == "section"


def test_page_zero_is_not_a_checkable_page():
    # Page 0 cannot be checked by a reader; never emit it as a page locator.
    loc = locator_for(_doc("pdf"), _pdf_sec("--- PAGE 0 ---\ntext"))
    assert loc is not None
    assert loc["kind"] != "page", "page 0 must not be handed out as a locator"
