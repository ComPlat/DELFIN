"""Cited retrieval: search the local literature/manual index WITH a locator.

This is the citation half of LS1. ``cite(query)`` runs the same TF-IDF
search ``search_docs`` runs (delfin.doc_server.search.DocSearchEngine,
index loaded from delfin.doc_server.indexer.get_default_index_path) but
returns each passage together with the thing a reader can check: the page
number for a PDF section, or the section id and title otherwise. A passage
an answer would quote must carry a locator or it is not handed out.

Why the locator travels with the passage (delfin/doc_server/search.py:178):
the model reads the result, and a citation it has to fetch separately is
one it will not fetch. Same principle, one step out of the raw engine: the
engine returns ``section_id`` and ``snippet``; this module derives the
human-checkable page/section from the section's text and metadata.

Model-agnostic, read-only, never raises.
"""

from __future__ import annotations

import re
from typing import Any

__all__ = ["cite", "locator_for"]

#: The marker the indexer embeds at the top of each PDF page's text
#: (delfin/doc_server/indexer.py, ``_chunk_pdf_into_sections``). Its page
#: number is the reader's locator for a passage that came from a PDF.
_PAGE_MARKER = re.compile(r"--- PAGE (\d+) ---")


def load_index() -> dict | None:
    """The current default doc index, or None when none exists.

    Loads from the same path the ``search_docs`` gate uses
    (get_default_index_path) so a cite and a search never disagree about
    which corpus is authoritative.
    """
    try:
        from delfin.doc_server.indexer import get_default_index_path
        path = get_default_index_path()
        if not path.exists():
            return None
        import json
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def locator_for(doc: dict, section: dict) -> dict | None:
    """The human-checkable locator of one section, or None.

    A PDF section's text carries a ``--- PAGE N ---`` marker; the first one
    is where the passage starts, and that page is the locator. Any other
    document type uses the section id + title. When neither is available
    (no page marker, and the section id and title are both empty) there is
    nothing a reader can be pointed to, and None says "do not cite this".
    """
    text = section.get("text") or ""
    page = _PAGE_MARKER.search(text)
    if page:
        return {"kind": "page", "value": int(page.group(1))}
    section_id = (section.get("section_id") or "").strip()
    title = (section.get("title") or "").strip()
    if section_id or title:
        return {"kind": "section", "value": section_id or title}
    return None


def _run_engine(index: dict, query: str, doc_filter: str = "",
                max_results: int = 10) -> tuple[Any, list]:
    """Build the engine (if derivable) and search; engine, results."""
    from delfin.doc_server.search import DocSearchEngine
    engine = DocSearchEngine(index)
    try:
        payload = engine.search(query, doc_filter=doc_filter,
                                max_results=max_results)
    except Exception:
        payload = {"results": [], "index": {}}
    return engine, payload.get("results", [])


def cite(
    query: str,
    *,
    index: dict | None = None,
    max_per_doc: int = 3,
    max_docs: int = 5,
) -> dict[str, Any]:
    """Search the local literature/manual index for citable passages on
    *query*.

    Parameters
    ----------
    query : str
        Free-text search query.
    index : dict, optional
        The index dict to search. When omitted, the current default index
        (get_default_index_path) is loaded; an absent index is an empty
        corpus, reported as such, never invented.
    max_per_doc, max_docs : int
        Caps on how many passages per document, and how many documents,
        are returned — a cite should support an answer, not flood it.

    Returns
    -------
    dict
        ``{"sources": [...], "note": str}``. Each source is a passage with
        everything a reader needs to check it:
        ``{doc_id, title, locator: {"kind": "page"|"section",
        "value": int|str}, snippet, reference}``. A passage with NO locator
        is never returned. ``note`` carries why the result is small or
        empty ("no source found") so the caller never guesses.
    """
    if index is None:
        index = load_index() or {}

    docs = (index.get("documents")) or {}
    if not docs:
        return {"sources": [], "note": "no source found"}

    # Search across the whole corpus once; the engine ranks globally.
    _engine, results = _run_engine(index, query, max_results=max_docs * 4)
    if not results:
        return {"sources": [], "note": "no source found"}

    sources: list[dict[str, Any]] = []
    seen_docs: dict[str, int] = {}
    for hit in results:
        doc_id = hit.get("doc_id", "")
        if doc_id not in docs:
            continue
        if seen_docs.get(doc_id, 0) >= max_per_doc:
            continue
        section_id = hit.get("section_id", "")
        section = docs[doc_id].get("sections", {}).get(section_id, {})
        locator = locator_for(docs[doc_id], section)
        if locator is None:
            # The one guarantee this module exists for: never quote a
            # passage the reader cannot be pointed at.
            continue
        seen_docs[doc_id] = seen_docs.get(doc_id, 0) + 1
        sources.append({
            "doc_id": doc_id,
            "title": docs[doc_id].get("title", ""),
            "section_id": section_id,
            "locator": locator,
            "snippet": hit.get("snippet", ""),
            "reference": hit.get("reference", ""),
        })
        if len(sources) >= max_docs:
            break

    if not sources:
        return {"sources": [], "note": "no source found"}
    return {"sources": sources, "note": ""}
