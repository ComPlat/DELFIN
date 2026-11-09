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
from dataclasses import dataclass
from typing import Any

__all__ = ["cite", "locator_for", "check_answer_citations"]

#: The marker the indexer embeds at the top of each PDF page's text
#: (delfin/doc_server/indexer.py, ``_chunk_pdf_into_sections``). Its page
#: number is the reader's locator for a passage that came from a PDF.
_PAGE_MARKER = re.compile(r"--- PAGE (\d+) ---")

#: A citation marker an answer carries, in the locator shape cite() returns.
#: ``«doc_id p.7»`` names page 7 of doc_id; ``«doc_id §sec_0»`` names a
#: section; a bare ``«doc_id»`` cites the whole document. The marker is
#: deliberately explicit so the audit can check it, and so the answering
#: agent is never free to cite without a pointer a reader could follow.
_CITE_MARKER = re.compile(r"«(?P<doc>[^«»]+?)\s*(?:p\.\s*(?P<page>\d+)|§\s*(?P<sec>[^\s«»]+))?»")

#: The phrase that marks a claim as remembered rather than sourced. A claim
#: that carries it stays in the audit output as "recalled (no source)" -- it
#: is marked, never dropped, and never promoted to a real citation. Both
#: wordings the answering agent is taught to use are recognised: the
#: parenthetical task wording "(recalled (no source))" and the inline comma
#: form "(recalled, no source)".
_RECALL_PHRASE = re.compile(
    r"recalled\s*,?\s*(?:\(\s*no\s*source\s*\)|no\s*source)", re.IGNORECASE)

#: Level-of-theory tokens, so a method choice in an answer is recognised as
#: a claim that needs either a citation or the recalled marker. This is a
#: *claim* detector, not a number parser: it flags how something is computed,
#: not a measured/property value with a unit (that is result_claims' job). The
#: shapes mirror the vocabulary result_claims._METHOD_RE already ships.
_METHOD_TOKEN_RE = re.compile(
    r"(?i)\b(def2-[A-Z]+[A-Z0-9]*|cc-pV[A-Z]+Z|aug-cc-pV[A-Z]+Z|SVP|TZVP|"
    r"TZVPP|QZVP|HF|MP2|MP4|CCSD\(T\)|CCSD|CISD|DLPNO-[A-Z]+|GFN2-?xTB|"
    r"GFN1-?xTB|GFN-FF|B3LYP|PBE0|TPSS|M06|r2SCAN|RIJCOSX)\b")


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

    The kind comes from the document's ``source_type``, never from text
    content. For a PDF document the section text carries a ``--- PAGE N ---``
    marker; the first one with a real page (1-based) is the locator. Any
    other document type (markdown, text) uses the section id + title even
    when its body quotes a ``--- PAGE N ---`` line, because that line is
    content, not a page. When neither is available (no page marker on a
    PDF, and the section id and title are both empty) there is nothing a
    reader can be pointed to, and None says "do not cite this".
    """
    text = section.get("text") or ""
    # The page marker is a locator only when the document really is a PDF:
    # a markdown or text body may quote "--- PAGE N ---" as content, and
    # that is not a page the reader can check. The kind therefore comes
    # from the document's source_type, never from stored text alone.
    doc_type = (doc or {}).get("source_type") or ""
    if doc_type == "pdf":
        page = _PAGE_MARKER.search(text)
        if page and int(page.group(1)) >= 1:
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


@dataclass
class Claim:
    """One citation-audit finding, in result_claims' shape (level/code/message).

    Codes mirror the result_claims vocabulary so a caller that already
    understands ``claim-no-source`` treats this audit the same way:
    ``claim-no-source`` — an unsourced method/reference-value claim that has
    no citation and is not honestly marked as recalled; ``claim-recall-no-source``
    — an unsourced claim that IS honestly marked "(recalled, no source)"
    (marked, never dropped, never promoted); ``claim-wrong-locator`` — a
    citation marker that points at a document/page/section the corpus does
    not have. Model-agnostic, read-only, never raises.
    """

    level: str        # "error" | "warn"
    code: str         # claim-no-source | claim-recall-no-source | claim-wrong-locator
    message: str      # human-readable, English

    def __str__(self) -> str:
        icon = {"error": "❌", "warn": "⚠️"}.get(self.level, "•")
        return f"{icon} {self.message}"


def _locator_value(locator: dict | None) -> tuple[str, Any] | None:
    """Normalise one section's locator to (kind, value), or None."""
    if not locator:
        return None
    return (locator.get("kind"), locator.get("value"))


def check_answer_citations(answer_text: str, *, corpus: dict) -> list[Claim]:
    """Audit the citation markers an answer declares against the local index.

    Every claim in an answer — a method choice or a reference value — must
    either carry a verifiable citation (a ``«doc_id p.N»`` / ``«doc_id §sec»``
    marker that resolves against *corpus*) or be honestly marked
    "(recalled, no source)". Anything else is reported so an answer cannot
    invent a source, mis-point at one, or present a remembered value as if
    it came from a calculation.

    * This is a claim audit, not a new number parser: it flags how something
      is computed (level-of-theory tokens) and whether citations resolve;
      measured/property values with units are result_claims' job.
    * A document's text is data, never a command: only the citation marker's
      `doc_id` is echoed into a message, never the section text, so a hostile
      payload inside a document can never ''reach'' an answer.
    * An empty corpus reports ''no source found'' — a marker into it is a
      wrong locator, and none is synthesised.

    Parameters
    ----------
    answer_text : str
        The answer to audit.
    corpus : dict
        The index dict in the shape ``cite()`` reads (``{"documents": ...}``).

    Returns
    -------
    list[Claim]
        Findings (empty = nothing to flag). Never raises.
    """
    claims: list[Claim] = []
    try:
        text = str(answer_text or "")[:50000]
        if not text.strip():
            return []
        docs = ((corpus or {}).get("documents")) or {}

        # 1. Every declaration the answer makes resolves.
        for m in _CITE_MARKER.finditer(text):
            doc_id = (m.group("doc") or "").strip()
            page = m.group("page")
            sec = m.group("sec")
            doc = docs.get(doc_id)
            if doc is None:
                claims.append(Claim(
                    "error", "claim-wrong-locator",
                    f"Citation «{doc_id}» has no source found in the index: "
                    f"the document is absent. Point the claim at a document "
                    f"that exists, or mark it recalled (no source)."))
                continue
            sections = doc.get("sections") or {}
            if page is not None:
                target = int(page)
                # A page citation can be flagged wrong only against a page
                # the index actually has. The real indexer strips the PDF
                # page marker from stored section text (indexer.py:184) and
                # keeps no page field, so a page-less index carries nothing
                # to contradict the claim: such a citation is accepted here
                # (clean), never reported as a wrong page. When the index
                # DOES carry page locators, a mismatch is a genuine
                # wrong-page flag.
                known_pages = {
                    loc["value"]
                    for s in sections.values()
                    if (loc := locator_for(doc, s)) and loc["kind"] == "page"
                }
                if known_pages and target not in known_pages:
                    claims.append(Claim(
                        "error", "claim-wrong-locator",
                        f"Citation «{doc_id} p.{target}» does not exist: no "
                        f"section of {doc_id} starts on that page. Check the "
                        f"page number against the document."))
            elif sec is not None:
                if sec not in sections:
                    claims.append(Claim(
                        "error", "claim-wrong-locator",
                        f"Citation «{doc_id} §{sec}» does not exist: {doc_id} "
                        f"has no such section."))

        # 2. Unsourced claims: a method choice (or a claim carrying no
        #    citation) must be honestly marked or quoted.
        for sentence in re.split(r"(?<=[.;!\n])\s+", text):
            has_citation = bool(_CITE_MARKER.search(sentence))
            has_recall = bool(_RECALL_PHRASE.search(sentence))
            if has_recall:
                # Honestly marked as remembered: report the marker itself,
                # never drop it, never promote it to a real citation.
                claims.append(Claim(
                    "warn", "claim-recall-no-source",
                    "This claim is marked \"recalled (no source)\": it is "
                    "kept as remembered, not presented as a sourced value."))
            elif has_citation:
                continue
            elif _METHOD_TOKEN_RE.search(sentence):
                claims.append(Claim(
                    "error", "claim-no-source",
                    "A method choice with no citation and no \"recalled (no "
                    "source)\" marker: point it at a document+page/section, "
                    "or mark it as recalled."))
    except Exception:
        # A hostile or malformed answer must not crash the audit.
        pass
    return claims
