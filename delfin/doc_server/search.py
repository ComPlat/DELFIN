"""TF-IDF search engine over the document index."""

from __future__ import annotations

import re
from typing import Any


class DocSearchEngine:
    """Search engine using TF-IDF with chemistry-aware tokenization.

    The engine is lazily initialized on the first query so that the MCP
    server starts up fast.
    """

    def __init__(self, index: dict) -> None:
        self._index = index
        self._vectorizer = None
        self._tfidf_matrix = None
        self._corpus_keys: list[tuple[str, str]] = []  # (doc_id, section_id)

    def _ensure_built(self) -> None:
        """Build the TF-IDF matrix if not already done."""
        if self._vectorizer is not None:
            return

        from sklearn.feature_extraction.text import TfidfVectorizer  # type: ignore

        corpus: list[str] = []
        keys: list[tuple[str, str]] = []

        for doc_id, doc in self._index.get("documents", {}).items():
            for section_id, section in doc.get("sections", {}).items():
                text = section.get("text", "")
                title = section.get("title", "")
                # Prepend title for boosted matching
                corpus.append(f"{title}\n{text}")
                keys.append((doc_id, section_id))

        if not corpus:
            self._corpus_keys = []
            return

        self._corpus_keys = keys

        # ``max_df=0.95`` is a FRACTION of the corpus: with a single section
        # the ceiling floors to zero documents, below ``min_df=1``, and the
        # vectoriser raises "max_df corresponds to < documents than min_df"
        # — so the first PDF a user indexed broke search_docs outright, on
        # every call, with an exception rather than an empty result. A
        # one-section corpus has no common terms to drop anyway.
        max_df = 0.95 if len(corpus) > 1 else 1.0

        # Custom token pattern that preserves hyphenated chemistry terms
        # like def2-TZVP, wB97X-D3, RIJCOSX, SARC/J, etc.
        self._vectorizer = TfidfVectorizer(
            token_pattern=r"(?u)\b[\w\-/\.]+\b",
            ngram_range=(1, 2),
            max_features=50000,
            sublinear_tf=True,
            min_df=1,
            max_df=max_df,
        )
        self._tfidf_matrix = self._vectorizer.fit_transform(corpus)

    def provenance(self) -> dict[str, Any]:
        """When this index was built, how big it is, what failed, what is stale."""
        from .indexer import index_provenance
        prov = index_provenance(self._index)
        prov["corpus_sections"] = len(self._corpus_keys)
        return prov

    def search(
        self,
        query: str,
        doc_filter: str = "",
        max_results: int = 10,
    ) -> dict[str, Any]:
        """Search the index for sections matching the query.

        Parameters
        ----------
        query : str
            Free-text search query.
        doc_filter : str, optional
            Restrict results to a specific ``doc_id``.
        max_results : int
            Maximum number of results to return.

        Returns
        -------
        dict
            ``{"results": [...], "index": {...}}``. Each result:
            ``{doc_id, section_id, title, doc_title, score, snippet}``.

            The payload carries the index provenance — build time, corpus
            size, documents that could not be extracted, documents whose
            source has changed since — because a bare list made "nothing
            matched your query" and "there is nothing here to match" the
            same answer, and made a hit from a section deleted from the
            source last month indistinguishable from a fresh one.
        """
        return {"results": self._ranked(query, doc_filter, max_results),
                "index": self.provenance()}

    def _ranked(
        self,
        query: str,
        doc_filter: str = "",
        max_results: int = 10,
    ) -> list[dict[str, Any]]:
        """The ranked hits alone, without the provenance wrapper."""
        self._ensure_built()

        if self._vectorizer is None or self._tfidf_matrix is None or not self._corpus_keys:
            return []

        from sklearn.metrics.pairwise import cosine_similarity  # type: ignore
        import numpy as np  # type: ignore

        query_vec = self._vectorizer.transform([query])
        scores = cosine_similarity(query_vec, self._tfidf_matrix).flatten()

        # Apply doc_filter
        if doc_filter:
            for i, (doc_id, _) in enumerate(self._corpus_keys):
                if doc_id != doc_filter:
                    scores[i] = 0.0

        # Get top results
        top_indices = np.argsort(scores)[::-1][:max_results]

        results: list[dict[str, Any]] = []
        docs = self._index.get("documents", {})
        for idx in top_indices:
            score = float(scores[idx])
            if score < 1e-6:
                break

            doc_id, section_id = self._corpus_keys[idx]
            doc = docs.get(doc_id, {})
            section = doc.get("sections", {}).get(section_id, {})
            text = section.get("text", "")

            # Build a snippet (first 300 chars)
            snippet = text[:300].replace("\n", " ").strip()
            if len(text) > 300:
                snippet += "..."

            results.append({
                "doc_id": doc_id,
                "section_id": section_id,
                "title": section.get("title", ""),
                "doc_title": doc.get("title", ""),
                "score": round(score, 4),
                "snippet": snippet,
            })

        return results
