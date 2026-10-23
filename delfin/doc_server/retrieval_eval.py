"""Measuring retrieval before changing it.

Input: a built index. Output: recall@k and MRR over a set of queries
DERIVED FROM THE CORPUS, so it exists for any user's documents and needs
nobody to write it.

How a query is made: each section contributes its own title AND a run of
words from its middle, and the document it came from is the answer. The
middle matters -- a truncating indexer keeps the beginning, so a set that
asks only about beginnings cannot see what such an indexer lost. Section titles are the authors' or
the typesetter's words, not a judgement of what matters, and a section
whose title is a bare number or a stock word ("Introduction", "Part 3")
contributes nothing -- it names no document.

What this measures, and what it does not. It measures RETRIEVABILITY: can
a document be reached at all by language taken from it. That is exactly
the failure it was built to catch -- three of four papers reached the
index with 24-49% of their text, so no wording could reach the rest of
them. It does NOT measure semantic generalisation: every query here
occurs in the corpus, so a lexical method is flattered. A set that
measured paraphrase would have to be written by a person or a model, and
would then measure whoever wrote it.

It is also DOCUMENT-level, and that is forgiving in a way worth knowing.
Measured against the chunking change that made every document reach the
index whole (24-49% -> 100% for the papers):

    before (truncating)   recall@1 0.918   recall@10 0.984   MRR 0.950
    after  (windowing)    recall@1 0.975   recall@10 1.000   MRR 0.988

Real, and modest -- because a document indexed at a quarter of its length
is still FOUND, by whichever quarter survived. What that change fixed was
what the agent can READ once it has arrived, which no document-level
figure sees. A passage-level measure would, and does not exist here yet.

Read the numbers accordingly: a fall is evidence, a rise is a hypothesis
until something outside this file agrees.
"""

from __future__ import annotations

import re
from typing import Any

#: Titles that name no document and make no query.
_EMPTY_TITLE = re.compile(
    r"^(part|front matter|full document|chapter|section|appendix|"
    r"introduction|conclusions?|references|abstract|methods?|results?|"
    r"discussion|acknowledge?ments?|supporting information)\b[\s\d.()]*$",
    re.I)


def _usable(title: str) -> bool:
    text = " ".join(str(title or "").split())
    if len(text) < 12 or _EMPTY_TITLE.match(text):
        return False
    # A heading that is mostly digits is a numbering, not a name.
    letters = sum(c.isalpha() for c in text)
    return letters >= max(8, len(text) // 2)


def _body_phrase(text: str, words: int = 12) -> str:
    """A run of words from the MIDDLE of *text*, or "".

    The middle, not the start: the start of a section is what a truncating
    indexer keeps, so a query taken from there is answered by an index
    that dropped everything else. The first version of this file asked
    only section titles and reported a 0.927 -> 0.957 improvement for a
    change that added 706,000 characters of body text -- because titles
    survived truncation too, so the set could not see what changed.
    """
    parts = str(text or "").split()
    if len(parts) < words * 3:
        return ""
    start = len(parts) // 2
    run = [w for w in parts[start:start + words] if any(c.isalnum() for c in w)]
    return " ".join(run) if len(run) >= words - 3 else ""


def queries(index: dict, per_document: int = 12) -> list[tuple[str, str]]:
    """(query, the doc_id that answers it), derived from *index*.

    Capped per document so one very large document -- the ORCA manual has
    2373 sections against a paper's dozen -- cannot decide the score on
    its own. The cap is applied over the longest sections, because a
    heading with substance under it is a heading somebody wrote.
    """
    out: list[tuple[str, str]] = []
    for doc_id, doc in (index.get("documents") or {}).items():
        sections = list((doc.get("sections") or {}).values())
        sections.sort(key=lambda s: len(s.get("text", "")), reverse=True)
        taken = 0
        for section in sections:
            title = section.get("title", "")
            if not _usable(title):
                continue
            out.append((" ".join(str(title).split()), doc_id))
            # And one from inside the same section. A title is a label and
            # survives most indexing faults; the body is what a reader
            # actually asks about, and what an indexer drops.
            phrase = _body_phrase(section.get("text", ""))
            if phrase:
                out.append((phrase, doc_id))
            taken += 1
            if taken >= per_document:
                break
    return out


def score(engine: Any, cases: list[tuple[str, str]], *,
          depth: int = 10) -> dict:
    """Run *cases* against *engine*. Returns recall@1, recall@k and MRR.

    A query whose document does not appear within *depth* contributes 0 to
    every figure rather than being dropped: a document that cannot be
    found is the result, not a missing measurement.
    """
    if not cases:
        return {"queries": 0, "recall@1": 0.0, f"recall@{depth}": 0.0,
                "mrr": 0.0, "unreachable": []}
    hits1 = hitsk = 0
    reciprocal = 0.0
    unreachable: list[str] = []
    for query, want in cases:
        try:
            results = engine.search(query, max_results=depth).get("results", [])
        except Exception:
            results = []
        rank = next((i + 1 for i, r in enumerate(results)
                     if r.get("doc_id") == want), None)
        if rank is None:
            unreachable.append(query)
            continue
        hitsk += 1
        hits1 += (rank == 1)
        reciprocal += 1.0 / rank
    n = len(cases)
    return {
        "queries": n,
        "recall@1": round(hits1 / n, 4),
        f"recall@{depth}": round(hitsk / n, 4),
        "mrr": round(reciprocal / n, 4),
        # Named, not counted: the ones worth reading are the misses.
        "unreachable": unreachable[:10],
    }


def report(index: dict, *, depth: int = 10, per_document: int = 12) -> dict:
    """Build the queries from *index* and score *index*'s own engine."""
    from .search import DocSearchEngine

    cases = queries(index, per_document=per_document)
    out = score(DocSearchEngine(index), cases, depth=depth)
    out["documents"] = len(index.get("documents") or {})
    return out
