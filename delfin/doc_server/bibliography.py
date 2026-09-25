"""What a document is, in the form a reader can check.

Input: a document's own metadata and the text of its first page. Output:
``title``, ``authors``, ``year``, ``journal``, ``doi`` and a one-line
``reference`` built from whichever of those exist.

Why. The index already held papers, and a search result named them by
filename: ``1 s2.0 S0021979724009044 main``. An agent could therefore read
a paper and recommend what it said, but could not tell anyone WHERE it
read it, which is the difference between an opinion and a citable claim.

Where the fields come from, in the order they are trusted. Measured over
the corpus rather than assumed:

    /Subject   the richest source -- publishers put the citation itself
               there (journal, volume, year, doi). Present and useful in
               3 of 4 papers measured.
    /Title     the real title in every one of the 4, not the typesetter
               garbage this field often carries elsewhere.
    /Author    author list, as printed.
    page 1     searched for a DOI only, and only when the metadata has
               none: it carried one in 1 of the 4.

No network. A lookup against a bibliographic service would resolve more,
and would also send a list of what the user reads to a third party and
fail on a machine with no route out. Both are worse than a reference with
a field missing, so every field is derived from the file itself.

Nothing here raises, and nothing here is required: a document with no
usable metadata keeps the title the indexer already gave it, which is
where this started. The floor is the current behaviour.
"""

from __future__ import annotations

import re
from typing import Any

#: A DOI as registered: "10." then a registrant code, then anything until
#: whitespace or a character that only ever ends one. Trailing punctuation
#: is stripped separately -- a DOI may legitimately contain a dot, so the
#: pattern cannot simply forbid it at the end.
_DOI_RE = re.compile(r"\b10\.\d{4,9}/[^\s\"'<>,;]+", re.I)
_YEAR_RE = re.compile(r"\b(19\d{2}|20\d{2})\b")
_TRAILING = ".,;:)]}"

#: Titles some producers write into /Title. None of them names a document.
_JUNK_TITLE_RE = re.compile(
    r"^(untitled|microsoft word\s*-|document\d*|print job|\d+\.(pdf|docx?)$)",
    re.I)


def doi_in(text: str) -> str:
    """The first DOI in *text*, lower-cased and unpunctuated, or "".

    A DOI is compared and resolved case-insensitively, so one spelling is
    stored. The trailing strip matters because a DOI at the end of a
    sentence and a DOI in a URL both arrive with characters that are not
    part of it.
    """
    match = _DOI_RE.search(str(text or ""))
    if not match:
        return ""
    return match.group(0).rstrip(_TRAILING).lower()


def year_in(text: str) -> str:
    """The first plausible publication year in *text*, or ""."""
    match = _YEAR_RE.search(str(text or ""))
    return match.group(0) if match else ""


def journal_in(subject: str) -> str:
    """The journal name from a /Subject line, or "".

    Publishers write the citation there in one of two shapes -- a name
    followed by a comma, or a name followed by the year. Everything after
    the first of those is volume, pages and doi, which the other fields
    already hold.
    """
    text = " ".join(str(subject or "").split())
    if not text:
        return ""
    head = re.split(r"\s*[,;]\s*|\s+(?=(?:19|20)\d{2}[.:])", text, maxsplit=1)[0]
    head = _DOI_RE.sub("", head).strip(" .,;:")
    # A /Subject that is a sentence is an abstract, not a citation.
    return head if 0 < len(head) <= 120 and head.count(" ") <= 12 else ""


def _clean(value: Any) -> str:
    text = " ".join(str(value or "").split())
    return "" if _JUNK_TITLE_RE.match(text) else text


def describe(metadata: dict | None = None, first_page: str = "",
             fallback_title: str = "") -> dict:
    """Everything known about one document, as fields.

    *metadata* is the PDF's own dictionary (pypdf spelling, ``/Title`` and
    friends); missing keys and a missing dictionary are the same thing.
    *fallback_title* is what the indexer would have used on its own, and is
    what comes back when the file says nothing about itself.
    """
    meta = {str(k): v for k, v in (metadata or {}).items()}
    subject = _clean(meta.get("/Subject"))
    title = _clean(meta.get("/Title")) or _clean(fallback_title)
    authors = _clean(meta.get("/Author"))
    doi = doi_in(subject) or doi_in(meta.get("/Keywords", "")) or doi_in(
        first_page[:4000])
    return {
        "title": title,
        "authors": authors,
        "year": year_in(subject) or "",
        "journal": journal_in(subject),
        "doi": doi,
    }


def _short_authors(authors: str) -> str:
    """"Jin and Merz Jr." -> "Jin et al."; a single name is left alone.

    Shortened only when there IS a list. Taking the last word of a lone
    author turned "Max-Planck-Institut fuer Kohlenforschung" into
    "Kohlenforschung" -- measured on the corpus, not imagined. A surname
    is the right handle for one name out of several; for one name on its
    own there is nothing to abbreviate, and guessing which words are a
    family name is guessing about people and organisations at once.
    """
    text = " ".join(str(authors or "").split())
    if not text:
        return ""
    parts = [p for p in re.split(r"\s*(?:,| and | & |;)\s*", text) if p.strip()]
    if len(parts) < 2:
        return text
    words = parts[0].split()
    return f"{words[-1]} et al." if words else text


def reference(fields: dict | None) -> str:
    """One line a reader can act on, or "".

    Built from whatever is present, in the order a citation is read:
    author, year, title, journal, doi. A field that is missing is left out
    rather than filled with a placeholder -- "n.d." and "Anon." look like
    data and are not.
    """
    f = fields or {}
    author = _short_authors(f.get("authors", ""))
    year = str(f.get("year") or "")
    title = str(f.get("title") or "")
    journal = str(f.get("journal") or "")
    doi = str(f.get("doi") or "")
    head = " ".join(x for x in (author, f"({year})" if year else "") if x)
    parts = [x for x in (head, title, journal) if x]
    line = ". ".join(parts)
    if doi:
        line = f"{line}. doi:{doi}" if line else f"doi:{doi}"
    return line.strip(" .") + ("" if not line else "")
