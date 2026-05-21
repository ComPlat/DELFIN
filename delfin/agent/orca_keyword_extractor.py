"""Extract ORCA keyword ground-truth from the indexed manual.

The DELFIN doc-index (``~/.delfin/doc_index.json``) holds the ORCA
6.1.1 manual broken into ~1250 sections.  This module scans those
sections to build a structured keyword namespace:

    {
      "casscf": {
        "block": "%casscf",
        "keywords": ["nel", "norb", "mult", "nroots", "ptmethod", …],
        "sections": ["ch3_13_…", "ch3_18_…"],
      },
      "tddft":  { … },
      "mp2":    { … },
      …
    }

Used to:

1. Generate `fact_verify` benchmark tasks whose expected/forbidden
   patterns are sourced from the manual, NOT from author memory.
2. Provide a programmatic ground-truth for "is keyword X a real
   ORCA keyword?" checks at runtime.

Design choices:

- Extraction is regex-based, not LLM-based — deterministic + cheap.
- Only obvious patterns are extracted (``%block`` headers + indented
  ``key value`` lines within manual examples).  False-positives are
  filtered via a stop-list of generic words.
- Output is committed to the repo as a snapshot — the manual is
  fixed for a given ORCA version, so the JSON only needs re-extraction
  when the manual itself updates.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


_DEFAULT_INDEX = Path.home() / ".delfin" / "doc_index.json"

# Known ORCA computation blocks worth extracting.  Each entry:
#   block_name: aliases used in the manual to introduce the block
_KNOWN_BLOCKS: dict[str, tuple[str, ...]] = {
    "casscf":  ("%casscf",),
    "tddft":   ("%tddft", "%tda"),
    "mp2":     ("%mp2", "%mdci"),
    "ccsd":    ("%mdci", "%ccsd"),
    "scf":     ("%scf",),
    "method":  ("%method",),
    "basis":   ("%basis",),
    "output":  ("%output",),
    "freq":    ("%freq",),
    "geom":    ("%geom",),
    "irc":     ("%irc",),
    "neb":     ("%neb",),
    "eprnmr":  ("%eprnmr", "%nmr"),
    "elprop":  ("%elprop",),
    "rocis":   ("%rocis", "%cis"),
    "mrci":    ("%mrci",),
    "cipsi":   ("%cipsi",),
    "esd":     ("%esd",),
    "compound": ("%compound",),
    "pal":     ("%pal",),
}

# Words that LOOK like keywords but are actually prose noise.
_KEYWORD_STOPLIST = frozenset({
    "and", "or", "the", "for", "with", "is", "are", "this", "that", "to",
    "of", "in", "by", "on", "as", "if", "not", "be", "we", "you", "it",
    "see", "note", "also", "such", "any", "all", "one", "two", "three",
    "value", "default", "true", "false", "yes", "no", "end", "begin",
    "section", "table", "figure", "example", "user", "manual", "block",
    "keyword", "input", "output", "file", "line", "page", "chapter",
    "fig", "tab", "eq", "ref", "doc",
})


def _load_index(path: Path | None = None) -> dict[str, Any]:
    p = path or _DEFAULT_INDEX
    with p.open(encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list) and data:
        data = data[0]
    return data


def _orca_sections(index: dict[str, Any]) -> dict[str, dict[str, Any]]:
    docs = index.get("documents", {})
    for doc_id, doc in docs.items():
        if "orca" in doc_id.lower():
            return doc.get("sections", {})
    return {}


def _extract_block_keywords(
    block_aliases: tuple[str, ...],
    sections: dict[str, dict[str, Any]],
) -> tuple[set[str], list[str]]:
    """Return ``(keywords, section_ids)`` for a given block.

    Keywords are identified by:
    1. Sections whose text contains ``%blockname``.
    2. Within those sections, the lines following the ``%blockname``
       up to ``end`` are scanned for ``\\w+`` tokens that look like
       parameter names (lowercase, alphanumeric, not in stop-list).
    """
    keywords: set[str] = set()
    section_ids: list[str] = []
    alias_re = re.compile(
        "|".join(re.escape(a) for a in block_aliases),
        re.IGNORECASE,
    )
    for sid, sec in sections.items():
        text = sec.get("text", "") or ""
        if not text:
            continue
        if not alias_re.search(text):
            continue
        section_ids.append(sid)
        # Scan all instances of the block in this section
        for m in alias_re.finditer(text):
            # Take the slice after the block opening up to "end" (or
            # 600 chars max — most blocks are well under that).
            start = m.end()
            tail = text[start:start + 1500]
            # Try to clip at first standalone "end"
            end_m = re.search(r"\bend\b", tail, re.IGNORECASE)
            if end_m:
                tail = tail[:end_m.start()]
            # Extract candidate keywords: lowercase words at start of
            # a non-empty line, OR after a comma / pipe / semicolon.
            for line in tail.splitlines():
                line = line.strip()
                if not line or line.startswith(("#", "//")):
                    continue
                # First whitespace-delimited token of the line
                first = re.match(r"([A-Za-z][A-Za-z0-9_]{1,30})", line)
                if first:
                    kw = first.group(1).lower()
                    if (kw not in _KEYWORD_STOPLIST
                            and len(kw) >= 2
                            and not kw.startswith("end")
                            and not kw.isdigit()):
                        keywords.add(kw)
    return keywords, section_ids


def extract_keyword_namespace(
    index_path: Path | None = None,
) -> dict[str, dict[str, Any]]:
    """Walk the indexed ORCA manual and return a ground-truth keyword
    namespace keyed by block name."""
    index = _load_index(index_path)
    sections = _orca_sections(index)
    out: dict[str, dict[str, Any]] = {}
    for block_name, aliases in _KNOWN_BLOCKS.items():
        keywords, section_ids = _extract_block_keywords(aliases, sections)
        if not keywords and not section_ids:
            continue
        out[block_name] = {
            "block": aliases[0],
            "aliases": list(aliases),
            "keywords": sorted(keywords),
            "sections": sorted(set(section_ids))[:20],
            # Snapshot-time count — useful for "did the manual change?"
            "n_keywords": len(keywords),
            "n_sections": len(set(section_ids)),
        }
    return out


def is_real_keyword(block: str, keyword: str,
                    namespace: dict[str, dict[str, Any]] | None = None,
                    path: Path | None = None) -> bool:
    """True if ``keyword`` appears in the ORCA manual for ``block``.

    Use to runtime-check whether a model's claimed keyword is real
    or hallucinated.
    """
    if namespace is None:
        namespace = extract_keyword_namespace(path)
    info = namespace.get(block.lower())
    if not info:
        return False
    return keyword.lower() in {k.lower() for k in info.get("keywords", [])}


__all__ = [
    "extract_keyword_namespace",
    "is_real_keyword",
]
