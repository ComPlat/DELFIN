"""Source ledger: records claims with their sources and flags weak research.

Pure library, no I/O, no network. W1 package -- web research that is
"researched with sources, never steered by a page". A research answer
carries claims; each claim needs at least two INDEPENDENT sources, each
with a retrieval date, and claims that two independent sources state
differently are conflicts. The ledger records the claim -> sources mapping
and reviews it; it does not judge truth, only sourcing.

Source record shape (a plain dict):
    {"url": str, "kind": "primary" | "secondary",
     "retrieved": "YYYY-MM-DD" ("" = unknown), "quote": str (optional span)}
"""
from __future__ import annotations

import re
from typing import Iterable

_KINDS = ("primary", "secondary")

# A claim key: normalised enough that "The SCS overlap is negative" and
# "the SCS overlap is negative." collide, loose enough that different
# claims stay different.
_PUNCT = re.compile(r"[^\w\s]")
_WS = re.compile(r"\s+")


def _claim_key(claim: str) -> str:
    text = _PUNCT.sub(" ", claim.lower())
    return _WS.sub(" ", text).strip()


def _unique_independent_sources(sources: Iterable[dict]) -> list[dict]:
    """Distinct-source view of a claim.

    Two sources are the SAME voice when they share a URL, or when one is
    a secondary: a secondary restates a primary's voice, it does not add
    one. Independence therefore needs >= 2 distinct PRIMARIES WITH a
    non-empty URL: a primary without a URL cannot be deduplicated
    against anything, so it can never serve as the second voice (two
    URL-less primaries are one paper cited twice, unverifiable).
    """
    primaries = 0
    any_source = False
    seen_urls: set[str] = set()
    for src in sources:
        any_source = True
        url = (src.get("url") or "").strip().lower()
        kind = (src.get("kind") or "primary").strip().lower()
        if kind != "primary":
            continue
        if not url:
            continue  # URL-less primary: not deduplicable, never a voice
        if url in seen_urls:
            continue
        seen_urls.add(url)
        primaries += 1
    if not any_source:
        return []
    if primaries >= 2:
        return [{"_voice": i} for i in range(primaries)]
    return [{"_voice": 0}]


def review_sources(claim: str, sources: list[dict]) -> list[dict]:
    """Issues for ONE claim. Issue dicts: {"claim", "type", ...detail}."""
    issues: list[dict] = []
    indep = _unique_independent_sources(sources)
    if len(indep) < 2:
        issues.append({"claim": claim, "type": "few_sources",
                       "n_sources": len(indep)})
    for src in sources:
        if not (src.get("retrieved") or "").strip():
            issues.append({"claim": claim, "type": "no_retrieved_date",
                           "url": (src.get("url") or "").strip()})
    return issues


def review_conflicts(records: list[tuple[str, list[dict]]]) -> list[dict]:
    """Claims that independent sources state differently.

    Two claims whose keys differ in polarity count as a conflict when
    each has >=1 source. Polarity is read from a word list of negations
    AND from antonym pairs (negative/positive), the two shapes claim
    contradictions actually take. Deliberately conservative: only
    polarity flips, not every wording difference, is a conflict.
    """
    _NEG = (" not ", " no ", " never ", "does not", "is not", "are not",
            "cannot", "can not", "negative")
    _POS = ("positive",)
    _STEM_SWAPS = [
        (" not ", " "), (" no ", " "), (" never ", " "),
        ("does not", "does"), ("is not", "is"), ("are not", "are"),
        ("cannot", "can"), ("can not", "can"),
        ("negative", "polarity-stem"),
        ("positive", "polarity-stem"),
    ]
    by_stem: dict[str, list[tuple[str, bool]]] = {}
    for claim, sources in records:
        key = _claim_key(claim)
        text = f" {key} "
        polarity = None
        if any(word in text for word in _NEG):
            polarity = "neg"
        elif any(word in text for word in _POS):
            polarity = "pos"
        stem = text
        for old, new in _STEM_SWAPS:
            stem = stem.replace(old, new)
        by_stem.setdefault(stem, []).append((claim, polarity))
    conflicts: list[dict] = []
    for stem, variants in by_stem.items():
        pols = {p for _, p in variants}
        # Polarity sets differing (incl. neg-vs-None: one source says a
        # method does not proceed, the other says it does) is a conflict;
        # two claims with the SAME polarity are not.
        if len(variants) > 1 and len(pols) > 1:
            conflicts.append({"type": "conflict",
                              "claims": [c for c, _ in variants],
                              "stem": stem.strip()})
    return conflicts


def review(records: list[tuple[str, list[dict]]]) -> list[dict]:
    """Full ledger review: per-claim issues + conflicts."""
    issues: list[dict] = []
    for claim, sources in records:
        issues.extend(review_sources(claim, sources))
    issues.extend(review_conflicts(records))
    return issues


def check_ledger(records: list[tuple[str, list[dict]]]) -> dict:
    """One-shot report over raw (claim, sources) pairs."""
    return {"claims": len(records), "issues": review(records)}


class SourceLedger:
    """Stateful wrapper: record() claims, review() the whole ledger."""

    def __init__(self) -> None:
        self._records: list[tuple[str, list[dict]]] = []

    def record(self, claim: str, sources: list[dict]) -> None:
        self._records.append((claim, list(sources)))

    def review(self) -> list[dict]:
        return review(self._records)

    def claims(self) -> int:
        return len(self._records)
