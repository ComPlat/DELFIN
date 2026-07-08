"""Runtime verify-enforcement guard against hallucinated ORCA keywords.

Phase 3 shipped the *documentation* and *benchmark* side of anti-
hallucination (Pattern 5 in the agent docs, the ``fact_verify`` task
class).  This module is the missing *runtime* piece: it scans the
agent's final answer for keyword claims that are NOT backed by the
ORCA manual ground-truth, so the dashboard can warn the user and — if
the turn never grounded itself via a doc-search — force exactly one
self-correction turn.

Design principle (same as ``generate_fact_tasks._validate_against_manual``):
the ground-truth is the extracted manual namespace, never author memory.
Both detectors are data-driven:

* **fake-keyword detector** (high precision): the union of every
  ``forbid`` entry in ``generate_fact_tasks._PROGRAM_BLOCK_TESTS`` that
  is confirmed absent from the manual namespace.  These are observed
  production hallucinations (``nactel``, ``nactorb``, ``multiplicity``
  used as a %casscf keyword, …).
* **unknown-keyword detector** (conservative): tokens written in the
  canonical ORCA block form ``keyword = value`` — and ONLY when the text
  actually shows ORCA input syntax (a ``%block`` marker) — that are not
  in the 1600+ keyword namespace and are not ordinary words.  It fires
  only when the answer is genuinely about ORCA keywords: DELFIN is a
  multi-tool agent, so backtick spans routinely quote CLI flags, file
  names, xTB methods (``gfn2``, ``gfnff``) and MANTA configs
  (``champion``, ``builder``) that are NOT ORCA keywords — those must
  never be judged against the ORCA namespace.

Nothing here touches the network or RDKit; it only reads the committed
``keywords_groundtruth_orca.json`` and is fully unit-testable.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


_HERE = Path(__file__).resolve().parent
_GROUNDTRUTH_PATH = _HERE / "pack" / "benchmark" / "keywords_groundtruth_orca.json"

# Tokens that pass the keyword-shape filter but are ordinary words a
# correct answer legitimately uses — never flag these as "unknown
# keywords" even when they appear backtick-quoted.
_KEYWORD_STOPWORDS = frozenset({
    "block", "blocks", "keyword", "keywords", "input", "output", "value",
    "values", "true", "false", "none", "default", "auto", "method",
    "methods", "option", "options", "orca", "manual", "section", "example",
    "string", "integer", "float", "bool", "list", "type", "name", "file",
    "calc", "calculation", "energy", "geometry", "basis", "functional",
})

# Forbid-list entries that are also ordinary words.  An answer may use
# these legitimately in prose ("the multiplicity is 3", "excited
# states"), so we only flag them as fake keywords when the answer
# *presents them as a keyword* (backtick-quoted or ``key = value``),
# never bare in prose.
_AMBIGUOUS_FAKES = frozenset({
    "multiplicity", "states", "density", "nuclei", "shifts", "tensor",
    "restart", "scaling", "increment", "guess", "level", "shift", "trust",
})


@dataclass(frozen=True)
class VerifyFlag:
    """One unverified keyword claim found in an answer."""

    keyword: str
    reason: str        # machine tag: "fake" | "unknown"
    suggestion: str    # user-facing nudge

    def message(self) -> str:
        return f"⚠️ Verify: '{self.keyword}' {self.suggestion}"


@lru_cache(maxsize=1)
def load_orca_namespace() -> frozenset[str]:
    """Return the set of all real ORCA keywords (lower-cased) across
    every block in the committed ground-truth.  Empty set if the file
    is missing — callers then degrade to the fake-keyword detector
    only (which has its own hard-coded list)."""
    if not _GROUNDTRUTH_PATH.exists():
        return frozenset()
    try:
        data = json.loads(_GROUNDTRUTH_PATH.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return frozenset()
    out: set[str] = set()
    for info in (data.get("blocks") or {}).values():
        for kw in info.get("keywords", []):
            out.add(str(kw).lower())
    return frozenset(out)


@lru_cache(maxsize=1)
def known_fake_keywords() -> frozenset[str]:
    """Union of every ``forbid`` entry across the curated block tests,
    filtered to those CONFIRMED absent from the manual namespace.

    Mirrors ``_validate_against_manual``: we never flag a token that is
    actually a real keyword, even if some curated forbid-list named it.
    """
    try:
        from delfin.agent.generate_fact_tasks import _PROGRAM_BLOCK_TESTS
    except Exception:
        return frozenset()
    namespace = load_orca_namespace()
    fakes: set[str] = set()
    for blocks in _PROGRAM_BLOCK_TESTS.values():
        for cfg in blocks.values():
            for kw in cfg.get("forbid", []) or []:
                low = str(kw).lower()
                # If the manual actually has it, it's not a fake.
                if namespace and low in namespace:
                    continue
                fakes.add(low)
    return frozenset(fakes)


# A token that "looks like an ORCA keyword": lower/mixed-case alnum,
# 3-24 chars, no spaces.  Used only inside keyword-presenting contexts.
_KEYWORD_SHAPE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{2,23}$")

# Inline-code spans: `nactel`, ``foo``.
_BACKTICK_SPAN = re.compile(r"`+([^`\n]{1,60})`+")

# "keyword = value" assignments, the canonical ORCA-block form.
_ASSIGN = re.compile(r"\b([A-Za-z][A-Za-z0-9_]{2,23})\s*=", re.MULTILINE)


def _word_present(text: str, word: str) -> bool:
    return re.search(rf"(?i)\b{re.escape(word)}\b", text) is not None


def scan_for_unverified_keywords(
    text: str,
    *,
    detect_unknown: bool = True,
) -> list[VerifyFlag]:
    """Scan ``text`` for keyword claims not backed by the ORCA manual.

    Returns a de-duplicated, order-stable list of :class:`VerifyFlag`.

    ``detect_unknown=False`` restricts the scan to the high-precision
    fake-keyword blocklist (skips the contextual unknown-keyword pass).
    """
    if not text or not text.strip():
        return []

    flags: list[VerifyFlag] = []
    seen: set[str] = set()

    # Tokens the answer *presents as* keywords: backtick-quoted spans AND
    # ``key = value`` assignments gate the ambiguous-fake detector.  The
    # unknown-keyword detector, by contrast, only considers ``assigned``
    # tokens (canonical ORCA block form) — a backtick span alone is NOT
    # evidence of an ORCA keyword claim (it quotes CLI flags, file names,
    # xTB methods, MANTA configs, ...).
    presented: list[str] = []
    assigned: list[str] = []
    for m in _BACKTICK_SPAN.finditer(text):
        presented.append(m.group(1).strip())
    for m in _ASSIGN.finditer(text):
        tok = m.group(1).strip()
        presented.append(tok)
        assigned.append(tok)
    presented_low = {t.lower() for t in presented}

    # 1. Fake-keyword blocklist — high precision, always on.  Ordinary-
    #    word fakes (_AMBIGUOUS_FAKES) only fire in keyword context.
    for fake in sorted(known_fake_keywords()):
        if fake in seen:
            continue
        if fake in _AMBIGUOUS_FAKES:
            hit = fake in presented_low
        else:
            hit = _word_present(text, fake)
        if hit:
            flags.append(VerifyFlag(
                keyword=fake,
                reason="fake",
                suggestion=("is not a real ORCA keyword (not in the manual) "
                            "— verify via search_docs and replace it."),
            ))
            seen.add(fake)

    if not detect_unknown:
        return flags

    namespace = load_orca_namespace()
    if not namespace:
        # No ground-truth loaded: can't tell unknown from known, so we
        # stay silent rather than risk false positives.
        return flags

    # 2. Unknown-keyword detector — fire ONLY in genuine ORCA-input
    #    context.  "About ORCA" means the text actually shows ORCA block
    #    syntax: a ``%block`` marker (``%scf``, ``%casscf``, ...), NOT the
    #    bare word "ORCA" (a downstream program mention) and NOT a plain
    #    "%" (e.g. "50 %").  Only ``keyword = value`` tokens are judged —
    #    never arbitrary backtick spans.
    has_orca_block = re.search(r"%[A-Za-z]", text) is not None
    if not has_orca_block:
        return flags

    for tok in assigned:
        low = tok.lower()
        if low in seen:
            continue
        if not _KEYWORD_SHAPE.match(tok):
            continue
        if low in _KEYWORD_STOPWORDS:
            continue
        if low in namespace:
            continue
        flags.append(VerifyFlag(
            keyword=tok,
            reason="unknown",
            suggestion=("is not in the ORCA manual namespace — verify via "
                        "search_docs before presenting it as a keyword."),
        ))
        seen.add(low)

    return flags


def correction_feedback(flags: list[VerifyFlag]) -> str:
    """Build the feedback message for the forced self-correction turn."""
    kws = ", ".join(f"'{f.keyword}'" for f in flags)
    return (
        f"The following keywords are not backed by the ORCA manual: {kws}. "
        "Look them up via search_docs and correct or remove them. Only "
        "mention keywords that actually exist in the manual."
    )
