"""Generic keyword-namespace extractor for quantum-chemistry program manuals.

Different programs use different input-block conventions:

  ORCA       — ``%blockname  keyword value  end``  (in `.inp` files)
  Turbomole  — ``$blockname \\n  keyword value \\n $end`` (control file)
  Gaussian   — ``# keywords`` route + `--Link1--` sections
  NWChem     — ``BLOCKNAME ...  END`` blocks
  Q-Chem     — ``$blockname ... $end`` (similar to Turbomole)
  Psi4       — ``set BLOCK { key value }`` blocks
  ADF/AMS    — ``Engine BLOCK ... End`` blocks

This module abstracts the program-specific bits into a ``ProgramConfig``
dataclass so the same extraction pipeline works for any quantum-chemistry
suite.  Preset configs ship for ORCA + Turbomole; adding a new program
is one new ``ProgramConfig`` value.

The output is a structured keyword namespace per block.  Section
identifiers are sanitized to bare chapter-number prefixes — no
copyrightable manual prose is persisted.

License-safety: the extractor reads from the user-local
``~/.delfin/doc_index.json`` and persists ONLY:
  - block markers (e.g. ``%casscf``, ``$dft``) — public syntax
  - keyword names (e.g. ``nel``, ``norb``) — factual terminology
  - bare chapter numbers (e.g. ``ch3_13``) — structural reference
No verbatim manual sentences, no section titles, no equations.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


_DEFAULT_INDEX = Path.home() / ".delfin" / "doc_index.json"


# Words that LOOK like keywords but are actually prose noise.
_DEFAULT_STOPLIST = frozenset({
    "and", "or", "the", "for", "with", "is", "are", "this", "that", "to",
    "of", "in", "by", "on", "as", "if", "not", "be", "we", "you", "it",
    "see", "note", "also", "such", "any", "all", "one", "two", "three",
    "value", "default", "true", "false", "yes", "no", "end", "begin",
    "section", "table", "figure", "example", "user", "manual", "block",
    "keyword", "input", "output", "file", "line", "page", "chapter",
    "fig", "tab", "eq", "ref", "doc",
})


@dataclass(frozen=True)
class ProgramConfig:
    """Per-program extraction config.

    The pipeline:
    1. Find documents whose id matches ``doc_id_pattern`` in the
       indexed doc-index.
    2. For each section in those docs, find ``block_open_re`` matches.
    3. For each match, scan the following text up to ``block_close_re``
       (or ``max_scan_chars``) for keyword candidates.
    4. A candidate is the first ``[A-Za-z][A-Za-z0-9_]+`` token on each
       non-empty, non-comment line within the block window.
    5. Filter via ``stoplist`` + minimum-length rule.
    """

    program: str                                # short id ("orca", "turbomole")
    doc_id_pattern: str                         # regex on doc_id
    block_open_re: str                          # regex with one group = block name
    block_close_re: str                         # regex marking block end
    known_blocks: dict[str, tuple[str, ...]]    # canonical_name → aliases (raw text)
    max_scan_chars: int = 1500
    min_keyword_len: int = 2
    stoplist: frozenset[str] = field(default_factory=lambda: _DEFAULT_STOPLIST)
    comment_starts: tuple[str, ...] = ("#", "//", "!", "*")
    # Section-ID sanitizer: regex that yields the safe prefix to keep.
    # Default = ``ch<digits>(_<digits>)?`` (chapter-section number only).
    chapter_re: str = r"^(ch\d+(?:_\d+)?)"


# ---------------------------------------------------------------------------
# Built-in presets
# ---------------------------------------------------------------------------


ORCA_CONFIG = ProgramConfig(
    program="orca",
    doc_id_pattern=r"orca",
    block_open_re=r"%(\w+)",
    block_close_re=r"\bend\b",
    known_blocks={
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
    },
)


TURBOMOLE_CONFIG = ProgramConfig(
    program="turbomole",
    doc_id_pattern=r"turbomole",
    block_open_re=r"\$(\w+)",
    block_close_re=r"\$end\b",
    known_blocks={
        # Curate as Turbomole manual gets indexed.  Common blocks:
        "dft":         ("$dft",),
        "ricc2":       ("$ricc2",),
        "rij":         ("$rij",),
        "marij":       ("$marij",),
        "ridft":       ("$ridft",),
        "scfconv":     ("$scfconv",),
        "scfdamp":     ("$scfdamp",),
        "drvopt":      ("$drvopt",),
        "soes":        ("$soes",),
        "ex_opt":      ("$ex_opt",),
        "embed":       ("$embed",),
        "cosmo":       ("$cosmo",),
        "freeze":      ("$freeze",),
    },
)


GAUSSIAN_CONFIG = ProgramConfig(
    program="gaussian",
    doc_id_pattern=r"gaussian|g16|g09",
    # Gaussian route options aren't block-delimited but inline on `#`-line;
    # the simple block pattern catches `--Link1--` sections + `# keywords`.
    block_open_re=r"#\s*(\w+)",
    block_close_re=r"--Link1--",
    known_blocks={
        # To be filled when a Gaussian manual is indexed.
    },
)


_KNOWN_PROGRAMS: dict[str, ProgramConfig] = {
    "orca":      ORCA_CONFIG,
    "turbomole": TURBOMOLE_CONFIG,
    "gaussian":  GAUSSIAN_CONFIG,
}


# ---------------------------------------------------------------------------
# Extraction pipeline (generic over ProgramConfig)
# ---------------------------------------------------------------------------


def _load_index(path: Path | None = None) -> dict[str, Any]:
    p = path or _DEFAULT_INDEX
    with p.open(encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list) and data:
        data = data[0]
    return data


def _matching_doc_sections(
    index: dict[str, Any], doc_id_pattern: str,
) -> dict[str, dict[str, Any]]:
    """Return the union of sections from all docs whose ID matches."""
    docs = index.get("documents", {})
    rx = re.compile(doc_id_pattern, re.IGNORECASE)
    out: dict[str, dict[str, Any]] = {}
    for doc_id, doc in docs.items():
        if rx.search(doc_id):
            out.update(doc.get("sections", {}))
    return out


def _extract_block_keywords(
    block_aliases: tuple[str, ...],
    sections: dict[str, dict[str, Any]],
    config: ProgramConfig,
) -> tuple[set[str], list[str]]:
    """Return ``(keywords, section_ids)`` for a single block."""
    keywords: set[str] = set()
    section_ids: list[str] = []
    alias_re = re.compile(
        "|".join(re.escape(a) for a in block_aliases),
        re.IGNORECASE,
    )
    close_re = re.compile(config.block_close_re, re.IGNORECASE)
    for sid, sec in sections.items():
        text = sec.get("text", "") or ""
        if not text or not alias_re.search(text):
            continue
        section_ids.append(sid)
        for m in alias_re.finditer(text):
            start = m.end()
            tail = text[start:start + config.max_scan_chars]
            end_m = close_re.search(tail)
            if end_m:
                tail = tail[:end_m.start()]
            for line in tail.splitlines():
                line = line.strip()
                if not line or line.startswith(config.comment_starts):
                    continue
                first = re.match(r"([A-Za-z][A-Za-z0-9_]{1,30})", line)
                if not first:
                    continue
                kw = first.group(1).lower()
                if (kw in config.stoplist
                        or len(kw) < config.min_keyword_len
                        or kw.isdigit()
                        or kw.startswith("end")):
                    continue
                keywords.add(kw)
    return keywords, section_ids


def _sanitize_section_id(sid: str, config: ProgramConfig) -> str:
    """Strip title-slug suffixes to bare chapter prefix."""
    m = re.match(config.chapter_re, sid, re.IGNORECASE)
    return m.group(1) if m else sid.split("_")[0]


def extract_namespace(
    config: ProgramConfig,
    *,
    index_path: Path | None = None,
    sanitize_section_ids: bool = True,
) -> dict[str, dict[str, Any]]:
    """Extract a structured keyword namespace for the given program."""
    index = _load_index(index_path)
    sections = _matching_doc_sections(index, config.doc_id_pattern)
    out: dict[str, dict[str, Any]] = {}
    for block_name, aliases in config.known_blocks.items():
        keywords, section_ids = _extract_block_keywords(
            aliases, sections, config,
        )
        if not keywords and not section_ids:
            continue
        if sanitize_section_ids:
            refs = sorted({
                _sanitize_section_id(s, config) for s in section_ids
            })
        else:
            refs = sorted(set(section_ids))
        out[block_name] = {
            "block": aliases[0],
            "aliases": list(aliases),
            "keywords": sorted(keywords),
            "chapter_refs": refs[:30],
            "n_keywords": len(keywords),
            "n_chapter_refs": len(refs),
        }
    return out


def is_real_keyword(
    program: str,
    block: str,
    keyword: str,
    namespace: dict[str, dict[str, Any]] | None = None,
    *,
    path: Path | None = None,
) -> bool:
    """True if ``keyword`` is in the manual's namespace for
    ``program`` / ``block``.  Case-insensitive."""
    if namespace is None:
        cfg = _KNOWN_PROGRAMS.get(program.lower())
        if cfg is None:
            return False
        namespace = extract_namespace(cfg, index_path=path)
    info = namespace.get(block.lower())
    if not info:
        return False
    return keyword.lower() in {k.lower() for k in info.get("keywords", [])}


def get_config(program: str) -> ProgramConfig | None:
    return _KNOWN_PROGRAMS.get(program.lower())


__all__ = [
    "ProgramConfig",
    "ORCA_CONFIG",
    "TURBOMOLE_CONFIG",
    "GAUSSIAN_CONFIG",
    "extract_namespace",
    "is_real_keyword",
    "get_config",
]
