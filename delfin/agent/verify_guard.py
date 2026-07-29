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
from typing import Optional


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


# ---------------------------------------------------------------------------
# Code-claim grounding — verify statements about code against reality
# ---------------------------------------------------------------------------
#
# The honesty addendum demands file:line citations, and the framework checks
# test claims (evidence ledger) and ORCA keywords (above) mechanically — but
# nothing checked claims ABOUT CODE. This closes that gap: citations in a
# final answer are cross-checked against (a) the filesystem and (b) the set
# of files the agent actually observed (read/grepped) this session.

# Extensions that make a dotted token a plausible FILE citation. Guards the
# path regex against dotted module names ("delfin.agent.memory_store") and
# version numbers, which must never be flagged.
_CODE_FILE_EXTS = frozenset({
    "py", "md", "json", "yaml", "yml", "toml", "cfg", "ini", "txt", "sh",
    "js", "ts", "html", "css", "csv", "xyz", "inp", "out", "mol", "hess",
    "gbw", "ipynb", "rst", "lock", "sql", "c", "h", "cpp", "hpp", "rs",
})


@dataclass(frozen=True)
class CodeClaimFlag:
    """One ungrounded code citation found in an answer."""

    path: str
    line: Optional[int]
    kind: str          # "nonexistent" | "unread"

    def message(self) -> str:
        ref = f"{self.path}:{self.line}" if self.line else self.path
        if self.kind == "nonexistent":
            return (f"⚠️ Verify: '{ref}' is cited but does not exist in the "
                    f"workspace — never invent file paths.")
        return (f"⚠️ Verify: '{ref}' is cited but was not read or grepped "
                f"this session — open it before describing it.")


def _iter_path_citations(text: str):
    """Yield (path, line|None) citations that plausibly reference files."""
    from delfin.agent.memory_store import _PATH_REF_RE
    for m in _PATH_REF_RE.finditer(text or ""):
        raw = m.group(1).strip().rstrip(".,;:")
        line = int(m.group(2)) if m.group(2) else None
        ext = raw.rsplit(".", 1)[-1].lower() if "." in raw else ""
        if "/" not in raw and ext not in _CODE_FILE_EXTS:
            continue
        if "/" in raw and "." in raw.rsplit("/", 1)[-1]:
            tail_ext = raw.rsplit(".", 1)[-1].lower()
            if tail_ext not in _CODE_FILE_EXTS:
                continue
        if raw.startswith(("http", "www.")) or "://" in raw:
            continue
        yield raw, line


def _is_observed(path: str, observed: frozenset[str]) -> bool:
    """True when a cited path matches any observed file (suffix-tolerant:
    observed entries may be absolute while the answer cites repo-relative,
    or vice versa)."""
    if not observed:
        return False
    norm = path.lstrip("./")
    for o in observed:
        on = str(o).replace("\\", "/").lstrip("./")
        if on == norm or on.endswith("/" + norm) or norm.endswith("/" + on):
            return True
    return False


def scan_for_ungrounded_code_claims(
    text: str,
    *,
    repo_root: Path | str | None = None,
    observed_files: Optional[frozenset[str] | set[str]] = None,
    max_flags: int = 8,
) -> list[CodeClaimFlag]:
    """Cross-check file citations in ``text`` against reality.

    - a cited path that exists nowhere in the workspace  -> "nonexistent"
      (hard flag: drives the forced self-correction turn)
    - a cited path that exists but was never observed via read/grep this
      session -> "unread" (soft flag only: injected context such as the
      repo map legitimately surfaces paths the model may mention)

    Deterministic, order-stable, capped at ``max_flags``. Never raises.
    """
    flags: list[CodeClaimFlag] = []
    seen: set[str] = set()
    obs = frozenset(str(p) for p in (observed_files or ()))
    root = Path(repo_root) if repo_root else None
    try:
        for path, line in _iter_path_citations(text):
            key = f"{path}:{line}"
            if key in seen or len(flags) >= max_flags:
                continue
            seen.add(key)
            if _is_observed(path, obs):
                continue
            exists = False
            try:
                p = Path(path)
                if p.is_absolute():
                    exists = p.is_file()
                elif root is not None:
                    exists = (root / path).is_file()
                else:
                    exists = Path(path).is_file()
            except OSError:
                exists = False
            if not exists:
                flags.append(CodeClaimFlag(path=path, line=line,
                                           kind="nonexistent"))
            else:
                flags.append(CodeClaimFlag(path=path, line=line,
                                           kind="unread"))
    except Exception:
        return flags
    return flags


def code_claim_feedback(flags: list[CodeClaimFlag]) -> str:
    """Feedback message for the forced self-correction turn (nonexistent
    citations only — 'unread' stays a soft warning)."""
    bad = [f"'{f.path}:{f.line}'" if f.line else f"'{f.path}'"
           for f in flags if f.kind == "nonexistent"]
    return (
        f"The following cited paths do not exist in the workspace: "
        f"{', '.join(bad)}. Read or grep the actual files and correct the "
        "answer — cite only paths you have verified."
    )


# ---------------------------------------------------------------------------
# Physical-quantity grounding — numbers with units need an evidence act
# ---------------------------------------------------------------------------
#
# The citation ledger above covers claims about CODE; nothing covered a
# chemistry agent's most consequential claim type: PHYSICAL QUANTITIES.
# Without this scanner an answer could state a numeric result ("the S1
# energy is 2.31 eV") in a turn that never read any calculation output and
# never performed a docs/calc lookup. This scanner detects unit-anchored
# numeric claims and flags them when the turn shows NO evidence act at all.
#
# Scope: this is a TURN-LEVEL gate, deliberately. Attributing individual
# numbers to individual files (did 2.31 eV really come from tddft.out?) is
# out of scope — once at least one evidence act exists this turn (any
# calculation-output-like file observed, or any lookup tool used), nothing
# is flagged. The gate only catches the zero-evidence case, which is the
# unambiguous hallucination signature.

# File extensions that mark an observed file as calculation output —
# reading one counts as an evidence act for numeric claims.
_CALC_OUTPUT_EXTS = frozenset({
    "out", "log", "json", "xyz", "hess", "gbw",
})

# Tool names whose use counts as an evidence act (docs / calc / file /
# history lookups). Compared against bare names; mcp__ns__name forms are
# reduced to their tail before comparison.
_EVIDENCE_TOOLS = frozenset({
    "search_docs", "read_section", "search_calcs", "get_calc",
    "read_file", "grep_file", "history_search", "web_fetch",
})

# Fallback shape match so lookup tools count as evidence regardless of the
# backend's naming scheme (CLI backends report e.g. "Read"/"Grep"/"WebFetch"
# instead of the snake_case names above). Editing/exec tools never match.
_EVIDENCE_TOOL_SHAPE = re.compile(r"(?i)(read|grep|search|fetch|docs|glob)")

# Number token for unit-anchored claims. The lookbehind stops mid-token
# matches (dotted version strings like 6.0.1 can never contribute their
# tail digits); the sign class includes the Unicode minus.
_QTY_NUM = r"(?<![\w.])[-+−]?\d+(?:\.\d+)?"

# Unit-anchored claim patterns: a number IMMEDIATELY before the unit
# (at most one whitespace char between them). Percentages and bare
# version numbers never match because every pattern requires a unit
# token. Order is significant only for readability — matches are sorted
# by text position afterwards; bare 'kcal' excludes 'kcal/mol' via a
# lookahead so the two patterns never double-flag one claim.
_QUANTITY_PATTERNS: tuple[tuple[str, "re.Pattern[str]"], ...] = tuple(
    (unit, re.compile(pattern, re.MULTILINE))
    for unit, pattern in (
        ("eV",       rf"{_QTY_NUM}\s?eV\b"),
        ("kcal/mol", rf"{_QTY_NUM}\s?kcal\s?/\s?mol\b"),
        ("kJ/mol",   rf"{_QTY_NUM}\s?kJ\s?/\s?mol\b"),
        ("Hartree",  rf"{_QTY_NUM}\s?(?:[Hh]artrees?|Eh)\b"),
        ("nm",       rf"{_QTY_NUM}\s?nm\b"),
        ("cm-1",     rf"{_QTY_NUM}\s?cm(?:\^-1|-1|⁻¹|−1)(?!\d)"),
        ("Debye",    rf"{_QTY_NUM}\s?[Dd]ebye\b"),
        ("Angstrom", rf"{_QTY_NUM}\s?(?:[Åå](?:ngstr(?:ö|o)ms?)?"
                     rf"|[Aa]ngstr(?:ö|o)ms?)\b"),
        ("ps",       rf"{_QTY_NUM}\s?ps\b"),
        ("fs",       rf"{_QTY_NUM}\s?fs\b"),
        # Kelvin: only the standalone-uppercase-K form directly after a
        # number, terminated by punctuation/whitespace/end — the loosest
        # anchor in the set, so it gets the strictest boundary.
        ("K",        rf"{_QTY_NUM}\s?K(?=$|[\s.,;:)\]!?])"),
        ("kcal",     rf"{_QTY_NUM}\s?kcal\b(?!\s?/)"),
        ("GHz",      rf"{_QTY_NUM}\s?GHz\b"),
    )
)

# Fenced code blocks (``` ... ```): never claims.
_FENCED_BLOCK = re.compile(r"```.*?```", re.DOTALL)


@dataclass(frozen=True)
class QuantityClaimFlag:
    """One physical-quantity claim stated without any evidence act."""

    quantity: str      # matched text, whitespace-normalized ("2.31 eV")
    unit: str          # canonical unit tag ("eV", "kcal/mol", ...)

    def message(self) -> str:
        return (f"⚠️ Verify: quantity '{self.quantity}' stated without any "
                f"calculation output or docs lookup this turn — state the "
                f"source or verify first.")


def _strip_non_claim_regions(text: str) -> str:
    """Blank out regions that are not the agent's own claims: fenced code
    blocks, inline backtick spans, and blockquoted lines ('> ...')."""
    t = _FENCED_BLOCK.sub(" ", text)
    t = _BACKTICK_SPAN.sub(" ", t)
    return "\n".join(
        "" if ln.lstrip().startswith(">") else ln
        for ln in t.split("\n")
    )


def _observed_calc_output(observed) -> bool:
    """True when any observed file looks like calculation output (.out,
    .log, .json, .xyz, .hess, .gbw, or an ORCA property file)."""
    for p in observed or ():
        name = str(p).replace("\\", "/").rsplit("/", 1)[-1].lower()
        ext = name.rsplit(".", 1)[-1] if "." in name else ""
        if ext in _CALC_OUTPUT_EXTS or "property" in name:
            return True
    return False


def _used_evidence_tool(tools) -> bool:
    """True when any tool name (bare or mcp__ns__name) is a lookup tool."""
    for t in tools or ():
        name = str(t).strip().lower()
        if name in _EVIDENCE_TOOLS or name.rsplit("__", 1)[-1] in _EVIDENCE_TOOLS:
            return True
        if _EVIDENCE_TOOL_SHAPE.search(name):
            return True
    return False


def scan_for_unsourced_quantities(
    text: str,
    *,
    observed_files: Optional[frozenset[str] | set[str]] = None,
    evidence_tools_used: Optional[frozenset[str] | set[str]] = None,
    max_flags: int = 6,
) -> list[QuantityClaimFlag]:
    """Scan ``text`` for physical-quantity claims made in a turn with no
    evidence act.

    A claim is a number immediately followed by a physical unit (eV,
    kcal/mol, kJ/mol, Hartree/Eh, nm, cm-1, Debye, Angstrom/Å, ps, fs,
    K, kcal, GHz). Backticked spans, fenced code blocks, blockquoted
    lines, percentages and version strings never count.

    Turn-level gate (documented design decision): when ``observed_files``
    contains any calculation-output-like file OR ``evidence_tools_used``
    contains any lookup tool, the turn had an evidence act and NOTHING is
    flagged — per-number attribution to per-file sources is out of scope.
    Only the zero-evidence turn is flagged.

    Deterministic, order-stable (text position), de-duplicated, capped at
    ``max_flags``. Never raises.
    """
    flags: list[QuantityClaimFlag] = []
    try:
        if not text or not text.strip() or max_flags <= 0:
            return []
        if _observed_calc_output(observed_files):
            return []
        if _used_evidence_tool(evidence_tools_used):
            return []
        scrubbed = _strip_non_claim_regions(text)
        matches: list[tuple[int, str, str]] = []
        for unit, rx in _QUANTITY_PATTERNS:
            for m in rx.finditer(scrubbed):
                if _is_hedged(scrubbed, m.start(), m.end()):
                    continue
                qty = " ".join(m.group(0).split())
                matches.append((m.start(), qty, unit))
        matches.sort(key=lambda t: (t[0], t[2]))
        seen: set[str] = set()
        for _pos, qty, unit in matches:
            key = qty.lower()
            if key in seen:
                continue
            seen.add(key)
            flags.append(QuantityClaimFlag(quantity=qty, unit=unit))
            if len(flags) >= max_flags:
                break
    except Exception:
        return flags
    return flags


def quantity_claim_feedback(flags: list[QuantityClaimFlag]) -> str:
    """Feedback message for the forced self-correction turn."""
    qtys = ", ".join(f"'{f.quantity}'" for f in flags)
    return (
        f"The following physical quantities were stated without any "
        f"evidence act this turn: {qtys}. Either read the source now "
        "(calculation output via get_calc/read_file, or the docs via "
        "search_docs) and cite it next to each value, or restate each "
        "value as unverified and say what would confirm it."
    )


# ---------------------------------------------------------------------------
# Hedge detection — explicitly uncertain claims are never enforcement targets
# ---------------------------------------------------------------------------
#
# Every claim scanner below (and the quantity scanner above) skips claims
# whose containing line carries an explicit uncertainty marker: enforcement
# exists to stop CONFIDENT fabrication, and an answer that says it is
# guessing has already disclosed its grounding status.

_HEDGE_MARKERS = re.compile(
    r"(?i)\b(?:"
    # German
    r"vermutlich|wahrscheinlich|m(?:ö|oe)glicherweise|vielleicht|"
    r"ungef(?:ä|ae)hr|etwa|circa|ca\.|sch(?:ä|ae)tzungsweise|grob|"
    r"unsicher|ungepr(?:ü|ue)ft|unverifiziert|"
    r"nicht\s+(?:gepr(?:ü|ue)ft|verifiziert|(?:ü|ue)berpr(?:ü|ue)ft|sicher|"
    r"nachgesehen|nachgeschaut)|"
    r"k(?:ö|oe)nnte|d(?:ü|ue)rfte|scheint|offenbar|"
    r"aus\s+dem\s+(?:ged(?:ä|ae)chtnis|kopf)|"
    # English
    r"probably|likely|perhaps|maybe|possibly|approximately|roughly|around|"
    r"unsure|uncertain|unverified|presumably|"
    r"not\s+(?:sure|certain|verified|checked|confirmed)|"
    r"i\s+(?:think|believe|guess|assume|recall|suspect|expect)|"
    r"have\s*n(?:o|')t\s+(?:checked|verified|looked|read|confirmed)|"
    r"did\s*n(?:o|')t\s+(?:check|verify|look|read|confirm)|"
    r"without\s+(?:checking|looking|verifying)|"
    r"from\s+memory|iirc|if\s+i\s+recall|"
    r"might|may\s+be|should\s+be|seems?|estimated?"
    r")\b"
)


def _is_hedged(text: str, start: int, end: int) -> bool:
    """True when the line containing ``text[start:end]`` carries an explicit
    uncertainty marker. Line-scoped on purpose: a hedge in one sentence must
    not blanket-immunize confident claims elsewhere in the answer."""
    lo = text.rfind("\n", 0, start) + 1
    hi = text.find("\n", end)
    if hi == -1:
        hi = len(text)
    return _HEDGE_MARKERS.search(text[lo:hi]) is not None


# ---------------------------------------------------------------------------
# Code-location claims — specific file/line assertions need observed evidence
# ---------------------------------------------------------------------------
#
# The citation scanner above checks that cited PATHS exist / were read, but a
# fabricated answer can assert a specific LOCATION without ever naming a
# path in citable form ("Zeile 26: class AgentEngine") — no scanner caught
# that. This one detects three claim shapes:
#
#   file_line  — path plus a concrete line number ("engine.py:281",
#                "engine.py, Zeile 281", "line 281 of engine.py")
#   bare_line  — a concrete line number next to a code symbol, no path
#                ("Zeile 26: class AgentEngine")
#   defined_in — a symbol-is-defined-in-file assertion without a line
#                ("AgentEngine is defined in delfin/agent/engine.py")
#
# Grounding rules (evidence = the observed-files ledger, i.e. files the
# model demonstrably read/grepped this session):
#
#   file_line  — ungrounded unless THAT path is in the ledger. Injected
#                context (repo map) surfaces paths but never line numbers,
#                so a line-number claim about an unread file is the
#                fabrication signature.
#   bare_line / defined_in — no attributable path, so these use the
#                turn-level rule: ungrounded only when the ledger is EMPTY
#                (zero files observed — the unambiguous zero-evidence case).
#
# When the backend exposes no ledger at all (``ledger_available=False``)
# the scanner stays silent: absence of bookkeeping is not evidence of
# fabrication. Hedged claims (see ``_is_hedged``) are never flagged, and
# statements without a concrete number ("near the top of the file") never
# match. No task-specific patterns, no filename allowlists.

# Path token with a code-file extension. Reuses the citation extension set
# so dotted module names and version strings never match.
_LOC_PATH = (
    r"(?:[\w\-][\w.\-]*/)*[\w\-][\w.\-]*"
    r"\.(?:" + "|".join(sorted(_CODE_FILE_EXTS)) + r")\b"
)

# Line-reference word (English + German).
_LINE_WORD = r"(?:lines?|zeilen?|ln\.?)"

# Same-line gap that must not skip over another path token — keeps the
# claim attributed to the NEAREST path (grounding is judged per path).
_LOC_GAP = rf"(?:(?!{_LOC_PATH})[^\n]){{0,40}}?"

# file_line claim shapes. One path group + one line group each.
_LOC_FILE_LINE_PATTERNS: tuple["re.Pattern[str]", ...] = (
    # path:N
    re.compile(rf"({_LOC_PATH}):(\d{{1,6}})\b", re.IGNORECASE),
    # path ... line N (short same-line gap)
    re.compile(rf"({_LOC_PATH}){_LOC_GAP}\b{_LINE_WORD}\s*(\d{{1,6}})\b",
               re.IGNORECASE),
    # line N ... of/in path (short same-line gap)
    re.compile(rf"\b{_LINE_WORD}\s*(\d{{1,6}})\b{_LOC_GAP}"
               rf"\b(?:of|in|von|im|der|aus)\b{_LOC_GAP}({_LOC_PATH})",
               re.IGNORECASE),
)

# bare_line claim shape: a line-number reference with no path on the line.
_LOC_BARE_LINE = re.compile(rf"\b{_LINE_WORD}\s*(\d{{1,6}})\b", re.IGNORECASE)
_LOC_PATH_RE = re.compile(_LOC_PATH, re.IGNORECASE)

# defined_in claim shapes (no line number): definition verb near a path.
_LOC_DEFINED_PATTERNS: tuple["re.Pattern[str]", ...] = (
    re.compile(rf"\b(?:defined|declared|implemented|located|lives)\b"
               rf"[^\n]{{0,60}}?({_LOC_PATH})", re.IGNORECASE),
    re.compile(rf"({_LOC_PATH})[^\n]{{0,60}}?"
               rf"\b(?:definiert|deklariert|implementiert)\b", re.IGNORECASE),
)

# Code-symbol signals gating the bare_line shape (both regexes checked
# against the claim's line; keyword one is case-insensitive, shape one is
# case-sensitive so CamelCase detection stays meaningful).
_CODE_KEYWORD_RE = re.compile(
    r"(?i)\b(?:class|def|function|method|klasse|funktion|methode|"
    r"defin(?:ed|iert)|declared|deklariert)\b")
_CODE_SHAPE_RE = re.compile(
    r"\b[A-Z][A-Za-z0-9]*[a-z][A-Z][A-Za-z0-9]*\b"   # CamelCase
    r"|\b[a-z][a-z0-9]*_[a-z0-9_]+\b"                 # snake_case
    r"|\b\w+\(\)")                                    # call form


@dataclass(frozen=True)
class LocationClaimFlag:
    """One code-location claim without matching observed evidence."""

    claim: str          # matched text, whitespace-normalized, truncated
    path: str           # cited path ("" for bare_line)
    line: Optional[int]
    kind: str           # "file_line" | "bare_line" | "defined_in"

    def message(self) -> str:
        return (f"⚠️ Verify: location claim '{self.claim}' is not backed by "
                f"any file read or grep this session — open the source or "
                f"state the location as unverified.")


def _scrub_keep_offsets(text: str) -> str:
    """Blank non-claim regions (fenced code blocks, blockquoted lines) and
    backtick characters while PRESERVING offsets, so hedge lookups on match
    positions stay valid. Inline-backtick content is kept: backticked
    ``path:line`` spans are exactly the claims this scanner targets."""
    t = _FENCED_BLOCK.sub(lambda m: " " * (m.end() - m.start()), text)
    t = "\n".join(
        " " * len(ln) if ln.lstrip().startswith(">") else ln
        for ln in t.split("\n")
    )
    return t.replace("`", " ")


def _line_of(text: str, start: int, end: int) -> str:
    lo = text.rfind("\n", 0, start) + 1
    hi = text.find("\n", end)
    return text[lo:hi if hi != -1 else len(text)]


def scan_for_ungrounded_location_claims(
    text: str,
    *,
    observed_files: Optional[frozenset[str] | set[str]] = None,
    ledger_available: bool = True,
    max_flags: int = 6,
) -> list[LocationClaimFlag]:
    """Scan ``text`` for specific code-location claims lacking observed
    evidence. See the section comment for claim shapes and grounding rules.

    Deterministic, order-stable, de-duplicated, capped at ``max_flags``.
    Never raises.
    """
    flags: list[LocationClaimFlag] = []
    try:
        if not text or not text.strip() or max_flags <= 0:
            return []
        if not ledger_available:
            return []
        obs = frozenset(str(p) for p in (observed_files or ()))
        any_observed = bool(obs)
        scrubbed = _scrub_keep_offsets(text)
        seen: set[str] = set()
        consumed: list[tuple[int, int]] = []   # spans claimed by file_line
        claimed_paths: set[str] = set()

        def _add(claim: str, path: str, line: Optional[int],
                 kind: str) -> None:
            key = f"{kind}:{path.lower()}:{line}"
            if key in seen or len(flags) >= max_flags:
                return
            seen.add(key)
            flags.append(LocationClaimFlag(
                claim=" ".join(claim.split())[:80],
                path=path, line=line, kind=kind))

        # 1. file_line — path + concrete line number.
        for rx in _LOC_FILE_LINE_PATTERNS:
            for m in rx.finditer(scrubbed):
                g1, g2 = m.group(1), m.group(2)
                path, num = (g1, g2) if not g1.isdigit() else (g2, g1)
                path = path.strip().rstrip(".,;:")
                consumed.append((m.start(), m.end()))
                claimed_paths.add(path.lower())
                if _is_hedged(scrubbed, m.start(), m.end()):
                    continue
                if _is_observed(path, obs):
                    continue
                _add(m.group(0), path, int(num), "file_line")

        # 2. bare_line — line number + code symbol, no path involved.
        for m in _LOC_BARE_LINE.finditer(scrubbed):
            if any(s <= m.start() < e for s, e in consumed):
                continue
            line_text = _line_of(scrubbed, m.start(), m.end())
            if _LOC_PATH_RE.search(line_text):
                continue   # path on the same line: file_line territory
            if not (_CODE_KEYWORD_RE.search(line_text)
                    or _CODE_SHAPE_RE.search(line_text)):
                continue   # no code symbol: not a code-location claim
            if _is_hedged(scrubbed, m.start(), m.end()):
                continue
            if any_observed:
                continue   # turn-level rule: any observation grounds these
            _add(m.group(0), "", int(m.group(1)), "bare_line")

        # 3. defined_in — symbol-defined-in-file, no line number.
        for rx in _LOC_DEFINED_PATTERNS:
            for m in rx.finditer(scrubbed):
                path = m.group(1).strip().rstrip(".,;:")
                if path.lower() in claimed_paths:
                    continue   # already covered by a file_line claim
                if _is_hedged(scrubbed, m.start(), m.end()):
                    continue
                if any_observed or _is_observed(path, obs):
                    continue   # turn-level rule (see section comment)
                _add(m.group(0), path, None, "defined_in")
    except Exception:
        return flags
    return flags


def location_claim_feedback(flags: list[LocationClaimFlag]) -> str:
    """Feedback message for the forced self-correction turn."""
    refs = ", ".join(f"'{f.claim}'" for f in flags)
    return (
        f"The following code-location claims are not backed by any file "
        f"read or grep in this session: {refs}. Verify now — read or grep "
        "the referenced source and cite the confirmed file and line — or "
        "restate the answer marking the location as unverified."
    )


def grounding_caveat(
    location_flags: list[LocationClaimFlag],
    quantity_flags: list[QuantityClaimFlag],
) -> str:
    """Visible caveat appended when the single correction turn still lacks
    grounding — the turn is never failed, the reader is warned instead."""
    items: list[str] = []
    items.extend(f"'{f.claim}'" for f in location_flags)
    items.extend(f"'{f.quantity}'" for f in quantity_flags)
    if not items:
        return ""
    return (
        "\n\n[verify] Caveat: the following claims remain unverified — no "
        "file read or lookup in this session backs them: "
        + ", ".join(items) + ". Treat them as unconfirmed."
    )
