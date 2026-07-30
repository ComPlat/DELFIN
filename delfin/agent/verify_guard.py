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
            # Prose alternations ("ORCA/xTB", "APIs/Bibliotheken",
            # "Input/Output") match the path shape but are not file
            # citations. An extensionless slash token only counts when
            # its FIRST segment is a real directory here — fabricated
            # paths imitate existing structure; word pairs do not.
            tail = path.rsplit("/", 1)[-1]
            if "/" in path and "." not in tail and not line:
                head = path.lstrip("./").split("/", 1)[0]
                try:
                    if root is None or not (root / head).is_dir():
                        continue
                except OSError:
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


# ---------------------------------------------------------------------------
# Functional claims — "it works now" requires the artifact to have been run
# ---------------------------------------------------------------------------
#
# The scanners above judge claims about code LOCATIONS and about measured
# QUANTITIES. A third class slipped past both: the claim that something the
# session BUILT now works — "das Skript läuft fehlerfrei", "the app is ready
# to use", "beide Spiele funktionieren im Browser". Passing tests over a
# module's logic and a server that answers on its port say nothing about
# whether the delivered artifact behaves as described, and interactive or
# browser behavior cannot be observed headlessly at all.
#
# Two claim kinds, judged against what the session actually EXECUTED:
#
#   interactive — the claim is about browser / UI / keyboard / mouse /
#                 playability. No headless check available here can observe
#                 that, so it is flagged whenever asserted plainly — unless
#                 the session drove a real UI (browser-automation tool or
#                 command).
#   runtime     — the claim is that something runs / starts / is usable.
#                 Grounded when the artifact named in the claim appears in a
#                 command the session actually ran. Launching a server is
#                 explicitly NOT evidence: it shows the server starts, never
#                 that what it serves works — such commands are recorded as
#                 start-only. When the claim names no artifact, the
#                 turn-level rule applies (mirroring bare_line above): it is
#                 flagged only when the session ran nothing at all.
#
# Never flagged: hedged claims (``_is_hedged``), explicit non-verification
# disclosures ("I could not verify ..."), negated statements, conditionals
# ("damit es funktioniert ..."), questions, and statements attributed to the
# user's own report. Multilingual by construction (German + English), since
# the agent answers in the user's language.

# Execution tools whose calls are recorded as "this session ran X".
_EXEC_TOOLS = frozenset({
    "bash", "bash_background", "run_tests", "run_command", "run_script",
    "shell", "exec", "execute", "terminal", "pytest",
})

# Read-only companions of the execution tools: they inspect a job, they do
# not run anything, so they must never count as an execution act.
_EXEC_TOOL_DENY = frozenset({
    "bash_output", "bash_status", "bash_kill", "watch_job", "job_monitor",
    "run_status",
})

# Fallback shape match so execution tools count regardless of the backend's
# naming scheme (CLI backends report e.g. "Bash").
_EXEC_TOOL_SHAPE = re.compile(
    r"(?i)(?:^|_)(?:bash|shell|exec|run|pytest|command|terminal)(?:_|$)")

# Tool-input keys that carry the executed command / target.
_EXEC_INPUT_KEYS = ("command", "cmd", "script", "target", "pytest_args",
                    "args", "file", "path")

# Commands that START something without exercising it. A server launch (and
# any backgrounded start) is evidence that the process starts — never that
# the artifact it serves behaves as claimed.
_START_ONLY_CMD_RE = re.compile(
    r"(?i)\b(?:bash_background|voila|voil(?:à|a)|streamlit|uvicorn|gunicorn|"
    r"flask|http\.server|serve|runserver|nodemon|vite|webpack|"
    r"jupyter\s+(?:notebook|lab|server)|ng\s+serve|next\s+dev|"
    r"npm\s+(?:run\s+)?(?:start|dev|serve)|yarn\s+(?:start|dev|serve)|"
    r"php\s+-S|rails\s+s(?:erver)?)\b")

# Test-runner commands, plus the shape of a test file: a test run grounds
# claims about the test files it discovers even when the command names only
# a directory.
_TEST_CMD_RE = re.compile(
    r"(?i)\b(?:pytest|run_tests|unittest|nose2?|tox|jest|vitest|"
    r"npm\s+(?:run\s+)?test)\b")
_TEST_FILE_RE = re.compile(r"(?i)^(?:test_[\w\-.]+|[\w\-.]+_test)\.(?:py|js|ts)$")

# A tool or command that actually drives a user interface. If one was used,
# interactive claims have real evidence behind them and are not flagged.
_UI_EXERCISE_TOOL_SHAPE = re.compile(
    r"(?i)(playwright|selenium|puppeteer|webdriver|chromedriver|browser_|"
    r"_browser|screenshot|chrome|devtools|click_element)")
_UI_EXERCISE_CMD_RE = re.compile(
    r"(?i)\b(?:playwright|cypress|selenium|puppeteer|webdriver|chromedriver|"
    r"xvfb-run)\b")

# Playability / hands-on predicates: asserting these IS an interactive claim,
# no extra marker needed.
_FUNC_PLAY_PAT = (
    r"(?:"
    # German
    r"\bspielbar\b|\bbedienbar\b|\bsteuerbar\b|"
    # both word orders: "kannst du ... spielen" / "du kannst ... spielen"
    r"\b(?:(?:kannst|k(?:ö|oe)nnen|kann)\s+(?:du|sie|ihr|man)|"
    r"(?:du|sie|ihr|man)\s+(?:kannst|k(?:ö|oe)nnen|kann))\b[^\n]{0,48}?"
    r"\b(?:spielen|steuern|bedienen|klicken|dr(?:ü|ue)cken)\b|"
    r"\bzum\s+Spielen\s+bereit\b|"
    # English
    r"\bplayable\b|"
    r"\byou\s+can\s+(?:now\s+)?(?:play|control|click|press|drag|move)\b|"
    r"\bready\s+to\s+play\b"
    r")"
)
_FUNC_PLAY_RE = re.compile(_FUNC_PLAY_PAT, re.IGNORECASE)

# General "it works / runs / is usable" predicates.
_FUNC_WORK_PAT = (
    r"(?:"
    # German
    r"\bfunktionier(?:t|en)\b|\bfunktionsf(?:ä|ae)hig\b|"
    r"\bfunktionst(?:ü|ue)chtig\b|\beinsatzbereit\b|\bbetriebsbereit\b|"
    r"\blauff(?:ä|ae)hig\b|\bbenutzbar\b|\bnutzbar\b|\bverwendbar\b|"
    r"\b(?:l(?:ä|ae)uft|laufen)\s+(?:jetzt\s+|nun\s+|wieder\s+)?"
    r"(?:fehlerfrei|problemlos|stabil|einwandfrei|durch|erfolgreich|sauber|"
    r"korrekt|wie\s+erwartet|jetzt|nun|wieder)\b|"
    r"\bstartet\s+(?:jetzt\s+|nun\s+)?"
    r"(?:erfolgreich|problemlos|fehlerfrei|sauber|jetzt|nun)\b|"
    r"\b(?:(?:kannst|k(?:ö|oe)nnen|kann)\s+(?:du|sie|ihr|man)|"
    r"(?:du|sie|ihr|man)\s+(?:kannst|k(?:ö|oe)nnen|kann))\s+"
    r"(?:(?:es|ihn|sie|das|die|den|jetzt|nun|direkt|sofort)\s+){0,3}"
    r"(?:nutzen|benutzen|verwenden|starten|(?:ö|oe)ffnen|ausf(?:ü|ue)hren)\b|"
    # English
    r"\bworks\b|\bwork\s+(?:correctly|fine|properly|now|as\s+expected|"
    r"in\s+the\s+browser)\b|\b(?:is|are)\s+working\b|"
    r"\bruns\s+(?:fine|successfully|correctly|cleanly|smoothly|now|"
    r"without\s+(?:errors|issues|problems))\b|"
    r"\b(?:starts|launches)\s+(?:successfully|fine|cleanly|now|"
    r"without\s+errors)\b|"
    r"\b(?:is|are)\s+(?:now\s+)?(?:fully\s+)?"
    r"(?:functional|operational|usable|ready\s+to\s+(?:use|run)|"
    r"ready\s+for\s+use)\b|"
    r"\byou\s+can\s+(?:now\s+)?(?:use|run|open|start|try)\b"
    r")"
)
_FUNC_WORK_RE = re.compile(_FUNC_WORK_PAT, re.IGNORECASE)

# An absolute completeness claim about verification ("vollständig
# getestet", "alles geprüft", "fully tested") asserts the ABSENCE of
# untested parts. No amount of executed commands can establish that — a
# green test run says what it covered, never what it did not. Observed in
# the field: a package whose e-mail path was never exercised was handed
# over as "vollständig getestet" while a real test run made the ordinary
# runtime check pass.
_FUNC_COMPLETENESS_RE = re.compile(
    r"(?i)\b(?:vollst[äa]ndig|komplett|l[üu]ckenlos|alles|alle\s+\w+|"
    r"fully|completely|thoroughly|end[- ]?to[- ]?end)\b[^.!?\n]{0,40}?"
    r"\b(?:getestet|gepr[üu]ft|verifiziert|abgedeckt|tested|verified|"
    r"validated|covered)\b"
)


_FUNC_PREDICATE_RE = re.compile(
    "(?:" + _FUNC_PLAY_PAT + "|" + _FUNC_WORK_PAT + ")", re.IGNORECASE)

# User-interface context markers. Combined with any functional predicate
# they make the claim an interactive one.
_FUNC_UI_MARKER_RE = re.compile(
    r"(?i)(?:"
    r"\bbrowser\b|\bwebseite\b|\bweb\s?page\b|\bweb-?app\b|\bwebapp\b|"
    r"\btastatur\b|\bpfeiltasten\b|\bleertaste\b|\btasten\b|\bkeyboard\b|"
    r"\barrow\s+keys\b|\bkey\s?press\b|\bwasd\b|\bmaus\b|\bmouse\b|"
    r"\bklick\w*\b|\bclick\w*\b|\bdrag\b|\bbutton\b|\bschaltfl(?:ä|ae)che\b|"
    r"\bgui\b|\bui\b|\boberfl(?:ä|ae)che\b|\bfrontend\b|\bfront-end\b|"
    r"\bwidget\w*\b|\bcanvas\b|\bipyevents\b|\bipywidgets\b|"
    r"\binteraktiv\w*\b|\binteractive\b|\banimation\b|\bon\s+screen\b|"
    r"\bim\s+browser\b|\bin\s+the\s+browser\b"
    r")")

# Explicit disclosure that the thing was NOT verified. Distinct from the
# generic hedge markers: an answer that says "I could not verify it in a
# browser" has already told the truth this guard exists to enforce.
_FUNC_DISCLOSURE_RE = re.compile(
    r"(?i)(?:"
    # German
    r"nicht\s+(?:verifiziert|verifizierbar|getestet|(?:ü|ue)berpr(?:ü|ue)ft|"
    r"gepr(?:ü|ue)ft|best(?:ä|ae)tigt|nachgewiesen|belegt)|"
    r"(?:konnte|kann|konnten|k(?:ö|oe)nnen)\s+(?:ich|wir)?\s*nicht\s+"
    r"[^\n]{0,40}?(?:verifizieren|(?:ü|ue)berpr(?:ü|ue)fen|pr(?:ü|ue)fen|"
    r"testen|best(?:ä|ae)tigen|nachweisen)|"
    r"ungetestet|ungepr(?:ü|ue)ft|unverifiziert|ohne\s+Gew(?:ä|ae)hr|"
    r"nicht\s+im\s+Browser\s+getestet|headless\s+nicht|"
    # English
    r"(?:could|can|was|were|am|is|are)\s*n(?:o|')?t\s+"
    r"[^\n]{0,40}?(?:verif(?:y|ied)|check(?:ed)?|confirm(?:ed)?|test(?:ed)?|"
    r"validat(?:e|ed))|"
    r"unable\s+to\s+(?:verify|check|confirm|test)|"
    r"\bunverified\b|\buntested\b|\bnot\s+verified\b|\bnot\s+tested\b|"
    r"no\s+way\s+to\s+(?:verify|check|confirm)|"
    r"have\s*n(?:o|')?t\s+(?:verified|tested|checked|confirmed)"
    r")")

# Negation near the predicate: "funktioniert nicht", "does not work".
_FUNC_NEGATION_RE = re.compile(
    r"(?i)(?:\bnicht\b|\bkein(?:e|en|er|em)?\b|\bnie\b|\bnot\b|n't\b|"
    r"\bcannot\b|\bfails?\b|\bfailing\b|\bbroken\b|\bkaputt\b|"
    r"\bfehlerhaft\b|\bnoch\s+nicht\b)")

# Conditional / subordinate / instructional framing: a requirement, an open
# question or a recipe — not an assertion about present state.
_FUNC_CONDITIONAL_RE = re.compile(
    r"(?i)(?:\bdamit\b|\bsoll\w*\b|\bfalls\b|\bwenn\b|\bsobald\b|\bob\b|"
    r"\bbevor\b|\bum\s+[^\n]{0,30}?\s+zu\s+\w+en\b|"
    r"\bshould\b|\bwould\b|\bonce\b|\bif\b|\bwhether\b|\bin\s+order\s+to\b|"
    r"\bto\s+make\s+it\b|\bmust\b|\bneeds?\s+to\b)")

# Explanatory usage: describing HOW something works is not a claim THAT it
# works ("so funktioniert der Parser", "here is how it works").
_FUNC_EXPLANATORY_RE = re.compile(
    r"(?i)(?:\bwie\s+funktionier(?:t|en)\b|\bso\s+funktionier(?:t|en)\b|"
    r"funktionier(?:t|en)\s+(?:wie\s+folgt|folgenderma(?:ß|ss)en|so\b)|"
    r"\bhow\s+(?:it|this|that|they|the\s+\w+)\s+works?\b|"
    r"\bhere\s+is\s+how\b|\bthe\s+way\s+it\s+works\b)")

# Software artifacts a functional claim can be about. Required for the
# turn-level kind (no artifact file named): it keeps the guard on claims
# about produced software and off general prose ("that works for me").
_ARTIFACT_NOUN_RE = re.compile(
    r"(?i)\b(?:skript|script|app|anwendung|programm|program|code|tool|"
    r"server|notebook|modul|module|cli|pipeline|setup|build|installation|"
    r"paket|package|datei|file|spiel\w*|game|dashboard|widget|"
    r"web-?app|befehl|command|kommando|demo|prototyp\w*|prototype)\b")

# The claim is reported back from the user, not asserted by the agent.
_FUNC_USER_SOURCE_RE = re.compile(
    r"(?i)(?:wie\s+du\s+(?:best(?:ä|ae)tigt|geschrieben|gesagt|berichtet)|"
    r"laut\s+(?:deiner|deinem|ihrer)|du\s+hast\s+(?:best(?:ä|ae)tigt|"
    r"berichtet|geschrieben)|dein(?:er|em)?\s+R(?:ü|ue)ckmeldung|"
    r"you\s+(?:confirmed|reported|said|told\s+me)|"
    r"as\s+you\s+(?:confirmed|reported|said)|per\s+your\s+report)")

# Artifact tokens a functional claim can be ABOUT (runnable files).
_ARTIFACT_EXTS = (
    "py", "sh", "bash", "ipynb", "js", "mjs", "cjs", "ts", "tsx", "jsx",
    "html", "htm", "jl", "rb", "go", "rs", "pl", "exe", "bat", "ps1",
)
_ARTIFACT_RE = re.compile(
    r"\b((?:[\w\-.]+/)*[\w\-]+\.(?:" + "|".join(_ARTIFACT_EXTS) + r"))\b",
    re.IGNORECASE)

# Sentence boundary: a terminator followed by whitespace/end, or a newline.
# The whitespace lookahead keeps "app.py" from splitting a sentence.
_SENT_BOUND_RE = re.compile(r"[.!?;:](?=\s|$)|\n")

# Work bound: how many predicate matches are examined at most, so a very
# long answer cannot turn the per-match sentence lookup into a hot loop.
_FUNC_MAX_MATCHES = 200


@dataclass(frozen=True)
class FunctionalClaimFlag:
    """One claim that something works, without the session exercising it."""

    claim: str          # normalized claim sentence, truncated
    subject: str        # artifact the claim is about ("" when unnamed)
    kind: str           # "interactive" | "unexercised" | "no_execution"

    def message(self) -> str:
        if self.kind == "interactive":
            return (f"⚠️ Verify: '{self.claim}' asserts interactive or browser "
                    f"behavior that nothing in this session exercised — it "
                    f"cannot be checked headlessly, so state it as unverified.")
        if self.kind == "unexercised":
            return (f"⚠️ Verify: '{self.claim}' — '{self.subject}' was never "
                    f"executed in this session (starting a server does not "
                    f"count) — run it or state the claim as unverified.")
        return (f"⚠️ Verify: '{self.claim}' claims working software, but this "
                f"session executed nothing — run it or state the claim as "
                f"unverified.")


def extract_exec_command(tool_name, tool_input) -> str:
    """Return a normalized command string when ``tool_name`` is an execution
    tool, else ``""``.

    The tool name is kept as the first token so command classification (see
    ``_START_ONLY_CMD_RE``) can tell a foreground run from a backgrounded
    start. Accepts the raw JSON string backends pass around, a decoded dict,
    or a plain command string. Never raises.
    """
    try:
        name = str(tool_name or "").strip().lower().rsplit("__", 1)[-1]
        if not name or name in _EXEC_TOOL_DENY:
            return ""
        if name not in _EXEC_TOOLS and not _EXEC_TOOL_SHAPE.search(name):
            return ""
        payload = tool_input
        if isinstance(payload, str) and payload.strip().startswith("{"):
            try:
                payload = json.loads(payload)
            except (json.JSONDecodeError, ValueError):
                pass
        parts: list[str] = [name]
        if isinstance(payload, dict):
            for key in _EXEC_INPUT_KEYS:
                val = payload.get(key)
                if isinstance(val, str) and val.strip():
                    parts.append(val.strip())
                elif isinstance(val, (list, tuple)):
                    parts.extend(str(v) for v in val
                                 if isinstance(v, (str, int, float)))
        elif isinstance(payload, str) and payload.strip():
            parts.append(payload.strip())
        return " ".join(" ".join(parts).split())[:400]
    except Exception:
        return ""


def _run_commands(exec_commands) -> list[str]:
    """Commands that actually EXERCISED something, i.e. everything except
    start-only launches (servers, backgrounded starts)."""
    out: list[str] = []
    for c in exec_commands or ():
        s = str(c)
        if _START_ONLY_CMD_RE.search(s):
            continue
        out.append(s)
    return out


def _artifact_exercised(subject: str, runs: list[str]) -> bool:
    """True when a foreground command named this artifact — or when a test
    runner ran and the artifact is a test file (runners discover those from
    a directory argument)."""
    base = subject.replace("\\", "/").rsplit("/", 1)[-1].lower()
    if not base:
        return False
    for c in runs:
        low = c.lower()
        if base in low:
            return True
        if _TEST_FILE_RE.match(base) and _TEST_CMD_RE.search(low):
            return True
    return False


def _ui_was_exercised(tools_used, exec_commands) -> bool:
    """True when the session drove a real user interface (browser-automation
    tool or command) — then interactive claims have evidence."""
    for t in tools_used or ():
        if _UI_EXERCISE_TOOL_SHAPE.search(str(t)):
            return True
    for c in exec_commands or ():
        if _UI_EXERCISE_CMD_RE.search(str(c)):
            return True
    return False


def _sentence_span(text: str, start: int, end: int) -> tuple[int, int]:
    """Span of the sentence containing ``text[start:end]``."""
    lo = 0
    for m in _SENT_BOUND_RE.finditer(text, 0, start):
        lo = m.end()
    hi = len(text)
    m = _SENT_BOUND_RE.search(text, end)
    if m is not None:
        hi = m.start()
    return lo, hi


def scan_for_unexercised_functional_claims(
    text: str,
    *,
    exec_commands: Optional[list[str]] = None,
    exec_ledger_available: bool = True,
    tools_used: Optional[frozenset[str] | set[str]] = None,
    max_flags: int = 4,
) -> list[FunctionalClaimFlag]:
    """Scan ``text`` for claims that something the session produced now
    works, and flag those the session never exercised.

    See the section comment for the claim kinds and the exemption rules.
    ``exec_ledger_available=False`` (the caller keeps no record of what ran)
    silences the runtime kind — absence of bookkeeping is not evidence of
    fabrication. Interactive claims are judged independently of that ledger:
    no amount of shell evidence makes browser behavior observable here.

    Deterministic, order-stable (text position), de-duplicated, capped at
    ``max_flags``. Never raises.
    """
    flags: list[FunctionalClaimFlag] = []
    try:
        if not text or not text.strip() or max_flags <= 0:
            return []
        scrubbed = _scrub_keep_offsets(text)
        cmds = [str(c) for c in (exec_commands or ())]
        runs = _run_commands(cmds)
        ui_ok = _ui_was_exercised(tools_used, cmds)
        seen_spans: set[tuple[int, int]] = set()
        seen_keys: set[str] = set()
        examined = 0
        # Completeness claims stand on their own: "alles getestet" needs no
        # functional predicate beside it, and no executed command can
        # establish the absence of untested parts.
        for m in _FUNC_COMPLETENESS_RE.finditer(scrubbed):
            if len(flags) >= max_flags:
                break
            span = _sentence_span(scrubbed, m.start(), m.end())
            sentence = scrubbed[span[0]:span[1]]
            if not sentence.strip():
                continue
            if _is_hedged(scrubbed, m.start(), m.end()):
                continue
            if (_FUNC_DISCLOSURE_RE.search(sentence)
                    or _FUNC_NEGATION_RE.search(sentence)
                    or _FUNC_CONDITIONAL_RE.search(sentence)
                    or _FUNC_EXPLANATORY_RE.search(sentence)
                    or _FUNC_USER_SOURCE_RE.search(sentence)):
                continue
            if scrubbed[span[1]:span[1] + 1] == "?":
                continue
            claim = " ".join(sentence.split())[:100]
            key = f"completeness::{claim.lower()}"
            if key in seen_keys:
                continue
            seen_keys.add(key)
            seen_spans.add(span)
            flags.append(FunctionalClaimFlag(
                claim=claim, subject="", kind="completeness"))

        for m in _FUNC_PREDICATE_RE.finditer(scrubbed):
            if len(flags) >= max_flags or examined >= _FUNC_MAX_MATCHES:
                break
            examined += 1
            span = _sentence_span(scrubbed, m.start(), m.end())
            if span in seen_spans:
                continue
            seen_spans.add(span)
            sentence = scrubbed[span[0]:span[1]]
            if not sentence.strip():
                continue
            # Not an assertion about present state, or already disclosed as
            # unverified — never an enforcement target.
            if _is_hedged(scrubbed, m.start(), m.end()):
                continue
            if (_FUNC_DISCLOSURE_RE.search(sentence)
                    or _FUNC_NEGATION_RE.search(sentence)
                    or _FUNC_CONDITIONAL_RE.search(sentence)
                    or _FUNC_EXPLANATORY_RE.search(sentence)
                    or _FUNC_USER_SOURCE_RE.search(sentence)):
                continue
            tail = scrubbed[span[1]:span[1] + 1]
            if tail == "?":
                continue                      # a question, not a claim
            claim = " ".join(sentence.split())[:100]
            interactive = bool(_FUNC_PLAY_RE.search(sentence)
                               or _FUNC_UI_MARKER_RE.search(sentence))
            if interactive:
                if ui_ok:
                    continue
                kind, subject = "interactive", ""
            else:
                if not exec_ledger_available:
                    continue
                subjects = _ARTIFACT_RE.findall(sentence)
                if subjects:
                    unrun = [s for s in subjects
                             if not _artifact_exercised(s, runs)]
                    if not unrun:
                        continue
                    kind, subject = "unexercised", unrun[0]
                else:
                    if runs:
                        continue      # turn-level rule: something did run
                    if not _ARTIFACT_NOUN_RE.search(sentence):
                        continue      # not a claim about produced software
                    kind, subject = "no_execution", ""
            key = f"{kind}:{subject.lower()}:{claim.lower()}"
            if key in seen_keys:
                continue
            seen_keys.add(key)
            flags.append(FunctionalClaimFlag(
                claim=claim, subject=subject, kind=kind))
    except Exception:
        return flags
    return flags


def functional_claim_caveat(flags: list[FunctionalClaimFlag]) -> str:
    """Visible caveat naming exactly what was claimed to work but never
    exercised.

    This class gets a caveat rather than a forced correction turn: the model
    cannot verify browser or interactive behavior headlessly either, so a
    retry would only invite a second confident assertion. Naming the gap is
    the honest outcome.
    """
    if not flags:
        return ""
    items: list[str] = []
    for f in flags:
        if f.kind == "completeness":
            items.append(
                f"'{f.claim}' — a completeness claim cannot be established "
                "by a test run: it says what was covered, never what was "
                "left out. Name the parts you did NOT exercise")
        elif f.kind == "interactive":
            items.append(f"'{f.claim}' — interactive/browser behavior was "
                         f"never exercised in this session")
        elif f.kind == "unexercised":
            items.append(f"'{f.claim}' — '{f.subject}' was never executed in "
                         f"this session (starting a server is not evidence "
                         f"that it works)")
        else:
            items.append(f"'{f.claim}' — this session executed nothing")
    note = ""
    if any(f.kind == "interactive" for f in flags):
        note = (" Interactive and browser behavior cannot be checked "
                "headlessly here.")
    return (
        "\n\n[verify] Caveat: the following was NOT verified in this session: "
        + "; ".join(items) + "." + note + " Treat it as unconfirmed."
    )
