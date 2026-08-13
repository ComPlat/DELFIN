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

import bisect
import contextvars as _contextvars
import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

from . import german as _german


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
            # A cited DIRECTORY is not a fabrication. Only files were
            # checked here, so naming a folder that exists -- "the forms are
            # in office_analysis/Reisekostenantraege/" -- was classified the
            # same as inventing a path, which is the hardest flag and forces
            # a correction turn. Observed in the field on a folder that held
            # the 31 files the sentence was about, and the answer that came
            # back from the forced retry was worse than the one it replaced.
            exists = False
            is_dir = False
            try:
                p = Path(path)
                if p.is_absolute():
                    target = p
                elif root is not None:
                    target = root / path
                else:
                    target = Path(path)
                exists = target.is_file()
                if not exists:
                    is_dir = target.is_dir()
            except OSError:
                exists = False
            if is_dir:
                # It is there. Whether the agent looked inside is the
                # "unread" question, never the "you made this up" one.
                flags.append(CodeClaimFlag(path=path, line=line,
                                           kind="unread"))
            elif not exists:
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


# ---------------------------------------------------------------------------
# What the tools actually returned this turn
# ---------------------------------------------------------------------------
#
# The turn-level gate below used to read: any lookup tool ran, or any
# calculation-output-like file was observed, therefore nothing is flagged.
# Measured on the shipped scanner:
#
#     3 flags | Die Messung ergab 2.31 eV für S1, 3.90 eV für S2 …
#     0 flags | ... with evidence_tools_used={"read_document"}
#
# One read of one unrelated spreadsheet exempted every energy in the
# answer. That is the same defect the office ledger was built to close on
# the other side of the house, and it has the same repair: a claim is
# grounded when the FIGURE appears in what a tool returned, not when some
# tool ran at all.
#
# The office ledger takes its figures from result dicts on purpose --
# never from rendered prose, so two readings of the same facts cannot
# disagree. Here there is no dict: the evidence for "2.31 eV" is the text
# of an ORCA output the agent read, and that text IS the primary source
# rather than a rendering of one.
#
# Bounded twice over. Only the head of a tool result is scanned, and only
# so many numbers are kept, because an .out file is megabytes and this
# module has already had one incident where a single long line made the
# grounding scan take minutes.
_TOOL_OUTPUT_SCAN_CHARS = 200_000
MAX_OBSERVED_NUMBERS = 5000

_OBSERVED_NUMBER_RE = re.compile(r"(?<![\w.])[-+−]?\d+(?:\.\d+)?(?![\w.])")

_observed_numbers: "_contextvars.ContextVar[Optional[list[float]]]" = (
    _contextvars.ContextVar("delfin_observed_numbers", default=None))
# Whether this pool can answer "the tools did NOT return that number".
# A tool result reaches the recorder as a HEAD slice, so a truncated one
# leaves the pool with a hole in it of unknown size -- and absence from a
# pool with a hole is not evidence of anything. This is the same rule the
# rest of the framework already runs on: unobserved is not zero.
_observed_complete: "_contextvars.ContextVar[bool]" = (
    _contextvars.ContextVar("delfin_observed_complete", default=True))


def reset_observed_numbers() -> None:
    """Forget what the previous turn's tools returned.

    Per TURN, like the other evidence ledgers here: the question is
    whether THIS answer states a figure THIS turn produced. Carried over,
    a stale energy would ground a later, unrelated one -- which is the
    failure this exists to catch, one turn late.

    Cleared to None and NOT to an empty list. The two say different
    things and the scanner turns on the difference: None is "no result
    reached this pool, so there is nothing to judge against", an empty
    list is "results came back and carried no number at all". Setting []
    here made every turn look like the second, and a turn whose evidence
    this pool cannot see -- one where the results never arrive through
    the recording path -- would have had its every quantity flagged.
    """
    _observed_numbers.set(None)
    _observed_complete.set(True)


def record_tool_numbers(output: Any, *, truncated: bool = False) -> int:
    """Take the numbers out of one tool result. Never raises.

    Returns how many were added, so a caller can tell "nothing numeric
    came back" from "this was never called" -- the distinction the gate
    below turns on.

    ``truncated`` says the result reached here as a HEAD slice, so the
    numbers past the cut were never seen. One such result puts a hole of
    unknown size in the pool, and a pool with a hole cannot answer "the
    tools did not return that number" -- so from then on the turn reports
    itself as unable to say. Getting this wrong would flag a chemist's
    correctly-quoted energy from deep inside an .out file.
    """
    try:
        if truncated:
            _observed_complete.set(False)
        pool = _observed_numbers.get()
        if pool is None:
            pool = []
            _observed_numbers.set(pool)
        if len(pool) >= MAX_OBSERVED_NUMBERS:
            return 0
        body = str(output or "")[:_TOOL_OUTPUT_SCAN_CHARS]
        added = 0
        for match in _OBSERVED_NUMBER_RE.finditer(body):
            try:
                pool.append(float(match.group(0).replace("−", "-")))
            except ValueError:
                continue
            added += 1
            if len(pool) >= MAX_OBSERVED_NUMBERS:
                break
        return added
    except Exception:
        return 0


def observations_are_complete() -> bool:
    """False once any result arrived truncated. Said out loud rather than
    folded into the pool, so a caller can tell the two apart."""
    return bool(_observed_complete.get())


def observed_numbers() -> Optional[list[float]]:
    """Every number this turn's tools returned, or None if it cannot say.

    None and [] are different answers and the scanner treats them so:
    None means there is nothing to judge against -- no result reached
    this pool, or one arrived cut short and the pool has a hole in it.
    An empty list means results came back and carried no number at all,
    which is a fact about the turn and gets checked like any other.
    """
    if not _observed_complete.get():
        return None
    return _observed_numbers.get()


def _grounded_in_observations(value: float, pool: list[float]) -> bool:
    """Is this value one the tools returned, or a difference of two?

    The difference is in because an energy gap is the most ordinary
    derived quantity there is: an answer that reads 2.31 and 3.90 out of
    an output and reports a 1.59 eV gap has done arithmetic, not
    invention, and flagging it would teach the model to stop reporting
    gaps. Bounded to the same base size the office ledger uses, because
    the derived set grows with its square and a large enough base makes
    every number derivable.
    """
    tolerance = max(abs(value) * 1e-4, 5e-3)
    for known in pool:
        if abs(value - known) <= tolerance:
            return True
    base = pool[:MAX_DERIVATION_BASE]
    for i, a in enumerate(base):
        for b in base[i + 1:]:
            if abs(value - abs(a - b)) <= tolerance:
                return True
    return False


# How many observed numbers take part in deriving a difference.
MAX_DERIVATION_BASE = 24

_CLAIM_VALUE_RE = re.compile(r"[-+−]?\d+(?:\.\d+)?")


def _claim_is_observed(quantity: str, pool: list[float]) -> bool:
    """Does the number inside a matched claim come from the tools?"""
    match = _CLAIM_VALUE_RE.search(str(quantity or ""))
    if not match:
        return False
    try:
        value = float(match.group(0).replace("−", "-"))
    except ValueError:
        return False
    return _grounded_in_observations(value, pool)


def scan_for_unsourced_quantities(
    text: str,
    *,
    observed_files: Optional[frozenset[str] | set[str]] = None,
    evidence_tools_used: Optional[frozenset[str] | set[str]] = None,
    numbers: Optional[list[float]] = None,
    max_flags: int = 6,
) -> list[QuantityClaimFlag]:
    """Scan ``text`` for physical-quantity claims made in a turn with no
    evidence act.

    A claim is a number immediately followed by a physical unit (eV,
    kcal/mol, kJ/mol, Hartree/Eh, nm, cm-1, Debye, Angstrom/Å, ps, fs,
    K, kcal, GHz). Backticked spans, fenced code blocks, blockquoted
    lines, percentages and version strings never count.

    Evidence is PER CLAIM. ``numbers`` is what this turn's tools actually
    returned (``observed_numbers()``); a quantity is grounded when its
    value is in there, or is the difference of two that are. What is not
    in there is flagged even in a turn full of tool calls -- which is the
    whole point, because the old rule exempted every energy in an answer
    on the strength of one read of one unrelated spreadsheet.

    ``numbers=None`` means no tool returned anything to judge against. The
    turn-level shape survives only for that case: a calculation output
    observed or a lookup tool used still exempts the turn, because
    without a pool there is nothing to check a value against and flagging
    an answer whose evidence this function simply cannot see would be a
    guard punishing work it is blind to.

    Deterministic, order-stable (text position), de-duplicated, capped at
    ``max_flags``. Never raises.
    """
    flags: list[QuantityClaimFlag] = []
    try:
        if not text or not text.strip() or max_flags <= 0:
            return []
        pool = numbers
        if pool is None:
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
            if pool is not None and _claim_is_observed(qty, pool):
                continue
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
# Claim scoping — a disarm excuses the claim it stands next to, nothing else
# ---------------------------------------------------------------------------
#
# Every exemption in this module (a hedge, a question, a conditional, a
# negation) has to be LOCAL to the claim it excuses. Scoped to the line, a
# markdown paragraph is one line, so one adverb three sentences later
# silenced a confident number at the top of it. Measured on the shipped
# code: "The S1 energy is 2.31 eV. The geometry is around the minimum."
# was not flagged, and nothing about the energy had changed.
#
# The boundary is deliberately CONSERVATIVE — it splits only where the next
# sentence unmistakably begins, so an abbreviation ("ca. 2.31 eV"), a file
# name ("app.py läuft"), a decimal and a version string never separate a
# disarm from the claim it belongs to. Erring here costs a missed flag, and
# the opposite error costs a correct answer being caveated, which is the
# one this project has been bitten by twice.
#
# There is deliberately NO "terminator at end of text" alternative: the
# backward scan runs over a bounded region, and inside one a `$` anchors to
# the region's end, so every abbreviation immediately left of the claim
# read as a sentence end. "Der Wert liegt bei ca. 2.31 eV" then split right
# between the hedge and the number it hedges. A missing boundary at the end
# of the text costs nothing — the window already stops there.
_CLAIM_SENT_BOUND_RE = re.compile(
    r"[.!?]+(?=\s+[\"'“(\[]?[A-ZÄÖÜ])"   # ". Next" / ". „Next"
    r"|\n"                                 # a line break always separates
)

# Clause boundary INSIDE a sentence. Used where the order matters: a
# subordinate clause that follows a claim does not retract it.
_CLAUSE_BOUND_RE = re.compile(r"[,;:]|\s[-–—]\s")


@lru_cache(maxsize=8)
def _sentence_bounds(text: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """(starts, ends) of every sentence boundary in ``text``.

    Computed once per answer and reused by every claim in it. The scanners
    ask for a sentence per match, and a long answer has thousands of
    matches — rescanning the text around each of them cost 15 s on 100 kB.
    Python caches a string's hash after the first use, so the lookup here
    is a dict probe plus an identity check.
    """
    found = [(m.start(), m.end()) for m in _CLAIM_SENT_BOUND_RE.finditer(text)]
    return tuple(s for s, _ in found), tuple(e for _, e in found)


def _claim_sentence_span(text: str, start: int, end: int) -> tuple[int, int]:
    """Span of the sentence containing ``text[start:end]``."""
    starts, ends = _sentence_bounds(text)
    i = bisect.bisect_right(ends, start) - 1
    lo = ends[i] if i >= 0 else 0
    j = bisect.bisect_left(starts, end)
    return lo, (starts[j] if j < len(starts) else len(text))


def _claim_sentence(text: str, start: int, end: int) -> str:
    lo, hi = _claim_sentence_span(text, start, end)
    return text[lo:hi]


def _claim_scope_with_lead_in(text: str, start: int, end: int) -> str:
    """The claim's sentence plus the one BEFORE it.

    Some disclosures are written as a lead-in — "Die Spalte ist mehrdeutig
    — welche Lesart gilt? Summe wäre 25.136 EUR." — and that order is the
    honest one. A disclosure that follows the figure is not: it is the
    note that gets left behind when somebody copies the number out, which
    is the failure the guard exists for. So the window reaches backwards
    only.
    """
    lo, hi = _claim_sentence_span(text, start, end)
    if lo > 0:
        pos = max(lo - 1, 0)
        lo = _claim_sentence_span(text, pos, pos)[0]
    return text[lo:hi]


def _clause_index(sentence: str, pos: int) -> int:
    """How many clause boundaries precede ``pos`` inside ``sentence``."""
    return sum(1 for _ in _CLAUSE_BOUND_RE.finditer(sentence, 0, max(pos, 0)))


def _governs(sentence: str, rx: "re.Pattern[str]", pos: int) -> bool:
    """True when ``rx`` matches in the clause holding ``pos`` or in one that
    PRECEDES it.

    Order carries meaning: "Wenn Python 3.11 installiert ist, läuft das
    Skript" states a condition, while "Das Skript läuft jetzt fehlerfrei,
    wenn Python 3.11 installiert ist" asserts the present state and then
    names a precondition. Both contain the same word; only the first is a
    non-assertion. A whole-sentence search could not tell them apart.
    """
    here = _clause_index(sentence, pos)
    for m in rx.finditer(sentence):
        if _clause_index(sentence, m.start()) <= here:
            return True
    return False


# ---------------------------------------------------------------------------
# Hedge detection — explicitly uncertain claims are never enforcement targets
# ---------------------------------------------------------------------------
#
# Every claim scanner below (and the quantity scanner above) skips claims
# whose containing SENTENCE carries an explicit uncertainty marker:
# enforcement exists to stop CONFIDENT fabrication, and an answer that says
# it is guessing has already disclosed its grounding status.

# THE SAME list the office figure guard uses — see delfin/agent/german.py.
# They used to be two, and the German half here was adverbs only while the
# English half had a full first-person vocabulary ("i think", "from
# memory", "iirc"). Measured on the shipped list, none of these was a
# hedge:
#
#     Die Summe liegt bei rund 45.000 EUR.
#     Ich schätze die Summe auf 45.000 EUR.
#     Ich glaube, es sind 45.000 EUR.
#     Soweit ich weiß sind es 45.000 EUR.
#
# and the first of them WAS a hedge to the other guard, in the same
# answer.
_HEDGE_MARKERS = _german.HEDGE_RE


@lru_cache(maxsize=8)
def _hedge_positions(text: str) -> tuple[int, ...]:
    """Where every uncertainty marker in ``text`` starts.

    Collected once per answer, like the sentence bounds: an answer with
    thousands of numbers in it asks this question thousands of times, and
    each answer used to cost a fresh scan of the surrounding text."""
    return tuple(m.start() for m in _HEDGE_MARKERS.finditer(text))


def _is_hedged(text: str, start: int, end: int) -> bool:
    """True when the SENTENCE containing ``text[start:end]`` carries an
    explicit uncertainty marker.

    Sentence-scoped, not line-scoped: a markdown paragraph is a single
    line, so "around" in the third sentence used to clear a bare number in
    the first. The hedge has to stand next to the claim it excuses.
    """
    lo, hi = _claim_sentence_span(text, start, end)
    marks = _hedge_positions(text)
    i = bisect.bisect_left(marks, lo)
    return i < len(marks) and marks[i] < hi


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
#
# The component repetitions are LENGTH-BOUNDED, and that is what keeps this
# cheap. Written unbounded, the final component is greedy and then has to
# find the extension dot, so a long run of word characters containing no
# dot makes the engine give back one character at a time from every
# starting position -- quadratic. Measured on one unbroken 20 kB line: 2.2s
# per pattern, six patterns, on every answer and every delegate report.
# Ordinary multi-line text of the same size costs 16ms, so it only ever bit
# on a pasted blob or a minified dump. That is precisely the input nobody
# anticipates, and the symptom is a turn that stops responding.
#
# 200 characters per component is far past any real path and bounds the
# work per position at a constant.
_LOC_COMPONENT = r"[\w\-][\w.\-]{0,200}"
_LOC_PATH = (
    r"(?:" + _LOC_COMPONENT + r"/){0,20}" + _LOC_COMPONENT +
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


def verification_marker(new_files) -> str:
    """The line printed when a forced correction turn actually verified.

    Derived from the RE-SCAN of the correction, never from a flag: the
    flag that used to drive the dashboard's "✓ Self-verification …" was
    assigned BEFORE the correction turn was even attempted, so the
    checkmark also appeared after a correction that raised a NEW claim,
    directly above the caveat saying the claims stayed unverified, and
    when the model had done nothing but add hedges.

    Emitted by the engine rather than a UI layer, so the CLI and headless
    paths — where a guard-forced retry was previously invisible — see it
    too.
    """
    names = sorted({str(p).replace("\\", "/").rsplit("/", 1)[-1]
                    for p in (new_files or ()) if str(p).strip()})
    if not names:
        return ""
    shown = ", ".join(names[:3])
    more = f" (+{len(names) - 3})" if len(names) > 3 else ""
    return (
        f"\n\n[verify] Self-check: die Angaben oben wurden gegen "
        f"{shown}{more} geprüft, gelesen im Korrekturzug."
    )


# ---------------------------------------------------------------------------
# The language a caveat is written in
# ---------------------------------------------------------------------------
#
# German, all of them. The line drawn here is not "user-facing text" —
# it is who READS the string and whether anything translates it:
#
#   * A tool result, an error, a note the model gets back — ENGLISH. The
#     model reads it and answers the user in the user's language, so it
#     is translated on the way out. Every message in office.py is that
#     kind and stays English.
#   * A caveat — GERMAN. It is appended to the FINISHED answer, after the
#     model's last turn; nothing comes after it to translate it, and it
#     asks the reader to do something ("bitte nachrechnen, bevor die Zahl
#     weitergegeben wird"). An instruction the reader skips is not a
#     warning.
#
# The ``[verify] Caveat:`` / ``[verify] Self-check:`` tags stay as they
# are: they are machine markers the engine and the tests key on, not
# prose, and translating a marker breaks what reads it.
#
# It used to be decided per function rather than by a rule, so one German
# answer about one spreadsheet came back with three caveats — German,
# English, German — two of them from this file, and the English one
# quoted the German noun back at the reader ("This answer states 31
# Belege"). The users of this framework write German; a warning they
# stop reading is worse than no warning.
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
        "\n\n[verify] Caveat: die folgenden Angaben sind unbelegt — kein "
        "Datei-Zugriff und keine Recherche in dieser Sitzung deckt sie ab: "
        + ", ".join(items) + ". Bitte als unbestätigt behandeln."
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
    # …and the same disclosure with German verb-second word order:
    # "getestet habe ich es allerdings nicht". The negation sits at the end
    # of the clause, which is where German puts it, so the "nicht <verb>"
    # form above cannot see it.
    r"(?:verifiziert|getestet|gepr(?:ü|ue)ft|(?:ü|ue)berpr(?:ü|ue)ft|"
    r"best(?:ä|ae)tigt|nachgewiesen|belegt)\s+"
    r"(?:habe|hab|haben|hat|hatte|wurde|wurden|ist|sind|war|waren)\b"
    r"[^\n]{0,40}?\bnicht\b|"
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


# Exemptions that hold wherever they stand in the sentence, because they
# are the honest disclosure this guard exists to elicit ("I could not
# verify it") or an attribution to the user's own report. Retracting those
# for standing after the predicate would punish exactly the phrasing the
# framework asks for.
_FUNC_SENTENCE_WIDE = (_FUNC_DISCLOSURE_RE, _FUNC_USER_SOURCE_RE)

# Exemptions that only hold when they GOVERN the predicate — same clause,
# or a clause before it. A trailing "…, wenn Python 3.11 installiert ist"
# names a precondition; it does not withdraw the assertion in front of it.
_FUNC_CLAUSE_LOCAL = (_FUNC_NEGATION_RE, _FUNC_CONDITIONAL_RE,
                      _FUNC_EXPLANATORY_RE)


def _func_disarmed(sentence: str, offset: int) -> bool:
    """True when ``sentence`` retracts the functional claim at ``offset``
    (an offset relative to the start of the sentence)."""
    for rx in _FUNC_SENTENCE_WIDE:
        if rx.search(sentence):
            return True
    for rx in _FUNC_CLAUSE_LOCAL:
        if _governs(sentence, rx, offset):
            return True
    return False


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
            if _func_disarmed(sentence, m.start() - span[0]):
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
            if _func_disarmed(sentence, m.start() - span[0]):
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
                f"'{f.claim}' — eine Vollständigkeitsaussage lässt sich "
                "durch einen Testlauf nicht belegen: er sagt, was geprüft "
                "wurde, nie was ausgelassen wurde. Bitte benennen, was NICHT "
                "ausgeführt wurde")
        elif f.kind == "interactive":
            items.append(f"'{f.claim}' — interaktives Verhalten bzw. das "
                         f"Verhalten im Browser wurde in dieser Sitzung nie "
                         f"ausgeführt")
        elif f.kind == "unexercised":
            items.append(f"'{f.claim}' — '{f.subject}' wurde in dieser "
                         f"Sitzung nie ausgeführt (einen Server zu starten "
                         f"belegt nicht, dass er funktioniert)")
        else:
            items.append(f"'{f.claim}' — in dieser Sitzung wurde nichts "
                         f"ausgeführt")
    note = ""
    if any(f.kind == "interactive" for f in flags):
        note = (" Interaktives Verhalten und Browser-Verhalten lassen sich "
                "hier ohne Anzeige nicht prüfen.")
    return (
        "\n\n[verify] Caveat: das Folgende wurde in dieser Sitzung NICHT "
        "geprüft: " + "; ".join(items) + "." + note
        + " Bitte als unbestätigt behandeln."
    )


# ---------------------------------------------------------------------------
# Totals over a column the reader could not read
# ---------------------------------------------------------------------------
#
# A spreadsheet column of values like "8.986" reads as 8986 or as 8.986
# depending on the convention, and when nothing in the column settles it the
# reader says so in its notes. Observed on the benchmark: the model was told
# exactly that, chose a reading anyway, and answered with a single total —
# with its own note admitting the assumption. In a report handed on to
# somebody else that figure is indistinguishable from a measured one.
#
# The prompt already forbids it and the tool result already says it. Neither
# binds, which is why this is a scanner.

# The gap between the column name and the phrase carries the example
# value — "values like '8.986'" — so it must allow the dot that made the
# value ambiguous in the first place. A character class excluding '.' put
# the one character the note is about outside the pattern.
_AMBIGUOUS_COLUMN_RE = re.compile(
    r"column '([^']{1,60})':.{0,120}?read as two different numbers")

# A stated total: a total word, then a number within the same clause.
_TOTAL_CLAIM_RE = re.compile(
    r"(?i)(gesamt(?:wert|summe|betrag)?|summe|insgesamt|total|zusammen)"
    r"[^.\n]{0,40}?(\d[\d.,]{2,})")

# Signs the answer is offering the reading rather than asserting it. These
# are the GOOD answer, and the scanner must not fire on them.
#
# TYPED, and read in the total's own sentence. A bare "?" used to be in
# this set and was searched over the WHOLE answer, so
#
#     "Gesamtsumme der Spalte Anschaffungswert: 45.231,50 EUR.
#      Soll ich noch nach Kostenstelle gruppieren?"
#
# passed while the same text without the trailing offer was flagged.
# Nothing about the total had changed — only the closing pleasantry. A
# question mark alone therefore no longer disarms anything; the sentence
# has to name the ambiguity, or state both readings (see
# ``_states_both_readings``).
# Naming the assumption is deliberately NOT in here: a note saying which
# reading was picked travels with the answer, not with the number, and the
# reader who copies the figure out leaves it behind. The answer has to
# present the ambiguity or ask.
_OFFERS_THE_CHOICE_RE = re.compile(
    r"(?i)(mehrdeutig|uneindeutig|nicht eindeutig|"
    r"zwei(?:erlei)? (?:lesart|deutung)|beide lesarten|"
    r"welche[srmn]? (?:lesart|konvention|format|schreibweise)|"
    r"ambiguous|which reading|which convention|both readings)")

# Number tokens inside one sentence, used to spot an answer that puts BOTH
# readings of the same figure side by side ("8.986 = 8986", "25.136 EUR /
# 25,136 EUR"). Two spellings of one digit sequence are the ambiguity made
# visible, which is the answer this guard is asking for.
_NUMBER_TOKEN_RE = re.compile(r"\d[\d.,]*\d|\d")


def _states_both_readings(sentence: str) -> bool:
    """True when the sentence spells one digit sequence two different ways."""
    seen: dict[str, str] = {}
    for m in _NUMBER_TOKEN_RE.finditer(sentence):
        raw = m.group(0).rstrip(".,")
        digits = raw.replace(".", "").replace(",", "")
        if not digits:
            continue
        previous = seen.get(digits)
        if previous is not None and previous != raw:
            return True
        seen.setdefault(digits, raw)
    return False


def extract_ambiguous_columns(tool_output: str) -> list[str]:
    """Column names a document reader reported as undecidable.

    Reads the reader's own note rather than re-deriving the judgement, so
    the two cannot disagree about which columns are in question.
    """
    if not tool_output:
        return []
    return [m.group(1) for m in _AMBIGUOUS_COLUMN_RE.finditer(str(tool_output))]


def scan_for_totals_over_ambiguous_columns(
    text: str, ambiguous_columns: list[str] | tuple[str, ...],
) -> list[str]:
    """Columns whose ambiguity the answer states a total over anyway.

    Returns the column names to warn about, empty when the answer either
    states no total or presents the ambiguity with it. Presenting both
    readings and asking which applies is the correct answer and must pass:
    the failure is a single figure offered as the value.
    """
    if not text or not ambiguous_columns:
        return []
    # A total whose OWN sentence presents the ambiguity is the good answer.
    # One that does not is asserted, whatever the rest of the reply says.
    asserted = False
    for m in _TOTAL_CLAIM_RE.finditer(text):
        scope = _claim_scope_with_lead_in(text, m.start(), m.end())
        if _OFFERS_THE_CHOICE_RE.search(scope):
            continue
        if _states_both_readings(scope):
            continue
        asserted = True
        break
    if not asserted:
        return []
    hit = []
    lowered = text.lower()
    for column in ambiguous_columns:
        name = str(column or "").strip()
        if not name:
            continue
        # Only when the answer is actually about that column: a total in a
        # reply that never mentions it is some other figure.
        if name.lower() in lowered and name not in hit:
            hit.append(name)
    return hit


# A count in an answer: "31 PDF-Dateien", "31 files", "29 rows".
#
# A CLOSED noun list was the whole guard, and it decided the outcome on a
# word the user never chose. Over one and the same 29-item listing, "31
# PDF-Dateien" was caught and "31 Rechnungen" was not; "Belege",
# "Formulare", "Datensätze" and "Positionen" walked through as well. So the
# shape is general — a two-plus-digit integer followed by a plural noun —
# and only the things that are MEASURES rather than items are excluded.
_COUNT_CLAIM_RE = re.compile(
    r"(?<![\w.,])(\d{2,6})\s+"
    r"([A-Za-zÄÖÜäöüß][\wÄÖÜäöüß]*(?:-[\wÄÖÜäöüß]+)*)\b"
)

# What a counted thing looks like, without a list of which things.
#
# German nouns are capitalised, and a count of them is plural: "31 Belege",
# never "31 Beleg". English plurals end in "s". Both halves are needed —
# on the German rule alone, "8.986 liest sich als 8986 oder als 8,986"
# produced the count claim "8986 oder als", because "oder" ends in a German
# plural suffix and "als" in an English one. Function words are short and
# lowercase; the two rules together exclude them without naming any.
_PLURAL_SUFFIXES_DE = ("en", "er", "se", "ze", "e", "n")
_ENGLISH_NON_PLURAL_ENDINGS = ("ss", "us", "is")

# Words that pass the shape but count no ITEMS: durations, sizes, money,
# physical units, and the chars/lines vocabulary the truncation markers
# themselves use. Compared against the last hyphen segment, lower-cased.
_COUNT_STOPWORDS = frozenset({
    # time
    "sekunden", "minuten", "stunden", "tage", "wochen", "monate", "jahre",
    "millisekunden", "seconds", "minutes", "hours", "days", "weeks",
    "months", "years", "milliseconds",
    # size / measure
    "zeichen", "chars", "characters", "bytes", "kilobytes", "megabytes",
    "gigabytes", "pixel", "pixels", "meter", "metern", "kilometer",
    "zentimeter", "millimeter", "meters", "metres", "kilometers", "miles",
    "gramm", "kilogramm", "grams", "kilograms", "tonnen", "liter", "litres",
    "liters", "grade", "degrees", "prozente", "percent",
    # money
    "euros", "euro", "cents", "dollars", "franken", "pfund",
    # frequency / vague
    "male", "times", "mal", "prozentpunkte",
})


# A second qualifying word after the noun, so an adjective in front of the
# real noun is reported whole ("31 neue Dateien", not "31 neue").
_COUNT_NEXT_WORD_RE = re.compile(
    r"[ \t]+([A-Za-zÄÖÜäöüß][\wÄÖÜäöüß]*(?:-[\wÄÖÜäöüß]+)*)\b")


def _is_countable_noun(word: str) -> bool:
    """True when ``word`` names things one can count."""
    head = word.rsplit("-", 1)[-1]
    tail = head.lower()
    if len(tail) < 3 or not tail.isalpha() or tail in _COUNT_STOPWORDS:
        return False
    if head[:1].isupper():                      # German noun
        return tail.endswith(_PLURAL_SUFFIXES_DE) or tail.endswith("s")
    return (len(tail) >= 4 and tail.endswith("s")   # English plural
            and not tail.endswith(_ENGLISH_NON_PLURAL_ENDINGS))


def _count_claims(text: str) -> list[tuple[int, str]]:
    """(number, claim text) for every counted-things claim in ``text``.

    Shared by both count guards so the two can never disagree about what
    counts as a count. A modifier is allowed between the number and the
    noun, so "31 neue Dateien" reads back as itself and not as "31 neue".
    """
    out: list[tuple[int, str]] = []
    body = text or ""
    for m in _COUNT_CLAIM_RE.finditer(body):
        first = m.group(2)
        if first.rsplit("-", 1)[-1].lower() in _COUNT_STOPWORDS:
            continue                    # "45 Minuten Videos" counts minutes
        nxt = _COUNT_NEXT_WORD_RE.match(body, m.end())
        if _is_countable_noun(first):
            claim = m.group(0)
            if nxt is not None and _is_countable_noun(nxt.group(1)):
                claim = body[m.start():nxt.end()]
        elif nxt is not None and _is_countable_noun(nxt.group(1)):
            claim = body[m.start():nxt.end()]
        else:
            continue
        try:
            out.append((int(m.group(1)), " ".join(claim.split())))
        except (TypeError, ValueError):
            continue
    return out


def scan_for_counts_over_truncated_output(
    text: str, truncated_tools: list[str] | tuple[str, ...],
) -> list[str]:
    """Counts an answer states after its only source was cut short.

    The failure this exists for: fill_series reported 29 complete of 31,
    the listing that would have settled it was cut at "... 2,000 chars
    total", and the answer asserted "31 PDF-Dateien verifiziert" while
    naming 29. The truncation marker was present and simply not honoured.

    Returns the counts to warn about, empty when nothing was truncated or
    the answer states no count. Deliberately does not try to decide
    whether the number is right -- only that its source could not have
    supported it.
    """
    if not text or not truncated_tools:
        return []
    try:
        return [claim for _n, claim in _count_claims(text)][:4]
    except Exception:
        return []


# A numbered or bulleted list item at the start of a line.
_LIST_ITEM_RE = re.compile(r"(?m)^\s*(?:\d{1,3}[.)]\s+|[-*•]\s+)\S")


def scan_for_count_vs_enumeration(text: str) -> list[tuple[int, int]]:
    """(claimed, listed) pairs where an answer states N and then lists M.

    No ledger and no tool trace: this reads the answer against itself.
    The field case is exact -- "ich habe 31 PDF-Dateien verifiziert",
    then items 1 to 29, then a table restating 31. Three inconsistent
    counts in one message, after both existing guard layers had run,
    because both compare the answer to EVIDENCE and neither compares it
    to itself.

    Only fires when the list is long enough to be the enumeration the
    number refers to (a count of 31 followed by two examples is not a
    contradiction) and when the gap is real rather than an off-by-one
    from a header row.
    """
    if not text:
        return []
    try:
        listed = len(_LIST_ITEM_RE.findall(text))
        if listed < 3:
            return []
        out: list[tuple[int, int]] = []
        for claimed, _claim in _count_claims(text):
            if claimed <= listed or claimed - listed <= 1:
                continue
            # A count far larger than the list is a summary, not an
            # enumeration ("3,000 rows" above a 10-row sample).
            if claimed > listed * 3:
                continue
            if (claimed, listed) not in out:
                out.append((claimed, listed))
        return out[:3]
    except Exception:
        return []


def count_vs_enumeration_caveat(pairs: list[tuple[int, int]]) -> str:
    """The note appended to an answer that contradicts its own list."""
    if not pairs:
        return ""
    claimed, listed = pairs[0]
    return (
        f"\n\n> ⚠️ Diese Antwort nennt {claimed}, führt aber {listed} "
        "Einträge auf. Eines von beiden ist falsch — bitte nachzählen, "
        "bevor die Zahl weitergegeben wird."
    )


def truncated_output_caveat(counts: list[str], tools: list[str]) -> str:
    """The note appended to an answer counting from a cut-short result."""
    if not counts:
        return ""
    named = ", ".join(counts[:3])
    where = ", ".join(sorted(set(tools))[:3]) or "ein Werkzeug-Ergebnis"
    return (
        "\n\n> ⚠️ Diese Antwort nennt " + named + ", aber die Ausgabe von "
        + where + " wurde in diesem Zug abgeschnitten. Eine Zahl, deren "
        "einzige Quelle abgeschnitten war, ist geschätzt und nicht gezählt "
        "— bitte gegen die vollständige Liste prüfen, bevor sie "
        "weitergegeben wird."
    )


def ambiguous_column_caveat(columns: list[str]) -> str:
    """The note appended to an answer that totalled an unreadable column."""
    if not columns:
        return ""
    named = ", ".join(f"'{c}'" for c in columns[:4])
    return (
        "\n\n> ⚠️ Die Spalte " + named + " ist nicht eindeutig lesbar: "
        "ein Wert wie `8.986` bedeutet 8986 oder 8,986, und nichts in der "
        "Spalte entscheidet das. Die Zahl oben beruht daher auf einer "
        "Annahme und ist nicht gemessen — bitte die Lesart bestätigen, "
        "bevor sie weitergegeben wird."
    )
