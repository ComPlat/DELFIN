"""ORCA simple-input checks DELFIN runs before an input reaches SLURM.

Consumed by ``delfin.api.validate_orca_input`` (delfin/api.py, the
function ``delfin.agent.calc_policy.check_orca_text`` calls before every
submit and the ops server exposes). Every finding is returned as a plain
``Issue`` here and mapped to ``ValidationIssue`` by the caller.

Sources of truth:
- keyword/block vocabulary: ``delfin.orca_knowledge`` (OPI 2.0.0 catalog
  + the words ORCA's own parser accepted),
- ``delfin.orca_knowledge.probe.probe_input`` when the static check
  cannot decide.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass

from delfin.orca_knowledge.probe import probe_input


@dataclass
class Issue:  # noqa: D401 - mirrors api.ValidationIssue fields
    severity: str   # "error" | "warning" | "info"
    code: str
    message: str
    suggestion: str = ""


#: Blocks ORCA knows but the OPI catalog does not list (block NAMES only;
#: their options are still checked against _GAP_BLOCK_OPTIONS below).
#: All are real ORCA blocks DELFIN itself writes. %pal/%maxcore are
#: range-checked in delfin/api.py's validate_orca_input %pal/%maxcore
#: branches; %scan is parsed by ORCA's own scan grammar and checked in
#: check_scan_block here.
_CATALOG_GAP_BLOCKS = {"maxcore", "pal", "scan"}

#: Options of the gap blocks ORCA's manual documents (ORCA 6 manual,
#: chapter 7: %pal and %maxcore are utility blocks with a fixed, tiny
#: option set -- they are not in the OPI typed catalog because OPI does
#: not model utility blocks).
_GAP_BLOCK_OPTIONS = {
    "pal": {"nprocs", "end"},
    "maxcore": {"end"},
}

_SOLVENT_FORMS = ("CPCM", "SMD", "CPCMC")


def _normalize_line(text: str) -> str:
    """Fullwidth look-alikes to ASCII before scanning (U+FF01 -> '!')."""
    return text.replace("\uff01", "!")


def _closest(words, word):
    """The known keyword closest to ``word`` (case-blind), else "".

    Token-aware: an invented keyword is often a real keyword's tokens in
    the wrong order (OptTight vs TightOpt). The CATALOG word's own best
    split (tight|opt matched against the full invented word) picks the
    closest real keyword. The suggestion carries the CATALOG's casing
    ('tightopt', straight from opi_keywords.json) -- no hand-written
    title-case map, it would drift from the catalog (QS ruling); ORCA
    accepts keywords case-insensitively, so the catalog form is exactly
    what the user can type. Documented in check_keyword_line too.
    """
    from difflib import SequenceMatcher
    pool = sorted(words)
    w = word.lower()

    def _split_ratio(k, w):
        # split the CATALOG word; each half is matched against the FULL
        # invented word. The winning split gives the token boundaries of
        # the real keyword (tight|opt).
        best = 0.0
        for i in range(1, len(k)):
            h1, h2 = k[:i], k[i:]
            score = min(SequenceMatcher(None, h1, w).ratio(),
                        SequenceMatcher(None, h2, w).ratio())
            if score > best:
                best = score
        return best

    ranked = []
    for k in pool:
        full = SequenceMatcher(None, w, k).ratio()
        ranked.append((-(full + _split_ratio(k, w)), k))
    ranked.sort()
    if ranked and -ranked[0][0] >= 0.5:
        return ranked[0][1]
    return ""


def _known_solvent(name: str) -> bool:
    """Whether ORCA knows ``name`` as a solvent (case/space-blind).

    Thin wrapper over orca_knowledge's own data (QS ruling: no patch to
    orca_knowledge): the same normalization as
    delfin/orca_knowledge/__init__.py:119 -- lower + stripped spaces --
    over the catalog's solvent list. Preferred source: solvents_from_opi
    (the OPI enum, what the pinned OPI knows); without the opi package
    installed it falls back to the same opi_keywords.json "solvents" key
    _keywords reads, so the check works everywhere DELFIN runs.
    """
    from delfin.orca_knowledge import _catalog
    try:
        try:
            from delfin.orca_knowledge import solvents_from_opi
            source = solvents_from_opi()
        except ImportError:
            source = _catalog().get("solvents", [])
        known = {s.lower().replace(" ", "") for s in source}
    except Exception:  # catalog unreadable: decide nothing, not wrongly
        return True
    return str(name or "").strip().lower().replace(" ", "") in known


def _solvent_pool() -> list[str]:
    """Solvent names for suggestions, same source order as _known_solvent."""
    try:
        from delfin.orca_knowledge import solvents_from_opi
        return solvents_from_opi()
    except ImportError:
        from delfin.orca_knowledge import _catalog
        return list(_catalog().get("solvents", []))


def check_keyword_line(text: str) -> list[Issue]:
    """Check every word on the ``!`` line(s) against ORCA's vocabulary."""
    from delfin.orca_knowledge import known_keyword
    issues: list[Issue] = []
    # A line that LOOKS like a keyword line but uses a non-ASCII
    # look-alike of `!` (U+FF01 etc.) is silently ignored by ORCA's
    # parser -- the whole line's keywords are dropped. Flag the line
    # BEFORE normalisation, naming the look-alike character (QS ruling);
    # normalize-then-scan alone can never see it, because no ASCII '!'
    # exists on the line.
    for line in text.splitlines():
        stripped = line.lstrip()
        if not stripped or stripped[0].isascii():
            continue
        import unicodedata
        norm0 = unicodedata.normalize("NFKD", stripped[0])
        if norm0 == "!":
            issues.append(Issue(
                severity="error", code="lookalike_keyword_line",
                message=(f"The line uses {stripped[0]!r} "
                         f"(U+{ord(stripped[0]):04X}) instead of the "
                         "ASCII `!` keyword marker; ORCA ignores this "
                         "line entirely, so its keywords are not set."),
                suggestion="Replace it with the ASCII `!` keyword marker.",
            ))
            break
    text = _normalize_line(text)
    for line in re.findall(r"^\s*!(.*)$", text, re.MULTILINE):
        for word in line.split():
            w = word.strip()
            if not w:
                continue
            if w.startswith("%"):
                issues.append(Issue(
                    severity="error", code="block_on_keyword_line",
                    message=f"`{w}` is a %-block written on the `!` line.",
                    suggestion="Blocks belong on their own line: "
                               f"`%{w.lstrip('%')} ... end`.",
                ))
                continue
            m = re.match(r"([A-Za-z]+)\(([^)]*)\)$", w)
            if m and m.group(1).upper() in _SOLVENT_FORMS:
                solvent = m.group(2)
                if not _known_solvent(solvent):
                    close = _closest(_solvent_pool(), solvent)
                    close = close.capitalize() if close else close
                    issues.append(Issue(
                        severity="error", code="unknown_solvent",
                        message=(f"`{m.group(1)}({solvent})` names a "
                                 "solvent ORCA does not know."),
                        suggestion=(f"Did you mean `{close}`?" if close
                                    else "Check the ORCA manual "
                                         "appendix for solvent names."),
                    ))
                continue
            if known_keyword(w):
                continue
            from delfin.orca_knowledge import _keywords
            kws, _solvents, aliases, _blocks = _keywords()
            close = _closest(kws | aliases, w)
            issues.append(Issue(
                severity="error", code="unknown_keyword",
                message=f"`{w}` is not an ORCA keyword.",
                suggestion=(f"Did you mean `{close}`?" if close
                            else "Check the ORCA manual chapter 6."),
            ))
    return issues


def check_keyword_conflicts(text: str) -> list[Issue]:
    """Keyword combinations ORCA's manual names as exclusive."""
    kw_words = []
    for line in re.findall(r"^\s*!(.*)$", _normalize_line(text), re.MULTILINE):
        kw_words += line.split()
    conflicts = [("NUMFREQ", "FREQ")]
    # The message quotes each token AS THE USER WROTE IT (as-written form
    # wins, uppercase only as fallback) -- a user reading only the message
    # must see their own words back (QS ruling).
    upper_words = [w.upper() for w in kw_words]
    def _as_written(key):
        for w in kw_words:
            if w.upper() == key:
                return w
        return key
    issues: list[Issue] = []
    for a, b in conflicts:
        if a in upper_words and b in upper_words:
            ca, cb = _as_written(a), _as_written(b)
            issues.append(Issue(
                severity="error", code="conflicting_keywords",
                message=f"{ca} and {cb} exclude each other.",
                suggestion=f"Pick one: {ca} for a numerical Hessian, "
                           f"{cb} for the analytic one.",
            ))
    return issues


def check_ri_aux_basis(text: str) -> list[Issue]:
    """RI without a matching auxiliary (/C or /J) basis is an error.

    ORCA 6 manual, section on RI approximations: RIJCOSX/RI need an
    auxiliary basis (def2-SVP/C etc.); without one ORCA aborts.
    """
    kw_words = []
    for line in re.findall(r"^\s*!(.*)$", _normalize_line(text), re.MULTILINE):
        kw_words += [w.upper() for w in line.split()]
    if "RI" not in kw_words:
        return []
    has_aux = any(w.endswith("/C") or w.endswith("/J") or w.endswith("/JK")
                  for w in kw_words)
    if has_aux:
        return []
    return [Issue(
        severity="error", code="ri_without_aux_basis",
        message="RI needs an auxiliary basis (/C or /J) in the keyword line.",
        suggestion="Add e.g. def2-SVP/C matching your orbital basis.",
    )]


def check_blocks(text: str) -> list[Issue]:
    """Check ``%block ... end`` names and options against the catalog."""
    from delfin.orca_knowledge import known_block_option
    issues: list[Issue] = []
    text = _normalize_line(text)
    for m in re.finditer(r"^\s*%\s*([A-Za-z0-9_]+)(.*)$", text, re.MULTILINE):
        block = m.group(1)
        block_low = block.lower()
        rest = m.group(2)
        # Inline utility blocks (%maxcore N, %pal nprocs N end) carry no
        # multi-line body: their options are the rest of the SAME line.
        # A search to the next `end` would swallow the following text.
        if block_low in ("maxcore",):
            body = rest
        elif re.search(r"\bend\b", rest, re.IGNORECASE):
            body = rest[:re.search(r"\bend\b", rest, re.IGNORECASE).start()]
        else:
            tail = text[m.end():]
            stop = re.search(r"^\s*end\s*$", tail, re.MULTILINE
                             | re.IGNORECASE)
            body = tail[:stop.start()] if stop else tail
        if block_low not in _CATALOG_GAP_BLOCKS and \
                not known_block_option(block):
            issues.append(Issue(
                severity="error", code="unknown_block",
                message=f"`%{block}` is not an ORCA block.",
                suggestion="Check the ORCA manual chapter 7 for block names.",
            ))
            continue
        # options until the block's end, skipping nested % sub-blocks and
        # pure data lines (numbers, coordinates, point charges).
        for line in body.splitlines():
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            if re.match(r"^[-+0-9.eE\s]+$", line):
                continue
            token = re.match(r"([A-Za-z0-9_]+)", line)
            if not token:
                continue
            option = token.group(1)
            if option.lower() == "end":
                break
            known = (known_block_option(block, option)
                     if block_low not in _CATALOG_GAP_BLOCKS
                     else option.lower() in _GAP_BLOCK_OPTIONS.get(
                         block_low, frozenset()))
            if not known:
                # ORCA warns and continues on unknown block options; the
                # parser probe (phase 3) decides inputs this cannot.
                issues.append(Issue(
                    severity="warning", code="unknown_block_option",
                    message=(f"`{option}` is not a known option of "
                             f"`%{block}`."),
                    suggestion="Verify against the ORCA manual or use "
                               "the parser probe for this input.",
                ))
    return issues


def check_scan_block(text: str) -> list[Issue]:
    """Scan-block checks (wave-23 incidents).

    - a constraint named twice (duplicate atoms) is an error;
    - a start value far from the geometry's current value is a warning;
    - an end point duplicating the start across the +-180 wrap is an
      error.
    """
    text = _normalize_line(text)
    m_scan = re.search(r"^\s*%scan\b(.*)$", text, re.MULTILINE | re.IGNORECASE)
    if not m_scan:
        return []
    tail = text[m_scan.end():]
    stop = re.search(r"^\s*end\s*$", tail, re.MULTILINE | re.IGNORECASE)
    body = tail[:stop.start()] if stop else tail

    atoms = _geometry_atoms(text)

    issues: list[Issue] = []
    seen: dict[tuple, int] = {}
    for entry in re.findall(r"\{([^}]*)\}\s*[SA]?", body):
        tokens = entry.split()
        if not tokens or tokens[0].upper() not in ("B", "A", "D", "S"):
            continue
        kind = tokens[0].upper()
        idx: list[int] = []
        rest = tokens[1:]
        while rest and re.fullmatch(r"\d+", rest[0]) and len(idx) < 4:
            idx.append(int(rest[0]))
            rest = rest[1:]
        if len(idx) < 2 or len(rest) < 3:
            continue
        start, end = float(rest[0]), float(rest[1])
        key = (kind, tuple(idx))
        if key in seen:
            issues.append(Issue(
                severity="error", code="duplicate_scan_constraint",
                message=(f"Scan constraint {kind} {idx} appears twice "
                         f"(lines for points {seen[key]} and later)."),
                suggestion="Each constraint may be listed once.",
            ))
        else:
            seen[key] = 1
        current = _geometry_value(kind, idx, atoms)
        if current is not None and abs(start - current) > 0.5:
            issues.append(Issue(
                severity="warning", code="scan_start_far_from_geometry",
                message=(f"Scan {kind} {idx} starts at {start} but the "
                         f"input geometry has {current:.2f}."),
                suggestion="Start the scan at the geometry's value.",
            ))
        if kind == "D" and _wraps_to_start(start, end):
            issues.append(Issue(
                severity="error", code="scan_end_wraps_to_start",
                message=(f"Dihedral scan {idx} ends at {end}, which is "
                         f"the start value {start} across the +-180 "
                         "wrap; the last point repeats the first."),
                suggestion="End the scan inside one branch, e.g. "
                           "179 -> 175 instead of 179 -> -179.",
            ))
    return issues


def _wraps_to_start(start: float, end: float) -> bool:
    """Whether ``end`` equals ``start`` through the +-180 dihedral wrap.

    True when the two values are the same angle represented on opposite
    sides of the +-180 branch (|start|, |end| near 180, opposite signs),
    i.e. the scan's last point repeats its first.
    """
    return (abs(abs(start) - 180.0) < 5.0 and abs(abs(end) - 180.0) < 5.0
            and (start > 0) != (end > 0))


def _geometry_atoms(text: str) -> list[tuple]:
    """(element, x, y, z) from the first * xyz block, if readable."""
    m = re.search(r"^\s*\*\s+xyz\s+[-+0-9]+\s+[-+0-9]+\s*$", text,
                  re.MULTILINE | re.IGNORECASE)
    if not m:
        return []
    tail = text[m.end():]
    stop = re.search(r"^\s*\*\s*$", tail, re.MULTILINE)
    out = []
    for line in (tail[:stop.start()] if stop else tail).splitlines():
        parts = line.split()
        if len(parts) >= 4:
            try:
                out.append((parts[0], float(parts[1]), float(parts[2]),
                            float(parts[3])))
            except ValueError:
                return out
    return out


def _geometry_value(kind: str, idx: tuple, atoms):
    """The geometry's bond (B), angle (A) or dihedral (D) value."""
    import math
    if not atoms or len(idx) < 2 or any(i >= len(atoms) for i in idx):
        return None
    p = [atoms[i][1:4] for i in idx]

    def dist(a, b):
        return math.dist(a, b)

    def angle(a, b, c):
        v1 = [b[i] - a[i] for i in range(3)]
        v2 = [c[i] - b[i] for i in range(3)]
        cos = (sum(v1[i] * v2[i] for i in range(3))
               / (math.dist(a, b) * math.dist(b, c) or 1e-12))
        return math.degrees(math.acos(max(-1.0, min(1.0, cos))))

    def dihedral(a, b, c, d):
        b1 = [b[i] - a[i] for i in range(3)]
        b2 = [c[i] - b[i] for i in range(3)]
        b3 = [d[i] - c[i] for i in range(3)]
        n1 = [b1[(i + 1) % 3] * b2[(i + 2) % 3]
              - b1[(i + 2) % 3] * b2[(i + 1) % 3] for i in range(3)]
        n2 = [b2[(i + 1) % 3] * b3[(i + 2) % 3]
              - b2[(i + 2) % 3] * b3[(i + 1) % 3] for i in range(3)]
        n1len = math.sqrt(sum(x * x for x in n1)) or 1e-12
        n2len = math.sqrt(sum(x * x for x in n2)) or 1e-12
        n1 = [x / n1len for x in n1]
        n2 = [x / n2len for x in n2]
        m1 = [n1[(i + 1) % 3] * b2[(i + 2) % 3]
              - n1[(i + 2) % 3] * b2[(i + 1) % 3] for i in range(3)]
        x = sum(n1[i] * n2[i] for i in range(3))
        y = sum(m1[i] * n2[i] for i in range(3))
        return math.degrees(math.atan2(y, x))

    if kind == "B" and len(idx) == 2:
        return dist(p[0], p[1])
    if kind == "A" and len(idx) == 3:
        return angle(p[0], p[1], p[2])
    if kind == "D" and len(idx) == 4:
        return dihedral(p[0], p[1], p[2], p[3])
    return None


def check_all(text: str, use_probe: bool = False) -> list[Issue]:
    """All static simple-input checks of this module.

    With ``use_probe=True`` an input the static checks cannot fully
    decide (their findings include warnings but no errors) is sent to
    ORCA's own parser (``delfin.orca_knowledge.probe.probe_input``) and
    a ``rejected`` verdict becomes an ``orca_parser_rejects`` error.
    Verdicts are cached by input hash, so the same input is never
    probed twice.
    """
    issues = (check_keyword_line(text) + check_keyword_conflicts(text)
              + check_ri_aux_basis(text) + check_blocks(text)
              + check_scan_block(text))
    if (use_probe and issues
            and not any(i.severity == "error" for i in issues)):
        verdict = _probe_verdict(text)
        if verdict["verdict"] == "rejected":
            issues.append(Issue(
                severity="error", code="orca_parser_rejects",
                message="ORCA's own parser rejects this input: "
                        + "; ".join(verdict["messages"])[:400],
            ))
        # "parsed"/"finished": no error. "unavailable"/"timeout"/"refused":
        # the static warnings stand on their own (recorded above); the
        # caller decides whether a timeout is worth blocking a submit.
    return issues


_PROBE_CACHE: dict[str, dict] = {}


def _probe_cache_key(text: str) -> str:
    """The cache key: sha256 of the input text itself."""
    return hashlib.sha256(text.encode()).hexdigest()


def _probe_verdict(text: str) -> dict:
    """probe_input(text) with a module-level cache keyed by input hash."""
    key = _probe_cache_key(text)
    if key not in _PROBE_CACHE:
        _PROBE_CACHE[key] = probe_input(text)
    return _PROBE_CACHE[key]


def clear_probe_cache() -> None:
    """Empty the probe-verdict cache (for tests and long sessions)."""
    _PROBE_CACHE.clear()
