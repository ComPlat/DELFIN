"""test_honesty — a pure scanner that proves a change's tests are honest.

Wave 17 measured what this module now checks mechanically: a feature
behind a switch was only tested with the switch OFF; probe files
``tests/test_zz_*`` were committed; tests loaded ``.gate/px`` copies
and passed only in the builder's tree; two real bugs (an undefined
name, a missing import) stayed invisible because new code sat inside a
bare ``except: pass`` — the tests were green throughout.

``check`` answers one question about a change: do its TESTS claim more
than they measure?  It is PURE — AST on the new code plus the paths and
hunk ranges it is handed; it never executes the code under test, never
touches the filesystem, and raises nothing.

Input shape: a list of change records, one per changed file::

    (path, new_start, new_end, source)

``source`` is the NEW text of the file (``None`` reads nothing — for a
path-only record the scanners that need text stay silent).  Hunk
boundaries say which lines are NEW: an old bare ``except`` in the same
file is outside the range and stays quiet.

Every finding carries ``path``, ``line`` (1-based) and a one-line
``reason`` naming the kind.  All four kinds are DIFF-SCOPED (QS ruling,
wave 18): the record set IS the change — kind 1 judges only handlers
inside the record's hunk range; kinds 2, 3 and 4 judge only files in
the record set.  Pre-existing violations elsewhere are history, not
this change's defect.  Kind 4 counts only tests PRESENT IN the record
set: a finding claims "no test in the change records sets it True /
leaves it off", never an absolute non-existence claim — hand the
change's test files over with the records.
"""

import ast
import os
import re
from dataclasses import dataclass

__all__ = ["HonestyFinding", "check"]


@dataclass(frozen=True)
class HonestyFinding:
    """One dishonesty pattern found in a change's new code or tests."""

    path: str
    line: int
    reason: str

    def message(self) -> str:
        return f"⚠️ test_honesty: {self.path}:{self.line} — {self.reason}"


# A bare ``except:`` with no exception class at all.
_BARE_EXCEPT = "bare"

# ``except Exception:`` with a body that only ``pass``es or only
# ``return None``s — swallows everything while pretending to handle it.
_SWALLOWING_EXCEPTION_BODIES = frozenset({"pass", "return_none"})

# Test-probe naming: committed ``tests/test_zz_*`` is a wave-rule
# violation (probes are uncommitted).
_TEST_ZZ_PATH = re.compile(r"^tests/test_zz_.+\.py$")

# Reading anything under .gate/ from a test: the gate's policy area,
# and a test that reads a .gate copy passes only in the builder's tree
# (the wave-17 failure).
_GATE_REFS = re.compile(r"[\"\']\.gate(/|[\"\'])")

# Settings-switch shape: an UPPER_SNAKE module constant assigned a
# bool literal (True/False), OR an env/settings lookup of a flag.
_SWITCH_ASSIGN = re.compile(
    r"^([A-Z][A-Z0-9_]{2,})\s*=\s*(True|False)\s*(?:#.*)?$", re.MULTILINE)

# How a test turns a switch on or off: attribute assignment in a test
# body or an env-var set in a test function.
_SWITCH_ON_TEST = re.compile(
    r"\.([A-Za-z][A-Za-z0-9_]*)\s*=\s*True"
    r"|setenv?\(\s*['\"]([A-Za-z0-9_]+)['\"]",
    re.MULTILINE)

_SWITCH_OFF_TEST = re.compile(
    r"\.([A-Za-z][A-Za-z0-9_]*)\s*=\s*False"
    r"|delenv\(\s*['\"]([A-Za-z0-9_]+)['\"]",
    re.MULTILINE)

# Only paths that hold product code can introduce a switch.
_SWITCH_OWNER_DIR = re.compile(r"^delfin/")


# ---------------------------------------------------------------------------
# kind 1: swallowing exception handlers around new code
# ---------------------------------------------------------------------------

def source_line_count(tree) -> int:
    """The last line the parsed text reaches (max of any node's end)."""
    last = 0
    for node in ast.walk(tree):
        end = getattr(node, "end_lineno", None)
        if end and end > last:
            last = end
    return last


def _except_findings(path: str, tree, start: int, end: int) -> list:
    findings = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Try):
            continue
        for handler in node.handlers:
            # The record's text is either the file's FULL new text (its
            # line count reaches the hunk end) or a hunk SNIPPET whose
            # first line corresponds to file line `start`.  AST lines
            # count the given text; map them into the file accordingly.
            if end and end > start and source_line_count(tree) < end:
                mapped = handler.lineno + max(start - 1, 0)
            else:
                mapped = handler.lineno
            if not (start <= mapped and (end == 0 or mapped <= end)):
                # Pre-existing handler outside this change's hunks.
                continue
            bare = handler.type is None
            body = handler.body
            swallowing = False
            if bare:
                swallowing = True
            else:
                names = _exception_names(handler.type)
                if ("exception" in names
                        or "baseexception" in names):
                    swallowing = _body_is_swallow(body)
            if not swallowing:
                continue
            if bare:
                reason = ("new bare `except:` swallows every error — "
                          "catch narrowly and make the failure visible "
                          "to a test (wave-17 lesson: this hid an "
                          "undefined name and a missing import)")
            else:
                reason = ("new `except Exception:` only "
                          + _body_label(body)
                          + " — swallows every error; catch narrowly "
                          "and re-raise, log, or let it fail")
            findings.append(HonestyFinding(
                path=path, line=mapped, reason=reason))
    return findings


def _exception_names(type_node) -> list:
    """Exception class(es) an ``except`` clause names, lowercased."""
    if type_node is None:
        return []
    if isinstance(type_node, ast.Name):
        return [type_node.id.lower()]
    if isinstance(type_node, ast.Attribute):
        return [type_node.attr.lower()]
    if isinstance(type_node, ast.Tuple):
        out = []
        for elt in type_node.elts:
            out.extend(_exception_names(elt))
        return out
    return []


def _body_is_swallow(body) -> bool:
    """Is the handler's body just ``pass`` / ``return None``?"""
    label = _body_label(body)
    if label == "pass":
        return True
    if label == "return_none":
        return True
    return False


def _body_label(body) -> str:
    """``pass`` | ``return_none`` | ``other`` for a handler body."""
    if len(body) == 1:
        stmt = body[0]
        if isinstance(stmt, ast.Pass):
            return "pass"
        if isinstance(stmt, ast.Return):
            val = stmt.value
            if val is None or (isinstance(val, ast.Constant)
                               and val.value is None):
                return "return_none"
    return "other"


# ---------------------------------------------------------------------------
# kind 2: a committed tests/test_zz_* probe
# ---------------------------------------------------------------------------

def _test_zz_findings(path: str, source: str) -> list:
    if not _TEST_ZZ_PATH.match(path):
        return []
    # Only committed content reaches check() as a record; a record
    # with no source is a path-only mention and not a commit.
    if source is None:
        return []
    return [HonestyFinding(
        path=path, line=1,
        reason=("tests/test_zz_* files are uncommitted probes; this one "
                "is in the change — move it out of the commit "
                "(mv to .gate/ per the wave rules)"))]


# ---------------------------------------------------------------------------
# kind 3: a test that reads from .gate/
# ---------------------------------------------------------------------------

def _gate_reader_findings(path: str, source: str, tree) -> list:
    if not path.startswith("tests/"):
        return []
    if not source or not _GATE_REFS.search(source):
        return []
    line = 1
    for i, ln in enumerate(source.splitlines(), 1):
        if _GATE_REFS.search(ln):
            line = i
            break
    # Prefer the AST line of the open()/read that mentions it.
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") in (
                "open", "read_text", "read_bytes"):
            try:
                src = ast.unparse(node)
            except Exception:
                continue
            if ".gate" in src:
                line = node.lineno
                break
    return [HonestyFinding(
        path=path, line=line,
        reason=("a test reads from .gate/ — that copy exists only in the "
                "builder's tree, so the test passes there and nowhere "
                "else (wave-17 failure); tests drive the real committed "
                "path"))]


# ---------------------------------------------------------------------------
# kind 4: a switch introduced with tests in only one state
# ---------------------------------------------------------------------------

def _switch_findings(path: str, source: str, all_paths) -> list:
    if not source or not _SWITCH_OWNER_DIR.match(path):
        return []
    switches = []
    for m in _SWITCH_ASSIGN.finditer(source):
        line = source.count("\n", 0, m.start()) + 1
        switches.append((m.group(1), line))
    if not switches:
        return []
    findings = []
    test_sources = [src for p, _s, _e, src in all_paths
                    if p.startswith("tests/") and src]
    for name, line in switches:
        low = name.lower()
        on = off = False
        for tsrc in test_sources:
            for m in _SWITCH_ON_TEST.finditer(tsrc):
                hit = (m.group(1) or m.group(2) or "").lower()
                if hit == low or hit.endswith(low) or low.endswith(hit):
                    on = True
            for m in _SWITCH_OFF_TEST.finditer(tsrc):
                hit = (m.group(1) or m.group(2) or "").lower()
                if hit == low or hit.endswith(low) or low.endswith(hit):
                    off = True
        if on and off:
            continue
        if not on and not off:
            reason = (f"switch {name} is introduced but no test sets it "
                      "True and none leaves it off — the feature behind "
                      "a switch is tested with the switch ON and OFF "
                      "(wave-17 failure)")
        elif not off:
            reason = (f"switch {name} is only ever tested ON — add the "
                      "test that leaves it off (wave-17 failure: a "
                      "feature behind a switch was only tested with the "
                      "switch off; both states belong in the suite)")
        else:
            reason = (f"switch {name} is only ever tested OFF — add the "
                      "test that sets it True")
        findings.append(HonestyFinding(path=path, line=line, reason=reason))
    return findings


# ---------------------------------------------------------------------------
# public entry
# ---------------------------------------------------------------------------

def check(diff_or_files) -> list:
    """Scan a change for tests that claim more than they measure.

    ``diff_or_files`` is a list of change records
    ``(path, new_start, new_end, source)``; ``source`` is the file's
    NEW text (``None`` = path-only record).  Returns a list of
    :class:`HonestyFinding` with file:line and a one-line reason.
    Pure: AST + the text given; no execution, no filesystem access.
    """
    records = list(diff_or_files or [])
    findings = []
    for path, start, end, source in records:
        tree = None
        if source and path.endswith(".py"):
            try:
                tree = ast.parse(source)
            except SyntaxError:
                tree = None  # syntax errors are the author's, not ours
        if source:
            findings.extend(_test_zz_findings(path, source))
            findings.extend(_gate_reader_findings(path, source, tree))
        if tree is not None:
            findings.extend(_except_findings(path, tree, start, end))
        if source:
            findings.extend(_switch_findings(path, source, records))
    return findings
