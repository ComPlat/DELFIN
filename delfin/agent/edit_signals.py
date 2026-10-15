"""What a turn's file changes did to the measure, and to the object.

Input: every file change the write path journals (the file tools through
``_capture_change`` with old and new text, the shell through
``_capture_raw_change`` with the pre-image bytes), every test run the
turn made (``note_command``), and at the end the user's message and the
model's answer.

Output, per session and per turn:

* :func:`measure_edits` -- changes to the measure: a test's assertion or
  expected value changed or removed, a skip/xfail/ignore added, a
  tolerance loosened, a test or a test file deleted, a grader file (CI
  configuration, test-runner configuration, benchmark or acceptance
  files, snapshots, thresholds) changed, rows removed from a data file;
* :func:`symptom_fixes` -- changes to the object that silence a failure
  instead of removing its cause: a broad except that does nothing, an
  empty catch, ``|| true`` / ``set +e``, an added sleep or retry, a raised
  timeout, a silenced warning, a rerun plugin or flaky marker, and a test
  re-run unchanged until it passed;
* :func:`turn_caveat` -- the one ``[verify]`` line naming them, or "".

Semantics.

* Only what a change ADDS counts: a pattern is counted in the added lines
  and in the removed lines, and only the excess is reported. Code that
  already had it is not flagged.
* A handler that does nothing counts only around a statement the old
  text had (existing code wrapped); new code written with a best-effort
  handler of its own is ordinary. A sleep or retry counts in a small
  insertion (at most four lines) or one that wraps existing lines.
* A file this session created is its own new work: correcting it is no
  measure edit, except a skip/xfail, a loosened tolerance or a widened
  assertion (equality became membership, ``or``, an inequality or a
  type check), which weaken a check whoever wrote it.
* The caveat is due only in a turn that set out to make a failing check
  pass: a test run of this session failed before the change, or the
  user's message asks for a fix (a fix verb, or a failure word in a
  message of at most 400 characters).
* It is withheld when the user's message asks for exactly that change
  (to the test, the expectation, the threshold, the snapshot, the rows;
  or a workaround), and when the answer gives the reason in a sentence
  that names the file.
* A fix turn that changed the object and whose answer names no cause
  gets one line asking for the root cause in one sentence, or the
  statement that it was not found.

Measured (offline replay of 200 recorded sessions, 531 user turns, 48 of
them with file edits, 10 fix turns with edits): 3 caveats. One is a
false alarm (a refactor that moved an existing callback into a
best-effort guard); two name test assertions edited after a red run in
files whose creation the replay cannot see (created through the shell,
which the live write path journals and then exempts). Before the
wrap rule, the session-created exemption and the short-message rule,
the same corpus gave 6 caveats with 4 false alarms. 117 of 245
recorded edit calls in that archive keep no arguments, so the replay
sees about half of the edits. Classification alone, over 363 real edit
pairs from the undo journals (consecutive pre-images of one file): 10
measure findings, every one a real change of the kind it names (skips
added for a missing program, removed tests, a changed expected value, a
changed fixture), and no symptom finding.

Rejected: refusing the edit (a test that encodes a wrong expectation must
be correctable, and the gate cannot tell which one it is); a model call
to classify the diff (a second instrument that grades the first); a
caveat on every test edit (writing and correcting tests is ordinary
work; the alarm case is the change that makes exactly the failing run
pass).

Limits: a change outside the write path's journal (a file another
process writes) is not seen. Patterns are lexical: a skip spelled in a
way the tables below do not know is not counted; a handler around
rewritten (not moved) code is not counted.

Cost: one difflib pass per journalled text change of at most 400 kB; no
prompt text unless a caveat is due.
"""

from __future__ import annotations

import difflib
import os
import re
import threading
from dataclasses import dataclass, field
from pathlib import PurePosixPath
from typing import Iterable, Optional

_LOCK = threading.Lock()
_MAX_TEXT = 400_000
_MAX_EVENTS = 400


@dataclass
class Signal:
    """One finding: which file, what kind, and the line that shows it."""

    path: str
    kind: str
    detail: str = ""
    seq: int = 0


@dataclass
class _Turn:
    seq: int = 0
    measure: list = field(default_factory=list)
    symptom: list = field(default_factory=list)
    object_changed: bool = False
    red_at: Optional[int] = None
    runs: list = field(default_factory=list)   # (seq, command, red)
    last_change_seq: int = 0


_TURNS: dict[str, _Turn] = {}
#: Paths each session created: its own new work, not a measure it was given.
_CREATED: dict[str, set] = {}
#: Kinds that weaken a check whoever wrote it. Correcting one's own new
#: test is ordinary work (recorded: "my test bug -- I tested membership in
#: an int"); making it skip on this host, accept more values or tolerate
#: more is not (recorded: `== MISSING` became `in (MISSING, OK)` after a
#: red run on one host).
_WEAKENING_EVEN_IN_OWN_WORK = frozenset({
    "skip/xfail added", "tolerance loosened", "assertion widened"})


def _turn(session: str) -> _Turn:
    key = str(session or "")
    t = _TURNS.get(key)
    if t is None:
        t = _TURNS[key] = _Turn()
    return t


def reset_turn(session: str) -> None:
    """Forget the previous turn's changes; a red run of the session stays
    known (the turn may be the one that fixes it)."""
    with _LOCK:
        old = _TURNS.get(str(session or ""))
        fresh = _Turn()
        if old is not None and old.red_at is not None:
            fresh.red_at = 0
        _TURNS[str(session or "")] = fresh


# ---------------------------------------------------------------------------
# what a path is
# ---------------------------------------------------------------------------

_TEST_NAME_RE = re.compile(
    r"(?:^test_.*\.py$|_test\.(?:py|go)$|\.(?:test|spec)\.[cm]?[jt]sx?$|"
    r"^conftest\.py$|Test\.java$|Tests?\.cs$|_spec\.rb$)")
_TEST_DIRS = frozenset({"tests", "test", "__tests__", "spec", "testing"})
_GRADER_NAME_RE = re.compile(
    r"(?:^\.gitlab-ci\.ya?ml$|^tox\.ini$|^pytest\.ini$|^\.coveragerc$|"
    r"^noxfile\.py$|^jest\.config\.[cm]?[jt]s$|^vitest\.config\.[cm]?[jt]s$|"
    r"^\.?(?:flake8|pylintrc)$|^ruff\.toml$|^\.pre-commit-config\.yaml$|"
    r"threshold|baseline.*\.(?:json|ya?ml|toml)$|^expected|golden|"
    r"\.snap$|accept(?:ance)?[_.-]|_accept(?:ance)?\.)", re.IGNORECASE)
_GRADER_DIR_RE = re.compile(
    r"(?:^|/)(?:\.github/workflows|\.circleci|__snapshots__|snapshots|"
    r"benchmarks?|acceptance|accept|golden|expected)(?:/|$)", re.IGNORECASE)
_CONFIG_NAMES = frozenset({"pyproject.toml", "setup.cfg", "package.json"})
_CONFIG_TEST_LINE_RE = re.compile(
    r"\[tool\.(?:pytest|coverage|ruff|mypy)|addopts|testpaths|"
    r"filterwarnings|fail_under|\bomit\b|exclude|--ignore|--deselect|"
    r"\"test\"\s*:|norecursedirs|python_files", re.IGNORECASE)
_DATA_EXT = frozenset({".csv", ".tsv", ".jsonl", ".dat"})
_CODE_EXT = frozenset({".py", ".js", ".jsx", ".ts", ".tsx", ".mjs", ".cjs",
                       ".go", ".rs", ".java", ".kt", ".rb", ".cs", ".c",
                       ".cc", ".cpp", ".h", ".hpp", ".jl", ".r", ".sh", ""})


def role(path: str) -> str:
    """``test``, ``grader``, ``data``, ``config`` or ``object``."""
    p = PurePosixPath(str(path or "").replace("\\", "/"))
    name = p.name
    parts = set(p.parts[:-1])
    if _TEST_NAME_RE.search(name) or parts & _TEST_DIRS:
        return "test"
    if _GRADER_NAME_RE.search(name) or _GRADER_DIR_RE.search(str(p.parent) + "/"):
        return "grader"
    if name in _CONFIG_NAMES:
        return "config"
    if p.suffix.lower() in _DATA_EXT:
        return "data"
    return "object"


# ---------------------------------------------------------------------------
# what a change adds
# ---------------------------------------------------------------------------

_SKIP_RE = re.compile(
    r"@pytest\.mark\.(?:skip|skipif|xfail)\b|\bpytest\.(?:skip|xfail)\(|"
    r"@unittest\.(?:skip|expectedFailure)|\bself\.skipTest\(|"
    r"\b(?:it|test|describe)\.(?:skip|todo)\(|\bx(?:it|describe|test)\(|"
    r"#\[ignore\]|\bt\.Skip(?:Now|f)?\(|@Disabled\b|@Ignore\b|"
    r"\bskip\s*=\s*True\b")
_ASSERT_RE = re.compile(
    r"^\s*(?:assert\b|self\.assert\w+\(|expect\(|assert_eq!|assert!|"
    r"assert_\w+\(|np\.testing\.assert|pd\.testing\.assert|"
    r"\w+\.should\b|t\.(?:Error|Fatal)f?\()")
_TEST_DEF_RE = re.compile(
    r"^\s*(?:async\s+)?def\s+(test_\w+)|^\s*(?:it|test)\(\s*['\"]([^'\"]+)|"
    r"^\s*func\s+(Test\w+)\(|^\s*#\[test\]")
_TOL_RE = re.compile(
    r"\b(rel|abs|atol|rtol|tol|tolerance|delta|places|eps|epsilon|"
    r"threshold|max_diff|max_error|decimal)\s*[=:]\s*"
    r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", re.IGNORECASE)
_TIMEOUT_RE = re.compile(
    r"\b(timeout|time_limit|deadline|max_wait|wait_s|timeout_s|"
    r"timeout_seconds)\s*[=:]\s*(\d+(?:\.\d+)?)", re.IGNORECASE)
_SWALLOW_RES = (
    ("broad except that does nothing",
     re.compile(r"^\s*except(?:\s+(?:Exception|BaseException))?(?:\s+as\s+\w+)?"
                r"\s*:\s*(?:pass|\.\.\.|continue|return(?:\s+None)?)?\s*(?:#.*)?$")),
    ("empty catch", re.compile(r"catch\s*(?:\([^)]*\))?\s*\{\s*\}")),
    ("`|| true` on a command", re.compile(r"\|\|\s*(?:true|:)\s*(?:$|[;)#])")),
    ("`set +e`", re.compile(r"^\s*set\s+\+e\b")),
    ("continue-on-error / allow_failure",
     re.compile(r"\b(?:continue-on-error|allow_failure)\s*:\s*true\b",
                re.IGNORECASE)),
)
_DOES_NOTHING_RE = re.compile(r"^\s*(?:pass|\.\.\.|continue|return(?:\s+None)?)"
                              r"\s*(?:#.*)?$")
_RETRY_RE = re.compile(
    r"\btime\.sleep\(|\basyncio\.sleep\(|^\s*sleep\s+\d|\bsetTimeout\(|"
    r"@retry\b|\btenacity\b|\bbackoff\.|\bfor\s+(?:attempt|retry|_?try|tries)"
    r"\w*\s+in\s+range\(|\bmax_retries\s*=|\bretries\s*=", re.IGNORECASE)
_SILENCE_RE = re.compile(
    r"warnings\.(?:filterwarnings|simplefilter)\(\s*['\"]ignore|"
    r"#\s*noqa\b|#\s*type:\s*ignore|eslint-disable|@SuppressWarnings|"
    r"#\s*pragma:\s*no\s*cover|-W\s*ignore|logging\.disable\(|"
    r"2>\s*/dev/null")
_RERUN_RE = re.compile(
    r"--reruns\b|@pytest\.mark\.flaky|@flaky\b|pytest-rerunfailures|"
    r"\bjest\.retryTimes\(|\bretries\s*:\s*\d")


def _lines(text: Optional[str]) -> list[str]:
    return (text or "")[:_MAX_TEXT].splitlines()


def _diff(old: Optional[str], new: Optional[str]) -> tuple[list[str], list[str]]:
    """(removed lines, added lines) of a line diff."""
    a, b = _lines(old), _lines(new)
    removed, added = [], []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(
            None, a, b, autojunk=False).get_opcodes():
        if tag in ("replace", "delete"):
            removed.extend(a[i1:i2])
        if tag in ("replace", "insert"):
            added.extend(b[j1:j2])
    return removed, added


def _excess(rx: re.Pattern, removed: list[str], added: list[str]) -> list[str]:
    """Added lines matching *rx* beyond the removed lines that match it."""
    hits = [ln for ln in added if rx.search(ln)]
    gone = sum(1 for ln in removed if rx.search(ln))
    return hits[gone:] if len(hits) > gone else []


def _short(line: str) -> str:
    return " ".join(line.split())[:90]


def _numbers(rx: re.Pattern, lines: list[str]) -> dict:
    out: dict = {}
    for ln in lines:
        for m in rx.finditer(ln):
            try:
                out.setdefault(m.group(1).lower(), []).append(float(m.group(2)))
            except ValueError:
                pass
    return out


def _loosened(removed: list[str], added: list[str]) -> list[str]:
    """Tolerance keywords whose value moved in the loose direction: up,
    except ``places``/``decimal`` (digits compared), where down is loose."""
    before, after = _numbers(_TOL_RE, removed), _numbers(_TOL_RE, added)
    out = []
    for key, olds in before.items():
        news = after.get(key)
        if not news:
            continue
        o, n = max(olds), max(news)
        loose = n < o if key in ("places", "decimal") else n > o
        if loose:
            out.append(f"{key} {o:g} -> {n:g}")
    return out


def _raised_timeouts(removed: list[str], added: list[str]) -> list[str]:
    before, after = _numbers(_TIMEOUT_RE, removed), _numbers(_TIMEOUT_RE, added)
    return [f"{k} {max(v):g} -> {max(after[k]):g}" for k, v in before.items()
            if after.get(k) and max(after[k]) > max(v)]


_WIDE_RE = re.compile(r"\bin\s*[(\[{]|\bor\b|[<>]=?|!=|\bany\(|"
                      r"\bis\s+not\s+None\b|\bisinstance\(")


def _widened(old: list[str], new: list[str]) -> Optional[str]:
    """The first assertion that now accepts more, as "old -> new": a
    membership, an ``or``, an inequality or a type check where the old
    line, its nearest counterpart among the added ones (difflib ratio at
    least 0.6), compared for equality."""
    for o in old:
        if "==" not in o:
            continue
        best = max(new, key=lambda n: difflib.SequenceMatcher(
            None, o, n).ratio(), default="")
        if not best or difflib.SequenceMatcher(None, o, best).ratio() < 0.6:
            continue
        if len(_WIDE_RE.findall(best)) > len(_WIDE_RE.findall(o)):
            return f"{_short(o)} -> {_short(best)}"
    return None


def _test_names(lines: list[str]) -> set:
    out = set()
    for ln in lines:
        m = _TEST_DEF_RE.search(ln)
        if m:
            out.add(next((g for g in m.groups() if g), "#[test]"))
    return out


def _hunks(old: Optional[str], new: Optional[str]) -> list:
    """(removed lines, added lines) per changed hunk."""
    a, b = _lines(old), _lines(new)
    out = []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(
            None, a, b, autojunk=False).get_opcodes():
        if tag != "equal":
            out.append((a[i1:i2], b[j1:j2]))
    return out


def _strip_set(lines: Iterable[str]) -> set:
    return {" ".join(ln.split()) for ln in lines if len(ln.strip()) >= 3}


def _wraps_existing(hunk_removed: list[str], body: list[str]) -> bool:
    """Whether *body* (the code the new handler guards) holds a statement
    the old text already had: existing code was wrapped, not new code
    written with its own handling."""
    old = _strip_set(hunk_removed)
    if not old:
        return False
    for ln in body:
        flat = " ".join(ln.split())
        if flat in old or any(o in flat for o in old if len(o) >= 6):
            return True
    return False


def _swallows(hunks: list) -> list[str]:
    """Handlers that do nothing, added around code the old text had.

    New code written with a best-effort handler of its own is ordinary;
    the symptom fix is the handler put around a statement that was there
    before -- the call that failed. Measured on recorded sessions: without
    this condition, two of six firings were new helper functions."""
    out = []
    for removed, added in hunks:
        for label, rx in _SWALLOW_RES:
            hits = []
            for i, ln in enumerate(added):
                if not rx.search(ln):
                    continue
                if label.startswith("broad except"):
                    if not re.search(r":\s*(?:pass|\.\.\.|continue|return)", ln):
                        nxt = next((x for x in added[i + 1:i + 3] if x.strip()),
                                   "")
                        if not _DOES_NOTHING_RE.match(nxt):
                            continue
                    start = next((k for k in range(i - 1, -1, -1)
                                  if re.match(r"^\s*try\s*:", added[k])), None)
                    body = added[start + 1:i] if start is not None else []
                    if start is not None and not _wraps_existing(removed, body):
                        continue
                elif label == "empty catch" or label.startswith("`|| true`"):
                    if not _wraps_existing(removed, [rx.sub("", ln)]):
                        continue
                hits.append(ln)
            gone = sum(1 for ln in removed if rx.search(ln))
            for ln in hits[gone:]:
                out.append(f"{label}: {_short(ln)}")
    return out


def _retries(hunks: list) -> list[str]:
    """A sleep or retry put into existing code: a small insertion (at most
    four lines) or one that re-indents lines the old text had."""
    out = []
    for removed, added in hunks:
        hits = _excess(_RETRY_RE, removed, added)
        if not hits:
            continue
        if len(added) <= 4 or _wraps_existing(removed, added):
            out.extend(_short(ln) for ln in hits)
    return out


def classify(path: str, old: Optional[str], new: Optional[str], *,
             created: bool = False, deleted: bool = False) -> dict:
    """The findings of one change: ``{"measure": [(kind, detail)],
    "symptom": [(kind, detail)], "role": ...}``."""
    r = role(path)
    measure: list = []
    symptom: list = []
    if deleted:
        if r == "test" and _test_names(_lines(old)):
            measure.append(("test file deleted", ""))
        elif r == "grader":
            measure.append(("grader file deleted", ""))
        return {"role": r, "measure": measure, "symptom": symptom}
    removed, added = _diff(old, new)
    if created:
        return {"role": r, "measure": [], "symptom": []}
    if r == "test" and PurePosixPath(path).suffix.lower() not in _CODE_EXT:
        # A snapshot, a golden file or a fixture under tests/: its content
        # is the expectation, or the input the expectation was made for.
        changed = [ln for ln in removed + added if ln.strip()]
        if changed:
            measure.append(("test data or snapshot changed",
                            _short(changed[0])))
    elif r == "test":
        for ln in _excess(_SKIP_RE, removed, added):
            measure.append(("skip/xfail added", _short(ln)))
        old_asserts = [ln for ln in removed if _ASSERT_RE.search(ln)]
        new_asserts = [ln for ln in added if _ASSERT_RE.search(ln)]
        if old_asserts and new_asserts:
            wide = _widened(old_asserts, new_asserts)
            if wide:
                measure.append(("assertion widened", wide))
            else:
                measure.append(("assertion or expected value changed",
                                f"{_short(old_asserts[0])} -> "
                                f"{_short(new_asserts[0])}"))
        elif old_asserts:
            measure.append(("assertion removed", _short(old_asserts[0])))
        gone = _test_names(removed) - _test_names(added)
        for name in sorted(gone)[:3]:
            measure.append(("test removed", name))
        for item in _loosened(removed, added):
            measure.append(("tolerance loosened", item))
        if any("pytest.raises" in ln for ln in removed) and not any(
                "pytest.raises" in ln for ln in added):
            measure.append(("expected exception no longer checked", ""))
    elif r == "grader" or (r == "config" and any(
            _CONFIG_TEST_LINE_RE.search(ln) for ln in removed + added)):
        changed = [ln for ln in added + removed if ln.strip()]
        if changed:
            measure.append(("grader changed", _short(changed[0])))
        for item in _loosened(removed, added):
            measure.append(("threshold loosened", item))
    elif r == "data":
        lost = len([ln for ln in removed if ln.strip()]) - len(
            [ln for ln in added if ln.strip()])
        if lost > 0:
            measure.append(("rows removed", f"{lost} row(s)"))
    hunks = _hunks(old, new)
    for item in _swallows(hunks):
        symptom.append(("failure swallowed", item))
    if r in ("object", "config", "grader"):
        for ln in _retries(hunks):
            symptom.append(("sleep or retry added", ln))
        for item in _raised_timeouts(removed, added):
            symptom.append(("timeout raised", item))
        for ln in _excess(_SILENCE_RE, removed, added):
            symptom.append(("warning or error silenced", _short(ln)))
    for ln in _excess(_RERUN_RE, removed, added):
        symptom.append(("rerun of failing tests added", _short(ln)))
    return {"role": r, "measure": measure, "symptom": symptom}


# ---------------------------------------------------------------------------
# recording
# ---------------------------------------------------------------------------

def note_change(session: str, path: str, old: Optional[str],
                new: Optional[str], *, created: bool = False,
                deleted: bool = False) -> dict:
    """Classify and record one journalled change. Never raises; returns
    the classification (empty on failure)."""
    try:
        if old is not None and new is not None and old == new:
            return {}
        if created is False and old is None and not deleted:
            created = True
        found = classify(path, old, new, created=created, deleted=deleted)
    except Exception:
        return {}
    with _LOCK:
        t = _turn(session)
        t.seq += 1
        t.last_change_seq = t.seq
        own = _CREATED.setdefault(str(session or ""), set())
        if created:
            own.add(os.path.normpath(str(path)))
        mine = os.path.normpath(str(path)) in own
        if found.get("role") == "object" or found.get("role") == "config":
            t.object_changed = True
        for kind, detail in (found.get("measure") or ()):
            if mine and kind not in _WEAKENING_EVEN_IN_OWN_WORK:
                continue
            if len(t.measure) < _MAX_EVENTS:
                t.measure.append(Signal(path, kind, detail, t.seq))
        for kind, detail in found.get("symptom") or ():
            if len(t.symptom) < _MAX_EVENTS:
                t.symptom.append(Signal(path, kind, detail, t.seq))
    return found


def note_command(session: str, command: str, output: str) -> None:
    """Record a test run and whether it was red. A test command re-run
    unchanged, red then green with no file change in between, is a rerun
    until green. Never raises."""
    try:
        from . import counter_test as _ct
        cmd = " ".join(str(command or "").split())
        if not _ct._is_test_command(cmd):
            return
        passed, failed, errors = _ct._counts(output or "")
        red = failed > 0
        if not red and not passed:
            return
    except Exception:
        return
    with _LOCK:
        t = _turn(session)
        t.seq += 1
        if red:
            t.red_at = t.seq
        elif t.runs:
            last = next((r for r in reversed(t.runs) if r[1] == cmd), None)
            if last and last[2] and last[0] > t.last_change_seq:
                t.symptom.append(Signal(cmd[:80], "flaky test re-run until "
                                        "it passed", "", t.seq))
        t.runs.append((t.seq, cmd, red))
        if len(t.runs) > 60:
            del t.runs[:-60]


def note_red_check(session: str) -> None:
    """Record that a check failed by some other route (CI, a benchmark)."""
    with _LOCK:
        t = _turn(session)
        t.seq += 1
        t.red_at = t.seq


def measure_edits(session: str, *, after_red_only: bool = False) -> list:
    """The turn's measure edits; with *after_red_only*, those made after
    a red run (a red run of an earlier turn counts as before all)."""
    with _LOCK:
        t = _TURNS.get(str(session or ""))
        if t is None:
            return []
        if after_red_only:
            if t.red_at is None:
                return []
            return [s for s in t.measure if s.seq > t.red_at]
        return list(t.measure)


def symptom_fixes(session: str) -> list:
    with _LOCK:
        t = _TURNS.get(str(session or ""))
        return list(t.symptom) if t else []


def red_seen(session: str) -> bool:
    with _LOCK:
        t = _TURNS.get(str(session or ""))
        return bool(t and t.red_at is not None)


def object_changed(session: str) -> bool:
    with _LOCK:
        t = _TURNS.get(str(session or ""))
        return bool(t and t.object_changed)


# ---------------------------------------------------------------------------
# the turn's caveat
# ---------------------------------------------------------------------------

_FIX_VERB_RE = re.compile(
    r"\b(?:fix|fixe|fixen|repair|debug|resolve|make\s+(?:\w+\s+){0,2}"
    r"(?:pass|green|work)|behebe?n?|beheb\w*|reparier\w*|grün\s+machen)\b",
    re.IGNORECASE)
_FAILURE_WORD_RE = re.compile(
    r"\b(?:broken|fails?|failing|failed|failure|errors?|crash\w*|bugs?|"
    r"regression|kaputt|fehler\w*|schl[äa]gt\s+fehl|fehlschl\w*|"
    r"läuft\s+nicht|stürzt\w*)\b", re.IGNORECASE)
#: A failure word alone asks for a fix only in a short message ("the test
#: fails"); a long brief mentions failures as its subject matter. Measured
#: on recorded briefs: a 1,561-character build brief matched "fail" twice
#: and "error" once without asking for any fix.
_SHORT_REQUEST = 400


def _asks_for_fix(text: str) -> bool:
    t = text or ""
    return bool(_FIX_VERB_RE.search(t)) or (
        len(t) <= _SHORT_REQUEST and bool(_FAILURE_WORD_RE.search(t)))


_ASKS_MEASURE_EDIT_RE = re.compile(
    r"\b(?:update|change|adjust|relax|loosen|tighten|rewrite|"
    r"edit|regenerate|refresh|skip|remove|delete|drop|pass\w*\s+an|"
    r"anpass\w*|änder\w*|aktualisier\w*|lösch\w*|entfern\w*)"
    r"\b[^.!?\n]{0,60}\b(?:tests?|expect\w*|assert\w*|threshold\w*|"
    r"toleranc\w*|snapshots?|golden|baseline|ci|workflow|benchmark|rows?|"
    r"outliers?|schwell\w*|toleranz\w*|erwart\w*|zeilen|ausreißer)\b|"
    # "fix the tests" asks for them to pass; "fix the expected value"
    # asks for the expectation itself to change.
    r"\b(?:fix|correct|korrigier\w*)\b[^.!?\n]{0,30}\b(?:expect\w*|"
    r"assert\w*|threshold\w*|toleranc\w*|snapshots?|golden|erwart\w*|"
    r"schwell\w*)\b|"
    r"\b(?:the\s+)?tests?\s+(?:is|are|was|were)\s+(?:wrong|outdated|"
    r"incorrect|broken)\b|\btests?\s+(?:ist|sind|war)\s+(?:falsch|veraltet)\b",
    re.IGNORECASE)
_ASKS_WORKAROUND_RE = re.compile(
    r"\bwork-?around\b|\bquick\s+fix\b|\bhot-?fix\b|\btemporar\w*|\bfor\s+now\b|"
    r"\bstopgap\b|\bvorübergehend\w*|\bübergangs\w*|\bnotlösung\w*|"
    r"\bprovisorisch\w*|\b(?:add|with|mit)\s+(?:a\s+)?retr(?:y|ies)\b|"
    r"\b(?:raise|increase|erhöh\w*)\s+(?:the\s+)?timeout\b|"
    r"\bsilence\b|\bsuppress\w*|\bunterdrück\w*", re.IGNORECASE)
_REASON_RE = re.compile(
    r"\bbecause\b|\bsince\b|\bas\s+the\b|\bwas\s+(?:wrong|incorrect|outdated|"
    r"stale|based\s+on)|\bdue\s+to\b|\bcaused\s+by\b|\broot\s+cause\b|"
    r"\bthe\s+(?:cause|reason)\b|\bspec(?:ification)?\b|\bweil\b|\bda\s+|"
    r"\bwar\s+(?:falsch|veraltet)|\bursache\b|\bgrund\b|\bliegt\s+an\b|"
    r"\blag\s+an\b", re.IGNORECASE)
_ROOT_CAUSE_RE = re.compile(
    r"\broot\s+cause\b|\bcaused\s+by\b|\bbecause\b|\bdue\s+to\b|"
    r"\bthe\s+(?:cause|reason|problem|issue|bug|defect|error)\s+"
    r"(?:was|is|lay|lies)\b|\b(?:cause|reason)\s*:|\bwas\s+(?:not|never|"
    r"missing|wrong|reading|using|called|passing|returning)\b|"
    r"\bdid\s+not\b|\bdidn't\b|\bwasn't\b|\bwhich\s+(?:made|meant|caused)\b|"
    r"\bnot\s+found\b|\bunknown\b|\bunclear\b|"
    r"\bursache\b|\bweil\b|\bliegt\s+an\b|\blag\s+an\b|\bdadurch\b|"
    r"\b(?:das|der)\s+(?:problem|fehler)\s+(?:war|ist|lag)\b|"
    r"\bnicht\s+gefunden\b|\bunklar\b", re.IGNORECASE)
_SENT_RE = re.compile(r"(?<=[.!?])\s+|\n+")


def _sentences(text: str) -> list[str]:
    return [s for s in _SENT_RE.split(text or "") if s.strip()]


def _reasoned(answer: str, path: str, extra: Iterable[str] = ()) -> bool:
    """The answer gives a reason in a sentence that names the file, or one
    of *extra* words."""
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
    names = [n for n in (path, base, stem, *extra) if n]
    for s in _sentences(answer):
        low = s.lower()
        if any(n.lower() in low for n in names) and _REASON_RE.search(s):
            return True
    return False


def fix_turn(user_message: str, session: str) -> bool:
    """The turn set out to make a failing check pass."""
    return red_seen(session) or _asks_for_fix(user_message)


def turn_caveat(session: str, user_message: str, answer: str) -> str:
    """The one ``[verify]`` caveat for the turn, or "". Never raises."""
    try:
        if not fix_turn(user_message, session):
            return ""
        parts: list[str] = []
        measure = []
        if not _ASKS_MEASURE_EDIT_RE.search(user_message or ""):
            seen = set()
            for s in measure_edits(session):
                key = (s.path, s.kind)
                if key in seen:
                    continue
                seen.add(key)
                if _reasoned(answer, s.path, ("test", "expected", "threshold",
                                               "tolerance", "erwart")):
                    continue
                measure.append(s)
        if measure:
            items = "; ".join(
                f"{_rel(s.path)}: {s.kind}" + (f" ({s.detail})" if s.detail
                                               else "")
                for s in measure[:4])
            parts.append(
                "the check was changed in a turn that set out to make it "
                f"pass — {items}. Say why the old check was wrong, or fix "
                "the code it checks instead")
        symptoms = []
        if not _ASKS_WORKAROUND_RE.search(user_message or ""):
            seen = set()
            for s in symptom_fixes(session):
                key = (s.path, s.kind)
                if key in seen:
                    continue
                seen.add(key)
                if _ROOT_CAUSE_RE.search(answer or "") and _reasoned(
                        answer, s.path, ()):
                    continue
                symptoms.append(s)
        if symptoms:
            items = "; ".join(
                f"{_rel(s.path)}: {s.kind}" + (f" ({s.detail})" if s.detail
                                               else "")
                for s in symptoms[:4])
            parts.append(
                f"a failure is silenced rather than removed — {items}. Name "
                "the cause it hides, or say why this is the right fix there "
                "and that it is a workaround")
        if (object_changed(session) and _asks_for_fix(user_message)
                and not _ROOT_CAUSE_RE.search(answer or "")):
            parts.append("this fix names no root cause: state it in one "
                         "sentence, or say it was not found")
        if not parts:
            return ""
        return "\n\n[verify] " + " ".join(
            p[0].upper() + p[1:] + "." for p in parts)
    except Exception:
        return ""


def _rel(path: str) -> str:
    try:
        rel = os.path.relpath(path)
        return rel if not rel.startswith("..") else path
    except ValueError:
        return path


__all__ = ["Signal", "role", "classify", "note_change", "note_command",
           "note_red_check", "reset_turn", "measure_edits", "symptom_fixes",
           "red_seen", "object_changed", "fix_turn", "turn_caveat"]
