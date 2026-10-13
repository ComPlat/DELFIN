"""Whether a fix's own test was seen failing before it was seen passing.

Input: the files a session changed (change-journal records
``{"path", "ts"}`` or bare paths) and its test runs, rebuilt from the
durable tool trace by :func:`durable_test_runs` and adapted by
:func:`normalize_test_runs`. The live per-turn test ledger cannot serve:
it is replaced at every turn start, so a red run from an earlier turn is
not in it.

Output of :func:`assess`: ``None`` when no product file changed,
``"shown"``, or one sentence starting ``"no control run: "`` that names
what is missing. Pure: no I/O, never raises on malformed records.

Semantics.

* A target is a test node (``tests/test_x.py::test_y``), a whole test
  file run without selection, or a file run under a selection (``-k``,
  ``-m``, ``--deselect``, ``--lf``, ``--sw``), keyed by that selection.
* A control is a target seen red and later green. A red run's
  ``FAILED <node>`` lines name nodes, and then the nodes, not the file,
  are what was red; such a node is controlled by a later green of the
  same node, of its file run without selection, of the whole suite run
  without selection, or of the same selection. A green of another node
  or another selection did not necessarily run it.
* A run counts when a test runner is on the command line (pytest,
  ``python -m pytest``, the gate wrapper) or it went through the
  run_tests tool, and its output reports a count or a non-zero exit. A
  run whose only problems are errors (a collection ImportError, a
  fixture setup error) never ran the test's assertions and is no red.
* When a test file changed, only targets in a changed test file decide.
  A product-only change is covered by any controlled node: there is no
  coverage data to say which test exercises it.
* A test file changed after the target's last red run and before its
  green is not the test that was seen red.

Order uses the trace's float timestamps; a record without one (hand-made
input) has unknown order, and then a red plus a green count. The journal
stamps whole seconds: a stamp ``c`` is a time in ``[c, c+1)``, so a
change is after a red at ``r`` only when ``c > r``.

Known limits.

* A node added to an existing file that never ran is in neither input; a
  controlled older node of the same file stands for it.
* The trace keeps the first 2000 characters of a shell run's output. A
  long red run loses its summary and FAILED lines and is red by its exit
  status alone, naming no node; a long green run whose count line was cut
  is not counted.
* A command line the live ledger's pattern does not recognise
  (``timeout 600 python -m pytest``, an interpreter by full path) is not
  a run here either; the two must agree.
* Output is what the command printed: ``pytest ...; echo "1 failed"`` is
  read as printed.

Rejected: requiring a journalled product change between red and green. A
control against a baseline copy (the red run in a checkout without the
fix) has no change in between and is the stronger control.
"""

from __future__ import annotations

import json
import re
import shlex
from typing import Dict, Iterable, List, Optional, Tuple

# Mirrors api_client._TEST_BASH_CMD_RE, plus the gate wrapper
# (``.delfin/gate-tools/gate TESTS...``).
_TEST_CMD_RE = re.compile(
    r"(?:^|[;&|(])\s*(?:python[\d.]*\s+-m\s+(?:pytest|unittest)\b"
    r"|pytest\b|py\.test\b)")
_GATE_CMD_RE = re.compile(
    r"(?:^|[;&|(])\s*(?:[\w./\\-]*gate-tools/gate|gate)\b")
_PYTEST_PASSED_RE = re.compile(r"(\d+)\s+passed\b")
_PYTEST_FAILED_RE = re.compile(r"(\d+)\s+failed\b")
_PYTEST_ERROR_RE = re.compile(r"(\d+)\s+errors?\b")
_UNITTEST_FAIL_RE = re.compile(
    r"FAILED\s*\((?:failures=(\d+))?(?:,\s*)?(?:errors=(\d+))?\)")
_FAILED_NODE_RE = re.compile(r"^FAILED\s+([\w./\\-]+\.py::\S+)",
                             re.MULTILINE)
_COLLECTION_ERROR_RE = re.compile(
    r"error during collection|^ERROR collecting\b", re.MULTILINE)
_ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_SHELL_SEPARATOR_RE = re.compile(r"\|\||&&|[|;&]")
_PYTHON_RE = re.compile(r"python[\d.]*")

# Options whose value is the next token, so that token is no target.
_VALUE_OPTIONS = frozenset({
    "-k", "-m", "-p", "-c", "-o", "-W", "-n", "--deselect", "--ignore",
    "--ignore-glob", "--rootdir", "--confcutdir", "--basetemp",
    "--maxfail", "--tb", "--timeout", "--durations", "--junitxml",
    "--import-mode", "--log-level",
})
# Options that run a subset of the named targets.
_SELECT_OPTIONS = frozenset({"-k", "-m", "--deselect"})
_SELECT_FLAGS = frozenset({"--lf", "--last-failed", "--sw", "--stepwise"})

_DOC_STEMS = frozenset({
    "readme", "changelog", "news", "contributing", "authors", "license",
    "notice",
})
_DOC_EXTENSIONS = frozenset({".md", ".rst", ".txt", ".html", ".htm", ".typ"})
_FAIL_STATES = ("failed", "error", "gave_up")


def _norm(path: str) -> str:
    """One spelling for comparison: forward slashes, no leading ``./``."""
    p = str(path or "").replace("\\", "/").strip()
    while p.startswith("./"):
        p = p[2:]
    return p


def _is_test_file(path: str) -> bool:
    """``test_*.py`` or ``*_test.py``, as api_client._looks_like_test_file."""
    base = _norm(path).rsplit("/", 1)[-1]
    return ((base.startswith("test_") and base.endswith(".py"))
            or base.endswith("_test.py"))


def _is_product_file(path: str) -> bool:
    """A changed file whose fix owes a control run: not a test file and
    not documentation (README, .md, .rst, ...)."""
    base = _norm(path).rsplit("/", 1)[-1]
    if _is_test_file(base):
        return False
    stem, dot, ext = base.rpartition(".")
    if not dot:
        stem, ext = base, ""
    if ext and "." + ext.lower() in _DOC_EXTENSIONS:
        return False
    return stem.lower() not in _DOC_STEMS


def _file_part(target: str) -> str:
    """The test file a target names (before any ``::``)."""
    return _norm(target.split("::", 1)[0])


def _same_test_file(a: str, b: str) -> bool:
    """Two spellings of one test file (suffix-tolerant, as
    api_client._same_test_file)."""
    an, bn = _norm(a), _norm(b)
    if not an or not bn:
        return False
    return an == bn or an.endswith("/" + bn) or bn.endswith("/" + an)


def _split(text: str) -> List[str]:
    try:
        return shlex.split(text)
    except ValueError:
        return text.split()


def _runner_args(tokens: List[str]) -> Optional[List[str]]:
    """The arguments after the test runner in one shell segment, or None
    when the segment runs no test runner."""
    for i, tok in enumerate(tokens):
        base = tok.rsplit("/", 1)[-1]
        if base in ("pytest", "py.test", "gate") \
                or tok.endswith("gate-tools/gate"):
            return tokens[i + 1:]
        if (_PYTHON_RE.fullmatch(base) and tokens[i + 1:i + 2] == ["-m"]
                and tokens[i + 2:i + 3] in (["pytest"], ["unittest"])):
            return tokens[i + 3:]
    return None


def _collect(args: List[str], targets: List[str],
             selection: List[str]) -> None:
    """Targets (``.py`` paths, ``::`` node ids) and selection options of
    one runner's arguments."""
    pending = ""
    for tok in args:
        if pending:
            if pending in _SELECT_OPTIONS:
                selection.append(f"{pending} {tok}")
            pending = ""
            continue
        if tok.startswith("-"):
            name, eq, value = tok.partition("=")
            if name in _SELECT_FLAGS:
                selection.append(name)
            elif name in _VALUE_OPTIONS:
                if eq:
                    if name in _SELECT_OPTIONS:
                        selection.append(f"{name} {value}")
                else:
                    pending = name
            elif not tok.startswith("--") and tok[:2] in ("-k", "-m") \
                    and len(tok) > 2:
                selection.append(f"{tok[:2]} {tok[2:]}")
            continue
        if ".py" in tok or "::" in tok:
            targets.append(tok)


def _parse_test_command(command: str) -> Tuple[List[str], str]:
    """``(targets, selection)`` of a test command.

    Every shell segment (split at ``|``, ``&&``, ``||``, ``;``, ``&``)
    that runs a test runner contributes its targets; ``selection`` joins
    the selecting options, ``""`` when the run is unselected. A command
    with no runner (the live ledger stores a run_tests target as its
    command) is read as bare runner arguments.
    """
    targets: List[str] = []
    selection: List[str] = []
    found = False
    for segment in _SHELL_SEPARATOR_RE.split(command or ""):
        args = _runner_args(_split(segment))
        if args is not None:
            found = True
            _collect(args, targets, selection)
    if not found:
        _collect(_split(command or ""), targets, selection)
    return targets, " ".join(selection)


def _failing_node_targets(output: str) -> List[str]:
    """Node ids named on a run's ``FAILED`` summary lines, colours
    stripped. ``ERROR`` lines (setup, collection) are not failures."""
    text = _ANSI_ESCAPE_RE.sub("", output or "")
    return _FAILED_NODE_RE.findall(text)


def _selected_key(target: str, selection: str) -> str:
    """The target key of a selected run: its file (or node) with the
    selection, e.g. ``tests/test_x.py::(-k test_a)``."""
    sep = "" if "::" in target else "::"
    return f"{_norm(target)}{sep}({selection})"


def _index_runs(test_runs: Iterable,
                reds: Optional[Dict[str, List[float]]] = None
                ) -> Dict[str, Tuple]:
    """Per target key: ``(red, green)``.

    Each side is ``None`` (never seen), ``"*"`` (seen, no timestamp) or
    the earliest float timestamp it was seen at. ``reds``, when given,
    receives every float red timestamp per key.
    """
    red_ts: Dict[str, object] = {}
    green_ts: Dict[str, object] = {}
    file_green: Dict[str, float] = {}
    generic_green: Optional[float] = None
    red_from: Dict[str, set] = {}
    reds = {} if reds is None else reds

    def _track(store: dict, key: str, ts) -> None:
        cur = store.get(key)
        if cur is None:
            store[key] = ts if ts is not None else "*"
        elif ts is not None and (cur == "*" or ts < cur):
            store[key] = ts

    def _red(key: str, ts) -> None:
        _track(red_ts, key, ts)
        if ts is not None:
            reds.setdefault(key, []).append(ts)

    for run in test_runs or ():
        if not isinstance(run, dict):
            continue
        try:
            ts = float(run["ts"]) if run.get("ts") is not None else None
        except (TypeError, ValueError):
            ts = None
        target = str(run.get("target") or run.get("file") or "").strip()
        selection = str(run.get("selection") or "").strip()
        key = _selected_key(target, selection) if selection else _norm(target)
        if run.get("red"):
            nodes = [_norm(n) for n in
                     _failing_node_targets(run.get("output") or "")]
            # A red run that names its failing nodes is red on those
            # nodes; the file it ran is red only when no node is named.
            # A selection stays red as well: its own later green is the
            # one that covers those nodes.
            if key and (selection or not nodes):
                _red(key, ts)
            for node in nodes:
                _red(node, ts)
                if key:
                    red_from.setdefault(node, set()).add(key)
        if run.get("green"):
            if key:
                _track(green_ts, key, ts)
                if not selection and "::" not in target and ts is not None:
                    f = _norm(target)
                    file_green[f] = min(file_green.get(f, ts), ts)
            elif not selection and ts is not None:
                generic_green = (ts if generic_green is None
                                 else min(generic_green, ts))

    out: Dict[str, Tuple] = {}
    for key in set(red_ts) | set(green_ts):
        if not _file_part(key):
            continue          # a selected whole-suite run: only folds below
        rt, gt = red_ts.get(key), green_ts.get(key)
        if gt is None and isinstance(rt, float):
            covering = [file_green.get(_file_part(key)), generic_green]
            covering += [green_ts.get(k) for k in red_from.get(key, ())]
            later = [g for g in covering if isinstance(g, float) and g > rt]
            if later:
                gt = min(later)
        out[key] = (rt, gt)
    return out


def normalize_test_runs(test_runs: Iterable) -> list:
    """Adapt per-run ledger records to :func:`assess` input.

    Input: ``{"command", "status", "failed", "ts", "output", ...}`` as
    :func:`durable_test_runs` and the live ledger write them. Output: one
    ``{"target", "red", "green", "selection"[, "ts", "output"]}`` per
    target the command names; a run naming no target yields one record
    with ``target=""`` (a whole-suite run) when it is green or names
    failing nodes. Green: no failures and a status that is no failure.
    Records already in assess shape (``target``/``red``/``green``) pass
    through.
    """
    out = []
    if not isinstance(test_runs, (list, tuple, set)):
        return out
    for r in test_runs:
        if not isinstance(r, dict):
            continue
        if "target" in r or "red" in r or "green" in r:
            out.append(r)
            continue
        status = str(r.get("status", "") or "")
        green = int(r.get("failed", 0) or 0) == 0 \
            and status not in _FAIL_STATES
        output = r.get("output")
        targets, selection = _parse_test_command(
            str(r.get("command") or "").strip())
        if not targets:
            if green or _failing_node_targets(output or ""):
                targets = [""]
        for target in targets:
            item = {"target": target, "red": not green, "green": green,
                    "selection": selection}
            if r.get("ts") is not None:
                item["ts"] = r["ts"]
            if output:
                item["output"] = output
            out.append(item)
    return out


def _call_args(entry: dict) -> dict:
    """A trace entry's call arguments (stored as JSON text), or {}."""
    raw = entry.get("input")
    if isinstance(raw, dict):
        return raw
    try:
        parsed = json.loads(raw) if isinstance(raw, str) else None
    except ValueError:
        parsed = None
    return parsed if isinstance(parsed, dict) else {}


def _is_test_command(command: str) -> bool:
    """A test runner is on the command line."""
    return bool(command) and bool(_TEST_CMD_RE.search(command)
                                  or _GATE_CMD_RE.search(command))


def _counts(output: str) -> Tuple[int, int, int]:
    """``(passed, failed, errors)`` from pytest or unittest output."""
    def _n(rx) -> int:
        m = rx.search(output or "")
        return int(m.group(1)) if m else 0
    passed, failed, errors = (_n(_PYTEST_PASSED_RE), _n(_PYTEST_FAILED_RE),
                              _n(_PYTEST_ERROR_RE))
    m = _UNITTEST_FAIL_RE.search(output or "")
    if m:
        failed += int(m.group(1) or 0)
        errors += int(m.group(2) or 0)
    return passed, failed, errors


def _run_tests_run(entry: dict):
    """``(command, output, passed, failed, exit_code)`` of a run_tests
    call, or None when its result holds no verdict.

    The call is recorded with its own arguments (``target``,
    ``pytest_args``) and a JSON result (``status``, ``summary``,
    ``failures[].node_id``); a timeout, a stopped or a background run, a
    usage error and an errors-only run hold none.
    """
    try:
        data = json.loads(str(entry.get("output") or ""))
    except ValueError:
        return None
    if not isinstance(data, dict):
        return None
    status = str(data.get("status") or "")
    summary = data.get("summary") or {}
    try:
        passed = int(summary.get("passed", 0) or 0)
        failed = int(summary.get("failed", 0) or 0)
        errors = int(summary.get("errors", 0) or 0)
    except (TypeError, ValueError, AttributeError):
        return None
    if status == "ok" and failed == 0 and errors == 0 and passed > 0:
        exit_code = 0
    elif status == "failed" and not (failed == 0 and errors):
        exit_code = 1
    else:
        return None
    args = _call_args(entry)
    extra = args.get("pytest_args")
    words = [str(args.get("target") or "")]
    words += [str(a) for a in extra] if isinstance(extra, list) else []
    nodes = [str(f.get("node_id") or "") for f in data.get("failures") or []
             if isinstance(f, dict)]
    output = "\n".join(f"FAILED {n}" for n in nodes if "::" in n)
    return (shlex.join([w for w in words if w]), output, passed,
            failed + errors, exit_code)


def _shell_run(entry: dict):
    """``(command, output, passed, failed, exit_code)`` of a shell test
    run, or None when it is no test run or shows no executed test."""
    cmd = str(_call_args(entry).get("command") or "")
    raw = entry.get("input")
    if not cmd and isinstance(raw, str) and not raw.lstrip().startswith("{"):
        cmd = raw            # an older recorder stored the bare command
    if not _is_test_command(cmd):
        return None
    output = str(entry.get("output") or "")
    ok = entry.get("ok")
    exit_code = (0 if ok else 1) if isinstance(ok, bool) else None
    passed, failed, errors = _counts(output)
    if failed == 0 and (errors or _COLLECTION_ERROR_RE.search(output)):
        return None          # errors only: the assertions never ran
    if failed == 0 and exit_code == 1:
        failed = 1           # red by exit status; the summary was cut
    if failed == 0 and passed == 0:
        return None          # --version, --collect-only: no test ran
    return cmd, output, passed, failed, exit_code


def durable_test_runs(trace_entries, *, now=None) -> list:
    """The session's test runs, across turns, from its tool trace.

    Input: tool_trace entries ``{tool, input, output, ok, ts}``. Output:
    one record per test run, ``{"tool", "command", "output",
    "exit_code", "status", "passed", "failed", "ts"}`` (the live ledger's
    shape plus the output), for :func:`normalize_test_runs`. ``now``
    replaces every timestamp (deterministic tests).
    """
    out = []
    if not isinstance(trace_entries, (list, tuple, set)):
        return out
    for e in trace_entries:
        if not isinstance(e, dict):
            continue
        tool = str(e.get("tool", "") or "")
        if tool.rsplit("__", 1)[-1] == "run_tests":
            run = _run_tests_run(e)
        else:
            run = _shell_run(e)
        if run is None:
            continue
        cmd, output, passed, failed, exit_code = run
        out.append({
            "tool": tool or "bash",
            "command": cmd[:300],
            "output": output[:2000],
            "exit_code": exit_code,
            "status": "failed" if failed else "ok",
            "passed": passed,
            "failed": failed,
            "ts": float(e.get("ts") or 0) if now is None else float(now),
        })
    return out


def _changes(changed_files: Iterable) -> List[Tuple[str, Optional[float]]]:
    """``(path, ts)`` per changed file; ``ts`` None when not stamped."""
    out = []
    for f in changed_files or ():
        if isinstance(f, dict):
            path, ts = f.get("path"), f.get("ts")
        else:
            path, ts = f, None
        path = _norm(path)
        if not path:
            continue
        try:
            ts = float(ts) if ts else None
        except (TypeError, ValueError):
            ts = None
        out.append((path, ts))
    return out


def assess(changed_files: Iterable, test_runs: Iterable) -> Optional[str]:
    """The report line for the session's control-run evidence.

    ``changed_files``: changed paths, or change-journal records
    ``{"path", "ts"}``. ``test_runs``: :func:`normalize_test_runs` output
    (``{"target" | "file", "red", "green", "selection", "ts", "output"}``).

    Precedence:
    1. a changed test file no run covers: changed but never run;
    2. in a changed test file: a target green but never red;
    3. ... green before red;
    4. ... red but never green afterwards;
    5. ... controlled, but the test file changed after its last red;
    6. ... a controlled target: ``"shown"``;
    7. product-only change: any controlled node is ``"shown"``; else the
       first missing half; else no run recorded.
    """
    changes = _changes(changed_files)
    changed = [p for p, _ in changes]
    if not any(_is_product_file(f) for f in changed):
        return None
    test_changed = [f for f in changed if _is_test_file(f)]
    reds: Dict[str, List[float]] = {}
    runs = _index_runs(test_runs, reds)
    run_files = {_file_part(t) for t in runs}

    for t in sorted(set(test_changed)):
        if not any(_same_test_file(t, f) for f in run_files):
            return f"no control run: test {t} was changed but never run"

    def _says(t: str, what: str) -> str:
        return f"no control run: test {t} {what}"

    def _floats(entry) -> bool:
        return all(isinstance(v, float) for v in entry)

    def _controlled(entry) -> bool:
        rt, gt = entry
        if rt is None or gt is None:
            return False
        return rt < gt if _floats(entry) else True

    red_files = {_file_part(t) for t in runs if "::" in t and runs[t][0]}

    def _never_red(t: str) -> bool:
        """Green, never red; a whole-file green whose file has a red node
        is that node's green, not a target of its own."""
        rt, gt = runs[t]
        if "::" not in t and t in red_files:
            return False
        return rt is None and gt is not None

    def _changed_after_red(t: str) -> bool:
        gt = runs[t][1]
        before = [r for r in reds.get(t, ()) if isinstance(gt, float) and r < gt]
        if not before:
            return False
        last_red = max(before)
        return any(ts is not None and last_red < ts < gt
                   and _is_test_file(p) and _same_test_file(p, _file_part(t))
                   for p, ts in changes)

    if test_changed:
        scope = [t for t in sorted(runs)
                 if any(_same_test_file(c, _file_part(t))
                        for c in test_changed)]
        for t in scope:
            if _never_red(t):
                return _says(t, "was never seen failing")
        for t in scope:
            if _floats(runs[t]) and runs[t][1] <= runs[t][0]:
                return _says(t, "was seen green before it was seen red")
        for t in scope:
            if runs[t][0] is not None and not _controlled(runs[t]):
                return _says(
                    t, "was seen failing but never passing after the change")
        for t in scope:
            if _changed_after_red(t):
                return _says(t, "was changed after it was seen failing")
        if any(_controlled(runs[t]) for t in scope):
            return "shown"
        return "no control run: no test run was recorded for the change"

    controlled = [t for t in sorted(runs) if _controlled(runs[t])]
    if any("::" in t for t in controlled):
        return "shown"
    if controlled:
        return ("no control run: the run on "
                f"{controlled[0]} does not name a node that was "
                "seen red and green")
    for t in sorted(runs):
        if _never_red(t):
            return _says(t, "was never seen failing")
    for t in sorted(runs):
        if runs[t][0] is not None:
            return _says(
                t, "was seen failing but never passing after the change")
    return "no control run: no test run was recorded for the change"
