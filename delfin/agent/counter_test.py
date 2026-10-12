"""Does a fix name its counter-test?

Four waves of reviews shared one recurring finding: a test that passes
with and without the fix, or a guard tested only in the direction "it
fires". A review step that asks "did the tests run green?" is answered
by a green run — it does not ask whether the same test was EVER seen
red before the change. A control run is the only evidence a fix actually
fixed something: the test, failing before the change (or against a copy
without it), passing after.

This module decides, from what a turn recorded, whether such a control
run is on record, and when it is not, gives the report the one sentence
it must carry. It is pure: both arguments are derived values, it does no
I/O, and it never blocks a turn — it only decides one report line
(``None`` for "add nothing").

WHICH LEDGER CARRIES THE RED-BEFORE (reviewers s27/s28, a correctness
crux resolved here): the verify step's in-memory test ledger
(``_observe_test_evidence``, api_client.py:11126) is PER-TURN — the
client replaces it at every turn start (api_client.py:6026-6030), and
``red_files`` is a live set cleared on a green run
(``_clear_green_targets``, api_client.py:10811). So the "seen red before
the change" observation may live in a PREVIOUS turn and is NOT in the
current turn's ``tests`` argument. The durable, cross-turn source with
per-node identity is ``delfin.agent.tool_trace`` (record :73, read :120):
every ``run_tests`` call, with its JSON result including the failing
NODE ids (``test_x.py::test_y``) — the same source report_flakes.py:127-139
reads. The phase-3 wiring therefore builds ``test_runs`` from tool_trace
(or a persisted failing-node record), not from the ephemeral per-turn
ledger, so a red from a previous turn is inside what ``assess`` sees.
"If the red half would not be visible, the sentence is wrong" — the
caller (task_evidence.check_completion_claim, task_evidence.py:327) is
told this in .gate/HANDOFF.md.

A control is per TEST TARGET, not per file: the SAME node
(``tests/test_x.py::test_y``) seen red and green. A red run of one node
and a green run of a different node in the same FILE do not combine into
a control — that is exactly the "test passes with and without the fix"
failure this package exists to catch (reviewer s28, probe #1).

KNOWN LIMITATION (QS ruling, condition 1): a NEW test NODE added to an
existing file is INVISIBLE to ``assess`` — ``changed_files`` carries
paths, ``test_runs`` carries recorded runs, and neither records a node
that was never run. If that brand-new node never runs, a red+green on an
unrelated, already-controlled node in the SAME file can still read as
``"shown"``. The inputs cannot represent the shape, so no rule-4 wording
can catch it; the client wiring must carry node-level change info (the
added test defs the git diff gives, and one run entry per failing/changed
node — report_flakes.node ids) to make the never-run-new-node case
representable. This module documents the gap and does not guess.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple


def _norm(path: str) -> str:
    """One canonical spelling for path comparison."""
    return str(path or "").replace("\\", "/").lstrip("./")


def _is_test_file(path: str) -> bool:
    """A file is a test file when its basename looks like pytest's.
    Mirrors api_client._looks_like_test_file (api_client.py:10805):
    ``test_*.py`` or ``*_test.py``."""
    base = str(path or "").replace("\\", "/").rsplit("/", 1)[-1]
    return ((base.startswith("test_") and base.endswith(".py"))
            or base.endswith("_test.py"))


# Documentation / plain-text artifacts are not "the fix": no control run
# is owed for a docs-only change, and reporting one for it would claim a
# test verified a README edit. Only source or behavior-changing files are
# product files here. Keep the set small, by content kind, not an allow
# list of every code extension.
_DOC_STEMS = frozenset({
    "readme", "changelog", "news", "contributing", "authors", "license",
    "notice",
})
_DOC_EXTENSIONS = frozenset({".md", ".rst", ".txt", ".html", ".htm", ".typ"})


def _is_product_file(path: str) -> bool:
    """A changed file whose fix owes a control run.

    A test file is never product. A documentation/plain-text artifact
    (README, .md, .rst, …) is treated as non-product: it is not the
    code a test would exercise, so it needs no counter-test. Everything
    else — source, config that changes behaviour — is a product file.
    """
    base = str(path or "").replace("\\", "/").rsplit("/", 1)[-1]
    if _is_test_file(base):
        return False
    ext = "." + base.rsplit(".", 1)[1].lower() if "." in base else ""
    if ext in _DOC_EXTENSIONS:
        return False
    stem = base.rsplit(".", 1)[0].lower() if "." in base else base.lower()
    return stem not in _DOC_STEMS


def _target_of(run: dict) -> str:
    """The test TARGET a recorded run exercised.

    Precision over per-file evidence: the durable ledger reports failing
    NODE ids ("tests/test_x.py::test_y") and the run's pytest target
    ("tests/test_x.py" or a node id), so the identity of a control is
    the node, not just the file. Accepts either the precise ``target``
    key or a plain ``file`` key (a whole-file run then has target=file).
    """
    tgt = str((run or {}).get("target") or "").strip()
    if tgt:
        return tgt
    f = str((run or {}).get("file") or "").strip()
    return _norm(f)


def _file_part(target: str) -> str:
    """The test FILE a target names (before any ``::node``)."""
    return _norm(target.split("::", 1)[0])


def _same_test_file(a: str, b: str) -> bool:
    """Whether two test-FILE spellings name the same file (suffix-tolerant,
    as api_client._same_test_file :10795)."""
    an, bn = _norm(a), _norm(b)
    if not an or not bn:
        return False
    return an == bn or an.endswith("/" + bn) or bn.endswith("/" + an)


def _index_runs(test_runs: Iterable) -> Dict[str, Tuple[bool, bool]]:
    """Aggregate runs per TEST TARGET -> (seen_red, seen_green)."""
    out: Dict[str, Tuple[bool, bool]] = {}
    for run in test_runs or ():
        if not isinstance(run, dict):
            continue
        tgt = _target_of(run)
        if not tgt:
            continue
        red = bool(run.get("red"))
        green = bool(run.get("green"))
        prev_red, prev_green = out.get(tgt, (False, False))
        out[tgt] = (prev_red or red, prev_green or green)
    return out


def normalize_test_runs(test_runs: Iterable) -> list:
    """Adapt the verify ledger's run records to :func:`assess` input shape.

    The verify step's ledger (``_observe_test_evidence``, api_client.py:11126)
    stores each run as ``{"tool", "command", "exit_code", "status", "passed",
    "failed", "ts"}`` — a per-RUN shape. ``assess`` expects per-run ``{"target",
    "red", "green"}``. This pure adapter maps one onto the other: the run's
    ``command`` (the pytest target / command) becomes ``target``; red and green
    follow the same rule ``check_completion_claim`` already uses for "green"
    (task_evidence.py:506-507): a run is green when neither ``failed`` nor
    ``status`` reports failure, red otherwise. Runs that say nothing about
    either (no ``command``/``target``, no status) are dropped — they carry no
    target to judge a control by.
    """
    _FAIL_STATES = ("failed", "error", "gave_up")
    out = []
    for r in (test_runs or ()) if isinstance(test_runs, (list, tuple, set)) else []:
        if not isinstance(r, dict):
            continue
        if "target" in r or "red" in r or "green" in r:
            out.append(r)  # already in assess shape
            continue
        target = str(r.get("command") or r.get("target") or "").strip()
        if not target:
            continue
        status = str(r.get("status", "") or "")
        failed = int(r.get("failed", 0) or 0)
        green = failed == 0 and status not in _FAIL_STATES
        out.append({"target": target, "red": not green, "green": green})
    return out


# PS # -- the durable cross-turn reconstruction --------------------------
# The parsing regexes below mirror delfin/agent/api_client.py:10732-10782
# (`_TEST_BASH_CMD_RE`, `_parse_test_counts`, `_failing_test_files`,
# `_clear_green_targets`) so the durable reconstruction in
# `durable_test_runs` classifies a run exactly the way the live per-turn
# ledger (`_observe_test_evidence`, api_client.py:11126) does. This module
# stays dependency-free; the regexes are copied so a pure, gated test can
# run without importing the api client.

import re as _re

_TEST_CMD_RE = _re.compile(
    r"(?:^|[;&|(])\s*(?:python[\d.]*\s+-m\s+(?:pytest|unittest)\b"
    r"|pytest\b|py\.test\b)"
)
_PYTEST_PASSED_RE = _re.compile(r"(\d+)\s+passed\b")
_PYTEST_FAILED_RE = _re.compile(r"(\d+)\s+failed\b")
_PYTEST_ERROR_RE = _re.compile(r"(\d+)\s+errors?\b")
_UNITTEST_FAIL_RE = _re.compile(
    r"FAILED\s*\((?:failures=(\d+))?(?:,\s*)?(?:errors=(\d+))?\)")
_PYTEST_FAILED_LINE_RE = _re.compile(
    r"^(?:FAILED|ERROR)\s+([\w./\\-]+\.py)", _re.MULTILINE)
_ANSI_ESCAPE_RE = _re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_FAIL_STATES = ("failed", "error", "gave_up")


def _parse_test_counts(output: str):
    """(passed, failed+errors) from pytest/unittest console output.
    (0, 0) when nothing matched — combine with the exit code/admitted
    failure, mirroring api_client._parse_test_counts :10750."""
    text = output or ""
    passed = failed = errors = 0
    m = _PYTEST_PASSED_RE.search(text)
    if m:
        passed = int(m.group(1))
    m = _PYTEST_FAILED_RE.search(text)
    if m:
        failed = int(m.group(1))
    m = _PYTEST_ERROR_RE.search(text)
    if m:
        errors = int(m.group(1))
    m = _UNITTEST_FAIL_RE.search(text)
    if m:
        failed += int(m.group(1) or 0)
        errors += int(m.group(2) or 0)
    return passed, failed + errors


def _failing_test_files(output: str):
    """Test files named in 'FAILED/ERROR tests/x.py' summary lines, colors
    stripped — mirror api_client._failing_test_files :10774."""
    clean = _ANSI_ESCAPE_RE.sub("", output or "")
    return {m.group(1) for m in _PYTEST_FAILED_LINE_RE.finditer(clean)}


def _trace_entry_command(entry: dict) -> str:
    """The shell command of a durable tool_trace entry, or ''.

    A trace entry carries its call input as JSON (``tool_trace.call_args``
    :145); a bash-like call's target test command is the ``command``
    argument when present, else the entry's own ``tool`` name cannot say
    which tests it ran. Mirrors tool_trace.command_of :171 but also reads
    the raw input string so the pure reconstruction needs no import."""
    raw = entry.get("input")
    if isinstance(raw, dict):
        return str(raw.get("command") or "")
    if isinstance(raw, str):
        try:
            parsed = __import__("json").loads(raw)
        except (ValueError, TypeError):
            parsed = None
        if isinstance(parsed, dict):
            return str(parsed.get("command") or "")
        return raw[:400] if _TEST_CMD_RE.search(raw) else ""
    return ""


def _run_status_from(failed: int, exit_code) -> str:
    """The same status word ``_observe_test_evidence`` would assign."""
    if failed:
        return "failed"
    if exit_code == 0:
        return "ok"
    if exit_code is None:
        return "unknown"
    return "failed" if exit_code != 0 else "ok"


def durable_test_runs(trace_entries, *, now=None) -> list:
    """Reconstruct the per-turn test-run ledger from the DURABLE tool_trace.

    ``assess`` needs a red from an EARLIER turn to decide a control, but
    the live verify ledger (``_observe_test_evidence``) is per-turn and is
    replaced at every turn start — so the "seen red before the change"
    observation is not in the current turn's ``tests`` argument. This pure
    function rebuilds the run records across turns from the durable
    ``tool_trace`` entries (``{tool, input, output, ok, ts}``): every call
    whose command targets pytest becomes one record in the same shape the
    per-turn ledger emits — ``{"tool", "command", "exit_code", "status",
    "passed", "failed", "ts"}`` — which ``normalize_test_runs`` then maps
    to ``{"target", "red", "green"}`` for ``assess``.

    A green run is one that demonstrably executed tests and failed none: a
    bare exit-0 from ``--version`` is not a run (mirrors
    ``_observe_test_evidence`` :11189-11192). An entry that carries no
    test command or no readable evidence is dropped.

    ``now`` pins the timestamp for determinism in tests (defaults to each
    entry's recorded ``ts``).

    Returns a list of run records suitable for ``normalize_test_runs``.
    """
    out = []
    for e in (trace_entries or ()) if isinstance(trace_entries, (list, tuple, set)) else []:
        if not isinstance(e, dict):
            continue
        cmd = _trace_entry_command(e)
        if not _TEST_CMD_RE.search(cmd):
            continue
        tool = str(e.get("tool", "") or "")
        output = str(e.get("output", "") or "")
        exit_code = e.get("ok")
        if isinstance(exit_code, bool):
            # ``ok`` (True/False) is the durable record's success flag; the
            # per-turn shape wants the integer exit code we infer from it.
            exit_code = 0 if exit_code else 1
        try:
            exit_code = int(exit_code)
        except (TypeError, ValueError):
            exit_code = None
        passed, failed = _parse_test_counts(output)
        # A non-zero exit beats a missing summary line (api_client :11179):
        # a failed run whose output lacks the "N failed" line is still red.
        if failed == 0 and isinstance(exit_code, int) and exit_code != 0:
            failed = 1
        # A run that demonstrated no tests (no "N passed", no "N failed")
        # is not evidence of anything — a bare `--version`, `--help` or
        # `--collect-only` — and must not masquerade as a green run
        # (api_client :11189-11192: only a run that demonstrably executed
        # tests clears red state).
        if failed == 0 and passed == 0:
            continue
        status = _run_status_from(failed, exit_code)
        ts = float(e.get("ts") or 0) if now is None else float(now)
        out.append({
            "tool": (tool or "bash"),
            "command": cmd[:300],
            "exit_code": exit_code,
            "status": status,
            "passed": passed,
            "failed": failed,
            "ts": ts,
        })
    return out


def _in_fix_scope(target: str, test_changed: List[str]) -> bool:
    """Whether a run ``target`` belongs to the fix's own counter-test.

    The fix's counter-test lives in the test file(s) the turn changed;
    a target is in scope when its test FILE is one of them. This is the
    guard that keeps an unrelated target's red+green from being read as
    a control for the fix (QS findings A/B/C): only evidence on the
    fix's own changed test records a control.
    """
    ft = _file_part(target)
    return any(_same_test_file(t, ft) for t in test_changed)


def assess(changed_files: Iterable, test_runs: Iterable) -> Optional[str]:
    """The report line for the turn's control-run evidence.

    ``changed_files``: paths the turn changed (product AND test files).
    ``test_runs``: recorded runs from the DURABLE ledger — every run
    ``{"target": <pytest target or node id>, "red": bool,
    "green": bool}`` (a whole-file run may use ``file``).

    A control is **per node, on the fix's own counter-test**: the SAME
    test node (``tests/test_x.py::test_y``) inside a changed test file
    seen red and green. A red+green on an unrelated node or an unrelated
    file is evidence for that node, not for the fix, and never proves a
    control for it (reviewers s27/s28 probes). A whole-file target
    (no ``::node``) is never a completed control: it cannot say which
    node carried the red (QS finding B).

    Returns:
    * ``None`` — the turn changed no product file; nothing is added.
    * ``"shown"`` — a control run is on record on the fix's own
      counter-test: a node inside a changed test file seen both red and
      green. When NO test file changed (a product-only fix), any
      controlled node stands in for the (pre-existing) counter-test.
    * one sentence naming what is missing otherwise — and a missing
      control on the fix's own test is not masked by a completed control
      on an unrelated file or an unrelated node (exception: a NEW node
      added to an existing file that never runs is invisible to both
      inputs, so a same-file controlled old node can still read as
      "shown" for it — see the module's KNOWN LIMITATION).

    Precedence (documented, not guessed):
    1. a changed test file that no run target covers — the fix's own
       counter-test sitting unrun;
    2. the fix's own scope: a node (or file) seen green but never red —
       the classic "passes with and without the fix";
    3. the fix's own scope: a node red but never green after the change;
    4. the fix's own scope: a node seen both red and green -> "shown";
    5. a product-only fix (no test file changed): a controlled node
       anywhere -> "shown"; an uncontrolled one names what it lacks;
    6. otherwise — a product change with no recorded run at all.
    """
    if not changed_files:
        return None
    changed = [_norm(f) for f in changed_files if _norm(f)]
    product_changed = [f for f in changed if _is_product_file(f)]
    if not product_changed:
        return None

    test_changed = [f for f in changed if _is_test_file(f)]
    runs_by_target = _index_runs(test_runs)
    run_files = {_file_part(t) for t in runs_by_target}

    # 1. A changed test file whose file is never covered by any run
    #    target. This is the fix's own counter-test unrun.
    for t in sorted(test_changed):
        if not any(_same_test_file(t, f) for f in run_files):
            return f"no control run: test {t} was changed but never run"

    def _node(t: str) -> bool:
        """A target that names a specific test node (not a whole file)."""
        return "::" in t

    def _target_sentence(t: str, verdict: str) -> str:
        return f"no control run: test {t} {verdict}"

    if test_changed:
        # 2/3/4. The fix's own counter-test decides. Evaluate its scope
        # in strict precedence, so a green-only node in the fix's own
        # file is never masked by a red+green on a different node.
        scope = [t for t in sorted(runs_by_target)
                 if _in_fix_scope(t, test_changed)]
        # 2. green but never red — the passing-without-the-fix marker.
        never_red = [t for t in scope if not runs_by_target[t][0]
                     and runs_by_target[t][1]]
        if never_red:
            return _target_sentence(never_red[0], "was never seen failing")
        # 3. red but never green after the change.
        red_not_green = [t for t in scope
                         if runs_by_target[t][0] and not runs_by_target[t][1]]
        if red_not_green:
            return _target_sentence(
                red_not_green[0],
                "was seen failing but never passing after the change")
        # 4. a completed control on a NODE inside the fix's own file.
        controlled = [t for t in scope
                      if runs_by_target[t][0] and runs_by_target[t][1]]
        if any(_node(t) for t in controlled):
            return "shown"
        # The fix's scope ran (whole file only) but no node was pinned.
        if controlled and not any(_node(t) for t in controlled):
            return ("no control run: the run on "
                    f"{controlled[0]} does not name a node that was "
                    "seen red and green")
        # The scope has no green and no red — no usable run for it.
        return "no control run: no test run was recorded for the change"

    # 5. A product-only fix (no test file changed): any controlled node
    #    stands in for the counter-test. A whole-file target that claims
    #    red AND green is ambiguous (which node failed?) and never reads
    #    as a completed control (QS finding B).
    controlled = [t for t in runs_by_target
                  if runs_by_target[t][0] and runs_by_target[t][1]]
    if any(_node(t) for t in controlled):
        return "shown"
    if controlled:
        return ("no control run: the run on "
                f"{controlled[0]} does not name a node that was "
                "seen red and green")
    never_red = [t for t in sorted(runs_by_target)
                 if not runs_by_target[t][0] and runs_by_target[t][1]]
    if never_red:
        return _target_sentence(never_red[0], "was never seen failing")
    red_not_green = [t for t in sorted(runs_by_target)
                     if runs_by_target[t][0] and not runs_by_target[t][1]]
    if red_not_green:
        return _target_sentence(
            red_not_green[0],
            "was seen failing but never passing after the change")

    # 6. A product change with no recorded run at all.
    return "no control run: no test run was recorded for the change"
