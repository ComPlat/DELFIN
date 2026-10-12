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


def _parse_test_command(command: str) -> list:
    """The concrete pytest targets a test command runs.

    Accepts the command forms DELFIN actually uses -- ``pytest\b``,
    ``py.test\b``, ``python -m pytest\b``/``python3 -m pytest\b``, and the
    DELFIN gate wrapper (``.../gate-tools/gate TESTS...`` / ``gate\b``) -- the
    F2 form every session uses -- and returns the target arguments: each
    token that names a ``.py`` file or a ``file::node`` id (several per
    command), with option tokens and their values skipped.

    Returns [] when the command is not a recognized test invocation.
    """
    if not command:
        return []
    tokens = command.split()
    if not tokens:
        return []
    # F2: recognize the test runner (incl. the gate wrapper) up front.
    def _is_runner(tok: str) -> bool:
        t = tok.strip()
        return (bool(_TEST_CMD_RE.search(t))
                or _GATE_CMD_RE.search(t)
                or t == "gate"
                or t.endswith("/gate")
                or t.endswith("gate-tools/gate"))
    hit = any(_is_runner(t) for t in tokens)
    if not hit:
        # A bare test-target command (phase-2 durable form):
        #   "tests/test_x.py::test_y"  or  "tests/test_x.py -q"
        # has no runner token but IS a test invocation. This parser is only
        # fed entries already classified as test runs, so a command that is
        # itself a .py/:: target counts as a test command.
        if not any(".py" in t or "::" in t for t in tokens):
            return []
    out: list = []
    skip_next = False
    for tok in tokens:
        if skip_next:
            skip_next = False
            continue
        if _is_runner(tok):
            continue
        if tok.startswith("-"):
            # an option; its inline value is part of it (--timeout=300)
            if not tok.startswith("--") or "=" not in tok:
                # a bare option that consumes the next token (e.g. -k expr)
                skip_next = True
            continue
        if ".py" in tok or "::" in tok:
            out.append(tok)
    return out


def _failing_node_targets(output: str) -> list:
    """Node ids named in a red run's own FAILED/ERROR lines.

    A whole-file red run reports which NODES failed in its console summary
    (``FAILED tests/test_x.py::test_a``). Those concrete nodes are the
    control's subject (F3): a later green run on the same file is a
    completed control for each. Colors are stripped first.
    """
    text = _ANSI_ESCAPE_RE.sub("", output or "")
    return _PYTEST_FAILED_NODE_RE.findall(text)


def _earlier(a, b):
    """The earlier KNOWN timestamp of two side-recordings.

    Each side value is ``None`` (never), ``"*"`` (seen, unknown ts) or a
    float. ``None``/``"*"`` are never "earlier" than a real float; between
    two unknown values neither is earlier (returns the first non-None or
    None)."""
    if isinstance(a, float) and isinstance(b, float):
        return a if a < b else b
    if isinstance(a, float):
        return a
    if isinstance(b, float):
        return b
    return a if a is not None else b


def _index_runs(test_runs: Iterable) -> Dict[str, Tuple]:
    """Aggregate runs per TEST TARGET -> (red_ts, green_ts).

    Order carries the control (F4): a target is controlled only when it was
    seen RED strictly BEFORE it was seen GREEN. Each side is stored as ONE
    value:
      * ``None``        -- that side was NEVER seen,
      * ``"*"``         -- seen, but at an UNKNOWN timestamp (a hand-made
        assess run with no ``ts``); treated permissively by assess,
      * a float         -- seen at that timestamp (real trace run).
    A red run's FAILED node lines pin concrete node targets at the run's ts
    (F3); a green run on a file covers that FILE's nodes, so a node red at
    ts_r is controlled by a green on the same file at a later ts.
    """
    red_ts: Dict[str, object] = {}
    green_ts: Dict[str, object] = {}
    file_green: Dict[str, object] = {}
    # A green run naming no test file ("pytest -q", "gate" bare) covers
    # EVERY suite, so it provides green coverage for any file a red run's
    # FAILED lines named (F5e).
    generic_green: object = None

    def _track(store: dict, tgt: str, ts) -> None:
        """Record a side for ``tgt`` seen at ``ts`` (float or None-unknown)."""
        cur = store.get(tgt)
        if cur is None:                    # never seen -> now seen
            store[tgt] = ts if ts is not None else "*"
        elif ts is None:                   # an unknown-ts sighting: keep a
            pass                           # known-ts one if we have it
        elif cur == "*" or ts < cur:
            store[tgt] = ts

    for run in test_runs or ():
        if not isinstance(run, dict):
            continue
        try:
            ts = float(run["ts"]) if run.get("ts") is not None else None
        except (TypeError, ValueError):
            ts = None
        red = bool(run.get("red"))
        green = bool(run.get("green"))
        tgt = _target_of(run)
        if red:
            if tgt:
                _track(red_ts, tgt, ts)
            # FAILED lines pin concrete node targets even when the command
            # carried no target path (F5e: bare "pytest -q").
            for n in _failing_node_targets(run.get("output") or ""):
                nt = _norm(n)
                if nt:
                    _track(red_ts, nt, ts)
        if green:
            if tgt:
                f = _file_part(tgt)
                if f:
                    _track(file_green, f, ts)
                _track(green_ts, tgt, ts)
            else:
                generic_green = _earlier(generic_green, ts)
    out: Dict[str, tuple] = {}
    for t in set(red_ts) | set(green_ts):
        rt = red_ts.get(t)
        gt = green_ts.get(t)
        if gt is None and rt not in (None, "*") and "." in t:
            # F3 fold: a later green run on the same FILE (or a whole suite,
            # generic_green) controls a target red at an earlier ts — for a
            # node target when its file passes later, AND for a whole-file
            # red run on a changed test file (s27 / operator F3) when that
            # file passes later. Only on KNOWN timestamps both sides.
            f = _file_part(t)
            fg = _earlier(file_green.get(f), generic_green)
            if isinstance(rt, float) and isinstance(fg, float) and fg > rt:
                gt = fg   # later whole-file green controls the target
        out[t] = (rt, gt)
    return out


def normalize_test_runs(test_runs: Iterable) -> list:
    """Adapt the verify ledger's run records to :func:`assess` input shape.

    The verify step's ledger (``_observe_test_evidence``, api_client.py:11126)
    stores each run as ``{"tool", "command", "exit_code", "status", "passed",
    "failed", "ts"}`` — a per-RUN shape, and the DURABLE reconstruction
    (``durable_test_runs``) keeps the same per-run shape with ``ts``. ``assess``
    expects per-run ``{"target", "red", "green"}`` with an optional ``ts`` and
    ``output``. This pure adapter maps one onto the other (F1): the run's
    ``command`` is PARSED into its concrete pytest targets (paths and node
    ids, several per command; gate-tool runs included — F2) via
    :func:`_parse_test_command``; one adapted record is emitted per target.
    red and green follow the same rule ``check_completion_claim`` already uses
    for "green" (task_evidence.py:506-507): a run is green when neither
    ``failed`` nor ``status`` reports failure, red otherwise. ``ts`` and the
    run's ``output`` (for F3's FAILED-node lines) are carried through. Runs
    that parse no test target are dropped — they carry no target to judge a
    control by. Dicts already in assess shape (``target``/``red``/``green``)
    pass through unchanged.
    """
    _FAIL_STATES = ("failed", "error", "gave_up")
    out = []
    for r in (test_runs or ()) if isinstance(test_runs, (list, tuple, set)) else []:
        if not isinstance(r, dict):
            continue
        if "target" in r or "red" in r or "green" in r:
            out.append(r)  # already in assess shape
            continue
        command = str(r.get("command") or r.get("target") or "").strip()
        ts = r.get("ts")
        output = r.get("output")
        status = str(r.get("status", "") or "")
        failed = int(r.get("failed", 0) or 0)
        green = failed == 0 and status not in _FAIL_STATES
        red = not green
        targets = _parse_test_command(command)
        if not targets:
            # A run whose command carried no test target path: a bare
            # "pytest -q" / "gate" run. It still counts when it holds
            # real evidence — a RED one carries FAILED-line nodes (F5e),
            # a GREEN one is whole-suite coverage for those nodes.
            if red and _failing_node_targets(output or ""):
                item = {"target": "", "red": True, "green": False}
                if ts is not None:
                    item["ts"] = ts
                if output:
                    item["output"] = output
                out.append(item)
            elif green:
                out.append({"target": "", "red": False, "green": True,
                            **({"ts": ts} if ts is not None else {}),
                            **({"output": output} if output else {})})
            continue
        for target in targets:
            item = {"target": target, "red": red, "green": green}
            if ts is not None:
                item["ts"] = ts
            if output:
                item["output"] = output
            out.append(item)
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
# A RED run's console summary names the concrete failing NODE ids
# (``FAILED tests/test_x.py::test_a``) -- F3 uses these to pin a whole-file
# red run's failing nodes so a later green on the same file is a control.
_PYTEST_FAILED_NODE_RE = _re.compile(
    r"^(?:FAILED|ERROR)\s+([\w./\\-]+\.py::[^\s:]+)", _re.MULTILINE)
# The DELFIN gate wrapper that every session runs tests through (F2):
#   .delfin/gate-tools/gate tests/test_x.py::test_a -q
_GATE_CMD_RE = _re.compile(r"(?:^|[;&|(])\s*(?:[\w./\\-]*gate-tools/gate|gate)\b")
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


def _is_test_command(command: str) -> bool:
    """Whether a durable shell command is a test invocation (F2).

    Accepts pytest/py.test/``python -m pytest|unittest`` AND the DELFIN gate
    wrapper (``.../gate-tools/gate TESTS...``) that every session actually
    runs tests through — plus a bare test-target command.
    """
    if not command:
        return False
    if bool(_TEST_CMD_RE.search(command)) or bool(_GATE_CMD_RE.search(command)):
        return True
    tokens = command.split()
    return any(".py" in t or "::" in t for t in tokens)


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
        return raw[:400] if _is_test_command(raw) else ""
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
        if not _is_test_command(cmd):
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
            "output": output[:2000],
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

    def _out_of_order(entry) -> bool:
        """Whether a target saw green strictly before red (F4)."""
        rt, gt = entry
        return (rt is not None and gt is not None
                and isinstance(rt, float) and isinstance(gt, float)
                and gt <= rt)

    def _controlled(entry) -> bool:
        """A completed control: red seen, green seen, red strictly BEFORE
        green. When timestamps are unknown (hand-made assess input with no
        ts), both-present counts as controlled (permissive, backwards
        compatible); when both known, the order decides.
        """
        rt, gt = entry
        if rt is None or gt is None:
            return False
        if not isinstance(rt, float) or not isinstance(gt, float):
            return True  # order unknown — both present
        return rt < gt

    if test_changed:
        # 2/3/4. The fix's own counter-test decides. Evaluate its scope
        # in strict precedence, so a green-only node in the fix's own
        # file is never masked by a red+green on a different node.
        scope = [t for t in sorted(runs_by_target)
                 if _in_fix_scope(t, test_changed)]
        # 2. green but never red — the passing-without-the-fix marker.
        never_red = [t for t in scope
                     if runs_by_target[t][0] is None
                     and runs_by_target[t][1] is not None]
        if never_red:
            return _target_sentence(never_red[0], "was never seen failing")
        # F4. green seen BEFORE red — not a control (the change broke the
        # test) and must not read as "shown".
        broke = [t for t in scope if _out_of_order(runs_by_target[t])]
        if broke:
            return _target_sentence(
                broke[0],
                "was seen green before it was seen red")
        # 3. red but never green after the change.
        red_not_green = [t for t in scope
                         if runs_by_target[t][0] is not None
                         and not _controlled(runs_by_target[t])]
        if red_not_green:
            return _target_sentence(
                red_not_green[0],
                "was seen failing but never passing after the change")
        # 4. a completed control on the fix's own counter-test: a NODE
        #    inside the fix's own changed test file, or (operator F3 /
        #    s27) a whole-file red-then-green on that changed test file,
        #    reads as "shown".
        controlled = [t for t in scope if _controlled(runs_by_target[t])]
        if any(_node(t) for t in controlled):
            return "shown"
        # Whole-file control on the fix's own changed test file: the file
        # was seen red then passed — a completed control for its scope.
        whole_controls = [t for t in controlled
                          if not _node(t)
                          and any(_same_test_file(t, c) for c in test_changed)]
        if whole_controls:
            return "shown"
        # The fix's scope ran (whole file only) but no changed-file control.
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
    controlled = [t for t in runs_by_target if _controlled(runs_by_target[t])]
    if any(_node(t) for t in controlled):
        return "shown"
    if controlled:
        return ("no control run: the run on "
                f"{controlled[0]} does not name a node that was "
                "seen red and green")
    never_red = [t for t in sorted(runs_by_target)
                 if runs_by_target[t][0] is None
                 and runs_by_target[t][1] is not None]
    if never_red:
        return _target_sentence(never_red[0], "was never seen failing")
    red_not_green = [t for t in sorted(runs_by_target)
                     if runs_by_target[t][0] is not None
                     and not _controlled(runs_by_target[t])]
    if red_not_green:
        return _target_sentence(
            red_not_green[0],
            "was seen failing but never passing after the change")

    # 6. A product change with no recorded run at all.
    return "no control run: no test run was recorded for the change"
