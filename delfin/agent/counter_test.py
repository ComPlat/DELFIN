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
"""

from __future__ import annotations

from typing import Dict, Iterable, Optional, Tuple


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


def assess(changed_files: Iterable, test_runs: Iterable) -> Optional[str]:
    """The report line for the turn's control-run evidence.

    ``changed_files``: paths the turn changed (product AND test files).
    ``test_runs``: recorded runs from the DURABLE ledger — every run
    ``{"target": <pytest target or node id>, "red": bool,
    "green": bool}`` (a whole-file run may use ``file``). A control is
    on record when the SAME target was seen red AND green.

    Returns:
    * ``None`` — the turn changed no product file; nothing is added.
    * ``"shown"`` — a control run is on record: some test target was
      seen both red and green.
    * one sentence naming what is missing otherwise. A completed control
      wins over a green-only unrelated run, so a properly-controlled fix
      is not sentenced just because the turn also holds a stale
      green-only run.

    Precedence (documented, not guessed):
    1. a changed test file that no run target covers — the fix's own
       counter-test sitting unrun;
    2. a target seen green but never red (the classic "passes with and
       without the fix");
    3. a target seen red but never green after the change;
    4. otherwise — a product change with no recorded run at all.
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

    # A complete control on some TARGET is the strongest evidence.
    if any(red and green for red, green in runs_by_target.values()):
        return "shown"

    # 2/3. No complete control. Name the first offender in deterministic
    # (sorted) order — the relative order between the two halves is
    # explicit: a never-seen-red target first, then red-but-not-green.
    never_red = [t for t in sorted(runs_by_target)
                 if not runs_by_target[t][0]]
    if never_red:
        return (f"no control run: test {never_red[0]}"
                " was never seen failing")
    red_not_green = [t for t in sorted(runs_by_target)
                     if runs_by_target[t][0] and not runs_by_target[t][1]]
    if red_not_green:
        return (f"no control run: test {red_not_green[0]}"
                " was seen failing but never passing after the change")

    # 4. A product change with no recorded run at all.
    return "no control run: no test run was recorded for the change"
