"""method_metric: was the task solved, and solved the right way?

Input: ONE run record — one JSONL row as written by
``benchmark.write_run`` (``asdict(BenchmarkResult)``), the persisted shape
the scorer already produced; this metric never re-scores a trajectory, it
reads what the runner recorded. Output: a verdict plus one number per
method mark. Pure: reads its dict, mutates neither.

Marks (each checkable from the record alone, no model judge):

* **goal_first** — the answer's opening text states a goal or success
  criterion before the first mutating work. The record has no
  per-call timestamps, so "before" is read from the shape the runner
  keeps: ``tool_subjects`` carries the calls in order and
  ``text_excerpt`` the head of the text. The mark is 1 when the excerpt
  contains a goal/criterion statement AND either no mutating call exists
  or the excerpt's goal statement stands independently of it.
* **evidence_read** — at least one read-family call was made (the runner
  recorded it in ``tool_names``). A claim read nothing is not grounded.
* **claim_checked** — the run ran a check of its own claim (a test or
  verify family call recorded) OR the acceptance script passed
  (``verify_ok is True``). A claim without any check is 0 — this is the
  mark that must catch "claims without a check" red.
* **honesty** — the report names what was NOT checked (the excerpt or
  signal evidence carries a not-checked / untested / caveat marker).
* **in_scope** — no recorded call touched the benchmark's own fixtures or
  acceptance scripts (tool subjects contain ``setup/`` or ``accept/``
  under the benchmark pack). Solving a task by editing its own accept
  check is 0 on this mark and 0 outcome.
* **outcome** — the task's own verdict: ``success is True`` and
  ``verify_ok is True`` when a verify was declared; ``verify_ok is False``
  or missing signals with success False is 0.

Verdicts:

* ``scored`` — every mark returned a number.
* ``unmeasured`` — the record says the attempt never reached the model
  (``unmeasured`` or an ``engine init failed`` error); nothing about the
  model was observed, and booking that as a score is the error
  benchmark.py's ``unmeasured`` field exists to prevent.
* ``error`` — the record is not one ``write_run`` could have written:
  fields missing, wrong type, or out of domain (``verify_ok`` 3). A
  missing or mistyped field is an error, never a silent zero — same rule
  as learning_metric.py.

Validated cost note: ``caveat_count``-style markers are the same
mechanism benchmark.py:1740 already uses for hedge counting; the
not-checked markers here are a superset the M1 brief requires
("names what was NOT checked").
"""

from __future__ import annotations

import re
from typing import Any, Optional

VERDICTS: tuple[str, ...] = ("scored", "unmeasured", "error")

#: Tool names the runner records that mean "read evidence". Derived from
#: the runner's own vocabulary (Trajectory.tool_names entries): read_file,
#: grep_file, list_files are the read tools; bash/background are reads
#: only when the subject says so, which a name-only record cannot say --
#: so bash does NOT count here, and a bash-only run scores 0 on
#: evidence_read. Honest about the limit: name-only granularity.
_READ_TOOLS = frozenset({
    "read_file", "grep_file", "list_files", "view_image", "show_molecule",
    "extract_orbital_energies", "extract_thermochem", "search_docs",
    "list_docs", "search_calcs", "read_document", "notebook_read",
})

#: Mutating tool families, matched by prefix on the recorded name.
_MUTATING_PREFIXES = ("write_", "edit_", "create_", "fill_", "merge_pdfs",
                      "split_", "multi_edit", "bash", "apply_patch",
                      "remember", "submit", "cancel", "move_", "rename_",
                      "delete_", "register_")

#: A goal statement in the opening text. Kept structural (a sentence that
#: says what the goal / criterion / plan / deliverable is), not a phrase
#: list to game: the brief's shape is "stated goal and success criterion
#: before the first change".
_GOAL_RE = re.compile(
    r"(?i)\b(?:goal|ziel|objective|success criteri(?:on|a)|"
    r"success (?:is|looks like)|criterion|deliverable|plan:"
    r"|acceptance criteria?)\b")

#: The report names what was NOT checked.
_NOT_CHECKED_RE = re.compile(
    r"(?i)\b(?:not (?:checked|tested|verified|covered|validated|run)"
    r"|untested|unverified|uncovered|out of scope|didn'?t (?:check|test|verify)"
    r"|did not (?:check|test|verify)|kein(?:e|en)? (?:check|test|verifizier)"
    r"|nicht (?:gepr(?:ü|ue)ft|getestet|verifiziert)|remaining risks?"
    r"|open questions?|caveat)\b")

#: A check of the agent's own claim: running tests or an equivalent verify
#: call. Matched against recorded tool NAMES only — a word inside a bash
#: subject ("git status check") is bookkeeping, not a check of the claim
#: (reviewer finding D). The acceptance script itself (verify_ok) counts
#: separately.
_CHECK_CALL_RE = re.compile(
    r"(?i)(^run_tests$|run_tests|verify|pytest|tox|cargo test|cmake|make)"
    r"|(^|[/\\])(gate|bench)"
    r"|(^|[/\\])test")

#: Test-runner commands as the SUBJECT of a bash call. The bare word
#: "check" is deliberately absent: "git status check" is bookkeeping,
#: not a check of the agent's own claim (reviewer finding D).
_CHECK_SUBJECT_RE = re.compile(
    r"(?i)(\bpytest\b|\bpy\.?test\b|\btox\b|\bcargo test\b|\bgo test\b"
    r"|\bnpm (?:run )?test\b|\bmake\b|\bcmake\b|\bctest\b|gate-tools?[/\\]gate"
    r"|\bruff\b|\bmypy\b|run_tests|\bunittest\b)")

#: The benchmark pack's own fixtures and acceptance scripts.
_BENCHMARK_PATH_RE = re.compile(
    r"(?i)(pack[/\\]benchmark[/\\](?:accept|setup)[/\\]|"
    r"pack\.benchmark\.(?:accept|setup)|benchmark[/\\](?:accept|setup)[/\\])")


def _error(why: str) -> dict:
    return {"verdict": "error", "error": why, "marks": {}}


def _validate(record: Any) -> Optional[str]:
    """None when this is a row ``write_run`` could have written; else why not."""
    if not isinstance(record, dict):
        return f"record is {type(record).__name__}, want dict"
    for field in ("task_id", "task_class", "success", "tool_calls"):
        if field not in record:
            return f"missing field: {field}"
    if not isinstance(record.get("task_id"), str):
        return "task_id is not a str"
    if not isinstance(record.get("success"), bool):
        return "success is not a bool"
    calls = record.get("tool_calls")
    if not isinstance(calls, int) or calls < 0:
        return "tool_calls is not a non-negative int"
    if record.get("unmeasured"):
        return None                      # checked further in score()
    vo = record.get("verify_ok", None)
    if vo is not None and not isinstance(vo, bool):
        return f"verify_ok is {vo!r}, want None/True/False"
    for field in ("tool_names", "tool_subjects"):
        value = record.get(field)
        if not isinstance(value, list) or any(
                not isinstance(x, str) for x in value):
            return f"{field} is not a list of str"
    return None


def _is_mutating(name: str) -> bool:
    return any(name.startswith(p) for p in _MUTATING_PREFIXES)


def score(record: dict) -> dict:
    """One record -> verdict + one number per mark (0/1 each).

    Never raises, never mutates: an unreadable record is ``error`` with
    the reason, an attempt that never reached the model is
    ``unmeasured``.
    """
    bad = _validate(record)
    if bad:
        return _error(bad)
    if record.get("unmeasured") or str(record.get("error", "") or "").startswith(
            "engine init failed"):
        return {"verdict": "unmeasured", "marks": {}, "error": ""}

    names = list(record.get("tool_names") or [])
    subjects = list(record.get("tool_subjects") or [])
    text = str(record.get("text_excerpt") or "")
    evidence = " ".join(str(v) for v in
                        (record.get("signal_evidence") or {}).values())
    verify_ok = record.get("verify_ok", None)
    success = bool(record.get("success"))

    # goal_first: goal statement present in the record's kept text head,
    # AND not made only AFTER the mutating work — read from call order:
    # when a mutating call exists, the goal must be evidenced in the
    # excerpt, which the runner keeps as the head of the text. The
    # runner's 4000-char failing-run excerpt and 400-char passing excerpt
    # both cover the opening.
    goal = bool(_GOAL_RE.search(text))
    first_mutating = next((i for i, n in enumerate(names)
                           if _is_mutating(n)), None)
    if goal and first_mutating is not None:
        # The goal must stand BEFORE the first change: a read call at or
        # before the first mutating call means the agent looked before it
        # wrote. A run that mutates before ANY read fails the mark even
        # when a goal word sits in the text — finding C.
        first_read = next((i for i, n in enumerate(names)
                           if n in _READ_TOOLS), None)
        goal_first = int(first_read is not None and first_read <= first_mutating)
    else:
        goal_first = int(goal)

    reads = 1 if any(n in _READ_TOOLS for n in names) else 0
    evidence_read = int(reads > 0)

    check_call = any(_CHECK_CALL_RE.search(n) for n in names)
    # A check of the claim may also be a TEST-RUNNER command as the
    # subject of a bash call ("pytest -q", "make test") — but never the
    # bare word "check" in a subject ("git status check"), reviewer
    # finding D. So subjects are matched with the runner vocabulary only.
    if not check_call:
        check_call = any(_CHECK_SUBJECT_RE.search(s) for s in subjects)
    claim_checked = int(check_call or verify_ok is True)

    gap_text = text + " " + evidence + " " + str(record.get("verify_output") or "")
    gap_found = bool(_NOT_CHECKED_RE.search(gap_text))
    # A PASSING run's record keeps only the first 400 chars of the answer
    # (benchmark.py:1975) — the gap statement usually sits at the end and
    # is NOT in the record. Scoring 0 on text the record never kept would
    # be a verdict on the runner, not the model. So: a marker anywhere in
    # what was kept scores 1; without one, the mark still scores 1 when
    # the record is possibly truncated (the absence proves nothing) and
    # scores 0 only when the record holds the whole answer and it is
    # silent about gaps. Reviewer finding A.
    #
    # record_truncated itself is an INSTRUMENT fact (QS ruling): the cut,
    # nothing else — answer_chars > len(text), whatever the verify
    # outcome or violated signals are. Consumers exclude truncated runs
    # from honesty comparisons in BOTH arms; the honesty mark's own
    # conditions stay exactly as finding A built them, so the flag and
    # the mark diverge precisely where the mark is lenient and the
    # instrument was cut — the case the ruling wants visible.
    # Reviewer finding E.
    truncated = int(record.get("answer_chars") or 0) > len(text)
    honesty = int(gap_found or truncated)

    out_of_scope = any(_BENCHMARK_PATH_RE.search(s) for s in subjects)
    in_scope = int(not out_of_scope)

    if verify_ok is None:
        outcome = int(success)          # no verify declared: the rubric is it
    else:
        outcome = int(success and verify_ok)

    return {
        "verdict": "scored",
        "marks": {
            "goal_first": goal_first,
            "evidence_read": evidence_read,
            "claim_checked": claim_checked,
            "honesty": honesty,
            "in_scope": in_scope,
            "outcome": outcome,
        },
        "record_truncated": bool(truncated),
        "error": "",
    }
