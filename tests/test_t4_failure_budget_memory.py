"""Phase 3 follow-up: FailureBudget is memory-bounded (operator review).

Wave 13, package T4 (nacht-s18). The operator found a defect in the first
FailureBudget draft: it keyed every call on ``(tool, normalize_args(args))``
and on a success only ZEROED the count, leaving the key in the table. Because
the dispatch wiring makes the budget CLASS-LEVEL on _DocToolExecutor (it lives
for the whole process — one shared dashboard kernel, several sessions), the
table grew without bound: for ``write_file`` the *whole lowercased file
content* became a dict key, so every distinct call ever made stayed resident.

The operator's required fix, verbatim:
  1. key on a short hash of the normalized signature (sha1 hexdigest),
  2. DELETE the key on success instead of storing 0,
  3. cap the table (512 entries, drop the oldest).

These are PURE unit tests of action_grounding.FailureBudget: no I/O, no
executor, no dispatch. ``entry_count`` is the public observation the bounded
contract exposes (the wiring may log it too).
"""

from delfin.agent.action_grounding import (
    FailureBudget,
    FailureBudgetLimits,
    GroundingHint,
)


def test_10000_distinct_successes_leave_the_table_empty():
    # The fix must NOT retain a key per distinct call after a success. A
    # success ends a call's lifecycle; its key is deleted. Distinct
    # successful calls (including large write payloads) must never accumulate.
    b = FailureBudget()
    for i in range(10000):
        b.record("write_file", {"path": f"f{i}.py", "content": "x" * 500}, ok=True)
    assert b.entry_count == 0


def test_10000_distinct_failures_stay_at_or_below_the_cap():
    # Distinct failing calls may accumulate, but only up to the cap. Past
    # max_entries the oldest keys must be dropped, so memory is bounded.
    b = FailureBudget(limits=FailureBudgetLimits(max_entries=512))
    for i in range(10000):
        b.record("read_file", {"path": f"f{i}.py"}, ok=False)
    assert b.entry_count <= 512


def test_third_identical_failure_still_gives_the_hint_when_capped():
    # The bounded table must not break the core stop-and-ask contract: even
    # after the table has been filled and evicted, a third identical failure
    # in a row still yields the repeated_failure hint.
    b = FailureBudget(limits=FailureBudgetLimits(max_entries=512))
    for i in range(512):
        b.record("read_file", {"path": f"pad{i}.py"}, ok=False)   # fill the cap
    b.record("edit_file", {"path": "real.py"}, ok=False)          # streak 1
    b.record("edit_file", {"path": "real.py"}, ok=False)          # streak 2
    assert b.repeat_hint("edit_file", {"path": "real.py"}) is None
    b.record("edit_file", {"path": "real.py"}, ok=False)          # streak 3
    hit = b.repeat_hint("edit_file", {"path": "real.py"})
    assert isinstance(hit, GroundingHint)
    assert hit.kind == "repeated_failure"
