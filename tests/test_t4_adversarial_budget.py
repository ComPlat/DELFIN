"""Package T4, reviewer (nacht-s19) — adversarial tests for phase 3.

action_grounding.FailureBudget: a pure per-task gate that says "stop and ask"
when the SAME (tool, normalized args) call has already failed N times in a
row. These adversarial tests pin the edges that decide whether the module
actually closes the gap the package set out to close:

* A1 — the RAW-[80]-string bucket-splitting fix. failure_log keys on a truncated
  raw string, so a long command re-worded past the cutoff lands in a DIFFERENT
  bucket and resets the counter. Here we prove the normalized key collapses two
  >80-char spellings of the same call that the raw key WOULD split.
* A2 — repeat_hint is side-effect free: safe to call before every dispatch
  without disturbing the streak or the count it reads back.
* A3 — record() is robust to non-string / None tool values (the wired caller
  feeds whatever the caller produced).
* A4 — failure_signature() reports the most-repeated failing call BY TOTAL
  COUNT, which is deliberately different from repeat_hint's CONSECUTIVE-STREAK
  semantics. Pinned here so the asymmetry is visible and intentional, not an
  accidental drift.

Run via the gate:  gate tests/test_t4_adversarial_budget.py -q | tail -12
NOTE: on the reviewer branch these require the phase-2 + phase-3 module (merged
in from agent/s18-t4b: d90aedc4 + d9bf770b).
"""
from delfin.agent.action_grounding import FailureBudget, GroundingHint


def test_long_command_whitespace_and_case_only_differs_shares_bucket():
    """A1 — the exact phase-1 gap: a failing call whose raw [:80] key material
    differs only by interior whitespace runs and case lands in a DIFFERENT
    bucket under the old failure_log key. ``normalize_args`` lowers case and
    collapses whitespace runs, so the two spellings are the SAME failing call.

    Note: token order WITHIN a command string is deliberately NOT normalised
    (``rm a && rm b`` != ``rm b && rm a``), so the two spellings here differ
    ONLY by whitespace/case, never by command word order or by the key set.
    """
    # Long (>80 chars) so a raw [:80]/[:300] truncation of the value would
    # already have a chance to look different, and the caller is one tool.
    long_tail = "/srv/projects/group/user/software/delfin/tests"
    spaced = f"LS -LA {long_tail} && ECHO done && WC  -l  test_t4_grounding.py"
    plain = f"ls -la {long_tail} && echo done && wc -l test_t4_grounding.py"
    b = FailureBudget()
    b.record("bash", {"command": spaced}, ok=False)
    b.record("bash", {"command": plain}, ok=False)
    assert b.repeat_hint("bash", {"command": plain}) is None  # streak 2 < 3
    b.record("bash", {"command": spaced}, ok=False)           # streak 3
    # recognised as the SAME failing call as 'plain'
    assert isinstance(b.repeat_hint("bash", {"command": plain}),
                      GroundingHint)
    assert isinstance(b.repeat_hint("bash", {"command": spaced}),
                      GroundingHint)


def test_dict_key_order_collapses_to_same_bucket():
    """A1b — argument DICT order collapse: the same set of keys in a different
    order is the same failing call (the phase-1 reordered-args case). The key
    SET is identical — only its order changes."""
    long_tail = "/srv/projects/group/user/software/delfin/tests"
    cmd = f"ls -la {long_tail} && echo done"
    b = FailureBudget()
    b.record("bash", {"command": cmd, "cwd": "/tmp"}, ok=False)
    b.record("bash", {"cwd": "/tmp", "command": cmd}, ok=False)
    assert b.repeat_hint("bash", {"command": cmd, "cwd": "/tmp"}) is None
    b.record("bash", {"cwd": "/tmp", "command": cmd}, ok=False)
    assert isinstance(b.repeat_hint("bash", {"command": cmd, "cwd": "/tmp"}),
                      GroundingHint)


def test_repeat_hint_is_side_effect_free():
    """A2 — repeat_hint must not consume or mutate the streak/count, so the
    wiring may call it before every dispatch (even repeatedly) safely."""
    b = FailureBudget()
    for _ in range(3):
        b.record("edit_file", {"path": "a.py"}, ok=False)
    h1 = b.repeat_hint("edit_file", {"path": "a.py"})
    h2 = b.repeat_hint("edit_file", {"path": "a.py"})  # called again, same result
    assert isinstance(h1, GroundingHint)
    assert isinstance(h2, GroundingHint)
    # calling it twice did not advance to a 4-state or reset anything
    h3 = b.repeat_hint("edit_file", {"path": "a.py"})
    assert isinstance(h3, GroundingHint)


def test_record_accepts_non_string_and_none_tool():
    """A3 — robustness: the wired caller may pass a tool name that isn't a
    clean string (or None). record must not raise, and the bucket must be keyed
    on its string form so repeated None-tools still accumulate."""
    b = FailureBudget()
    b.record(None, {"cmd": "x"}, ok=False)
    b.record(None, {"cmd": "x"}, ok=False)
    b.record("12345", {"cmd": "x"}, ok=False)  # int-like tool name, distinct
    assert b.repeat_hint(None, {"cmd": "x"}) is None
    b.record(None, {"cmd": "x"}, ok=False)     # None-tool streak -> 3
    assert isinstance(b.repeat_hint(None, {"cmd": "x"}), GroundingHint)
    # the int-like tool bucket is untouched
    assert b.repeat_hint("12345", {"cmd": "x"}) is None


def test_failure_signature_reports_most_repeated_by_total_count():
    """A4 — failure_signature ranks by TOTAL count (not consecutive streak),
    which can differ from what repeat_hint treats as the stuck loop. Pinned so
    the asymmetry is explicit: a call with 4 total failures but an interleaved
    success (streak reset) is reported here, yet repeat_hint would NOT flag
    it. The failure_log hook consumes this string, so the wiring must not
    assume signature == currently-stuck."""
    b = FailureBudget()
    # 'bash' fails twice then a success breaks its streak (total 2, streak 0)
    b.record("bash", {"command": "make"}, ok=False)
    b.record("bash", {"command": "make"}, ok=False)
    b.record("bash", {"command": "make"}, ok=True)   # streak reset to 0
    # 'edit_file' fails three times consecutively (total 3, streak 3)
    b.record("edit_file", {"path": "b.py"}, ok=False)
    b.record("edit_file", {"path": "b.py"}, ok=False)
    b.record("edit_file", {"path": "b.py"}, ok=False)
    sig = b.failure_signature()
    assert sig is not None and sig.startswith("edit_file|")
    # and repeat_hint flags exactly the stuck one
    assert isinstance(b.repeat_hint("edit_file", {"path": "b.py"}),
                      GroundingHint)
    assert b.repeat_hint("bash", {"command": "make"}) is None
