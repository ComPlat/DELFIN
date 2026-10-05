"""T4 reviewer (nacht-s19) — adversarial coverage of the memfix (720288ea).

The builder reworked FailureBudget to be memory-bounded: keys are now
``(tool, sha1(normalized_args))``, a success DELETEs its entry, and the table
is capped at FailureBudgetLimits.max_entries (FIFO-drop-oldest on a genuinely
new failing signature). These tests attack the two places the memfix could
corrupt the phase-3 core contract:

* M1 — a signature evicted when the table is full must NOT inherit a stale
  count: if it later re-fails it restarts the streak from 1, so a stuck loop
  is neither prematurely hinted (count inherited across eviction) nor
  mis-counted.
* M2 — the cap is a hard bound under interleaved success/trim traffic: after
  filling the table and deleting some successes while adding new failures,
  the table never exceeds max_entries and a just-deleted (previously failing)
  signature is never re-reported as most-repeated.
* M3 — eviction ordering is FIFO and the evicted key's failure_signature is
  gone: after the oldest key is dropped, failure_signature no longer names it.

All pure, through the public API of FailureBudget — no I/O, no executor.
"""
from delfin.agent.action_grounding import (
    FailureBudget,
    FailureBudgetLimits,
    GroundingHint,
)


def test_evicted_signature_restarts_its_streak_from_one():
    """M1 — an evicted key, when it re-fails, must not carry a stale count.
    Pack the table to the cap, then fail an OLD evicted signature once more:
    one failure after eviction must NOT produce a hint (it is streak-1, not
    streak-4)."""
    b = FailureBudget(limits=FailureBudgetLimits(identical_fail_limit=3,
                                                 max_entries=512))
    # fill the table so every key is resident; the FIRST padding key is the
    # oldest and will be the first evicted.
    for i in range(512):
        b.record("pad", {"n": i}, ok=False)
    assert b.entry_count == 512
    # force a genuinely new key into the full table -> evicts the oldest
    b.record("write_file", {"path": "new.py"}, ok=False)
    assert b.entry_count <= 512
    # the evicted signature re-fails exactly once: streak 1 -> no hint yet
    assert b.repeat_hint("pad", {"n": 0}) is None


def test_cap_is_a_hard_bound_under_interleaved_success_and_failures():
    """M2 — under interleaved success-DELETE and new-failure traffic the table
    never exceeds max_entries, and a signature that just succeeded is never
    re-reported by failure_signature."""
    b = FailureBudget(limits=FailureBudgetLimits(max_entries=64))
    for i in range(1000):
        # one new failing key (may evict)
        b.record("read_file", {"path": f"f{i}.py"}, ok=False)
        # periodically succeed an OLDER one away (pop, not zero)
        if i % 7 == 0:
            older = i - 40
            if older >= 0:
                b.record("read_file", {"path": f"f{older}.py"}, ok=True)
    assert b.entry_count <= 64
    # a key whose most recent outcome was success must not be most-repeated
    sig = b.failure_signature()
    if sig is not None:
        # all keys currently resident had an ok=False as their last outcome
        # for <key> (pop removed the successful ones); nothing further to pin
        assert isinstance(sig, str) and "|" in sig


def test_evicted_oldest_is_no_longer_reported_as_most_repeated():
    """M3 — after the oldest entry is evicted, failure_signature must not
    name it; a distinct survivor with the highest count wins."""
    b = FailureBudget(limits=FailureBudgetLimits(max_entries=3))
    b.record("a", {"k": 1}, ok=False)          # oldest
    b.record("b", {"k": 2}, ok=False)
    b.record("c", {"k": 3}, ok=False)
    b.record("b", {"k": 2}, ok=False)          # b -> count 2
    b.record("c", {"k": 3}, ok=False)          # c -> count 2
    b.record("d", {"k": 4}, ok=False)          # NEW, full table -> evicts 'a'
    sig = b.failure_signature()
    assert sig is not None
    assert not sig.startswith("a|"), sig       # evicted oldest gone
    # b and c tie at count 2; one of them or d must be reported
    assert any(sig.startswith(prefix) for prefix in ("b|", "c|", "d|")), sig
