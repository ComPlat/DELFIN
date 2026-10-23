"""Prompt sections were deleted on a test that cannot detect their use.

``ContextUsageTracker`` decides "the agent referenced section S" by
substring-matching the model's REPLY against a fixed keyword list, and
after 8 injections below a 15% hit rate the loader stops sending that
section at all.

The signal does not measure what it claims. Read off this machine's own
ledger, 500 entries:

    honesty_addendum      500 injected    0 hits     0%
    external_memory       438 injected    0 hits     0%
    live_state            265 injected    0 hits     0%
    briefing / profile /
    playbook                8 injected    0 hits     0%   ← being skipped
    repo_map              182 injected   48 hits    26%   ← kept

The honesty addendum is the proof. It scored zero over five hundred
turns, and the agent demonstrably follows it -- it refuses, it grounds
its claims, the benchmark measures a 92% verification rate. It simply
never writes the word "honesty". A section that works SILENTLY is
indistinguishable, in this counter, from one that is ignored.

The one section that survives does so for the opposite reason:
``repo_map``'s markers include the literal ``delfin/``, which appears in
almost every reply in this repository. So the mechanism reliably keeps
the largest block and deletes the ones whose value is invisible.

Nine section names have no entry in the marker table at all, so
``any(m in reply for m in [])`` is False by construction -- they are
logged as a miss every single turn without any check being performed.

And ``memory`` is not safe, it is merely early: the three memory blocks
(project memory, external memory, episodic recall) are all gated on that
one name, so they vanish together the moment the counter reaches eight.

Skipping is therefore off. The recording stays -- the data is worth
having, and whoever builds a real usage signal will need it -- but no
section is deleted from the prompt on the strength of a proxy that
cannot tell use from silence.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import context_tracker as CT


@pytest.fixture
def tracker(tmp_path):
    return CT.ContextUsageTracker(path=tmp_path / "usage.jsonl")


def _saturate(tracker, section: str, n: int = 40) -> None:
    """More than enough samples, none of them a hit."""
    for _ in range(n):
        tracker.record_usage([section], "an answer that never says the word")


# ---------------------------------------------------------------------------
# Nothing is deleted on this signal
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("section", [
    "memory", "external_memory", "episodes", "live_state",
    "briefing", "profile", "playbook", "repo_map",
    "honesty_addendum", "refusal_addendum",
])
def test_no_section_is_skipped_however_low_it_scores(tracker, section):
    _saturate(tracker, section)
    assert not tracker.should_skip(section)


def test_the_three_memory_blocks_share_one_gate(tracker):
    """They are gated on the single name "memory", so a decision about it
    is a decision about all three."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "prompt_loader.py").read_text(encoding="utf-8")
    assert src.count('_should_skip_section("memory"') >= 2


# ---------------------------------------------------------------------------
# ...and a section nobody can measure is not recorded as a failure
# ---------------------------------------------------------------------------

def test_a_section_without_markers_is_not_counted_as_a_miss(tracker):
    """`any(m in reply for m in [])` is False by construction. Recording
    that as evidence of non-use is recording nothing as something."""
    hits = tracker.record_usage(["live_state"], "any reply at all")
    assert hits.get("live_state") is None


def test_a_section_with_markers_is_still_measured(tracker):
    hits = tracker.record_usage(["playbook"], "according to the playbook, step 1")
    assert hits.get("playbook") is True


def test_a_measurable_section_that_misses_is_recorded_as_a_miss(tracker):
    hits = tracker.record_usage(["playbook"], "nothing relevant here")
    assert hits.get("playbook") is False


def test_the_unmeasurable_section_leaves_no_miss_in_the_ledger(tracker):
    tracker.record_usage(["live_state", "playbook"], "nothing relevant")
    entry = json.loads((tracker._path).read_text(encoding="utf-8").splitlines()[0])
    assert "live_state" not in entry["hits"]
    assert entry["hits"]["playbook"] is False


def test_which_sections_are_measurable_is_stated_not_guessed():
    measurable = {s for s, m in CT._SECTION_MARKERS.items() if m}
    assert "playbook" in measurable
    assert "live_state" not in measurable
    assert "honesty_addendum" not in measurable


# ---------------------------------------------------------------------------
# The recording itself is untouched — the data stays worth having
# ---------------------------------------------------------------------------

def test_usage_is_still_recorded(tracker):
    tracker.record_usage(["repo_map"], "see delfin/agent/engine.py")
    assert tracker._path.exists()


def test_the_rate_report_still_works(tracker):
    tracker.record_usage(["playbook"], "according to the playbook")
    tracker.record_usage(["playbook"], "nothing")
    rates = tracker.get_hit_rates()
    assert rates["playbook"] == pytest.approx(0.5)


def test_a_reenable_switch_exists_and_defaults_off():
    """Whoever builds a real usage signal should not have to re-derive
    this; the decision is a constant, not a deletion."""
    assert CT.SKIP_UNUSED_SECTIONS is False
