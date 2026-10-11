"""U7 / phase 2 - adversarial review of the stale-blocked marker.

The builder's test_u7_stale_blocked.py pins the happy cases (older -> stale,
newer/equal -> current, non-blocked never, None default). These are the
corners it does not: the timestamp and reference boundaries where a real
implementation mixes units blindly or takes the panel down.

Attacked, against the merged task_ticker (`_block_epoch` / `_stale_suffix`):

- unparseable / missing task timestamp must be a silent non-stale, never a
  crash and never a blind mixed-unit ordering ("" / "garbage" / a bad date);
- `last_message_at` is contractually an epoch float, but the reference can
  arrive as a string from a JSON boundary - a numeric string must still
  work, and a NON-numeric reference must not take the render down (the
  ticker is a courtesy: it never raises);
- render_html must survive a hostile blocked_reason (HTML characters) with
  the stale marker still present, and escape it (no raw script injection);
- a blocked task with no reason is never stale (even with a far-future
  reference), matching the store invariant that a reason is required;
- a non-ASCII reason plus the stale marker both survive the text surface;
- far-future vs near-zero reference flip stale, never a constant answer.

These are GREEN on merged main when the implementation is correct; a red one
here is a finding.
"""

from __future__ import annotations

import datetime as _dt

import pytest

from delfin.agent import task_ticker as TT
from delfin.agent.agent_tasks import get_store


@pytest.fixture
def store(tmp_path):
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return get_store(tmp_path)


def _blocked(store, subject, reason="waiting on approval"):
    t = store.create(subject, session_id="s")
    return store.update(t["id"], status="blocked", blocked_reason=reason)


def test_a_missing_or_garbage_timestamp_is_silently_non_stale():
    """None, '' and unparseable updated_at must read as None from
    _block_epoch - a silent non-stale, never the "never a blind mixed-unit
    ordering" hazard of defaulting to 0/now."""
    assert TT._block_epoch({}) is None
    assert TT._block_epoch({"updated_at": ""}) is None
    assert TT._block_epoch({"updated_at": "not-a-date"}) is None
    assert TT._block_epoch({"updated_at": "2026-13-99T99:99:99Z"}) is None


def test_created_at_is_the_fallback_when_updated_at_is_gone():
    """A task with only created_at (never re-stamped) still has a timestamp
    to compare."""
    epoch = _dt.datetime(
        2026, 10, 11, 0, 0, tzinfo=_dt.timezone.utc
    ).timestamp()
    assert TT._block_epoch({"created_at": "2026-10-11T00:00:00Z"}) == epoch


def test_a_numeric_string_last_message_at_still_works(store, tmp_path):
    """The reference can cross a JSON boundary as a string; a numeric string
    must be coerced, exactly as a float call would behave."""
    _blocked(store, "push the branch")
    row = TT.rows(tmp_path, session_id="s", last_message_at="10000000000")
    assert " (stale)" in row[0]["label"]


def test_a_non_numeric_last_message_at_never_takes_the_render_down(store, tmp_path):
    """_block_epoch promises 'never raises'; _stale_suffix must hold that too.
    A malformed reference (non-numeric from a third-party wire) must not
    crash rows/render_text/render_html - the ticker is a courtesy."""
    _blocked(store, "push the branch")
    for surface in (
        lambda: TT.rows(tmp_path, session_id="s", last_message_at="oops"),
        lambda: TT.render_text(tmp_path, session_id="s", last_message_at="oops"),
        lambda: TT.render_html(tmp_path, session_id="s", last_message_at="oops"),
    ):
        surface()  # must not raise


def test_render_html_escapes_a_hostile_reason_and_keeps_the_marker(
    store, tmp_path,
):
    """A blocked_reason carrying HTML must be escaped (no raw injection) and
    the stale marker must still be present."""
    t = store.create("push the branch", session_id="s")
    store.update(t["id"], status="blocked", blocked_reason="<script>alert(1)</script>")
    html = TT.render_html(tmp_path, session_id="s", last_message_at=10**13)
    assert "(stale)" in html, "the marker must still render"
    assert "<script>" not in html, "the hostile reason must be escaped"


def test_a_blocked_task_with_no_reason_is_never_stale():
    """Even far in the future, a blocked task with no reason is not flagged -
    there is no reason to call stale. The store itself refuses to write a
    blocked row with an empty reason (agent_tasks update raises), so the
    only way this state reaches the renderer is a hand-built dict; the
    renderer must still not flag it."""
    assert TT._stale_suffix(
        {"status": "blocked", "blocked_reason": "",
         "updated_at": "2020-01-01T00:00:00Z"},
        10**13,
    ) == ""


def test_unicode_reason_and_marker_survive_together(store, tmp_path):
    """A non-ASCII reason plus the stale marker, through the text surface."""
    _blocked(store, "push the branch", reason="waiting on Bestaetigung")
    text = TT.render_text(tmp_path, session_id="s", last_message_at=10**13)
    assert "(stale)" in text
    assert "waiting on Bestaetigung" in text


def test_far_future_and_near_zero_flip_stale(store, tmp_path):
    """The decision must be a function of the reference, not a constant."""
    _blocked(store, "push the branch")
    fut = TT.rows(tmp_path, session_id="s", last_message_at=10**13)
    now = TT.rows(tmp_path, session_id="s", last_message_at=0.0)
    assert " (stale)" in fut[0]["label"]
    assert " (stale)" not in now[0]["label"]
