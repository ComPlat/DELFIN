"""A finished delegate left the chip list and could not be reached again.

Every delegation added a row below the message box and the row stayed for
the rest of the session, so the list grew without bound and the running
delegates were pushed down by rows with nothing left to report. The list
also carried a session-start cutoff, which meant a delegate from an
earlier session was unreachable from the UI although its saved session was
still on disk -- indistinguishable, from the reader's side, from deleted.

Finished delegates are now rows in one archive dropdown: still openable in
the same drill-in transcript view, no longer growing the live list, and no
longer scoped to the session that happened to start them.

What is asserted here: the selection rule (``archive_rows``) as a pure
function, and -- for the parts that only exist as widget wiring -- that the
tick cannot pull a reader out of an entry they opened.
"""

from __future__ import annotations

import inspect

from delfin.dashboard import tab_agent as T


def _finished(sa_id, *, kind="Explore", desc="", error=""):
    return {"sa_id": sa_id, "subagent_type": kind,
            "description": desc, "error": error}


# ---------------------------------------------------------------------------
# The selection rule
# ---------------------------------------------------------------------------

def test_a_finished_delegate_becomes_an_archive_row():
    rows = T.archive_rows({}, [_finished("a1", kind="Explore", desc="read the pack")])
    assert [sa_id for _, sa_id in rows] == ["a1"]
    label = rows[0][0]
    assert "Explore" in label and "read the pack" in label


def test_the_outcome_is_visible_without_opening_the_row():
    ok, bad = T.archive_rows(
        {}, [_finished("a1"), _finished("a2", error="boom")])
    assert ok[0].startswith("✅")
    assert bad[0].startswith("❌")


def test_a_running_delegate_is_not_also_filed():
    """The live chip list carries it; two rows would offer a stale copy."""
    rows = T.archive_rows({"a1": {"type": "Explore"}},
                          [_finished("a1"), _finished("a2")])
    assert [sa_id for _, sa_id in rows] == ["a2"]


def test_the_store_order_is_kept():
    """list_finished() returns newest first, and the dropdown shows that."""
    rows = T.archive_rows({}, [_finished("new"), _finished("mid"),
                               _finished("old")])
    assert [sa_id for _, sa_id in rows] == ["new", "mid", "old"]


def test_a_delegate_from_an_earlier_session_is_still_reachable():
    """No session-start cutoff: age bounds the list, not session identity.

    The record carries a finish time long before any plausible session
    start, and it must still produce a row -- the saved session is on disk
    and the reader asked to be able to study it.
    """
    old = _finished("a1")
    old["finished_at"] = 1.0
    assert [sa_id for _, sa_id in T.archive_rows({}, [old])] == ["a1"]


def test_a_record_without_an_id_is_skipped_not_raised():
    """The store is read defensively: a half-written record must not take
    the whole control down with it."""
    rows = T.archive_rows({}, [{}, {"sa_id": ""}, _finished("a1")])
    assert [sa_id for _, sa_id in rows] == ["a1"]


def test_nothing_finished_is_an_empty_list_not_a_placeholder_row():
    """The placeholder belongs to the widget, so the row builder stays
    honest about there being nothing to show -- that is what hides it."""
    assert T.archive_rows({}, []) == []
    assert T.archive_rows({}, None) == []


# ---------------------------------------------------------------------------
# The wiring the rule cannot cover
# ---------------------------------------------------------------------------

def _refresh_source() -> str:
    src = inspect.getsource(T)
    start = src.index("def _refresh_agent_view_chips():")
    return src[start:src.index("agent_archive_dropdown.observe")]


def test_the_tick_cannot_pull_the_reader_out_of_an_open_entry():
    """Reassigning ``options`` resets ``value`` and fires the observer.

    The 1.5 s live tick rebuilds this control, so without suppression the
    rebuild would re-enter the pick handler and send the reader back to
    Main while they were reading.
    """
    body = _refresh_source()
    assert "_view_archive_programmatic" in body
    i = body.index("agent_archive_dropdown.options =")
    guard = body[:i]
    assert "_view_archive_programmatic\"] = True" in guard, (
        "options are reassigned before the observer is suppressed")
    handler = inspect.getsource(T).split("def _on_archive_pick")[1]
    assert "_view_archive_programmatic" in handler.split("\n\n")[0]


def test_an_open_archive_entry_counts_as_a_valid_view():
    """Otherwise the next tick finds the selected id absent from the chip
    values and resets the view to Main."""
    body = _refresh_source()
    assert "vals = [v for _, v, _ in entries] + [v for _, v in archive]" in body


def test_the_archive_row_is_rebuilt_only_on_a_real_change():
    body = _refresh_source()
    assert "_view_archive_sig" in body


def test_the_control_is_in_the_layout_under_the_chips():
    src = inspect.getsource(T)
    i = src.index("         agent_view_chips,\n")
    assert "agent_archive_dropdown," in src[i:i + 400]
