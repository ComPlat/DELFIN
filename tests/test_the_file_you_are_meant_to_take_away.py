"""The two artifacts most worth handing over produced no card.

The chat surfaces what a turn produced: a plot, a report, a workbook, a
zip — rendered inline with a download control, so the user does not have
to leave the conversation to get the thing they asked for.

Two producers were missing from the list that triggers it, and they are
the two whose whole purpose is that the user takes the file away.
`draft_email` writes a .eml that exists precisely so it is opened and
sent; `fill_series` writes one document per row of a table. Neither
appeared. On top of that the renderer had no branch for .eml at all, so
even naming it in the answer produced nothing.

Both were invisible for the same reason five defects were earlier the
same day: neither tool had ever been called in a benchmark run, so nobody
ever saw the card that did not appear.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.dashboard.tab_agent import (
    _FILE_CREATING_TOOLS, _is_file_creating, _render_artifact_body,
    _render_artifact_inline,
)


@pytest.fixture
def draft(tmp_path):
    from delfin.agent import office

    path = tmp_path / "anfrage.eml"
    office.draft_email(path, to="einkauf@example.org",
                       cc="chef@example.org",
                       subject="Rückfrage zu Beleg R-014",
                       body="Guten Tag,\n\nzu R-014 fehlt der Betrag.")
    return path


# ---------------------------------------------------------------------------
# The trigger
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("tool", ["draft_email", "fill_series"])
def test_the_producers_that_hand_something_over_trigger_a_card(tool):
    assert _is_file_creating(tool), (
        f"{tool} writes a file for the user to keep and surfaces nothing")


@pytest.mark.parametrize("tool", ["draft_email", "fill_series"])
def test_the_trigger_works_under_both_tool_vocabularies(tool):
    """A set holding one spelling silently does nothing on the other
    backend; _is_file_creating strips the transport namespace so the
    concept is what matters."""
    assert _is_file_creating(f"mcp__kit-coding__{tool}")


def test_every_document_producer_in_the_catalogue_is_in_the_set():
    """The general form. A producer added later and forgotten here is a
    file the user is never offered."""
    from delfin.agent import api_client as A

    producers = set()
    for entry in A._DOC_TOOLS_OPENAI:
        fn = entry["function"]
        props = set(fn.get("parameters", {}).get("properties", {}))
        if {"output", "output_dir"} & props:
            producers.add(fn["name"])
    producers |= {"create_pdf", "create_docx", "draft_email", "edit_sheet",
                  "fill_series", "publish_report"}
    missing = sorted(n for n in producers if not _is_file_creating(n))
    assert not missing, (
        "these leave a file behind and surface no card: " + ", ".join(missing))


# ---------------------------------------------------------------------------
# The card
# ---------------------------------------------------------------------------

def test_a_draft_is_rendered_at_all(draft):
    assert _render_artifact_body(draft) is not None


def test_the_card_says_who_it_is_to_and_what_it_is_about(draft):
    html = _render_artifact_inline(draft)
    assert html
    assert "einkauf@example.org" in html
    assert "chef@example.org" in html


def test_an_encoded_subject_is_decoded_for_the_reader(draft):
    """A German subject travels RFC-2047 encoded, and
    =?utf-8?q?R=C3=BCckfrage?= tells nobody anything."""
    raw = draft.read_text(encoding="utf-8", errors="replace")
    assert "=?utf-8?" in raw, "premise gone: the subject is not encoded"
    html = _render_artifact_inline(draft)
    assert "Rückfrage zu Beleg R-014" in html
    assert "=?utf-8?" not in html


def test_the_card_says_it_will_not_be_sent(draft):
    """The tool never sends. A card that looked like an outbox would be
    worse than no card."""
    html = _render_artifact_inline(draft)
    assert re.search(r"(?i)nicht versendet|selbst senden", html)


def test_the_draft_can_actually_be_downloaded(draft):
    html = _render_artifact_inline(draft)
    assert "download=" in html and "anfrage.eml" in html


def test_a_series_document_still_renders(tmp_path):
    """fill_series writes .docx per row; the trigger is only worth
    anything if the renderer draws what it produces."""
    from delfin.agent import office

    out = tmp_path / "brief.docx"
    office.create_docx(out, [{"heading": "Zahlungserinnerung"},
                             {"paragraph": "Sehr geehrte Meier GmbH,"}])
    html = _render_artifact_inline(out)
    assert html and "download=" in html


def test_an_unreadable_draft_does_not_take_the_chat_down(tmp_path):
    broken = tmp_path / "broken.eml"
    broken.write_bytes(b"\x00\x01\x02 not a message")
    _render_artifact_body(broken)  # must not raise
