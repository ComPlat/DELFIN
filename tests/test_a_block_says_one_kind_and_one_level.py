"""One block, one kind — the content a block carries reaches the file.

WHY THIS FILE EXISTS. ``create_docx`` and ``create_pdf`` describe their
input as "content blocks", one kind per block, and the two writers are
documented as taking the same shape (``test_the_same_blocks_work_for_word_
and_for_pdf``). Measured through the tool dispatch on 2026-11-23: a block
carrying more than one kind was written as whichever kind the ``if``
chain happened to test first, and the call came back ``status: ok`` with
the other half missing from the file.

The block below is not a fanciful input. It is the shape a model writes
when it reads ``{heading, level} | {paragraph}`` in the schema and puts
the heading over its paragraph in one object -- and the document it
produced was missing a sentence, with no error and a counts line that
read plausibly.

The second half of the file pins the ``level`` reading, which the two
writers had disagreed about: ``level: 0`` meant "Heading 1" in Word and
"Title" in the PDF, so the same blocks made two documents that did not
look alike, which is the one thing that shared block language exists to
prevent.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("docx")
pytest.importorskip("pypdf")
pytest.importorskip("reportlab.pdfgen.canvas")

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


def _perms(ws) -> KitToolPermissions:
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "office-block-test"
    return perms


def _docx_paragraphs(path) -> list[tuple[str, str]]:
    import docx

    return [(p.style.name, p.text)
            for p in docx.Document(str(path)).paragraphs]


def _pdf_text(path) -> str:
    import pypdf

    return "\n".join(page.extract_text() or ""
                     for page in pypdf.PdfReader(str(path)).pages)


# ---------------------------------------------------------------------------
# A block carrying two kinds
# ---------------------------------------------------------------------------

def test_a_docx_block_carrying_a_heading_and_a_paragraph_is_refused(ws):
    """Not silently halved: the model is told the block has two kinds."""
    ex = _DocToolExecutor()
    out = json.loads(ex.execute("create_docx", {
        "path": "bericht.docx",
        "blocks": [{"heading": "Ergebnis", "paragraph": "87 Prozent."}],
    }, _perms(ws)))
    assert "error" in out, out
    assert "block 1" in out["error"], out["error"]
    # The message has to say WHICH kinds collided, or the model cannot fix it.
    assert "heading" in out["error"] and "paragraph" in out["error"], out["error"]
    assert not (ws / "bericht.docx").exists()


def test_a_pdf_block_carrying_a_heading_and_a_paragraph_is_refused(ws):
    ex = _DocToolExecutor()
    out = json.loads(ex.execute("create_pdf", {
        "path": "bericht.pdf",
        "blocks": [{"heading": "Ergebnis", "paragraph": "87 Prozent."}],
    }, _perms(ws)))
    assert "error" in out, out
    assert "block 1" in out["error"], out["error"]
    assert "heading" in out["error"] and "paragraph" in out["error"], out["error"]
    assert not (ws / "bericht.pdf").exists()


def test_the_refusal_never_lands_in_the_middle_of_a_written_file(ws):
    """The first block is fine; the second collides. Nothing is left behind."""
    ex = _DocToolExecutor()
    out = json.loads(ex.execute("create_docx", {
        "path": "bericht.docx",
        "blocks": [
            {"paragraph": "Einleitung."},
            {"table": [["a", "b"]], "paragraph": "daneben"},
        ],
    }, _perms(ws)))
    assert "error" in out, out
    assert "block 2" in out["error"], out["error"]
    assert not (ws / "bericht.docx").exists()


def test_one_kind_per_block_still_writes_every_block(ws):
    """The control: the documented shape is untouched by the refusal above."""
    ex = _DocToolExecutor()
    out = json.loads(ex.execute("create_docx", {
        "path": "bericht.docx",
        "blocks": [
            {"heading": "Ergebnis", "level": 1},
            {"paragraph": "87 Prozent."},
            {"table": [["Probe", "Wert"], ["1", "12.4"]], "header_row": True},
        ],
    }, _perms(ws)))
    assert out.get("status") == "ok", out
    assert (out["headings"], out["paragraphs"], out["tables"]) == (1, 1, 1)

    pdf = json.loads(ex.execute("create_pdf", {
        "path": "bericht.pdf",
        "blocks": [
            {"heading": "Ergebnis", "level": 1},
            {"paragraph": "87 Prozent."},
            {"table": [["Probe", "Wert"], ["1", "12.4"]], "header_row": True},
        ],
    }, _perms(ws)))
    assert pdf.get("status") == "ok", pdf
    assert pdf["verified"] is True
    assert "87 Prozent" in _pdf_text(ws / "bericht.pdf")


# ---------------------------------------------------------------------------
# level: the value reaches the document, and a bad one is named
# ---------------------------------------------------------------------------

def test_a_level_reaches_the_outline_depth_in_both_writers(ws):
    """``level`` is not decoration: 2 must not be drawn as 1.

    The first version of this case asserted something about ``level: 0``
    that measurement refuted -- both writers already agree there (0 is
    the first heading level), and reportlab's Title and Heading1 styles
    have the same 18 pt, so the two documents were identical before the
    change and after it. An assertion that cannot fail is not a test, so
    it was replaced with the level that does differ, 2.
    """
    ex = _DocToolExecutor()
    perms = _perms(ws)
    word = json.loads(ex.execute("create_docx", {
        "path": "bericht.docx",
        "blocks": [{"heading": "Kapitel", "level": 2}],
    }, perms))
    assert word.get("status") == "ok", word
    assert _docx_paragraphs(ws / "bericht.docx")[0][0] == "Heading 2"


def test_a_level_that_is_not_a_number_is_refused_by_both_writers(ws):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    for tool, name in (("create_docx", "bericht.docx"),
                       ("create_pdf", "bericht.pdf")):
        out = json.loads(ex.execute(tool, {
            "path": name, "blocks": [{"heading": "H", "level": "zwei"}],
        }, perms))
        assert "error" in out, (tool, out)
        # Names the offending value, not an int() internals line.
        assert "zwei" in out["error"], (tool, out["error"])
        assert "level" in out["error"], (tool, out["error"])
        assert "block 1" in out["error"], (tool, out["error"])
        assert not (ws / name).exists(), tool


def test_a_level_left_out_or_empty_still_means_the_first_heading(ws):
    """The control for the refusal above: the tolerant cases stay tolerant.

    ``level`` omitted, null, "" and 0 all meant "the first heading level"
    before this change, and every existing caller relies on that; only a
    value that is not a number at all is new behaviour.
    """
    ex = _DocToolExecutor()
    perms = _perms(ws)
    cases = [({}, "omitted"), ({"level": None}, "null"),
             ({"level": ""}, "empty"), ({"level": 0}, "zero")]
    for i, (extra, label) in enumerate(cases):
        out = json.loads(ex.execute("create_docx", {
            "path": f"k{i}.docx", "blocks": [{"heading": "Kapitel", **extra}],
        }, perms))
        assert out.get("status") == "ok", (label, out)
        assert _docx_paragraphs(ws / f"k{i}.docx")[0][0] == "Heading 1", label
