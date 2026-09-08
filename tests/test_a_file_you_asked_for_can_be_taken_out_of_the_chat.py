"""An artifact you can see, and an artifact you can keep.

The chat drew the types it knew how to draw — images, SVG, CSV, JSON,
structures — and returned nothing for everything else. So a PDF, a Word
file, a workbook or a zip, which is exactly what ``create_pdf``,
``create_docx``, ``merge_pdfs`` and ``publish_report`` produce, appeared
in the conversation as nothing at all: the agent said the report was
ready and the report was nowhere.

And nothing that WAS drawn could be saved. Voila serves the notebook,
not the workspace, so a ``file://`` or relative link resolves to nothing
in the browser — the reason images are embedded rather than linked. The
download follows the same reasoning: the bytes ride inside the page and
the browser's own ``download`` attribute writes the file, so there is no
route to serve and nothing leaves the machine that was not already on
screen.

What is deliberately still silent: source and text files. They are in
the workspace panel and quoted in the answer, and a download button under
every scratch .py the agent writes is noise, not a feature.
"""

from __future__ import annotations

import base64
import re
import zipfile

import pytest

from delfin.dashboard.tab_agent import (_ARTIFACT_DOWNLOAD_LIMIT,
                                        _render_artifact_inline)

_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQ"
    "DwAEhQGAhKmMIQAAAABJRU5ErkJggg==")


def _href(html: str) -> str:
    m = re.search(r'href="data:([^;]+);base64,([^"]+)"', html)
    assert m, f"no download control in: {html[:400]}"
    return m.group(1), base64.b64decode(m.group(2))


def test_an_image_can_now_be_saved_as_well_as_seen(tmp_path):
    p = tmp_path / "plot.png"
    p.write_bytes(_PNG)
    html = _render_artifact_inline(p)
    assert 'download="plot.png"' in html
    mime, data = _href(html)
    assert data == _PNG, "the download would write different bytes"


def test_the_control_sits_inside_the_card_it_belongs_to(tmp_path):
    """A stray button between two artifacts belongs to neither."""
    p = tmp_path / "plot.png"
    p.write_bytes(_PNG)
    html = _render_artifact_inline(p).rstrip()
    assert html.endswith("</div>")
    assert html.index("download=") < html.rindex("</div>")


def test_a_pdf_is_shown_and_not_merely_announced(tmp_path):
    fitz = pytest.importorskip("fitz")
    p = tmp_path / "report.pdf"
    doc = fitz.open()
    page = doc.new_page()
    page.insert_text((72, 72), "Ergebnisbericht")
    doc.save(str(p))
    doc.close()

    html = _render_artifact_inline(p)
    assert html is not None, "a PDF produced no card at all"
    assert "page 1 of 1" in html
    assert 'download="report.pdf"' in html
    mime, data = _href(html)
    assert mime == "application/pdf"
    assert data.startswith(b"%PDF")


def test_a_workbook_says_what_is_in_it(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    p = tmp_path / "messwerte.xlsx"
    wb = openpyxl.Workbook()
    wb.active.title = "Proben"
    wb.active.append(["id", "wert"])
    wb.active.append([1, 2.5])
    wb.save(str(p))
    html = _render_artifact_inline(p)
    assert html is not None
    assert "Proben" in html
    assert 'download="messwerte.xlsx"' in html


def test_a_document_says_what_is_in_it(tmp_path):
    docx = pytest.importorskip("docx")
    p = tmp_path / "protokoll.docx"
    d = docx.Document()
    d.add_paragraph("Versuchsprotokoll vom Montag")
    d.save(str(p))
    html = _render_artifact_inline(p)
    assert html is not None
    assert "Versuchsprotokoll" in html
    assert 'download="protokoll.docx"' in html


def test_an_archive_lists_itself_before_you_open_it(tmp_path):
    p = tmp_path / "ergebnisse.zip"
    with zipfile.ZipFile(p, "w") as zf:
        zf.writestr("a/spektrum.csv", "x,y\n1,2\n")
        zf.writestr("a/notiz.txt", "hallo")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "spektrum.csv" in html and "notiz.txt" in html
    assert 'download="ergebnisse.zip"' in html


def test_a_scratch_source_file_stays_silent(tmp_path):
    """The noise floor. Everything the agent writes while working would
    otherwise get a card of its own."""
    for name in ("tagreport.py", "notes.md", "out.log", "run.sh"):
        p = tmp_path / name
        p.write_text("print('hi')\n")
        assert _render_artifact_inline(p) is None, name


def test_a_file_too_large_to_carry_is_named_instead_of_truncated(tmp_path):
    """A partial download is worse than none, and the path is what the
    user needs in order to reach it another way."""
    p = tmp_path / "huge.zip"
    with p.open("wb") as fh:
        fh.write(b"PK\x05\x06" + b"\0" * (_ARTIFACT_DOWNLOAD_LIMIT + 1))
    html = _render_artifact_inline(p)
    assert html is not None
    assert "download=" not in html
    assert str(p) in html


def test_a_broken_document_still_offers_the_file(tmp_path):
    """A preview that cannot be built must not swallow the artifact —
    that is the failure this whole card exists to end."""
    p = tmp_path / "kaputt.docx"
    p.write_bytes(b"not really a document")
    html = _render_artifact_inline(p)
    assert html is not None
    assert 'download="kaputt.docx"' in html


# ---------------------------------------------------------------------------
# ...and the scan has to run at all
# ---------------------------------------------------------------------------

def test_the_scan_fires_on_the_names_this_backend_actually_uses():
    """The gate read `tool_name in ("Write", "Edit", "Bash",
    "NotebookEdit")`. Those are the Anthropic spellings; the KIT and
    Ollama surfaces call the same tools `write_file`, `edit_file` and
    `bash`. So on the surface these models are served from, the artifact
    scan never ran — the agent wrote the file and the chat showed
    nothing, whatever the renderer could have drawn."""
    from delfin.dashboard.tab_agent import _is_file_creating

    for name in ("write_file", "edit_file", "multi_edit", "apply_patch",
                 "bash", "Write", "Edit", "Bash", "NotebookEdit"):
        assert _is_file_creating(name), name


def test_the_tools_that_produce_the_formats_a_user_asks_for_count_too():
    """create_pdf, create_docx, merge_pdfs and publish_report were on
    neither list, so the one case where the user explicitly asked to be
    GIVEN a file was the case that stayed silent."""
    from delfin.dashboard.tab_agent import _is_file_creating

    for name in ("create_pdf", "create_docx", "merge_pdfs", "split_pdf",
                 "fill_docx_template", "fill_pdf_form", "edit_sheet",
                 "publish_report"):
        assert _is_file_creating(name), name


def test_an_mcp_namespaced_tool_is_recognised_by_its_bare_name():
    from delfin.dashboard.tab_agent import _is_file_creating

    assert _is_file_creating("mcp__delfin-ops__plot_uvvis_spectrum")
    assert _is_file_creating("mcp__kit-coding__write_file")


def test_a_read_only_tool_does_not_trigger_a_directory_diff():
    """The gate exists so the scan is not run after every tool call."""
    from delfin.dashboard.tab_agent import _is_file_creating

    for name in ("read_file", "grep_file", "list_files", "web_search",
                 "search_docs", "mcp__delfin-docs__read_section"):
        assert not _is_file_creating(name), name


def test_the_same_report_written_again_is_a_new_artifact(tmp_path):
    """The turn-start snapshot keyed on the file NAME, so a report was
    shown the first time it was written and never again — and "mach den
    Bericht nochmal mit den korrigierten Zahlen" is the second time,
    which is the one being waited for."""
    from delfin.dashboard.tab_agent import _artifact_stamp

    p = tmp_path / "bericht.pdf"
    p.write_bytes(b"%PDF-1.4 first")
    before = _artifact_stamp(p)
    p.write_bytes(b"%PDF-1.4 second, with the corrected numbers")
    assert _artifact_stamp(p) != before


def test_a_file_that_vanished_does_not_take_the_scan_down(tmp_path):
    from delfin.dashboard.tab_agent import _artifact_stamp

    assert _artifact_stamp(tmp_path / "gone.pdf") == (0, -1)
