"""Enter makes a paragraph, and the size and the font reach the file.

Editing was block by block: pressing Enter dropped out of the paragraph, so a
document could be corrected but not written.  And the size control changed the
page and nothing else -- the run's size and font were dropped on the way to the
kernel, so the file came back in the size it had.
"""

from __future__ import annotations

import json

import pytest

from delfin.dashboard import docx_view
from delfin.dashboard.context import DashboardContext

docx = pytest.importorskip("docx")


@pytest.fixture
def opened(tmp_path):
    from delfin.dashboard import tab_calculations_browser as browser

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    doc = docx.Document()
    doc.add_paragraph("First paragraph.")
    doc.add_paragraph("Second paragraph.")
    path = tmp_path / "calc" / "report.docx"
    doc.save(path)

    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = sent.append
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()
    state = refs["sheet_state"]
    state["docx_path"] = str(path)
    state["docx_doc"] = docx_view.read_document(path)
    state["docx_edit"] = {"active": True}
    return refs, state, path


def _message(refs, **payload):
    refs["calc_sheet_payload_input"].value = json.dumps({"kind": "docx", **payload})
    refs["calc_sheet_action_btn"].click()


def _texts(path):
    return [p.text for p in docx.Document(str(path)).paragraphs]


def test_enter_adds_a_paragraph_where_it_was_pressed(opened):
    refs, state, path = opened

    _message(refs, address="new:1", after="p:0",
             runs=[{"t": "paragraph.", "b": 0, "i": 0, "u": 0}])
    _message(refs, address="p:0",
             runs=[{"t": "First ", "b": 0, "i": 0, "u": 0}])
    refs["calc_docx_save"]()

    assert _texts(path) == ["First ", "paragraph.", "Second paragraph."]


def test_enter_twice_keeps_the_order_it_was_typed_in(opened):
    refs, state, path = opened

    _message(refs, address="new:1", after="p:0",
             runs=[{"t": "middle", "b": 0, "i": 0, "u": 0}])
    _message(refs, address="new:2", after="new:1",
             runs=[{"t": "last", "b": 0, "i": 0, "u": 0}])
    refs["calc_docx_save"]()

    assert _texts(path) == ["First paragraph.", "middle", "last",
                            "Second paragraph."]


def test_typing_in_a_paragraph_that_is_not_in_the_file_yet(opened):
    refs, state, path = opened

    _message(refs, address="new:1", after="p:0",
             runs=[{"t": "", "b": 0, "i": 0, "u": 0}])
    _message(refs, address="new:1",
             runs=[{"t": "typed after Enter", "b": 0, "i": 0, "u": 0}])
    refs["calc_docx_save"]()

    assert _texts(path)[1] == "typed after Enter"


def test_backspace_at_the_start_takes_the_paragraph_away(opened):
    refs, state, path = opened

    _message(refs, address="p:0",
             runs=[{"t": "First paragraph.Second paragraph.", "b": 0,
                    "i": 0, "u": 0}])
    _message(refs, address="p:1", drop=1)
    refs["calc_docx_save"]()

    assert _texts(path) == ["First paragraph.Second paragraph."]


def test_a_paragraph_made_and_unmade_again_never_reaches_the_file(opened):
    refs, state, path = opened

    _message(refs, address="new:1", after="p:0",
             runs=[{"t": "oops", "b": 0, "i": 0, "u": 0}])
    _message(refs, address="new:1", drop=1)

    assert state.get("docx_inserts") == []


def test_the_size_and_the_font_reach_the_file(opened):
    refs, state, path = opened

    _message(refs, address="p:0", runs=[
        {"t": "First paragraph.", "b": 0, "i": 0, "u": 0,
         "s": 18, "f": "Georgia"},
    ])
    refs["calc_docx_save"]()

    run = docx.Document(str(path)).paragraphs[0].runs[0]
    assert run.font.size.pt == 18
    assert run.font.name == "Georgia"


def test_a_font_the_view_does_not_offer_is_refused():
    with pytest.raises(docx_view.DocxError):
        docx_view.check_font("Comic Sans MS")
    assert docx_view.check_font("georgia") == "Georgia"
    assert docx_view.check_font("") == ""


def test_the_bold_button_is_not_the_paragraph_class():
    """'dw-b' is the editable paragraph block, and it carries a bottom margin.
    The Bold button wearing the same class sat above its neighbours."""
    markup = docx_view.toolbar_html()
    assert 'dw-bold' in markup
    assert 'class="dw-btn dw-b"' not in markup
    assert "'.dw-bold': 'bold'" in docx_view.edit_js('calc-scope-1')
