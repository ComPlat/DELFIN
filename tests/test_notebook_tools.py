"""Cell-level notebook reading and editing.

``delfin/agent/notebook_tools.py`` is a write path over the user's own
files -- ``apply_edit`` rewrites a whole ``.ipynb`` in place -- and no
test file imported it. What existed went through the tool dispatcher and
asserted on the JSON the tool returns, so the module's own contract (what
survives an edit, what is rejected, what a truncated cell looks like) was
never stated anywhere.

Path containment is deliberately NOT this module's job: the docstring
says the caller resolves the path and runs the sandbox check first. The
last test here pins that division so the guarantee is not assumed to live
in two places and end up in neither.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import notebook_tools as NT


def _nb(*cells) -> dict:
    return {"cells": list(cells), "metadata": {}, "nbformat": 4,
            "nbformat_minor": 5}


def _code(source, *, cell_id=None, outputs=None, execution_count=None):
    cell = {"cell_type": "code", "metadata": {}, "source": source,
            "outputs": outputs or [], "execution_count": execution_count}
    if cell_id:
        cell["id"] = cell_id
    return cell


def _write(tmp_path, nb, name="analyse.ipynb"):
    p = tmp_path / name
    p.write_text(json.dumps(nb, indent=1), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def test_a_source_written_as_lines_reads_as_one_string(tmp_path):
    """nbformat allows either; the agent must not have to know which."""
    p = _write(tmp_path, _nb(_code(["import pandas as pd\n", "pd\n"]),
                             _code("x = 1\n")))
    cells = NT.read_cells(p)
    assert [c.source for c in cells] == ["import pandas as pd\npd\n", "x = 1\n"]
    assert [c.idx for c in cells] == [0, 1]
    assert [c.cell_type for c in cells] == ["code", "code"]


def test_a_markdown_cell_keeps_its_type(tmp_path):
    p = _write(tmp_path, _nb({"cell_type": "markdown", "metadata": {},
                              "source": "# Titel\n"}))
    assert [c.cell_type for c in NT.read_cells(p)] == ["markdown"]


def test_an_image_output_is_summarised_not_pasted(tmp_path):
    """A base64 plot in the tool result costs tokens and tells the agent
    nothing it can act on."""
    blob = "iVBORw0KGgo" + "A" * 5000
    p = _write(tmp_path, _nb(_code("plot()\n", outputs=[
        {"output_type": "display_data", "data": {"image/png": blob}}])))
    summary = NT.read_cells(p)[0].output_summary
    assert blob not in summary
    assert "display_data" in summary and "image/png" in summary


def test_a_stream_output_reports_how_much_text_there_was(tmp_path):
    p = _write(tmp_path, _nb(_code("print(x)\n", outputs=[
        {"output_type": "stream", "name": "stdout",
         "text": ["a" * 120, "b" * 80]}])))
    summary = NT.read_cells(p)[0].output_summary
    assert "stream/stdout" in summary
    assert "200 chars" in summary


def test_an_error_output_names_the_exception(tmp_path):
    """The one output the agent has to react to."""
    p = _write(tmp_path, _nb(_code("1/0\n", outputs=[
        {"output_type": "error", "ename": "ZeroDivisionError",
         "evalue": "division by zero", "traceback": ["..."]}])))
    assert "error(ZeroDivisionError)" in NT.read_cells(p)[0].output_summary


def test_a_cell_without_output_says_nothing_rather_than_zero(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    assert NT.read_cells(p)[0].output_summary == ""


def test_a_huge_cell_is_cut_in_the_middle_and_says_so(tmp_path):
    """Head and tail, because the imports are at the top and the failing
    line is usually at the bottom."""
    source = "HEAD\n" + ("f\n" * 4000) + "TAIL"
    p = _write(tmp_path, _nb(_code(source)))
    out = NT.read_cells(p, max_source_chars=200)[0].source
    assert out.startswith("HEAD")
    assert out.endswith("TAIL")
    assert "chars from middle of cell omitted" in out
    assert len(out) < len(source)


def test_a_short_cell_is_not_touched(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    assert NT.read_cells(p, max_source_chars=200)[0].source == "x = 1\n"


def test_a_notebook_with_no_cells_reads_as_no_cells(tmp_path):
    assert NT.read_cells(_write(tmp_path, _nb())) == []


def test_a_malformed_cell_is_skipped_not_crashed_on(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n"), "not a cell",
                             _code("y = 2\n")))
    cells = NT.read_cells(p)
    assert [c.source for c in cells] == ["x = 1\n", "y = 2\n"]


# ---------------------------------------------------------------------------
# Editing
# ---------------------------------------------------------------------------

def test_replacing_a_cell_leaves_its_neighbours_alone(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n"), _code("y = 2\n"),
                             _code("z = 3\n")))
    before, after = NT.apply_edit(p, 1, "replace", source="y = 99\n")
    assert (before, after) == (3, 3)
    assert [c.source for c in NT.read_cells(p)] == [
        "x = 1\n", "y = 99\n", "z = 3\n"]


def test_a_replaced_cell_keeps_its_identity_and_loses_its_stale_output(
        tmp_path):
    """The id is what a notebook's own history and links are keyed on; the
    output belongs to the code that is gone."""
    p = _write(tmp_path, _nb(_code(
        "x = 1\n", cell_id="abc123", execution_count=7,
        outputs=[{"output_type": "stream", "name": "stdout",
                  "text": ["1\n"]}])))
    NT.apply_edit(p, 0, "replace", source="x = 2\n")
    cell = json.loads(p.read_text(encoding="utf-8"))["cells"][0]
    assert cell["id"] == "abc123"
    assert cell["outputs"] == []
    assert cell["execution_count"] is None


def test_a_replaced_cell_keeps_its_metadata(tmp_path):
    p = _write(tmp_path, _nb({"cell_type": "code", "source": "x = 1\n",
                              "outputs": [], "execution_count": None,
                              "metadata": {"tags": ["parameters"]}}))
    NT.apply_edit(p, 0, "replace", source="x = 2\n")
    cell = json.loads(p.read_text(encoding="utf-8"))["cells"][0]
    assert cell["metadata"] == {"tags": ["parameters"]}


@pytest.mark.parametrize("mode,expected", [
    ("insert_before", ["neu\n", "x = 1\n", "y = 2\n"]),
    ("insert_after", ["x = 1\n", "neu\n", "y = 2\n"]),
])
def test_an_insert_lands_on_the_side_it_names(tmp_path, mode, expected):
    p = _write(tmp_path, _nb(_code("x = 1\n"), _code("y = 2\n")))
    before, after = NT.apply_edit(p, 0, mode, source="neu\n")
    assert (before, after) == (2, 3)
    assert [c.source for c in NT.read_cells(p)] == expected


def test_an_inserted_cell_gets_an_id_of_its_own(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n", cell_id="abc123")))
    NT.apply_edit(p, 0, "insert_after", source="y = 2\n")
    ids = [c.get("id") for c in
           json.loads(p.read_text(encoding="utf-8"))["cells"]]
    assert ids[0] == "abc123"
    assert ids[1] and ids[1] != "abc123"


def test_an_inserted_markdown_cell_carries_no_code_fields(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    NT.apply_edit(p, 0, "insert_before", source="# Titel\n",
                  cell_type="markdown")
    cell = json.loads(p.read_text(encoding="utf-8"))["cells"][0]
    assert cell["cell_type"] == "markdown"
    assert "outputs" not in cell and "execution_count" not in cell


def test_appending_at_the_end_is_allowed(tmp_path):
    """``insert_before`` at ``n`` is how a cell is appended; the bound for
    an insert is one past the last cell, unlike replace and delete."""
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    before, after = NT.apply_edit(p, 1, "insert_before", source="y = 2\n")
    assert (before, after) == (1, 2)
    assert [c.source for c in NT.read_cells(p)] == ["x = 1\n", "y = 2\n"]


def test_deleting_a_cell_removes_exactly_one(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n"), _code("y = 2\n"),
                             _code("z = 3\n")))
    before, after = NT.apply_edit(p, 1, "delete")
    assert (before, after) == (3, 2)
    assert [c.source for c in NT.read_cells(p)] == ["x = 1\n", "z = 3\n"]


def test_a_delete_needs_no_source(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    assert NT.apply_edit(p, 0, "delete") == (1, 0)


@pytest.mark.parametrize("mode,idx", [
    ("replace", 2), ("replace", -1),
    ("delete", 2), ("delete", -1),
    ("insert_before", 3), ("insert_after", -1),
])
def test_an_index_outside_the_notebook_is_refused(tmp_path, mode, idx):
    p = _write(tmp_path, _nb(_code("x = 1\n"), _code("y = 2\n")))
    original = p.read_text(encoding="utf-8")
    with pytest.raises(ValueError, match="out of range"):
        NT.apply_edit(p, idx, mode, source="neu\n")
    assert p.read_text(encoding="utf-8") == original, (
        "the notebook was rewritten by a call that was rejected")


def test_an_unknown_mode_is_refused(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    with pytest.raises(ValueError, match="mode must be one of"):
        NT.apply_edit(p, 0, "append", source="neu\n")


def test_an_unknown_cell_type_is_refused(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    with pytest.raises(ValueError, match="cell_type must be one of"):
        NT.apply_edit(p, 0, "replace", source="neu\n", cell_type="python")


def test_an_edit_without_a_source_is_refused(tmp_path):
    """Silently writing an empty cell would delete the code the model
    meant to keep."""
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    with pytest.raises(ValueError, match="source is required"):
        NT.apply_edit(p, 0, "replace")


def test_the_rest_of_the_notebook_survives_an_edit(tmp_path):
    """Kernel and language info live outside ``cells`` and are what makes
    the file openable."""
    nb = _nb(_code("x = 1\n"))
    nb["metadata"] = {"kernelspec": {"name": "python3",
                                     "display_name": "Python 3"},
                      "language_info": {"name": "python", "version": "3.11"}}
    p = _write(tmp_path, nb)
    NT.apply_edit(p, 0, "replace", source="x = 2\n")
    after = json.loads(p.read_text(encoding="utf-8"))
    assert after["metadata"] == nb["metadata"]
    assert after["nbformat"] == 4 and after["nbformat_minor"] == 5


def test_the_file_stays_valid_json_a_notebook_reader_accepts(tmp_path):
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    NT.apply_edit(p, 0, "insert_after", source="print('hallo')\n")
    text = p.read_text(encoding="utf-8")
    assert text.endswith("\n")
    nb = json.loads(text)
    assert [c["cell_type"] for c in nb["cells"]] == ["code", "code"]


def test_a_written_source_is_a_list_of_lines(tmp_path):
    """Jupyter writes list-of-lines; emitting one long string churns the
    diff of every cell the agent touches."""
    p = _write(tmp_path, _nb(_code("x = 1\n")))
    NT.apply_edit(p, 0, "replace", source="a = 1\nb = 2\n")
    assert json.loads(p.read_text(encoding="utf-8"))["cells"][0]["source"] == [
        "a = 1\n", "b = 2\n"]


# ---------------------------------------------------------------------------
# Where containment lives
# ---------------------------------------------------------------------------

def test_the_module_itself_does_not_police_the_path(tmp_path):
    """Stated, not assumed: this module takes an ALREADY-RESOLVED path and
    the tool layer does the sandbox check. Writing that down here is what
    keeps the guarantee from being expected in two places and implemented
    in neither -- the next test is the half that does enforce it."""
    outside = tmp_path / "outside.ipynb"
    outside.write_text(json.dumps(_nb(_code("x = 1\n"))), encoding="utf-8")
    NT.apply_edit(outside, 0, "replace", source="x = 2\n")
    assert [c.source for c in NT.read_cells(outside)] == ["x = 2\n"]


def test_the_tool_layer_refuses_a_notebook_outside_the_workspace(tmp_path):
    from delfin.agent import api_client as A

    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "outside.ipynb"
    outside.write_text(json.dumps(_nb(_code("x = 1\n"))), encoding="utf-8")

    perms = A.KitToolPermissions(workspace=ws, mode="acceptEdits")
    perms.confirm_callback = lambda n, a, prev: True
    out = A._doc_executor.execute(
        "notebook_edit",
        {"path": str(outside), "cell_idx": 0, "mode": "replace",
         "source": "x = 99\n"},
        perms)

    assert "error" in out, out
    assert [c.source for c in NT.read_cells(outside)] == ["x = 1\n"], (
        "a notebook outside the workspace was rewritten")
