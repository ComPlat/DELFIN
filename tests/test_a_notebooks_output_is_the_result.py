"""notebook_read named the exception class and dropped the message.

Found by driving the tool the way a model drives it -- through
``executor.execute("notebook_read", …)`` -- on 2026-09-09. It is one of
the 32 of 72 advertised tools that no benchmark run has ever called, so
until then it had only ever been exercised by its own author's tests.

A notebook's outputs are not decoration; for a scientist they are the
result. The two questions you open a failed notebook with are "what did
this compute" and "why did that cell fail", and the tool answered:

    "output_summary": "1 output(s): error(ValueError)"
    "output_summary": "1 output(s): stream/stdout (~35 chars)"
    "output_summary": "1 output(s): execute_result(text/plain)"

-- the exception CLASS with no message, no line and no stack; a count of
characters instead of the characters; the existence of a value instead of
the value. The module's reasoning was sound and its conclusion too broad:
keeping a base64 PNG out of the prompt is worth doing, and it was applied
to text as well.

What makes it a defect rather than a taste is the tool's own description,
which tells the agent to use it INSTEAD of read_file. read_file would
have dumped the traceback -- as raw JSON, ugly, but there. The
specialised tool returned strictly less of what mattered than the generic
one it replaces.

So text comes back and bytes do not, capped per cell like source already
was, and truncated from the FRONT because that is the end an output is
informative at.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import notebook_tools as NB
from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _nb(cells: list[dict]) -> dict:
    return {"cells": cells, "metadata": {}, "nbformat": 4, "nbformat_minor": 5}


def _code(source: str, outputs: list[dict]) -> dict:
    return {"cell_type": "code", "execution_count": 1, "metadata": {},
            "outputs": outputs, "source": [source]}


@pytest.fixture
def ws(tmp_path):
    return tmp_path


def _read(ws: Path, nb: dict, **kw) -> dict:
    (ws / "n.ipynb").write_text(json.dumps(nb), encoding="utf-8")
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    return json.loads(_doc_executor.execute(
        "notebook_read", {"path": "n.ipynb", **kw}, perms))


# ---------------------------------------------------------------------------
# The three things that were dropped
# ---------------------------------------------------------------------------

def test_a_failed_cell_reports_why(ws):
    out = _read(ws, _nb([_code("E = float(row['energy'])\n", [{
        "output_type": "error", "ename": "ValueError",
        "evalue": "could not convert string to float: 'n/a'",
        "traceback": ["ValueError  Traceback (most recent call last)",
                      "Cell In[3], line 2",
                      "----> 2 E = float(row['energy'])"],
    }])]))
    text = out["cells"][0]["output"]
    assert "ValueError" in text
    assert "could not convert string to float: 'n/a'" in text   # the message
    assert "line 2" in text                                     # the place
    # ...and the summary still says what KIND of output it is.
    assert "error(ValueError)" in out["cells"][0]["output_summary"]


def test_what_a_cell_printed_comes_back(ws):
    out = _read(ws, _nb([_code("print(G)\n", [{
        "output_type": "stream", "name": "stdout",
        "text": ["G = -76.4123 Eh\n", "dG = 2.31 kcal/mol\n"]}])]))
    assert "-76.4123" in out["cells"][0]["output"]
    assert "2.31 kcal/mol" in out["cells"][0]["output"]


def test_a_result_value_comes_back(ws):
    out = _read(ws, _nb([_code("w / w.sum()\n", [{
        "output_type": "execute_result", "execution_count": 2, "metadata": {},
        "data": {"text/plain": ["array([0.61, 0.21, 0.11, 0.07])"]}}])]))
    assert "0.61" in out["cells"][0]["output"]


# ---------------------------------------------------------------------------
# ...without putting back what the summary existed to keep out
# ---------------------------------------------------------------------------

def test_a_plot_does_not_arrive_as_base64(ws):
    """The reason the summary was written. A figure's PNG is thousands of
    tokens of nothing readable, and it stays out."""
    blob = "iVBORw0KGgo" + "A" * 40_000
    out = _read(ws, _nb([_code("plt.plot(x, y)\n", [{
        "output_type": "display_data", "metadata": {},
        "data": {"image/png": blob,
                 "text/plain": ["<Figure size 640x480>"]}}])]))
    cell = out["cells"][0]
    assert blob[:40] not in json.dumps(out)
    # The readable half of the same output is kept.
    assert "<Figure size 640x480>" in cell["output"]
    assert "image/png" in cell["output_summary"]


def test_a_long_output_is_capped_from_the_front(ws):
    """A loop that printed for ten minutes ends with the answer, and a
    traceback ends with the line that failed. Source truncates in the
    middle; output must not."""
    text = "".join(f"step {i}\n" for i in range(2000)) + "FINAL RESULT 42\n"
    out = _read(ws, _nb([_code("run()\n", [{
        "output_type": "stream", "name": "stdout", "text": [text]}])]),
        max_output_chars=300)
    got = out["cells"][0]["output"]
    assert len(got) <= 400
    assert "FINAL RESULT 42" in got
    assert "step 0" not in got
    assert "omitted" in got


def test_a_cell_with_no_output_has_no_output_key(ws):
    """A key per cell saying "nothing here" is prompt spent on nothing."""
    out = _read(ws, _nb([
        {"cell_type": "markdown", "metadata": {}, "source": ["# Title\n"]},
        _code("x = 1\n", []),
    ]))
    assert "output" not in out["cells"][0]
    assert "output" not in out["cells"][1]


def test_ansi_colour_codes_are_stripped(ws):
    """Jupyter writes tracebacks with terminal colour in them; the escape
    sequences are noise in a prompt and can break a diff."""
    out = _read(ws, _nb([_code("boom()\n", [{
        "output_type": "error", "ename": "KeyError", "evalue": "'x'",
        "traceback": ["\x1b[0;31mKeyError\x1b[0m: 'x'"]}])]))
    assert "\x1b" not in out["cells"][0]["output"]
    assert "KeyError" in out["cells"][0]["output"]


# ---------------------------------------------------------------------------
# The primitive, directly
# ---------------------------------------------------------------------------

def test_read_cells_never_raises_on_a_malformed_output(tmp_path):
    p = tmp_path / "odd.ipynb"
    p.write_text(json.dumps(_nb([
        {"cell_type": "code", "metadata": {}, "source": "x=1",
         "outputs": ["not a dict", {"output_type": "stream"},
                     {"output_type": "execute_result", "data": None}]},
    ])), encoding="utf-8")
    cells = NB.read_cells(p)
    assert cells[0].output_text == ""


def test_source_truncation_is_unchanged(tmp_path):
    """The output cap is a second, separate cap -- it must not have moved
    the one that was already there."""
    p = tmp_path / "big.ipynb"
    p.write_text(json.dumps(_nb([_code("y = 1\n" * 5000, [])])),
                 encoding="utf-8")
    cells = NB.read_cells(p, max_source_chars=200)
    assert "middle of cell omitted" in cells[0].source
