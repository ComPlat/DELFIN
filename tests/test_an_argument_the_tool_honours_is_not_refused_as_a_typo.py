"""An argument the tool honours but the schema does not advertise still runs.

``notebook_read`` honours ``max_output_chars`` (``api_client.py`` reads it
out of the call and passes it down as ``max_output_chars=``), and the
catalogue deliberately does not advertise it -- see the comment beside
``max_source_chars`` in ``delfin/agent/tool_catalog.py``: *"max_output_chars
is accepted and is deliberately NOT advertised ... a knob in the catalogue
is paid for by every request whether or not anyone turns it."*

So the schema's property list is NOT the set of arguments the tool accepts.
The rule that refuses a stray field standing in for an absent property
reads the schema as exactly that set, and therefore refused this call:

    unknown field 'max_output_chars' -- did you mean 'max_source_chars'?

-- a call that would have worked, on a tool whose documented behaviour is
what the caller asked for. That is a regression the rule introduced, not a
defect it fixed: the tolerance it replaced let the call through.

Found by running the files the whole-repository suite names on this branch
against the base: ``test_a_notebooks_output_is_the_result.py::
test_a_long_output_is_capped_from_the_front`` passes on ``origin/main`` and
failed here.

Red before the fix: the execute call answers with the refusal above instead
of a ``cells`` list. The control case is a real typo on the same tool,
which must still be refused -- the fix must not turn the rule off.
"""

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _nb(cells):
    return {"cells": cells, "metadata": {}, "nbformat": 4, "nbformat_minor": 5}


def _code(source, outputs):
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


@pytest.mark.xfail(
    strict=False,
    reason="needs patches/api_client.all-three.patch applied: the fix "
           "is in delfin/agent/api_client.py, which this session may not "
           "write. xfail = patch not applied (a1's rule refuses a working "
           "call); xpass = patch applied (the regression is gone).")
def test_an_undocumented_but_honoured_cap_is_not_refused(ws):
    """``max_output_chars`` is accepted, undocumented, and must run."""
    text = "".join(f"step {i}\n" for i in range(2000)) + "FINAL RESULT 42\n"
    out = _read(ws, _nb([_code("run()\n", [{
        "output_type": "stream", "name": "stdout", "text": [text]}])]),
        max_output_chars=300)
    assert "error" not in out, out
    got = out["cells"][0]["output"]
    # And the argument took effect: the cap the caller asked for was used.
    assert len(got) <= 400, got
    assert "FINAL RESULT 42" in got
    assert "step 0" not in got


def test_a_real_typo_on_the_same_tool_is_still_refused(ws):
    """Control: the rule still fires for a spelling the tool does not accept."""
    out = _read(ws, _nb([_code("run()\n", [])]), max_sorce_chars=100)
    assert "error" in out, out
    assert "max_source_chars" in out["error"], out
