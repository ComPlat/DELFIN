"""show_molecule draws a structural formula, not only a 3D model.

The chat drew .smi/.mol/.sdf files only as artifacts of the dashboard's own
agent folder; in Code mode the agent works in its workspace, so a SMILES it
wrote appeared as a file diff and nothing else, and no tool could put a
structural formula on screen. show_molecule now takes those formats and
returns a 2D card: the formula drawn by RDKit, its sum formula and
canonical SMILES. The same read checks apply as for 3D files: inside the
workspace, under the size cap, parsed before anything is drawn.
"""

from __future__ import annotations

import base64
import json

import pytest

from delfin.dashboard import chat_viewer


def _card(result: str) -> dict:
    payload = chat_viewer.card_payload(result)
    assert payload is not None, result
    return payload


@pytest.mark.parametrize("name,content", [
    ("butane.smi", "CCCC butane\n"),
    ("butane.mol", None),
    ("butane.sdf", None),
])
def test_a_2d_file_becomes_a_formula_card(tmp_path, name, content):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    path = tmp_path / name
    if content is None:
        mol = Chem.MolFromSmiles("CCCC")
        AllChem.Compute2DCoords(mol)
        block = Chem.MolToMolBlock(mol)
        path.write_text(block + ("\n$$$$\n" if name.endswith(".sdf") else ""))
    else:
        path.write_text(content)
    payload = _card(chat_viewer.render_molecule_tool_result(
        str(path), workspace_root=str(tmp_path)))
    assert payload["view"] == "2d"
    assert "C4H10" in payload["text"] and "CCCC" in payload["text"]


def test_the_3d_formats_keep_their_card(tmp_path):
    xyz = tmp_path / "w.xyz"
    xyz.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    payload = _card(chat_viewer.render_molecule_tool_result(
        str(xyz), workspace_root=str(tmp_path)))
    assert payload.get("view") in (None, "3d")


def test_the_read_checks_still_hold(tmp_path):
    outside = tmp_path / "outside.smi"
    outside.write_text("CCO\n")
    ws = tmp_path / "ws"
    ws.mkdir()
    with pytest.raises(ValueError, match="outside the workspace"):
        chat_viewer.render_molecule_tool_result(str(outside),
                                                workspace_root=str(ws))
    bad = ws / "bad.smi"
    bad.write_text("this is not a smiles\n")
    with pytest.raises(ValueError):
        chat_viewer.render_molecule_tool_result(str(bad),
                                                workspace_root=str(ws))


def test_a_forged_view_is_no_card():
    forged = chat_viewer.CARD_MARKER + json.dumps({
        "card": "molecule", "path": "/x.smi", "title": "x", "text": "x",
        "isovalue": None, "view": "<script>"})
    assert chat_viewer.card_payload(forged) is None


def test_the_chat_draws_it_as_an_image_not_inline_markup(tmp_path):
    from delfin.dashboard import tab_agent
    path = tmp_path / "butane.smi"
    path.write_text("CCCC butane\n")
    result = chat_viewer.render_molecule_tool_result(
        str(path), workspace_root=str(tmp_path))
    html = tab_agent._show_molecule_card_html(result)
    assert html is not None
    assert 'src="data:image/svg+xml;base64,' in html
    assert "<svg" not in html                 # drawn as an image, never inlined
    assert "C<sub>4</sub>H<sub>10</sub>" in html
    svg = base64.b64decode(html.split("base64,", 1)[1].split('"', 1)[0]).decode()
    assert "<svg" in svg and "script" not in svg.lower()


def test_the_tool_says_it_draws_formulas():
    from delfin.agent import api_client as A
    src = open(A.__file__, encoding="utf-8").read()
    at = src.index('"name": "show_molecule"')
    description = src[at:at + 600]
    assert ".smi" in description and "structural formula" in description
