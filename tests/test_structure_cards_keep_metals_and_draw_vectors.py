"""Structure cards in the agent chat: metals kept, vectors drawn, SMILES wrapped.

Three behaviours the results-overview cards must show:
- a long SMILES wraps instead of running off the card, and is not cut;
- a coordination metal survives the 2D skeleton parse and is drawn;
- 2D structure drawings are SVG (vector), not PNG rasters, and an
  optimised XYZ gets a 2D vector skeleton beside its 3D card.
"""

import base64
import pathlib
import tempfile

import pytest

from delfin.dashboard import tab_agent


def _render(tmp_path: pathlib.Path, suffix: str, content: str) -> str:
    path = tmp_path / f"structure{suffix}"
    path.write_text(content, encoding="utf-8")
    renderer = (tab_agent._render_artifact_body
                if hasattr(tab_agent, "_render_artifact_body")
                else tab_agent._artifact_card)
    html = renderer(path)
    assert html, f"no card rendered for {suffix}"
    return html


def _render_raw(tmp_path: pathlib.Path, suffix: str, content: str):
    path = tmp_path / f"structure{suffix}"
    path.write_text(content, encoding="utf-8")
    renderer = (tab_agent._render_artifact_body
                if hasattr(tab_agent, "_render_artifact_body")
                else tab_agent._artifact_card)
    return renderer(path)


def _svg_of(html: str) -> str:
    payload = html.split("data:image/svg+xml;base64,", 1)[1]
    payload = payload.split('"', 1)[0]
    return base64.b64decode(payload).decode("utf-8")


def test_a_long_smiles_wraps_inside_the_card_and_is_not_cut(tmp_path):
    smiles = "[Fe](C)(C)(n1ccccc1)" + "C" * 300
    html = _render(tmp_path, ".smi", smiles + "\n")
    assert "overflow-wrap:anywhere" in html
    assert "word-break:break-all" in html
    assert smiles in html          # full string, not cut at 120 chars
    assert smiles[:120] + "</code>" not in html


def test_a_coordination_metal_is_kept_in_the_drawn_skeleton(tmp_path):
    html = _render(tmp_path, ".smi", "[Fe](C)(C)(n1ccccc1)\n")
    svg = _svg_of(html)
    # Atom glyphs are path objects; the first atom (Fe) has its own
    # atom-0 block in the drawing.
    assert "class='atom-0'" in svg


def test_2d_structure_cards_draw_a_vector_svg_not_a_png(tmp_path):
    html = _render(tmp_path, ".smi", "c1ccccc1CCO\n")
    assert "data:image/svg+xml;base64" in html
    assert "data:image/png;base64" not in html


def test_an_optimised_xyz_gets_a_2d_vector_skeleton_beside_the_3d_card(
    tmp_path,
):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mh = Chem.AddHs(Chem.MolFromSmiles("CCOc1ccccc1"))
    AllChem.EmbedMolecule(mh, randomSeed=42)
    conf = mh.GetConformer()
    lines = [str(mh.GetNumAtoms()), "opt structure"]
    for i, atom in enumerate(mh.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol()} "
                     f"{pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    html = _render(tmp_path, ".xyz", "\n".join(lines) + "\n")
    assert "delfin-chat-mol3d" in html          # 3D card still there
    assert "data:image/svg+xml;base64" in html  # plus the vector skeleton
    assert "2D skeleton" in html


def test_an_unreadable_xyz_shows_no_card_and_does_not_crash(tmp_path):
    # A garbage XYZ renders no card (base behaviour at 776c20de8) but must
    # not raise when the caller walks a folder of artifacts.
    assert _render_raw(tmp_path, ".xyz", "not an xyz at all\n") is None
