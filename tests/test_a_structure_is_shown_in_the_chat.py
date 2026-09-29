"""A structure, a trajectory and an orbital are interactive in the chat.

Input: a file the agent produced in the workspace.
Output: the chat card the artifact renderer returns for it.

Before this, `.smi`, `.mol` and `.sdf` were drawn as a 2D picture and
`.xyz` -- the format DELFIN itself produces most, from smiles_to_xyz and
from every optimisation -- got a line of text: atom count and formula.
A trajectory and a cube got nothing at all, because neither suffix was
in the renderable set.

The viewer is 3Dmol, the build vendored with the package and installed
at start-up, so a login node with no outbound network draws the same
picture. The same library and the same two isosurface values the
Calculations tab uses, so a molecule does not change appearance
depending on which tab it is opened in.

These tests use written files and no binary, no network and no
installed viewer: the JS is asserted as markup, which is what the
renderer returns.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_agent as T


_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.76 0.59 0.0\nH -0.76 0.59 0.0\n"


def _xyz(tmp_path, frames: int = 1, name: str = "struct.xyz"):
    p = tmp_path / name
    p.write_text(_WATER * frames, encoding="utf-8")
    return p


# -- what counts as a frame ------------------------------------------------

def test_a_single_structure_is_one_frame(tmp_path):
    assert T._xyz_frame_count(_WATER) == 1


def test_a_trajectory_counts_its_frames(tmp_path):
    assert T._xyz_frame_count(_WATER * 12) == 12


def test_a_cut_off_last_frame_is_not_counted():
    """A killed job leaves half a frame. Counting by dividing the line
    total would report a fraction and hand the viewer a frame that is
    not there."""
    truncated = _WATER * 3 + "3\nwater\nO 0.0 0.0 0.0\n"
    assert T._xyz_frame_count(truncated) == 3


def test_something_that_is_not_an_xyz_is_no_frames():
    assert T._xyz_frame_count("Gaussian output\n  SCF Done: -76.4\n") == 0
    assert T._xyz_frame_count("") == 0


# -- which files get a viewer ---------------------------------------------

def test_a_structure_a_trajectory_and_an_orbital_are_offered(tmp_path):
    single = T._mol3d_payload(_xyz(tmp_path))
    assert single and single["fmt"] == "xyz" and single["frames"] == 1

    traj = T._mol3d_payload(_xyz(tmp_path, frames=8, name="opt.trj"))
    assert traj and traj["frames"] == 8

    cube = tmp_path / "homo.cube"
    cube.write_text("orbital\n density\n 3 0.0 0.0 0.0\n", encoding="utf-8")
    assert (T._mol3d_payload(cube) or {}).get("fmt") == "cube"


def test_a_file_the_viewer_cannot_read_gets_no_card(tmp_path):
    other = tmp_path / "notes.txt"
    other.write_text("just text", encoding="utf-8")
    assert T._mol3d_payload(other) is None

    lying = tmp_path / "output.xyz"          # the suffix says xyz, the body does not
    lying.write_text("SCF Done: -76.4\n", encoding="utf-8")
    assert T._mol3d_payload(lying) is None


def test_a_file_too_large_to_carry_in_the_transcript_is_refused(tmp_path):
    """The payload travels inside the page's HTML and the transcript
    keeps every card. A long run belongs in the Calculations tab."""
    big = tmp_path / "long.xyz"
    big.write_text("x" * (T._MOL3D_MAX_BYTES + 1), encoding="utf-8")
    assert T._mol3d_payload(big) is None

    many = _xyz(tmp_path, frames=T._MOL3D_MAX_FRAMES + 1, name="many.xyz")
    assert T._mol3d_payload(many) is None


# -- the card itself -------------------------------------------------------

def test_the_card_carries_the_coordinates_and_a_way_to_start(tmp_path):
    card = T._mol3d_card_html(_xyz(tmp_path), "struct.xyz")
    assert card is not None
    assert 'data-mol3d="' in card, "the coordinates travel in an attribute"
    assert "O 0.0 0.0 0.0" in card
    assert "__delfinMol3D" in card, "nothing would start the viewer"
    assert 'class="delfin-mol3d-stage"' in card


def test_the_card_keeps_the_text_line_under_the_viewer(tmp_path):
    """3Dmol may be absent and the element may never become visible. A
    grey box says less than the formula did before this existed."""
    card = T._mol3d_card_html(_xyz(tmp_path), "struct.xyz")
    assert "3 atoms" in card and "H2 O" in card


def test_a_trajectory_says_so_and_names_its_frames(tmp_path):
    card = T._mol3d_card_html(_xyz(tmp_path, frames=8), "opt.xyz")
    assert 'data-mol3d-frames="8"' in card
    assert "8 frames" in card
    assert "trajectory" in card


def test_the_data_cannot_break_out_of_its_attribute(tmp_path):
    """The comment line of an XYZ is free text and reaches the browser.
    Unescaped it would end the attribute and the rest would be parsed as
    markup -- with an event handler in it if somebody wanted one."""
    hostile = tmp_path / "evil.xyz"
    hostile.write_text(
        '1\n" onmouseover="alert(1)" x="<img src=x onerror=alert(2)>\n'
        "O 0.0 0.0 0.0\n", encoding="utf-8")
    card = T._mol3d_card_html(hostile, "evil.xyz")
    assert card is not None
    assert 'onmouseover="alert(1)"' not in card
    assert "<img src=x onerror=alert(2)>" not in card
    assert "&lt;img" in card or "&amp;lt;img" in card


def test_two_cards_do_not_share_a_stage(tmp_path):
    a = T._mol3d_card_html(_xyz(tmp_path, name="a.xyz"), "a.xyz")
    b = T._mol3d_card_html(_xyz(tmp_path, name="b.xyz"), "b.xyz")
    import re
    ids = [re.search(r'id="(delfin_chat_mol3d_\d+)"', c).group(1) for c in (a, b)]
    assert ids[0] != ids[1], "two viewers in one stage draw over each other"


def test_the_renderer_returns_the_card_for_every_read_format(tmp_path):
    """The whole chain, as the chat calls it."""
    for name, frames in (("struct.xyz", 1), ("opt.trj", 5)):
        body = T._render_artifact_body(_xyz(tmp_path, frames=frames, name=name))
        assert body and "delfin-chat-mol3d" in body, name


def test_a_viewer_switched_off_in_the_settings_is_not_drawn(tmp_path, monkeypatch):
    """The setting is a user's decision about every viewer in DELFIN;
    the chat is not an exception to it."""
    from delfin.dashboard import molecule_viewer

    monkeypatch.setattr(molecule_viewer, "get_viewer_profile",
                        lambda: {"enabled": False})
    assert T._mol3d_card_html(_xyz(tmp_path), "struct.xyz") is None

    # ... and the text card still answers what the molecule is.
    body = T._render_artifact_body(_xyz(tmp_path))
    assert body and "H2 O" in body and "delfin-chat-mol3d" not in body


def test_the_page_script_is_installed_once_and_starts_the_viewer():
    """The card calls window.__delfinMol3D; something has to define it."""
    import pathlib
    src = pathlib.Path(T.__file__).read_text(encoding="utf-8")
    assert "window.__delfinMol3D = function" in src
    assert "__delfinMol3DInstalled" in src, "the guard against a second install"
    assert "addModelsAsFrames" in src, "a trajectory would not animate"
    assert "addVolumetricData" in src, "an orbital would have no isosurface"
