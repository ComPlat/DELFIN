"""The Review tab for a long blinded rating: many categories, an unsure mark, a window-sized view.

Synthetic data only: invented file and category names, a water molecule
moved a little per frame.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("ipywidgets")

from delfin import review as rv  # noqa: E402
from delfin.dashboard import tab_review  # noqa: E402
from delfin.dashboard.tab_review import ReviewPanel  # noqa: E402

CATS20 = [f"C{k:02d}" for k in range(1, 21)]


def _xyz(stem: str, n_frames: int) -> str:
    return "".join(
        f"3\n{stem} frame{k} label-{k}\n"
        f"O 0.000 0.000 {0.117 + 0.01 * k:.3f}\n"
        f"H 0.000 0.757 {-0.467 + 0.01 * k:.3f}\n"
        f"H 0.000 -0.757 {-0.467 + 0.01 * k:.3f}\n"
        for k in range(n_frames))


@pytest.fixture
def review_dir(tmp_path: Path) -> Path:
    d = tmp_path / "set" / "structures"
    d.mkdir(parents=True)
    for k in range(1, 5):
        (d / f"T{k:04d}.xyz").write_text(_xyz(f"T{k:04d}", 1))
    return d


def _panel(where: Path, cats=None, **kw) -> ReviewPanel:
    panel = ReviewPanel(where)
    if cats is not None:
        panel.categories_input.value = "\n".join(cats)
    for name, value in kw.items():
        getattr(panel, name).value = value
    assert panel.open(str(where)) is not None, panel.status_html.value
    return panel


# ---------------------------------------------------------------- keys 0 and Shift+digit
def test_key_zero_is_the_tenth_and_shift_digits_the_eleventh_to_nineteenth(review_dir):
    panel = _panel(review_dir, CATS20)
    labels = [b.description for b in panel.category_buttons]
    assert labels[0].startswith("1  ") and labels[8].startswith("9  ")
    assert labels[9].startswith("0  ")
    assert labels[10].startswith("⇧1  ") and labels[18].startswith("⇧9  ")
    assert labels[19] == "C20"                                # mouse only, no key
    for key in ("0", "s1", "s9"):
        panel.handle_key(key)
    assert panel.selected_categories() == ["C10", "C11", "C19"]
    panel.handle_key("s1")                                    # toggles off again
    panel.handle_key("b")
    rec = panel.session.record("T0001.xyz", 0)
    assert rec["categories"] == ["C10", "C19"]


def test_category_keys_map_both_ways():
    for i in range(19):
        assert tab_review.category_index(tab_review.category_key(i)) == i
    assert tab_review.category_key(19) == ""
    assert tab_review.category_index("s0") is None
    assert tab_review.category_index("x") is None


def test_a_category_key_for_a_frame_no_longer_shown_is_dropped(review_dir):
    panel = _panel(review_dir, CATS20)
    shown = panel.session.position
    panel.handle_key("p")
    panel.keys_input.value = f"5 0 {shown}"
    panel.keys_input.value = f"6 s2 {shown}"
    assert panel.selected_categories() == []
    panel.keys_input.value = f"7 s2 {panel.session.position}"
    assert panel.selected_categories() == ["C12"]


def test_the_key_script_reads_digits_from_their_place_and_sends_shift_as_s():
    js = tab_review._keys_js()
    assert "ev.code" in js and "Digit|Numpad" in js
    assert "ev.shiftKey" in js and "'s'+dm[2]" in js
    for key in ("u:'u'", "h:'h'", "r:'r'"):
        assert key in js


# ---------------------------------------------------------------- appending categories
def test_appending_categories_keeps_the_ratings_and_the_keys(review_dir):
    panel = _panel(review_dir, ["A", "B", "C"])
    panel.handle_key("2")
    panel.handle_key("b")
    rated = dict(panel.session.record("T0001.xyz", 0))
    panel.categories_input.value = "C\nA\nNEW"                 # reordered, B left out
    panel.apply_categories()
    assert panel.session.categories == ["A", "B", "C", "NEW"]
    assert [b._review_category for b in panel.category_buttons] == ["A", "B", "C", "NEW"]
    assert panel.category_buttons[1].description.startswith("2  ")
    assert panel.category_buttons[3].description.startswith("4  ")
    added = panel.add_categories("X, Y,NEW")
    assert added == ["X", "Y"]
    on_disk = json.loads(panel.session.review_path.read_text())
    assert on_disk["categories"] == ["A", "B", "C", "NEW", "X", "Y"]
    assert panel.session.record("T0001.xyz", 0) == rated
    panel.handle_key("left")
    assert panel.selected_categories() == ["B"]


def test_a_choice_on_the_shown_frame_survives_adding_a_category(review_dir):
    panel = _panel(review_dir, ["A", "B"])
    panel.handle_key("p")                                     # now a review with ratings
    panel.handle_key("2")
    panel.add_categories(["Z"])
    assert panel.selected_categories() == ["B"]
    panel.handle_key("3")
    panel.handle_key("b")
    assert panel.session.record("T0002.xyz", 0)["categories"] == ["B", "Z"]


def test_before_any_rating_the_field_replaces_the_list(review_dir):
    panel = _panel(review_dir, ["A", "B"])
    panel.categories_input.value = "Q\nR"
    panel.apply_categories()
    assert panel.session.categories == ["Q", "R"]


# ---------------------------------------------------------------- kategorien.txt
def test_a_kategorien_file_one_folder_up_fills_a_new_review(review_dir):
    (review_dir.parent / rv.CATEGORY_FILE_NAME).write_text("# one per line\nTOPO\n\nCLASH\nH\n")
    panel = ReviewPanel(review_dir)                           # start dir is the folder
    assert panel.categories_input.value.splitlines() == ["TOPO", "CLASH", "H"]
    panel.open(str(review_dir))
    assert panel.session.categories == ["TOPO", "CLASH", "H"]
    assert [b._review_category for b in panel.category_buttons] == ["TOPO", "CLASH", "H"]


def test_the_kategorien_file_follows_the_path_and_never_replaces_a_typed_list(review_dir,
                                                                               tmp_path):
    (review_dir / rv.CATEGORY_FILE_NAME).write_text("K1\nK2\n")
    panel = ReviewPanel(tmp_path)
    assert panel.categories_input.value.splitlines() == list(rv.DEFAULT_CATEGORIES)
    panel.path_input.value = str(review_dir)
    assert panel.categories_input.value.splitlines() == ["K1", "K2"]
    panel.path_input.value = str(tmp_path)                    # an earlier fill is not "typed"
    assert panel.categories_input.value.splitlines() == list(rv.DEFAULT_CATEGORIES)
    panel.categories_input.value = "MINE"
    panel.open(str(review_dir))
    assert panel.session.categories == ["MINE"]


def test_a_resumed_review_keeps_its_own_categories_over_a_kategorien_file(review_dir):
    first = _panel(review_dir, ["OWN1", "OWN2"])
    first.handle_key("p")
    (review_dir / rv.CATEGORY_FILE_NAME).write_text("K1\nK2\n")
    again = ReviewPanel(review_dir)
    again.open(str(review_dir))
    assert again.session.categories == ["OWN1", "OWN2"]
    assert again.session.progress()["rated"] == 1


# ---------------------------------------------------------------- the unsure mark
def test_u_adds_and_removes_unsicher_in_the_comment(review_dir):
    panel = _panel(review_dir)
    panel.handle_key("u")
    assert panel.note_input.value == "unsicher" and panel.unsure_btn.value
    on_disk = json.loads(panel.session.review_path.read_text())
    assert on_disk["notes"]["T0001.xyz"]["0"]["note"] == "unsicher"   # kept before a rating
    panel.handle_key("b")
    assert panel.session.record("T0001.xyz", 0)["note"] == "unsicher"
    assert not panel.unsure_btn.value                         # next frame, no mark
    panel.handle_key("left")
    assert panel.unsure_btn.value
    panel.handle_key("u")
    assert panel.note_input.value == ""
    assert panel.session.record("T0001.xyz", 0)["note"] == ""
    assert json.loads(panel.session.review_path.read_text())[
        "reviews"]["T0001.xyz"]["0"]["note"] == ""


def test_the_unsure_mark_keeps_the_rest_of_the_comment_and_the_button_does_the_same(review_dir):
    panel = _panel(review_dir)
    panel.note_input.value = "ring odd"
    panel.unsure_btn.value = True
    assert panel.note_input.value == "ring odd unsicher"
    panel.handle_key("u")
    assert panel.note_input.value == "ring odd"
    panel.note_input.value = "unsure about the ring"          # typed by hand
    assert panel.unsure_btn.value
    shown = panel.session.position
    panel.handle_key("p")
    panel.keys_input.value = f"9 u {shown}"                   # stale: dropped
    assert panel.note_input.value == ""


def test_the_unsure_word_helpers():
    assert rv.with_unsure("", True) == "unsicher"
    assert rv.with_unsure("x unsicher", True) == "x unsicher"
    assert rv.with_unsure("x, unsicher", False) == "x"
    assert rv.with_unsure("ring unsure", False) == "ring"
    assert rv.is_unsure("UNSICHER") and rv.is_unsure("bit unsure")
    assert not rv.is_unsure("Unsicherheit") and not rv.is_unsure("")


# ---------------------------------------------------------------- the window-sized view
def test_the_viewer_height_follows_the_window_through_the_shared_mechanism(review_dir,
                                                                          monkeypatch):
    """No fixed pixel height: a viewport length (as the tab bodies use) on the
    viewer box, and the dashboard's viewer factory resizes the canvas with it."""
    from delfin.dashboard import molecule_viewer
    assert "100vh" in tab_review.VIEWER_HEIGHT
    seen = {}
    monkeypatch.setattr(molecule_viewer, "render_xyz_in_output",
                        lambda out, xyz, width=None, height=None:
                        seen.update(width=width, height=height))
    panel = _panel(review_dir)
    assert seen == {"width": "100%", "height": tab_review.VIEWER_HEIGHT}
    assert panel.sidebar.layout.height == tab_review.VIEWER_HEIGHT
    assert not panel.rating_box.layout.max_height       # no inner scroll bar
    view = molecule_viewer.py3Dmol.view(width=tab_review.VIEWER_WIDTH,
                                        height=tab_review.VIEWER_HEIGHT)
    molecule_viewer.apply_molecule_view_style(view)
    assert f"height: {tab_review.VIEWER_HEIGHT};" in view.startjs
    assert "window.__delfinCreateViewer(" in view.startjs
    assert "window.__delfinKeepSized(viewer, element)" in \
        molecule_viewer.RIGHT_MOUSE_TRANSLATE_PATCH_JS


def test_the_list_is_folded_away_in_the_blinded_mode_and_the_setup_closes(review_dir):
    panel = _panel(review_dir, blind_cb=True)
    assert panel.sidebar.layout.display == "none"
    assert panel.setup.selected_index is None
    assert panel.more.selected_index is None
    assert "0 / 4 rated" in panel.progress_html.value
    panel.handle_key("p")
    panel.handle_key("b")
    assert "2 / 4 rated" in panel.progress_html.value
    assert "1</b> pass" in panel.progress_html.value
    assert "1</b> block" in panel.progress_html.value


def test_h_atoms_can_be_hidden_and_r_draws_the_view_again(review_dir, monkeypatch):
    from delfin.dashboard import molecule_viewer
    drawn = []
    monkeypatch.setattr(molecule_viewer, "render_xyz_in_output",
                        lambda out, xyz, width=None, height=None: drawn.append(xyz))
    panel = _panel(review_dir)
    assert drawn[-1].splitlines()[0] == "3"                   # all atoms by default
    panel.handle_key("h")
    assert drawn[-1].splitlines()[0] == "1"
    assert all(not ln.startswith("H ") for ln in drawn[-1].splitlines()[2:])
    panel.handle_key("right")
    assert drawn[-1].splitlines()[0] == "1"                   # stays hidden on the next frame
    n = len(drawn)
    panel.handle_key("r")
    assert len(drawn) == n + 1
    panel.handle_key("h")
    assert drawn[-1].splitlines()[0] == "3"


# ---------------------------------------------------------------- one list, a go-to field
def test_a_folder_of_one_frame_files_shows_one_list_with_the_verdicts(review_dir):
    """Four one-frame files: the frame list (always one line) is folded away
    and each file carries its own verdict mark."""
    panel = _panel(review_dir)
    assert panel.frames_box.layout.display == "none"
    assert panel.files_box.layout.display == ""
    panel.handle_key("p")
    labels = [label for label, _ in panel.file_list.options]
    assert labels[0] == "✓ T0001.xyz" and labels[1] == "· T0002.xyz"


def test_go_to_jumps_by_part_of_a_name_or_by_number_when_blinded(review_dir):
    panel = _panel(review_dir)
    panel.goto_input.value = "t0003"
    assert panel.session.current.file == "T0003.xyz"
    assert panel.goto_input.value == ""
    panel.goto_input.value = "0002"
    assert panel.session.current.file == "T0002.xyz"
    panel.goto_input.value = "nothing"
    assert panel.session.current.file == "T0002.xyz"
    assert "nothing matches" in panel.status_html.value
    blind = _panel(review_dir, blind_cb=True)
    assert blind.go_to("#3") and blind.session.position == 2
    assert not blind.go_to("T0001")


def test_a_long_file_list_shows_a_window_around_the_current_file(tmp_path):
    d = tmp_path / "many"
    d.mkdir()
    for k in range(1, 2 * tab_review.FRAME_WINDOW + 51):
        (d / f"T{k:04d}.xyz").write_text(_xyz(f"T{k:04d}", 1))
    panel = _panel(d)
    assert len(panel.file_list.options) == 2 * tab_review.FRAME_WINDOW
    panel.go_to(f"T{2 * tab_review.FRAME_WINDOW + 50:04d}")
    assert panel.file_list.value == f"T{2 * tab_review.FRAME_WINDOW + 50:04d}.xyz"
