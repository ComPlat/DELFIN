"""The Review tab: its keys reach it in the dashboard, and its screen is laid out for rating.

Synthetic data only: invented file names, a water molecule moved a little per frame.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("ipywidgets")

from delfin import review as rv  # noqa: E402
from delfin.dashboard import tab_review  # noqa: E402
from delfin.dashboard.tab_review import ReviewPanel  # noqa: E402


def _xyz(stem: str, n_frames: int) -> str:
    return "".join(
        f"3\n{stem} frame{k} label-{k}\n"
        f"O 0.000 0.000 {0.117 + 0.01 * k:.3f}\n"
        f"H 0.000 0.757 {-0.467 + 0.01 * k:.3f}\n"
        f"H 0.000 -0.757 {-0.467 + 0.01 * k:.3f}\n"
        for k in range(n_frames))


@pytest.fixture
def folder(tmp_path: Path) -> Path:
    d = tmp_path / "structures"
    d.mkdir()
    (d / "a_01.xyz").write_text(_xyz("a_01", 2))
    (d / "a_02.xyz").write_text(_xyz("a_02", 3))
    return d


class _Ctx:
    """The two parts of the dashboard context the tab uses for its script."""

    def __init__(self):
        import ipywidgets as widgets
        self.js_output = widgets.Output()
        self.init_js_parts: list[str] = []

    def add_init_js(self, script):
        self.init_js_parts.append(script)


# ---------------------------------------------------------------- keys in the dashboard
def test_the_key_script_goes_out_with_the_dashboard_start_up_scripts():
    """Appended to the shared output on its own, the dashboard's one start-up send
    cleared it again before the browser ran it, and no key reached the tab."""
    ctx = _Ctx()
    ReviewPanel(Path.cwd(), ctx=ctx)
    assert any("__delfinReviewKeys" in s for s in ctx.init_js_parts)
    assert not ctx.js_output.outputs


def test_escape_leaves_a_text_field_of_the_tab():
    js = tab_review._keys_js()
    assert "Escape" in js and "blur()" in js


# ---------------------------------------------------------------- layout
def test_the_setup_closes_and_the_rating_screen_opens_once_a_review_is_open(folder):
    panel = ReviewPanel(folder)
    assert panel.setup.selected_index == 0
    assert panel.work.layout.display == "none"
    assert all(b.disabled for b in (panel.pass_btn, panel.block_btn, *panel.category_buttons))
    assert panel.open(str(folder)) is not None
    assert panel.setup.selected_index is None
    assert panel.work.layout.display == ""
    assert not panel.pass_btn.disabled and not panel.category_buttons[0].disabled
    assert "0 / 5" in panel.progress_html.value
    assert "a_01.xyz" in panel.position_html.value


def test_the_list_folds_away_and_comes_back(folder):
    panel = ReviewPanel(folder)
    panel.open(str(folder))
    panel.list_toggle.value = False
    assert panel.sidebar.layout.display == "none"
    panel.list_toggle.value = True
    assert panel.sidebar.layout.display == ""


def test_a_rating_says_that_it_was_saved_and_a_chosen_category_stands_out(folder):
    panel = ReviewPanel(folder)
    panel.open(str(folder))
    panel.handle_key("2")
    assert panel.category_buttons[1].button_style == "warning"
    panel.handle_key("b")
    assert "saved" in panel.saved_html.value and "block" in panel.saved_html.value
    assert "a_01.xyz" in panel.saved_html.value
    assert panel.category_buttons[1].button_style == ""       # next frame, nothing chosen
    panel.handle_key("left")
    assert panel.category_buttons[1].button_style == "warning"
    assert "block" in panel.verdict_html.value


def test_categories_set_before_a_new_review_are_kept_and_a_resumed_one_keeps_its_own(folder):
    panel = ReviewPanel(folder)
    panel.categories_input.value = "TOPO\nCLASH\nH"
    panel.apply_categories()
    panel.open(str(folder))
    assert panel.session.categories == ["TOPO", "CLASH", "H"]
    assert [b._review_category for b in panel.category_buttons] == ["TOPO", "CLASH", "H"]
    panel.handle_key("3")
    panel.handle_key("b")
    again = ReviewPanel(folder)                       # default list in the field
    again.open(str(folder))
    assert again.session.categories == ["TOPO", "CLASH", "H"]
    assert again.session.record("a_01.xyz", 0)["categories"] == ["H"]


def test_the_blinded_screen_names_no_file_and_no_label(folder):
    panel = ReviewPanel(folder)
    panel.blind_cb.value = True
    panel.seed_input.value = 7
    panel.open(str(folder))
    panel.handle_key("p")
    shown = " ".join([panel.position_html.value, panel.saved_html.value,
                      panel.verdict_html.value, panel.context_html.value,
                      *(label for label, _ in panel.frame_list.options)])
    assert "structure #" in shown
    for leak in ("a_01", "a_02", ".xyz", "label-"):
        assert leak not in shown
    assert panel.files_box.layout.display == "none"
    assert panel.file_comment_box.layout.display == "none"


def test_a_session_opened_without_a_list_uses_the_default_categories(folder):
    s = rv.ReviewSession.open(folder)
    assert s.categories == list(rv.DEFAULT_CATEGORIES)
    t = rv.ReviewSession.open(folder, name="other", default_categories=["A", " ", "B"])
    assert t.categories == ["A", "B"]


# ---------------------------------------------------------------- keys that arrive late
def test_a_rating_key_pressed_for_a_frame_no_longer_shown_rates_nothing(folder):
    """The browser sends the position it showed with every key; the kernel has
    moved on after a rating, so a second P sent before the screen caught up
    (a double press, a held key) would otherwise rate a frame never seen."""
    panel = ReviewPanel(folder)
    panel.open(str(folder))
    shown = panel.session.position
    panel.keys_input.value = f"1 p {shown}"
    assert panel.session.record("a_01.xyz", 0)["verdict"] == "pass"
    assert panel.session.position != shown
    panel.keys_input.value = f"2 b {shown}"                 # still the old position
    assert panel.session.record("a_01.xyz", 1) is None
    assert f"delfin-review-at-{panel.session.position}" in panel.keys_input._dom_classes
    panel.keys_input.value = f"3 b {panel.session.position}"
    assert panel.session.record("a_01.xyz", 1)["verdict"] == "block"


def test_the_key_script_ignores_a_held_key_and_sends_the_shown_position():
    js = tab_review._keys_js()
    assert "ev.repeat" in js
    assert "delfin-review-at-" in js


def test_a_focused_frame_list_or_checkbox_does_not_swallow_the_keys():
    """Only a field one types into keeps the keys; after a click into the frame
    list or on a checkbox of the tab, P and B still rate."""
    js = tab_review._keys_js()
    assert "tag==='SELECT'||" not in js
    assert "checkbox" in js


# ---------------------------------------------------------------- blinded, after a switch
def test_switching_to_blinded_leaves_no_file_name_in_the_hidden_file_list(folder):
    panel = ReviewPanel(folder)
    panel.open(str(folder))
    panel.blind_cb.value = True
    labels = " ".join(str(label) for label, _ in panel.file_list.options)
    assert "a_01" not in labels and "a_02" not in labels


# ---------------------------------------------------------------- narrow screens
def test_the_viewer_follows_the_width_of_its_column(folder, monkeypatch):
    """A fixed 720 px canvas ran under the rating column below about 1250 px."""
    from delfin.dashboard import molecule_viewer
    seen = {}
    monkeypatch.setattr(molecule_viewer, "render_xyz_in_output",
                        lambda out, xyz, width=None, height=None: seen.update(width=width))
    panel = ReviewPanel(folder)
    panel.open(str(folder))
    assert seen["width"] == "100%"
    for field in (panel.path_input, panel.findings_input, panel.name_input):
        assert field.layout.max_width == "100%"
