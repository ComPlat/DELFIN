"""The Review tab: frames of multi-frame XYZ files are rated, saved at once and resumed.

Everything here is synthetic: invented file names, a water molecule moved a little
from frame to frame.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from delfin import review as rv


def _xyz(stem: str, n_frames: int) -> str:
    blocks = []
    for k in range(n_frames):
        dz = 0.01 * k
        blocks.append(
            f"3\n{stem} frame{k} label-{k}\n"
            f"O 0.000 0.000 {0.117 + dz:.3f}\n"
            f"H 0.000 0.757 {-0.467 + dz:.3f}\n"
            f"H 0.000 -0.757 {-0.467 + dz:.3f}\n"
        )
    return "".join(blocks)


@pytest.fixture
def folder(tmp_path: Path) -> Path:
    d = tmp_path / "structures"
    d.mkdir()
    (d / "a_01.xyz").write_text(_xyz("a_01", 3))
    (d / "a_02.xyz").write_text(_xyz("a_02", 2))
    (d / "a_03.xyz").write_text(_xyz("a_03", 4))
    (tmp_path / "index.tsv").write_text("id\tsmiles\tn_atoms\na_01\tO\t3\na_02\tO\t3\n")
    return d


# ---------------------------------------------------------------- loading
def test_a_multi_frame_xyz_is_split_into_its_frames():
    frames = rv.parse_frames(_xyz("x", 3) + "\n", "x.xyz")
    assert [f.index for f in frames] == [0, 1, 2]
    assert frames[1].comment == "x frame1 label-1"
    assert frames[1].label == "label-1"
    assert frames[2].n_atoms == 3 and frames[2].coords.count("\n") == 2
    assert frames[0].xyz().startswith("3\nx frame0 label-0\nO ")


def test_a_folder_opens_every_xyz_in_file_and_frame_order(folder):
    s = rv.ReviewSession.open(folder)
    assert len(s) == 9
    assert s.order[:4] == [("a_01.xyz", 0), ("a_01.xyz", 1), ("a_01.xyz", 2), ("a_02.xyz", 0)]
    assert s.review_path == folder / "review_structures.json"
    assert s.data["files"]["a_01.xyz"]["index"]["smiles"] == "O"
    assert "frame 2/3" in s.display_name(1)


def test_a_single_file_is_reviewed_into_its_own_folder(folder):
    s = rv.ReviewSession.open(folder / "a_03.xyz", name="alice")
    assert len(s) == 4
    assert s.review_path == folder / "review_alice.json"


# ---------------------------------------------------------------- navigation
def test_navigation_stops_at_both_ends_and_finds_the_next_unrated(folder):
    s = rv.ReviewSession.open(folder)
    assert s.prev().key == ("a_01.xyz", 0)
    s.go(100)
    assert s.position == len(s) - 1
    s.go(0)
    s.rate("pass", advance=False)
    s.go(1)
    s.rate("pass", advance=False)
    s.go(0)
    assert s.next_unrated().key == ("a_01.xyz", 2)
    assert s.go_to_file("a_03.xyz").key == ("a_03.xyz", 0)


def test_rating_moves_to_the_next_unrated_frame(folder):
    s = rv.ReviewSession.open(folder)
    s.go(2)
    s.rate("block", ["bond length"])
    s.go(1)
    s.rate("pass")
    assert s.position == 3


# ---------------------------------------------------------------- saving and resuming
def test_every_rating_is_on_disk_at_once_with_everything_that_links_it(folder):
    s = rv.ReviewSession.open(folder)
    s.rate("block", ["contact/clash", "other"], "H too close")
    data = json.loads(s.review_path.read_text())
    rec = data["reviews"]["a_01.xyz"]["0"]
    assert rec["verdict"] == "block"
    assert rec["categories"] == ["contact/clash", "other"]
    assert rec["note"] == "H too close"
    assert rec["file"] == "a_01.xyz" and rec["frame"] == 0 and rec["label"] == "label-0"
    assert rec["sha256"] == rv.file_sha256(folder / "a_01.xyz")
    assert rec["smiles"] == "O"
    assert rec["timestamp"]


def test_reopening_resumes_and_keeps_comments(folder):
    s = rv.ReviewSession.open(folder)
    s.rate("pass")
    s.rate("block", ["wrong isomer"])
    s.set_file_comment("a_01.xyz", "first file looks fine")
    s.set_session_comment("first pass")
    again = rv.ReviewSession.open(folder)
    assert again.progress() == {"rated": 2, "total": 9, "pass": 1, "block": 1}
    assert again.record("a_01.xyz", 1)["categories"] == ["wrong isomer"]
    assert again.data["files"]["a_01.xyz"]["comment"] == "first file looks fine"
    assert again.data["comment"] == "first pass"


def test_a_frame_rated_again_keeps_the_earlier_rating(folder):
    s = rv.ReviewSession.open(folder)
    s.rate("pass", advance=False)
    s.rate("block", ["other"], advance=False)
    rec = s.record("a_01.xyz", 0)
    assert rec["verdict"] == "block"
    assert [h["verdict"] for h in rec["history"]] == ["pass"]


def test_a_changed_file_no_longer_counts_as_rated(folder):
    s = rv.ReviewSession.open(folder)
    s.rate("pass")
    (folder / "a_01.xyz").write_text(_xyz("a_01", 3).replace("0.757", "0.758"))
    again = rv.ReviewSession.open(folder)
    assert again.record("a_01.xyz", 0) is None
    assert any("changed" in w for w in again.warnings)
    assert again.data["reviews"]["a_01.xyz"]["0"]["verdict"] == "pass"


def test_an_unknown_verdict_is_refused(folder):
    s = rv.ReviewSession.open(folder)
    with pytest.raises(ValueError):
        s.rate("maybe")
    assert not s.review_path.exists()


# ---------------------------------------------------------------- blinded mode
def test_the_blinded_order_is_shuffled_the_same_way_every_time(folder):
    a = rv.ReviewSession.open(folder, blind=True, seed=7)
    b = rv.ReviewSession.open(folder, blind=True, seed=7)
    c = rv.ReviewSession.open(folder, blind=True, seed=8)
    plain = rv.ReviewSession.open(folder)
    assert a.order == b.order
    assert sorted(a.order) == sorted(plain.order)
    assert a.order != plain.order or c.order != plain.order
    assert rv.blinded_order(list(reversed(plain.order)), 7) == a.order


def test_the_blinded_mode_hides_names_and_findings_but_keeps_the_mapping(folder):
    (folder / "findings.json").write_text(json.dumps({"a_01.xyz": {"0": ["long bond"]}}))
    s = rv.ReviewSession.open(folder, blind=True, seed=3)
    assert s.display_name() == "structure #1"
    hidden = s.current
    s.go_to_file("a_01.xyz")
    assert s.findings_for(s.current) == []
    s.go(0)
    s.rate("block", ["other"])
    data = json.loads(s.review_path.read_text())
    rec = data["reviews"][hidden.file][str(hidden.index)]
    assert rec["blind_number"] == 1 and rec["seed"] == 3
    assert data["order"][0] == [hidden.file, hidden.index]
    assert "findings_shown" not in rec


def test_findings_are_shown_next_to_their_frame_in_the_normal_mode(folder):
    (folder / "findings.json").write_text(json.dumps({"a_02": {"1": ["short contact"]}}))
    s = rv.ReviewSession.open(folder)
    s.go_to_file("a_02.xyz")
    assert s.findings_for(s.current) == []
    s.next()
    assert s.findings_for(s.current) == ["short contact"]


# ---------------------------------------------------------------- CSV and CLI
def test_the_summary_counts_categories_and_writes_csv(folder, capsys):
    s = rv.ReviewSession.open(folder)
    s.rate("pass")
    s.rate("block", ["bond length", "contact/clash"], "two problems")
    s.rate("block", ["bond length"])
    s.export_csv()
    rows = list(csv.DictReader(open(s.review_path.with_suffix(".csv"))))
    assert [r["verdict"] for r in rows] == ["pass", "block", "block"]
    assert rows[1]["categories"] == "bond length;contact/clash"
    assert rows[1]["note"] == "two problems"

    from delfin.cli import main as delfin_main
    out_csv = folder / "out.csv"
    assert delfin_main(["review", "summary", str(s.review_path), "--csv", str(out_csv)]) == 0
    printed = capsys.readouterr().out
    assert "rated 3 / 9" in printed and "pass 1" in printed and "block 2" in printed
    line = next(row for row in printed.splitlines() if row.startswith("bond length"))
    assert line.split()[-2:] == ["2", "2"]
    assert len(list(csv.DictReader(open(out_csv)))) == 3


def test_the_summary_of_a_missing_file_fails_cleanly(tmp_path, capsys):
    assert rv.main(["summary", str(tmp_path / "none.json")]) == 1


# ---------------------------------------------------------------- the tab, headless
def test_the_tab_rates_by_key_and_by_button(folder):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_registry
    from delfin.dashboard.tab_review import ReviewPanel

    assert "delfin.dashboard.tab_review" in tab_registry._BUILTIN_DYNAMIC_TABS
    panel = ReviewPanel(folder)
    assert panel.open(str(folder)) is not None
    panel.handle_key("1")
    panel.note_input.value = "too long"
    panel.handle_key("b")
    rec = panel.session.record("a_01.xyz", 0)
    assert rec["verdict"] == "block" and rec["categories"] == ["bond length"]
    assert rec["note"] == "too long"
    assert panel.session.position == 1
    assert panel.note_input.value == ""
    panel.handle_key("p")
    panel.handle_key("left")
    assert panel.session.position == 1
    assert panel.frame_list.options[0][0].startswith("✗")
    panel.multi_cb.value = False
    panel.category_buttons[0].value = True
    panel.category_buttons[1].value = True
    assert panel.selected_categories() == [rv.DEFAULT_CATEGORIES[1]]
    assert "2 / 9" in panel.progress_html.value
    panel.blind_cb.value = True
    assert panel.files_box.layout.display == "none"
    assert all(label.split()[-1].startswith("#") for label, _ in panel.frame_list.options)
    assert panel.export().is_file()
