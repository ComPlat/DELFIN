"""Red first: a damaged branch_list.json must survive a pin() call.

The write path is ``BL._update`` -> ``_load_all`` -> ``_save_all``: the
lenient read turns a damaged store into ``{}``, so every pinned and
hidden branch is gone after the next pin, with no copy anywhere.
"""
import json

import pytest

from delfin.dashboard import branch_list as BL


@pytest.fixture
def prefs_file(tmp_path, monkeypatch):
    path = tmp_path / "branch_list.json"
    path.write_text(json.dumps(
        {"repos": {"r1": {"pinned": ["feature-a"],
                          "hidden": ["wip-b"]}}}), encoding="utf-8")
    monkeypatch.setattr(BL, "_prefs_path", lambda: path)
    return path


def test_pin_after_damage_keeps_the_old_lists(prefs_file, monkeypatch):
    prefs_file.write_text("{not json", encoding="utf-8")
    assert BL.load_prefs("r1") == {"pinned": [], "hidden": []}  # reader lenient

    BL.pin("r1", "feature-c")          # the writer must not lose the store

    kept = sorted(p for p in prefs_file.parent.iterdir()
                  if p.name.startswith(prefs_file.name + ".unreadable-"))
    assert kept, "the damaged file was overwritten with no copy anywhere"
    assert "not json" in kept[0].read_text(encoding="utf-8")
    data = json.loads(prefs_file.read_text(encoding="utf-8"))
    key = BL.repo_key("r1")
    assert key in data["repos"]
    row = data["repos"][key]
    assert row["pinned"] == ["feature-c"]


def test_untouched_store_is_not_quarantined_by_a_read(prefs_file):
    BL.pin("r1", "feature-a")
    assert not list(prefs_file.parent.glob("*.unreadable-*"))
