"""Red first: a damaged dashboard_header.json must survive save_prefs.

``HU.save_prefs`` reads through the lenient ``load_prefs`` and writes the
whole store back, so one damaged byte dropped the header's remembered
state -- measured shape, same as the branch list.
"""
import json

import pytest

from delfin.dashboard import header_ui as HU


@pytest.fixture
def prefs_file(tmp_path, monkeypatch):
    path = tmp_path / "dashboard_header.json"
    path.write_text(json.dumps({"hint_compact": True}), encoding="utf-8")
    monkeypatch.setattr(HU, "_prefs_path", lambda: path)
    return path


def test_save_prefs_after_damage_keeps_the_old_keys(prefs_file, capsys):
    prefs_file.write_text("{not json", encoding="utf-8")
    assert HU.load_prefs() == {}                       # reader stays lenient

    HU.save_prefs(ring_size=7)                         # the writer must not lose the store

    kept = sorted(p for p in prefs_file.parent.iterdir()
                  if p.name.startswith(prefs_file.name + ".unreadable-"))
    assert kept, "the damaged file was overwritten with no copy anywhere"
    assert "not json" in kept[0].read_text(encoding="utf-8")
    assert json.loads(prefs_file.read_text(encoding="utf-8")) == {"ring_size": 7}
    out = capsys.readouterr().out
    assert "dashboard_header.json" in out, "the loss must be said, not silent"
    assert "unreadable" in out or "moved" in out
