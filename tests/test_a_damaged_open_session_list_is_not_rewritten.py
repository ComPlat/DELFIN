"""Red first: a damaged agent_open_sessions.json must survive a save.

``AS.save_open_sessions`` builds the whole file from the rows it is
given; it never reads. The loss the lenient ``load_open_sessions``
enabled was indirect: the startup read answers ``[]`` for a damaged
file, the restore drops every open session, and ``_persist`` (line
~3021) then saves that empty list over the store. This test pins the
direct writer path: rows handed in after a damaged read must land in a
file that still holds the damaged original somewhere.
"""
import json

import pytest

from delfin.dashboard import agent_sessions as AS


@pytest.fixture
def open_file(tmp_path, monkeypatch):
    path = tmp_path / "agent_open_sessions.json"
    path.write_text(json.dumps(
        {"sessions": [{"session_id": "s1", "workspace": "/ws",
                       "host": "h1"}]}), encoding="utf-8")
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", path)
    return path


def test_save_open_sessions_after_damage_keeps_the_old_rows(open_file, capsys):
    open_file.write_text("{not json", encoding="utf-8")
    assert AS.load_open_sessions() == []               # reader stays lenient

    AS.save_open_sessions([{"session_id": "s2", "workspace": "/ws2"}])

    kept = sorted(p for p in open_file.parent.iterdir()
                  if p.name.startswith(open_file.name + ".unreadable-"))
    assert kept, "the damaged file was overwritten with no copy anywhere"
    assert "not json" in kept[0].read_text(encoding="utf-8")
    data = json.loads(open_file.read_text(encoding="utf-8"))
    ids = [row["session_id"] for row in data["sessions"]]
    assert ids == ["s2"]
    out = capsys.readouterr().out
    assert "agent_open_sessions.json" in out, "the loss must be said, not silent"
