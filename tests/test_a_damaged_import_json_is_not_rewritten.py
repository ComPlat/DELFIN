"""Red first: a damaged delfin_import.json must survive a written command.

``_append_to_json`` reads the job's ``delfin_import.json`` leniently and
writes the whole structure back: one damaged byte drops every recorded
define command the moment the next command is appended. The file is the
only record of what was typed into Define, and the replay button plays
exactly what it holds.
"""
import json

from delfin.dashboard import tab_turbomole_builder as TTB


def test_append_after_damage_keeps_the_old_commands(tmp_path, capsys):
    json_file = tmp_path / "delfin_import.json"
    json_file.write_text(json.dumps(
        {"define_commands": ["a coord", "scf"], "settings": {}}),
        encoding="utf-8")
    json_file.write_text("{not json", encoding="utf-8")

    data = TTB._read_import_json(json_file)
    data.setdefault('define_commands', []).append("browzig")
    json_file.write_text(json.dumps(data, indent=2))

    kept = sorted(p for p in tmp_path.iterdir()
                  if p.name.startswith(json_file.name + ".unreadable-"))
    assert kept, "the damaged file was overwritten with no copy anywhere"
    assert "not json" in kept[0].read_text(encoding="utf-8")
    data = json.loads(json_file.read_text(encoding="utf-8"))
    assert data["define_commands"] == ["browzig"]
    out = capsys.readouterr().out
    assert "delfin_import.json" in out, "the loss must be said, not silent"
