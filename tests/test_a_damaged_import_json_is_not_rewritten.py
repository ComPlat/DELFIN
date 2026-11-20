"""Red first: a damaged delfin_import.json must survive a written command.

``_append_to_json`` reads the job's ``delfin_import.json`` leniently and
writes the whole structure back: one damaged byte drops every recorded
define command the moment the next command is appended. The file is the
only record of what was typed into Define, and the replay button plays
exactly what it holds.
"""
import json

from delfin.dashboard import tab_turbomole_builder as TTB


def test_append_after_damage_keeps_the_old_commands(tmp_path, caplog):
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
    # The report is a log warning, not stdout: in a served dashboard
    # nothing reads the kernel's stdout, and the Settings tab's
    # diagnosis is the route a person sees (it lists the quarantine
    # file). The record is asserted here, where its wording lives.
    said = "\n".join(r.getMessage() for r in caplog.records)
    assert "delfin_import.json" in said, (
        "the loss must be said, not silent")
    assert "moved to" in said and "kept aside" in said
