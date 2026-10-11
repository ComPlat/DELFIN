"""Package L4 (phase 2): the dashboard user-model panel.

The panel is pure functions over the ``UserModel`` store: visible rows,
an escaped text render, "forget this" and "forget everything", all through
the store so a forgotten entry is gone from what the next session loads.
"""

from delfin.agent import user_model
from delfin.dashboard import user_model_panel


def _populate(tmp_path, *, n=2):
    m = user_model.UserModel(path=tmp_path / "user_model.json")
    m.update_from_user_message(
        "I am a quantum chemist working on ORCA", source_role="user"
    )
    m.update_from_user_message(
        "I prefer pytest over unittest", source_role="user"
    )
    return str(m.path())


def test_panel_lists_stored_rows(tmp_path):
    path = _populate(tmp_path)
    rows = user_model_panel.collect_rows(path)
    fields = {r["field"] for r in rows}
    assert "role" in fields
    assert "preferences" in fields
    for row in rows:
        assert row["value"]
        assert row["source_role"] == "user"
        assert row["date"]


def test_panel_render_is_html_safe(tmp_path):
    path = _populate(tmp_path)
    # Store a value that would break an unescaped <tag> in the panel.
    m = user_model.UserModel(path=path)
    m.update_from_user_message(
        "I want answers like <b>bold</b>", source_role="user"
    )
    out = user_model_panel.render(path)
    assert "<b>bold</b>" not in out  # raw tag must not reach the panel


def test_panel_empty_when_nothing_stored(tmp_path):
    out = user_model_panel.render(str(tmp_path / "user_model.json"))
    assert "(no user model" in out


def test_panel_forget_entry_removes_from_next_load(tmp_path):
    path = _populate(tmp_path)
    # First entry is the role (index 0 of role bucket).
    res = user_model_panel.forget_entry("role", 0, path)
    assert res["removed"] is True
    m2 = user_model.UserModel(path=path)
    # The forgotten role entry is gone from what a new instance loads ...
    assert m2.entries().get("role", []) == []
    # ... while the untouched field survives.
    assert m2.entries().get("preferences", []) != []


def test_panel_forget_entry_absent_reports_false(tmp_path):
    path = _populate(tmp_path)
    res = user_model_panel.forget_entry("expertise", 0, path)
    assert res["removed"] is False


def test_panel_forget_everything(tmp_path):
    path = _populate(tmp_path)
    res = user_model_panel.forget_everything(path)
    assert res["removed"] >= 2
    m2 = user_model.UserModel(path=path)
    for field in user_model.FIELDS:
        assert m2.entries().get(field, []) == []
