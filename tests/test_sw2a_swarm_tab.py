"""SW2a phase 1: the swarm tab exists, registers, and is pure data.

The tab is a VIEW over delfin/agent/swarm (SW1) functions: no orchestration
logic of its own. The testable logic (spec -> display rows, status report ->
cards, triage queue -> panel) lives in module-level helpers that take plain
dicts and return plain dicts/strings, so these tests run without ipywidgets
and without a browser, like tests/test_tab_agent_activity.py.

RED CONTROL (phase 1): on the commit before the tab module exists every test
here fails with ModuleNotFoundError/AttributeError — nothing fakes the module.
"""
from __future__ import annotations

import pytest

from delfin.agent.swarm.status import STATUS_EXAMPLE, STATUS_SCHEMA_VERSION


def _import_tab():
    from delfin.dashboard import tab_swarm as tab
    return tab


# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------


def test_tab_module_imports_and_has_create_tab():
    tab = _import_tab()
    assert callable(tab.create_tab)
    assert hasattr(tab, "HAS_WIDGETS")


def test_swarm_tab_registered_in_the_dashboard_tab_list():
    """Pin: the tab list on main is a SUBSET of the tab list after SW2a —
    the swarm tab is added next to the agent tabs, nothing else moves."""
    from delfin.dashboard import tab_registry as R

    _import_tab()  # registration side effect
    ids = {t.id for t in R.registered_tabs()}
    assert "swarm" in ids
    titles = {t.title for t in R.registered_tabs()}
    assert "Swarm" in titles


# ---------------------------------------------------------------------------
# Status report -> board cards (pure data, no widgets)
# ---------------------------------------------------------------------------


def test_status_example_renders_to_card_rows():
    """SW1's own STATUS_EXAMPLE document must render without loss: one card
    per package, phase and member health carried through."""
    tab = _import_tab()
    rows = tab.board_cards(STATUS_EXAMPLE)
    assert [r["name"] for r in rows] == ["core"]
    core = rows[0]
    assert core["phase"] == "in_review"
    assert [m["name"] for m in core["members"]] == ["s6", "s12"]
    assert [m["health"] for m in core["members"]] == ["active", "active"]
    assert core["members"][1]["waits_on"] == "builder commit"


def test_board_cards_pass_through_every_phase_vocabulary():
    """No phase from SW1's schema is renamed or dropped by the view."""
    tab = _import_tab()
    from delfin.agent.swarm.status import STATUS_SCHEMA
    phases = STATUS_SCHEMA["properties"]["packages"]["items"]["properties"][
        "phase"]["enum"]
    for phase in phases:
        doc = {
            "schema_version": STATUS_SCHEMA_VERSION,
            "goal": "g",
            "packages": [{"name": "p", "phase": phase, "members": []}],
            "waiting_for_user": [],
        }
        rows = tab.board_cards(doc)
        assert rows[0]["phase"] == phase


def test_board_cards_rejects_a_foreign_shape_loudly():
    """A document that is not a swarm status report is not rendered into a
    quiet empty board: the view fails loudly so the user sees the mismatch."""
    tab = _import_tab()
    with pytest.raises((KeyError, TypeError, ValueError)):
        tab.board_cards({"goal": "not a status doc"})


# ---------------------------------------------------------------------------
# Guided flow: spec text -> plan proposal (pure data)
# ---------------------------------------------------------------------------


def test_spec_errors_come_back_as_inline_messages():
    """A validation failure (SW1 SpecError) becomes a message dict, not a
    traceback — the guided flow shows it inline next to the editable table."""
    tab = _import_tab()
    result = tab.validate_spec_text("spec_version: 1\npackages: []\n")
    assert result["ok"] is False
    assert result["error"]  # prose from SW1
    assert result["packages"] == []


def test_a_valid_spec_text_validates_and_lists_packages():
    tab = _import_tab()
    spec_text = (
        "spec_version: 1\n"
        'goal: "demo"\n'
        "packages:\n"
        "  - name: core\n"
        '    owned: ["delfin/agent/swarm/**"]\n'
        "    roles:\n"
        "      builder: true\n"
        "      reviewer: true\n"
    )
    result = tab.validate_spec_text(spec_text)
    assert result["ok"] is True
    assert result["error"] == ""
    assert [p["name"] for p in result["packages"]] == ["core"]
    assert [p["roles"] for p in result["packages"]] == [
        {"builder": True, "reviewer": True}]


def test_plan_preview_is_the_dry_run_text():
    """The dry-run preview shows exactly what would start: the same text the
    CLI renders via SW1's render_dry_run, plus the goal."""
    tab = _import_tab()
    spec_text = (
        "spec_version: 1\n"
        'goal: "demo goal"\n'
        "packages:\n"
        "  - name: core\n"
        '    owned: ["delfin/agent/swarm/**"]\n'
        "    roles:\n"
        "      builder: true\n"
    )
    result = tab.plan_preview(spec_text)
    assert result["ok"] is True
    assert "demo goal" in result["preview"]
    assert "Nothing started (dry run)" in result["preview"]


def test_protected_scope_error_is_shown_inline():
    """A glob over a protected path is refused by SW1 and must appear as an
    inline validation message in the guided flow (target item 1)."""
    tab = _import_tab()
    spec_text = (
        "spec_version: 1\n"
        'goal: "bad"\n'
        "packages:\n"
        "  - name: evil\n"
        '    owned: ["delfin/agent/sandbox.py"]\n'
        "    roles:\n"
        "      builder: true\n"
    )
    result = tab.validate_spec_text(spec_text)
    assert result["ok"] is False
    assert "protected" in result["error"].lower()


# ---------------------------------------------------------------------------
# Triage queue -> panel rows (pure data)
# ---------------------------------------------------------------------------


def test_triage_list_sorts_escalations_first():
    tab = _import_tab()
    roster = {"a-builder", "a-reviewer"}
    queue = [
        {"id": "m1", "kind": "session_message", "from": "a-reviewer",
         "member": "a-builder", "text": "hi"},
        {"id": "t1", "kind": "push_request", "from": "a-builder",
         "member": "a-builder", "text": "push me"},
    ]
    rows = tab.triage_rows(queue, roster=roster)
    assert [r["id"] for r in rows] == ["t1", "m1"]
    assert rows[0]["action"] == "escalate"
    assert rows[1]["action"] == "auto_answer"


def test_triage_rows_carry_the_dialog_text():
    """The user decides on the item's real text, not a summary."""
    tab = _import_tab()
    item = {"id": "m2", "kind": "session_message", "from": "a-reviewer",
            "member": "a-builder", "text": "Please approve branch merge."}
    rows = tab.triage_rows([item], roster={"a-builder", "a-reviewer"})
    assert rows[0]["text"] == "Please approve branch merge."


def test_pause_resume_absent_without_sw1_verbs():
    """Fail-closed (finding sent to QS): SW1 ships no pause/resume module.
    The view must NOT render buttons that would lie about what they do.
    If SW1 adds the verbs, this test is rewritten along its intent."""
    tab = _import_tab()
    verbs = tab.available_actions()
    assert "pause" not in verbs
    assert "resume" not in verbs
    assert "start" in verbs          # via SW1 start_members
    assert "merge_plan" in verbs     # via SW1 merge_plan (shown as command)
