"""SW2a adversarial tests, round 2 (reviewer nacht-s16).

Attack the guided flow's start path and the pure helpers' edges.
Expected RED against builder phase 1 (agent/s15-sw2ab @ a281e616a):

F1 (TOCTOU) The Start handler reloads the spec from the TEXTAREA, not
    from the previewed document (tab_swarm.py:329-340). Text edited
    between preview and Start launches a spec the user never saw. Fix
    shape: plan_preview records the previewed spec text (module-level
    last_previewed_spec_text) and Start must launch exactly that; a
    failed/absent preview must keep start disabled.
F3 Edge cases of the pure helpers: None, empty, unicode, huge inputs
    must fail loudly (raise ValueError) or return empty structures --
    never a silent wrong render.
"""

from __future__ import annotations

import pytest

from delfin.agent.swarm.status import STATUS_EXAMPLE


VALID_SPEC = (
    "spec_version: 1\n"
    'goal: "demo"\n'
    "packages:\n"
    "  - name: core\n"
    '    owned: ["delfin/agent/swarm/**"]\n'
    "    roles:\n"
    "      builder: true\n"
)

EVIL_SPEC = (
    "spec_version: 1\n"
    'goal: "demo"\n'
    "packages:\n"
    "  - name: evil\n"
    '    owned: ["delfin/agent/sandbox.py"]\n'
    "    roles:\n"
    "      builder: true\n"
)


# ---------------------------------------------------------------------------
# F1: Start must launch the previewed spec, not the textarea (TOCTOU)
# ---------------------------------------------------------------------------


def test_plan_preview_records_the_previewed_spec_text():
    """F1 part 1: the preview must record exactly what it previewed, so
    the Start handler can be bound to it. RED on phase 1: nothing is
    recorded (module has no last_previewed_spec_text)."""
    from delfin.dashboard import tab_swarm

    result = tab_swarm.plan_preview(VALID_SPEC)
    assert result["ok"] is True
    assert tab_swarm.last_previewed_spec_text == VALID_SPEC, (
        "plan_preview does not record the previewed spec text, so Start "
        "cannot be bound to what the user saw (TOCTOU, "
        "tab_swarm.py:329-340)")


def test_failed_preview_does_not_overwrite_a_good_recorded_spec():
    """F1 part 2: a FAILED preview must not leave a recorded spec that
    Start could launch -- the record must be cleared (or never set), so
    Start re-validates before every launch. RED on phase 1: the record
    does not exist yet, so nothing is cleared."""
    from delfin.dashboard import tab_swarm

    assert tab_swarm.plan_preview(VALID_SPEC)["ok"] is True
    good = tab_swarm.last_previewed_spec_text
    bad = tab_swarm.plan_preview(EVIL_SPEC)
    assert bad["ok"] is False
    assert tab_swarm.last_previewed_spec_text != good, (
        "a failed preview left the previously recorded spec in place; "
        "Start must not be able to launch a spec whose last preview "
        "failed")


def test_evil_spec_never_validates_through_the_tab_path():
    """A protected-path spec must fail validation through plan_preview
    exactly as through validate_spec_text (same SW1 rules, one path)."""
    from delfin.dashboard import tab_swarm

    result = tab_swarm.plan_preview(EVIL_SPEC)
    assert result["ok"] is False
    assert "protected" in result["error"].lower()
    assert result["preview"] == ""


# ---------------------------------------------------------------------------
# F3: helper edge cases
# ---------------------------------------------------------------------------


def test_board_cards_on_empty_and_none_reports():
    """None and an empty dict are not status reports: raise loudly, never
    render a silent empty board."""
    from delfin.dashboard import tab_swarm

    with pytest.raises(ValueError):
        tab_swarm.board_cards(None)
    with pytest.raises(ValueError):
        tab_swarm.board_cards({})


def test_board_cards_on_huge_and_unicode_documents():
    """A report with a 10_000-package list and unicode/hostile names must
    pass through the pure helper without crashing; escaping happens at
    the widget layer with _html.escape (regression pin)."""
    from delfin.dashboard import tab_swarm

    pkg = dict(STATUS_EXAMPLE["packages"][0])
    doc = {
        "schema_version": STATUS_EXAMPLE["schema_version"],
        "goal": "göαl \U0001f512 with <script>alert(1)</script>",
        "packages": [dict(pkg, name=f"p{i}") for i in range(10_000)],
        "waiting_for_user": [],
    }
    rows = tab_swarm.board_cards(doc)
    assert len(rows) == 10_000
    assert rows[0]["name"] == "p0"
    assert rows[-1]["name"] == "p9999"


def test_waiting_panel_passes_waiting_rows_through_verbatim():
    """A member waiting on the user reaches the panel with its real text
    (the user decides on what was said, not a summary)."""
    from delfin.dashboard import tab_swarm

    doc = dict(STATUS_EXAMPLE, waiting_for_user=[
        {"member": "s12", "waits_on": "builder commit <b>bold</b>"}])
    rows = tab_swarm.waiting_panel(doc)
    assert rows[0]["waits_on"] == "builder commit <b>bold</b>"


def test_triage_rows_rejects_a_non_list_queue():
    from delfin.dashboard import tab_swarm

    with pytest.raises(ValueError):
        tab_swarm.triage_rows("not a list")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        tab_swarm.triage_rows(None)  # type: ignore[arg-type]
