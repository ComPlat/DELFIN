"""Tests for module-level helper functions in tab_agent.py.

Only pure functions are exercised here — anything that touches ipywidgets
state needs a notebook kernel.
"""
from __future__ import annotations

from delfin.dashboard.tab_agent import (
    _SLASH_COMMANDS,
    _filter_slash_commands,
    _format_solo_domain_state,
    _render_slash_palette_html,
    _render_subagent_pane_html,
    _render_todo_pane_html,
)


# ---------------------------------------------------------------------------
# TodoWrite plan pane renderer
# ---------------------------------------------------------------------------

def test_todo_pane_empty_returns_empty_string():
    """No todos → empty HTML so the caller can hide the widget."""
    assert _render_todo_pane_html([]) == ""


def test_todo_pane_renders_progress_counter():
    todos = [
        {"content": "A", "status": "completed", "activeForm": "Doing A"},
        {"content": "B", "status": "in_progress", "activeForm": "Doing B"},
        {"content": "C", "status": "pending", "activeForm": "Doing C"},
    ]
    html = _render_todo_pane_html(todos)
    assert "1/3 done" in html
    assert "1 active" in html


def test_todo_pane_uses_active_form_for_in_progress():
    todos = [
        {"content": "Run pytest", "status": "in_progress",
         "activeForm": "Running pytest"},
    ]
    html = _render_todo_pane_html(todos)
    # in_progress → activeForm shown
    assert "Running pytest" in html
    # plain content not present
    assert "Run pytest</span>" not in html


def test_todo_pane_uses_content_for_completed_and_pending():
    todos = [
        {"content": "Done thing", "status": "completed",
         "activeForm": "Doing thing"},
        {"content": "Future thing", "status": "pending",
         "activeForm": "Doing future"},
    ]
    html = _render_todo_pane_html(todos)
    assert "Done thing" in html
    assert "Future thing" in html
    # activeForm strings should NOT appear (they are only for in_progress)
    assert "Doing thing" not in html
    assert "Doing future" not in html


def test_todo_pane_includes_progress_bar_width():
    todos = [
        {"content": "x", "status": "completed", "activeForm": "x"},
        {"content": "y", "status": "completed", "activeForm": "y"},
        {"content": "z", "status": "pending", "activeForm": "z"},
        {"content": "w", "status": "pending", "activeForm": "w"},
    ]
    html = _render_todo_pane_html(todos)
    # 2 of 4 = 50%
    assert "width:50%" in html


def test_todo_pane_escapes_content():
    """Content from the agent must be HTML-escaped to prevent injection."""
    todos = [
        {"content": "<script>alert('x')</script>", "status": "pending",
         "activeForm": "boom"},
    ]
    html = _render_todo_pane_html(todos)
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_todo_pane_truncates_long_content():
    long_content = "x" * 500
    todos = [
        {"content": long_content, "status": "pending", "activeForm": "do"},
    ]
    html = _render_todo_pane_html(todos)
    # Renderer caps at 160 chars per todo
    assert "x" * 161 not in html
    assert "x" * 100 in html  # but reasonable length still present


def test_todo_pane_handles_missing_keys_gracefully():
    """Robust against partial todos (missing status/content)."""
    todos = [
        {"content": "no status"},                     # missing status
        {"status": "pending"},                         # missing content
        {"status": "in_progress"},                     # missing both content + activeForm
    ]
    html = _render_todo_pane_html(todos)
    # Should not raise; should still render header
    assert "0/3 done" in html
    assert "no status" in html


def test_todo_pane_no_active_means_no_active_badge():
    """Header shows '· N active' only when in_progress > 0."""
    todos = [
        {"content": "a", "status": "completed", "activeForm": "a"},
        {"content": "b", "status": "pending", "activeForm": "b"},
    ]
    html = _render_todo_pane_html(todos)
    assert "active" not in html


def test_todo_pane_status_glyphs_present():
    todos = [
        {"content": "a", "status": "completed", "activeForm": "a"},
        {"content": "b", "status": "in_progress", "activeForm": "b"},
        {"content": "c", "status": "pending", "activeForm": "c"},
    ]
    html = _render_todo_pane_html(todos)
    assert "✓" in html  # completed
    assert "→" in html  # in_progress
    assert "○" in html  # pending


# ---------------------------------------------------------------------------
# Slash-command catalog + palette helpers
# ---------------------------------------------------------------------------

def test_slash_catalog_is_nonempty():
    assert len(_SLASH_COMMANDS) > 30
    # Every entry follows the (category, command, summary, has_args) shape
    for entry in _SLASH_COMMANDS:
        assert len(entry) == 4
        cat, cmd, summary, has_args = entry
        assert isinstance(cat, str) and cat
        assert cmd.startswith("/")
        assert isinstance(summary, str) and summary
        assert isinstance(has_args, bool)


def test_slash_catalog_no_duplicate_commands():
    """Each command literal should be unique to avoid ambiguous insertions."""
    seen = [cmd for _cat, cmd, _summary, _has in _SLASH_COMMANDS]
    assert len(seen) == len(set(seen)), (
        f"Duplicate commands in catalog: {sorted(seen)}"
    )


def test_filter_empty_returns_all():
    out = _filter_slash_commands(_SLASH_COMMANDS, "")
    assert len(out) == len(_SLASH_COMMANDS)


def test_filter_blank_returns_all():
    """Whitespace-only query is treated as empty."""
    out = _filter_slash_commands(_SLASH_COMMANDS, "   ")
    assert len(out) == len(_SLASH_COMMANDS)


def test_filter_matches_command_literal():
    out = _filter_slash_commands(_SLASH_COMMANDS, "/help")
    assert any(cmd == "/help" for _cat, cmd, _s, _h in out)


def test_filter_matches_summary_text():
    out = _filter_slash_commands(_SLASH_COMMANDS, "permission profile")
    assert any("perm" in cmd for _cat, cmd, _s, _h in out)


def test_filter_matches_category():
    out = _filter_slash_commands(_SLASH_COMMANDS, "git")
    cmds = {cmd for _cat, cmd, _s, _h in out}
    assert "/git status" in cmds


def test_filter_case_insensitive():
    out_lower = _filter_slash_commands(_SLASH_COMMANDS, "control")
    out_upper = _filter_slash_commands(_SLASH_COMMANDS, "CONTROL")
    assert {c for _cat, c, _s, _h in out_lower} == {c for _cat, c, _s, _h in out_upper}


def test_filter_no_match_returns_empty():
    out = _filter_slash_commands(_SLASH_COMMANDS, "thiswillnotmatchanything12345")
    assert out == []


def test_render_palette_groups_by_category():
    sample = [
        ("Git", "/git status", "show status", False),
        ("Git", "/git log", "show log", False),
        ("Session", "/help", "help", False),
    ]
    html = _render_slash_palette_html(sample)
    assert "Git" in html
    assert "Session" in html
    assert "/git status" in html


def test_render_palette_empty_shows_no_match_message():
    html = _render_slash_palette_html([], query="xyz")
    assert "No commands match" in html


def test_render_palette_marks_args_with_trailing_space():
    """Commands with has_args=True should insert with a trailing space."""
    sample = [("Calc", "/calc ls", "list", True)]
    html = _render_slash_palette_html(sample)
    # data-command attribute carries the insert-text
    assert 'data-command="/calc ls "' in html


def test_render_palette_marks_no_args_without_trailing_space():
    sample = [("Session", "/help", "help", False)]
    html = _render_slash_palette_html(sample)
    assert 'data-command="/help"' in html


def test_render_palette_escapes_query():
    html = _render_slash_palette_html([], query="<script>")
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


# ---------------------------------------------------------------------------
# Subagent pane renderer
# ---------------------------------------------------------------------------

def test_subagent_pane_empty_returns_empty_string():
    assert _render_subagent_pane_html([]) == ""


def test_subagent_pane_running_call_has_spinner_glyph():
    calls = [{
        "subagent_type": "general-purpose",
        "description": "research X",
        "prompt": "Find papers on Y",
        "status": "running",
        "output": "",
    }]
    html = _render_subagent_pane_html(calls)
    assert "↻" in html  # running glyph
    assert "1 running" in html
    assert "general-purpose" in html
    assert "research X" in html


def test_subagent_pane_done_shows_collapsible_output():
    calls = [{
        "subagent_type": "Explore",
        "description": "scan repo",
        "prompt": "find files",
        "status": "done",
        "output": "Found 12 files",
    }]
    html = _render_subagent_pane_html(calls)
    assert "✓" in html  # done glyph
    assert "<details" in html
    assert "Found 12 files" in html


def test_subagent_pane_no_running_omits_active_badge():
    calls = [{
        "subagent_type": "general-purpose",
        "description": "x", "prompt": "p",
        "status": "done", "output": "ok",
    }]
    html = _render_subagent_pane_html(calls)
    assert "running" not in html


def test_subagent_pane_escapes_user_text():
    calls = [{
        "subagent_type": "general-purpose",
        "description": "<script>alert('x')</script>",
        "prompt": "<b>bold</b>",
        "status": "running",
        "output": "<script>",
    }]
    html = _render_subagent_pane_html(calls)
    assert "<script>alert" not in html
    assert "&lt;script&gt;" in html


def test_subagent_pane_truncates_long_output_in_preview():
    long_output = "x" * 5000
    calls = [{
        "subagent_type": "Explore", "description": "d", "prompt": "p",
        "status": "done", "output": long_output,
    }]
    html = _render_subagent_pane_html(calls)
    # Renderer caps preview at 600 chars
    assert "x" * 700 not in html
    # But the chars-total label reflects the real size
    assert "5000 chars" in html


def test_subagent_pane_handles_multiple_calls():
    calls = [
        {"subagent_type": "Explore", "description": "a", "prompt": "p1",
         "status": "done", "output": "result a"},
        {"subagent_type": "general-purpose", "description": "b", "prompt": "p2",
         "status": "running", "output": ""},
    ]
    html = _render_subagent_pane_html(calls)
    assert "2 call(s)" in html
    assert "1 running" in html
    assert "result a" in html


def test_subagent_pane_handles_missing_keys_gracefully():
    """Robust against partial dicts."""
    calls = [{}]  # totally empty
    html = _render_subagent_pane_html(calls)
    assert "1 call(s)" in html
    # Default subagent_type is shown
    assert "agent" in html.lower()


# ---------------------------------------------------------------------------
# Solo-mode domain-state formatter
# ---------------------------------------------------------------------------

def test_solo_domain_state_empty_returns_empty_string():
    assert _format_solo_domain_state({}) == ""


def test_solo_domain_state_skips_empty_values():
    """Snapshot with only empty/None values yields no block."""
    snap = {
        "calc_dir": "",
        "selected": None,
        "control": {},
        "orca_builder": {},
        "job_summary": "",
        "workspace_files": [],
        "active_tab": "",
        "perm_profile": "",
    }
    assert _format_solo_domain_state(snap) == ""


def test_solo_domain_state_starts_with_header():
    out = _format_solo_domain_state({"calc_dir": "/c"})
    assert out.startswith("--- Domain State ---")
    assert "calc_dir: /c" in out


def test_solo_domain_state_renders_calc_and_selected():
    out = _format_solo_domain_state({
        "calc_dir": "/home/user/calc",
        "selected": "Jerome/complexes_1_200/job1",
    })
    assert "calc_dir: /home/user/calc" in out
    assert "selected: Jerome/complexes_1_200/job1" in out


def test_solo_domain_state_renders_control_compact():
    out = _format_solo_domain_state({
        "control": {
            "functional": "PBE0",
            "main_basisset": "def2-TZVP",
            "PAL": 40,
            "charge": 0,
            "multiplicity": 1,
            "solvent": "DMSO",
        },
    })
    # Compact "FUNC/BASIS" plus key=value
    assert "control: PBE0/def2-TZVP" in out
    assert "PAL=40" in out
    assert "multiplicity=1" in out
    assert "solvent=DMSO" in out
    # charge=0 is the default and should be SKIPPED (treated as empty)
    assert "charge=0" not in out


def test_solo_domain_state_renders_orca_builder():
    out = _format_solo_domain_state({
        "orca_builder": {"method": "BP86", "basis": "def2-SVP", "mult": 3},
    })
    assert "orca_builder: BP86/def2-SVP" in out
    assert "mult=3" in out


def test_solo_domain_state_method_only_no_basis():
    out = _format_solo_domain_state({
        "control": {"functional": "PBE0"},
    })
    assert "control: PBE0" in out
    # No slash because no basis
    assert "PBE0/" not in out


def test_solo_domain_state_renders_job_summary():
    out = _format_solo_domain_state({"job_summary": "2 RUNNING, 5 PENDING"})
    assert "jobs: 2 RUNNING, 5 PENDING" in out


def test_solo_domain_state_workspace_files_truncated_at_eight():
    files = [f"file_{i}.csv" for i in range(15)]
    out = _format_solo_domain_state({"workspace_files": files})
    assert "file_0.csv" in out
    assert "file_7.csv" in out
    assert "file_8.csv" not in out
    assert "(+7 more)" in out


def test_solo_domain_state_active_tab_and_perms():
    out = _format_solo_domain_state({
        "active_tab": "DELFIN Agent",
        "perm_profile": "ask_all",
    })
    assert "active_tab: DELFIN Agent" in out
    assert "perms: ask_all" in out


def test_solo_domain_state_handles_non_dict_control_gracefully():
    """If control isn't a dict (corrupt state), the formatter must not crash."""
    out = _format_solo_domain_state({
        "calc_dir": "/c", "control": "not a dict",
    })
    assert "calc_dir: /c" in out
    # No control line emitted
    assert "control:" not in out


def test_solo_domain_state_full_snapshot_well_formed():
    out = _format_solo_domain_state({
        "calc_dir": "/calc",
        "selected": "job_x/orca.out",
        "control": {"functional": "PBE0", "main_basisset": "def2-TZVP", "PAL": 8},
        "orca_builder": {"method": "BP86", "basis": "def2-SVP"},
        "job_summary": "1 RUNNING",
        "workspace_files": ["a.csv", "b.png"],
        "active_tab": "Calculations",
        "perm_profile": "repo_free",
    })
    # All eight lines are present in order
    expected_keys = [
        "calc_dir:", "selected:", "control:", "orca_builder:",
        "jobs:", "workspace:", "active_tab:", "perms:",
    ]
    last_index = -1
    for key in expected_keys:
        idx = out.find(key)
        assert idx > last_index, f"missing or out-of-order: {key}"
        last_index = idx
