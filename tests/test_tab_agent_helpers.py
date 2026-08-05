"""Tests for module-level helper functions in tab_agent.py.

Only pure functions are exercised here — anything that touches ipywidgets
state needs a notebook kernel.
"""
from __future__ import annotations

from delfin.dashboard.tab_agent import (
    _SLASH_COMMANDS,
    _TAB_SUGGESTIONS,
    _build_full_transcript_handoff,
    _compute_turn_verdict,
    _extract_action_commands,
    _extract_denied_tool_name,
    _filter_slash_commands,
    _format_session_boot,
    _format_solo_domain_state,
    _is_structurally_blocked,
    _record_solo_turn_outcome,  # backwards-compat alias
    _record_turn_outcome,
    _render_action_confirmation_html,
    _render_artifact_inline,
    _render_molecule_to_png_b64,
    _render_slash_palette_html,
    _render_subagent_pane_html,
    _render_todo_pane_html,
    _render_xyz_summary,
    _should_show_action_confirmation,
    _suggestion_for_tab,
    _text_requests_confirmation,
    _tool_in_allowlist,
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


# ---------------------------------------------------------------------------
# D1 — confirmation extraction helpers
# ---------------------------------------------------------------------------

def test_extract_no_actions_returns_empty():
    assert _extract_action_commands("") == []
    assert _extract_action_commands("just plain text") == []


def test_extract_single_action():
    text = "Klar, ich mache das.\nACTION: /control key functional BP86\nFertig."
    assert _extract_action_commands(text) == ["/control key functional BP86"]


def test_extract_multiple_actions_in_order():
    text = (
        "Plan:\n"
        "ACTION: /control key functional PBE0\n"
        "ACTION: /control key main_basisset def2-TZVP\n"
        "ACTION: /control key PAL 8\n"
        "Done."
    )
    out = _extract_action_commands(text)
    assert out == [
        "/control key functional PBE0",
        "/control key main_basisset def2-TZVP",
        "/control key PAL 8",
    ]


def test_extract_strips_whitespace():
    text = "ACTION:    /submit   "
    assert _extract_action_commands(text) == ["/submit"]


def test_extract_ignores_lines_without_slash_command():
    text = "ACTION: not a command\nACTION: /jobs"
    assert _extract_action_commands(text) == ["/jobs"]


def test_extract_only_top_of_line_triggers():
    """Inline 'ACTION:' inside prose should NOT be picked up."""
    text = "I would do ACTION: /jobs but actually no."
    assert _extract_action_commands(text) == []


def test_confirmation_detection_german():
    assert _text_requests_confirmation("Soll ich das submitten?")
    assert _text_requests_confirmation("Möchtest du weitermachen?")
    assert _text_requests_confirmation("Darf ich den Job abschicken?")


def test_confirmation_detection_english():
    assert _text_requests_confirmation("Should I proceed with the recalc?")
    assert _text_requests_confirmation("Shall we go ahead?")
    assert _text_requests_confirmation("Do you want me to continue?")


def test_confirmation_detection_case_insensitive():
    assert _text_requests_confirmation("SHOULD I CONTINUE?")
    assert _text_requests_confirmation("Soll Ich Weitermachen?")


def test_no_confirmation_for_plain_status():
    assert not _text_requests_confirmation("")
    assert not _text_requests_confirmation("Done. Job 123 submitted.")
    assert not _text_requests_confirmation("Reading the file now.")


def test_no_confirmation_for_question_without_intent():
    """A plain '?' or rhetorical question is not a confirmation request."""
    assert not _text_requests_confirmation("Why does this happen?")


def test_should_show_confirmation_when_actions_and_question():
    text = (
        "Plan: setze BP86.\n"
        "ACTION: /control key functional BP86\n"
        "Soll ich das ausführen?"
    )
    assert _should_show_action_confirmation(text) is True


def test_should_show_confirmation_skipped_without_question():
    text = (
        "Setze BP86.\n"
        "ACTION: /control key functional BP86\n"
        "Fertig."
    )
    # No confirmation phrase → auto-exec runs, no buttons shown
    assert _should_show_action_confirmation(text) is False


def test_should_show_confirmation_skipped_without_actions():
    """Confirmation phrase alone — without ACTIONs — never shows buttons."""
    text = "Soll ich weitermachen?"
    assert _should_show_action_confirmation(text) is False


def test_should_show_confirmation_empty_text():
    assert _should_show_action_confirmation("") is False


def test_should_show_confirmation_multiple_actions():
    text = (
        "Mein Plan:\n"
        "ACTION: /control key functional PBE0\n"
        "ACTION: /control key main_basisset def2-TZVP\n"
        "Should I proceed?"
    )
    assert _should_show_action_confirmation(text) is True


def test_render_action_confirmation_lists_commands():
    html = _render_action_confirmation_html([
        "/control key functional BP86",
        "/orca submit",
    ])
    assert "2 action" in html
    assert "/control key functional BP86" in html
    assert "/orca submit" in html


def test_render_action_confirmation_empty_returns_empty_string():
    assert _render_action_confirmation_html([]) == ""


def test_render_action_confirmation_escapes_user_text():
    html = _render_action_confirmation_html([
        "/control set <script>alert('x')</script>",
    ])
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_render_action_confirmation_truncates_very_long_commands():
    """Very long commands are capped so the row stays scannable."""
    long_cmd = "/control set " + "x" * 500
    html = _render_action_confirmation_html([long_cmd])
    assert "x" * 200 not in html  # cap is at 160 chars


# ---------------------------------------------------------------------------
# D2 — tab-change suggestions
# ---------------------------------------------------------------------------

def test_suggestion_known_actionable_tab():
    out = _suggestion_for_tab("Submit Job")
    assert out is not None
    assert "CONTROL" in out


def test_suggestion_recalc_tab():
    out = _suggestion_for_tab("Recalc")
    assert out is not None
    assert "recalc" in out.lower()


def test_suggestion_job_status_tab():
    out = _suggestion_for_tab("Job Status")
    assert out is not None
    assert "/jobs check" in out or "Job-Events" in out


def test_suggestion_calculations_tab_offers_skill():
    out = _suggestion_for_tab("Calculations")
    assert out is not None
    assert "/skill energy-table" in out


def test_suggestion_orca_builder_tab():
    out = _suggestion_for_tab("ORCA Builder")
    assert out is not None
    # Case-insensitive: the suggestion is a sentence now, so the words the
    # test cares about sit mid-sentence rather than capitalised.
    lowered = out.lower()
    assert "input" in lowered or "builder" in lowered


def test_suggestion_silent_tabs_return_none():
    # Tabs intentionally mapped to "" → no suggestion noise
    for tab in ("Settings", "Archive", "Remote Archive",
                "DELFIN Agent", "Agent Activity"):
        assert _suggestion_for_tab(tab) is None, f"{tab} should be silent"


def test_suggestion_unknown_tab_returns_none():
    assert _suggestion_for_tab("Some Future Tab") is None


def test_suggestion_empty_or_falsy_returns_none():
    assert _suggestion_for_tab("") is None
    assert _suggestion_for_tab(None) is None  # type: ignore[arg-type]


def test_tab_suggestions_table_is_complete():
    """Every dashboard tab the user can land on must be either in the
    suggestion table or explicitly silenced (empty string)."""
    # If new tabs appear in the dashboard, they should be added to
    # _TAB_SUGGESTIONS so this test stays meaningful.
    must_have = {
        "Submit Job", "Recalc", "Job Status", "Calculations",
        "Literature", "ORCA Builder", "Settings",
        "Archive", "DELFIN Agent", "Agent Activity",
    }
    missing = must_have - _TAB_SUGGESTIONS.keys()
    assert missing == set(), f"missing suggestions for: {missing}"


# ---------------------------------------------------------------------------
# Workspace artifact renderer (D4 — inline visualisation)
# ---------------------------------------------------------------------------

# A 1x1 transparent PNG (smallest valid PNG) for image tests
_TINY_PNG_BYTES = bytes.fromhex(
    "89504e470d0a1a0a0000000d49484452"
    "0000000100000001080600000000"  # IHDR — 1x1 RGBA
    "1f15c4890000000a49444154789c63000100000500010d0a2db40000000049454e44ae426082"
)


def test_artifact_unsupported_returns_none(tmp_path):
    p = tmp_path / "thing.bin"
    p.write_bytes(b"\x00\x01\x02")
    assert _render_artifact_inline(p) is None


def test_artifact_missing_path_returns_none(tmp_path):
    assert _render_artifact_inline(tmp_path / "nope.png") is None


def test_artifact_png_renders_base64_img(tmp_path):
    p = tmp_path / "plot.png"
    p.write_bytes(_TINY_PNG_BYTES)
    html = _render_artifact_inline(p)
    assert html is not None
    assert "<img" in html
    assert "data:image/png;base64," in html
    assert "plot.png" in html


def test_artifact_jpg_uses_jpeg_mime(tmp_path):
    """JPG and JPEG suffix both map to image/jpeg mime."""
    p = tmp_path / "x.jpg"
    p.write_bytes(b"\xff\xd8\xff\xe0fake")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "data:image/jpeg;base64," in html


def test_artifact_huge_png_skipped(tmp_path, monkeypatch):
    """Files over the 2 MB cap render a placeholder instead of inlining."""
    p = tmp_path / "big.png"
    p.write_bytes(b"x" * 2_500_000)
    html = _render_artifact_inline(p)
    assert html is not None
    assert "too large" in html
    # Real bytes must NOT be embedded
    assert "data:image/png;base64," not in html


def test_artifact_svg_inline(tmp_path):
    p = tmp_path / "fig.svg"
    p.write_text(
        "<?xml version='1.0'?>"
        "<svg xmlns='http://www.w3.org/2000/svg'><rect width='10' height='10'/></svg>"
    )
    html = _render_artifact_inline(p)
    assert html is not None
    assert "<svg" in html
    assert "<?xml" not in html  # XML declaration is stripped
    assert "fig.svg" in html


def test_artifact_csv_table_with_header(tmp_path):
    p = tmp_path / "energies.csv"
    p.write_text("name,energy\nA,-1.23\nB,-2.34\nC,-3.45\n")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "<table" in html
    assert "<thead" in html
    assert "name" in html and "energy" in html
    assert "-1.23" in html
    assert "energies.csv" in html


def test_artifact_csv_truncates_long_files(tmp_path):
    rows = ["col1,col2"] + [f"r{i},v{i}" for i in range(50)]
    p = tmp_path / "big.csv"
    p.write_text("\n".join(rows))
    html = _render_artifact_inline(p)
    assert html is not None
    # Caps at 12 rows
    assert "r0,v0" not in html.replace(" ", "")
    assert "r10" in html  # within cap
    assert "more rows" in html


def test_artifact_tsv_uses_tab_separator(tmp_path):
    p = tmp_path / "x.tsv"
    p.write_text("a\tb\nc\td\n")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "a" in html and "b" in html
    # Must not have split on commas
    assert "<td" in html


def test_artifact_csv_escapes_html():
    """CSV cells with HTML special chars must be escaped."""
    import tempfile, os
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", delete=False, encoding="utf-8",
    ) as f:
        f.write("name,value\n<script>,&hack\n")
        path = f.name
    try:
        html = _render_artifact_inline(path)
        assert html is not None
        assert "<script>" not in html
        assert "&lt;script&gt;" in html
    finally:
        os.unlink(path)


def test_artifact_json_pretty(tmp_path):
    p = tmp_path / "data.json"
    p.write_text('{"key": "value", "n": 42}')
    html = _render_artifact_inline(p)
    assert html is not None
    assert "<pre" in html
    assert "key" in html
    assert "data.json" in html


def test_artifact_json_huge_skipped(tmp_path):
    p = tmp_path / "big.json"
    p.write_text("x" * 7000)
    html = _render_artifact_inline(p)
    assert html is not None
    assert "too large" in html
    assert "<pre" not in html


# ---------------------------------------------------------------------------
# Molecule artifacts (SMI / MOL / SDF / XYZ)
# ---------------------------------------------------------------------------

import pytest as _pt

_HAS_RDKIT = True
try:
    from rdkit import Chem  # noqa: F401
except Exception:
    _HAS_RDKIT = False


@_pt.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")
def test_molecule_to_png_smiles_valid():
    b64 = _render_molecule_to_png_b64("c1ccccc1")
    assert b64 is not None
    assert len(b64) > 100  # actual PNG bytes


@_pt.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")
def test_molecule_to_png_invalid_smiles_returns_none():
    assert _render_molecule_to_png_b64("not_a_smiles_!!!") is None


@_pt.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")
def test_artifact_smi_renders_with_image(tmp_path):
    p = tmp_path / "benzene.smi"
    p.write_text("c1ccccc1 benzene")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "data:image/png;base64," in html
    assert "c1ccccc1" in html
    assert "benzene.smi" in html


@_pt.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")
def test_artifact_smi_invalid_falls_back_to_text(tmp_path):
    """Invalid SMILES still produce a labelled stub instead of None."""
    p = tmp_path / "bad.smi"
    p.write_text("nonsense_!!!")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "bad.smi" in html
    # No PNG, since rendering failed
    assert "data:image/png;base64," not in html


def test_artifact_smi_empty_returns_none(tmp_path):
    p = tmp_path / "empty.smi"
    p.write_text("")
    assert _render_artifact_inline(p) is None


@_pt.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")
def test_artifact_mol_renders_image(tmp_path):
    """Generate a real .mol file and verify it renders."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")  # ethanol
    p = tmp_path / "ethanol.mol"
    Chem.MolToMolFile(mol, str(p))
    html = _render_artifact_inline(p)
    assert html is not None
    assert "data:image/png;base64," in html
    assert "atoms" in html  # atom-count label


@_pt.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")
def test_artifact_sdf_renders_first_molecule(tmp_path):
    from rdkit import Chem
    writer = Chem.SDWriter(str(tmp_path / "set.sdf"))
    for smi in ("CCO", "c1ccccc1"):
        writer.write(Chem.MolFromSmiles(smi))
    writer.close()
    html = _render_artifact_inline(tmp_path / "set.sdf")
    assert html is not None
    assert "data:image/png;base64," in html


def test_artifact_mol_invalid_returns_none(tmp_path):
    p = tmp_path / "bad.mol"
    p.write_text("not a valid mol file")
    assert _render_artifact_inline(p) is None


def test_xyz_summary_basic(tmp_path):
    p = tmp_path / "water.xyz"
    p.write_text("3\nwater molecule\nO 0.0 0.0 0.0\nH 0.0 0.7 0.6\nH 0.0 -0.7 0.6\n")
    out = _render_xyz_summary(p)
    assert out is not None
    n_atoms, formula = out
    assert n_atoms == 3
    assert "H2" in formula
    assert "O" in formula


def test_xyz_summary_invalid_first_line(tmp_path):
    p = tmp_path / "bad.xyz"
    p.write_text("not a number\ncomment\nO 0 0 0\n")
    assert _render_xyz_summary(p) is None


def test_xyz_summary_empty(tmp_path):
    p = tmp_path / "empty.xyz"
    p.write_text("")
    assert _render_xyz_summary(p) is None


def test_artifact_xyz_renders_summary_block(tmp_path):
    p = tmp_path / "h2o.xyz"
    p.write_text("3\nwater\nO 0 0 0\nH 0 1 0\nH 1 0 0\n")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "h2o.xyz" in html
    assert "3 atoms" in html
    assert "H2" in html


def test_artifact_xyz_malformed_returns_none(tmp_path):
    p = tmp_path / "bad.xyz"
    p.write_text("garbage\n")
    assert _render_artifact_inline(p) is None


def test_artifact_filename_escaped(tmp_path):
    p = tmp_path / "weird<name>.png"
    # Some filesystems may reject this; if so skip the test
    try:
        p.write_bytes(_TINY_PNG_BYTES)
    except OSError:
        import pytest as _pytest
        _pytest.skip("filesystem rejected weird filename")
    html = _render_artifact_inline(p)
    assert html is not None
    assert "<name>" not in html
    assert "&lt;name&gt;" in html


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


# ---------------------------------------------------------------------------
# Step 5 — structural-block detection for permission denials
# ---------------------------------------------------------------------------

def test_extract_denied_tool_name_from_dict_string():
    """Real CLI denial payload is a literal dict string."""
    raw = "{'tool_name': 'Edit', 'tool_input': {'file_path': '/tmp/x.py'}}"
    assert _extract_denied_tool_name(raw) == "Edit"


def test_extract_denied_tool_name_from_dict_object():
    raw = {"tool_name": "Bash", "tool_input": {"command": "rm -rf /"}}
    assert _extract_denied_tool_name(raw) == "Bash"


def test_extract_denied_tool_name_returns_empty_when_missing():
    assert _extract_denied_tool_name("") == ""
    assert _extract_denied_tool_name("not a dict") == ""
    assert _extract_denied_tool_name("{'no_tool_key': 1}") == ""


def test_tool_in_allowlist_none_means_unrestricted():
    """A None allow-list = CLI default = every tool allowed."""
    assert _tool_in_allowlist("Edit", None) is True
    assert _tool_in_allowlist("Anything", None) is True


def test_tool_in_allowlist_empty_list_blocks_everything():
    assert _tool_in_allowlist("Edit", []) is False
    assert _tool_in_allowlist("Read", []) is False


def test_tool_in_allowlist_pattern_entries_count_as_base_tool():
    """``Bash(git *)`` enables the Bash tool."""
    cli = ["Read", "Glob", "Bash(git *)"]
    assert _tool_in_allowlist("Bash", cli) is True
    assert _tool_in_allowlist("Read", cli) is True
    assert _tool_in_allowlist("Edit", cli) is False


def test_is_structurally_blocked_dashboard_edit():
    """Real-world case: dashboard mode lacks Edit, agent tries to Edit."""
    dashboard_tools = [
        "Read", "Grep", "Glob", "Write", "Bash", "WebSearch", "WebFetch",
    ]
    raw = "{'tool_name': 'Edit', 'tool_input': {'file_path': 'foo.py'}}"
    assert _is_structurally_blocked(raw, dashboard_tools) is True


def test_is_structurally_blocked_solo_edit_is_allowed():
    solo_tools = [
        "Read", "Grep", "Glob", "Bash", "Edit", "Write",
        "WebSearch", "WebFetch",
    ]
    raw = "{'tool_name': 'Edit', 'tool_input': {'file_path': 'foo.py'}}"
    assert _is_structurally_blocked(raw, solo_tools) is False


def test_is_structurally_blocked_returns_false_when_unrestricted():
    """``allowed_tools=None`` means CLI default — never structurally blocked."""
    raw = "{'tool_name': 'Edit', 'tool_input': {}}"
    assert _is_structurally_blocked(raw, None) is False


def test_is_structurally_blocked_returns_false_when_tool_unknown():
    """Empty tool name → can't determine, treat as not structurally blocked."""
    assert _is_structurally_blocked("", ["Read"]) is False
    assert _is_structurally_blocked("garbage", ["Read"]) is False


# ---------------------------------------------------------------------------
# Step 3 — full-transcript handoff for live mode-switch
# ---------------------------------------------------------------------------

def test_handoff_empty_chat_returns_empty():
    assert _build_full_transcript_handoff([], "dashboard", "solo") == ""


def test_handoff_chat_with_only_empty_content_returns_empty():
    """No real content → no handoff (avoids empty header-only blob)."""
    msgs = [
        {"role": "user", "content": ""},
        {"role": "assistant", "content": None},
    ]
    assert _build_full_transcript_handoff(msgs, "dashboard", "solo") == ""


def test_handoff_includes_mode_pair_in_header():
    msgs = [{"role": "user", "content": "hi"}]
    out = _build_full_transcript_handoff(msgs, "dashboard", "solo")
    assert "[Mode switch: dashboard → solo]" in out
    assert "operating as the **solo** agent" in out


def test_handoff_renders_all_known_roles():
    msgs = [
        {"role": "user", "content": "fix bug"},
        {"role": "assistant", "content": "located the issue"},
        {"role": "system", "content": "Cycle complete"},
        {"role": "tool", "content": "diff output", "tool_name": "Edit"},
    ]
    out = _build_full_transcript_handoff(msgs, "dashboard", "solo")
    assert "### User\nfix bug" in out
    assert "### Assistant\nlocated the issue" in out
    assert "### System\nCycle complete" in out
    assert "### Tool (Edit)\ndiff output" in out


def test_handoff_unknown_role_rendered_as_note():
    """Unrecognised role keys are surfaced under "Note", never silently dropped."""
    msgs = [{"role": "guest", "content": "external note"}]
    out = _build_full_transcript_handoff(msgs, "a", "b")
    assert "### Note (guest)\nexternal note" in out


def test_handoff_truncates_long_messages():
    long_text = "x" * 12000
    msgs = [{"role": "user", "content": long_text}]
    out = _build_full_transcript_handoff(
        msgs, "dashboard", "solo", max_chars_per_msg=4000
    )
    # Cap respected (header + "[truncated ... chars]" footer)
    assert "[truncated" in out
    assert "x" * 5000 not in out


def test_handoff_preserves_message_order():
    msgs = [
        {"role": "user", "content": "first"},
        {"role": "assistant", "content": "reply1"},
        {"role": "user", "content": "second"},
        {"role": "assistant", "content": "reply2"},
    ]
    out = _build_full_transcript_handoff(msgs, "x", "y")
    i_first = out.find("first")
    i_reply1 = out.find("reply1")
    i_second = out.find("second")
    i_reply2 = out.find("reply2")
    assert 0 < i_first < i_reply1 < i_second < i_reply2


def test_handoff_skips_messages_with_empty_content():
    """Empty-content messages are silently dropped from the transcript."""
    msgs = [
        {"role": "user", "content": "real"},
        {"role": "assistant", "content": ""},
        {"role": "user", "content": None},
    ]
    out = _build_full_transcript_handoff(msgs, "x", "y")
    assert "### User\nreal" in out
    assert "### Assistant" not in out  # empty assistant skipped


def test_handoff_includes_context_framing():
    """The framing tells the new agent it's a continuation, not a fresh start."""
    msgs = [{"role": "user", "content": "hi"}]
    out = _build_full_transcript_handoff(msgs, "dashboard", "solo")
    assert "Prior conversation" in out
    assert "End of prior context" in out


# ---------------------------------------------------------------------------
# S4 — Session boot context (once per fresh dashboard session)
# ---------------------------------------------------------------------------

def test_session_boot_empty_inputs_returns_empty():
    """No data → no header — caller must be free to skip the block entirely."""
    assert _format_session_boot() == ""
    assert _format_session_boot(outcomes=[], jobs=[], commits=[]) == ""


def test_session_boot_renders_outcome_block():
    outcomes = [
        {"verdict": "PASS", "mode": "solo", "task": "fix bug", "cost_usd": 0.012},
        {"verdict": "FAIL", "mode": "dashboard", "task": "recalc job_X", "cost_usd": 0.234},
    ]
    out = _format_session_boot(outcomes=outcomes)
    assert "[Session boot context]" in out
    assert "Recent outcomes" in out
    assert "[PASS]" in out and "[FAIL]" in out
    assert "fix bug" in out
    assert "$0.234" in out


def test_session_boot_outcomes_truncated_to_last_5():
    outcomes = [
        {"verdict": "PASS", "mode": "x", "task": f"t{i}", "cost_usd": 0.01}
        for i in range(20)
    ]
    out = _format_session_boot(outcomes=outcomes)
    # Most recent 5 are t15..t19 — t14 must NOT show
    assert "'t19'" in out
    assert "'t14'" not in out


def test_session_boot_renders_job_block():
    jobs = [
        {"name": "Cas_red_1", "status": "RUNNING"},
        {"name": "job_2", "status": "RUNNING"},
        {"name": "job_3", "status": "PENDING"},
    ]
    out = _format_session_boot(jobs=jobs)
    assert "Active SLURM jobs" in out
    assert "RUNNING:2" in out
    assert "PENDING:1" in out
    assert "Cas_red_1(R)" in out


def test_session_boot_renders_commits_and_branch():
    out = _format_session_boot(
        commits=["abc1234 fix scf", "def5678 add NEB", "999aaaa upd docs"],
        branch="GUPPY",
    )
    assert "Recent commits:" in out
    assert "abc1234" in out
    assert "Active branch: GUPPY" in out


def test_session_boot_renders_calc_path():
    out = _format_session_boot(calc_path="/home/u/calc/Jerome/c1")
    assert "Active calc folder: /home/u/calc/Jerome/c1" in out


def test_session_boot_truncates_above_max_chars():
    """Cap respected even with absurd outcome list."""
    outcomes = [
        {"verdict": "PASS", "mode": "x", "task": "x" * 200, "cost_usd": 0.0}
        for _ in range(50)
    ]
    out = _format_session_boot(outcomes=outcomes, max_chars=400)
    assert "[truncated]" in out
    assert len(out) <= 600  # header + 400 + "[truncated]" footer


def test_session_boot_handles_dataclass_like_objects():
    """Outcomes can be CycleOutcome dataclasses (attrs) instead of dicts."""
    class _Fake:
        verdict = "PASS"
        mode = "solo"
        task = "real task"
        cost_usd = 0.05
    out = _format_session_boot(outcomes=[_Fake(), _Fake()])
    assert "real task" in out and "[PASS]" in out


def test_session_boot_skips_missing_blocks():
    """A boot context with only commits must NOT leak empty section headers."""
    out = _format_session_boot(commits=["abc1234 fix"])
    assert "Recent commits:" in out
    assert "Recent outcomes" not in out
    assert "Active SLURM jobs" not in out
    assert "Active calc folder" not in out


# ---------------------------------------------------------------------------
# B2 — _compute_turn_verdict (real PASS/FAIL/PARTIAL detection)
# ---------------------------------------------------------------------------

class _FakeEngine:
    """Minimal stand-in for AgentEngine — only the attrs the verdict uses."""
    def __init__(self, *, stop_requested: bool = False):
        self._stop_requested = stop_requested


def test_verdict_pass_default():
    """A clean response without denied commands or errors → PASS."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="All set. Done.",
        denied_commands=[],
        state={},
    )
    assert v == "PASS"


def test_verdict_fail_on_explicit_error_type():
    """Engine bubbled up an exception → FAIL regardless of text content."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="Successfully did the thing.",  # text claims success
        denied_commands=[],
        state={},
        error_type="exception",
    )
    assert v == "FAIL"


def test_verdict_fail_on_stop_requested_with_short_response():
    """User hit ESC and almost nothing came out → FAIL."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(stop_requested=True),
        response_text="Working on it",  # < 80 chars
        denied_commands=[],
        state={},
    )
    assert v == "FAIL"


def test_verdict_pass_on_stop_requested_with_long_response():
    """Stop requested AFTER a substantial answer is not a fail."""
    long = "x" * 200
    v = _compute_turn_verdict(
        engine=_FakeEngine(stop_requested=True),
        response_text=long,
        denied_commands=[],
        state={},
    )
    assert v == "PASS"


def test_verdict_fail_on_denied_plus_giveup_phrase():
    """Tool calls blocked + agent says it cannot continue → FAIL."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="I'm unable to complete this without the Edit tool.",
        denied_commands=["Edit"],
        state={},
    )
    assert v == "FAIL"


def test_verdict_partial_on_question_tag():
    """Agent paused with QUESTION: → PARTIAL (waiting on user)."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="Did some work. QUESTION: which file should I touch next?",
        denied_commands=[],
        state={},
    )
    assert v == "PARTIAL"


def test_verdict_partial_on_incomplete_phrase():
    """Agent self-reports incomplete work → PARTIAL."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="Couldn't fully complete the task — maxiter ran out.",
        denied_commands=[],
        state={},
    )
    assert v == "PARTIAL"


def test_verdict_partial_when_denied_but_progress_made():
    """Tool calls blocked but the agent still produced a substantive reply
    (no give-up phrase) → PARTIAL — work continued under restrictions."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="Skipped the Edit; here's the patch you can apply manually.",
        denied_commands=["Edit"],
        state={},
    )
    assert v == "PARTIAL"


def test_verdict_records_german_giveup_phrase():
    """The fail-phrase list covers German fallbacks too."""
    v = _compute_turn_verdict(
        engine=_FakeEngine(),
        response_text="Ich konnte nicht weitermachen, alle Edits sind blockiert.",
        denied_commands=["Edit"],
        state={},
    )
    assert v == "FAIL"


# ---------------------------------------------------------------------------
# B1 — _record_turn_outcome (renamed, tracks all modes)
# ---------------------------------------------------------------------------

def test_record_turn_outcome_skips_empty_task():
    """An empty user task is a no-op (no outcome write)."""

    class _NoCallEngine:
        mode = "solo"
        def record_cycle_outcome(self, **kw):
            raise AssertionError("must NOT be called for empty task")

    out = _record_turn_outcome(
        _NoCallEngine(),
        user_task="",
        response_text="some reply",
        state={},
        start_time=None,
    )
    assert out == {}


def test_record_turn_outcome_skips_empty_response():
    """An empty assistant response is also a no-op."""

    class _NoCallEngine:
        mode = "solo"
        def record_cycle_outcome(self, **kw):
            raise AssertionError("must NOT be called for empty response")

    out = _record_turn_outcome(
        _NoCallEngine(),
        user_task="real task",
        response_text="",
        state={},
        start_time=None,
    )
    assert out == {}


def test_record_turn_outcome_passes_computed_verdict():
    """The verdict produced by _compute_turn_verdict reaches the engine call."""
    captured: dict = {}

    class _Engine:
        mode = "dashboard"
        _stop_requested = False
        def record_cycle_outcome(self, **kw):
            captured.update(kw)
            return {}

    _record_turn_outcome(
        _Engine(),
        user_task="fix bug",
        response_text="I'm unable to do that, Edit is blocked.",
        state={"_denied_commands": ["Edit"]},
        start_time=None,
    )
    assert captured["verdict"] == "FAIL"
    assert captured["denied_commands"] == ["Edit"]


def test_record_turn_outcome_runs_for_all_modes():
    """B1 — dashboard / quick / reviewed all reach record_cycle_outcome."""
    captured: list[str] = []

    class _Engine:
        _stop_requested = False
        def __init__(self, mode):
            self.mode = mode
        def record_cycle_outcome(self, **kw):
            captured.append(self.mode)
            return {}

    for mode in ("solo", "dashboard", "quick", "reviewed", "tdd"):
        _record_turn_outcome(
            _Engine(mode),
            user_task="t",
            response_text="ok done",
            state={},
            start_time=None,
        )
    assert captured == ["solo", "dashboard", "quick", "reviewed", "tdd"]


def test_record_solo_turn_outcome_alias_still_works():
    """The old name is kept as an alias so external callers don't break."""
    assert _record_solo_turn_outcome is _record_turn_outcome


# ---------------------------------------------------------------------------
# Plan-accept arming heuristic + header wait chip
# ---------------------------------------------------------------------------


def test_greeting_does_not_look_like_a_plan():
    from delfin.dashboard.tab_agent import _looks_like_plan_response
    assert not _looks_like_plan_response("Hallo! Ich bin bereit zu helfen.")
    assert not _looks_like_plan_response("")
    assert not _looks_like_plan_response("Kurze Antwort ohne Struktur." * 3)


def test_numbered_plan_is_recognised():
    from delfin.dashboard.tab_agent import _looks_like_plan_response
    plan = (
        "Hier ist mein Plan für die Umsetzung des Features:\n"
        "1. Zuerst die bestehende Struktur analysieren und die "
        "betroffenen Module identifizieren.\n"
        "2. Danach die neue Funktion implementieren, mit Tests "
        "für jeden Randfall.\n"
        "3. Anschließend die Dokumentation aktualisieren.\n"
        "4. Zum Schluss den vollen Testlauf ausführen.\n"
    )
    assert _looks_like_plan_response(plan)


def test_bulleted_plan_is_recognised():
    from delfin.dashboard.tab_agent import _looks_like_plan_response
    plan = (
        "Vorgehen für die Migration der Datenbankschicht im Detail:\n"
        "- Backup der bestehenden Daten anlegen und verifizieren\n"
        "- Schema-Migration mit Rollback-Pfad implementieren\n"
        "- Integrationstests gegen die migrierte Datenbank laufen lassen\n"
    )
    assert _looks_like_plan_response(plan)


def test_short_list_is_not_a_plan():
    from delfin.dashboard.tab_agent import _looks_like_plan_response
    assert not _looks_like_plan_response("1. a\n2. b\n")


def test_wait_chip_renders_and_clears():
    from delfin.dashboard.tab_agent import _wait_chip_html
    assert _wait_chip_html("") == ""
    # The label is passed IN, so the test should assert the chip renders
    # whatever it was given rather than pinning one particular label.
    out = _wait_chip_html("Plan approval")
    assert "waiting:" in out and "Plan approval" in out
    long = _wait_chip_html("x" * 100)
    assert "…" in long
    # HTML in the wait text must be escaped, not rendered.
    assert "<script>" not in _wait_chip_html("<script>alert(1)</script>")


# ---------------------------------------------------------------------------
# Honest restore: a "pending" plan that already executed must say so
# ---------------------------------------------------------------------------


def _trace(*tools):
    return [{"tool": t} for t in tools]


def test_post_plan_executions_counted_from_field_case():
    """Shape of archive case 20260729-112242: plan submitted, then the
    execution turn ran writes/bash before the host died."""
    from delfin.dashboard.tab_agent import _count_post_plan_executions
    entries = _trace(
        "mcp__kit-coding__project_introspect",
        "mcp__kit-coding__task_create",
        "mcp__kit-coding__exit_plan_mode",
        "mcp__kit-coding__task_update",
        "mcp__kit-coding__bash",
        "mcp__kit-coding__write_file",
        "mcp__kit-coding__bash_background",
        "mcp__kit-coding__bash_kill",
    )
    assert _count_post_plan_executions(entries) == 4


def test_no_execution_after_plan_submission_counts_zero():
    from delfin.dashboard.tab_agent import _count_post_plan_executions
    assert _count_post_plan_executions(_trace(
        "mcp__kit-coding__task_create",
        "mcp__kit-coding__exit_plan_mode",
        "mcp__kit-coding__task_update",
    )) == 0
    # Writes BEFORE the submission belong to an earlier phase.
    assert _count_post_plan_executions(_trace(
        "mcp__kit-coding__write_file",
        "mcp__kit-coding__exit_plan_mode",
    )) == 0
    assert _count_post_plan_executions([]) == 0
    assert _count_post_plan_executions([{}, {"name": ""}]) == 0


def test_resubmitted_plan_resets_the_counter():
    from delfin.dashboard.tab_agent import _count_post_plan_executions
    assert _count_post_plan_executions(_trace(
        "mcp__kit-coding__exit_plan_mode",
        "mcp__kit-coding__bash",
        "mcp__kit-coding__exit_plan_mode",
    )) == 0


# ---------------------------------------------------------------------------
# Stale watchdog: prefill is not a stall
# ---------------------------------------------------------------------------


def _watchdog_source() -> str:
    from pathlib import Path as _P
    return (_P(__file__).resolve().parent.parent / "delfin" / "dashboard"
            / "tab_agent.py").read_text(encoding="utf-8")


def test_first_token_budget_is_separate_and_larger():
    """Field case 20260729-122058: a dashboard turn was killed after 121 s
    of silence although the provider had simply not started yet (measured
    ~96 s time-to-first-token on the same endpoint for a one-word turn).
    Waiting for the FIRST token must use its own, much larger budget."""
    src = _watchdog_source()
    assert "first_token_kill_after_s" in src
    assert "waiting_for_first = not state.get(\"_stream_saw_output\")" in src
    assert "budget = first_token_kill if waiting_for_first else kill_after" in src
    assert "max(\n            600.0, kill_after * 4.0)" in src


def test_stream_output_marks_first_token_for_every_channel():
    src = _watchdog_source()
    # text, thinking and tool use all count as "the provider started".
    assert src.count('state["_stream_saw_output"] = True') >= 3
    # ... and the flag is reset when a turn starts.
    assert 'state["_stream_saw_output"] = False' in src


def test_watchdog_kill_is_not_reported_as_a_backend_failure():
    src = _watchdog_source()
    assert 'state["_watchdog_stopped"] = (' in src
    assert 'if not state.pop("_watchdog_stopped", "")' in src
    assert "Turn ended by DELFIN's watchdog" in src
    # The old wording blamed the CLI for our own stop.
    assert "The CLI ended the turn" not in src


def test_kill_watch_rearms_instead_of_expiring():
    """The first-token budget outlives the initial timer, so a not-yet-due
    check must re-arm rather than silently stop watching."""
    src = _watchdog_source()
    idx = src.find("def _check_kill():")
    assert idx > 0
    block = src[idx:idx + 2000]
    assert "again = _threading.Timer(" in block
    assert 'state["_stale_kill_timer"] = again' in block


# ---------------------------------------------------------------------------
# Mid-turn persistence + plan-mode verify exemption
# ---------------------------------------------------------------------------


def test_session_is_checkpointed_before_and_during_a_turn():
    """Field complaint: after a reload only the plan was there, not the
    conversation that produced it — auto-save ran only at turn end."""
    src = _watchdog_source()
    assert "def _checkpoint_session(" in src
    # before the turn runs, on every finalized message, and before the
    # plan-approval block parks the worker
    assert src.count("_checkpoint_session(") >= 4
    assert "_checkpoint_session(min_interval_s=0.0)" in src


def test_checkpoint_is_throttled_and_never_raises():
    src = _watchdog_source()
    idx = src.find("def _checkpoint_session(")
    block = src[idx:idx + 900]
    assert "min_interval_s" in block
    assert "_last_checkpoint_ts" in block
    assert "except Exception:" in block


def test_plan_turns_are_exempt_from_claim_grounding():
    """A plan names files it INTENDS to create; grounding them against the
    workspace made every plan a false alarm (field case 20260729-125618:
    3 flags on tetris_game.py / snake_game.py / game_dashboard.ipynb)."""
    from pathlib import Path as _P
    eng = (_P(__file__).resolve().parent.parent / "delfin" / "agent"
           / "engine.py").read_text(encoding="utf-8")
    assert "def _turn_describes_intent(" in eng
    assert 'getattr(self.kit_permissions, "mode", "") or "") == "plan"' in eng
    assert "exit_plan_mode" in eng
    assert "and not self._turn_describes_intent()" in eng
    src = _watchdog_source()
    assert "_describes_intent = engine._turn_describes_intent()" in src
    assert "and not _describes_intent" in src


# ---------------------------------------------------------------------------
# ACTION continuation budget — progress vs. repetition
# ---------------------------------------------------------------------------

def _run_rounds(rounds, ceiling=12, repeat_limit=2):
    """Drive the pure decision over a list of per-round ACTION sets.

    Mirrors the bookkeeping the continuation loop does (signature history,
    executed-command set, round counter) and returns
    ``(rounds_completed, stop_reason)``.
    """
    from delfin.dashboard.tab_agent import (
        _decide_action_round, _normalize_action_command)
    sigs: list[str] = []
    seen: set[str] = set()
    used = 0
    stale = 0
    for cmds in rounds:
        used += 1
        dec = _decide_action_round(
            cmds, sigs, seen, rounds_used=used, stale_rounds=stale,
            ceiling=ceiling, repeat_limit=repeat_limit)
        sigs.append(dec.signature)
        if dec.stale:
            stale += 1
        for c in cmds:
            n = _normalize_action_command(c)
            if n:
                seen.add(n)
        if not dec.proceed:
            return used, dec.reason
    return used, ""


def test_distinct_action_sets_keep_going_up_to_the_ceiling():
    """Genuine multi-step work — a different command every round — must
    not be cut short at the old flat cap of 3. The only thing that ends
    such a sequence is the absolute ceiling."""
    rounds = [[f"/tab step{i}"] for i in range(20)]
    used, reason = _run_rounds(rounds, ceiling=12)
    assert reason == "ceiling"
    assert used == 12


def test_four_distinct_rounds_are_progress_not_a_loop():
    """The field case: read the error file, open the folder, inspect the
    input, answer. Four distinct rounds must all be allowed."""
    used, reason = _run_rounds([
        ["/tab calc"],
        ["/calc show job.err"],
        ["/calc show job.inp"],
        ["/tab orca"],
    ], ceiling=12)
    assert (used, reason) == (4, "")


def test_same_action_set_twice_stops():
    """Two occurrences of the identical set is a loop — stop immediately,
    do not bleed to the ceiling."""
    used, reason = _run_rounds([
        ["/tab orca"],
        ["/tab orca"],
        ["/tab orca"],
    ], ceiling=12)
    assert reason == "repeat"
    assert used == 2


def test_repeat_detection_ignores_order_and_case_and_spacing():
    """Same work in a different order (or with cosmetic differences) is
    still the same work."""
    used, reason = _run_rounds([
        ["/tab orca", "/jobs"],
        ["/JOBS", "/tab   orca"],
    ], ceiling=12)
    assert (used, reason) == (2, "repeat")


def test_alternating_two_sets_stops_on_the_first_repeat():
    """A/B/A alternation dodged a previous-round-only comparison; the full
    per-turn signature history catches it."""
    used, reason = _run_rounds([
        ["/tab orca"],
        ["/jobs"],
        ["/tab orca"],
    ], ceiling=12)
    assert (used, reason) == (3, "repeat")


def test_empty_round_stops():
    """A round that executed no ACTION has nothing to feed back — another
    model turn would be paid for nothing."""
    from delfin.dashboard.tab_agent import _decide_action_round
    dec = _decide_action_round([], [], set(), rounds_used=1)
    assert dec.proceed is False
    assert dec.reason == "no_actions"
    assert dec.new_commands == ()


def test_round_of_only_already_executed_commands_is_no_progress():
    """A reshuffled/partial re-emission is not an identical set, but it
    still adds nothing new."""
    from delfin.dashboard.tab_agent import _decide_action_round
    dec = _decide_action_round(
        ["/tab orca"], ["/jobs\x1f/tab orca"], {"/tab orca", "/jobs"},
        rounds_used=2)
    assert dec.proceed is False
    assert dec.reason == "no_progress"


def test_one_new_command_among_old_ones_is_progress():
    from delfin.dashboard.tab_agent import _decide_action_round
    dec = _decide_action_round(
        ["/tab orca", "/calc show new.err"], ["/tab orca"], {"/tab orca"},
        rounds_used=1)
    assert dec.proceed is True
    assert dec.reason == "progress"
    assert dec.new_commands == ("/calc show new.err",)


def test_repeat_limit_is_clamped_to_at_least_two():
    """A limit of 1 (or 0) would stop on the very first round and make the
    loop useless — the helper clamps it."""
    from delfin.dashboard.tab_agent import _decide_action_round
    dec = _decide_action_round(
        ["/tab orca"], [], set(), rounds_used=1, repeat_limit=1)
    assert dec.proceed is True
    dec = _decide_action_round(
        ["/tab orca"], [], set(), rounds_used=1, repeat_limit=0)
    assert dec.proceed is True


def test_higher_repeat_limit_tolerates_one_more_occurrence():
    """The setting has to actually do something — with a limit of 3 the
    same set may occur three times before the turn ends."""
    used, reason = _run_rounds([
        ["/tab orca"], ["/tab orca"], ["/tab orca"],
    ], ceiling=12, repeat_limit=3)
    assert (used, reason) == (3, "repeat")


def test_action_round_signature_is_order_insensitive():
    from delfin.dashboard.tab_agent import _action_round_signature
    assert (_action_round_signature(["/a", "/b"])
            == _action_round_signature(["/b", "/a"]))
    assert _action_round_signature([]) == ""
    assert _action_round_signature(["/a"]) != _action_round_signature(["/b"])


def test_decide_action_round_survives_garbage_limits():
    """Settings come from a user-editable file — bad values must not raise
    inside the streaming worker."""
    from delfin.dashboard.tab_agent import _decide_action_round
    dec = _decide_action_round(
        ["/tab orca"], [], set(), rounds_used="x",
        ceiling="nope", repeat_limit=None)
    assert dec.proceed is True


def test_resolve_action_round_limits_defaults(monkeypatch):
    """No setting → generous ceiling, repeat limit 2."""
    from delfin import user_settings
    from delfin.dashboard import tab_agent as ta
    monkeypatch.setattr(user_settings, "load_settings", lambda: {})
    assert ta._resolve_action_round_limits() == (
        ta._ACTION_ROUND_CEILING_DEFAULT, 2)


def test_resolve_action_round_limits_reads_agent_settings(monkeypatch):
    from delfin import user_settings
    from delfin.dashboard import tab_agent as ta
    monkeypatch.setattr(
        user_settings, "load_settings",
        lambda: {"agent": {"max_action_rounds": 5,
                           "action_repeat_limit": 4}})
    assert ta._resolve_action_round_limits() == (5, 4)


def test_resolve_action_round_limits_zero_disables_the_ceiling(monkeypatch):
    from delfin import user_settings
    from delfin.dashboard import tab_agent as ta
    monkeypatch.setattr(
        user_settings, "load_settings",
        lambda: {"agent": {"max_action_rounds": 0}})
    ceiling, repeat = ta._resolve_action_round_limits()
    assert ceiling >= 10_000
    assert repeat == 2


def test_resolve_action_round_limits_never_raises(monkeypatch):
    from delfin import user_settings
    from delfin.dashboard import tab_agent as ta

    def _boom():
        raise OSError("settings file corrupt")

    monkeypatch.setattr(user_settings, "load_settings", _boom)
    assert ta._resolve_action_round_limits() == (
        ta._ACTION_ROUND_CEILING_DEFAULT, 2)
    monkeypatch.setattr(
        user_settings, "load_settings",
        lambda: {"agent": {"max_action_rounds": "many",
                           "action_repeat_limit": "lots"}})
    assert ta._resolve_action_round_limits() == (
        ta._ACTION_ROUND_CEILING_DEFAULT, 2)


# -- the stop note ----------------------------------------------------------

def test_ceiling_stop_note_is_honest_and_actionable():
    """It must name what ran, say the ROUND LIMIT ended the turn (not that
    the work is done or that the agent misbehaved), and say how to go on."""
    from delfin.dashboard.tab_agent import _format_action_stop_note
    note = _format_action_stop_note(
        "ceiling", ["/tab calc", "/calc show job.err"], 12, [])
    assert "/tab calc" in note and "/calc show job.err" in note
    assert "12" in note
    assert "agent.max_action_rounds" in note
    assert "continue" in note.lower()
    # Honesty: no claim that the agent misbehaved, no claim of completion
    assert "kept emitting" not in note
    assert "completed" not in note.lower()


def test_ceiling_stop_note_lists_unexecuted_commands():
    from delfin.dashboard.tab_agent import _format_action_stop_note
    note = _format_action_stop_note(
        "ceiling", ["/tab calc"], 12, ["/jobs", "/tab orca"])
    assert "Not executed" in note
    assert "/jobs" in note and "/tab orca" in note


def test_repeat_stop_note_says_the_agent_repeated_itself():
    from delfin.dashboard.tab_agent import _format_action_stop_note
    note = _format_action_stop_note("repeat", ["/tab orca"], 12, [])
    assert "/tab orca" in note
    assert "already run" in note
    assert "follow-up" in note
    for reason in ("repeat", "no_progress"):
        assert _format_action_stop_note(reason, ["/tab orca"], 12, [])


def test_stop_note_is_empty_for_a_turn_that_finished_on_its_own():
    from delfin.dashboard.tab_agent import _format_action_stop_note
    assert _format_action_stop_note("progress", ["/tab orca"], 12, []) == ""
    assert _format_action_stop_note("no_actions", [], 12, []) == ""


def test_stop_note_truncates_long_command_lists():
    from delfin.dashboard.tab_agent import _format_action_stop_note
    note = _format_action_stop_note(
        "ceiling", [f"/tab t{i}" for i in range(9)], 12, [])
    assert "+4 more" in note


# ---------------------------------------------------------------------------
# An exhausted account budget is not a model problem
# ---------------------------------------------------------------------------


def test_exceeded_budget_is_classified_as_quota():
    """Field case 2026-07-30: the gateway answered 'Error code: 400 -
    ExceededBudget ... Spend=5.55, Budget=5.0'. The agent read it as a
    per-call failure and kept improvising against a hard wall."""
    from delfin.dashboard.tab_agent import _classify_model_error_text as cls
    for text in (
        "Error code: 400 - {'detail': 'ExceededBudget: End User=x@kit.edu "
        "over budget. Spend=5.548, Budget=5.0'}",
        "insufficient_quota: you exceeded your current quota",
        "quota exceeded for this project",
    ):
        assert cls(text) == "quota", text


def test_quota_wins_over_the_other_classes():
    """A budget message must never be read as a flapping backend or a bad
    key — those invite a retry or a model switch, which cannot help."""
    from delfin.dashboard.tab_agent import _classify_model_error_text as cls
    assert cls("ExceededBudget ... service unavailable, model group=x") == "quota"


def test_ordinary_errors_are_not_quota():
    from delfin.dashboard.tab_agent import _classify_model_error_text as cls
    assert cls("Error code: 503 - temporarily unavailable") == "temp"
    assert cls("AuthenticationError: invalid subscription") == "auth"
    assert cls("Error code: 400 - maximum context length exceeded") == ""


def test_spend_amounts_are_extracted_for_the_message():
    from delfin.dashboard.tab_agent import _QUOTA_SPEND_RE
    m = _QUOTA_SPEND_RE.search("over budget. Spend=5.548838, Budget=5.0")
    assert m and m.group(1).startswith("5.54") and m.group(2) == "5.0"


def test_quota_message_states_that_switching_models_will_not_help():
    from pathlib import Path as _P
    src = (_P(__file__).resolve().parent.parent / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    assert "is_quota_exhausted" in src
    assert "cap sits on the ACCOUNT" in src
    assert "does not help either" in src
