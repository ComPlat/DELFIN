"""Remote Archive tab: browse a configured remote archive over SSH."""

from __future__ import annotations

import base64
import csv as _csv
import html
import io
import json
import posixpath
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path, PurePosixPath

import ipywidgets as widgets
from IPython.display import HTML, clear_output, display

from delfin.remote_archive import (
    TEXT_PREVIEW_MAX_BYTES,
    fetch_remote_file,
    list_remote_entries,
    normalize_remote_relative_path,
    read_remote_text_preview,
    remote_delete,
    remote_duplicate,
    remote_mkdir,
    remote_rename,
)
from delfin.ssh_transfer_jobs import (
    create_download_job,
    ensure_jobs_dir,
    launch_transfer_job,
    list_transfer_jobs,
)
from delfin.user_settings import load_transfer_settings

from .constants import CALC_SEARCH_OPTIONS

try:
    from delfin.remote_archive import list_remote_folder_files as _list_remote_folder_files
except ImportError:
    _list_remote_folder_files = None

from .helpers import disable_spellcheck
from .input_processing import is_smiles, smiles_to_xyz_quick
from .molecule_viewer import (
    DEFAULT_3DMOL_STYLE_JS,
    DEFAULT_3DMOL_ZOOM,
    coord_to_xyz,
    parse_xyz_frames,
    patch_viewer_mouse_controls_js,
)

REMOTE_FULL_FETCH_MAX_BYTES = 128 * 1024 * 1024
REMOTE_TEXT_RENDER_MAX_CHARS = 400_000
REMOTE_MOL_SIZE = 450
REMOTE_MOL_DYNAMIC_SCALE = 0.9725
REMOTE_LEFT_DEFAULT = 375
REMOTE_LEFT_MIN = 375
REMOTE_LEFT_MAX = 520
REMOTE_XYZ_LARGE_TRAJ_FRAMES = 2000
REMOTE_XYZ_PLAY_FPS_DEFAULT = 10
REMOTE_XYZ_PLAY_FPS_MIN = 1
REMOTE_XYZ_PLAY_FPS_MAX = 60
REMOTE_SEARCH_MAX_MATCHES = 2000
REMOTE_HIGHLIGHT_MAX_CHARS = 400_000
REMOTE_DOWNLOAD_MAX_BYTES = 25 * 1024 * 1024


def create_tab(ctx):
    """Create the Remote Archive browser tab."""
    state = {
        "config": None,
        "current_relative_path": "",
        "entries": [],
        "filtered_entries": [],
        "selected_entry": None,
        "selected_remote_path": "",
        "file_content": "",
        "file_preview_note": "",
        "visualize_kind": "",
        "visualize_payload": "",
        "visualize_enabled": False,
        "current_xyz_frames": [],
        "current_xyz_index": 0,
        "traj_playing": False,
        "traj_play_toggle_guard": False,
        "traj_viewer_ready": False,
        # Search state
        "search_spans": [],
        "current_match": -1,
        "search_truncated": False,
        "selected_file_path": None,
        "selected_file_size": 0,
        # Table extraction state
        "table_col_defs": [
            {"name": "Value 1", "type": "regex", "pattern": "", "occ": "last"},
        ],
        "table_col_widgets": [],
        "table_csv_data": "",
        "table_panel_active": False,
    }
    scope_id = f"remote-archive-scope-{abs(id(state))}"
    remote_resize_mol_fn = f"remoteArchiveResizeMolViewer_{abs(id(state))}"
    viewer_mouse_patch_js = patch_viewer_mouse_controls_js("viewer", "el")
    mol3d_counter = [0]

    title = widgets.HTML("<h3>📂 Remote Archive</h3>")
    info_html = widgets.HTML(value="")
    status_html = widgets.HTML(value="")
    path_html = widgets.HTML(value="<b>📂 Path:</b> /")

    up_btn = widgets.Button(
        description="⬆ Up",
        button_style="warning",
        layout=widgets.Layout(width="62px", height="28px"),
    )
    home_btn = widgets.Button(
        description="🏠",
        button_style="info",
        layout=widgets.Layout(width="62px", height="28px"),
    )
    refresh_btn = widgets.Button(
        description="🔄",
        layout=widgets.Layout(width="62px", height="28px"),
    )
    open_btn = widgets.Button(
        description="Open",
        button_style="primary",
        layout=widgets.Layout(width="88px", height="28px"),
    )
    filter_input = widgets.Text(
        placeholder="Filter remote files...",
        layout=widgets.Layout(flex="1 1 auto", min_width="160px", height="28px"),
    )
    filter_input.add_class("remote-archive-filter")
    sort_dropdown = widgets.Dropdown(
        options=[("A-Z", "name"), ("Newest", "date_desc"), ("Oldest", "date_asc")],
        value="name",
        layout=widgets.Layout(width="98px", height="28px"),
    )
    file_list = widgets.Select(
        options=[],
        rows=22,
        layout=widgets.Layout(width="100%", flex="1 1 0", min_height="0", margin="-4px 0 0 0"),
    )
    file_list.add_class("remote-file-list")
    keyboard_action_input = widgets.Text(
        value="",
        layout=widgets.Layout(width="1px", height="1px", display="none"),
    )
    keyboard_action_input.add_class("remote-cmd-keyboard-action")
    file_info_html = widgets.HTML(value="")
    selected_path_html = widgets.HTML(value="", layout=widgets.Layout(display="none"))
    transfer_back_btn = widgets.Button(
        description="→ Calculations",
        button_style="info",
        layout=widgets.Layout(width="126px", min_width="126px", height="26px"),
        disabled=True,
    )
    transfer_to_archive_btn = widgets.Button(
        description="→ Archive",
        button_style="info",
        layout=widgets.Layout(width="104px", min_width="104px", height="26px"),
        disabled=True,
    )
    transfer_jobs_btn = widgets.Button(
        description="Running Transfers",
        layout=widgets.Layout(width="142px", min_width="142px", height="26px"),
    )
    transfer_jobs_refresh_btn = widgets.Button(
        description="Refresh Jobs",
        layout=widgets.Layout(width="112px", min_width="112px", height="26px"),
    )
    copy_path_btn = widgets.Button(
        description="PATH",
        layout=widgets.Layout(width="70px", min_width="70px", height="26px"),
        disabled=True,
    )
    copy_btn = widgets.Button(
        description="Copy",
        button_style="info",
        layout=widgets.Layout(width="80px", min_width="80px", height="26px"),
        disabled=True,
    )
    view_toggle = widgets.ToggleButton(
        description="Visualize",
        value=False,
        disabled=True,
        button_style="warning",
        layout=widgets.Layout(width="110px", min_width="110px", height="26px"),
    )

    # -- File operation buttons -------------------------------------------------
    new_folder_btn = widgets.Button(
        description="\U0001f4c1 New Folder",
        layout=widgets.Layout(width="110px", min_width="110px", height="26px"),
        disabled=True,
    )
    rename_btn = widgets.Button(
        description="\u270f Rename",
        layout=widgets.Layout(width="90px", min_width="90px", height="26px"),
        disabled=True,
    )
    duplicate_btn = widgets.Button(
        description="\U0001f4cb Duplicate",
        layout=widgets.Layout(width="106px", min_width="106px", height="26px"),
        disabled=True,
    )
    delete_btn = widgets.Button(
        description="\U0001f5d1 Delete",
        button_style="danger",
        layout=widgets.Layout(width="90px", min_width="90px", height="26px"),
        disabled=True,
    )
    # Confirmation/input dialogs
    ops_status_html = widgets.HTML(value="")
    rename_input = widgets.Text(
        placeholder="New name...",
        layout=widgets.Layout(width="200px", height="28px"),
    )
    new_folder_input = widgets.Text(
        placeholder="Folder name...",
        layout=widgets.Layout(width="200px", height="28px"),
    )
    confirm_panel = widgets.HBox(
        [], layout=widgets.Layout(display="none", gap="6px", align_items="center", width="100%"),
    )

    # -- Search widgets ---------------------------------------------------------
    search_input = widgets.Text(
        placeholder="Search in file...",
        continuous_update=True,
        layout=widgets.Layout(width="140px", min_width="140px", height="26px"),
    )
    search_input.add_class("remote-search-input")
    search_suggest = widgets.Dropdown(
        options=["(Select)"] + CALC_SEARCH_OPTIONS,
        value="(Select)",
        layout=widgets.Layout(width="200px", min_width="200px", height="26px"),
    )
    search_btn = widgets.Button(
        description="🔍",
        layout=widgets.Layout(width="85px", min_width="85px", height="26px"),
    )
    search_prev_btn = widgets.Button(
        description="◀",
        layout=widgets.Layout(width="58px", min_width="58px", height="26px"),
        disabled=True,
    )
    search_next_btn = widgets.Button(
        description="▶",
        layout=widgets.Layout(width="58px", min_width="58px", height="26px"),
        disabled=True,
    )
    search_result_html = widgets.HTML(
        value="",
        layout=widgets.Layout(width="180px", min_width="180px"),
    )

    # -- Content navigation widgets ---------------------------------------------
    top_btn = widgets.Button(
        description="⬆ Top",
        layout=widgets.Layout(width="78px", min_width="78px", height="26px"),
    )
    bottom_btn = widgets.Button(
        description="⬇ End",
        layout=widgets.Layout(width="78px", min_width="78px", height="26px"),
    )

    # -- Download widgets -------------------------------------------------------
    download_btn = widgets.Button(
        description="Download",
        layout=widgets.Layout(width="100px", min_width="100px", height="26px"),
        disabled=True,
    )
    download_status_html = widgets.HTML(
        value="",
        layout=widgets.Layout(width="100%", overflow_x="hidden"),
    )

    # -- Editable path input ----------------------------------------------------
    path_input = widgets.Text(
        value="/",
        continuous_update=False,
        layout=widgets.Layout(
            flex="1 1 0", min_width="0", width="1px", max_width="100%",
            height="24px", overflow_x="hidden", margin="0", padding="0",
        ),
    )
    path_input.add_class("remote-path-input")
    path_prefix_html = widgets.HTML(
        value="<b>📂 Path:</b>",
        layout=widgets.Layout(width="100%"),
    )
    path_input_box = widgets.VBox(
        [path_prefix_html, path_input],
        layout=widgets.Layout(
            width="100%", overflow_x="hidden",
            align_items="stretch", gap="2px",
        ),
    )

    # -- Table presets (same as Calculations Browser) ---------------------------
    _TABLE_PRESETS = [
        ('— Preset —',                '',         'regex', '',                                                         'last'),
        ('ORCA Total Energy (Eh)',    'Energy',   'regex', r'FINAL SINGLE POINT ENERGY\s+(-?[\d.]+)',                 'last'),
        ('ORCA SCF Energy (Eh)',      'SCF E',    'regex', r'E\(SCF\)\s*=\s*(-?[\d.]+)',                             'last'),
        ('ORCA HOMO (eV)',            'HOMO',     'regex', r'(?i)homo.*?(-?[\d.]+)\s*eV',                            'last'),
        ('ORCA LUMO (eV)',            'LUMO',     'regex', r'(?i)lumo.*?(-?[\d.]+)\s*eV',                            'last'),
        ('ORCA Dipole total (D)',     'Dipole',   'regex', r'(?i)total dipole moment.*?([\d.]+)',                     'last'),
        ('ORCA S\u00b2 (before)',     'S2',       'regex', r'<S\*\*2> Expectation value\s+([\d.]+)',                 'last'),
        ('ORCA Net Charge',           'Charge',   'regex', r'Total Charge\s*\.*\s*(-?\d+)',                           'first'),
        ('ORCA Nr. Atoms',            'Atoms',    'regex', r'Number of atoms\s*\.\.\.\s*(\d+)',                       'first'),
        ('ORCA Nr. Basis Fns',        'Basis',    'regex', r'Number of basis functions\s*\.\.\.\s*(\d+)',             'first'),
        ('TM SCF Energy (Eh)',        'E(TM)',    'regex', r'SCF energy\s+:\s+(-?[\d.]+)',                           'last'),
        ('TM HOMO (eV)',              'HOMO(TM)', 'regex', r'HOMO.*?(-?[\d.]+)\s*eV',                                'last'),
        ('DELFIN TDDFT to_state',     'to_state', 'json',  'tddft_absorption.transitions[*].to_state',                'all'),
        ('DELFIN TDDFT energy_eV',    'energy_eV','json',  'tddft_absorption.transitions[*].energy_eV',               'all'),
        ('DELFIN TDDFT wavelength',   'nm',       'json',  'tddft_absorption.transitions[*].wavelength_nm',           'all'),
        ('DELFIN TDDFT fosc',         'fosc',     'json',  'tddft_absorption.transitions[*].oscillator_strength',     'all'),
        ('DELFIN S0→S1 wavelength',   'S0S1_nm',  'json',  'tddft_absorption.transitions[?(@.from_state=="S0" && @.to_state=="S1")].wavelength_nm', 'first'),
        ('DELFIN S0→S1 energy_eV',    'S0S1_eV',  'json',  'tddft_absorption.transitions[?(@.from_state=="S0" && @.to_state=="S1")].energy_eV', 'first'),
        ('JSON key path',             'JSON',     'json',  'results.energy.final',                                     'last'),
    ]
    _tp_labels = [p[0] for p in _TABLE_PRESETS]
    _tp_names  = [p[1] for p in _TABLE_PRESETS]
    _tp_types  = [p[2] for p in _TABLE_PRESETS]
    _tp_pats   = [p[3] for p in _TABLE_PRESETS]
    _tp_occs   = [p[4] for p in _TABLE_PRESETS]
    _table_number_re = re.compile(
        r'[-+]?(?:(?:\d+[.,]\d*)|(?:[.,]\d+)|\d+)(?:[eEdD][-+]?\d+)?'
    )
    _table_number_full_re = re.compile(
        r'^[-+]?(?:(?:\d+[.,]\d*)|(?:[.,]\d+)|\d+)(?:[eEdD][-+]?\d+)?$'
    )

    # -- Table extraction widgets -----------------------------------------------
    table_btn = widgets.Button(
        description="📊 Table",
        layout=widgets.Layout(width="84px", height="26px"),
        disabled=(_list_remote_folder_files is None),
    )
    table_file_input = widgets.Text(
        placeholder="e.g. orca.out or result.json",
        layout=widgets.Layout(flex="1 1 auto", min_width="120px", height="26px"),
    )
    table_scope_dd = widgets.Dropdown(
        options=[("All folders", "all")],
        value="all",
        layout=widgets.Layout(width="130px", height="26px"),
    )
    table_recursive_cb = widgets.ToggleButton(
        value=False,
        description="Recurse",
        layout=widgets.Layout(width="86px", height="26px"),
    )
    table_decimal_comma_btn = widgets.ToggleButton(
        value=False,
        description="Dot to Comma",
        tooltip="Convert decimal point to decimal comma in table/CSV values",
        layout=widgets.Layout(width="126px", height="26px"),
    )
    table_add_col_btn = widgets.Button(
        description="+ Column",
        layout=widgets.Layout(width="80px", height="26px"),
    )
    table_run_btn = widgets.Button(
        description="▶ Run",
        button_style="primary",
        layout=widgets.Layout(width="68px", height="26px"),
    )
    table_close_btn = widgets.Button(
        description="✕",
        layout=widgets.Layout(width="32px", height="26px"),
    )
    table_csv_btn = widgets.Button(
        description="⬇ CSV",
        layout=widgets.Layout(width="68px", height="26px", display="none"),
    )
    table_status_html = widgets.HTML(value="")
    table_output_html = widgets.HTML(
        value="",
        layout=widgets.Layout(
            width="100%", overflow_x="auto", overflow_y="auto",
            flex="1 1 0", min_height="0",
        ),
    )
    table_cols_box = widgets.VBox(
        [],
        layout=widgets.Layout(width="100%", gap="4px"),
    )
    table_panel = widgets.VBox(
        [
            widgets.HBox(
                [widgets.HTML("<b>📊 Extract Table</b>"), table_close_btn],
                layout=widgets.Layout(
                    justify_content="space-between", align_items="center", width="100%",
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<span style="white-space:nowrap"><b>File:</b></span>'),
                    table_file_input,
                    table_scope_dd,
                    table_recursive_cb,
                    table_decimal_comma_btn,
                ],
                layout=widgets.Layout(gap="6px", align_items="center", width="100%"),
            ),
            table_cols_box,
            widgets.HBox(
                [table_add_col_btn, table_run_btn, table_csv_btn],
                layout=widgets.Layout(gap="6px", align_items="center"),
            ),
            table_status_html,
            table_output_html,
        ],
        layout=widgets.Layout(
            display="none", width="100%", padding="8px", gap="6px",
            border="1px solid #e0e0e0", border_radius="4px", overflow_x="hidden",
            flex="1 1 0", min_height="0",
        ),
    )

    content_label = widgets.HTML(
        "<div style='height:26px; line-height:26px; margin:0 0 8px 0;'>"
        "<b>📄 File Content:</b></div>"
    )
    viewer_label = widgets.HTML(
        "<div style='height:26px; line-height:26px; margin:0 0 8px 0;'>"
        "<b>🔬 Molecule Preview:</b></div>"
    )
    frame_label_html = widgets.HTML(value="", layout=widgets.Layout(display="none"))
    frame_input = widgets.BoundedIntText(
        value=1,
        min=1,
        max=1,
        step=1,
        layout=widgets.Layout(width="74px", height="28px", display="none"),
    )
    frame_input.add_class("remote-xyz-frame-input")
    frame_total_html = widgets.HTML(value="", layout=widgets.Layout(display="none"))
    xyz_loop_checkbox = widgets.ToggleButton(
        value=True,
        description="Loop",
        button_style="info",
        layout=widgets.Layout(width="82px", height="32px"),
    )
    xyz_fps_input = widgets.BoundedIntText(
        value=REMOTE_XYZ_PLAY_FPS_DEFAULT,
        min=REMOTE_XYZ_PLAY_FPS_MIN,
        max=REMOTE_XYZ_PLAY_FPS_MAX,
        step=1,
        layout=widgets.Layout(width="72px", height="28px"),
    )
    xyz_play_btn = widgets.ToggleButton(
        value=False,
        description="Play",
        icon="play",
        button_style="success",
        layout=widgets.Layout(width="86px", height="32px"),
    )
    xyz_play_btn.add_class("remote-xyz-play-btn")
    xyz_copy_btn = widgets.Button(
        description="📋 Copy Coordinates",
        button_style="success",
        layout=widgets.Layout(width="176px", min_width="176px", height="32px"),
    )
    viewer_output = widgets.Output(
        layout=widgets.Layout(
            width="100%",
            border="2px solid #1976d2",
            min_height="0",
            overflow="hidden",
        )
    )
    viewer_output.add_class("remote-mol-viewer")
    preview_html = widgets.HTML(
        value="",
        layout=widgets.Layout(width="100%", flex="1 1 0", min_height="0", overflow="auto"),
    )
    transfer_jobs_html = widgets.HTML(
        value="",
        layout=widgets.Layout(width="100%", overflow_x="hidden", overflow_y="auto", flex="1 1 0", min_height="0"),
    )
    transfer_jobs_panel = widgets.VBox(
        [
            widgets.HBox(
                [widgets.HTML("<b>Running Transfers</b>"), transfer_jobs_refresh_btn],
                layout=widgets.Layout(gap="6px", align_items="center", flex_flow="row wrap"),
            ),
            transfer_jobs_html,
        ],
        layout=widgets.Layout(
            display="none",
            gap="6px",
            width="100%",
            flex="1 1 0",
            min_height="0",
            overflow="hidden",
        ),
    )

    def _set_status(message, color="#455a64"):
        status_html.value = f'<span style="color:{color};">{message}</span>'

    def _run_js(script):
        ctx.run_js(script)

    def _copy_to_clipboard(text, label="content"):
        text_payload = json.dumps(str(text or ""))
        label_payload = json.dumps(str(label or "content"))
        _run_js(
            "(function(){"
            f"const text={text_payload};"
            f"const label={label_payload};"
            "function _manualPrompt(){"
            "try{window.prompt('Copy to clipboard (Cmd+C/Ctrl+C, Enter):', text);}catch(_e){}"
            "}"
            "function _legacyCopy(){"
            "try{"
            "const ta=document.createElement('textarea');"
            "ta.value=text;"
            "ta.setAttribute('readonly','readonly');"
            "ta.style.position='fixed';"
            "ta.style.top='-1000px';"
            "ta.style.left='-1000px';"
            "ta.style.opacity='0';"
            "document.body.appendChild(ta);"
            "ta.focus();"
            "ta.select();"
            "ta.setSelectionRange(0, ta.value.length);"
            "const ok=document.execCommand('copy');"
            "document.body.removeChild(ta);"
            "return !!ok;"
            "}catch(_e){return false;}"
            "}"
            "if(navigator.clipboard && navigator.clipboard.writeText){"
            "navigator.clipboard.writeText(text).catch(function(){"
            "if(!_legacyCopy()) _manualPrompt();"
            "});"
            "}else{"
            "if(!_legacyCopy()) _manualPrompt();"
            "}"
            "})();"
        )

    def _traj_can_play():
        frame_count = len(state.get("current_xyz_frames") or [])
        return 1 < frame_count <= REMOTE_XYZ_LARGE_TRAJ_FRAMES

    def _update_loop_button_style():
        xyz_loop_checkbox.button_style = "info" if xyz_loop_checkbox.value else ""

    def _set_play_button_state(active, sync_value=True):
        active = bool(active)
        state["traj_playing"] = active
        xyz_play_btn.description = "Stop" if active else "Play"
        xyz_play_btn.icon = "stop" if active else "play"
        xyz_play_btn.button_style = "danger" if active else "success"
        if sync_value and xyz_play_btn.value != active:
            state["traj_play_toggle_guard"] = True
            try:
                xyz_play_btn.value = active
            finally:
                state["traj_play_toggle_guard"] = False

    def _stop_xyz_playback(update_button=True):
        scope_key_json = json.dumps(scope_id)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                if (window._remoteTrajPlayTimerByScope && window._remoteTrajPlayTimerByScope[scopeKey]) {{
                    clearInterval(window._remoteTrajPlayTimerByScope[scopeKey]);
                    delete window._remoteTrajPlayTimerByScope[scopeKey];
                }}
            }})();
            """
        )
        if update_button:
            _set_play_button_state(False, sync_value=True)

    def _update_traj_control_state():
        frames = state.get("current_xyz_frames") or []
        has_xyz = bool(frames) and view_toggle.value and state.get("visualize_kind") == "xyz"
        xyz_controls.layout.display = "flex" if has_xyz else "none"
        xyz_playback_row.layout.display = "flex" if has_xyz else "none"
        xyz_tray_controls.layout.display = "flex" if has_xyz else "none"
        can_play = has_xyz and _traj_can_play()
        xyz_loop_checkbox.disabled = not can_play
        xyz_fps_input.disabled = not can_play
        xyz_play_btn.disabled = not can_play
        xyz_copy_btn.disabled = not has_xyz
        _update_loop_button_style()
        if not can_play:
            _stop_xyz_playback(update_button=True)

    def _set_view_toggle(value, disabled=None):
        try:
            view_toggle.unobserve(_on_view_toggle, names="value")
        except Exception:
            pass
        if disabled is not None:
            view_toggle.disabled = disabled
        view_toggle.value = value
        view_toggle.observe(_on_view_toggle, names="value")

    def _transfer_jobs_dir():
        return ensure_jobs_dir()

    def _transfer_status_color(status):
        mapping = {
            "queued": "#1976d2",
            "running": "#1976d2",
            "retrying": "#ef6c00",
            "success": "#2e7d32",
            "warning": "#ef6c00",
            "failed": "#d32f2f",
        }
        return mapping.get(str(status or "").lower(), "#555")

    def _transfer_items_html(entry, limit=4):
        raw_sources = entry.get("sources", []) or []
        rendered = []
        direction = str(entry.get("direction") or "push").lower()
        for raw_source in raw_sources[:limit]:
            source_text = str(raw_source or "").strip()
            if not source_text:
                continue
            label = source_text
            if direction == "push":
                source_path = Path(source_text)
                label = source_path.name or source_text
            else:
                label = PurePosixPath(source_text).name or source_text
            rendered.append(f'<code title="{html.escape(source_text)}">{html.escape(label)}</code>')
        if not rendered:
            return "n/a"
        remaining = len(raw_sources) - len(rendered)
        summary = ", ".join(rendered)
        if remaining > 0:
            summary = f"{summary} + {remaining} more"
        return summary

    def _render_transfer_jobs(limit=8):
        try:
            entries = list_transfer_jobs(jobs_dir=_transfer_jobs_dir(), limit=limit)
        except Exception as exc:
            transfer_jobs_html.value = (
                f'<span style="color:#d32f2f;">Could not load transfer jobs: '
                f'{html.escape(str(exc))}</span>'
            )
            return
        if not entries:
            transfer_jobs_html.value = '<span style="color:#555;">No background transfer jobs yet.</span>'
            return

        blocks = []
        for entry in entries:
            status = str(entry.get("status") or "unknown")
            color = _transfer_status_color(status)
            direction = str(entry.get("direction") or "push").lower()
            remote = (
                f'{entry.get("user", "?")}@{entry.get("host", "?")}:'
                f'{entry.get("remote_path", "?")}'
            )
            local_target = str(entry.get("local_target") or "").strip()
            endpoint = remote if direction == "push" else f"{remote} -> {local_target or '?'}"
            summary = str(entry.get("last_summary") or entry.get("last_error") or "").strip()
            retry_note = ""
            retry_in = entry.get("retry_in_seconds", 0) or 0
            if status == "retrying" and retry_in:
                retry_note = f" Retry in about {int(retry_in)} s."
            attempts = entry.get("attempt", 0)
            max_retries = entry.get("max_retries", 0)
            log_path = entry.get("log_path", "")
            items_html = _transfer_items_html(entry)
            updated = str(entry.get("updated_at") or "").replace("T", " ")
            if updated.endswith("+00:00"):
                updated = updated[:-6] + " UTC"
            blocks.append(
                '<div style="border:1px solid #d9dee3;border-radius:6px;padding:8px 10px;'
                'margin:0 0 8px 0;background:#fafbfc;">'
                f'<div><b>{html.escape(entry.get("job_id", "job"))}</b> '
                f'<span style="color:{color};font-weight:600;">{html.escape(status.upper())}</span></div>'
                f'<div><code>{html.escape(endpoint)}</code></div>'
                f'<div>{int(entry.get("source_count", 0) or 0)} item(s), '
                f'attempt {int(attempts or 0)}/{int(max_retries or 0) + 1}</div>'
                f'<div>Items: {items_html}</div>'
                f'<div>{html.escape(summary or "No status message.")}{html.escape(retry_note)}</div>'
                f'<div style="color:#555;">Updated: {html.escape(updated or "-")}</div>'
                f'<div style="color:#555;">Log: <code>{html.escape(str(log_path))}</code></div>'
                '</div>'
            )
        transfer_jobs_html.value = "".join(blocks)

    def _update_transfer_jobs_visibility():
        jobs_visible = transfer_jobs_panel.layout.display != "none"
        transfer_jobs_btn.button_style = "primary" if jobs_visible else ""
        if jobs_visible:
            left_panel.add_class("remote-transfer-jobs-mode")
            filter_row.layout.display = "none"
            file_list.layout.display = "none"
        else:
            left_panel.remove_class("remote-transfer-jobs-mode")
            filter_row.layout.display = "flex"
            file_list.layout.display = ""

    def _settings_summary(config):
        if not config:
            return "No remote archive configured."
        remote_root = str(config.get("remote_path") or "/")
        return (
            f'<b>Remote target:</b> '
            f'<code>{html.escape(str(config.get("user") or "?"))}@'
            f'{html.escape(str(config.get("host") or "?"))}:{html.escape(remote_root)}</code> '
            f'<span style="color:#616161;">(browse inside this root; selected items can be moved back into local Calculations or Archive)</span>'
        )

    def _format_size(size_bytes):
        size = int(size_bytes or 0)
        if size >= 1024 ** 3:
            return f"{size / (1024 ** 3):.2f} GB"
        if size >= 1024 ** 2:
            return f"{size / (1024 ** 2):.2f} MB"
        if size >= 1024:
            return f"{size / 1024:.2f} KB"
        return f"{size} B"

    def _entry_icon(entry):
        if entry.get("is_dir"):
            return "📂"
        suffix = str(entry.get("suffix") or "").lower()
        name = str(entry.get("name") or "")
        if name.lower() == "coord":
            return "🔬"
        if suffix == ".xyz":
            return "🔬"
        if suffix == ".png":
            return "🖼"
        if suffix in {".inp", ".sh"}:
            return "📝"
        if suffix in {".cube", ".cub"}:
            return "🧊"
        if suffix in {".gbw", ".cis", ".densities"}:
            return "💾"
        if suffix in {".doc", ".docx"}:
            return "📃"
        return "📄"

    def _entry_label(entry):
        return f'{_entry_icon(entry)} {entry.get("name", "")}'

    def _entry_title(entry):
        if entry.get("is_dir"):
            return entry.get("name", "")
        return f'{entry.get("name", "")} ({_format_size(entry.get("size", 0))})'

    def _current_remote_path(config=None):
        cfg = config or state.get("config") or {}
        root_path = str(cfg.get("remote_path") or "/").rstrip("/") or "/"
        rel = normalize_remote_relative_path(state.get("current_relative_path", ""))
        return root_path if not rel else f"{root_path}/{rel}"

    def _entry_by_relative_path(relative_path):
        rel = normalize_remote_relative_path(relative_path)
        for entry in state.get("filtered_entries", []):
            if normalize_remote_relative_path(entry.get("relative_path", "")) == rel:
                return entry
        return None

    _path_syncing = [False]

    def _update_path_html():
        config = state.get("config")
        _path_syncing[0] = True
        try:
            if not config:
                path_input.value = "/"
            else:
                path_input.value = _current_remote_path(config)
        finally:
            _path_syncing[0] = False

    def _set_viewer_visible(is_visible):
        viewer_container.layout.display = "flex" if is_visible else "none"

    def _set_selected_path_display(path_value):
        selected_path_html.value = ""

    def _clear_viewer():
        viewer_output.clear_output()
        _set_viewer_visible(False)
        state["traj_viewer_ready"] = False
        scope_key_json = json.dumps(scope_id)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                if (window._remoteMolViewerByScope) delete window._remoteMolViewerByScope[scopeKey];
                if (window._remoteTrajViewerByScope) delete window._remoteTrajViewerByScope[scopeKey];
            }})();
            """
        )

    def _hide_frame_controls():
        state["current_xyz_frames"] = []
        state["current_xyz_index"] = 0
        state["traj_viewer_ready"] = False
        frame_label_html.layout.display = "none"
        frame_input.layout.display = "none"
        frame_total_html.layout.display = "none"
        frame_label_html.value = ""
        frame_total_html.value = ""
        frame_input.max = 1
        frame_input.value = 1
        _update_traj_control_state()

    def _reset_visualization_state():
        state["visualize_kind"] = ""
        state["visualize_payload"] = ""
        state["visualize_enabled"] = False
        _set_view_toggle(False, disabled=True)
        _clear_viewer()
        _hide_frame_controls()

    def _set_visualization(kind="", payload="", frames=None):
        state["visualize_kind"] = str(kind or "")
        state["visualize_payload"] = str(payload or "")
        state["visualize_enabled"] = bool(kind)
        state["current_xyz_frames"] = list(frames or [])
        state["current_xyz_index"] = 0
        state["traj_viewer_ready"] = False
        if state["visualize_enabled"]:
            _set_view_toggle(bool(view_toggle.value), disabled=False)
        else:
            _set_view_toggle(False, disabled=True)
            _clear_viewer()
        _update_traj_control_state()

    def _clear_preview(message="Select a remote file to preview."):
        file_info_html.value = ""
        selected_path_html.value = ""
        state["file_content"] = ""
        state["file_preview_note"] = ""
        state["selected_file_path"] = None
        state["selected_file_size"] = 0
        preview_html.value = (
            "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
            "padding:12px; background:#fafafa;'>"
            f"{html.escape(message)}</div>"
        )
        copy_btn.disabled = True
        copy_path_btn.disabled = True
        download_btn.disabled = True
        download_status_html.value = ""
        content_toolbar.layout.display = "none"
        _reset_search_state()
        _reset_visualization_state()

    def _render_text_preview(text, *, note=""):
        content = str(text or "")
        truncated = False
        if len(content) > REMOTE_TEXT_RENDER_MAX_CHARS:
            content = content[:REMOTE_TEXT_RENDER_MAX_CHARS]
            truncated = True
        note_parts = []
        if note:
            note_parts.append(note)
        if truncated:
            note_parts.append(f"display limited to {REMOTE_TEXT_RENDER_MAX_CHARS:,} characters")
        note_html = ""
        if note_parts:
            note_html = (
                "<div style='color:#616161; margin:0 0 8px 0;'>"
                + " | ".join(html.escape(part) for part in note_parts)
                + "</div>"
            )
        preview_html.value = (
            "<style>"
            ".remote-match { background: #fff59d; padding: 0 2px; }"
            ".remote-match.current { background: #ffcc80; }"
            "</style>"
            "<div id='remote-content-box' style='height:100%; overflow-y:auto; overflow-x:hidden;"
            " border:1px solid #ddd; padding:6px; background:#fafafa; width:100%; box-sizing:border-box;'>"
            f"{note_html}"
            "<div id='remote-content-text' style='white-space:pre-wrap; overflow-wrap:anywhere;"
            " word-break:break-word; font-family:monospace; font-size:12px; line-height:1.3;'>"
            f"{html.escape(content)}"
            "</div></div>"
        )

    def _render_image_preview(local_path):
        data = Path(local_path).read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        preview_html.value = (
            "<div style='height:100%; width:100%; border:1px solid #ddd; padding:8px; "
            "background:#fafafa; box-sizing:border-box; display:flex; align-items:center; "
            "justify-content:center; overflow:auto;'>"
            f"<img src='data:image/png;base64,{b64}' style='max-width:100%; max-height:100%; "
            "object-fit:contain; display:block;' /></div>"
        )

    def _render_docx_preview(local_path):
        try:
            import mammoth
        except ImportError:
            _render_text_preview("mammoth not installed.\n\nRun: pip install mammoth")
            return

        def _convert_image_to_base64(image):
            try:
                with image.open() as handle:
                    payload = handle.read()
                content_type = image.content_type or "image/png"
                b64 = base64.b64encode(payload).decode("ascii")
                return {"src": f"data:{content_type};base64,{b64}"}
            except Exception:
                return {"src": ""}

        with open(str(local_path), "rb") as docx_handle:
            result = mammoth.convert_to_html(
                docx_handle,
                convert_image=mammoth.images.img_element(_convert_image_to_base64),
            )
        preview_html.value = (
            "<div style='height:100%; overflow:auto; border:1px solid #ddd; background:white; "
            "padding:14px; box-sizing:border-box;'>"
            "<style>"
            ".remote-docx img { max-width:100%; height:auto; margin:10px 0; }"
            ".remote-docx table { border-collapse:collapse; width:100%; margin:10px 0; }"
            ".remote-docx td, .remote-docx th { border:1px solid #ddd; padding:8px; }"
            "</style>"
            f"<div class='remote-docx'>{result.value}</div>"
            "</div>"
        )

    def _render_3dmol(data, fmt="xyz", volumetric=False):
        mol3d_counter[0] += 1
        viewer_id = f"remote_mol3d_{mol3d_counter[0]}"
        wrapper_id = f"remote_mol_wrap_{mol3d_counter[0]}"
        data_json = json.dumps(str(data or ""))
        scope_key_json = json.dumps(scope_id)
        view_scope_json = json.dumps(f"{scope_id}:{state.get('current_relative_path') or '/'}")
        volumetric_js = ""
        if volumetric:
            volumetric_js = (
                "viewer.addVolumetricData(molData,'cube',{isoval:0.02,color:'#0026ff',opacity:0.85});"
                "viewer.addVolumetricData(molData,'cube',{isoval:-0.02,color:'#b00010',opacity:0.85});"
            )
        with viewer_output:
            clear_output()
            display(
                HTML(
                    f"""
                    <div id="{wrapper_id}" class="remote-mol-stage-wrapper" style="width:100%;">
                        <div id="{viewer_id}" style="width:100%;height:{REMOTE_MOL_SIZE}px;position:relative;"></div>
                    </div>
                    <script>
                    if (typeof $3Dmol === "undefined") {{
                        var _s = document.createElement("script");
                        _s.src = "https://3Dmol.org/build/3Dmol-min.js";
                        document.head.appendChild(_s);
                    }}
                    (function() {{
                        var tries = 0;
                        function initViewer() {{
                            var el = document.getElementById("{viewer_id}");
                            var mv = el ? el.closest('.remote-mol-viewer') : null;
                            var scopeRoot = el ? el.closest('.{scope_id}') : null;
                            if (!scopeRoot) scopeRoot = document.querySelector('.{scope_id}');
                            if (!el || typeof $3Dmol === "undefined" || !mv || mv.offsetParent === null) {{
                                tries += 1;
                                if (tries < 80) setTimeout(initViewer, 50);
                                return;
                            }}
                            var rightPanel = scopeRoot ? scopeRoot.querySelector('.remote-right') : null;
                            if (rightPanel) {{
                                var mvRect = mv.getBoundingClientRect();
                                if (mvRect.top > 0 || mvRect.height > 0) {{
                                    var rightRect = rightPanel.getBoundingClientRect();
                                    var topChildren = Array.prototype.slice.call(rightPanel.children || []);
                                    var host = null;
                                    for (var i = 0; i < topChildren.length; i++) {{
                                        if (topChildren[i].contains(mv)) {{
                                            host = topChildren[i];
                                            break;
                                        }}
                                    }}
                                    var reservedBelow = 0;
                                    if (host) {{
                                        var passed = false;
                                        for (var j = 0; j < topChildren.length; j++) {{
                                            var child = topChildren[j];
                                            if (child === host) {{
                                                passed = true;
                                                continue;
                                            }}
                                            if (!passed) continue;
                                            var style = window.getComputedStyle(child);
                                            if (!style || style.display === 'none' || style.visibility === 'hidden') continue;
                                            var childRect = child.getBoundingClientRect();
                                            if (childRect.height > 0) reservedBelow += childRect.height;
                                        }}
                                    }}
                                    var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
                                    var row = mv.closest('.remote-mol-view-row');
                                    var rowRect = row ? row.getBoundingClientRect() : rightRect;
                                    var tray = scopeRoot ? scopeRoot.querySelector('.remote-xyz-tray-controls') : null;
                                    var trayStyle = tray ? window.getComputedStyle(tray) : null;
                                    var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
                                    var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
                                    var availW = Math.max(120, rowRect.width - trayWidth - 16);
                                    var h = Math.floor(availH * {REMOTE_MOL_DYNAMIC_SCALE});
                                    var w = Math.floor(Math.min(h * 1.2, availW));
                                    if (h >= 80 && w >= 120) {{
                                        mv.style.width = w + 'px';
                                        mv.style.height = h + 'px';
                                        el.style.width = w + 'px';
                                        el.style.height = h + 'px';
                                    }}
                                }}
                            }}
                            window._remoteMolViewStateByScope = window._remoteMolViewStateByScope || {{}};
                            window._remoteMolViewerByScope = window._remoteMolViewerByScope || {{}};
                            window._remoteTrajViewerByScope = window._remoteTrajViewerByScope || {{}};
                            window._remoteMolViewScopeKeyByScope = window._remoteMolViewScopeKeyByScope || {{}};
                            var scopeKey = {scope_key_json};
                            var viewScope = {view_scope_json};
                            var previousViewer =
                                window._remoteMolViewerByScope[scopeKey]
                                || window._remoteTrajViewerByScope[scopeKey]
                                || null;
                            var previousScope =
                                window._remoteMolViewScopeKeyByScope[scopeKey] || viewScope;
                            if (previousViewer && typeof previousViewer.getView === 'function') {{
                                try {{
                                    window._remoteMolViewStateByScope[previousScope] = previousViewer.getView();
                                }} catch (_e) {{}}
                            }}
                            var savedView = window._remoteMolViewStateByScope[viewScope] || null;
                            var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                            {viewer_mouse_patch_js}
                            var molData = {data_json};
                            viewer.addModel(molData, "{fmt}");
                            viewer.setStyle({{}}, {DEFAULT_3DMOL_STYLE_JS});
                            {volumetric_js}
                            if (savedView && typeof viewer.setView === 'function') {{
                                try {{
                                    viewer.setView(savedView);
                                }} catch (_e) {{
                                    viewer.zoomTo();
                                    viewer.center();
                                    viewer.zoom({DEFAULT_3DMOL_ZOOM});
                                }}
                            }} else {{
                                viewer.zoomTo();
                                viewer.center();
                                viewer.zoom({DEFAULT_3DMOL_ZOOM});
                            }}
                            viewer.render();
                            window._remoteMolViewerByScope[scopeKey] = viewer;
                            window._remoteMolViewScopeKeyByScope[scopeKey] = viewScope;
                            var scopeRoot2 = document.querySelector('.{scope_id}');
                            var wrappers = scopeRoot2
                                ? scopeRoot2.querySelectorAll('.remote-mol-stage-wrapper')
                                : document.querySelectorAll('.remote-mol-stage-wrapper');
                            wrappers.forEach(function(w) {{
                                if (w.id !== "{wrapper_id}") w.remove();
                            }});
                            if (window["{remote_resize_mol_fn}"]) {{
                                setTimeout(window["{remote_resize_mol_fn}"], 200);
                            }}
                        }}
                        setTimeout(initViewer, 0);
                    }})();
                    </script>
                    """
                )
            )
        _set_viewer_visible(True)

    def _render_xyz_in_viewer(xyz_text):
        _render_3dmol(xyz_text, fmt="xyz", volumetric=False)

    def _render_xyz_trajectory_viewer(initial_load=False):
        frames = state.get("current_xyz_frames") or []
        if not frames:
            return
        idx = max(0, min(len(frames) - 1, int(state.get("current_xyz_index", 0))))
        state["current_xyz_index"] = idx
        if len(frames) > REMOTE_XYZ_LARGE_TRAJ_FRAMES:
            _render_xyz_in_viewer(_frame_to_xyz(frames[idx]))
            return
        if not initial_load and state.get("traj_viewer_ready"):
            scope_key_json = json.dumps(scope_id)
            _run_js(
                f"""
                setTimeout(function(){{
                    var scopeKey = {scope_key_json};
                    var viewer = window._remoteTrajViewerByScope
                        ? window._remoteTrajViewerByScope[scopeKey]
                        : null;
                    if (viewer) {{
                        viewer.setFrame({idx});
                        viewer.render();
                    }}
                }}, 0);
                """
            )
            return

        full_xyz = "".join(_frame_to_xyz(frame) for frame in frames)
        mol3d_counter[0] += 1
        viewer_id = f"remote_trj_viewer_{mol3d_counter[0]}"
        wrapper_id = f"remote_mol_wrap_{mol3d_counter[0]}"
        scope_key_json = json.dumps(scope_id)
        with viewer_output:
            clear_output()
            display(
                HTML(
                    f"""
                    <div id="{wrapper_id}" class="remote-mol-stage-wrapper" style="width:100%;">
                        <div id="{viewer_id}" style="width:100%;height:{REMOTE_MOL_SIZE}px;position:relative;"></div>
                    </div>
                    <script>
                    if (typeof $3Dmol === "undefined") {{
                        var _s = document.createElement("script");
                        _s.src = "https://3Dmol.org/build/3Dmol-min.js";
                        document.head.appendChild(_s);
                    }}
                    (function() {{
                        var tries = 0;
                        function initViewer() {{
                            var el = document.getElementById("{viewer_id}");
                            var mv = el ? el.closest('.remote-mol-viewer') : null;
                            var scopeRoot = el ? el.closest('.{scope_id}') : null;
                            if (!scopeRoot) scopeRoot = document.querySelector('.{scope_id}');
                            if (!el || typeof $3Dmol === "undefined" || !mv || mv.offsetParent === null) {{
                                tries += 1;
                                if (tries < 80) setTimeout(initViewer, 50);
                                return;
                            }}
                            var rightPanel = scopeRoot ? scopeRoot.querySelector('.remote-right') : null;
                            if (rightPanel) {{
                                var mvRect = mv.getBoundingClientRect();
                                if (mvRect.top > 0 || mvRect.height > 0) {{
                                    var rightRect = rightPanel.getBoundingClientRect();
                                    var topChildren = Array.prototype.slice.call(rightPanel.children || []);
                                    var host = null;
                                    for (var i = 0; i < topChildren.length; i++) {{
                                        if (topChildren[i].contains(mv)) {{
                                            host = topChildren[i];
                                            break;
                                        }}
                                    }}
                                    var reservedBelow = 0;
                                    if (host) {{
                                        var passed = false;
                                        for (var j = 0; j < topChildren.length; j++) {{
                                            var child = topChildren[j];
                                            if (child === host) {{
                                                passed = true;
                                                continue;
                                            }}
                                            if (!passed) continue;
                                            var style = window.getComputedStyle(child);
                                            if (!style || style.display === 'none' || style.visibility === 'hidden') continue;
                                            var childRect = child.getBoundingClientRect();
                                            if (childRect.height > 0) reservedBelow += childRect.height;
                                        }}
                                    }}
                                    var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
                                    var row = mv.closest('.remote-mol-view-row');
                                    var rowRect = row ? row.getBoundingClientRect() : rightRect;
                                    var tray = scopeRoot ? scopeRoot.querySelector('.remote-xyz-tray-controls') : null;
                                    var trayStyle = tray ? window.getComputedStyle(tray) : null;
                                    var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
                                    var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
                                    var availW = Math.max(120, rowRect.width - trayWidth - 16);
                                    var h = Math.floor(availH * {REMOTE_MOL_DYNAMIC_SCALE});
                                    var w = Math.floor(Math.min(h * 1.2, availW));
                                    if (h >= 80 && w >= 120) {{
                                        mv.style.width = w + 'px';
                                        mv.style.height = h + 'px';
                                        el.style.width = w + 'px';
                                        el.style.height = h + 'px';
                                    }}
                                }}
                            }}
                            var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                            {viewer_mouse_patch_js}
                            viewer.addModelsAsFrames(`{full_xyz}`, "xyz");
                            viewer.setStyle({{}}, {DEFAULT_3DMOL_STYLE_JS});
                            viewer.zoomTo();
                            viewer.center();
                            viewer.zoom({DEFAULT_3DMOL_ZOOM});
                            viewer.setFrame({idx});
                            viewer.render();
                            window._remoteTrajViewerByScope = window._remoteTrajViewerByScope || {{}};
                            window._remoteMolViewerByScope = window._remoteMolViewerByScope || {{}};
                            window._remoteTrajViewerByScope[{scope_key_json}] = viewer;
                            window._remoteMolViewerByScope[{scope_key_json}] = viewer;
                            var scopeRoot2 = document.querySelector('.{scope_id}');
                            var wrappers = scopeRoot2
                                ? scopeRoot2.querySelectorAll('.remote-mol-stage-wrapper')
                                : document.querySelectorAll('.remote-mol-stage-wrapper');
                            wrappers.forEach(function(w) {{
                                if (w.id !== "{wrapper_id}") w.remove();
                            }});
                            if (window["{remote_resize_mol_fn}"]) {{
                                setTimeout(window["{remote_resize_mol_fn}"], 200);
                            }}
                        }}
                        setTimeout(initViewer, 0);
                    }})();
                    </script>
                    """
                )
            )
        state["traj_viewer_ready"] = True
        _set_viewer_visible(True)

    def _render_cube_in_viewer(cube_text):
        _render_3dmol(cube_text, fmt="cube", volumetric=True)

    def _frame_to_xyz(frame):
        comment, xyz_block, n_atoms = frame
        return f"{n_atoms}\n{comment}\n{xyz_block}\n"

    def _render_selected_frame():
        frames = state.get("current_xyz_frames") or []
        if not frames:
            return
        index = max(0, min(len(frames) - 1, int(state.get("current_xyz_index", 0))))
        state["current_xyz_index"] = index
        frame_input.unobserve(_on_frame_change, names="value")
        try:
            frame_input.value = index + 1
            frame_input.max = len(frames)
        finally:
            frame_input.observe(_on_frame_change, names="value")
        comment = str(frames[index][0] or "")
        large_traj_note = (
            f' <span style="color:#888;font-size:0.85em;">'
            f'(large trajectory, single-frame mode)</span>'
            if len(frames) > REMOTE_XYZ_LARGE_TRAJ_FRAMES else ""
        )
        frame_label_html.value = (
            f"{html.escape(comment[:100])}{'...' if len(comment) > 100 else ''}{large_traj_note}"
        )
        frame_total_html.value = f"<b>/ {len(frames)}</b>"
        frame_label_html.layout.display = "block"
        frame_input.layout.display = "inline-flex"
        frame_total_html.layout.display = "inline-flex"
        if len(frames) > 1 and view_toggle.value:
            _render_xyz_trajectory_viewer(initial_load=not state.get("traj_viewer_ready"))
        elif view_toggle.value:
            _render_xyz_in_viewer(_frame_to_xyz(frames[index]))

    def _join_remote_relative(base_relative, child_relative):
        base = normalize_remote_relative_path(base_relative)
        raw_child = str(child_relative or "").strip().replace("\\", "/")
        if not raw_child:
            return base
        if raw_child.startswith("/"):
            return normalize_remote_relative_path(raw_child)
        base_dir = posixpath.dirname(base) if base else ""
        combined = posixpath.join(base_dir, raw_child)
        return normalize_remote_relative_path(combined)

    def _extract_orca_xyz_block(text):
        lines = str(text or "").splitlines()
        start_idx = None
        for index, line in enumerate(lines):
            if re.match(r"^\s*\*\s*xyz\b", line, flags=re.IGNORECASE):
                start_idx = index + 1
                break
        if start_idx is None:
            return None
        coords = []
        for line in lines[start_idx:]:
            if re.match(r"^\s*\*", line):
                break
            if not line.strip():
                continue
            coords.append(line.rstrip())
        return coords if coords else None

    def _build_xyz_from_input(text, title, base_relative_path):
        coords = _extract_orca_xyz_block(text)
        if coords:
            return f"{len(coords)}\n{title}\n" + "\n".join(coords)

        match = re.search(
            r"^\s*\*\s*xyzfile\s+-?\d+\s+\d+\s+(\S+)",
            str(text or ""),
            flags=re.IGNORECASE | re.MULTILINE,
        )
        if match and state.get("config"):
            remote_ref = _join_remote_relative(base_relative_path, match.group(1))
            try:
                local_ref = fetch_remote_file(
                    state["config"]["host"],
                    state["config"]["user"],
                    state["config"]["remote_path"],
                    state["config"]["port"],
                    remote_ref,
                )
                if Path(local_ref).exists():
                    return Path(local_ref).read_text(errors="ignore")
            except Exception:
                pass

        lines = [line for line in str(text or "").strip().splitlines() if line.strip()]
        if not lines:
            return None
        if re.fullmatch(r"\d+", lines[0].strip()):
            return "\n".join(lines)
        if is_smiles(text):
            smiles_line = lines[0].strip()
            xyz_string, num_atoms, _method, error = smiles_to_xyz_quick(smiles_line)
            if not error and xyz_string:
                return f"{num_atoms}\n{title}\n{xyz_string}"
        return f"{len(lines)}\n{title}\n" + "\n".join(lines)

    def _show_file_info(entry, extra=""):
        name = str(entry.get("name") or "")
        size_str = _format_size(entry.get("size", 0))
        extra_html = f" {extra}" if extra else ""
        state["selected_remote_path"] = _current_entry_remote_path(entry)
        file_info_html.value = (
            f"<b><span style='word-break:break-all;'>{html.escape(name)}</span></b> "
            f"<span style='color:#616161;'>({html.escape(size_str)})</span>{extra_html}"
        )
        _set_selected_path_display(state["selected_remote_path"])
        copy_path_btn.disabled = False
        content_toolbar.layout.display = "flex"

    def _current_entry_remote_path(entry):
        config = state.get("config") or {}
        root_path = str(config.get("remote_path") or "/").rstrip("/") or "/"
        relative_path = normalize_remote_relative_path(entry.get("relative_path", ""))
        if not relative_path:
            return root_path
        return f"{root_path}/{relative_path}"

    def _preview_local_file(local_path, entry, *, note=""):
        path = Path(local_path)
        suffix = path.suffix.lower()
        lower_name = path.name.lower()

        state["file_preview_note"] = str(note or "")
        state["file_content"] = ""
        state["selected_file_path"] = str(path)
        try:
            state["selected_file_size"] = path.stat().st_size
        except Exception:
            state["selected_file_size"] = 0
        copy_btn.disabled = True
        download_btn.disabled = False
        _reset_search_state()
        _reset_visualization_state()

        if lower_name == "coord":
            content = path.read_text(errors="ignore")
            xyz_text = coord_to_xyz(content)
            frames = parse_xyz_frames(xyz_text) if xyz_text else []
            state["file_content"] = content
            _show_file_info(entry, "Turbomole coord")
            _render_text_preview(content, note=note)
            _set_visualization("xyz" if xyz_text else "", xyz_text or "", frames=frames)
            copy_btn.disabled = not bool(content)
            if view_toggle.value and xyz_text:
                _render_selected_frame()
            return

        if suffix == ".xyz":
            content = path.read_text(errors="ignore")
            frames = parse_xyz_frames(content)
            state["file_content"] = content
            _show_file_info(entry, f"{len(frames) or 1} frame(s)")
            _render_text_preview(content, note=note)
            _set_visualization("xyz" if frames else "", content if frames else "", frames=frames)
            copy_btn.disabled = not bool(content)
            if frames:
                _set_view_toggle(True, disabled=False)
                _update_view()
            return

        if suffix in {".png"}:
            state["file_content"] = ""
            _show_file_info(entry)
            _render_image_preview(path)
            return

        if suffix in {".cube", ".cub"}:
            content = path.read_text(errors="ignore")
            state["file_content"] = content
            _show_file_info(entry)
            preview_html.value = (
                "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
                "padding:12px; background:#fafafa;'>3D volumetric preview</div>"
            )
            _set_visualization("cube", content)
            copy_btn.disabled = not bool(content)
            _set_view_toggle(True, disabled=False)
            _update_view()
            return

        if suffix in {".doc", ".docx"}:
            state["file_content"] = ""
            _show_file_info(entry)
            _render_docx_preview(path)
            return

        if suffix in {".gbw", ".cis", ".densities", ".tmp"}:
            state["file_content"] = ""
            _show_file_info(entry)
            _render_text_preview("Binary file.\n\nPreview is not available for this file type.", note=note)
            return

        content = path.read_text(errors="ignore")
        state["file_content"] = content
        _show_file_info(entry)
        _render_text_preview(content, note=note)
        copy_btn.disabled = not bool(content)

        if suffix in {".out", ".log"}:
            coords = _extract_orca_xyz_block(content)
            xyz_text = f"{len(coords)}\n{path.name}\n" + "\n".join(coords) if coords else ""
            frames = parse_xyz_frames(xyz_text) if xyz_text else []
            _set_visualization("xyz" if xyz_text else "", xyz_text, frames=frames)
            if view_toggle.value and xyz_text:
                _render_selected_frame()
            return

        if suffix == ".inp":
            xyz_text = _build_xyz_from_input(content, path.name, entry.get("relative_path", ""))
            frames = parse_xyz_frames(xyz_text) if xyz_text else []
            _set_visualization("xyz" if xyz_text else "", xyz_text or "", frames=frames)
            if view_toggle.value and xyz_text:
                _render_selected_frame()
            return

    def _selected_entry():
        value = file_list.value
        return _entry_by_relative_path(value)

    # -- File operation handlers -----------------------------------------------

    def _on_new_folder_click(_button=None):
        new_folder_input.value = ""
        confirm_panel.children = [
            widgets.HTML("<b>New folder name:</b>"),
            new_folder_input,
            widgets.Button(description="Create", button_style="primary",
                          layout=widgets.Layout(width="70px", height="28px")),
            widgets.Button(description="Cancel",
                          layout=widgets.Layout(width="70px", height="28px")),
        ]
        confirm_panel.children[2].on_click(_on_new_folder_confirm)
        confirm_panel.children[3].on_click(_on_confirm_cancel)
        confirm_panel.layout.display = "flex"

    def _on_new_folder_confirm(_button=None):
        config = state.get("config")
        name = new_folder_input.value.strip()
        confirm_panel.layout.display = "none"
        if not config or not name:
            ops_status_html.value = '<span style="color:#d32f2f;">Invalid folder name.</span>'
            return
        try:
            remote_mkdir(config["host"], config["user"], config["remote_path"],
                        config["port"], state.get("current_relative_path", ""), name)
            ops_status_html.value = f'<span style="color:#2e7d32;">Created folder: {html.escape(name)}</span>'
            _refresh_listing(set_status=False)
        except Exception as exc:
            ops_status_html.value = f'<span style="color:#d32f2f;">Create failed: {html.escape(str(exc))}</span>'

    def _on_rename_click(_button=None):
        entry = _selected_entry()
        if not entry:
            return
        rename_input.value = entry.get("name", "")
        confirm_panel.children = [
            widgets.HTML(f"<b>Rename</b> {html.escape(entry.get('name', ''))}:"),
            rename_input,
            widgets.Button(description="Rename", button_style="primary",
                          layout=widgets.Layout(width="70px", height="28px")),
            widgets.Button(description="Cancel",
                          layout=widgets.Layout(width="70px", height="28px")),
        ]
        confirm_panel.children[2].on_click(_on_rename_confirm)
        confirm_panel.children[3].on_click(_on_confirm_cancel)
        confirm_panel.layout.display = "flex"

    def _on_rename_confirm(_button=None):
        config = state.get("config")
        entry = _selected_entry()
        new_name = rename_input.value.strip()
        confirm_panel.layout.display = "none"
        if not config or not entry or not new_name:
            ops_status_html.value = '<span style="color:#d32f2f;">Invalid input.</span>'
            return
        rel = normalize_remote_relative_path(entry.get("relative_path", ""))
        if not rel:
            return
        try:
            remote_rename(config["host"], config["user"], config["remote_path"],
                         config["port"], rel, new_name)
            ops_status_html.value = f'<span style="color:#2e7d32;">Renamed to: {html.escape(new_name)}</span>'
            _refresh_listing(set_status=False)
        except Exception as exc:
            ops_status_html.value = f'<span style="color:#d32f2f;">Rename failed: {html.escape(str(exc))}</span>'

    def _on_duplicate_click(_button=None):
        config = state.get("config")
        entry = _selected_entry()
        if not config or not entry:
            return
        rel = normalize_remote_relative_path(entry.get("relative_path", ""))
        if not rel:
            return
        try:
            new_name = remote_duplicate(config["host"], config["user"], config["remote_path"],
                                       config["port"], rel)
            ops_status_html.value = f'<span style="color:#2e7d32;">Duplicated as: {html.escape(new_name)}</span>'
            _refresh_listing(set_status=False)
        except Exception as exc:
            ops_status_html.value = f'<span style="color:#d32f2f;">Duplicate failed: {html.escape(str(exc))}</span>'

    def _on_delete_click(_button=None):
        entry = _selected_entry()
        if not entry:
            return
        name = entry.get("name", "")
        confirm_panel.children = [
            widgets.HTML(f"<b>Delete</b> <code>{html.escape(name)}</code>?"),
            widgets.Button(description="Yes, delete", button_style="danger",
                          layout=widgets.Layout(width="90px", height="28px")),
            widgets.Button(description="Cancel",
                          layout=widgets.Layout(width="70px", height="28px")),
        ]
        confirm_panel.children[1].on_click(_on_delete_confirm)
        confirm_panel.children[2].on_click(_on_confirm_cancel)
        confirm_panel.layout.display = "flex"

    def _on_delete_confirm(_button=None):
        config = state.get("config")
        entry = _selected_entry()
        confirm_panel.layout.display = "none"
        if not config or not entry:
            return
        rel = normalize_remote_relative_path(entry.get("relative_path", ""))
        if not rel:
            return
        try:
            result = remote_delete(config["host"], config["user"], config["remote_path"],
                                  config["port"], [rel])
            deleted = result.get("deleted", [])
            errors = result.get("errors", [])
            if deleted and not errors:
                ops_status_html.value = f'<span style="color:#2e7d32;">Deleted {len(deleted)} item(s).</span>'
            elif errors:
                ops_status_html.value = f'<span style="color:#d32f2f;">{html.escape("; ".join(errors))}</span>'
            _refresh_listing(set_status=False)
        except Exception as exc:
            ops_status_html.value = f'<span style="color:#d32f2f;">Delete failed: {html.escape(str(exc))}</span>'

    def _on_confirm_cancel(_button=None):
        confirm_panel.layout.display = "none"

    # --------------------------------------------------------------------------

    def _update_buttons():
        config_ready = state.get("config") is not None
        has_parent = bool(normalize_remote_relative_path(state.get("current_relative_path", "")))
        selected = _selected_entry()
        transfer_back_btn.disabled = not (config_ready and selected)
        transfer_to_archive_btn.disabled = not (config_ready and selected)
        open_btn.disabled = not bool(selected)
        up_btn.disabled = (not config_ready) or (not has_parent)
        home_btn.disabled = not config_ready
        refresh_btn.disabled = not config_ready
        file_list.disabled = not config_ready
        filter_input.disabled = not config_ready
        sort_dropdown.disabled = not config_ready
        if not selected:
            copy_path_btn.disabled = True
        copy_btn.disabled = not bool(state.get("file_content"))
        new_folder_btn.disabled = not config_ready
        rename_btn.disabled = not (config_ready and selected)
        duplicate_btn.disabled = not (config_ready and selected)
        delete_btn.disabled = not (config_ready and selected)
        if not state.get("visualize_enabled"):
            _set_view_toggle(False, disabled=True)

    def _apply_filter():
        query = str(filter_input.value or "").strip().lower()
        entries = list(state.get("entries", []))
        if query:
            entries = [entry for entry in entries if query in str(entry.get("name") or "").lower()]
        state["filtered_entries"] = entries
        if not entries:
            file_list.options = [("(No matches)", "")]
            file_list.index = 0
            state["selected_entry"] = None
            _update_buttons()
            return
        file_list.options = [(_entry_label(entry), entry.get("relative_path", "")) for entry in entries]
        previous = state.get("selected_entry") or {}
        previous_value = previous.get("relative_path", "")
        available_values = {entry.get("relative_path", "") for entry in entries}
        if previous_value in available_values:
            file_list.value = previous_value
        else:
            file_list.index = None
            state["selected_entry"] = None
        _update_buttons()

    def _load_config(set_status=True):
        try:
            config = load_transfer_settings()
        except Exception as exc:
            state["config"] = None
            info_html.value = _settings_summary(None)
            _clear_preview("Remote archive settings are invalid.")
            if set_status:
                _set_status(f"Remote archive settings invalid: {html.escape(str(exc))}", color="#d32f2f")
            _update_buttons()
            return None

        state["config"] = config
        info_html.value = _settings_summary(config)
        _update_path_html()
        if not config:
            _clear_preview("Configure a remote target in Settings first.")
            if set_status:
                _set_status("No remote target configured yet. Open Settings first.", color="#ef6c00")
        return config

    def _refresh_listing(set_status=True):
        config = state.get("config") or _load_config(set_status=False)
        if not config:
            state["entries"] = []
            state["filtered_entries"] = []
            file_list.options = [("(Configure Settings first)", "")]
            file_list.index = 0
            _update_buttons()
            return

        _update_path_html()
        _set_status("Loading remote directory ...", color="#455a64")
        ctx.set_busy(True)
        try:
            listing = list_remote_entries(
                config["host"],
                config["user"],
                config["remote_path"],
                config["port"],
                state.get("current_relative_path", ""),
                sort_mode=sort_dropdown.value,
            )
        except Exception as exc:
            state["entries"] = []
            state["filtered_entries"] = []
            file_list.options = [("(Remote directory unavailable)", "")]
            file_list.index = 0
            _clear_preview("Could not load remote directory.")
            _set_status(html.escape(str(exc)), color="#d32f2f")
            _update_buttons()
            ctx.set_busy(False)
            return
        finally:
            ctx.set_busy(False)

        state["current_relative_path"] = normalize_remote_relative_path(listing.get("relative_path", ""))
        state["entries"] = list(listing.get("entries", []))
        _apply_filter()
        if set_status:
            _set_status(
                f"Loaded {len(state['entries'])} item(s) from <code>{html.escape(_current_remote_path(config))}</code>.",
                color="#2e7d32",
            )

    def _clear_filter_for_navigation():
        if str(filter_input.value or "").strip():
            filter_input.value = ""

    def _navigate_home(_button=None):
        _clear_filter_for_navigation()
        state["current_relative_path"] = ""
        state["selected_entry"] = None
        _refresh_listing(set_status=True)

    def _navigate_up(_button=None):
        current = normalize_remote_relative_path(state.get("current_relative_path", ""))
        if not current:
            return
        _clear_filter_for_navigation()
        parent = str(PurePosixPath(current).parent)
        state["current_relative_path"] = "" if parent == "." else normalize_remote_relative_path(parent)
        state["selected_entry"] = None
        _refresh_listing(set_status=True)

    def _open_entry(entry):
        if not entry:
            return
        state["selected_entry"] = entry
        if entry.get("is_dir"):
            _clear_filter_for_navigation()
            state["current_relative_path"] = normalize_remote_relative_path(entry.get("relative_path", ""))
            _refresh_listing(set_status=True)
            return
        _preview_selected_file(entry)

    def _open_selection(_button=None):
        _open_entry(_selected_entry())

    def _preview_selected_file(entry):
        config = state.get("config")
        if not config:
            return
        state["selected_entry"] = entry
        name = str(entry.get("name") or "")
        size = int(entry.get("size") or 0)
        suffix = str(entry.get("suffix") or "").lower()
        lower_name = name.lower()
        always_fetch = lower_name == "coord" or suffix in {
            ".xyz",
            ".png",
            ".cube",
            ".cub",
            ".doc",
            ".docx",
            ".out",
            ".log",
            ".inp",
        }
        _set_status(f"Loading <code>{html.escape(name)}</code> ...", color="#455a64")
        ctx.set_busy(True)
        try:
            if always_fetch or size <= REMOTE_FULL_FETCH_MAX_BYTES:
                local_path = fetch_remote_file(
                    config["host"],
                    config["user"],
                    config["remote_path"],
                    config["port"],
                    entry.get("relative_path", ""),
                )
                note = "Remote file cached locally for preview."
                _preview_local_file(local_path, entry, note=note)
            else:
                preview = read_remote_text_preview(
                    config["host"],
                    config["user"],
                    config["remote_path"],
                    config["port"],
                    entry.get("relative_path", ""),
                    max_bytes=TEXT_PREVIEW_MAX_BYTES,
                )
                state["file_content"] = str(preview.get("text", "") or "")
                state["file_preview_note"] = "Large remote file preview."
                state["selected_file_path"] = None
                state["selected_file_size"] = size
                _reset_search_state()
                _reset_visualization_state()
                _show_file_info(entry)
                note = "Large remote file preview."
                if preview.get("truncated"):
                    note = f"{note} Only the first {TEXT_PREVIEW_MAX_BYTES:,} bytes are shown."
                _render_text_preview(preview.get("text", ""), note=note)
                copy_btn.disabled = not bool(state["file_content"])
                download_btn.disabled = True
        except Exception as exc:
            _clear_preview(f"Could not load {name}.")
            _set_status(html.escape(str(exc)), color="#d32f2f")
            return
        finally:
            ctx.set_busy(False)

        _set_status(
            f"Preview ready: <code>{html.escape(name)}</code> from <code>{html.escape(_current_remote_path(config))}</code>.",
            color="#2e7d32",
        )
        _update_buttons()
        _update_view()

    def _update_view():
        # When table panel is active, hide everything else
        if state.get("table_panel_active", False):
            content_label.layout.display = "none"
            preview_html.layout.display = "none"
            content_toolbar.layout.display = "none"
            _set_viewer_visible(False)
            return
        show_visualize = bool(view_toggle.value and state.get("visualize_enabled"))
        has_content = bool(state.get("file_content"))
        if show_visualize:
            content_label.layout.display = "none"
            preview_html.layout.display = "none"
            content_toolbar.layout.display = "none"
            _set_viewer_visible(True)
            if state.get("visualize_kind") == "cube":
                _render_cube_in_viewer(state.get("visualize_payload", ""))
            elif state.get("visualize_kind") == "xyz":
                _render_selected_frame()
            else:
                _clear_viewer()
        else:
            _stop_xyz_playback(update_button=True)
            _set_viewer_visible(False)
            content_label.layout.display = "block"
            preview_html.layout.display = "block"
            content_toolbar.layout.display = "flex" if has_content else "none"
        _update_traj_control_state()

    def _on_view_toggle(change=None):
        _update_view()

    def _on_copy_click(_button=None):
        content = str(state.get("file_content") or "")
        if not content:
            return
        _copy_to_clipboard(content, label="remote file content")

    def _on_copy_path_click(_button=None):
        remote_path = str(state.get("selected_remote_path") or "")
        if not remote_path:
            return
        _set_selected_path_display(remote_path)
        _copy_to_clipboard(remote_path, label="remote path")

    def _start_transfer_back(local_target, destination_label):
        config = state.get("config")
        entry = _selected_entry()
        if not config or not entry:
            _set_status("Select a remote file or directory first.", color="#d32f2f")
            return
        relative_path = normalize_remote_relative_path(entry.get("relative_path", ""))
        if not relative_path:
            _set_status("Select a remote file or directory first.", color="#d32f2f")
            return
        local_target = Path(local_target).resolve()
        transfer_back_btn.disabled = True
        transfer_to_archive_btn.disabled = True
        try:
            job = create_download_job(
                [relative_path],
                config["host"],
                config["user"],
                config["remote_path"],
                config["port"],
                local_target,
                delete_remote_on_success=True,
            )
            pid = launch_transfer_job(job, python_executable=sys.executable)
        except Exception as exc:
            _set_status(
                f"Transfer back failed to start: {html.escape(str(exc))}",
                color="#d32f2f",
            )
            _update_buttons()
            return
        transfer_jobs_panel.layout.display = "flex"
        _update_transfer_jobs_visibility()
        _render_transfer_jobs()
        destination = local_target / Path(relative_path)
        _set_status(
            f"Started transfer job to {html.escape(destination_label)} "
            f"<code>{html.escape(job['job_id'])}</code> for "
            f"<code>{html.escape(relative_path)}</code> into "
            f"<code>{html.escape(str(destination))}</code>. "
            "The remote source will be removed only after rsync finishes successfully. "
            f"PID: <code>{int(pid)}</code>.",
            color="#2e7d32",
        )
        _update_buttons()

    def _on_transfer_back_click(_button=None):
        _start_transfer_back(ctx.calc_dir, "Calculations")

    def _on_transfer_to_archive_click(_button=None):
        _start_transfer_back(ctx.archive_dir, "Archive")

    def _on_xyz_copy(_button=None):
        frames = state.get("current_xyz_frames") or []
        if not frames:
            return
        idx = max(0, min(len(frames) - 1, int(state.get("current_xyz_index", 0))))
        _copy_to_clipboard(_frame_to_xyz(frames[idx]).rstrip(), label=f"xyz frame {idx + 1}")

    def _start_xyz_playback():
        if not _traj_can_play():
            _set_play_button_state(False, sync_value=True)
            return
        if not state.get("traj_viewer_ready"):
            _render_selected_frame()
        _stop_xyz_playback(update_button=False)
        _set_play_button_state(True, sync_value=False)
        frame_count = len(state.get("current_xyz_frames") or [])
        try:
            fps = int(xyz_fps_input.value)
        except Exception:
            fps = REMOTE_XYZ_PLAY_FPS_DEFAULT
        fps = max(REMOTE_XYZ_PLAY_FPS_MIN, min(REMOTE_XYZ_PLAY_FPS_MAX, fps))
        delay_ms = max(16, int(round(1000.0 / float(fps))))
        loop_enabled = bool(xyz_loop_checkbox.value)
        start_frame = int(state.get("current_xyz_index", 0)) + 1
        scope_key_json = json.dumps(scope_id)
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var frameCount = {int(frame_count)};
                var loopEnabled = {str(loop_enabled).lower()};
                var delayMs = {int(delay_ms)};
                var startFrame = {int(start_frame)};
                window._remoteTrajPlayTimerByScope = window._remoteTrajPlayTimerByScope || {{}};
                if (window._remoteTrajPlayTimerByScope[scopeKey]) {{
                    clearInterval(window._remoteTrajPlayTimerByScope[scopeKey]);
                    delete window._remoteTrajPlayTimerByScope[scopeKey];
                }}
                function getViewer() {{
                    if (window._remoteTrajViewerByScope && window._remoteTrajViewerByScope[scopeKey]) {{
                        return window._remoteTrajViewerByScope[scopeKey];
                    }}
                    if (window._remoteMolViewerByScope && window._remoteMolViewerByScope[scopeKey]) {{
                        return window._remoteMolViewerByScope[scopeKey];
                    }}
                    return null;
                }}
                var scopeRoot = document.querySelector('.{scope_id}');
                if (!scopeRoot) return;
                var frameWidget = scopeRoot.querySelector('.remote-xyz-frame-input');
                var frameInput = frameWidget ? frameWidget.querySelector('input') : null;
                var playButton = scopeRoot.querySelector('.remote-xyz-play-btn button');
                if (!frameInput) return;
                var current = parseInt(frameInput.value || String(startFrame), 10);
                if (!isFinite(current) || current < 1 || current > frameCount) {{
                    current = Math.max(1, Math.min(frameCount, startFrame));
                }}
                function syncFrameValue(next) {{
                    frameInput.value = String(next);
                    frameInput.dispatchEvent(new Event('input', {{bubbles: true}}));
                    frameInput.dispatchEvent(new Event('change', {{bubbles: true}}));
                }}
                var timer = setInterval(function() {{
                    var nowVal = parseInt(frameInput.value || String(current), 10);
                    if (!isFinite(nowVal) || nowVal < 1 || nowVal > frameCount) nowVal = current;
                    var next = nowVal + 1;
                    if (next > frameCount) {{
                        if (loopEnabled) next = 1;
                        else {{
                            clearInterval(timer);
                            delete window._remoteTrajPlayTimerByScope[scopeKey];
                            if (playButton) {{
                                try {{ playButton.click(); }} catch (_e) {{}}
                            }}
                            return;
                        }}
                    }}
                    current = next;
                    syncFrameValue(next);
                    var viewer = getViewer();
                    if (viewer) {{
                        try {{
                            viewer.setFrame(next - 1);
                            viewer.render();
                        }} catch (_e) {{}}
                    }}
                }}, delayMs);
                window._remoteTrajPlayTimerByScope[scopeKey] = timer;
            }})();
            """
        )

    def _on_xyz_play_change(change):
        if change.get("name") != "value":
            return
        if state.get("traj_play_toggle_guard"):
            return
        if bool(change.get("new")):
            _start_xyz_playback()
        else:
            _stop_xyz_playback(update_button=True)

    def _on_xyz_loop_change(change):
        if change.get("name") != "value":
            return
        _update_loop_button_style()

    def _on_filter_change(change):
        _apply_filter()

    def _on_sort_change(change):
        if change.get("name") != "value":
            return
        _refresh_listing(set_status=False)

    def _on_transfer_jobs_toggle(_button=None):
        if transfer_jobs_panel.layout.display == "none":
            _render_transfer_jobs()
            transfer_jobs_panel.layout.display = "flex"
        else:
            transfer_jobs_panel.layout.display = "none"
        _update_transfer_jobs_visibility()

    def _on_transfer_jobs_refresh(_button=None):
        _render_transfer_jobs()

    def _on_selection_change(change):
        entry = _selected_entry()
        state["selected_entry"] = entry
        if not entry:
            _clear_preview()
            _update_buttons()
            return
        if entry.get("is_dir"):
            state["selected_remote_path"] = _current_entry_remote_path(entry)
            file_info_html.value = (
                f"<b><span style='word-break:break-all;'>{html.escape(str(entry.get('name') or 'Folder'))}</span></b> "
                "<span style='color:#616161;'>(directory)</span>"
            )
            _set_selected_path_display(state["selected_remote_path"])
            copy_path_btn.disabled = False
            copy_btn.disabled = True
            download_btn.disabled = True
            _reset_search_state()
            _reset_visualization_state()
            preview_html.value = (
                "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
                "padding:12px; background:#fafafa;'>"
                "Directory selected. Press <b>Enter</b> or click <b>Open</b> to enter it.</div>"
            )
        else:
            _preview_selected_file(entry)
        _update_buttons()

    def _on_keyboard_action(change):
        action = str(change.get("new") or "").strip()
        keyboard_action_input.value = ""
        if action != "open":
            return
        _open_selection()

    def _on_frame_change(change):
        if change.get("name") != "value":
            return
        frames = state.get("current_xyz_frames") or []
        if not frames:
            return
        state["current_xyz_index"] = max(0, min(len(frames) - 1, int(frame_input.value) - 1))
        _render_selected_frame()

    # -- Search helpers ---------------------------------------------------------

    def _pos_to_line_col(pos):
        if pos < 0:
            return 0, 0
        content = state.get("file_content", "")
        line = content.count("\n", 0, pos) + 1
        last_nl = content.rfind("\n", 0, pos)
        col = pos + 1 if last_nl == -1 else pos - last_nl
        return line, col

    def _update_search_nav_buttons():
        if not state["search_spans"]:
            search_prev_btn.disabled = True
            search_next_btn.disabled = True
            return
        search_prev_btn.disabled = state["current_match"] <= 0
        search_next_btn.disabled = state["current_match"] >= len(state["search_spans"]) - 1

    def _update_search_result():
        if not state["search_spans"]:
            search_result_html.value = ""
            return
        notes = []
        if state.get("search_truncated"):
            notes.append(f"limited to {REMOTE_SEARCH_MAX_MATCHES} matches")
        note_html = ""
        if notes:
            note_html = f' <span style="color:#777;">({", ".join(notes)})</span>'
        if state["current_match"] < 0:
            search_result_html.value = (
                f'<span style="color:green;">{len(state["search_spans"])} matches</span>{note_html}'
            )
            return
        start, _ = state["search_spans"][state["current_match"]]
        line, col = _pos_to_line_col(start)
        search_result_html.value = (
            f'<b>{state["current_match"] + 1}/{len(state["search_spans"])}</b> '
            f'<span style="color:#555;">(line {line}, col {col})</span>{note_html}'
        )

    def _find_matches(query):
        """Search in state['file_content'] (in-memory, case-insensitive)."""
        content = state.get("file_content", "")
        if not content or not query:
            return []
        # For very large files with a cached local path, try ripgrep
        local_path = state.get("selected_file_path")
        file_size = state.get("selected_file_size", 0) or 0
        if file_size > 5 * 1024 * 1024 and local_path and Path(local_path).exists():
            try:
                result = subprocess.run(
                    ["rg", "-a", "-i", "-F", "--byte-offset", "--only-matching", query, str(local_path)],
                    capture_output=True, text=True, timeout=30,
                )
                spans = []
                for line in result.stdout.splitlines():
                    if ":" in line:
                        parts = line.split(":", 2)
                        if len(parts) >= 2:
                            try:
                                offset = int(parts[0])
                                match_text = parts[1] if len(parts) == 2 else parts[2]
                                spans.append((offset, offset + len(match_text)))
                                if len(spans) >= REMOTE_SEARCH_MAX_MATCHES + 1:
                                    break
                            except ValueError:
                                continue
                return spans
            except (FileNotFoundError, subprocess.TimeoutExpired, Exception):
                pass
        # In-memory search
        q_lower = query.lower()
        q_len = len(query)
        spans = []
        hay = content.lower()
        pos = hay.find(q_lower)
        while pos >= 0:
            spans.append((pos, pos + q_len))
            if len(spans) >= REMOTE_SEARCH_MAX_MATCHES + 1:
                break
            pos = hay.find(q_lower, pos + 1)
        return spans

    def _should_highlight():
        return len(state.get("file_content", "")) <= REMOTE_HIGHLIGHT_MAX_CHARS

    def _apply_highlight(query, current_index):
        if query is None:
            query = ""
        js = r"""
        (function() {
            const box = document.getElementById('remote-content-box');
            const el = document.getElementById('remote-content-text');
            if (!box || !el) return;
            const text = el.textContent || '';
            const q = __QUERY__;
            if (!q) {
                el.textContent = text;
                return;
            }
            function esc(s) {
                return s.replace(/[&<>"]/g, function(c) {
                    return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];
                });
            }
            function escapeRegExp(s) {
                return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
            }
            const re = new RegExp(escapeRegExp(q), 'gi');
            let html = '';
            let last = 0;
            let i = 0;
            let m;
            while ((m = re.exec(text)) !== null) {
                html += esc(text.slice(last, m.index));
                const cls = (i === __INDEX__) ? 'remote-match current' : 'remote-match';
                const id = (i === __INDEX__) ? 'remote-current-match' : '';
                html += `<mark class="${cls}" ${id ? 'id="' + id + '"' : ''}>${esc(m[0])}</mark>`;
                last = m.index + m[0].length;
                i++;
            }
            html += esc(text.slice(last));
            el.innerHTML = html;
            if (__INDEX__ >= 0) {
                setTimeout(function() {
                    const mark = document.getElementById('remote-current-match');
                    if (mark) { mark.scrollIntoView({block: 'center'}); }
                }, 0);
            }
        })();
        """
        js = js.replace("__QUERY__", repr(query)).replace("__INDEX__", str(current_index))
        _run_js(js)

    def _scroll_to(target):
        if target == "top":
            _run_js("""
            setTimeout(function(){
                const box = document.getElementById('remote-content-box');
                if (box) { box.scrollTop = 0; }
            }, 0);
            """)
        elif target == "bottom":
            _run_js("""
            setTimeout(function(){
                const box = document.getElementById('remote-content-box');
                if (box) { box.scrollTop = box.scrollHeight; }
            }, 0);
            """)
        elif target == "match" and state["current_match"] >= 0:
            _run_js("""
            setTimeout(function(){
                const box = document.getElementById('remote-content-box');
                const el = document.getElementById('remote-current-match');
                if (!box || !el) return;
                const boxRect = box.getBoundingClientRect();
                const elRect = el.getBoundingClientRect();
                const delta = (elRect.top - boxRect.top) - (box.clientHeight / 2);
                box.scrollTop += delta;
            }, 0);
            """)

    def _reset_search_state():
        state["search_spans"] = []
        state["current_match"] = -1
        state["search_truncated"] = False
        search_result_html.value = ""
        _update_search_nav_buttons()

    def _do_search(_button=None):
        query = (search_input.value or "").strip()
        if not query and search_suggest.value and search_suggest.value != "(Select)":
            query = search_suggest.value
            search_input.value = query
        state["search_spans"] = []
        state["current_match"] = -1
        state["search_truncated"] = False

        if not query or not state.get("file_content"):
            search_result_html.value = ""
            _update_search_nav_buttons()
            if _should_highlight():
                _apply_highlight("", -1)
            return

        spans = _find_matches(query)
        if len(spans) > REMOTE_SEARCH_MAX_MATCHES:
            state["search_truncated"] = True
            spans = spans[:REMOTE_SEARCH_MAX_MATCHES]
        state["search_spans"] = spans

        if not state["search_spans"]:
            search_result_html.value = '<span style="color:red;">0 matches</span>'
            _update_search_nav_buttons()
            if _should_highlight():
                _apply_highlight("", -1)
            return

        state["current_match"] = 0
        _update_search_nav_buttons()
        _update_search_result()
        if _should_highlight():
            _apply_highlight(query, state["current_match"])
        _scroll_to("match")

    def _on_search_suggest(change):
        value = change.get("new")
        if not value or value == "(Select)":
            return
        search_input.value = value

    def _on_search_input_change(change):
        if len(state.get("file_content", "")) > REMOTE_HIGHLIGHT_MAX_CHARS:
            return
        _do_search()

    def _on_search_submit(_widget):
        _do_search()

    def _prev_match(_button=None):
        if not state["search_spans"] or state["current_match"] <= 0:
            return
        state["current_match"] -= 1
        _update_search_nav_buttons()
        _update_search_result()
        if _should_highlight():
            _apply_highlight(search_input.value.strip(), state["current_match"])
        _scroll_to("match")

    def _next_match(_button=None):
        if not state["search_spans"] or state["current_match"] >= len(state["search_spans"]) - 1:
            return
        state["current_match"] += 1
        _update_search_nav_buttons()
        _update_search_result()
        if _should_highlight():
            _apply_highlight(search_input.value.strip(), state["current_match"])
        _scroll_to("match")

    def _on_top_click(_button=None):
        _scroll_to("top")

    def _on_bottom_click(_button=None):
        _scroll_to("bottom")

    # -- Download helpers -------------------------------------------------------

    def _trigger_download(filename, payload, mime="application/octet-stream"):
        b64 = base64.b64encode(payload).decode("ascii")
        js = (
            "(function(){"
            f"const fileName={json.dumps(filename)};"
            f"const mimeType={json.dumps(mime)};"
            f"const b64={json.dumps(b64)};"
            "try{"
            "const bin=atob(b64);"
            "const len=bin.length;"
            "const bytes=new Uint8Array(len);"
            "for(let i=0;i<len;i++){bytes[i]=bin.charCodeAt(i);}"
            "const blob=new Blob([bytes],{type:mimeType});"
            "const url=URL.createObjectURL(blob);"
            "const a=document.createElement('a');"
            "a.href=url;a.download=fileName;document.body.appendChild(a);a.click();"
            "setTimeout(function(){URL.revokeObjectURL(url);if(a.parentNode){a.parentNode.removeChild(a);}},1500);"
            "}catch(err){console.error('Download failed:', err);}"
            "})();"
        )
        _run_js(js)

    def _on_download_click(_button=None):
        entry = state.get("selected_entry")
        config = state.get("config")
        if not entry or not config:
            download_status_html.value = '<span style="color:#d32f2f;">No file selected.</span>'
            return

        local_path = state.get("selected_file_path")
        if not local_path or not Path(local_path).exists():
            # Try to fetch
            try:
                local_path = fetch_remote_file(
                    config["host"], config["user"],
                    config["remote_path"], config["port"],
                    entry.get("relative_path", ""),
                )
            except Exception as exc:
                download_status_html.value = (
                    f'<span style="color:#d32f2f;">Fetch failed: {html.escape(str(exc))}</span>'
                )
                return

        try:
            target = Path(local_path)
            payload = b""
            filename = ""
            mime = "application/octet-stream"

            if target.is_dir():
                download_status_html.value = (
                    f'<span style="color:#1976d2;">Preparing ZIP for '
                    f'{html.escape(target.name or str(target))}...</span>'
                )
                with tempfile.TemporaryDirectory(prefix="delfin_download_") as tmpdir:
                    archive_name = target.name or "folder"
                    archive_base = Path(tmpdir) / archive_name
                    archive_path = Path(shutil.make_archive(
                        str(archive_base), "zip",
                        root_dir=str(target.parent),
                        base_dir=target.name,
                    ))
                    payload = archive_path.read_bytes()
                    filename = archive_path.name
                    mime = "application/zip"
            else:
                payload = target.read_bytes()
                filename = target.name

            size_bytes = len(payload)
            if size_bytes > REMOTE_DOWNLOAD_MAX_BYTES:
                download_status_html.value = (
                    f'<span style="color:#d32f2f;">Download too large '
                    f'({_format_size(size_bytes)}). Limit: {_format_size(REMOTE_DOWNLOAD_MAX_BYTES)}.</span>'
                )
                return

            _trigger_download(filename, payload, mime=mime)
            download_status_html.value = (
                f'<span style="color:#2e7d32;">Download started: '
                f'{html.escape(filename)} ({_format_size(size_bytes)})</span>'
            )
        except Exception as exc:
            download_status_html.value = (
                f'<span style="color:#d32f2f;">Download failed: {html.escape(str(exc))}</span>'
            )

    # -- Editable path input handler --------------------------------------------

    def _on_path_input_change(change):
        if _path_syncing[0]:
            return
        raw = str(change.get("new") or "").strip()
        if not raw:
            _update_path_html()
            return
        config = state.get("config")
        if not config:
            _update_path_html()
            return
        root = str(config.get("remote_path") or "/").rstrip("/") or "/"
        # Compute relative path from the input
        if raw.startswith(root):
            rel = raw[len(root):].strip("/")
        elif raw.startswith("/"):
            rel = raw.strip("/")
        else:
            rel = raw.strip("/")
        state["current_relative_path"] = normalize_remote_relative_path(rel)
        state["selected_entry"] = None
        _refresh_listing(set_status=True)

    # -- Table extraction helpers -----------------------------------------------

    def _normalize_table_number(token):
        value = str(token).strip().replace("D", "E").replace("d", "e")
        if "," in value and "." not in value:
            value = value.replace(",", ".")
        return value

    def _json_collect_key_values(node, key):
        matches = []
        stack = [node]
        while stack:
            current = stack.pop()
            if isinstance(current, dict):
                if key in current:
                    matches.append(current[key])
                for value in reversed(list(current.values())):
                    stack.append(value)
            elif isinstance(current, list):
                for item in reversed(current):
                    stack.append(item)
        return matches

    def _normalize_json_path(path):
        normalized = str(path or "").strip()
        if not normalized:
            return ""
        if normalized.startswith("$"):
            normalized = normalized[1:]
        return normalized.lstrip(".")

    def _json_split_path(path):
        parts = []
        buf = []
        bracket_depth = 0
        quote = None
        escaped = False
        for ch in path:
            if quote is not None:
                buf.append(ch)
                if escaped:
                    escaped = False
                    continue
                if ch == "\\":
                    escaped = True
                elif ch == quote:
                    quote = None
                continue
            if ch in ("'", '"') and bracket_depth > 0:
                quote = ch
                buf.append(ch)
                continue
            if ch == "[":
                bracket_depth += 1
                buf.append(ch)
                continue
            if ch == "]":
                if bracket_depth > 0:
                    bracket_depth -= 1
                buf.append(ch)
                continue
            if ch == "." and bracket_depth == 0:
                token = "".join(buf).strip()
                if token:
                    parts.append(token)
                buf = []
                continue
            buf.append(ch)
        tail = "".join(buf).strip()
        if tail:
            parts.append(tail)
        return parts

    def _json_parse_value(raw):
        token = str(raw).strip()
        if len(token) >= 2 and token[0] == token[-1] and token[0] in ('"', "'"):
            return token[1:-1]
        low = token.lower()
        if low == "true":
            return True
        if low == "false":
            return False
        if low == "null":
            return None
        try:
            if any(ch in token for ch in (".", "e", "E")):
                return float(token)
            return int(token)
        except Exception:
            return token

    def _json_compare_values(lhs, op, rhs):
        if op in ("==", "!="):
            equal = lhs == rhs
            if not equal:
                try:
                    equal = float(lhs) == float(rhs)
                except Exception:
                    equal = str(lhs) == str(rhs)
            return equal if op == "==" else (not equal)
        try:
            lval, rval = float(lhs), float(rhs)
        except Exception:
            lval, rval = str(lhs), str(rhs)
        if op == ">":
            return lval > rval
        if op == ">=":
            return lval >= rval
        if op == "<":
            return lval < rval
        if op == "<=":
            return lval <= rval
        return False

    def _json_item_matches_filter(item, expression):
        if not isinstance(item, dict):
            return False
        clauses = [c.strip() for c in re.split(r"\s*(?:&&|\band\b)\s*", expression) if c.strip()]
        if not clauses:
            return False
        for clause in clauses:
            m = re.match(r"^@\.(\w+)\s*(==|!=|>=|<=|>|<)\s*(.+)$", clause)
            if not m:
                return False
            key, op, raw_rhs = m.groups()
            rhs = _json_parse_value(raw_rhs)
            lhs = item.get(key)
            if not _json_compare_values(lhs, op, rhs):
                return False
        return True

    def _json_parse_part(part):
        text = part.strip()
        if not text:
            return "", None, None
        if "[" not in text or not text.endswith("]"):
            return text, None, None
        bracket_pos = text.find("[")
        key = text[:bracket_pos].strip()
        selector = text[bracket_pos + 1:-1].strip()
        if selector == "*":
            return key, "wildcard", None
        if re.fullmatch(r"-?\d+", selector):
            return key, "index", int(selector)
        if selector.startswith("?(") and selector.endswith(")"):
            return key, "filter", selector[2:-1].strip()
        return key, None, None

    def _json_apply_selector(values, selector_kind, selector_value):
        if selector_kind is None:
            return values
        selected = []
        for value in values:
            if selector_kind == "wildcard":
                if isinstance(value, list):
                    selected.extend(value)
                elif isinstance(value, dict):
                    selected.extend(value.values())
            elif selector_kind == "index":
                if isinstance(value, list):
                    try:
                        selected.append(value[selector_value])
                    except Exception:
                        pass
            elif selector_kind == "filter":
                if isinstance(value, list):
                    for item in value:
                        if _json_item_matches_filter(item, selector_value):
                            selected.append(item)
                elif _json_item_matches_filter(value, selector_value):
                    selected.append(value)
        return selected

    def _json_step(node, part):
        key, selector_kind, selector_value = _json_parse_part(part)
        if isinstance(node, dict):
            if key == "*":
                matches = list(node.values())
            elif key in node:
                matches = [node[key]]
            else:
                matches = _json_collect_key_values(node, key)
            return _json_apply_selector(matches, selector_kind, selector_value)
        if isinstance(node, list):
            lower_key = key.lower()
            if key in ("*", "[]"):
                matches = list(node)
            elif lower_key == "last":
                matches = [node[-1]] if node else []
            elif lower_key == "first":
                matches = [node[0]] if node else []
            elif key == "":
                matches = [node]
            else:
                try:
                    matches = [node[int(key)]]
                except Exception:
                    matches = []
                    for item in node:
                        if isinstance(item, dict) and key in item:
                            matches.append(item[key])
                    if not matches:
                        for item in node:
                            if isinstance(item, (dict, list)):
                                matches.extend(_json_step(item, key))
            return _json_apply_selector(matches, selector_kind, selector_value)
        return []

    def _json_extract_path_values(json_data, path):
        path = _normalize_json_path(path)
        parts = _json_split_path(path)
        if not parts:
            return []
        nodes = [json_data]
        for part in parts:
            next_nodes = []
            for node in nodes:
                next_nodes.extend(_json_step(node, part))
            if not next_nodes:
                return []
            nodes = next_nodes
        flat_nodes = []
        for node in nodes:
            if isinstance(node, list):
                flat_nodes.extend(node)
            else:
                flat_nodes.append(node)
        return flat_nodes

    def _extract_values(col_def, text, json_data):
        typ = col_def.get("type", "regex")
        pat = col_def.get("pattern", "").strip()
        occ = col_def.get("occ", "last")
        if not pat:
            return ["—"]
        try:
            if typ in ("regex", "text"):
                rx = re.compile(re.escape(pat) if typ == "text" else pat)
                if rx.groups > 0:
                    raw_matches = rx.finditer(text)
                    values = []
                    for m in raw_matches:
                        g = m.group(1)
                        val = str(g).strip() if g is not None else "—"
                        line_num = text[:m.start()].count("\n") + 1
                        values.append((val, line_num))
                    if not values:
                        return ["—"]
                else:
                    values = []
                    for m in rx.finditer(text):
                        line_num = text[:m.start()].count("\n") + 1
                        line_end = text.find("\n", m.end())
                        if line_end == -1:
                            line_end = len(text)
                        tail = text[m.end():line_end]
                        nums = [nm.group(0) for nm in _table_number_re.finditer(tail)]
                        if nums:
                            values.append((_normalize_table_number(nums[0]), line_num))
                            continue
                        near = text[m.end():min(len(text), m.end() + 240)]
                        near_nums = [nm.group(0) for nm in _table_number_re.finditer(near)]
                        if near_nums:
                            values.append((_normalize_table_number(near_nums[0]), line_num))
                            continue
                        values.append((m.group(0).strip(), line_num))
                    if not values:
                        return ["—"]
                if occ == "all":
                    return values
                if occ == "first":
                    return [values[0][0]]
                return [values[-1][0]]
            elif typ == "json":
                if json_data is None:
                    return ["—"]
                values = _json_extract_path_values(json_data, pat)
                if not values:
                    return ["—"]
                rendered = [str(v) for v in values]
                if occ == "all":
                    return rendered
                if occ == "first":
                    return [rendered[0]]
                return [rendered[-1]]
        except Exception as exc:
            return [f"err:{exc}"]
        return ["—"]

    def _format_table_output_value(raw_value):
        value = str(raw_value)
        if not table_decimal_comma_btn.value:
            return value
        if not _table_number_full_re.fullmatch(value.strip()):
            return value
        return value.replace(".", ",")

    def _render_extract_table_html(headers, rows):
        th = "".join(
            f'<th style="padding:4px 10px;border:1px solid #ddd;'
            f'background:#f0f4f8;white-space:nowrap;text-align:left">'
            f'{html.escape(str(h))}</th>'
            for h in headers
        )
        trs = ""
        for i, row in enumerate(rows):
            bg = "#ffffff" if i % 2 == 0 else "#f7f9fc"
            tds = "".join(
                f'<td style="padding:4px 10px;border:1px solid #ddd;'
                f'white-space:nowrap">{html.escape(str(v))}</td>'
                for v in row
            )
            trs += f'<tr style="background:{bg}">{tds}</tr>'
        return (
            '<div style="overflow-x:auto">'
            '<table style="border-collapse:collapse;font-size:12px;width:auto">'
            f'<thead><tr>{th}</tr></thead>'
            f'<tbody>{trs}</tbody>'
            '</table></div>'
        )

    def _collect_table_col_values():
        for i, rw in enumerate(state["table_col_widgets"]):
            if i < len(state["table_col_defs"]):
                state["table_col_defs"][i]["name"] = rw["name"].value
                state["table_col_defs"][i]["type"] = rw["type_dd"].value
                state["table_col_defs"][i]["occ"] = rw["occ_dd"].value
                state["table_col_defs"][i]["pattern"] = rw["pattern"].value

    def _rebuild_table_col_rows():
        rows = []
        state["table_col_widgets"] = []
        for i, col in enumerate(state["table_col_defs"]):
            name_w = widgets.Text(
                value=col.get("name", ""),
                placeholder="Column name",
                layout=widgets.Layout(width="100px", height="26px"),
            )
            type_dd = widgets.Dropdown(
                options=[("Text", "text"), ("Regex", "regex"), ("JSON path", "json")],
                value=col.get("type", "text"),
                layout=widgets.Layout(width="92px", height="26px"),
            )
            occ_dd = widgets.Dropdown(
                options=[("last", "last"), ("first", "first"), ("all", "all")],
                value=col.get("occ", "last"),
                layout=widgets.Layout(width="66px", height="26px"),
            )
            _pat_placeholders = {"regex": "regex pattern", "json": "key.path", "text": "literal text"}
            pat_w = widgets.Text(
                value=col.get("pattern", ""),
                placeholder=_pat_placeholders.get(col.get("type", "text"), "literal text"),
                layout=widgets.Layout(flex="1 1 auto", min_width="80px", height="26px"),
            )
            preset_dd = widgets.Dropdown(
                options=_tp_labels,
                value=_tp_labels[0],
                layout=widgets.Layout(width="180px", height="26px"),
            )
            rm_btn = widgets.Button(
                description="✕",
                layout=widgets.Layout(width="32px", height="26px"),
            )

            def _on_type_change(change, pw=pat_w):
                _ph = {"regex": "regex pattern", "json": "key.path", "text": "literal text"}
                pw.placeholder = _ph.get(change["new"], "literal text")

            def _on_preset(change, idx=i):
                label = change["new"]
                if label == _tp_labels[0]:
                    return
                pi = _tp_labels.index(label)
                _collect_table_col_values()
                state["table_col_defs"][idx].update({
                    "name": _tp_names[pi],
                    "type": _tp_types[pi],
                    "pattern": _tp_pats[pi],
                    "occ": _tp_occs[pi],
                })
                if _tp_types[pi] == "json" and not table_file_input.value.strip():
                    table_file_input.value = "DELFIN_Data.json"
                _rebuild_table_col_rows()

            def _on_remove(b, idx=i):
                _collect_table_col_values()
                if len(state["table_col_defs"]) > 1:
                    state["table_col_defs"].pop(idx)
                _rebuild_table_col_rows()

            type_dd.observe(_on_type_change, names="value")
            preset_dd.observe(_on_preset, names="value")
            rm_btn.on_click(_on_remove)

            rows.append(widgets.HBox(
                [name_w, type_dd, occ_dd, pat_w, preset_dd, rm_btn],
                layout=widgets.Layout(gap="4px", align_items="center", width="100%"),
            ))
            state["table_col_widgets"].append({
                "name": name_w, "type_dd": type_dd, "occ_dd": occ_dd,
                "pattern": pat_w, "preset_dd": preset_dd,
            })
        table_cols_box.children = tuple(rows)

    def _on_table_toggle(_button=None):
        if table_panel.layout.display == "none":
            table_panel.layout.display = ""
            state["table_panel_active"] = True
            _rebuild_table_col_rows()
        else:
            table_panel.layout.display = "none"
            state["table_panel_active"] = False
        _update_view()

    def _on_table_close(_button=None):
        table_panel.layout.display = "none"
        state["table_panel_active"] = False
        _update_view()

    def _on_table_add_col(_button=None):
        _collect_table_col_values()
        n = len(state["table_col_defs"]) + 1
        state["table_col_defs"].append(
            {"name": f"Value {n}", "type": "regex", "pattern": "", "occ": "last"}
        )
        _rebuild_table_col_rows()

    def _on_table_run(_button=None):
        if _list_remote_folder_files is None:
            table_status_html.value = (
                '<span style="color:#d32f2f;">Table extraction not available '
                '(list_remote_folder_files not found in remote_archive module).</span>'
            )
            return
        _collect_table_col_values()
        target_file = table_file_input.value.strip()
        if not target_file:
            table_status_html.value = '<span style="color:#d32f2f;">Please enter a filename.</span>'
            return
        config = state.get("config")
        if not config:
            table_status_html.value = '<span style="color:#d32f2f;">No remote archive configured.</span>'
            return

        current_rel = state.get("current_relative_path", "")
        table_status_html.value = '<span style="color:#1976d2;">Scanning remote folders...</span>'
        ctx.set_busy(True)
        try:
            # Use list_remote_folder_files to get all matching files in one SSH call
            file_map = _list_remote_folder_files(
                config["host"], config["user"],
                config["remote_path"], config["port"],
                current_rel, target_file,
                recursive=bool(table_recursive_cb.value),
            )
        except Exception as exc:
            table_status_html.value = (
                f'<span style="color:#d32f2f;">Remote scan failed: {html.escape(str(exc))}</span>'
            )
            ctx.set_busy(False)
            return

        if not file_map:
            table_status_html.value = '<span style="color:#d32f2f;">No matching files found.</span>'
            ctx.set_busy(False)
            return

        cols = state["table_col_defs"]
        headers = ["Job"]
        for i, c in enumerate(cols):
            headers.append(c.get("name", f"Col {i+1}"))
            if c.get("occ") == "all":
                headers.append(c.get("name", f"Col {i+1}") + " (Line)")
        _n_header_cols = len(headers) - 1
        rows = []
        missing = 0

        for folder_name, file_rel_paths in sorted(file_map.items()):
            if not file_rel_paths:
                rows.append([folder_name] + ["—"] * _n_header_cols)
                missing += 1
                continue
            for file_rel_path in file_rel_paths:
                try:
                    local_path = fetch_remote_file(
                        config["host"], config["user"],
                        config["remote_path"], config["port"],
                        file_rel_path,
                    )
                    content = Path(local_path).read_text(errors="replace")
                except Exception:
                    rows.append([folder_name] + ["err"] * _n_header_cols)
                    continue

                json_data = None
                if file_rel_path.lower().endswith(".json"):
                    try:
                        json_data = json.loads(content)
                    except Exception:
                        pass

                col_values = [_extract_values(c, content, json_data) for c in cols]
                row_count = max((len(vs) for vs in col_values), default=1)
                for row_idx in range(row_count):
                    row = [folder_name]
                    for c, values in zip(cols, col_values):
                        is_all = c.get("occ") == "all"
                        if not values or values == ["—"]:
                            row.append("—")
                            if is_all:
                                row.append("—")
                        elif row_idx < len(values):
                            v = values[row_idx]
                            if isinstance(v, tuple):
                                row.append(_format_table_output_value(v[0]))
                                if is_all:
                                    row.append(str(v[1]) if v[1] is not None else "—")
                            else:
                                row.append(_format_table_output_value(v))
                                if is_all:
                                    row.append("—")
                        else:
                            row.append("—")
                            if is_all:
                                row.append("—")
                    rows.append(row)

        ctx.set_busy(False)
        table_output_html.value = _render_extract_table_html(headers, rows)
        buf = io.StringIO()
        csv_rows = [["" if c == "—" else c for c in row] for row in rows]
        _csv.writer(buf, delimiter=";").writerows([headers] + csv_rows)
        state["table_csv_data"] = buf.getvalue()
        table_csv_btn.layout.display = "inline-flex"
        folder_count = len(file_map)
        row_count = len(rows)
        msg = f'<span style="color:#2e7d32;">{folder_count} folder(s) processed'
        if row_count != folder_count:
            msg += f", {row_count} row(s) generated"
        if missing:
            msg += f", {missing} without {html.escape(target_file)}"
        table_status_html.value = msg + ".</span>"

    def _on_table_csv_download(_button=None):
        data = state.get("table_csv_data", "")
        if not data:
            return
        b64 = base64.b64encode(data.encode("utf-8-sig")).decode()
        _run_js(
            'var a=document.createElement("a");'
            f'a.href="data:text/csv;charset=utf-8;base64,{b64}";'
            'a.download="extract_table.csv";'
            'document.body.appendChild(a);a.click();document.body.removeChild(a);'
        )

    # -- Layout -----------------------------------------------------------------

    controls_row = widgets.HBox(
        [up_btn, home_btn, refresh_btn, open_btn, new_folder_btn, rename_btn, duplicate_btn, delete_btn, table_btn],
        layout=widgets.Layout(width="100%", gap="6px", flex_flow="row wrap"),
    )
    filter_row = widgets.HBox(
        [filter_input, sort_dropdown],
        layout=widgets.Layout(width="100%", gap="6px", align_items="center"),
    )
    filter_row.add_class("remote-filter-row")
    left_panel = widgets.VBox(
        [info_html, path_input_box, controls_row, filter_row, file_list, transfer_jobs_panel, status_html],
        layout=widgets.Layout(
            flex=f"0 0 {REMOTE_LEFT_DEFAULT}px",
            min_width=f"{REMOTE_LEFT_MIN}px",
            max_width=f"{REMOTE_LEFT_MAX}px",
            padding="5px",
            gap="6px",
            overflow="hidden",
        ),
    )
    xyz_controls = widgets.HBox(
        [
            widgets.HBox(
                [widgets.HTML("<b>Frame:</b>"), frame_input, frame_total_html],
                layout=widgets.Layout(gap="10px", align_items="center", min_width="170px", flex="0 0 auto"),
            ),
            xyz_copy_btn,
        ],
        layout=widgets.Layout(
            display="none",
            gap="12px",
            margin="0 0 6px 0",
            align_items="center",
            justify_content="space-between",
            flex_flow="row nowrap",
            width="100%",
        ),
    )
    xyz_playback_row = widgets.HBox(
        [
            xyz_loop_checkbox,
            widgets.HBox(
                [widgets.HTML("<b>FPS:</b>"), xyz_fps_input],
                layout=widgets.Layout(gap="6px", align_items="center", flex="0 0 auto"),
            ),
            xyz_play_btn,
        ],
        layout=widgets.Layout(
            display="none",
            gap="12px",
            align_items="center",
            width="100%",
            justify_content="space-between",
        ),
    )
    viewer_wrap = widgets.Box(
        [viewer_output],
        layout=widgets.Layout(flex="0 0 auto", min_width="0", width="auto"),
    )
    viewer_wrap.add_class("remote-mol-view-wrap")
    xyz_tray_controls = widgets.VBox(
        [frame_label_html, xyz_controls, xyz_playback_row],
        layout=widgets.Layout(
            display="none",
            gap="14px",
            align_items="stretch",
            width="360px",
            min_width="340px",
            max_width="420px",
            margin="0",
        ),
    )
    xyz_tray_controls.add_class("remote-xyz-tray-controls")
    viewer_row = widgets.HBox(
        [viewer_wrap, xyz_tray_controls],
        layout=widgets.Layout(
            width="100%",
            gap="12px",
            align_items="flex-start",
            justify_content="flex-start",
            flex_flow="row nowrap",
        ),
    )
    viewer_row.add_class("remote-mol-view-row")
    viewer_container = widgets.VBox(
        [viewer_label, viewer_row],
        layout=widgets.Layout(display="none", margin="0 0 10px 0", width="100%", align_items="stretch"),
    )
    top_toolbar = widgets.HBox(
        [
            file_info_html,
            widgets.HBox(
                [transfer_jobs_btn, transfer_back_btn, transfer_to_archive_btn, copy_path_btn, copy_btn, download_btn, view_toggle],
                layout=widgets.Layout(
                    gap="10px",
                    flex_flow="row wrap",
                    justify_content="flex-end",
                    align_items="center",
                    width="100%",
                    overflow_x="hidden",
                ),
            ),
        ],
        layout=widgets.Layout(align_items="center", justify_content="space-between", width="100%"),
    )
    content_toolbar = widgets.HBox(
        [
            top_btn, bottom_btn,
            widgets.HTML("&nbsp;│&nbsp;"),
            search_input, search_suggest, search_btn,
            widgets.HTML("&nbsp;&nbsp;"),
            search_prev_btn, search_next_btn,
            search_result_html,
        ],
        layout=widgets.Layout(
            display="none", margin="5px 0", width="100%", overflow_x="hidden", gap="6px",
            flex_flow="row wrap", align_items="center",
        ),
    )
    download_status_row = widgets.HBox(
        [download_status_html],
        layout=widgets.Layout(width="100%"),
    )
    right_panel = widgets.VBox(
        [top_toolbar, content_toolbar, download_status_row, ops_status_html, confirm_panel, viewer_container, content_label, preview_html, table_panel],
        layout=widgets.Layout(
            flex="1 1 0",
            min_width="0",
            padding="5px",
            gap="6px",
            overflow="hidden",
        ),
    )
    splitter = widgets.HTML(
        "<div class='remote-splitter' title='Drag to resize'></div>",
        layout=widgets.Layout(height="100%", width="10px", margin="0", overflow="visible"),
    )

    css = widgets.HTML(
        "<style>"
        f".{scope_id}, .{scope_id} * {{ overflow-x:hidden !important; box-sizing:border-box; }}"
        f".{scope_id} {{ height:calc(100vh - 145px); max-height:calc(100vh - 145px); "
        "display:flex; flex-direction:column; overflow:hidden !important; }}"
        f".{scope_id} .remote-left {{ display:flex !important; flex-direction:column !important; }}"
        f".{scope_id} .remote-right {{ display:flex !important; flex-direction:column !important; overflow:hidden !important; }}"
        f".{scope_id} .widget-select select {{ height:100% !important; }}"
        f".{scope_id} .widget-select {{ flex:1 1 0 !important; min-height:0 !important; }}"
        f".{scope_id} .widget-output {{ overflow:hidden !important; }}"
        f".{scope_id} .widget-output .output_area, .{scope_id} .widget-output .output_subarea,"
        f" .{scope_id} .widget-output .output_wrapper, .{scope_id} .widget-output .jp-OutputArea-child,"
        f" .{scope_id} .widget-output .jp-OutputArea-output {{ overflow:hidden !important; margin:0 !important; padding:0 !important; width:100% !important; }}"
        f".{scope_id} .remote-left.remote-transfer-jobs-mode .remote-file-list,"
        f" .{scope_id} .remote-left.remote-transfer-jobs-mode .remote-filter-row {{ display:none !important; }}"
        f".{scope_id} .remote-splitter {{ width:8px; height:100%; cursor:col-resize;"
        " background:linear-gradient(to right, #d6d6d6, #f2f2f2, #d6d6d6);"
        " border-radius:4px; display:block; z-index:10; pointer-events:auto !important; position:relative; }"
        f".{scope_id} .remote-splitter:hover {{ background:linear-gradient("
        "to right, #b0b0b0, #e0e0e0, #b0b0b0); }"
        f".{scope_id} .remote-mol-view-row {{ width:100% !important; gap:12px !important; align-items:flex-start !important; flex-wrap:nowrap !important; }}"
        f".{scope_id} .remote-mol-view-wrap {{ flex:0 0 auto !important; min-width:0 !important; width:auto !important; }}"
        f".{scope_id} .remote-xyz-tray-controls {{ width:360px !important; min-width:340px !important; max-width:420px !important; }}"
        "</style>"
    )

    tab_widget = widgets.VBox(
        [
            css,
            title,
            keyboard_action_input,
            widgets.HBox(
                [left_panel, splitter, right_panel],
                layout=widgets.Layout(
                    width="100%",
                    align_items="stretch",
                    gap="12px",
                    flex="1 1 0",
                    min_height="0",
                ),
            ),
        ],
        layout=widgets.Layout(width="100%", max_width="100%", padding="10px", overflow="hidden"),
    )
    tab_widget.add_class(scope_id)
    tab_widget.add_class("remote-archive-tab")
    left_panel.add_class("remote-left")
    right_panel.add_class("remote-right")

    home_btn.on_click(_navigate_home)
    up_btn.on_click(_navigate_up)
    refresh_btn.on_click(lambda _button: _refresh_listing(set_status=True))
    open_btn.on_click(_open_selection)
    transfer_jobs_btn.on_click(_on_transfer_jobs_toggle)
    transfer_jobs_refresh_btn.on_click(_on_transfer_jobs_refresh)
    transfer_back_btn.on_click(_on_transfer_back_click)
    transfer_to_archive_btn.on_click(_on_transfer_to_archive_click)
    copy_btn.on_click(_on_copy_click)
    copy_path_btn.on_click(_on_copy_path_click)
    new_folder_btn.on_click(_on_new_folder_click)
    rename_btn.on_click(_on_rename_click)
    duplicate_btn.on_click(_on_duplicate_click)
    delete_btn.on_click(_on_delete_click)
    view_toggle.observe(_on_view_toggle, names="value")
    filter_input.observe(_on_filter_change, names="value")
    sort_dropdown.observe(_on_sort_change, names="value")
    file_list.observe(_on_selection_change, names="value")
    keyboard_action_input.observe(_on_keyboard_action, names="value")
    frame_input.observe(_on_frame_change, names="value")
    xyz_loop_checkbox.observe(_on_xyz_loop_change, names="value")
    xyz_play_btn.observe(_on_xyz_play_change, names="value")
    xyz_copy_btn.on_click(_on_xyz_copy)

    # Search event wiring
    search_btn.on_click(_do_search)
    search_input.observe(_on_search_input_change, names="value")
    try:
        search_input.on_submit(_on_search_submit)
    except AttributeError:
        pass
    search_suggest.observe(_on_search_suggest, names="value")
    search_prev_btn.on_click(_prev_match)
    search_next_btn.on_click(_next_match)

    # Navigation event wiring
    top_btn.on_click(_on_top_click)
    bottom_btn.on_click(_on_bottom_click)

    # Download event wiring
    download_btn.on_click(_on_download_click)

    # Path input event wiring
    path_input.observe(_on_path_input_change, names="value")

    # Table event wiring
    table_btn.on_click(_on_table_toggle)
    table_close_btn.on_click(_on_table_close)
    table_add_col_btn.on_click(_on_table_add_col)
    table_run_btn.on_click(_on_table_run)
    table_csv_btn.on_click(_on_table_csv_download)

    disable_spellcheck(ctx, class_name="remote-archive-filter")
    disable_spellcheck(ctx, class_name="remote-search-input")
    disable_spellcheck(ctx, class_name="remote-path-input")

    _clear_preview()
    _load_config(set_status=False)
    _refresh_listing(set_status=False)
    _update_path_html()
    _update_buttons()
    _update_transfer_jobs_visibility()

    init_js = f"""
    (function() {{
        function resizeRemoteArchiveViewer(scopeRoot) {{
            if (!scopeRoot || scopeRoot.offsetParent === null) return;
            var rightPanel = scopeRoot.querySelector('.remote-right');
            var mv = scopeRoot.querySelector('.remote-mol-viewer');
            if (!rightPanel || !mv || mv.offsetParent === null) return;
            var stage = mv.querySelector('[id^="remote_mol3d_"], [id^="remote_trj_viewer_"]');
            var container = mv.closest('.widget-vbox');
            if (!container || container.style.display === 'none') return;
            var mvRect = mv.getBoundingClientRect();
            if (mvRect.top === 0 && mvRect.height === 0) return;
            var rightRect = rightPanel.getBoundingClientRect();
            var topChildren = Array.prototype.slice.call(rightPanel.children || []);
            var host = null;
            for (var i = 0; i < topChildren.length; i++) {{
                if (topChildren[i].contains(mv)) {{
                    host = topChildren[i];
                    break;
                }}
            }}
            var reservedBelow = 0;
            if (host) {{
                var passed = false;
                for (var j = 0; j < topChildren.length; j++) {{
                    var child = topChildren[j];
                    if (child === host) {{
                        passed = true;
                        continue;
                    }}
                    if (!passed) continue;
                    var style = window.getComputedStyle(child);
                    if (!style || style.display === 'none' || style.visibility === 'hidden') continue;
                    var childRect = child.getBoundingClientRect();
                    if (childRect.height > 0) reservedBelow += childRect.height;
                }}
            }}
            var availH = rightRect.bottom - mvRect.top - reservedBelow - 10;
            var row = mv.closest('.remote-mol-view-row');
            var rowRect = row ? row.getBoundingClientRect() : rightRect;
            var tray = scopeRoot.querySelector('.remote-xyz-tray-controls');
            var trayStyle = tray ? window.getComputedStyle(tray) : null;
            var trayVisible = !!(tray && trayStyle && trayStyle.display !== 'none');
            var trayWidth = trayVisible ? tray.getBoundingClientRect().width : 0;
            var availW = Math.max(120, rowRect.width - trayWidth - 16);
            var h = Math.floor(availH * {REMOTE_MOL_DYNAMIC_SCALE});
            var w = Math.floor(Math.min(h * 1.2, availW));
            if (h < 80 || w < 120) return;
            mv.style.width = w + 'px';
            mv.style.height = h + 'px';
            if (stage) {{
                stage.style.width = w + 'px';
                stage.style.height = h + 'px';
            }}
            var scopeKey = {json.dumps(scope_id)};
            window._remoteMolViewStateByScope = window._remoteMolViewStateByScope || {{}};
            window._remoteMolViewScopeKeyByScope = window._remoteMolViewScopeKeyByScope || {{}};
            var viewer = null;
            if (window._remoteMolViewerByScope && window._remoteMolViewerByScope[scopeKey]) {{
                viewer = window._remoteMolViewerByScope[scopeKey];
            }} else if (window._remoteTrajViewerByScope && window._remoteTrajViewerByScope[scopeKey]) {{
                viewer = window._remoteTrajViewerByScope[scopeKey];
            }}
            var savedView = null;
            if (viewer && typeof viewer.getView === 'function') {{
                try {{
                    savedView = viewer.getView();
                    var vs = window._remoteMolViewScopeKeyByScope[scopeKey] || null;
                    if (savedView && vs) window._remoteMolViewStateByScope[vs] = savedView;
                }} catch (_e) {{}}
            }}
            if (viewer && typeof viewer.resize === 'function') {{
                try {{
                    viewer.resize();
                    viewer.render();
                    if (savedView && typeof viewer.setView === 'function') {{
                        viewer.setView(savedView);
                        viewer.render();
                    }}
                }} catch (_e) {{}}
            }}
        }}
        function _setWidgetField(root, cls, value) {{
            if (!root) return false;
            var node = root.querySelector('.' + cls);
            if (!node) return false;
            var field = null;
            if (node.matches && (node.matches('input') || node.matches('textarea'))) {{
                field = node;
            }} else if (node.querySelector) {{
                field = node.querySelector('input, textarea');
            }}
            if (!field) return false;
            var strVal = value == null ? '' : String(value);
            var nativeSet = Object.getOwnPropertyDescriptor(
                Object.getPrototypeOf(field), 'value'
            ) || Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')
              || Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value');
            if (nativeSet && nativeSet.set) {{
                nativeSet.set.call(field, strVal);
            }} else {{
                field.value = strVal;
            }}
            field.dispatchEvent(new Event('input', {{ bubbles: true }}));
            field.dispatchEvent(new Event('change', {{ bubbles: true }}));
            return true;
        }}
        function installRemoteArchiveEnter(root) {{
            if (!root) return;
            var selectEl = root.querySelector('.remote-file-list select');
            if (!selectEl) {{
                if (!root._delfinRemoteArchiveEnterRetryScheduled) {{
                    root._delfinRemoteArchiveEnterRetryScheduled = true;
                    setTimeout(function() {{
                        root._delfinRemoteArchiveEnterRetryScheduled = false;
                        installRemoteArchiveEnter(root);
                    }}, 150);
                }}
                return;
            }}
            if (selectEl.dataset.remoteEnterBound === '1') return;
            selectEl.dataset.remoteEnterBound = '1';
            selectEl.addEventListener('keydown', function(e) {{
                if (!e || e.key !== 'Enter') return;
                e.preventDefault();
                e.stopPropagation();
                _setWidgetField(root, 'remote-cmd-keyboard-action', 'open');
            }}, true);
            selectEl.addEventListener('dblclick', function(e) {{
                e.preventDefault();
                e.stopPropagation();
                _setWidgetField(root, 'remote-cmd-keyboard-action', 'open');
            }}, true);
        }}
        function installRemoteArchiveSplitter(root) {{
            if (!root) return;
            var left = root.querySelector('.remote-left');
            var right = root.querySelector('.remote-right');
            var splitter = root.querySelector('.remote-splitter');
            var scopeKey = {json.dumps(scope_id)};
            if (!left || !right || !splitter || splitter.dataset.bound === scopeKey) return;
            splitter.dataset.bound = scopeKey;
            var minW = {REMOTE_LEFT_MIN};
            var maxW = {REMOTE_LEFT_MAX};
            function onMove(e) {{
                var box = left.parentElement.getBoundingClientRect();
                var w = e.clientX - box.left;
                if (w < minW) w = minW;
                if (w > maxW) w = maxW;
                left.style.flex = '0 0 ' + w + 'px';
                left.style.minWidth = w + 'px';
                left.style.maxWidth = w + 'px';
            }}
            function onUp() {{
                document.removeEventListener('mousemove', onMove);
                document.removeEventListener('mouseup', onUp);
                if (window["{remote_resize_mol_fn}"]) {{
                    setTimeout(function() {{
                        window["{remote_resize_mol_fn}"]();
                    }}, 50);
                }}
            }}
            splitter.addEventListener('mousedown', function(e) {{
                e.preventDefault();
                document.addEventListener('mousemove', onMove);
                document.addEventListener('mouseup', onUp);
            }});
        }}
        function bootRemoteArchiveEnter() {{
            var root = document.querySelector('.{scope_id}');
            if (!root) return;
            window["{remote_resize_mol_fn}"] = function() {{
                resizeRemoteArchiveViewer(document.querySelector('.{scope_id}'));
            }};
            installRemoteArchiveEnter(root);
            installRemoteArchiveSplitter(root);
            setTimeout(window["{remote_resize_mol_fn}"], 150);
        }}
        if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', bootRemoteArchiveEnter, {{ once: true }});
        }}
        bootRemoteArchiveEnter();
        setTimeout(bootRemoteArchiveEnter, 200);
        setTimeout(bootRemoteArchiveEnter, 1000);
        if (!window._remoteArchiveResizeBoundByScope) {{
            window._remoteArchiveResizeBoundByScope = {{}};
        }}
        if (!window._remoteArchiveResizeBoundByScope[{json.dumps(scope_id)}]) {{
            window._remoteArchiveResizeBoundByScope[{json.dumps(scope_id)}] = true;
            window.addEventListener('resize', function() {{
                if (window["{remote_resize_mol_fn}"]) {{
                    window["{remote_resize_mol_fn}"]();
                }}
            }});
        }}
    }})();
    """

    return tab_widget, {"init_js": init_js}
