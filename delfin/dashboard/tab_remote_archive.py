"""Remote Archive tab: browse a configured remote archive over SSH."""

from __future__ import annotations

import base64
import html
import json
import posixpath
import re
import sys
from pathlib import Path, PurePosixPath

import ipywidgets as widgets
from IPython.display import HTML, clear_output, display

from delfin.remote_archive import (
    TEXT_PREVIEW_MAX_BYTES,
    fetch_remote_file,
    list_remote_entries,
    normalize_remote_relative_path,
    read_remote_text_preview,
)
from delfin.ssh_transfer_jobs import (
    create_download_job,
    ensure_jobs_dir,
    launch_transfer_job,
    list_transfer_jobs,
)
from delfin.user_settings import load_transfer_settings

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

    def _update_path_html():
        config = state.get("config")
        if not config:
            path_html.value = "<b>📂 Path:</b> /"
            return
        path_html.value = (
            f"<b>📂 Path:</b> "
            f"<code>{html.escape(_current_remote_path(config))}</code>"
        )

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
        preview_html.value = (
            "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
            "padding:12px; background:#fafafa;'>"
            f"{html.escape(message)}</div>"
        )
        copy_btn.disabled = True
        copy_path_btn.disabled = True
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
            "<div style='height:100%; overflow:auto; border:1px solid #ddd; background:#fafafa; "
            "padding:10px; box-sizing:border-box;'>"
            f"{note_html}<pre style='white-space:pre-wrap; margin:0; font-family:monospace;'>"
            f"{html.escape(content)}</pre></div>"
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
        copy_btn.disabled = True
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
            if view_toggle.value:
                _render_cube_in_viewer(content)
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

    def _navigate_home(_button=None):
        state["current_relative_path"] = ""
        state["selected_entry"] = None
        _refresh_listing(set_status=True)

    def _navigate_up(_button=None):
        current = normalize_remote_relative_path(state.get("current_relative_path", ""))
        if not current:
            return
        parent = str(PurePosixPath(current).parent)
        state["current_relative_path"] = "" if parent == "." else normalize_remote_relative_path(parent)
        state["selected_entry"] = None
        _refresh_listing(set_status=True)

    def _open_entry(entry):
        if not entry:
            return
        state["selected_entry"] = entry
        if entry.get("is_dir"):
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
                _reset_visualization_state()
                _show_file_info(entry)
                note = "Large remote file preview."
                if preview.get("truncated"):
                    note = f"{note} Only the first {TEXT_PREVIEW_MAX_BYTES:,} bytes are shown."
                _render_text_preview(preview.get("text", ""), note=note)
                copy_btn.disabled = not bool(state["file_content"])
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
        show_visualize = bool(view_toggle.value and state.get("visualize_enabled"))
        if show_visualize:
            content_label.layout.display = "none"
            preview_html.layout.display = "none"
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

    controls_row = widgets.HBox(
        [up_btn, home_btn, refresh_btn, open_btn],
        layout=widgets.Layout(width="100%", gap="6px", flex_flow="row wrap"),
    )
    filter_row = widgets.HBox(
        [filter_input, sort_dropdown],
        layout=widgets.Layout(width="100%", gap="6px", align_items="center"),
    )
    filter_row.add_class("remote-filter-row")
    left_panel = widgets.VBox(
        [info_html, path_html, controls_row, filter_row, file_list, transfer_jobs_panel, status_html],
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
                [transfer_jobs_btn, transfer_back_btn, transfer_to_archive_btn, copy_path_btn, copy_btn, view_toggle],
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
    right_panel = widgets.VBox(
        [top_toolbar, viewer_container, content_label, preview_html],
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
    view_toggle.observe(_on_view_toggle, names="value")
    filter_input.observe(_on_filter_change, names="value")
    sort_dropdown.observe(_on_sort_change, names="value")
    file_list.observe(_on_selection_change, names="value")
    keyboard_action_input.observe(_on_keyboard_action, names="value")
    frame_input.observe(_on_frame_change, names="value")
    xyz_loop_checkbox.observe(_on_xyz_loop_change, names="value")
    xyz_play_btn.observe(_on_xyz_play_change, names="value")
    xyz_copy_btn.on_click(_on_xyz_copy)

    disable_spellcheck(ctx, class_name="remote-archive-filter")

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
            var viewer = null;
            if (window._remoteMolViewerByScope && window._remoteMolViewerByScope[scopeKey]) {{
                viewer = window._remoteMolViewerByScope[scopeKey];
            }} else if (window._remoteTrajViewerByScope && window._remoteTrajViewerByScope[scopeKey]) {{
                viewer = window._remoteTrajViewerByScope[scopeKey];
            }}
            if (viewer && typeof viewer.resize === 'function') {{
                try {{
                    viewer.resize();
                    viewer.render();
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
