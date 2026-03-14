"""Remote Archive tab: read-only browsing inside the configured SSH target root."""

from __future__ import annotations

import base64
import html
import json
import posixpath
import re
from pathlib import Path, PurePosixPath

import ipywidgets as widgets
import py3Dmol
from IPython.display import HTML, clear_output, display

from delfin.remote_archive import (
    TEXT_PREVIEW_MAX_BYTES,
    fetch_remote_file,
    list_remote_entries,
    normalize_remote_relative_path,
    read_remote_text_preview,
)
from delfin.user_settings import load_transfer_settings

from .helpers import disable_spellcheck
from .input_processing import is_smiles, smiles_to_xyz_quick
from .molecule_viewer import (
    apply_molecule_view_style,
    coord_to_xyz,
    parse_xyz_frames,
)

REMOTE_FULL_FETCH_MAX_BYTES = 128 * 1024 * 1024
REMOTE_TEXT_RENDER_MAX_CHARS = 400_000
REMOTE_VIEWER_WIDTH = 700
REMOTE_VIEWER_HEIGHT = 420
REMOTE_LEFT_DEFAULT = 375
REMOTE_LEFT_MIN = 375
REMOTE_LEFT_MAX = 520


def create_tab(ctx):
    """Create the Remote Archive browser tab."""
    state = {
        "config": None,
        "current_relative_path": "",
        "entries": [],
        "filtered_entries": [],
        "selected_entry": None,
        "selected_remote_path": "",
        "current_xyz_frames": [],
        "current_xyz_index": 0,
    }
    scope_id = f"remote-archive-scope-{abs(id(state))}"

    title = widgets.HTML("<h3>📂 Remote Archive</h3>")
    info_html = widgets.HTML(value="")
    status_html = widgets.HTML(value="")
    path_html = widgets.HTML(value="<b>📂 Path:</b> /")

    open_settings_btn = widgets.Button(
        description="Settings",
        layout=widgets.Layout(width="92px", height="28px"),
    )
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
    remote_dblclick_input = widgets.Text(
        value="",
        layout=widgets.Layout(width="1px", height="1px", display="none"),
    )
    remote_dblclick_input.add_class("remote-archive-cmd-dblclick")
    file_info_html = widgets.HTML(value="")
    selected_path_html = widgets.HTML(value="")
    content_label = widgets.HTML(
        "<div style='height:26px; line-height:26px; margin:0 0 8px 0;'>"
        "<b>📄 File Preview:</b></div>"
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
    frame_total_html = widgets.HTML(value="", layout=widgets.Layout(display="none"))
    viewer_output = widgets.Output(
        layout=widgets.Layout(
            width="100%",
            border="2px solid #1976d2",
            min_height="0",
            overflow="hidden",
        )
    )
    preview_html = widgets.HTML(
        value="",
        layout=widgets.Layout(width="100%", flex="1 1 0", min_height="0", overflow="auto"),
    )

    def _set_status(message, color="#455a64"):
        status_html.value = f'<span style="color:{color};">{message}</span>'

    def _settings_summary(config):
        if not config:
            return "No remote archive configured."
        remote_root = str(config.get("remote_path") or "/")
        return (
            f'<b>Remote target:</b> '
            f'<code>{html.escape(str(config.get("user") or "?"))}@'
            f'{html.escape(str(config.get("host") or "?"))}:{html.escape(remote_root)}</code> '
            f'<span style="color:#616161;">(read-only browse inside this root)</span>'
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

    def _clear_viewer():
        viewer_output.clear_output()
        _set_viewer_visible(False)

    def _hide_frame_controls():
        state["current_xyz_frames"] = []
        state["current_xyz_index"] = 0
        frame_label_html.layout.display = "none"
        frame_input.layout.display = "none"
        frame_total_html.layout.display = "none"
        frame_label_html.value = ""
        frame_total_html.value = ""
        frame_input.max = 1
        frame_input.value = 1

    def _clear_preview(message="Select a remote file to preview."):
        file_info_html.value = ""
        selected_path_html.value = ""
        preview_html.value = (
            "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
            "padding:12px; background:#fafafa;'>"
            f"{html.escape(message)}</div>"
        )
        _clear_viewer()
        _hide_frame_controls()

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

    def _render_xyz_in_viewer(xyz_text):
        _set_viewer_visible(True)
        with viewer_output:
            clear_output()
            view = py3Dmol.view(width=REMOTE_VIEWER_WIDTH, height=REMOTE_VIEWER_HEIGHT)
            view.addModel(xyz_text, "xyz")
            apply_molecule_view_style(view)
            view.show()

    def _render_cube_in_viewer(cube_text):
        _set_viewer_visible(True)
        with viewer_output:
            clear_output()
            view = py3Dmol.view(width=REMOTE_VIEWER_WIDTH, height=REMOTE_VIEWER_HEIGHT)
            view.addModel(cube_text, "cube")
            view.addVolumetricData(cube_text, "cube", {"isoval": 0.02, "color": "#0026ff", "opacity": 0.85})
            view.addVolumetricData(cube_text, "cube", {"isoval": -0.02, "color": "#b00010", "opacity": 0.85})
            view.setBackgroundColor("white")
            view.zoomTo()
            view.center()
            view.render()
            view.show()

    def _frame_to_xyz(frame):
        comment, xyz_block, n_atoms = frame
        return f"{n_atoms}\n{comment}\n{xyz_block}\n"

    def _render_selected_frame():
        frames = state.get("current_xyz_frames") or []
        if not frames:
            return
        index = max(0, min(len(frames) - 1, int(state.get("current_xyz_index", 0))))
        state["current_xyz_index"] = index
        _render_xyz_in_viewer(_frame_to_xyz(frames[index]))
        frame_label_html.value = f"<b>Frame:</b>"
        frame_total_html.value = f"<b>/ {len(frames)}</b>"

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
        selected_path_html.value = (
            f"<span style='color:#616161; word-break:break-all;'>"
            f"<code>{html.escape(state['selected_remote_path'])}</code></span>"
        )

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

        _hide_frame_controls()
        _clear_viewer()

        if lower_name == "coord":
            content = path.read_text(errors="ignore")
            xyz_text = coord_to_xyz(content)
            _show_file_info(entry, "Turbomole coord")
            _render_text_preview(content, note=note)
            if xyz_text:
                _render_xyz_in_viewer(xyz_text)
            return

        if suffix == ".xyz":
            content = path.read_text(errors="ignore")
            frames = parse_xyz_frames(content)
            _show_file_info(entry, f"{len(frames) or 1} frame(s)")
            _render_text_preview(content, note=note)
            if frames:
                state["current_xyz_frames"] = frames
                state["current_xyz_index"] = 0
                if len(frames) > 1:
                    frame_label_html.layout.display = "inline-flex"
                    frame_input.layout.display = "inline-flex"
                    frame_total_html.layout.display = "inline-flex"
                    frame_input.max = len(frames)
                    frame_input.value = 1
                _render_selected_frame()
            return

        if suffix in {".png"}:
            _show_file_info(entry)
            _render_image_preview(path)
            return

        if suffix in {".cube", ".cub"}:
            content = path.read_text(errors="ignore")
            _show_file_info(entry)
            preview_html.value = (
                "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
                "padding:12px; background:#fafafa;'>3D volumetric preview</div>"
            )
            _render_cube_in_viewer(content)
            return

        if suffix in {".doc", ".docx"}:
            _show_file_info(entry)
            _render_docx_preview(path)
            return

        if suffix in {".gbw", ".cis", ".densities", ".tmp"}:
            _show_file_info(entry)
            _render_text_preview("Binary file.\n\nPreview is not available for this file type.", note=note)
            return

        content = path.read_text(errors="ignore")
        _show_file_info(entry)
        _render_text_preview(content, note=note)

        if suffix in {".out", ".log"}:
            coords = _extract_orca_xyz_block(content)
            if coords:
                _render_xyz_in_viewer(f"{len(coords)}\n{path.name}\n" + "\n".join(coords))
            return

        if suffix == ".inp":
            xyz_text = _build_xyz_from_input(content, path.name, entry.get("relative_path", ""))
            if xyz_text:
                _render_xyz_in_viewer(xyz_text)
            return

    def _selected_entry():
        value = file_list.value
        return _entry_by_relative_path(value)

    def _update_buttons():
        config_ready = state.get("config") is not None
        has_parent = bool(normalize_remote_relative_path(state.get("current_relative_path", "")))
        up_btn.disabled = (not config_ready) or (not has_parent)
        home_btn.disabled = not config_ready
        refresh_btn.disabled = not config_ready
        file_list.disabled = not config_ready
        filter_input.disabled = not config_ready
        sort_dropdown.disabled = not config_ready

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
                _show_file_info(entry)
                note = "Large remote file preview."
                if preview.get("truncated"):
                    note = f"{note} Only the first {TEXT_PREVIEW_MAX_BYTES:,} bytes are shown."
                _render_text_preview(preview.get("text", ""), note=note)
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

    def _on_filter_change(change):
        _apply_filter()

    def _on_sort_change(change):
        if change.get("name") != "value":
            return
        _refresh_listing(set_status=False)

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
            selected_path_html.value = (
                f"<span style='color:#616161; word-break:break-all;'>"
                f"<code>{html.escape(state['selected_remote_path'])}</code></span>"
            )
            preview_html.value = (
                "<div style='color:#616161; border:1px solid #e0e0e0; border-radius:6px; "
                "padding:12px; background:#fafafa;'>"
                "Directory selected. Double-click it or press <b>Enter</b> to enter it.</div>"
            )
            _clear_viewer()
            _hide_frame_controls()
        else:
            _preview_selected_file(entry)
        _update_buttons()

    def _on_dblclick(change):
        value = str(change.get("new") or "").strip()
        remote_dblclick_input.value = ""
        if not value:
            return
        _open_entry(_entry_by_relative_path(value))

    def _on_frame_change(change):
        if change.get("name") != "value":
            return
        frames = state.get("current_xyz_frames") or []
        if not frames:
            return
        state["current_xyz_index"] = max(0, min(len(frames) - 1, int(frame_input.value) - 1))
        _render_selected_frame()

    def _open_settings(_button=None):
        ctx.select_tab("Settings")

    controls_row = widgets.HBox(
        [up_btn, home_btn, refresh_btn, open_settings_btn],
        layout=widgets.Layout(width="100%", gap="6px", flex_flow="row wrap"),
    )
    filter_row = widgets.HBox(
        [filter_input, sort_dropdown],
        layout=widgets.Layout(width="100%", gap="6px", align_items="center"),
    )
    hidden_inputs = widgets.VBox(
        [remote_dblclick_input],
        layout=widgets.Layout(display="none"),
    )
    left_panel = widgets.VBox(
        [info_html, path_html, controls_row, hidden_inputs, filter_row, file_list, status_html],
        layout=widgets.Layout(
            flex=f"0 0 {REMOTE_LEFT_DEFAULT}px",
            min_width=f"{REMOTE_LEFT_MIN}px",
            max_width=f"{REMOTE_LEFT_MAX}px",
            padding="5px",
            gap="6px",
            overflow="hidden",
        ),
    )
    frame_row = widgets.HBox(
        [frame_label_html, frame_input, frame_total_html],
        layout=widgets.Layout(width="100%", gap="6px", align_items="center"),
    )
    viewer_container = widgets.VBox(
        [viewer_label, frame_row, viewer_output],
        layout=widgets.Layout(display="none", margin="0 0 10px 0", width="100%", align_items="stretch"),
    )
    right_panel = widgets.VBox(
        [file_info_html, selected_path_html, viewer_container, content_label, preview_html],
        layout=widgets.Layout(
            flex="1 1 0",
            min_width="0",
            padding="5px",
            gap="6px",
            overflow="hidden",
        ),
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
        "</style>"
    )

    init_js = f"""
    (function() {{
        function eventPointNode(e) {{
            if (document.elementFromPoint && e && typeof e.clientX === 'number' && typeof e.clientY === 'number') {{
                var pointed = document.elementFromPoint(e.clientX, e.clientY);
                if (pointed) return pointed;
            }}
            return e && e.target ? e.target : null;
        }}
        function optionIndex(selectEl, opt) {{
            if (!selectEl || !opt || !selectEl.options) return -1;
            return Array.prototype.indexOf.call(selectEl.options, opt);
        }}
        function optionFromIndex(selectEl, idx) {{
            if (!selectEl || !selectEl.options || idx < 0 || idx >= selectEl.options.length) return null;
            return selectEl.options[idx];
        }}
        function optionVisualHeight(selectEl) {{
            if (!selectEl || !selectEl.options || !selectEl.options.length) return 0;
            for (var i = 0; i < selectEl.options.length; i++) {{
                var rect = selectEl.options[i].getBoundingClientRect ? selectEl.options[i].getBoundingClientRect() : null;
                if (rect && rect.height > 0) return rect.height;
            }}
            return 0;
        }}
        function optionIndexAtPoint(selectEl, e) {{
            if (!selectEl || !selectEl.options || !selectEl.options.length) return -1;
            var node = eventPointNode(e);
            if (node && node.tagName === 'OPTION') return optionIndex(selectEl, node);
            if (node && node.closest) {{
                var optNode = node.closest('option');
                if (optNode) return optionIndex(selectEl, optNode);
            }}
            if (!e || typeof e.clientX !== 'number' || typeof e.clientY !== 'number') return -1;
            var rect = selectEl.getBoundingClientRect();
            if (e.clientX < rect.left || e.clientX > rect.right || e.clientY < rect.top || e.clientY > rect.bottom) return -1;
            var optionCount = selectEl.options.length;
            var optionHeight = optionVisualHeight(selectEl);
            if ((!optionHeight || !isFinite(optionHeight)) && selectEl.scrollHeight && selectEl.scrollHeight > 0) {{
                optionHeight = selectEl.scrollHeight / optionCount;
            }}
            if ((!optionHeight || !isFinite(optionHeight)) && rect.height > 0) {{
                var visibleRows = Math.max(1, Math.min(optionCount, Number(selectEl.size) || optionCount));
                optionHeight = rect.height / visibleRows;
            }}
            if (!optionHeight || !isFinite(optionHeight)) return -1;
            var localY = (e.clientY - rect.top) + (selectEl.scrollTop || 0);
            var contentHeight = optionHeight * optionCount;
            if (localY < 0 || localY >= contentHeight) return -1;
            var rawIdx = Math.floor(localY / optionHeight);
            if (rawIdx < 0) rawIdx = 0;
            if (rawIdx >= optionCount) rawIdx = optionCount - 1;
            return rawIdx;
        }}
        function optionAtPoint(selectEl, e) {{
            return optionFromIndex(selectEl, optionIndexAtPoint(selectEl, e));
        }}
        function setWidgetInput(root, cls, value) {{
            if (!root) return false;
            var field = root.querySelector('.' + cls + ' input, .' + cls + ' textarea');
            if (!field) return false;
            var strVal = String(value == null ? '' : value);
            var nativeSet = Object.getOwnPropertyDescriptor(Object.getPrototypeOf(field), 'value')
                || Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')
                || Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value');
            if (nativeSet && nativeSet.set) nativeSet.set.call(field, strVal);
            else field.value = strVal;
            field.dispatchEvent(new Event('input', {{ bubbles: true }}));
            field.dispatchEvent(new Event('change', {{ bubbles: true }}));
            return true;
        }}
        function initRemoteScope(attempt) {{
            var root = document.querySelector('.{scope_id}');
            if (!root) {{
                if ((attempt || 0) < 40) setTimeout(function() {{ initRemoteScope((attempt || 0) + 1); }}, 250);
                return;
            }}
            if (root._remoteArchiveExplorerReady) return;
            root._remoteArchiveExplorerReady = true;
            root._dblLastTime = 0;
            root._dblLastValue = '';
            root._dblLastX = 0;
            root._dblLastY = 0;

            root.addEventListener('mousedown', function(e) {{
                if (e.button != null && e.button !== 0) return;
                var selectEl = root.querySelector('.remote-file-list select');
                if (!selectEl) return;
                var currentOpt = optionAtPoint(selectEl, e);
                if (!currentOpt) return;
                var currentLabel = String(currentOpt.textContent || currentOpt.innerText || '').trim();
                var currentValue = String(currentOpt.value || '');
                if (!currentLabel || currentLabel.charAt(0) === '(' || !currentValue) return;
                var now = Date.now();
                var dx = Math.abs((typeof e.clientX === 'number' ? e.clientX : 0) - (Number(root._dblLastX) || 0));
                var dy = Math.abs((typeof e.clientY === 'number' ? e.clientY : 0) - (Number(root._dblLastY) || 0));
                if (currentValue === String(root._dblLastValue || '') && (now - root._dblLastTime) < 500 && dx <= 20 && dy <= 20) {{
                    root._dblLastTime = 0;
                    root._dblLastValue = '';
                    root._dblLastX = 0;
                    root._dblLastY = 0;
                    e.preventDefault();
                    e.stopPropagation();
                    setWidgetInput(root, 'remote-archive-cmd-dblclick', currentValue);
                    return;
                }}
                root._dblLastTime = now;
                root._dblLastValue = currentValue;
                root._dblLastX = (typeof e.clientX === 'number') ? e.clientX : 0;
                root._dblLastY = (typeof e.clientY === 'number') ? e.clientY : 0;
            }}, true);

            root.addEventListener('keydown', function(e) {{
                var selectEl = root.querySelector('.remote-file-list select');
                if (!selectEl || e.target !== selectEl) return;
                if (e.key === 'Enter') {{
                    var selectedIndex = typeof selectEl.selectedIndex === 'number' ? selectEl.selectedIndex : -1;
                    var selectedOpt = optionFromIndex(selectEl, selectedIndex);
                    if (!selectedOpt) return;
                    var currentLabel = String(selectedOpt.textContent || selectedOpt.innerText || '').trim();
                    var currentValue = String(selectedOpt.value || '');
                    if (!currentLabel || currentLabel.charAt(0) === '(' || !currentValue) return;
                    e.preventDefault();
                    e.stopPropagation();
                    setWidgetInput(root, 'remote-archive-cmd-dblclick', currentValue);
                }}
            }}, true);
        }}
        initRemoteScope(0);
    }})();
    """

    tab_widget = widgets.VBox(
        [
            css,
            title,
            widgets.HBox(
                [left_panel, right_panel],
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

    open_settings_btn.on_click(_open_settings)
    home_btn.on_click(_navigate_home)
    up_btn.on_click(_navigate_up)
    refresh_btn.on_click(lambda _button: _refresh_listing(set_status=True))
    filter_input.observe(_on_filter_change, names="value")
    sort_dropdown.observe(_on_sort_change, names="value")
    file_list.observe(_on_selection_change, names="value")
    remote_dblclick_input.observe(_on_dblclick, names="value")
    frame_input.observe(_on_frame_change, names="value")

    disable_spellcheck(ctx, class_name="remote-archive-filter")

    _clear_preview()
    _load_config(set_status=False)
    _refresh_listing(set_status=False)
    _update_path_html()
    _update_buttons()

    return tab_widget, {"init_js": init_js}
