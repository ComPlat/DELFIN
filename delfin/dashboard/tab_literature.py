"""Literature tab: browse, upload, view, and delete reference documents.

Uses the same two-panel layout, draggable splitter, and dynamic sizing as
the Calculations and Archive tabs — but without any workflow features.
Just browse, upload, view, delete, and re-index.
"""

from __future__ import annotations

import base64
import html as _html
import json
import shutil
from pathlib import Path

import ipywidgets as widgets


# ---------------------------------------------------------------------------
# File-type helpers
# ---------------------------------------------------------------------------

_ICON_MAP = {
    '.pdf': '\U0001F4D5',   # 📕
    '.md': '\U0001F4DD',    # 📝
    '.txt': '\U0001F4C4',   # 📄
    '.rst': '\U0001F4C4',   # 📄
    '.docx': '\U0001F4D8',  # 📘
    '.doc': '\U0001F4D8',   # 📘
    '.xlsx': '\U0001F4CA',  # 📊
    '.csv': '\U0001F4CA',   # 📊
    '.png': '\U0001F5BC',   # 🖼
    '.jpg': '\U0001F5BC',   # 🖼
    '.jpeg': '\U0001F5BC',  # 🖼
    '.bib': '\U0001F4DA',   # 📚
}

_TEXT_EXTS = {
    '.md', '.txt', '.rst', '.csv', '.log', '.ini', '.cfg',
    '.yaml', '.yml', '.json', '.bib', '.tex', '.sty',
}
_IMAGE_EXTS = {'.png', '.jpg', '.jpeg', '.gif', '.svg', '.webp'}

_MAX_TEXT_BYTES = 2 * 1024 * 1024       # 2 MB
_MAX_IMAGE_BYTES = 20 * 1024 * 1024     # 20 MB
_MAX_PDF_BYTES = 50 * 1024 * 1024       # 50 MB


def _file_icon(path: Path) -> str:
    if path.is_dir():
        return '\U0001F4C2'  # 📂
    return _ICON_MAP.get(path.suffix.lower(), '\U0001F4C4')


def _label_to_name(label: str) -> str:
    """Strip emoji prefix from file-list label."""
    for i in range(1, 4):
        if i < len(label) and label[i] == ' ':
            return label[i + 1:]
    return label


# ---------------------------------------------------------------------------
# Tab factory
# ---------------------------------------------------------------------------

def create_tab(ctx):
    """Create the Literature tab.  Returns ``(tab_widget, refs_dict)``."""

    # ── resolve literature directory ──────────────────────────────────
    lit_dir = (ctx.repo_dir / 'literature') if ctx.repo_dir else (Path.cwd() / 'literature')
    lit_dir.mkdir(parents=True, exist_ok=True)

    # ── layout constants (match Calculations Browser) ─────────────────
    LEFT_DEFAULT = 375
    LEFT_MIN = 375
    LEFT_MAX = 520
    CONTENT_HEIGHT = 760

    # ── closure state ─────────────────────────────────────────────────
    state = {
        'current_path': '',
        'all_items': [],
        'selected_file': None,
    }
    scope_id = f'lit-scope-{abs(id(state))}'

    # ==================================================================
    #  LEFT PANEL — navigation, filter, file list, upload
    # ==================================================================

    # Path display
    lit_path_prefix = widgets.HTML(
        value='<b>\U0001F4DA Path:</b>',
        layout=widgets.Layout(width='100%'),
    )
    lit_path_input = widgets.Text(
        value='/',
        continuous_update=False,
        layout=widgets.Layout(
            flex='1 1 0', min_width='0', width='1px', max_width='100%',
            height='24px', overflow_x='hidden', margin='0', padding='0',
        ),
    )
    lit_path_label = widgets.VBox(
        [lit_path_prefix, lit_path_input],
        layout=widgets.Layout(width='100%', overflow_x='hidden', align_items='stretch', gap='2px'),
    )

    # Navigation buttons (identical style to Calculations Browser)
    lit_back_btn = widgets.Button(
        description='\u2B06 Up', button_style='warning',
        layout=widgets.Layout(width='58px', height='26px'), disabled=True,
    )
    lit_home_btn = widgets.Button(
        description='\U0001F3E0', button_style='info',
        layout=widgets.Layout(width='58px', height='26px'),
    )
    lit_refresh_btn = widgets.Button(
        description='\U0001F504',
        layout=widgets.Layout(width='58px', height='26px'),
    )
    lit_delete_btn = widgets.Button(
        description='\U0001F5D1 Delete', button_style='danger',
        layout=widgets.Layout(width='80px', height='26px'),
    )
    lit_reindex_btn = widgets.Button(
        description='\U0001F50D Index',
        tooltip='Rebuild the doc server search index',
        layout=widgets.Layout(width='80px', height='26px'),
    )

    lit_rename_btn = widgets.Button(
        description='Rename',
        layout=widgets.Layout(width='78px', height='26px'),
    )
    lit_rename_input = widgets.Text(
        placeholder='New name',
        layout=widgets.Layout(flex='1 1 auto', min_width='120px', height='26px', display='none'),
    )
    lit_rename_confirm = widgets.Button(
        description='OK', button_style='success',
        layout=widgets.Layout(width='50px', height='26px', display='none'),
    )
    lit_rename_cancel = widgets.Button(
        description='X',
        layout=widgets.Layout(width='36px', height='26px', display='none'),
    )

    lit_nav_controls_row = widgets.HBox(
        [lit_back_btn, lit_home_btn, lit_refresh_btn, lit_delete_btn, lit_rename_btn, lit_reindex_btn],
        layout=widgets.Layout(
            width='100%', overflow_x='hidden',
            justify_content='flex-start', gap='6px',
        ),
    )
    lit_rename_row = widgets.HBox(
        [lit_rename_input, lit_rename_confirm, lit_rename_cancel],
        layout=widgets.Layout(
            width='100%', overflow_x='hidden', display='none',
            justify_content='flex-start', gap='4px',
        ),
    )

    lit_nav_bar = widgets.VBox([
        lit_path_label,
        lit_nav_controls_row,
        lit_rename_row,
    ], layout=widgets.Layout(width='100%', overflow_x='hidden'))

    # Filter & sort row
    lit_filter = widgets.Text(
        placeholder='Filter files...', continuous_update=True,
        layout=widgets.Layout(
            flex='1 1 auto', min_width='0', width='auto', max_width='100%',
            height='28px', overflow_x='hidden', margin='0', padding='0',
        ),
    )
    lit_sort = widgets.Dropdown(
        options=[('A-Z', 'name'), ('Newest', 'date_desc'), ('Oldest', 'date_asc')],
        value='name',
        layout=widgets.Layout(width='90px', min_width='90px', height='26px', margin='0 0 0 4px'),
    )
    lit_filter_sort_row = widgets.HBox(
        [lit_filter, lit_sort],
        layout=widgets.Layout(
            width='100%', margin='0 0 8px 0',
            align_items='center', justify_content='flex-start', overflow_x='hidden',
        ),
    )

    # File list
    lit_file_list = widgets.SelectMultiple(
        options=[], rows=22,
        layout=widgets.Layout(width='100%', flex='1 1 0', min_height='0', margin='-4px 0 0 0'),
    )
    lit_file_list.add_class('lit-file-list')

    # Upload: hidden FileUpload (1px, same as Calculations) + styled trigger button
    lit_upload = widgets.FileUpload(
        accept='',
        multiple=True,
        description='',
        layout=widgets.Layout(width='1px', height='1px', overflow='hidden'),
    )
    lit_upload_btn = widgets.Button(
        description='\U0001F4E4 Upload',
        layout=widgets.Layout(width='100%', height='30px', margin='6px 0 0 0'),
    )

    # Status
    lit_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )

    lit_left = widgets.VBox([
        lit_nav_bar,
        lit_filter_sort_row,
        lit_file_list,
        lit_upload_btn,
        lit_upload,  # hidden 1px widget
    ], layout=widgets.Layout(
        flex=f'0 0 {LEFT_DEFAULT}px',
        min_width=f'{LEFT_MIN}px',
        max_width=f'{LEFT_MAX}px',
        padding='5px', overflow_x='hidden', overflow_y='hidden',
    ))

    # ==================================================================
    #  RIGHT PANEL — file info, content viewer, delete confirm
    # ==================================================================

    lit_file_info = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )

    lit_path_display = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )

    # Delete confirmation row
    lit_del_label = widgets.HTML('')
    lit_del_yes = widgets.Button(
        description='Yes', button_style='warning',
        layout=widgets.Layout(width='60px', height='26px'),
    )
    lit_del_no = widgets.Button(
        description='No',
        layout=widgets.Layout(width='60px', height='26px'),
    )
    lit_del_confirm = widgets.HBox(
        [lit_del_label, lit_del_yes, lit_del_no],
        layout=widgets.Layout(display='none', gap='6px', align_items='center'),
    )
    lit_del_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )

    # Content area (main viewer — same frame style as Calculations Browser)
    _CONTENT_FRAME = (
        "height:100%; overflow-y:auto; overflow-x:hidden;"
        " border:1px solid #ddd; padding:6px;"
        " background:#fafafa; width:100%; box-sizing:border-box;"
    )
    _lit_dir_display = str(lit_dir)
    _CONTENT_DEFAULT = (
        f"<div style='{_CONTENT_FRAME} color:#666; font-family:sans-serif;"
        f" min-height:200px; padding:20px;'>"
        f"<p style='font-size:14px;margin:0 0 12px 0;'>Select a file to preview</p>"
        f"<p style='font-size:12px;color:#888;margin:0 0 6px 0;'>"
        f"Drop files here or click <b>Upload</b> for small files.</p>"
        f"<p style='font-size:12px;color:#888;margin:0 0 6px 0;'>"
        f"For large PDFs (ORCA Manual etc.) copy directly:</p>"
        f"<code style='font-size:11px;background:#eee;padding:4px 8px;"
        f"border-radius:3px;display:inline-block;margin:4px 0;'>"
        f"cp ~/Downloads/your_file.pdf {_html.escape(_lit_dir_display)}/</code>"
        f"<p style='font-size:12px;color:#888;margin:8px 0 0 0;'>"
        f"Then click <b>\U0001F504</b> to refresh.</p>"
        f"</div>"
    )
    lit_content = widgets.HTML(
        value=_CONTENT_DEFAULT,
        layout=widgets.Layout(
            width='100%', display='block', overflow_x='hidden',
            flex='1 1 0', min_height='0',
        ),
    )
    lit_content.add_class('lit-content-area')

    lit_right = widgets.VBox([
        widgets.HBox([lit_file_info], layout=widgets.Layout(
            align_items='center', justify_content='space-between', width='100%',
        )),
        lit_path_display,
        lit_del_confirm,
        lit_del_status,
        lit_status,
        lit_content,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='5px',
        overflow_x='hidden', overflow_y='hidden',
    ))

    # ==================================================================
    #  SPLITTER (draggable, same style as Calculations Browser)
    # ==================================================================

    lit_splitter = widgets.HTML(
        "<div class='lit-splitter' title='Drag to resize'></div>",
        layout=widgets.Layout(height='100%', width='10px', margin='0', overflow='visible'),
    )

    # ==================================================================
    #  CSS (matches Calculations Browser patterns)
    # ==================================================================

    lit_css = widgets.HTML(
        '<style>'
        f'.{scope_id}, .{scope_id} * {{ overflow-x:hidden !important; box-sizing:border-box; }}'
        f'.{scope_id} {{ overflow:hidden !important;'
        ' height:calc(100vh - 145px) !important;'
        ' max-height:calc(100vh - 145px) !important;'
        ' display:flex !important; flex-direction:column !important; }}'
        '.lit-right { overflow:hidden !important;'
        ' display:flex !important; flex-direction:column !important; }'
        '.lit-left { display:flex !important; flex-direction:column !important; }'
        '.lit-left .widget-select { flex:1 1 0 !important; min-height:0 !important; }'
        '.lit-left .widget-select select { height:100% !important; }'
        '.lit-left .widget-select-multiple { flex:1 1 0 !important; min-height:0 !important; }'
        '.lit-left .widget-select-multiple select { height:100% !important; }'
        '.lit-content-area { flex:1 1 0 !important; min-height:0 !important;'
        ' overflow-y:auto !important; overflow-x:hidden !important; }'
        '.lit-content-area .widget-html-content { height:100%; }'
        '.lit-left, .lit-right { overflow-x:hidden !important; }'
        '.lit-splitter { width:8px; height:100%; cursor:col-resize;'
        ' background:linear-gradient(to right, #d6d6d6, #f2f2f2, #d6d6d6);'
        ' border-radius:4px; display:block;'
        ' z-index:10; pointer-events:auto !important; position:relative; }'
        '.lit-splitter:hover { background:linear-gradient('
        'to right, #b0b0b0, #e0e0e0, #b0b0b0); }'
        f'.{scope_id} .widget-vbox, .{scope_id} .widget-hbox {{ overflow-y:hidden !important; }}'
        '.lit-left .widget-vbox { overflow:hidden !important; }'
        '.lit-left code { display:block !important; overflow:hidden !important;'
        ' text-overflow:ellipsis !important; white-space:nowrap !important; }'
        f'.{scope_id} input::-webkit-scrollbar {{ width:0; height:0; display:none; }}'
        f'.{scope_id} input {{ scrollbar-width:none; }}'
        '.lit-file-list select { font-family:monospace !important; font-size:13px !important; }'
        # Drop-zone highlight (activated by JS drag events)
        '.lit-drop-active .lit-file-list select {'
        ' border:2px dashed #1976d2 !important; background:#e3f2fd !important; }'
        '</style>'
    )

    # ==================================================================
    #  MAIN TAB ASSEMBLY
    # ==================================================================

    tab_widget = widgets.VBox([
        lit_css,
        widgets.HTML('<h3>\U0001F4DA Literature</h3>'),
        widgets.HBox(
            [lit_left, lit_splitter, lit_right],
            layout=widgets.Layout(
                width='100%', overflow_x='hidden',
                align_items='stretch', gap='12px',
                flex='1 1 0', min_height='0',
            ),
        ),
    ], layout=widgets.Layout(
        padding='10px', overflow_x='hidden',
        width='100%', max_width='100%',
    ))
    tab_widget.add_class(scope_id)
    lit_left.add_class('lit-left')
    lit_right.add_class('lit-right')

    # ==================================================================
    #  JAVASCRIPT — splitter drag + drop-zone + dblclick navigation
    # ==================================================================

    _init_js = f"""
    (function() {{
        function initLitSplitter(attempt) {{
            var root = document.querySelector('.{scope_id}');
            if (!root) {{
                if ((attempt || 0) < 40) {{
                    setTimeout(function() {{ initLitSplitter((attempt||0)+1); }}, 100);
                    return;
                }}
                return;
            }}
            /* --- Splitter drag --- */
            var left = root.querySelector('.lit-left');
            var right = root.querySelector('.lit-right');
            var splitter = root.querySelector('.lit-splitter');
            if (left && right && splitter && splitter.dataset.bound !== '{scope_id}') {{
                splitter.dataset.bound = '{scope_id}';
                var minW = {LEFT_MIN};
                var maxW = {LEFT_MAX};
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
                }}
                splitter.addEventListener('mousedown', function(e) {{
                    e.preventDefault();
                    document.addEventListener('mousemove', onMove);
                    document.addEventListener('mouseup', onUp);
                }});
            }}

            /* --- Find the hidden upload input --- */
            function getUploadInput() {{
                return root.querySelector('.lit-left .widget-upload input[type="file"]');
            }}

            /* --- Upload button triggers the hidden file input --- */
            var uploadBtn = root.querySelector('.lit-left .jupyter-button');
            if (uploadBtn) {{
                /* Find the last button in left panel (our styled upload btn) */
                var allBtns = root.querySelectorAll('.lit-left .jupyter-button');
                var triggerBtn = allBtns[allBtns.length - 1];
                if (triggerBtn && triggerBtn.textContent.indexOf('Upload') >= 0) {{
                    triggerBtn.addEventListener('click', function(e) {{
                        e.preventDefault();
                        e.stopPropagation();
                        var inp = getUploadInput();
                        if (inp) inp.click();
                    }});
                }}
            }}

            /* --- Drop-zone visual feedback on file list --- */
            var leftEl = root.querySelector('.lit-left');
            var fileListEl = root.querySelector('.lit-file-list');
            if (fileListEl && leftEl) {{
                fileListEl.addEventListener('dragover', function(e) {{
                    e.preventDefault();
                    leftEl.classList.add('lit-drop-active');
                }});
                fileListEl.addEventListener('dragleave', function(e) {{
                    leftEl.classList.remove('lit-drop-active');
                }});
                fileListEl.addEventListener('drop', function(e) {{
                    e.preventDefault();
                    leftEl.classList.remove('lit-drop-active');
                    var inp = getUploadInput();
                    if (inp && e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files.length) {{
                        try {{
                            inp.files = e.dataTransfer.files;
                            inp.dispatchEvent(new Event('change', {{bubbles: true}}));
                        }} catch (_err) {{}}
                    }}
                }});
            }}

            /* --- Drop on the whole tab (also triggers upload) --- */
            root.addEventListener('dragover', function(e) {{
                e.preventDefault();
                var l = root.querySelector('.lit-left');
                if (l) l.classList.add('lit-drop-active');
            }});
            root.addEventListener('dragleave', function(e) {{
                if (e.relatedTarget && root.contains(e.relatedTarget)) return;
                var l = root.querySelector('.lit-left');
                if (l) l.classList.remove('lit-drop-active');
            }});
            root.addEventListener('drop', function(e) {{
                e.preventDefault();
                var l = root.querySelector('.lit-left');
                if (l) l.classList.remove('lit-drop-active');
                var inp = getUploadInput();
                if (inp && e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files.length) {{
                    try {{
                        inp.files = e.dataTransfer.files;
                        inp.dispatchEvent(new Event('change', {{bubbles: true}}));
                    }} catch (_err) {{}}
                }}
            }});
        }}
        initLitSplitter(0);
    }})();
    """

    # ==================================================================
    #  LOGIC
    # ==================================================================

    def _cur_dir() -> Path:
        return (lit_dir / state['current_path']) if state['current_path'] else lit_dir

    def _set_status(msg: str, color: str = '#666') -> None:
        lit_status.value = f'<span style="font-size:12px;color:{color}">{msg}</span>'

    def _display_path(p: Path) -> str:
        try:
            rel = p.resolve().relative_to(lit_dir.resolve())
        except Exception:
            return _html.escape(p.name)
        t = rel.as_posix()
        return '/' if t in ('', '.') else f'/{_html.escape(t)}'

    # ── List directory ────────────────────────────────────────────────

    def list_directory() -> None:
        state['all_items'] = []
        state['selected_file'] = None
        lit_delete_btn.disabled = True
        lit_del_confirm.layout.display = 'none'
        lit_del_status.value = ''
        lit_file_info.value = ''
        lit_path_display.value = ''
        lit_content.value = _CONTENT_DEFAULT
        lit_filter.value = ''

        rel = state['current_path'] or ''
        lit_path_input.value = '/' if not rel else f'/{rel}'
        lit_back_btn.disabled = (rel == '')

        cur = _cur_dir()
        if not cur.exists():
            lit_file_list.options = ['(Folder not found)']
            return

        items = []
        try:
            sort_mode = lit_sort.value
            if sort_mode == 'date_desc':
                entries = sorted(
                    cur.iterdir(),
                    key=lambda x: (not x.is_dir(), -(x.stat().st_mtime if x.exists() else 0)),
                )
            elif sort_mode == 'date_asc':
                entries = sorted(
                    cur.iterdir(),
                    key=lambda x: (not x.is_dir(), x.stat().st_mtime if x.exists() else 0),
                )
            else:
                entries = sorted(
                    cur.iterdir(),
                    key=lambda x: (not x.is_dir(), x.name.lower()),
                )
            for entry in entries:
                if entry.name.startswith('.'):
                    continue
                icon = _file_icon(entry)
                items.append(f'{icon} {entry.name}')
        except PermissionError:
            items = ['(Permission denied)']

        state['all_items'] = items
        lit_file_list.options = items if items else ['(Empty — drop files here)']

    # ── Filter ────────────────────────────────────────────────────────

    def filter_list(change=None) -> None:
        query = lit_filter.value.strip().lower()
        if not query:
            lit_file_list.options = state['all_items'] or ['(Empty — drop files here)']
            return
        filtered = [it for it in state['all_items'] if query in it.lower()]
        lit_file_list.options = filtered if filtered else ['(No matches)']

    # ── Open / preview file ───────────────────────────────────────────

    def open_file(path: Path) -> None:
        state['selected_file'] = path
        lit_delete_btn.disabled = False
        lit_del_confirm.layout.display = 'none'
        lit_del_status.value = ''

        if not path.exists():
            lit_content.value = (
                f"<div style='{_CONTENT_FRAME} font-family:monospace;'>File not found</div>"
            )
            return

        size = path.stat().st_size
        if size > 1024 * 1024:
            size_str = f'{size / (1024 * 1024):.2f} MB'
        elif size > 1024:
            size_str = f'{size / 1024:.2f} KB'
        else:
            size_str = f'{size} bytes'

        lit_file_info.value = (
            f'<b>{_html.escape(path.name)}</b>'
            f'<span style="color:#888;margin-left:8px;">({size_str})</span>'
        )
        lit_path_display.value = (
            f'<code style="font-size:11px;color:#666;">'
            f'{_html.escape(str(path))}</code>'
        )

        suffix = path.suffix.lower()

        # PDF — inline via base64 data URI iframe (with border frame)
        if suffix == '.pdf':
            if size > _MAX_PDF_BYTES:
                lit_content.value = (
                    f"<div style='{_CONTENT_FRAME} font-family:monospace;'>"
                    f"PDF too large to preview ({size_str})</div>"
                )
                return
            try:
                data = base64.b64encode(path.read_bytes()).decode('ascii')
                lit_content.value = (
                    f"<div style='{_CONTENT_FRAME} padding:0;'>"
                    f'<iframe src="data:application/pdf;base64,{data}" '
                    f'style="width:100%;height:100%;min-height:600px;border:none;'
                    f'border-radius:2px;"></iframe></div>'
                )
            except Exception as exc:
                lit_content.value = (
                    f"<div style='{_CONTENT_FRAME} font-family:monospace;'>"
                    f"Cannot display PDF: {_html.escape(str(exc))}</div>"
                )
            return

        # Images (with border frame)
        if suffix in _IMAGE_EXTS:
            if size > _MAX_IMAGE_BYTES:
                lit_content.value = (
                    f"<div style='{_CONTENT_FRAME} font-family:monospace;'>"
                    f"Image too large ({size_str})</div>"
                )
                return
            try:
                mime = {'svg': 'image/svg+xml'}.get(
                    suffix.lstrip('.'), f'image/{suffix.lstrip(".")}',
                )
                data = base64.b64encode(path.read_bytes()).decode('ascii')
                lit_content.value = (
                    f"<div style='{_CONTENT_FRAME} text-align:center; padding:10px;'>"
                    f'<img src="data:{mime};base64,{data}" '
                    f'style="max-width:100%;max-height:600px;"/></div>'
                )
            except Exception as exc:
                lit_content.value = (
                    f"<div style='{_CONTENT_FRAME} font-family:monospace;'>"
                    f"Cannot display: {_html.escape(str(exc))}</div>"
                )
            return

        # Text-like files (identical style to Calculations Browser)
        if size > _MAX_TEXT_BYTES:
            lit_content.value = (
                f"<div style='{_CONTENT_FRAME} font-family:monospace;'>"
                f"File too large to preview ({size_str})</div>"
            )
            return
        try:
            text = path.read_text(encoding='utf-8', errors='replace')
            escaped = _html.escape(text)
            lit_content.value = (
                f"<div style='{_CONTENT_FRAME}'>"
                f"<div style='white-space:pre-wrap; overflow-wrap:anywhere;"
                f" word-break:break-word; font-family:monospace; font-size:12px;"
                f" line-height:1.3;'>{escaped}</div></div>"
            )
        except Exception:
            lit_content.value = (
                f"<div style='{_CONTENT_FRAME} font-family:monospace;'>"
                f"Cannot preview this file type ({suffix})</div>"
            )

    # ── Event handlers ────────────────────────────────────────────────

    def on_select(change) -> None:
        selected = change.get('new', ())
        if not selected:
            return
        label = selected[0] if isinstance(selected, (list, tuple)) else selected
        if not label or label.startswith('('):
            return
        name = _label_to_name(label)
        full = _cur_dir() / name
        if full.is_dir():
            state['current_path'] = str(
                Path(state['current_path']) / name
            ) if state['current_path'] else name
            list_directory()
        else:
            open_file(full)

    def on_back(btn) -> None:
        if not state['current_path']:
            return
        parent = str(Path(state['current_path']).parent)
        state['current_path'] = '' if parent == '.' else parent
        list_directory()

    def on_home(btn) -> None:
        state['current_path'] = ''
        list_directory()

    def on_refresh(btn) -> None:
        list_directory()
        _set_status('Refreshed', '#2e7d32')

    def on_path_submit(change) -> None:
        raw = str(change.get('new', '/')).strip()
        if raw in ('', '/'):
            state['current_path'] = ''
        else:
            cleaned = raw.lstrip('/')
            candidate = lit_dir / cleaned
            if candidate.is_dir():
                state['current_path'] = cleaned
            else:
                _set_status(f'Not a folder: {raw}', '#d32f2f')
                return
        list_directory()

    def on_sort_change(change) -> None:
        list_directory()

    def on_reindex(btn) -> None:
        _set_status('Building index...', '#1565c0')
        try:
            from delfin.doc_server.indexer import build_index, get_default_index_path
            index = build_index(lit_dir)
            out = get_default_index_path()
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(index, indent=2, ensure_ascii=False), encoding='utf-8')
            n_docs = index['document_count']
            n_secs = sum(d['section_count'] for d in index['documents'].values())
            _set_status(f'\u2705 Index: {n_docs} docs, {n_secs} sections', '#2e7d32')
        except Exception as exc:
            _set_status(f'Index error: {exc}', '#d32f2f')

    def on_upload(change) -> None:
        entries = change.get('new', ())
        if not entries:
            return
        target = _cur_dir()
        saved = 0
        for entry in entries:
            name = entry.get('name', '') if isinstance(entry, dict) else getattr(entry, 'name', '')
            content = entry.get('content', b'') if isinstance(entry, dict) else getattr(entry, 'content', b'')
            if not name:
                continue
            dest = target / name
            if dest.exists():
                stem, sfx = dest.stem, dest.suffix
                c = 2
                while dest.exists():
                    dest = target / f'{stem}_{c}{sfx}'
                    c += 1
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(bytes(content))
            saved += 1
        if saved:
            _set_status(f'Uploaded {saved} file(s)', '#2e7d32')
            list_directory()
        try:
            lit_upload.value = ()
        except Exception:
            pass

    def on_delete_click(btn) -> None:
        sel = state.get('selected_file')
        if not sel or not sel.exists():
            return
        kind = 'folder' if sel.is_dir() else 'file'
        lit_del_label.value = (
            f'<b style="color:#d32f2f;">Delete {kind} '
            f'<code>{_html.escape(sel.name)}</code>?</b>'
        )
        lit_del_confirm.layout.display = ''

    def on_del_yes(btn) -> None:
        sel = state.get('selected_file')
        lit_del_confirm.layout.display = 'none'
        if not sel or not sel.exists():
            return
        try:
            if sel.is_dir():
                shutil.rmtree(sel)
            else:
                sel.unlink()
            lit_del_status.value = (
                f'<span style="color:#d32f2f;font-size:12px;">'
                f'Deleted {_html.escape(sel.name)}</span>'
            )
            state['selected_file'] = None
            lit_delete_btn.disabled = True
            list_directory()
        except Exception as exc:
            lit_del_status.value = (
                f'<span style="color:#d32f2f;font-size:12px;">'
                f'Delete failed: {_html.escape(str(exc))}</span>'
            )

    def on_del_no(btn) -> None:
        lit_del_confirm.layout.display = 'none'

    # ── Rename ────────────────────────────────────────────────────────

    def on_rename_click(btn) -> None:
        sel = state.get('selected_file')
        if not sel or not sel.exists():
            _set_status('Select a file first', '#d32f2f')
            return
        lit_rename_input.value = sel.name
        lit_rename_row.layout.display = ''
        lit_rename_input.layout.display = ''
        lit_rename_confirm.layout.display = ''
        lit_rename_cancel.layout.display = ''

    def on_rename_confirm(btn) -> None:
        sel = state.get('selected_file')
        new_name = lit_rename_input.value.strip()
        lit_rename_row.layout.display = 'none'
        lit_rename_input.layout.display = 'none'
        lit_rename_confirm.layout.display = 'none'
        lit_rename_cancel.layout.display = 'none'
        if not sel or not sel.exists() or not new_name:
            return
        if new_name == sel.name:
            return
        dest = sel.parent / new_name
        if dest.exists():
            _set_status(f'{new_name} already exists', '#d32f2f')
            return
        try:
            sel.rename(dest)
            _set_status(f'Renamed to {new_name}', '#2e7d32')
            state['selected_file'] = dest
            list_directory()
        except Exception as exc:
            _set_status(f'Rename failed: {exc}', '#d32f2f')

    def on_rename_cancel(btn) -> None:
        lit_rename_row.layout.display = 'none'
        lit_rename_input.layout.display = 'none'
        lit_rename_confirm.layout.display = 'none'
        lit_rename_cancel.layout.display = 'none'

    # ── Wire events ───────────────────────────────────────────────────

    lit_file_list.observe(on_select, names='value')
    lit_back_btn.on_click(on_back)
    lit_home_btn.on_click(on_home)
    lit_refresh_btn.on_click(on_refresh)
    lit_path_input.observe(on_path_submit, names='value')
    lit_sort.observe(on_sort_change, names='value')
    lit_filter.observe(filter_list, names='value')
    lit_reindex_btn.on_click(on_reindex)
    lit_upload.observe(on_upload, names='value')
    lit_delete_btn.on_click(on_delete_click)
    lit_del_yes.on_click(on_del_yes)
    lit_del_no.on_click(on_del_no)
    lit_rename_btn.on_click(on_rename_click)
    lit_rename_confirm.on_click(on_rename_confirm)
    lit_rename_cancel.on_click(on_rename_cancel)

    # ── Initial load ──────────────────────────────────────────────────
    list_directory()

    return tab_widget, {'init_js': _init_js}
