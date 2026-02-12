"""Calculations Browser tab: file browser with viewer, search, recalc, and more."""

import base64
import html as _html
import json
import re
import shutil
import subprocess
import threading
from pathlib import Path

import ipywidgets as widgets
import py3Dmol
from IPython.display import HTML, clear_output, display

from .constants import CALC_SEARCH_OPTIONS
from .input_processing import parse_inp_resources, parse_resource_settings
from .helpers import disable_spellcheck
from .molecule_viewer import coord_to_xyz, parse_xyz_frames
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import MolToImage


def create_tab(ctx):
    """Create the Calculations Browser tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- layout constants ---------------------------------------------------
    CALC_CONTENT_HEIGHT = 760
    CALC_MOL_SIZE = 480
    CALC_MOL_ZOOM = 0.9
    CALC_LEFT_DEFAULT = 320
    CALC_LEFT_MIN = 320
    CALC_LEFT_MAX = 520
    # -- state (closure-captured) -------------------------------------------
    state = {
        'current_path': '',
        'file_content': '',
        'all_items': [],
        'search_spans': [],
        'current_match': -1,
        'selected_inp_path': None,
        'selected_inp_base': '',
        'recalc_active': False,
        'delete_current': False,
        'xyz_frames': [],
        'xyz_current_frame': [0],
        'report_running': {},
        'traj_viewer_ready': False,
        'preselect': {
            'active': False,
            'entries': [],
            'index': 0,
            'decisions': {},
            'csv_path': '',
            'job_dir': '',
        },
    }

    # -- widgets ------------------------------------------------------------
    calc_path_label = widgets.HTML(
        value='<b>📂 Path:</b> /',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_back_btn = widgets.Button(
        description='⬆ Up', button_style='warning',
        layout=widgets.Layout(width='58px', height='26px'), disabled=True,
    )
    calc_home_btn = widgets.Button(
        description='🏠', button_style='info',
        layout=widgets.Layout(width='58px', height='26px'),
    )
    calc_refresh_btn = widgets.Button(
        description='🔄',
        layout=widgets.Layout(width='58px', height='26px'),
    )
    calc_delete_btn = widgets.Button(
        description='🗑 Delete', button_style='danger',
        layout=widgets.Layout(width='80px', height='26px'),
    )

    # Delete confirmation
    calc_delete_yes_btn = widgets.Button(
        description='Yes', button_style='danger',
        layout=widgets.Layout(width='60px', height='26px'),
    )
    calc_delete_no_btn = widgets.Button(
        description='No',
        layout=widgets.Layout(width='60px', height='26px'),
    )
    calc_delete_label = widgets.HTML('<b>Delete selected item?</b>')
    calc_delete_confirm = widgets.HBox(
        [calc_delete_label, calc_delete_yes_btn, calc_delete_no_btn],
        layout=widgets.Layout(display='none', gap='6px', align_items='center'),
    )
    calc_delete_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )

    # Filter & sort
    calc_folder_search = widgets.Text(
        placeholder='Filter files...', continuous_update=True,
        layout=widgets.Layout(
            flex='1 1 auto', min_width='0', width='auto', max_width='100%',
            height='28px', overflow_x='hidden', margin='0', padding='0',
        ),
    )
    calc_sort_dropdown = widgets.Dropdown(
        options=[('A-Z', 'name'), ('Newest', 'date_desc'), ('Oldest', 'date_asc')],
        value='name',
        layout=widgets.Layout(width='90px', min_width='90px', height='26px', margin='0 0 0 4px'),
    )
    calc_filter_sort_row = widgets.HBox(
        [calc_folder_search, calc_sort_dropdown],
        layout=widgets.Layout(
            width='100%', margin='0 0 8px 0',
            align_items='center', justify_content='flex-start',
            overflow_x='hidden',
        ),
    )
    calc_filter_sort_row.add_class('calc-filter-row')

    # File list
    calc_file_list = widgets.Select(
        options=[], rows=22,
        layout=widgets.Layout(
            width='100%', flex='1 1 0', min_height='0', margin='-4px 0 0 0'
        ),
    )

    # Content area
    calc_content_area = widgets.HTML(
        value='',
        layout=widgets.Layout(
            width='100%', display='block', overflow_x='hidden',
            flex='1 1 0', min_height='0',
        ),
    )

    # Molecule viewer
    calc_mol_label = widgets.HTML(
        "<div style='height:26px; line-height:26px; margin:0 0 8px 0;'>"
        "<b>🔬 Molecule Preview:</b></div>"
    )
    calc_mol_viewer = widgets.Output(
        layout=widgets.Layout(
            width=f'{CALC_MOL_SIZE}px', height=f'{CALC_MOL_SIZE}px',
            border='2px solid #1976d2',
            overflow='hidden', padding='0', border_radius='0',
        ),
    )
    calc_mol_viewer.add_class('calc-mol-viewer')

    # XYZ trajectory controls
    calc_xyz_frame_label = widgets.HTML(
        value='', layout=widgets.Layout(margin='4px 0'),
    )
    calc_xyz_frame_input = widgets.BoundedIntText(
        value=1, min=1, max=1, step=1,
        layout=widgets.Layout(width='60px', height='28px'),
    )
    calc_xyz_frame_total = widgets.HTML(
        value='<b>/ 1</b>', layout=widgets.Layout(width='40px'),
    )
    calc_xyz_copy_btn = widgets.Button(
        description='📋 Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='160px', height='32px'),
    )
    calc_xyz_controls = widgets.HBox(
        [widgets.HTML('<b>Frame:</b>'), calc_xyz_frame_input, calc_xyz_frame_total, calc_xyz_copy_btn],
        layout=widgets.Layout(display='none', gap='6px', margin='6px 0', align_items='center'),
    )

    # Coord (Turbomole) copy button
    calc_coord_copy_btn = widgets.Button(
        description='📋 Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='160px', height='32px'),
    )
    calc_coord_controls = widgets.HBox(
        [calc_coord_copy_btn],
        layout=widgets.Layout(display='none', gap='6px', margin='6px 0', align_items='center'),
    )

    # Copy / path / report buttons
    calc_copy_path_btn = widgets.Button(
        description='PATH', button_style='',
        layout=widgets.Layout(width='70px', min_width='70px', height='26px'), disabled=True,
    )
    calc_copy_btn = widgets.Button(
        description='Copy', button_style='info',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px'), disabled=True,
    )
    calc_report_btn = widgets.Button(
        description='Report', button_style='success',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px'), disabled=True,
    )
    calc_report_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%'),
    )

    # File info / path display / content label
    calc_file_info = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_path_display = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_content_label = widgets.HTML(
        value='<b>📄 File Content:</b>',
        layout=widgets.Layout(width='100%', overflow_x='hidden', margin='8px 0 0 0'),
    )
    calc_view_toggle = widgets.ToggleButton(
        description='Visualize', value=False, disabled=True, button_style='warning',
        layout=widgets.Layout(width='110px', min_width='110px', height='26px'),
    )

    # Recalc widgets
    calc_recalc_btn = widgets.Button(
        description='Recalc', button_style='warning',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px'), disabled=True,
    )
    calc_submit_recalc_btn = widgets.Button(
        description='Submit Recalc', button_style='success',
        layout=widgets.Layout(width='130px', min_width='130px', height='26px'), disabled=True,
    )
    calc_recalc_time = widgets.Text(
        value='24:00:00', description='Time limit',
        layout=widgets.Layout(width='200px', height='26px'),
    )
    calc_recalc_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_edit_area = widgets.Textarea(
        value='',
        layout=widgets.Layout(
            width='100%', flex='1 1 0', min_height='0', display='none'
        ),
    )
    calc_edit_area.add_class('delfin-nospell')
    calc_edit_area.add_class('calc-edit-area')
    disable_spellcheck(ctx)

    # Content navigation
    calc_top_btn = widgets.Button(
        description='⬆ Top',
        layout=widgets.Layout(width='95px', min_width='95px', height='26px'),
    )
    calc_bottom_btn = widgets.Button(
        description='⬇ End',
        layout=widgets.Layout(width='95px', min_width='95px', height='26px'),
    )

    # Content search
    calc_search_input = widgets.Text(
        placeholder='Search in file...', continuous_update=False,
        layout=widgets.Layout(width='140px', min_width='140px', height='26px'),
    )
    calc_search_suggest = widgets.Dropdown(
        options=['(Select)'] + CALC_SEARCH_OPTIONS,
        value='(Select)',
        layout=widgets.Layout(width='200px', min_width='200px', height='26px'),
    )
    calc_search_btn = widgets.Button(
        description='🔍',
        layout=widgets.Layout(width='85px', min_width='85px', height='26px'),
    )
    calc_prev_btn = widgets.Button(
        description='◀',
        layout=widgets.Layout(width='58px', min_width='58px', height='26px'), disabled=True,
    )
    calc_next_btn = widgets.Button(
        description='▶',
        layout=widgets.Layout(width='58px', min_width='58px', height='26px'), disabled=True,
    )
    calc_search_result = widgets.HTML(
        value='', layout=widgets.Layout(width='180px', min_width='180px'),
    )

    # OCCUPIER Override widgets
    calc_options_dropdown = widgets.Dropdown(
        options=['(Options)'], value='(Options)',
        layout=widgets.Layout(width='150px', min_width='150px', height='26px', display='none'),
    )
    calc_options_dropdown.add_class('calc-options-dropdown')
    calc_override_input = widgets.Text(
        value='', placeholder='STAGE=INDEX',
        layout=widgets.Layout(width='140px', min_width='140px', height='26px', display='none'),
    )
    calc_override_time = widgets.Text(
        value='08:00:00', placeholder='HH:MM:SS',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px', display='none'),
    )
    calc_override_btn = widgets.Button(
        description='Submit', button_style='success',
        layout=widgets.Layout(width='70px', min_width='70px', height='26px', display='none'),
    )
    calc_override_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', display='none'),
    )

    # -- preselection widgets ----------------------------------------------
    calc_preselect_title = widgets.HTML('<b>Preselection</b>')
    calc_preselect_progress = widgets.HTML('')
    calc_preselect_info = widgets.HTML('')
    calc_preselect_img = widgets.Output()
    calc_preselect_yes = widgets.Button(
        description='Yes', button_style='success',
        layout=widgets.Layout(width='80px', height='30px'),
    )
    calc_preselect_no = widgets.Button(
        description='No', button_style='danger',
        layout=widgets.Layout(width='80px', height='30px'),
    )
    calc_preselect_status = widgets.HTML('')
    calc_preselect_container = widgets.VBox([
        calc_preselect_title,
        calc_preselect_progress,
        calc_preselect_info,
        calc_preselect_img,
        widgets.HBox([calc_preselect_yes, calc_preselect_no], layout=widgets.Layout(gap='10px')),
        calc_preselect_status,
    ], layout=widgets.Layout(display='none', margin='10px 0', width='100%'))

    # -- compound layout widgets (must be defined before closures that use them)
    calc_mol_container = widgets.VBox([
        calc_mol_label,
        calc_xyz_frame_label,
        calc_xyz_controls,
        calc_coord_controls,
        calc_mol_viewer,
    ], layout=widgets.Layout(display='none', margin='0 0 10px 0', width='100%', align_items='flex-end'))

    calc_content_toolbar = widgets.HBox([
        calc_top_btn, calc_bottom_btn,
        widgets.HTML('&nbsp;│&nbsp;'),
        calc_search_input, calc_search_suggest, calc_search_btn,
        widgets.HTML('&nbsp;&nbsp;'),
        calc_prev_btn, calc_next_btn,
        calc_options_dropdown, calc_override_input, calc_override_time, calc_override_btn,
        calc_search_result,
    ], layout=widgets.Layout(
        margin='5px 0', width='100%', overflow_x='hidden', gap='6px',
        flex_flow='row wrap', align_items='center',
    ))

    calc_recalc_toolbar = widgets.HBox([
        calc_recalc_time,
        calc_submit_recalc_btn,
        calc_recalc_status,
    ], layout=widgets.Layout(
        margin='5px 0', width='100%', overflow_x='hidden', gap='8px',
        align_items='center', display='none',
    ))

    # -- helper closures ----------------------------------------------------
    RDLogger.DisableLog('rdApp.*')

    def _calc_is_complete_mutation_csv(selected_name):
        return selected_name and selected_name.lower() == 'complete_mutation_space.csv'

    def _calc_get_selected_path():
        selected = calc_file_list.value
        if not selected or selected.startswith('('):
            return None
        name = selected[2:].strip()
        return (
            (_calc_dir() / state['current_path'] / name)
            if state['current_path']
            else (_calc_dir() / name)
        )

    def _calc_parse_mutation_csv(csv_path):
        entries = []
        try:
            lines = csv_path.read_text(encoding='utf-8', errors='ignore').splitlines()
        except Exception:
            return entries
        for line in lines:
            if not line.strip():
                continue
            if line.lower().startswith('label;smiles'):
                continue
            if ';' not in line:
                continue
            label, smiles = line.split(';', 1)
            label = label.strip()
            smiles = smiles.strip()
            if not smiles:
                continue
            entries.append((label, smiles))
        return entries

    def _calc_read_decision_file(path):
        entries = set()
        if not path.exists():
            return entries
        try:
            lines = path.read_text(encoding='utf-8', errors='ignore').splitlines()
        except Exception:
            return entries
        for line in lines:
            if not line.strip():
                continue
            if line.lower().startswith('label;smiles'):
                continue
            if ';' not in line:
                continue
            label, smiles = line.split(';', 1)
            label = label.strip()
            smiles = smiles.strip()
            if label or smiles:
                entries.add((label, smiles))
        return entries

    def _calc_next_preselect_index(start_idx):
        entries = state['preselect']['entries']
        decisions = state['preselect']['decisions']
        idx = max(0, start_idx)
        while idx < len(entries):
            if str(idx) not in decisions:
                return idx
            idx += 1
        return len(entries)

    def _calc_preselect_write_outputs():
        entries = state['preselect']['entries']
        decisions = state['preselect']['decisions']
        job_dir = Path(state['preselect']['job_dir'])
        if not job_dir:
            return
        preselected_path = job_dir / 'preselected_mutation_space.csv'
        rejected_path = job_dir / 'rejected_mutation_space.csv'
        yes_lines = ['Label;SMILES']
        no_lines = ['Label;SMILES']
        for i, (label, smiles) in enumerate(entries):
            decision = decisions.get(str(i))
            if decision == 'yes':
                yes_lines.append(f'{label};{smiles}')
            elif decision == 'no':
                no_lines.append(f'{label};{smiles}')
        preselected_path.write_text('\n'.join(yes_lines), encoding='utf-8')
        rejected_path.write_text('\n'.join(no_lines), encoding='utf-8')

    def _calc_preselect_save_state():
        job_dir = Path(state['preselect']['job_dir'])
        if not job_dir:
            return
        state_path = job_dir / 'selection_state.json'
        payload = {
            'csv_path': state['preselect']['csv_path'],
            'index': state['preselect']['index'],
            'decisions': state['preselect']['decisions'],
        }
        state_path.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    def _calc_preselect_load(csv_path):
        entries = _calc_parse_mutation_csv(csv_path)
        job_dir = csv_path.parent
        state_path = job_dir / 'selection_state.json'
        decisions = {}
        index = 0

        if state_path.exists():
            try:
                data = json.loads(state_path.read_text(encoding='utf-8', errors='ignore'))
                if data.get('csv_path') == str(csv_path):
                    decisions = data.get('decisions', {}) or {}
                    index = int(data.get('index', 0))
            except Exception:
                decisions = {}

        if not decisions:
            yes_set = _calc_read_decision_file(job_dir / 'preselected_mutation_space.csv')
            no_set = _calc_read_decision_file(job_dir / 'rejected_mutation_space.csv')
            for i, entry in enumerate(entries):
                if entry in yes_set:
                    decisions[str(i)] = 'yes'
                elif entry in no_set:
                    decisions[str(i)] = 'no'

        state['preselect']['entries'] = entries
        state['preselect']['decisions'] = decisions
        state['preselect']['csv_path'] = str(csv_path)
        state['preselect']['job_dir'] = str(job_dir)
        state['preselect']['index'] = _calc_next_preselect_index(index)
        calc_preselect_title.value = f'<b>Preselection:</b> {csv_path.name}'

    def _calc_preselect_render_current():
        entries = state['preselect']['entries']
        idx = state['preselect']['index']
        total = len(entries)
        if total == 0:
            calc_preselect_progress.value = '<b>No entries found.</b>'
            calc_preselect_info.value = ''
            calc_preselect_status.value = '<span style="color:#d32f2f;">Empty or invalid CSV.</span>'
            calc_preselect_yes.disabled = True
            calc_preselect_no.disabled = True
            with calc_preselect_img:
                clear_output(wait=True)
            return
        if idx >= total:
            calc_preselect_progress.value = f'<b>Done.</b> ({total} / {total})'
            calc_preselect_info.value = ''
            calc_preselect_status.value = '<span style="color:green;">All entries processed.</span>'
            calc_preselect_yes.disabled = True
            calc_preselect_no.disabled = True
            with calc_preselect_img:
                clear_output(wait=True)
            return

        label, smiles = entries[idx]
        calc_preselect_progress.value = f'<b>{idx + 1} / {total}</b>'
        calc_preselect_info.value = (
            f'<div><b>Label:</b> {_html.escape(label)}</div>'
            f'<div><b>SMILES:</b> <code>{_html.escape(smiles)}</code></div>'
        )
        calc_preselect_status.value = ''
        calc_preselect_yes.disabled = False
        calc_preselect_no.disabled = False
        with calc_preselect_img:
            clear_output(wait=True)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                display(HTML('<i>Invalid SMILES</i>'))
            else:
                try:
                    rdDepictor.Compute2DCoords(mol)
                except Exception:
                    pass
                img = MolToImage(mol, size=(360, 360))
                display(img)

    def _calc_preselect_show(show):
        if show:
            calc_preselect_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_content_label.layout.display = 'none'
            calc_mol_container.layout.display = 'none'
            calc_edit_area.layout.display = 'none'
        else:
            calc_preselect_container.layout.display = 'none'

    def _calc_preselect_decide(decision):
        entries = state['preselect']['entries']
        idx = state['preselect']['index']
        if idx >= len(entries):
            return
        state['preselect']['decisions'][str(idx)] = decision
        state['preselect']['index'] = _calc_next_preselect_index(idx + 1)
        _calc_preselect_write_outputs()
        _calc_preselect_save_state()
        _calc_preselect_render_current()

    def _run_js(script):
        ctx.run_js(script)

    def _render_3dmol(data, fmt='xyz', stick_radius=0.1, sphere_scale=0.25,
                      extra_fn=None):
        """Render a 3D molecule via py3Dmol with dynamic resize afterwards."""
        calc_mol_viewer.clear_output()
        with calc_mol_viewer:
            view = py3Dmol.view(width=CALC_MOL_SIZE, height=CALC_MOL_SIZE)
            view.addModel(data, fmt)
            view.setStyle({}, {'stick': {'radius': stick_radius}, 'sphere': {'scale': sphere_scale}})
            if extra_fn:
                extra_fn(view, data)
            view.zoomTo()
            view.zoom(CALC_MOL_ZOOM)
            view.show()

    def _set_view_toggle(value, disabled=None):
        """Set visualize toggle without triggering multiple UI refreshes."""
        try:
            calc_view_toggle.unobserve(_on_view_toggle, names='value')
        except Exception:
            pass
        if disabled is not None:
            calc_view_toggle.disabled = disabled
        calc_view_toggle.value = value
        calc_view_toggle.observe(_on_view_toggle, names='value')

    def _on_view_toggle(change=None):
        calc_update_view()

    def _calc_dir():
        return ctx.calc_dir

    # -- view / display helpers ---------------------------------------------
    def calc_update_view():
        show_mol = calc_view_toggle.value
        if show_mol:
            calc_mol_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_edit_area.layout.display = 'none'
            calc_content_label.layout.display = 'none'
            calc_content_toolbar.layout.display = 'none'
            calc_recalc_toolbar.layout.display = 'none'
        else:
            calc_mol_container.layout.display = 'none'
            calc_content_label.layout.display = 'block'
            if state['recalc_active']:
                calc_content_toolbar.layout.display = 'none'
                calc_content_area.layout.display = 'none'
                calc_edit_area.layout.display = 'block'
                calc_recalc_toolbar.layout.display = 'flex'
            else:
                calc_content_toolbar.layout.display = 'flex'
                calc_content_area.layout.display = 'block'
                calc_edit_area.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'none'

    def calc_set_message(message):
        calc_content_area.value = (
            "<div style='height:100%; overflow-y:auto;"
            " border:1px solid #ddd; padding:6px; font-family:monospace;"
            " white-space:pre; background:#fafafa;'>"
            f"{_html.escape(message)}"
            "</div>"
        )

    def calc_update_path_display():
        display_path = '/' + state['current_path'] if state['current_path'] else '/'
        calc_path_label.value = (
            f'<b>📂 Path:</b> <code style="background:#eee;padding:2px 5px;'
            f' white-space:nowrap; overflow:hidden; text-overflow:ellipsis;'
            f' display:block; max-width:100%;" title="{display_path}">{display_path}</code>'
        )
        calc_back_btn.disabled = (state['current_path'] == '')

    def calc_update_report_btn():
        current_dir = _calc_dir() / state['current_path'] if state['current_path'] else _calc_dir()
        calc_report_btn.disabled = not (current_dir / 'CONTROL.txt').exists()

    def calc_render_content(scroll_to=None):
        if not state['file_content']:
            calc_set_message('Select a file...')
            return
        text = state['file_content']
        calc_content_area.value = (
            "<style>"
            ".calc-match { background: #fff59d; padding: 0 2px; }"
            ".calc-match.current { background: #ffcc80; }"
            "</style>"
            "<div id='calc-content-box' style='height:100%;"
            " overflow-y:auto; overflow-x:hidden; border:1px solid #ddd; padding:6px;"
            " background:#fafafa; width:100%; box-sizing:border-box;'>"
            "<div id='calc-content-text' style='white-space:pre-wrap; overflow-wrap:anywhere;"
            " word-break:break-word; font-family:monospace; font-size:12px; line-height:1.3;'>"
            f"{_html.escape(text)}"
            "</div></div>"
        )
        if scroll_to:
            calc_scroll_to(scroll_to)

    def calc_scroll_to(target):
        if target == 'top':
            _run_js("""
            setTimeout(function(){
                const box = document.getElementById('calc-content-box');
                if (box) { box.scrollTop = 0; }
            }, 0);
            """)
        elif target == 'bottom':
            _run_js("""
            setTimeout(function(){
                const box = document.getElementById('calc-content-box');
                if (box) { box.scrollTop = box.scrollHeight; }
            }, 0);
            """)
        elif target == 'match' and state['current_match'] >= 0:
            _run_js("""
            setTimeout(function(){
                const box = document.getElementById('calc-content-box');
                const el = document.getElementById('calc-current-match');
                if (!box || !el) return;
                const boxRect = box.getBoundingClientRect();
                const elRect = el.getBoundingClientRect();
                const delta = (elRect.top - boxRect.top) - (box.clientHeight / 2);
                box.scrollTop += delta;
            }, 0);
            """)

    def calc_apply_highlight(query, current_index):
        if query is None:
            query = ''
        js = r"""
        (function() {
            const box = document.getElementById('calc-content-box');
            const el = document.getElementById('calc-content-text');
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
                return s.replace(/[.*+?^${}()|[\\]\\]/g, '\\$&');
            }
            const re = new RegExp(escapeRegExp(q), 'gi');
            let html = '';
            let last = 0;
            let i = 0;
            let m;
            while ((m = re.exec(text)) !== null) {
                html += esc(text.slice(last, m.index));
                const cls = (i === __INDEX__) ? 'calc-match current' : 'calc-match';
                const id = (i === __INDEX__) ? 'calc-current-match' : '';
                html += `<mark class="${cls}" ${id ? 'id="' + id + '"' : ''}>${esc(m[0])}</mark>`;
                last = m.index + m[0].length;
                i++;
            }
            html += esc(text.slice(last));
            el.innerHTML = html;
        })();
        """
        js = js.replace('__QUERY__', repr(query)).replace('__INDEX__', str(current_index))
        _run_js(js)

    # -- XYZ helpers --------------------------------------------------------
    def calc_extract_orca_xyz_block(text):
        lines = text.splitlines()
        start_idx = None
        for i, line in enumerate(lines):
            if re.match(r'^\s*\*\s*xyz\b', line, flags=re.IGNORECASE):
                start_idx = i + 1
                break
        if start_idx is None:
            return None
        coords = []
        for line in lines[start_idx:]:
            if re.match(r'^\s*\*', line):
                break
            if not line.strip():
                continue
            coords.append(line.rstrip())
        return coords if coords else None

    def calc_build_xyz_from_input(text, title='input.txt', base_dir=None):
        coords = calc_extract_orca_xyz_block(text)
        if coords:
            atom_count = len(coords)
            return f"{atom_count}\n{title}\n" + '\n'.join(coords)
        # ORCA: XYZFile reference
        m = re.search(r'^\s*\*\s*xyzfile\s+-?\d+\s+\d+\s+(\S+)',
                      text, flags=re.IGNORECASE | re.MULTILINE)
        if m:
            xyz_path = m.group(1)
            try:
                p = Path(xyz_path)
                if not p.is_absolute() and base_dir:
                    p = Path(base_dir) / p
                if p.exists():
                    return p.read_text(errors='ignore')
            except Exception:
                pass
        lines = [ln for ln in text.strip().splitlines() if ln.strip()]
        if not lines:
            return None
        if re.fullmatch(r'\d+', lines[0].strip()):
            return '\n'.join(lines)
        atom_count = len(lines)
        return f"{atom_count}\n{title}\n" + '\n'.join(lines)

    def calc_update_xyz_viewer(initial_load=False):
        frames = state['xyz_frames']
        if not frames:
            return
        # Ensure viewer area is visible before rendering to avoid flicker.
        calc_mol_container.layout.display = 'block'
        calc_content_area.layout.display = 'none'
        idx = state['xyz_current_frame'][0]
        if idx < 0 or idx >= len(frames):
            return
        comment, xyz_block, n_atoms = frames[idx]

        # Update frame input (suppress observer temporarily)
        calc_xyz_frame_input.unobserve(calc_on_xyz_input_change, names='value')
        calc_xyz_frame_input.value = idx + 1
        calc_xyz_frame_input.max = len(frames)
        calc_xyz_frame_input.observe(calc_on_xyz_input_change, names='value')
        calc_xyz_frame_total.value = f"<b>/ {len(frames)}</b>"
        calc_xyz_frame_label.value = (
            f"{_html.escape(comment[:100])}{'...' if len(comment) > 100 else ''}"
        )

        # Trajectory: use JS viewer and keep orientation when switching frames
        if len(frames) > 1:
            if initial_load or not state['traj_viewer_ready']:
                full_xyz = ""
                for comm, block, natoms in frames:
                    full_xyz += f"{natoms}\n{comm}\n{block}\n"
                calc_mol_viewer.clear_output()
                with calc_mol_viewer:
                    viewer_id = "calc_trj_viewer"
                    html_content = f"""
                    <div id="{viewer_id}" style="width:{CALC_MOL_SIZE}px;height:{CALC_MOL_SIZE}px;position:relative;"></div>
                    <script>
                    (function() {{
                        var tries = 0;
                        function initViewer() {{
                            var el = document.getElementById("{viewer_id}");
                            if (!el || typeof $3Dmol === "undefined") {{
                                tries += 1;
                                if (tries < 50) {{
                                    setTimeout(initViewer, 100);
                                }}
                                return;
                            }}
                            var viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor: "white"}});
                            var xyz = `{full_xyz}`;
                            viewer.addModelsAsFrames(xyz, "xyz");
                            viewer.setStyle({{}}, {{stick: {{radius: 0.1}}, sphere: {{scale: 0.25}}}});
                            viewer.zoomTo();
                            viewer.zoom({CALC_MOL_ZOOM});
                            viewer.setFrame(0);
                            viewer.render();
                            window.calc_trj_viewer = viewer;
                            window._calcMolViewer = viewer;
                            if (window.calcResizeMolViewer) setTimeout(window.calcResizeMolViewer, 200);
                        }}
                        setTimeout(initViewer, 0);
                    }})();
                    </script>
                    """
                    display(HTML(html_content))
                state['traj_viewer_ready'] = True
            else:
                _run_js(f"""
                setTimeout(function(){{
                    if (window.calc_trj_viewer) {{
                        window.calc_trj_viewer.setFrame({idx});
                        window.calc_trj_viewer.render();
                    }}
                }}, 0);
                """)
            return

        # Single frame: render with py3Dmol
        xyz_content = f"{n_atoms}\n{comment}\n{xyz_block}\n"
        _render_3dmol(xyz_content)

    # -- ORCA terminated check ----------------------------------------------
    def calc_orca_terminated_normally(path):
        try:
            with path.open('rb') as f:
                f.seek(0, 2)
                size = f.tell()
                f.seek(max(0, size - 20000))
                tail = f.read()
            text = tail.decode('utf-8', errors='ignore')
            return '****ORCA TERMINATED NORMALLY****' in text
        except Exception:
            return False

    # -- directory listing --------------------------------------------------
    def calc_list_directory():
        state['file_content'] = ''
        state['all_items'] = []
        state['search_spans'] = []
        state['current_match'] = -1
        calc_file_list.options = []
        calc_mol_viewer.clear_output()
        calc_mol_container.layout.display = 'none'
        calc_file_info.value = ''
        calc_path_display.value = ''
        calc_search_result.value = ''
        calc_folder_search.value = ''
        calc_prev_btn.disabled = True
        calc_next_btn.disabled = True
        _set_view_toggle(False, True)
        calc_copy_btn.disabled = True
        calc_copy_path_btn.disabled = True
        calc_update_view()
        calc_set_message('Select a file...')

        current_dir = _calc_dir() / state['current_path'] if state['current_path'] else _calc_dir()
        if not current_dir.exists():
            calc_file_list.options = ['(Folder not found)']
            return

        items = []
        try:
            sort_mode = calc_sort_dropdown.value
            if sort_mode == 'date_desc':
                entries = sorted(
                    current_dir.iterdir(),
                    key=lambda x: (not x.is_dir(), -x.stat().st_mtime if x.exists() else 0),
                )
            elif sort_mode == 'date_asc':
                entries = sorted(
                    current_dir.iterdir(),
                    key=lambda x: (not x.is_dir(), x.stat().st_mtime if x.exists() else 0),
                )
            else:  # name
                entries = sorted(
                    current_dir.iterdir(),
                    key=lambda x: (not x.is_dir(), x.name.lower()),
                )
            for entry in entries:
                if entry.is_dir():
                    items.append(f'📁 {entry.name}')
                else:
                    suffix = entry.suffix.lower()
                    if suffix == '.xyz':
                        items.append(f'🔬 {entry.name}')
                    elif suffix == '.png':
                        items.append(f'🖼 {entry.name}')
                    elif suffix in ['.out', '.log']:
                        is_delfin = entry.name.lower().startswith('delfin_')
                        ok = True if is_delfin else calc_orca_terminated_normally(entry)
                        icon = '📄' if ok else '❌'
                        items.append(f'{icon} {entry.name}')
                    elif suffix in ['.inp', '.sh']:
                        items.append(f'📝 {entry.name}')
                    elif suffix in ['.cube', '.cub']:
                        items.append(f'🧊 {entry.name}')
                    elif suffix in ['.gbw', '.cis', '.densities']:
                        items.append(f'💾 {entry.name}')
                    elif suffix in ['.doc', '.docx']:
                        items.append(f'📃 {entry.name}')
                    else:
                        items.append(f'📄 {entry.name}')
        except PermissionError:
            items = ['(Permission denied)']

        state['all_items'] = items if items else ['(Empty folder)']
        calc_file_list.options = state['all_items']
        calc_update_path_display()
        calc_update_report_btn()

    def calc_filter_file_list(change=None):
        query = calc_folder_search.value.strip().lower()
        if not query:
            calc_file_list.options = state['all_items']
        else:
            filtered = [item for item in state['all_items'] if query in item.lower()]
            calc_file_list.options = filtered if filtered else ['(No matches)']

    # -- recalc helpers -----------------------------------------------------
    def calc_reset_recalc_state():
        state['selected_inp_path'] = None
        state['selected_inp_base'] = ''
        state['recalc_active'] = False
        calc_recalc_btn.disabled = True
        calc_submit_recalc_btn.disabled = True
        calc_recalc_status.value = ''
        calc_edit_area.value = ''
        calc_edit_area.layout.display = 'none'
        calc_recalc_toolbar.layout.display = 'none'
        calc_content_label.value = '<b>📄 File Content:</b>'

    def calc_get_recalc_base_name(stem):
        match = re.match(r'^(.*)_recalc_(\d+)$', stem)
        return match.group(1) if match else stem

    def calc_next_recalc_index(job_dir, base_name):
        pattern = re.compile(rf'^{re.escape(base_name)}_recalc_(\d+)\.(?:inp|out)$')
        max_idx = 0
        for entry in job_dir.iterdir():
            m = pattern.match(entry.name)
            if m:
                max_idx = max(max_idx, int(m.group(1)))
        return max_idx + 1

    def calc_find_workspace_root(current_path):
        if not current_path:
            return None
        path_parts = current_path.split('/')
        for i in range(len(path_parts), 0, -1):
            test_path = '/'.join(path_parts[:i])
            if (_calc_dir() / test_path / 'CONTROL.txt').exists():
                return test_path
        return None

    # -- delete helpers -----------------------------------------------------
    def calc_delete_hide_confirm():
        state['delete_current'] = False
        calc_delete_confirm.layout.display = 'none'
        calc_delete_status.value = ''

    # -- search helpers -----------------------------------------------------
    def calc_pos_to_line_col(pos):
        if pos < 0:
            return 0, 0
        line = state['file_content'].count('\n', 0, pos) + 1
        last_nl = state['file_content'].rfind('\n', 0, pos)
        col = pos + 1 if last_nl == -1 else pos - last_nl
        return line, col

    def calc_update_nav_buttons():
        if not state['search_spans']:
            calc_prev_btn.disabled = True
            calc_next_btn.disabled = True
            return
        calc_prev_btn.disabled = (state['current_match'] <= 0)
        calc_next_btn.disabled = (state['current_match'] >= len(state['search_spans']) - 1)

    def calc_update_search_result():
        if not state['search_spans']:
            calc_search_result.value = ''
            return
        if state['current_match'] < 0:
            calc_search_result.value = (
                f'<span style="color:green;">{len(state["search_spans"])} matches</span>'
            )
            return
        start, _ = state['search_spans'][state['current_match']]
        line, col = calc_pos_to_line_col(start)
        calc_search_result.value = (
            f'<b>{state["current_match"] + 1}/{len(state["search_spans"])}</b> '
            f'<span style="color:#555;">(line {line}, col {col})</span>'
        )

    def calc_show_match():
        if not state['search_spans'] or state['current_match'] < 0:
            return
        calc_update_nav_buttons()
        calc_update_search_result()
        calc_apply_highlight(calc_search_input.value.strip(), state['current_match'])
        calc_scroll_to('match')

    # -- options dropdown ---------------------------------------------------
    def calc_update_options_dropdown():
        selected = calc_file_list.value
        sel_lower = selected.lower() if selected else ''
        if selected:
            name = selected[2:].strip()
            if _calc_is_complete_mutation_csv(name):
                calc_options_dropdown.options = ['(Options)', 'Preselection']
                calc_options_dropdown.value = '(Options)'
                calc_options_dropdown.layout.display = 'block'
                return
        if selected and 'OCCUPIER.txt' in selected:
            calc_options_dropdown.options = ['(Options)', 'Override']
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
        elif selected and ('CONTROL.txt' in selected or sel_lower.endswith('.inp')):
            calc_options_dropdown.options = ['(Options)', 'Recalc']
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'block'
        else:
            calc_options_dropdown.options = ['(Options)']
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.layout.display = 'none'
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            if not calc_view_toggle.value and not state['recalc_active']:
                calc_content_area.layout.display = 'block'

    # -- event handlers -----------------------------------------------------
    def calc_on_back(b):
        if state['current_path']:
            parts = state['current_path'].split('/')
            state['current_path'] = '/'.join(parts[:-1])
            calc_list_directory()

    def calc_on_home(b):
        state['current_path'] = ''
        calc_list_directory()

    def calc_on_refresh(b):
        calc_list_directory()

    def calc_go_top(b):
        if state['file_content']:
            calc_render_content(scroll_to='top')

    def calc_go_bottom(b):
        if state['file_content']:
            calc_render_content(scroll_to='bottom')

    def calc_do_search(b=None):
        query = calc_search_input.value.strip()
        if not query and calc_search_suggest.value and calc_search_suggest.value != '(Select)':
            query = calc_search_suggest.value
            calc_search_input.value = query
        state['search_spans'] = []
        state['current_match'] = -1

        if not query or not state['file_content']:
            calc_search_result.value = ''
            calc_update_nav_buttons()
            calc_apply_highlight('', -1)
            return

        pattern = re.compile(re.escape(query), re.IGNORECASE)
        state['search_spans'] = [m.span() for m in pattern.finditer(state['file_content'])]

        if not state['search_spans']:
            calc_search_result.value = '<span style="color:red;">0 matches</span>'
            calc_update_nav_buttons()
            calc_apply_highlight('', -1)
            return

        state['current_match'] = 0
        calc_update_nav_buttons()
        calc_update_search_result()
        calc_apply_highlight(query, state['current_match'])
        calc_scroll_to('match')

    def calc_on_suggest(change):
        value = change['new']
        if not value or value == '(Select)':
            return
        calc_search_input.value = value

    def calc_prev_match(b):
        if state['search_spans'] and state['current_match'] > 0:
            state['current_match'] -= 1
            calc_show_match()

    def calc_next_match(b):
        if state['search_spans'] and state['current_match'] < len(state['search_spans']) - 1:
            state['current_match'] += 1
            calc_show_match()

    def calc_on_xyz_input_change(change):
        new_val = change['new']
        if state['xyz_frames'] and 1 <= new_val <= len(state['xyz_frames']):
            state['xyz_current_frame'][0] = new_val - 1
            calc_update_xyz_viewer(initial_load=False)

    def calc_on_xyz_copy(button):
        frames = state['xyz_frames']
        if not frames:
            return
        idx = state['xyz_current_frame'][0]
        if idx < 0 or idx >= len(frames):
            return
        comment, xyz_block, n_atoms = frames[idx]
        xyz_content = f"{n_atoms}\n{comment}\n{xyz_block}"
        escaped = xyz_content.replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n')
        _run_js(
            f"navigator.clipboard.writeText('{escaped}')"
            f".then(() => console.log('Copied frame {idx+1}'))"
            f".catch(err => console.error('Copy failed:', err));"
        )

    def calc_on_coord_copy(button):
        if not state['file_content']:
            return
        escaped = state['file_content'].replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n')
        _run_js(
            f"navigator.clipboard.writeText('{escaped}')"
            f".then(() => console.log('Copied coord'))"
            f".catch(err => console.error('Copy failed:', err));"
        )

    def calc_on_content_copy(button):
        content = ''
        if state['recalc_active']:
            content = calc_edit_area.value or ''
        if not content:
            content = state['file_content'] or ''
        if not content:
            return
        escaped = content.replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n')
        _run_js(
            f"navigator.clipboard.writeText('{escaped}')"
            f".then(() => console.log('Copied file content'))"
            f".catch(err => console.error('Copy failed:', err));"
        )

    def calc_on_path_copy(button):
        selected = calc_file_list.value
        if not selected:
            return
        name = selected[2:].strip()
        if state['current_path']:
            full_path = str(_calc_dir() / state['current_path'] / name)
        else:
            full_path = str(_calc_dir() / name)
        calc_path_display.value = (
            f'<input type="text" value="{_html.escape(full_path)}" onclick="this.select()"'
            f' style="width:100%;font-family:monospace;font-size:12px;border:1px solid #aaa;'
            f'padding:2px;background:#f8f8f8" readonly>'
        )
        escaped = full_path.replace('\\', '\\\\').replace("'", "\\'")
        _run_js(
            f"navigator.clipboard.writeText('{escaped}')"
            f".then(() => console.log('Copied path'))"
            f".catch(err => console.error('Copy failed:', err));"
        )

    def calc_on_report(button):
        current_dir = _calc_dir() / state['current_path'] if state['current_path'] else _calc_dir()
        dir_key = str(current_dir)

        if dir_key in state['report_running']:
            calc_report_status.value = '<span style="color:orange;">Report already running for this folder...</span>'
            return

        state['report_running'][dir_key] = True
        calc_report_status.value = f'<span style="color:blue;">Generating report in {current_dir.name}...</span>'

        def run_report():
            try:
                result = subprocess.run(
                    ['delfin', '--report', 'docx'],
                    cwd=str(current_dir),
                    capture_output=True, text=True, timeout=600,
                )
                if result.returncode == 0:
                    calc_report_status.value = f'<span style="color:green;">Report generated in {current_dir.name}!</span>'
                else:
                    calc_report_status.value = (
                        f'<span style="color:red;">Error in {current_dir.name}: '
                        f'{_html.escape(result.stderr[:100])}</span>'
                    )
            except subprocess.TimeoutExpired:
                calc_report_status.value = f'<span style="color:red;">Timeout in {current_dir.name}</span>'
            except Exception as e:
                calc_report_status.value = f'<span style="color:red;">Error: {_html.escape(str(e)[:100])}</span>'
            finally:
                state['report_running'].pop(dir_key, None)

        threading.Thread(target=run_report, daemon=True).start()

    def calc_on_options_change(change):
        if change['new'] == 'Override':
            calc_override_input.layout.display = 'block'
            calc_override_time.layout.display = 'block'
            calc_override_btn.layout.display = 'block'
            calc_override_status.layout.display = 'block'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            calc_content_area.layout.display = 'block'
            selected = calc_file_list.value
            if selected:
                name = selected[2:].strip()
                stage_name = name.replace('_OCCUPIER.txt', '').replace('OCCUPIER.txt', '')
                calc_override_input.value = f'{stage_name}=' if stage_name else ''
        elif change['new'] == 'Recalc':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'block'
            calc_override_btn.layout.display = 'block'
            calc_override_status.layout.display = 'block'
            calc_override_status.value = ''
            calc_edit_area.value = state['file_content']
            calc_edit_area.layout.display = 'block'
            _calc_preselect_show(False)
            calc_content_area.layout.display = 'none'
        elif change['new'] == 'Preselection':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            csv_path = _calc_get_selected_path()
            if not csv_path or not csv_path.exists():
                calc_preselect_status.value = '<span style="color:#d32f2f;">File not found.</span>'
                _calc_preselect_show(True)
                return
            _calc_preselect_load(csv_path)
            _calc_preselect_show(True)
            _calc_preselect_render_current()
        else:
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            calc_content_area.layout.display = 'block'

    def calc_on_override_start(button):
        is_recalc = calc_options_dropdown.value == 'Recalc'
        selected = calc_file_list.value
        selected_name = selected[2:].strip() if selected and len(selected) > 2 else ''
        is_inp_recalc = bool(is_recalc and selected_name.lower().endswith('.inp'))

        if not is_recalc:
            override_value = calc_override_input.value.strip()
            if not override_value or '=' not in override_value:
                calc_override_status.value = (
                    '<span style="color:#d32f2f;">Enter STAGE=INDEX (e.g., red_step_2=1)</span>'
                )
                return

        # ORCA recalc for .inp selection
        if is_inp_recalc:
            inp_path = None
            if state['selected_inp_path'] and Path(state['selected_inp_path']).exists():
                inp_path = Path(state['selected_inp_path'])
            elif selected_name:
                inp_path = (
                    (_calc_dir() / state['current_path'] / selected_name)
                    if state['current_path']
                    else (_calc_dir() / selected_name)
                )
                if not inp_path.exists():
                    inp_path = None
            if inp_path is None:
                calc_override_status.value = '<span style="color:#d32f2f;">No .inp file selected.</span>'
                return

            inp_content = calc_edit_area.value or inp_path.read_text()
            if not inp_content.strip():
                calc_override_status.value = '<span style="color:#d32f2f;">Input content is empty.</span>'
                return

            job_dir = inp_path.parent
            base_stem = inp_path.stem
            base_name = calc_get_recalc_base_name(base_stem)
            next_idx = calc_next_recalc_index(job_dir, base_name)
            job_name = f"{base_name}_recalc_{next_idx}"

            try:
                new_inp = job_dir / f"{job_name}.inp"
                new_inp.write_text(inp_content)

                slurm_time = calc_override_time.value.strip() or '24:00:00'
                pal_used, maxcore_used = parse_inp_resources(inp_content)
                if pal_used is None and ctx.orca_pal_widget is not None:
                    pal_used = ctx.orca_pal_widget.value
                if maxcore_used is None and ctx.orca_maxcore_widget is not None:
                    maxcore_used = ctx.orca_maxcore_widget.value

                result = ctx.backend.submit_orca(
                    job_dir=job_dir,
                    job_name=job_name,
                    inp_file=f'{job_name}.inp',
                    time_limit=slurm_time,
                    pal=pal_used or 40,
                    maxcore=maxcore_used or 6000,
                )
            except Exception as exc:
                calc_override_status.value = (
                    f'<span style="color:#d32f2f;">Error: {_html.escape(str(exc))}</span>'
                )
                return

            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                calc_override_status.value = (
                    f'<span style="color:#2e7d32;">Submitted {_html.escape(job_name)}'
                    f' (ID: {_html.escape(job_id)})</span>'
                )
            else:
                msg = result.stderr or result.stdout or 'Unknown error'
                calc_override_status.value = f'<span style="color:#d32f2f;">{_html.escape(msg)}</span>'
            return

        # Find workspace root (where CONTROL.txt is located)
        if not state['current_path']:
            calc_override_status.value = '<span style="color:#d32f2f;">No folder selected</span>'
            return

        workspace_root = calc_find_workspace_root(state['current_path'])
        if not workspace_root:
            calc_override_status.value = (
                '<span style="color:#d32f2f;">No CONTROL.txt found in parent directories</span>'
            )
            return

        workspace_path = _calc_dir() / workspace_root

        # For Recalc: save edited CONTROL.txt first
        if is_recalc:
            try:
                control_path = workspace_path / 'CONTROL.txt'
                control_path.write_text(calc_edit_area.value)
                state['file_content'] = calc_edit_area.value
                calc_render_content(scroll_to='top')
            except Exception as e:
                calc_override_status.value = (
                    f'<span style="color:#d32f2f;">Error saving CONTROL.txt: {e}</span>'
                )
                return

        calc_override_status.value = '<span style="color:#1976d2;">Submitting job...</span>'

        try:
            # Derive PAL/maxcore
            pal_used = None
            maxcore_used = None
            control_path = workspace_path / 'CONTROL.txt'
            if control_path.exists():
                pal_used, maxcore_used = parse_resource_settings(control_path.read_text())

            if pal_used is None or maxcore_used is None:
                inp_candidate = None
                try:
                    if state['selected_inp_path'] and Path(state['selected_inp_path']).exists():
                        inp_candidate = Path(state['selected_inp_path'])
                except Exception:
                    inp_candidate = None

                if inp_candidate is None:
                    inp_files = sorted(workspace_path.glob('*.inp'))
                    inp_candidate = inp_files[0] if inp_files else None

                if inp_candidate is not None and inp_candidate.exists():
                    pal_inp, maxcore_inp = parse_inp_resources(inp_candidate.read_text())
                    if pal_used is None:
                        pal_used = pal_inp
                    if maxcore_used is None:
                        maxcore_used = maxcore_inp

            if pal_used is None or maxcore_used is None:
                calc_override_status.value = (
                    '<span style="color:#d32f2f;">PAL/maxcore not found in CONTROL.txt or .inp</span>'
                )
                return

            if is_recalc:
                job_name = f'recalc_{workspace_root.replace("/", "_")}'
                result = ctx.backend.submit_delfin(
                    job_dir=workspace_path,
                    job_name=job_name,
                    mode='delfin-recalc',
                    time_limit=calc_override_time.value,
                    pal=pal_used,
                    maxcore=maxcore_used,
                )
            else:
                override_value = calc_override_input.value.strip()
                job_name = f'occ_{override_value.replace("=", "_")}'
                result = ctx.backend.submit_delfin(
                    job_dir=workspace_path,
                    job_name=job_name,
                    mode='delfin-recalc-override',
                    time_limit=calc_override_time.value,
                    pal=pal_used,
                    maxcore=maxcore_used,
                    override=override_value,
                )

            if result.returncode == 0:
                calc_override_status.value = (
                    f'<span style="color:#388e3c;">Submitted from {workspace_root}: '
                    f'{result.stdout.strip()}</span>'
                )
            else:
                calc_override_status.value = (
                    f'<span style="color:#d32f2f;">Error: {result.stderr.strip()}</span>'
                )
        except Exception as e:
            calc_override_status.value = f'<span style="color:#d32f2f;">Error: {e}</span>'

    def calc_on_recalc_click(button):
        if not state['selected_inp_path'] or not state['selected_inp_path'].exists():
            calc_recalc_status.value = '<span style="color:#d32f2f;">No .inp file selected.</span>'
            return
        try:
            content = state['selected_inp_path'].read_text()
        except Exception as exc:
            calc_recalc_status.value = (
                f'<span style="color:#d32f2f;">Error reading file: {_html.escape(str(exc))}</span>'
            )
            return
        calc_edit_area.value = content
        state['recalc_active'] = True
        calc_submit_recalc_btn.disabled = False
        _set_view_toggle(False, False)
        calc_content_label.value = '<b>📄 File Content (editable):</b>'
        calc_update_view()

    def calc_on_submit_recalc(button):
        if not state['selected_inp_path'] or not state['selected_inp_path'].exists():
            calc_recalc_status.value = '<span style="color:#d32f2f;">No .inp file selected.</span>'
            return
        inp_content = calc_edit_area.value
        if not inp_content.strip():
            calc_recalc_status.value = '<span style="color:#d32f2f;">Input content is empty.</span>'
            return

        job_dir = state['selected_inp_path'].parent
        base_stem = state['selected_inp_path'].stem
        base_name = calc_get_recalc_base_name(base_stem)
        next_idx = calc_next_recalc_index(job_dir, base_name)
        job_name = f"{base_name}_recalc_{next_idx}"

        try:
            inp_path = job_dir / f"{job_name}.inp"
            inp_path.write_text(inp_content)

            slurm_time = calc_recalc_time.value.strip() or '24:00:00'
            pal_used, maxcore_used = parse_inp_resources(inp_content)
            if pal_used is None and ctx.orca_pal_widget is not None:
                pal_used = ctx.orca_pal_widget.value
            if maxcore_used is None and ctx.orca_maxcore_widget is not None:
                maxcore_used = ctx.orca_maxcore_widget.value

            result = ctx.backend.submit_orca(
                job_dir=job_dir,
                job_name=job_name,
                inp_file=f'{job_name}.inp',
                time_limit=slurm_time,
                pal=pal_used or 40,
                maxcore=maxcore_used or 6000,
            )
        except Exception as exc:
            calc_recalc_status.value = (
                f'<span style="color:#d32f2f;">Error: {_html.escape(str(exc))}</span>'
            )
            return

        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
            calc_recalc_status.value = (
                f'<span style="color:#2e7d32;">Submitted {_html.escape(job_name)}'
                f' (ID: {_html.escape(job_id)})</span>'
            )
        else:
            msg = result.stderr or result.stdout or 'Unknown error'
            calc_recalc_status.value = f'<span style="color:#d32f2f;">{_html.escape(msg)}</span>'

    def calc_on_delete_click(button):
        if not calc_file_list.value or calc_file_list.value.startswith('('):
            if state['current_path']:
                calc_delete_label.value = (
                    f'<b>Delete current folder?</b> '
                    f'<code>{_html.escape(state["current_path"])}</code>'
                )
                state['delete_current'] = True
                calc_delete_confirm.layout.display = 'flex'
            else:
                calc_delete_status.value = '<span style="color:#d32f2f;">No item selected.</span>'
            return
        calc_delete_label.value = '<b>Delete selected item?</b>'
        state['delete_current'] = False
        calc_delete_confirm.layout.display = 'flex'

    def calc_on_delete_yes(button):
        delete_current = state['delete_current']
        if delete_current:
            if not state['current_path']:
                calc_delete_status.value = '<span style="color:#d32f2f;">No current folder.</span>'
                calc_delete_confirm.layout.display = 'none'
                state['delete_current'] = False
                return
            full_path = _calc_dir() / state['current_path']
            name = state['current_path']
        else:
            selected = calc_file_list.value
            if not selected or selected.startswith('('):
                calc_delete_status.value = '<span style="color:#d32f2f;">No item selected.</span>'
                calc_delete_confirm.layout.display = 'none'
                state['delete_current'] = False
                return
            name = selected[2:].strip()
            full_path = (
                (_calc_dir() / state['current_path'] / name)
                if state['current_path']
                else (_calc_dir() / name)
            )
        try:
            if full_path.is_dir():
                shutil.rmtree(full_path)
            else:
                full_path.unlink()
            calc_delete_status.value = f'<span style="color:#2e7d32;">Deleted: {_html.escape(name)}</span>'
        except Exception as exc:
            calc_delete_status.value = (
                f'<span style="color:#d32f2f;">Delete failed: {_html.escape(str(exc))}</span>'
            )
        calc_delete_confirm.layout.display = 'none'
        state['delete_current'] = False
        if delete_current:
            cp = state['current_path']
            state['current_path'] = cp.rsplit('/', 1)[0] if '/' in cp else ''
        saved_filter = calc_folder_search.value
        calc_list_directory()
        if saved_filter:
            calc_folder_search.value = saved_filter
            calc_filter_file_list()
        calc_file_list.value = None

    def calc_on_delete_no(button):
        calc_delete_hide_confirm()

    def calc_on_sort_change(change):
        saved_filter = calc_folder_search.value
        calc_list_directory()
        if saved_filter:
            calc_folder_search.value = saved_filter
            calc_filter_file_list()

    # -- file selection handler (main) --------------------------------------
    def calc_on_select(change):
        selected = change['new']
        if not selected or selected.startswith('('):
            return

        name = selected[2:].strip()
        full_path = (
            (_calc_dir() / state['current_path'] / name)
            if state['current_path']
            else (_calc_dir() / name)
        )

        # Navigate into directory
        if full_path.is_dir():
            state['current_path'] = (
                f'{state["current_path"]}/{name}' if state['current_path'] else name
            )
            calc_list_directory()
            return

        # Reset search state
        state['search_spans'] = []
        state['current_match'] = -1
        calc_search_result.value = ''
        calc_search_input.value = ''
        calc_update_nav_buttons()
        calc_delete_hide_confirm()
        calc_update_options_dropdown()

        # Reset display (avoid flashing text area before viewer is ready)
        calc_mol_viewer.clear_output()
        calc_mol_container.layout.display = 'none'
        calc_content_area.layout.display = 'none'
        calc_content_label.layout.display = 'none'
        _set_view_toggle(False, True)
        calc_file_info.value = ''
        calc_path_display.value = ''
        state['file_content'] = ''
        calc_reset_recalc_state()
        _calc_preselect_show(False)

        # Reset xyz trajectory controls
        state['xyz_frames'].clear()
        state['xyz_current_frame'][0] = 0
        calc_xyz_controls.layout.display = 'none'
        calc_coord_controls.layout.display = 'none'
        calc_xyz_frame_label.value = ''
        calc_xyz_frame_input.value = 1
        calc_xyz_frame_input.max = 1
        calc_xyz_frame_total.value = '<b>/ 1</b>'

        if not full_path.exists():
            calc_set_message('File not found')
            return

        # File info
        size = full_path.stat().st_size
        if size > 1024 * 1024:
            size_str = f'{size / (1024 * 1024):.2f} MB'
        elif size > 1024:
            size_str = f'{size / 1024:.2f} KB'
        else:
            size_str = f'{size} B'

        suffix = full_path.suffix.lower()
        name_lower = full_path.name.lower()

        # --- Turbomole coord file ---
        if name_lower == 'coord':
            _set_view_toggle(True, False)
            calc_mol_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_update_view()
            try:
                content = full_path.read_text()
                state['file_content'] = content
                calc_render_content(scroll_to='top')
                calc_copy_btn.disabled = False
                calc_copy_path_btn.disabled = False
                xyz_data = coord_to_xyz(content)
                if xyz_data:
                    calc_file_info.value = (
                        f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                        f' ({size_str}, Turbomole coord)'
                    )
                    state['xyz_frames'].clear()
                    state['xyz_frames'].append(('TM coord', xyz_data.split('\n', 2)[2] if '\n' in xyz_data else '', xyz_data.count('\n') - 1))
                    state['xyz_current_frame'][0] = 0
                    calc_xyz_controls.layout.display = 'none'
                    calc_coord_controls.layout.display = 'flex'
                    _render_3dmol(xyz_data, stick_radius=0.15)
                else:
                    calc_file_info.value = f'<b>{_html.escape(name)}</b> (could not parse)'
            except Exception as e:
                calc_file_info.value = f'<b>{_html.escape(name)}</b> (error: {e})'
            return

        # --- XYZ file (with trajectory support) ---
        if suffix == '.xyz':
            _set_view_toggle(True, False)
            calc_mol_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_update_view()
            try:
                content = full_path.read_text()
                state['file_content'] = content
                calc_render_content(scroll_to='top')
                calc_copy_btn.disabled = False
                calc_copy_path_btn.disabled = False
                lines = content.split('\n')

                state['xyz_frames'].clear()
                state['xyz_frames'].extend(parse_xyz_frames(content))
                n_frames = len(state['xyz_frames'])

                if n_frames > 1:
                    calc_file_info.value = (
                        f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                        f' ({size_str}, {n_frames} frames)'
                    )
                    state['xyz_current_frame'][0] = 0
                    calc_xyz_frame_input.max = n_frames
                    calc_xyz_frame_input.value = 1
                    calc_xyz_frame_total.value = f"<b>/ {n_frames}</b>"
                    calc_xyz_controls.layout.display = 'flex'
                    calc_coord_controls.layout.display = 'none'
                    calc_update_xyz_viewer(initial_load=True)
                else:
                    calc_file_info.value = (
                        f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                        f' ({size_str}, {len(lines)} lines)'
                    )
                    if state['xyz_frames']:
                        state['xyz_current_frame'][0] = 0
                        calc_xyz_frame_input.max = len(state['xyz_frames'])
                        calc_xyz_frame_input.value = 1
                        calc_xyz_frame_total.value = f"<b>/ {len(state['xyz_frames'])}</b>"
                        calc_xyz_controls.layout.display = 'flex'
                        calc_coord_controls.layout.display = 'none'
                    else:
                        calc_xyz_controls.layout.display = 'none'
                    calc_xyz_frame_label.value = ''
                    _render_3dmol(content)
                calc_render_content(scroll_to='top')
            except Exception as e:
                calc_set_message(f'Error: {e}')
            return

        # --- PNG image ---
        if suffix == '.png':
            try:
                data = full_path.read_bytes()
                b64 = base64.b64encode(data).decode('ascii')
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str})'
                )
                calc_content_area.value = (
                    "<div style='border:1px solid #ddd; padding:6px; background:#fafafa;'>"
                    f"<img src='data:image/png;base64,{b64}' style='max-width:50%;"
                    f" max-height:{CALC_CONTENT_HEIGHT}px; height:auto; display:block;' />"
                    "</div>"
                )
                calc_update_view()
            except Exception as e:
                calc_set_message(f'Error: {e}')
            return

        # --- Cube/Cub volumetric data ---
        if suffix in ['.cube', '.cub']:
            _set_view_toggle(False, False)
            calc_update_view()
            try:
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str})'
                )
                lower_name = name.lower()
                if lower_name.endswith('.esp.cube') or lower_name.endswith('.esp.cub'):
                    from delfin.reporting.esp_report import generate_esp_png_for_state

                    workspace_root = None
                    search = full_path.parent
                    while True:
                        if (search / 'CONTROL.txt').exists():
                            workspace_root = search
                            break
                        if search.parent == search:
                            break
                        search = search.parent
                    if workspace_root is None:
                        workspace_root = full_path.parent

                    rendered = generate_esp_png_for_state(workspace_root, 'S0')
                    calc_set_message('ESP rendered via DELFIN report pipeline. Toggle Visualize to view.')
                    state['file_content'] = ''
                    with calc_mol_viewer:
                        if rendered and rendered.exists():
                            data = rendered.read_bytes()
                            b64 = base64.b64encode(data).decode('ascii')
                            display(HTML(
                                f"<img src='data:image/png;base64,{b64}'"
                                f" style='max-width:100%; max-height:{CALC_CONTENT_HEIGHT}px; height:auto;"
                                f" display:block;' />"
                            ))
                        else:
                            print('ESP render failed (no PNG generated).')
                else:
                    content = full_path.read_text()
                    state['file_content'] = ''
                    calc_set_message('Cube file loaded. Toggle Visualize to render isosurfaces.')
                    def _cube_extra(view, raw):
                        view.addVolumetricData(raw, 'cube', {'isoval': 0.02, 'color': '#0026ff', 'opacity': 0.85})
                        view.addVolumetricData(raw, 'cube', {'isoval': -0.02, 'color': '#b00010', 'opacity': 0.85})
                    _render_3dmol(content, fmt='cube', extra_fn=_cube_extra)
            except Exception as e:
                calc_set_message(f'Error: {e}')
            return

        # --- Binary files ---
        if suffix in ['.gbw', '.cis', '.densities', '.tmp']:
            calc_file_info.value = (
                f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                f' ({size_str})'
            )
            calc_set_message(f'Binary file ({size_str})\n\nCannot display binary content.')
            return

        # --- DOCX files ---
        if suffix in ['.doc', '.docx']:
            try:
                import mammoth

                def convert_image_to_base64(image):
                    try:
                        with image.open() as img_stream:
                            img_data = img_stream.read()
                            content_type = image.content_type if image.content_type else 'image/png'
                            b64 = base64.b64encode(img_data).decode('ascii')
                            return {'src': f'data:{content_type};base64,{b64}'}
                    except Exception:
                        return {'src': ''}

                with open(str(full_path), 'rb') as docx_file:
                    result = mammoth.convert_to_html(
                        docx_file,
                        convert_image=mammoth.images.img_element(convert_image_to_base64),
                    )
                    html_content = result.value

                styled_html = (
                    '<div style="font-family: Calibri, Arial, sans-serif; line-height: 1.5;'
                    ' padding: 15px; background: white;">'
                    '<style>'
                    '.docx-content img { max-width: 100%; height: auto; margin: 10px 0;'
                    ' border: 1px solid #ddd; border-radius: 4px; }'
                    '.docx-content p { margin: 8px 0; }'
                    '.docx-content h1 { color: #1976d2; font-size: 1.8em; margin: 20px 0 10px 0; }'
                    '.docx-content h2 { color: #1976d2; font-size: 1.5em; margin: 18px 0 8px 0; }'
                    '.docx-content h3 { color: #333; font-size: 1.2em; margin: 15px 0 6px 0; }'
                    '.docx-content table { border-collapse: collapse; margin: 10px 0; width: 100%; }'
                    '.docx-content td, .docx-content th { border: 1px solid #ddd; padding: 8px; }'
                    '.docx-content ul, .docx-content ol { margin: 10px 0; padding-left: 25px; }'
                    '</style>'
                    f'<div class="docx-content">{html_content}</div>'
                    '</div>'
                )

                text_only = re.sub(r'<[^>]+>', '', html_content)
                text_only = re.sub(r'\s+', ' ', text_only).strip()
                state['file_content'] = text_only

                num_images = html_content.count('<img ')
                img_info = f', {num_images} image{"s" if num_images != 1 else ""}' if num_images > 0 else ''

                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str}{img_info})'
                )
                calc_content_area.value = (
                    "<div id='calc-content-box' style='height:100%; overflow-y:auto;"
                    " overflow-x:hidden; border:1px solid #ddd; background:#fafafa;"
                    " width:100%; box-sizing:border-box;'>" + styled_html + "</div>"
                )
                calc_copy_btn.disabled = False
                calc_copy_path_btn.disabled = False
                calc_update_view()

            except ImportError:
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str})'
                )
                calc_set_message('mammoth not installed.\n\nRun: pip install mammoth')
            except Exception as e:
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str})'
                )
                calc_set_message(f'Error reading Word file: {e}')
            return

        # --- Text files (default) ---
        try:
            content = full_path.read_text()
            state['file_content'] = content
            calc_copy_btn.disabled = False
            calc_copy_path_btn.disabled = False
            lines = content.split('\n')
            calc_file_info.value = (
                f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                f' ({size_str}, {len(lines)} lines)'
            )
            calc_render_content(scroll_to='top')

            # input.txt: try to show molecule
            if name == 'input.txt':
                _set_view_toggle(True, False)
                calc_mol_container.layout.display = 'block'
                calc_content_area.layout.display = 'none'
                calc_update_view()
                try:
                    xyz_content = calc_build_xyz_from_input(
                        content, title=name, base_dir=full_path.parent
                    )
                    if xyz_content:
                        state['xyz_frames'].clear()
                        state['xyz_frames'].extend(parse_xyz_frames(xyz_content))
                        if state['xyz_frames']:
                            state['xyz_current_frame'][0] = 0
                            calc_xyz_frame_input.max = len(state['xyz_frames'])
                            calc_xyz_frame_input.value = 1
                            calc_xyz_frame_total.value = f"<b>/ {len(state['xyz_frames'])}</b>"
                            calc_xyz_controls.layout.display = 'flex'
                            calc_coord_controls.layout.display = 'none'
                        else:
                            calc_xyz_controls.layout.display = 'none'
                        calc_xyz_frame_label.value = ''
                        _render_3dmol(xyz_content)
                except Exception as exc:
                    calc_set_message(f'Error: {exc}')

            # .inp files: prepare visualization but do not auto-open it
            elif suffix == '.inp':
                try:
                    xyz_content = calc_build_xyz_from_input(
                        content, title=name, base_dir=full_path.parent
                    )
                    if xyz_content:
                        state['xyz_frames'].clear()
                        state['xyz_frames'].extend(parse_xyz_frames(xyz_content))
                        if state['xyz_frames']:
                            state['xyz_current_frame'][0] = 0
                            calc_xyz_frame_input.max = len(state['xyz_frames'])
                            calc_xyz_frame_input.value = 1
                            calc_xyz_frame_total.value = f"<b>/ {len(state['xyz_frames'])}</b>"
                            calc_xyz_controls.layout.display = 'flex'
                            calc_coord_controls.layout.display = 'none'
                        else:
                            calc_xyz_controls.layout.display = 'none'
                        calc_xyz_frame_label.value = ''
                        _render_3dmol(xyz_content)
                        _set_view_toggle(False, False)
                    else:
                        _set_view_toggle(False, True)
                except Exception as exc:
                    _set_view_toggle(False, True)
                    calc_set_message(f'Error: {exc}')
                calc_mol_container.layout.display = 'none'
                calc_content_area.layout.display = 'block'

            if suffix == '.inp':
                state['selected_inp_path'] = full_path
                state['selected_inp_base'] = full_path.stem
                calc_recalc_btn.disabled = False

            calc_update_view()
        except UnicodeDecodeError:
            calc_file_info.value = (
                f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                f' ({size_str})'
            )
            calc_set_message(f'Binary file ({size_str})\n\nCannot display.')
        except Exception as e:
            calc_set_message(f'Error: {e}')

    # -- wiring -------------------------------------------------------------
    calc_back_btn.on_click(calc_on_back)
    calc_home_btn.on_click(calc_on_home)
    calc_refresh_btn.on_click(calc_on_refresh)
    calc_top_btn.on_click(calc_go_top)
    calc_bottom_btn.on_click(calc_go_bottom)
    calc_search_btn.on_click(calc_do_search)
    calc_search_input.observe(lambda change: calc_do_search(), names='value')
    calc_search_suggest.observe(calc_on_suggest, names='value')
    calc_view_toggle.observe(_on_view_toggle, names='value')
    calc_prev_btn.on_click(calc_prev_match)
    calc_next_btn.on_click(calc_next_match)
    calc_recalc_btn.on_click(calc_on_recalc_click)
    calc_submit_recalc_btn.on_click(calc_on_submit_recalc)
    calc_delete_btn.on_click(calc_on_delete_click)
    calc_delete_yes_btn.on_click(calc_on_delete_yes)
    calc_delete_no_btn.on_click(calc_on_delete_no)
    calc_file_list.observe(calc_on_select, names='value')
    calc_folder_search.observe(calc_filter_file_list, names='value')
    calc_options_dropdown.observe(calc_on_options_change, names='value')
    calc_override_btn.on_click(calc_on_override_start)
    calc_preselect_yes.on_click(lambda _b: _calc_preselect_decide('yes'))
    calc_preselect_no.on_click(lambda _b: _calc_preselect_decide('no'))
    calc_sort_dropdown.observe(calc_on_sort_change, names='value')
    calc_xyz_frame_input.observe(calc_on_xyz_input_change, names='value')
    calc_xyz_copy_btn.on_click(calc_on_xyz_copy)
    calc_coord_copy_btn.on_click(calc_on_coord_copy)
    calc_copy_btn.on_click(calc_on_content_copy)
    calc_copy_path_btn.on_click(calc_on_path_copy)
    calc_report_btn.on_click(calc_on_report)

    # -- initialise ---------------------------------------------------------
    if ctx.calc_dir.exists():
        calc_list_directory()
    else:
        calc_file_list.options = ['(calc folder not found)']

    # -- layout -------------------------------------------------------------
    calc_nav_bar = widgets.VBox([
        calc_path_label,
        widgets.HBox(
            [calc_back_btn, calc_home_btn, calc_refresh_btn, calc_delete_btn],
            layout=widgets.Layout(
                width='100%', overflow_x='hidden',
                justify_content='flex-start', gap='6px',
            ),
        ),
    ], layout=widgets.Layout(width='100%', overflow_x='hidden'))

    calc_left = widgets.VBox([
        calc_nav_bar,
        calc_filter_sort_row,
        calc_file_list,
    ], layout=widgets.Layout(
        flex=f'0 0 {CALC_LEFT_DEFAULT}px',
        min_width=f'{CALC_LEFT_MIN}px',
        max_width=f'{CALC_LEFT_MAX}px',
        padding='5px', overflow_x='hidden', overflow_y='hidden',
    ))

    calc_css = widgets.HTML(
        '<style>'
        '#calc-content-box { overflow-x:hidden !important; }'
        '.calc-tab, .calc-tab * { overflow-x:hidden !important; box-sizing:border-box; }'
        '.calc-tab { overflow:hidden !important;'
        ' height:calc(100vh - 145px) !important;'
        ' max-height:calc(100vh - 145px) !important;'
        ' display:flex !important; flex-direction:column !important; }'
        '.calc-right { overflow:hidden !important;'
        ' display:flex !important; flex-direction:column !important; }'
        '.calc-left { display:flex !important; flex-direction:column !important; }'
        '.calc-left .widget-select { flex:1 1 0 !important; min-height:0 !important; }'
        '.calc-left .widget-select select { height:100% !important; }'
        '.calc-content-area { flex:1 1 0 !important; min-height:0 !important;'
        ' overflow-y:auto !important; overflow-x:hidden !important; }'
        '.calc-content-area .widget-html-content { height:100%; }'
        '.calc-edit-area textarea { height:100% !important; }'
        '.calc-tab .widget-vbox, .calc-tab .widget-hbox { overflow-y:hidden !important; }'
        '.calc-right .widget-output, .calc-right .jupyter-widgets-output-area { overflow:hidden !important; }'
        '.calc-right .widget-output .output_area { overflow:hidden !important; }'
        '.calc-right .widget-output .output { overflow:hidden !important; }'
        '.calc-right .output_scroll, .calc-right .output_subarea,'
        ' .calc-right .output_wrapper { overflow:hidden !important; }'
        '.calc-right .jp-OutputArea, .calc-right .jp-OutputArea-child,'
        ' .calc-right .jp-OutputArea-output { overflow:hidden !important; }'
        '.calc-right .p-Widget { overflow:hidden !important; }'
        '.calc-right .widget-textarea { overflow:visible !important; }'
        '.calc-right .widget-textarea textarea { overflow-y:auto !important; }'
        '.calc-left, .calc-right { overflow-x:hidden !important; }'
        '.calc-left * { max-width:100% !important; }'
        '.calc-right * { max-width:100% !important; }'
        '.calc-left .widget-select, .calc-left .widget-select select'
        ' { overflow-x:hidden !important; }'
        '.calc-left .widget-select select'
        ' { text-overflow: ellipsis; white-space: nowrap; }'
        '.calc-left .widget-text { flex:1 1 auto !important; min-width:0 !important; width:auto !important; }'
        '.calc-left .widget-text input { width:100% !important; }'
        '.calc-filter-row { width:100% !important; display:flex !important; }'
        '.calc-filter-row > .widget-text { flex:1 1 auto !important; min-width:0 !important; width:auto !important; }'
        '.calc-filter-row > .widget-text input { width:100% !important; min-width:0 !important; }'
        '.calc-filter-row > .widget-dropdown { flex:0 0 90px !important; margin-left:auto !important; }'
        '.calc-left .widget-text input'
        ' { overflow-x:hidden !important; overflow-y:hidden !important; text-overflow: ellipsis; }'
        '.calc-right .widget-text input'
        ' { overflow-x:hidden !important; overflow-y:hidden !important; }'
        '.calc-tab .widget-text { overflow:visible !important; }'
        '.calc-tab .widget-dropdown { height:26px !important; }'
        '.calc-tab .widget-dropdown select { '
        'height:26px !important; max-height:26px !important; line-height:26px !important; '
        'overflow:hidden !important; padding:0 6px !important; }'
        '.calc-tab .widget-dropdown, .calc-tab .widget-dropdown select { overflow:hidden !important; }'
        '.calc-options-dropdown { width:150px !important; min-width:150px !important; }'
        '.calc-options-dropdown select { width:150px !important; max-width:150px !important; }'
        '.calc-left .widget-text { width:auto !important; }'
        '.calc-left .widget-text input { width:100% !important; }'
        '.calc-filter-row > .widget-text { width:calc(100% - 94px) !important; }'
        '.calc-tab .widget-button { flex:0 0 auto !important; }'
        '.calc-tab .widget-dropdown, .calc-tab .widget-text { flex:0 0 auto !important; }'
        '.calc-tab input { overflow:hidden !important; height:26px !important;'
        ' line-height:26px !important; padding:0 6px !important; box-sizing:border-box !important; }'
        '.calc-tab input::-webkit-scrollbar { width:0; height:0; display:none; }'
        '.calc-tab input { scrollbar-width:none; }'
        '.calc-left .widget-vbox { overflow:hidden !important; }'
        '.calc-left code { display:block !important; overflow:hidden !important;'
        ' text-overflow:ellipsis !important; white-space:nowrap !important; }'
        '.calc-splitter { width:8px; height:100%; cursor:col-resize;'
        ' background:linear-gradient(to right, #d6d6d6, #f2f2f2, #d6d6d6);'
        ' border-radius:4px; display:block; }'
        '.calc-splitter:hover { background:linear-gradient('
        'to right, #d0d0d0, #f0f0f0, #d0d0d0); }'
        '.calc-mol-viewer { overflow:hidden !important; padding:0 !important; }'
        '.calc-mol-viewer .output_area, .calc-mol-viewer .output_subarea,'
        ' .calc-mol-viewer .output_wrapper, .calc-mol-viewer .jp-OutputArea-child,'
        ' .calc-mol-viewer .jp-OutputArea-output'
        ' { padding:0 !important; margin:0 !important; }'
        '</style>'
    )

    calc_right = widgets.VBox([
        widgets.HBox([
            calc_file_info,
            widgets.HBox(
                [calc_copy_path_btn, calc_copy_btn, calc_report_btn, calc_view_toggle],
                layout=widgets.Layout(
                    gap='10px',
                    flex_flow='row wrap',
                    justify_content='flex-end',
                    align_items='center',
                    width='100%',
                    overflow_x='hidden',
                ),
            ),
        ], layout=widgets.Layout(align_items='center', justify_content='space-between', width='100%')),
        calc_path_display,
        calc_mol_container,
        calc_content_label,
        calc_preselect_container,
        calc_delete_confirm,
        calc_delete_status,
        calc_recalc_toolbar,
        calc_content_toolbar,
        calc_override_status,
        calc_content_area,
        calc_edit_area,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='5px',
        overflow_x='hidden', overflow_y='hidden',
    ))

    calc_splitter = widgets.HTML(
        "<div class='calc-splitter' title='Drag to resize'></div>",
        layout=widgets.Layout(height='100%', width='6px', margin='0 2px'),
    )

    tab_widget = widgets.VBox([
        calc_css,
        widgets.HTML('<h3>📂 Calculations Browser</h3>'),
        widgets.HBox(
            [calc_left, calc_splitter, calc_right],
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
    tab_widget.add_class('calc-tab')
    calc_left.add_class('calc-left')
    calc_right.add_class('calc-right')
    calc_content_area.add_class('calc-content-area')

    # Enable draggable splitter (no visible slider)
    _run_js(f"""
    setTimeout(function() {{
        const left = document.querySelector('.calc-left');
        const right = document.querySelector('.calc-right');
        const splitter = document.querySelector('.calc-splitter');
        if (!left || !right || !splitter) return;
        if (splitter.dataset.bound === '1') return;
        splitter.dataset.bound = '1';

        const minW = {CALC_LEFT_MIN};
        const maxW = {CALC_LEFT_MAX};

        function onMove(e) {{
            const box = left.parentElement.getBoundingClientRect();
            let w = e.clientX - box.left;
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
    }}, 0);
    """)

    # Dynamic mol-viewer resize: keep square, flush with left panel bottom
    # Also monkey-patch $3Dmol.createViewer to capture viewer instances.
    _run_js("""
    setTimeout(function() {
        /* Monkey-patch $3Dmol.createViewer to store every new viewer */
        function patchCreateViewer() {
            if (typeof $3Dmol !== 'undefined' && !$3Dmol._patched) {
                var orig = $3Dmol.createViewer;
                $3Dmol.createViewer = function() {
                    var v = orig.apply(this, arguments);
                    window._calcMolViewer = v;
                    setTimeout(function() {
                        if (window.calcResizeMolViewer) window.calcResizeMolViewer();
                    }, 300);
                    return v;
                };
                $3Dmol._patched = true;
            }
        }
        patchCreateViewer();
        var patchInterval = setInterval(function() {
            patchCreateViewer();
            if (typeof $3Dmol !== 'undefined' && $3Dmol._patched) clearInterval(patchInterval);
        }, 500);

        window.calcResizeMolViewer = function() {
            var left = document.querySelector('.calc-left');
            var mv = document.querySelector('.calc-mol-viewer');
            if (!left || !mv || mv.offsetParent === null) return;
            var container = mv.closest('.widget-vbox');
            if (!container || container.style.display === 'none') return;
            var leftRect = left.getBoundingClientRect();
            var mvRect = mv.getBoundingClientRect();
            if (mvRect.top === 0 && mvRect.height === 0) return;
            var availH = leftRect.bottom - mvRect.top - 6;
            var availW = mv.parentElement.getBoundingClientRect().width - 4;
            var h = Math.floor(availH * 0.975);
            var w = Math.floor(Math.min(h * 1.2, availW));
            if (h < 200) h = 200;
            if (w < 240) w = 240;
            mv.style.width = w + 'px';
            mv.style.height = h + 'px';
            var inner = mv.querySelector('[id^="3dmolviewer"], [id^="calc_trj"]');
            if (inner) { inner.style.width = w + 'px'; inner.style.height = h + 'px'; }
            var v = window._calcMolViewer || window.calc_trj_viewer;
            if (v && typeof v.resize === 'function') { v.resize(); v.render(); }
        };
        window.addEventListener('resize', function() {
            setTimeout(window.calcResizeMolViewer, 100);
        });
        var tab = document.querySelector('.calc-tab');
        if (tab) {
            new MutationObserver(function() {
                setTimeout(window.calcResizeMolViewer, 200);
            }).observe(tab, {attributes: true, subtree: true, attributeFilter: ['style']});
        }
    }, 0);
    """)

    return tab_widget, {
        'calc_list_directory': calc_list_directory,
    }
