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
from IPython.display import HTML, clear_output, display

from .constants import CALC_SEARCH_OPTIONS
from .input_processing import (
    parse_inp_resources,
    parse_resource_settings,
    smiles_to_xyz,
    contains_metal,
    is_smiles,
)
from .helpers import disable_spellcheck
from .molecule_viewer import coord_to_xyz, parse_xyz_frames
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Chem.Draw import MolToImage


def create_tab(ctx):
    """Create the Calculations Browser tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- layout constants ---------------------------------------------------
    CALC_CONTENT_HEIGHT = 760
    CALC_MOL_SIZE = 460
    CALC_MOL_DYNAMIC_SCALE = 0.9725
    CALC_MOL_ZOOM = 0.9
    CALC_LEFT_DEFAULT = 320
    CALC_LEFT_MIN = 320
    CALC_LEFT_MAX = 520
    CALC_PRESELECT_VIZ_SIZE = 520
    # Large-file guardrails
    CALC_TEXT_FULL_READ_BYTES = 2 * 1024 * 1024
    CALC_TEXT_CHUNK_BYTES = 2 * 1024 * 1024
    CALC_SEARCH_MAX_MATCHES = 2000
    CALC_HIGHLIGHT_MAX_CHARS = 400_000
    # -- state (closure-captured) -------------------------------------------
    state = {
        'current_path': '',
        'file_content': '',
        'all_items': [],
        'search_spans': [],
        'current_match': -1,
        'search_truncated': False,
        'selected_file_path': None,
        'selected_file_size': 0,
        'selected_inp_path': None,
        'selected_inp_base': '',
        'recalc_active': False,
        'delete_current': False,
        'file_is_preview': False,
        'file_preview_note': '',
        'file_chunk_start': 0,
        'file_chunk_end': 0,
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
            'mode': 'complete_preselection',
            'molblock_cache': {},
            'regen_seed_counter': 0,
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
            width='100%',
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
    calc_chunk_prev_btn = widgets.Button(
        description='◀ Part',
        layout=widgets.Layout(width='85px', min_width='85px', height='26px'),
        disabled=True,
    )
    calc_chunk_prev_btn.add_class('calc-chunk-prev-trigger')
    calc_chunk_next_btn = widgets.Button(
        description='Part ▶',
        layout=widgets.Layout(width='85px', min_width='85px', height='26px'),
        disabled=True,
    )
    calc_chunk_next_btn.add_class('calc-chunk-next-trigger')
    calc_chunk_label = widgets.HTML(
        value='',
        layout=widgets.Layout(width='190px', min_width='190px'),
    )
    calc_chunk_request_start = widgets.IntText(
        value=0,
        layout=widgets.Layout(width='1px', height='1px'),
    )
    calc_chunk_request_start.add_class('calc-chunk-start-input')
    calc_chunk_request_ratio = widgets.FloatText(
        value=0.0,
        layout=widgets.Layout(width='1px', height='1px'),
    )
    calc_chunk_request_ratio.add_class('calc-chunk-ratio-input')
    calc_chunk_request_btn = widgets.Button(
        description='chunk-load',
        layout=widgets.Layout(width='1px', height='1px'),
    )
    calc_chunk_request_btn.add_class('calc-chunk-load-trigger')

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
    calc_preselect_img = widgets.Output(layout=widgets.Layout(
        width='100%',
        height='auto',
        margin='0', padding='0', border='1px solid #ccc', overflow='hidden',
        display='flex', align_items='center', justify_content='center',
    ))
    calc_preselect_img.add_class('calc-preselect-viz')
    calc_preselect_3d = widgets.Output(layout=widgets.Layout(
        width='100%',
        height='auto',
        margin='0', padding='0', border='1px solid #ccc', overflow='hidden',
        display='flex', align_items='center', justify_content='center',
    ))
    calc_preselect_3d.add_class('calc-preselect-viz')
    calc_preselect_img_wrap = widgets.Box(
        [calc_preselect_img],
        layout=widgets.Layout(
            flex='1 1 220px', min_width='180px', max_width=f'{CALC_PRESELECT_VIZ_SIZE}px',
            width='100%', overflow='hidden',
        ),
    )
    calc_preselect_img_wrap.add_class('calc-preselect-viz-wrap')
    calc_preselect_3d_wrap = widgets.Box(
        [calc_preselect_3d],
        layout=widgets.Layout(
            flex='1 1 220px', min_width='180px', max_width=f'{CALC_PRESELECT_VIZ_SIZE}px',
            width='100%', overflow='hidden',
        ),
    )
    calc_preselect_3d_wrap.add_class('calc-preselect-viz-wrap')
    calc_preselect_yes = widgets.Button(
        description='Yes', button_style='success',
        layout=widgets.Layout(width='80px', height='30px'),
    )
    calc_preselect_no = widgets.Button(
        description='No', button_style='danger',
        layout=widgets.Layout(width='80px', height='30px'),
    )
    calc_preselect_prev = widgets.Button(
        description='◀ Prev', button_style='',
        layout=widgets.Layout(width='90px', height='30px'),
    )
    calc_preselect_next = widgets.Button(
        description='Next ▶', button_style='',
        layout=widgets.Layout(width='90px', height='30px'),
    )
    calc_preselect_move = widgets.Button(
        description='Move', button_style='warning',
        layout=widgets.Layout(width='190px', height='30px', display='none'),
    )
    calc_preselect_close = widgets.Button(
        description='Close', button_style='',
        layout=widgets.Layout(width='90px', height='30px'),
    )
    calc_preselect_new3d = widgets.Button(
        description='New 3D Structure', button_style='info',
        layout=widgets.Layout(width='150px', height='30px'),
    )
    calc_preselect_status = widgets.HTML('')
    calc_preselect_container = widgets.VBox([
        calc_preselect_title,
        calc_preselect_progress,
        calc_preselect_info,
        widgets.HBox(
            [calc_preselect_img_wrap, calc_preselect_3d_wrap],
            layout=widgets.Layout(gap='20px', overflow_x='hidden', flex_flow='row wrap'),
        ),
        widgets.HBox(
            [
                calc_preselect_prev, calc_preselect_next, calc_preselect_close,
                calc_preselect_new3d, calc_preselect_move, calc_preselect_yes, calc_preselect_no,
            ],
            layout=widgets.Layout(gap='10px', flex_flow='row wrap'),
        ),
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
    calc_chunk_hidden_row = widgets.HBox(
        [
            calc_chunk_prev_btn, calc_chunk_next_btn, calc_chunk_label,
            calc_chunk_request_start, calc_chunk_request_ratio, calc_chunk_request_btn,
        ],
        layout=widgets.Layout(
            display='flex',
            height='0px',
            min_height='0px',
            overflow='hidden',
            visibility='hidden',
        ),
    )

    # -- helper closures ----------------------------------------------------
    RDLogger.DisableLog('rdApp.*')

    def _calc_parse_complex_mol(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol, sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            ))
            return mol
        except Exception:
            pass
        try:
            mol.UpdatePropertyCache(strict=False)
            return mol
        except Exception:
            pass
        return Chem.MolFromSmiles(smiles, sanitize=False)

    def _calc_preselect_smiles_to_3d_molblock(smiles, seed=42):
        try:
            seed = int(seed)
        except Exception:
            seed = 42
        is_complex = contains_metal(smiles)
        mol = None
        if is_complex:
            mol = _calc_parse_complex_mol(smiles)
            if mol is None:
                return None
            try:
                mol = Chem.AddHs(mol, addCoords=False)
            except Exception:
                pass
            try:
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                params.randomSeed = seed
                params.maxIterations = 500
                conf_id = AllChem.EmbedMolecule(mol, params)
                if conf_id < 0:
                    AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=seed)
            except Exception:
                try:
                    rdDepictor.Compute2DCoords(mol)
                    conf = mol.GetConformer()
                    for i in range(mol.GetNumAtoms()):
                        pos = conf.GetAtomPosition(i)
                        conf.SetAtomPosition(i, (pos.x, pos.y, 0.0))
                except Exception:
                    return None
        else:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is None:
                    return None
            try:
                mol = Chem.AddHs(mol)
                params = AllChem.ETKDGv3() if hasattr(AllChem, 'ETKDGv3') else AllChem.ETKDG()
                params.useRandomCoords = True
                params.randomSeed = seed
                conf_id = AllChem.EmbedMolecule(mol, params)
                if conf_id < 0:
                    AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=seed)
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass
        try:
            return Chem.MolToMolBlock(mol)
        except Exception:
            return None

    def _calc_is_complete_mutation_csv(selected_name):
        return selected_name and selected_name.lower() == 'complete_mutation_space.csv'

    def _calc_is_preselected_mutation_csv(selected_name):
        return selected_name and selected_name.lower() == 'preselected_mutation_space.csv'

    def _calc_is_rejected_mutation_csv(selected_name):
        return selected_name and selected_name.lower() == 'rejected_mutation_space.csv'

    def _calc_is_mutation_csv(selected_name):
        return (
            _calc_is_complete_mutation_csv(selected_name)
            or _calc_is_preselected_mutation_csv(selected_name)
            or _calc_is_rejected_mutation_csv(selected_name)
        )

    def _calc_mode_for_mutation_csv_name(selected_name):
        if _calc_is_complete_mutation_csv(selected_name):
            return 'complete_visualize'
        if _calc_is_preselected_mutation_csv(selected_name):
            return 'preselected_visualize'
        if _calc_is_rejected_mutation_csv(selected_name):
            return 'rejected_visualize'
        return None

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

    def _calc_write_mutation_csv(csv_path, entries):
        lines = ['Label;SMILES']
        for label, smiles in entries:
            lines.append(f'{label};{smiles}')
        csv_path.write_text('\n'.join(lines), encoding='utf-8')

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

    def _calc_clamp_preselect_index(index, total):
        if total <= 0:
            return 0
        try:
            idx = int(index)
        except Exception:
            idx = 0
        return max(0, min(idx, total - 1))

    def _calc_preselect_write_outputs():
        entries = state['preselect']['entries']
        decisions = state['preselect']['decisions']
        job_dir_raw = state['preselect'].get('job_dir', '')
        if not job_dir_raw:
            return
        job_dir = Path(job_dir_raw)
        preselected_path = job_dir / 'preselected_mutation_space.csv'
        rejected_path = job_dir / 'rejected_mutation_space.csv'
        yes_entries = []
        no_entries = []
        for i, entry in enumerate(entries):
            decision = decisions.get(str(i))
            if decision == 'yes':
                yes_entries.append(entry)
            elif decision == 'no':
                no_entries.append(entry)
        _calc_write_mutation_csv(preselected_path, yes_entries)
        _calc_write_mutation_csv(rejected_path, no_entries)

    def _calc_preselect_save_state():
        job_dir_raw = state['preselect'].get('job_dir', '')
        if not job_dir_raw:
            return
        job_dir = Path(job_dir_raw)
        state_path = job_dir / 'selection_state.json'
        payload = {
            'csv_path': state['preselect']['csv_path'],
            'index': state['preselect']['index'],
            'decisions': state['preselect']['decisions'],
        }
        state_path.write_text(json.dumps(payload, indent=2), encoding='utf-8')

    def _calc_sync_complete_decision_for_entry(job_dir, entry, decision):
        complete_path = job_dir / 'complete_mutation_space.csv'
        if not complete_path.exists():
            return
        complete_entries = _calc_parse_mutation_csv(complete_path)
        if not complete_entries:
            return
        matched_indices = [i for i, e in enumerate(complete_entries) if e == entry]
        if not matched_indices:
            return
        state_path = job_dir / 'selection_state.json'
        payload = {'csv_path': str(complete_path), 'index': 0, 'decisions': {}}
        if state_path.exists():
            try:
                data = json.loads(state_path.read_text(encoding='utf-8', errors='ignore'))
                if isinstance(data, dict):
                    payload.update(data)
            except Exception:
                pass
        decisions = payload.get('decisions', {})
        if not isinstance(decisions, dict):
            decisions = {}
        for i in matched_indices:
            decisions[str(i)] = decision
        payload['csv_path'] = str(complete_path)
        payload['decisions'] = decisions
        if not isinstance(payload.get('index'), int):
            payload['index'] = 0
        state_path.write_text(json.dumps(payload, indent=2), encoding='utf-8')

        if (
            state['preselect'].get('job_dir')
            and Path(state['preselect']['job_dir']) == job_dir
            and str(state['preselect'].get('mode', '')).startswith('complete_')
        ):
            for i in matched_indices:
                state['preselect']['decisions'][str(i)] = decision

    def _calc_preselect_title_for_mode(mode, csv_name):
        if mode == 'complete_preselection':
            return f'<b>Preselection:</b> {csv_name}'
        if mode == 'complete_visualize':
            return f'<b>Visualize:</b> {csv_name}'
        if mode == 'preselected_visualize':
            return f'<b>Visualize (Preselected):</b> {csv_name}'
        if mode == 'rejected_visualize':
            return f'<b>Visualize (Rejected):</b> {csv_name}'
        return f'<b>Visualize:</b> {csv_name}'

    def _calc_preselect_load(csv_path, mode='complete_preselection'):
        entries = _calc_parse_mutation_csv(csv_path)
        job_dir = csv_path.parent
        decisions = {}
        index = 0

        if mode in ('complete_preselection', 'complete_visualize'):
            state_path = job_dir / 'selection_state.json'
            if state_path.exists():
                try:
                    data = json.loads(state_path.read_text(encoding='utf-8', errors='ignore'))
                    if data.get('csv_path') == str(csv_path):
                        decisions = data.get('decisions', {}) or {}
                        index = int(data.get('index', 0))
                except Exception:
                    decisions = {}
                    index = 0

            if not decisions:
                yes_set = _calc_read_decision_file(job_dir / 'preselected_mutation_space.csv')
                no_set = _calc_read_decision_file(job_dir / 'rejected_mutation_space.csv')
                for i, entry in enumerate(entries):
                    if entry in yes_set:
                        decisions[str(i)] = 'yes'
                    elif entry in no_set:
                        decisions[str(i)] = 'no'

            if mode == 'complete_preselection':
                index = _calc_next_preselect_index(index)
            else:
                index = _calc_clamp_preselect_index(index, len(entries))
        else:
            index = 0

        state['preselect']['entries'] = entries
        state['preselect']['decisions'] = decisions
        state['preselect']['csv_path'] = str(csv_path)
        state['preselect']['job_dir'] = str(job_dir)
        state['preselect']['mode'] = mode
        state['preselect']['index'] = index
        state['preselect']['molblock_cache'] = {}
        state['preselect']['regen_seed_counter'] = 0
        calc_preselect_title.value = _calc_preselect_title_for_mode(mode, csv_path.name)

    def _calc_preselect_cache_key(index):
        return f"{state['preselect'].get('csv_path', '')}::{int(index)}"

    def _calc_preselect_get_or_build_molblock(index, smiles, force_new=False):
        cache = state['preselect'].setdefault('molblock_cache', {})
        key = _calc_preselect_cache_key(index)
        if not force_new and key in cache:
            return cache.get(key)
        seed = 42
        if force_new:
            state['preselect']['regen_seed_counter'] = state['preselect'].get('regen_seed_counter', 0) + 1
            seed += state['preselect']['regen_seed_counter']
        mol_block = _calc_preselect_smiles_to_3d_molblock(smiles, seed=seed)
        if mol_block is not None:
            cache[key] = mol_block
        return mol_block

    def _calc_preselect_render_3d_molblock(mol_block):
        with calc_preselect_3d:
            clear_output(wait=True)
            if mol_block is None:
                display(HTML('<i>3D conversion failed</i>'))
                return

            _mol3d_counter[0] += 1
            viewer_id = f"calc_preselect_3dmol_{_mol3d_counter[0]}"
            mol_json = json.dumps(mol_block)
            display(HTML(f"""
                <div id="{viewer_id}" style="width:100%;height:100%;position:relative;"></div>
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
                        var box = el ? el.closest('.calc-preselect-viz') : null;
                        if (!el || typeof $3Dmol === "undefined"
                            || !box || box.offsetParent === null) {{
                            tries += 1;
                            if (tries < 80) setTimeout(initViewer, 50);
                            return;
                        }}
                        var rect = box.getBoundingClientRect();
                        var side = Math.floor(Math.min(rect.width || 0, rect.height || 0));
                        if (side >= 100) {{
                            el.style.width = side + 'px';
                            el.style.height = side + 'px';
                        }}
                        var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                        var molData = {mol_json};
                        viewer.addModel(molData, "mol");
                        viewer.setStyle({{}}, {{
                            stick: {{
                                colorscheme: "Jmol",
                                radius: 0.11,
                                singleBonds: false,
                                doubleBondScaling: 0.65,
                                tripleBondScaling: 0.65
                            }},
                            sphere: {{colorscheme: "Jmol", scale: 0.28}}
                        }});
                        viewer.zoomTo();
                        viewer.center();
                        viewer.zoom(0.90);
                        viewer.render();
                        window._calcPreselect3dViewer = viewer;
                        window._calcPreselect3dElId = "{viewer_id}";
                        if (window.calcResizePreselect3D) {{
                            setTimeout(window.calcResizePreselect3D, 120);
                        }}
                    }}
                    setTimeout(initViewer, 0);
                }})();
                </script>
            """))

    def _calc_preselect_render_current():
        mode = state['preselect'].get('mode', 'complete_preselection')
        is_complete_preselection = (mode == 'complete_preselection')
        is_preselected_visualize = (mode == 'preselected_visualize')
        is_rejected_visualize = (mode == 'rejected_visualize')
        entries = state['preselect']['entries']
        idx = state['preselect']['index']
        total = len(entries)

        calc_preselect_yes.layout.display = 'block' if is_complete_preselection else 'none'
        calc_preselect_no.layout.display = 'block' if is_complete_preselection else 'none'
        calc_preselect_move.layout.display = (
            'block' if (is_preselected_visualize or is_rejected_visualize) else 'none'
        )

        if is_preselected_visualize:
            calc_preselect_move.description = 'Move to Rejected'
            calc_preselect_move.button_style = 'danger'
        elif is_rejected_visualize:
            calc_preselect_move.description = 'Revive to Preselected'
            calc_preselect_move.button_style = 'success'
        else:
            calc_preselect_move.description = 'Move'
            calc_preselect_move.button_style = 'warning'

        if total == 0:
            calc_preselect_progress.value = '<b>No entries found.</b>'
            calc_preselect_info.value = ''
            calc_preselect_status.value = '<span style="color:#d32f2f;">Empty or invalid CSV.</span>'
            calc_preselect_prev.disabled = True
            calc_preselect_next.disabled = True
            calc_preselect_yes.disabled = True
            calc_preselect_no.disabled = True
            calc_preselect_move.disabled = True
            calc_preselect_new3d.disabled = True
            with calc_preselect_img:
                clear_output(wait=True)
            with calc_preselect_3d:
                clear_output(wait=True)
            return

        if idx >= total:
            calc_preselect_progress.value = f'<b>Done.</b> ({total} / {total})'
            calc_preselect_info.value = ''
            calc_preselect_status.value = '<span style="color:green;">All entries processed.</span>'
            calc_preselect_prev.disabled = False
            calc_preselect_next.disabled = True
            calc_preselect_yes.disabled = True
            calc_preselect_no.disabled = True
            calc_preselect_move.disabled = True
            calc_preselect_new3d.disabled = True
            with calc_preselect_img:
                clear_output(wait=True)
            with calc_preselect_3d:
                clear_output(wait=True)
            return

        label, smiles = entries[idx]
        decision = state['preselect']['decisions'].get(str(idx))
        decision_html = ''
        if decision == 'yes':
            decision_html = '<div><b>Decision:</b> <span style="color:#2e7d32;">preselected</span></div>'
        elif decision == 'no':
            decision_html = '<div><b>Decision:</b> <span style="color:#c62828;">rejected</span></div>'

        calc_preselect_progress.value = f'<b>{idx + 1} / {total}</b>'
        calc_preselect_info.value = (
            f'<div><b>Label:</b> {_html.escape(label)}</div>'
            f'<div><b>SMILES:</b> <code>{_html.escape(smiles)}</code></div>'
            f'{decision_html}'
        )

        if is_complete_preselection:
            calc_preselect_status.value = '<span style="color:#555;">Use Yes/No to classify.</span>'
        elif mode == 'complete_visualize':
            calc_preselect_status.value = '<span style="color:#555;">Visualize-only mode.</span>'
        elif is_preselected_visualize:
            calc_preselect_status.value = '<span style="color:#555;">Move entries to rejected if needed.</span>'
        elif is_rejected_visualize:
            calc_preselect_status.value = '<span style="color:#555;">Revive entries back to preselected if needed.</span>'
        else:
            calc_preselect_status.value = ''

        calc_preselect_prev.disabled = (idx <= 0)
        calc_preselect_next.disabled = (idx >= total - 1)
        calc_preselect_yes.disabled = not is_complete_preselection
        calc_preselect_no.disabled = not is_complete_preselection
        calc_preselect_move.disabled = not (is_preselected_visualize or is_rejected_visualize)
        calc_preselect_new3d.disabled = False

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
                img = MolToImage(mol, size=(CALC_PRESELECT_VIZ_SIZE, CALC_PRESELECT_VIZ_SIZE))
                display(img)

        mol_block = _calc_preselect_get_or_build_molblock(idx, smiles, force_new=False)
        _calc_preselect_render_3d_molblock(mol_block)

    def _calc_preselect_new_3d_structure(_btn=None):
        entries = state['preselect']['entries']
        idx = state['preselect']['index']
        if idx < 0 or idx >= len(entries):
            return
        _label, smiles = entries[idx]
        mol_block = _calc_preselect_get_or_build_molblock(idx, smiles, force_new=True)
        if mol_block is None:
            calc_preselect_status.value = '<span style="color:#d32f2f;">New 3D generation failed.</span>'
            return
        _calc_preselect_render_3d_molblock(mol_block)
        calc_preselect_status.value = '<span style="color:#555;">New 3D structure generated.</span>'

    def _calc_preselect_show(show):
        if show:
            calc_preselect_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_content_label.layout.display = 'none'
            calc_mol_container.layout.display = 'none'
            calc_edit_area.layout.display = 'none'
            # In preselect/visualize mode the content toolbar (Top/End/Search) is redundant.
            calc_content_toolbar.layout.display = 'none'
            calc_recalc_toolbar.layout.display = 'none'
        else:
            calc_preselect_container.layout.display = 'none'
            # Restore toolbars based on current mode.
            if state['recalc_active']:
                calc_content_toolbar.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'flex'
            elif calc_view_toggle.value:
                calc_content_toolbar.layout.display = 'none'
                calc_recalc_toolbar.layout.display = 'none'
            else:
                calc_content_toolbar.layout.display = 'flex'
                calc_recalc_toolbar.layout.display = 'none'

    def _calc_preselect_prev_entry(_btn=None):
        entries = state['preselect']['entries']
        total = len(entries)
        if total == 0:
            return
        idx = state['preselect']['index']
        if idx >= total:
            state['preselect']['index'] = total - 1
        elif idx > 0:
            state['preselect']['index'] = idx - 1
        _calc_preselect_render_current()

    def _calc_preselect_next_entry(_btn=None):
        entries = state['preselect']['entries']
        total = len(entries)
        if total == 0:
            return
        idx = state['preselect']['index']
        if idx < total - 1:
            state['preselect']['index'] = idx + 1
        _calc_preselect_render_current()

    def _calc_preselect_move_current(_btn=None):
        mode = state['preselect'].get('mode', '')
        if mode not in ('preselected_visualize', 'rejected_visualize'):
            return
        entries = state['preselect']['entries']
        idx = state['preselect']['index']
        if idx < 0 or idx >= len(entries):
            return
        job_dir_raw = state['preselect'].get('job_dir', '')
        if not job_dir_raw:
            return
        job_dir = Path(job_dir_raw)
        entry = entries[idx]

        if mode == 'preselected_visualize':
            source_path = job_dir / 'preselected_mutation_space.csv'
            target_path = job_dir / 'rejected_mutation_space.csv'
            target_decision = 'no'
            status_msg = '<span style="color:#c62828;">Moved to rejected.</span>'
        else:
            source_path = job_dir / 'rejected_mutation_space.csv'
            target_path = job_dir / 'preselected_mutation_space.csv'
            target_decision = 'yes'
            status_msg = '<span style="color:#2e7d32;">Revived to preselected.</span>'

        source_entries = _calc_parse_mutation_csv(source_path)
        target_entries = _calc_parse_mutation_csv(target_path)
        removed = False
        new_source = []
        for item in source_entries:
            if not removed and item == entry:
                removed = True
                continue
            new_source.append(item)
        if not removed:
            new_source = [e for i, e in enumerate(entries) if i != idx]
        if entry not in target_entries:
            target_entries.append(entry)

        _calc_write_mutation_csv(source_path, new_source)
        _calc_write_mutation_csv(target_path, target_entries)
        _calc_sync_complete_decision_for_entry(job_dir, entry, target_decision)

        selected_path = _calc_get_selected_path()
        if selected_path and selected_path.exists():
            try:
                state['file_content'] = selected_path.read_text(errors='ignore')
                calc_render_content(scroll_to='top')
            except Exception:
                pass

        state['preselect']['entries'] = new_source
        if state['preselect']['index'] >= len(new_source):
            state['preselect']['index'] = max(0, len(new_source) - 1)
        _calc_preselect_render_current()
        calc_preselect_status.value = status_msg

    def _calc_preselect_close_view(_btn=None):
        if calc_options_dropdown.value != '(Options)':
            calc_options_dropdown.value = '(Options)'
            return
        _set_view_toggle(False, False)
        _calc_preselect_show(False)
        calc_update_view()

    def _calc_preselect_decide(decision):
        if state['preselect'].get('mode', '') != 'complete_preselection':
            return
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

    _mol3d_counter = [0]

    def _render_3dmol(data, fmt='xyz', stick_radius=0.1, sphere_scale=0.25,
                      extra_fn=None):
        """Render a 3D molecule via JS with correct initial sizing."""
        import json
        _mol3d_counter[0] += 1
        viewer_id = f"mol3d_{_mol3d_counter[0]}"
        wrapper_id = f"calc_mol_wrap_{_mol3d_counter[0]}"
        data_json = json.dumps(data)
        view_scope_json = json.dumps(state.get('current_path') or '/')

        volumetric_js = ""
        if fmt == 'cube':
            volumetric_js = (
                "viewer.addVolumetricData(molData,'cube',"
                "{isoval:0.02,color:'#0026ff',opacity:0.85});"
                "viewer.addVolumetricData(molData,'cube',"
                "{isoval:-0.02,color:'#b00010',opacity:0.85});"
            )

        html = f"""
        <div id="{wrapper_id}" class="calc-mol-stage-wrapper" style="width:100%;">
            <div id="{viewer_id}" style="width:100%;height:{CALC_MOL_SIZE}px;position:relative;"></div>
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
                var mv = el ? el.closest('.calc-mol-viewer') : null;
                if (!el || typeof $3Dmol === "undefined"
                    || !mv || mv.offsetParent === null) {{
                    tries += 1;
                    if (tries < 80) setTimeout(initViewer, 50);
                    return;
                }}
                /* Compute correct size from layout before creating viewer */
                var lft = document.querySelector('.calc-left');
                if (lft) {{
                    var leftRect = lft.getBoundingClientRect();
                    var mvRect = mv.getBoundingClientRect();
                    if (mvRect.top > 0 || mvRect.height > 0) {{
                        var availH = leftRect.bottom - mvRect.top - 6;
                        var availW = mv.parentElement
                            ? mv.parentElement.getBoundingClientRect().width - 4
                            : mvRect.width;
                        var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
                        var w = Math.floor(Math.min(h * 1.2, availW));
                        if (h >= 200 && w >= 240) {{
                            mv.style.width = w + 'px';
                            mv.style.height = h + 'px';
                            el.style.width = w + 'px';
                            el.style.height = h + 'px';
                        }}
                    }}
                }}
                window._calcMolViewStateByScope = window._calcMolViewStateByScope || {{}};
                var viewScope = {view_scope_json};
                var previousViewer = window._calcMolViewer || window.calc_trj_viewer;
                var previousScope = window._calcMolViewScopeKey || viewScope;
                if (previousViewer && typeof previousViewer.getView === 'function') {{
                    try {{
                        window._calcMolViewStateByScope[previousScope] = previousViewer.getView();
                    }} catch (_e) {{}}
                }}
                var savedView = window._calcMolViewStateByScope[viewScope] || null;
                var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                var molData = {data_json};
                viewer.addModel(molData, "{fmt}");
                viewer.setStyle({{}}, {{
                    stick: {{radius: {stick_radius}}},
                    sphere: {{scale: {sphere_scale}}}
                }});
                {volumetric_js}
                if (savedView && typeof viewer.setView === 'function') {{
                    try {{
                        viewer.setView(savedView);
                    }} catch (_e) {{
                        viewer.zoomTo();
                        viewer.zoom({CALC_MOL_ZOOM});
                    }}
                }} else {{
                    viewer.zoomTo();
                    viewer.zoom({CALC_MOL_ZOOM});
                }}
                viewer.render();
                window._calcMolViewer = viewer;
                window._calcMolViewScopeKey = viewScope;
                var wrappers = document.querySelectorAll('.calc-mol-stage-wrapper');
                wrappers.forEach(function(w) {{
                    if (w.id !== "{wrapper_id}") w.remove();
                }});
                if (window.calcResizeMolViewer) {{
                    setTimeout(window.calcResizeMolViewer, 200);
                }}
            }}
            setTimeout(initViewer, 0);
        }})();
        </script>
        """
        with calc_mol_viewer:
            display(HTML(html))

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

    def _calc_open_mutation_visualize_from_selected():
        selected = calc_file_list.value
        if not selected or selected.startswith('('):
            return False
        name = selected[2:].strip()
        mode = _calc_mode_for_mutation_csv_name(name)
        if not mode:
            if name.lower() != 'input.txt':
                return False
            input_path = _calc_get_selected_path()
            if not input_path or not input_path.exists():
                return False
            input_text = state.get('file_content', '')
            if not input_text:
                try:
                    input_text = input_path.read_text(errors='ignore')
                except Exception:
                    input_text = ''
            if not is_smiles(input_text):
                return False
            smiles_lines = [ln.strip() for ln in input_text.splitlines() if ln.strip()]
            if not smiles_lines:
                return False
            state['preselect']['entries'] = [('input_smiles', smiles_lines[0])]
            state['preselect']['decisions'] = {}
            state['preselect']['csv_path'] = str(input_path)
            state['preselect']['job_dir'] = str(input_path.parent)
            state['preselect']['mode'] = 'complete_visualize'
            state['preselect']['index'] = 0
            state['preselect']['molblock_cache'] = {}
            state['preselect']['regen_seed_counter'] = 0
            calc_preselect_title.value = '<b>Visualize:</b> input.txt (SMILES)'
            _calc_preselect_show(True)
            _calc_preselect_render_current()
            return True
        csv_path = _calc_get_selected_path()
        if not csv_path or not csv_path.exists():
            calc_preselect_status.value = '<span style="color:#d32f2f;">File not found.</span>'
            _calc_preselect_show(True)
            return True
        _calc_preselect_load(csv_path, mode=mode)
        _calc_preselect_show(True)
        _calc_preselect_render_current()
        return True

    def _on_view_toggle(change=None):
        if calc_view_toggle.value:
            if _calc_open_mutation_visualize_from_selected():
                return
        else:
            if calc_preselect_container.layout.display != 'none':
                _calc_preselect_show(False)
        calc_update_view()

    def _calc_dir():
        return ctx.calc_dir

    # -- view / display helpers ---------------------------------------------
    def _calc_format_bytes(n_bytes):
        if n_bytes >= 1024 * 1024:
            return f'{n_bytes / (1024 * 1024):.2f} MB'
        if n_bytes >= 1024:
            return f'{n_bytes / 1024:.2f} KB'
        return f'{n_bytes} B'

    def _calc_hide_chunk_controls():
        calc_chunk_prev_btn.disabled = True
        calc_chunk_next_btn.disabled = True
        calc_chunk_label.value = ''

    def _calc_update_chunk_controls():
        size = int(state.get('selected_file_size') or 0)
        if (
            not state.get('file_is_preview')
            or not state.get('selected_file_path')
            or size <= CALC_TEXT_FULL_READ_BYTES
        ):
            _calc_hide_chunk_controls()
            return
        start = int(state.get('file_chunk_start') or 0)
        end = int(state.get('file_chunk_end') or 0)
        total_parts = max(1, (size + CALC_TEXT_CHUNK_BYTES - 1) // CALC_TEXT_CHUNK_BYTES)
        part_idx = min(total_parts, max(1, (start // CALC_TEXT_CHUNK_BYTES) + 1))
        calc_chunk_prev_btn.disabled = (start <= 0)
        calc_chunk_next_btn.disabled = (end >= size or end <= start)
        calc_chunk_label.value = f"<span style='display:none;'>Part {part_idx}/{total_parts}</span>"

    def _calc_read_text_chunk(path, size_bytes, start_byte):
        if size_bytes <= 0:
            return '', 0, 0
        safe_start = max(0, min(int(start_byte), max(0, size_bytes - 1)))
        with path.open('rb') as f:
            f.seek(safe_start)
            raw = f.read(CALC_TEXT_CHUNK_BYTES)
        safe_end = safe_start + len(raw)
        return raw.decode('utf-8', errors='ignore'), safe_start, safe_end

    def _calc_load_text_preview_chunk(path, size_bytes, start_byte, rerun_search=False):
        content, chunk_start, chunk_end = _calc_read_text_chunk(path, size_bytes, start_byte)
        state['file_content'] = content
        state['file_is_preview'] = True
        state['file_chunk_start'] = chunk_start
        state['file_chunk_end'] = chunk_end
        state['file_preview_note'] = ''
        calc_render_content(scroll_to=None)
        _calc_update_chunk_controls()
        query = (calc_search_input.value or '').strip()
        if rerun_search and query:
            calc_do_search()
        elif not query:
            state['search_spans'] = []
            state['current_match'] = -1
            state['search_truncated'] = False
            calc_search_result.value = ''
            calc_update_nav_buttons()

    def _calc_should_highlight():
        return len(state.get('file_content', '')) <= CALC_HIGHLIGHT_MAX_CHARS

    def _calc_is_chunk_mode():
        return bool(
            state.get('file_is_preview')
            and int(state.get('selected_file_size') or 0) > CALC_TEXT_FULL_READ_BYTES
            and state.get('selected_file_path')
        )

    def _calc_find_global_matches(path, query, max_matches):
        q_bytes = query.encode('utf-8', errors='ignore')
        if not q_bytes:
            return []
        q_lower = q_bytes.lower()
        q_len = len(q_lower)
        overlap = max(q_len - 1, 256)
        spans = []
        carry = b''
        file_offset = 0
        with path.open('rb') as f:
            while True:
                block = f.read(CALC_TEXT_CHUNK_BYTES)
                if not block:
                    break
                data = carry + block
                lower = data.lower()
                search_from = 0
                while True:
                    idx = lower.find(q_lower, search_from)
                    if idx < 0:
                        break
                    abs_start = file_offset - len(carry) + idx
                    abs_end = abs_start + q_len
                    if abs_start >= 0:
                        spans.append((abs_start, abs_end))
                        if len(spans) >= max_matches:
                            return spans
                    search_from = idx + 1
                if len(data) > overlap:
                    carry = data[-overlap:]
                else:
                    carry = data
                file_offset += len(block)
        return spans

    def _calc_line_col_from_file_offset(path, byte_pos):
        target = max(0, int(byte_pos))
        line = 1
        last_nl = -1
        seen = 0
        with path.open('rb') as f:
            while seen < target:
                need = min(CALC_TEXT_CHUNK_BYTES, target - seen)
                if need <= 0:
                    break
                chunk = f.read(need)
                if not chunk:
                    break
                idx = 0
                while True:
                    nl = chunk.find(b'\n', idx)
                    if nl < 0:
                        break
                    line += 1
                    last_nl = seen + nl
                    idx = nl + 1
                seen += len(chunk)
        col = target + 1 if last_nl < 0 else target - last_nl
        return line, col

    def _calc_local_text_span_from_file_bytes(path, chunk_start, abs_start, abs_end):
        chunk_start_i = max(0, int(chunk_start))
        start_i = max(0, int(abs_start))
        end_i = max(start_i + 1, int(abs_end))
        local_b_start = max(0, start_i - chunk_start_i)
        local_b_end = max(local_b_start + 1, end_i - chunk_start_i)
        try:
            with path.open('rb') as f:
                f.seek(chunk_start_i)
                raw = f.read(local_b_end)
            if not raw:
                return 0, 1
            local_b_start = min(local_b_start, len(raw))
            local_b_end = min(max(local_b_start + 1, local_b_end), len(raw))
            local_start = len(raw[:local_b_start].decode('utf-8', errors='ignore'))
            local_end = len(raw[:local_b_end].decode('utf-8', errors='ignore'))
            if local_end <= local_start:
                local_end = local_start + 1
            return local_start, local_end
        except Exception:
            fallback_start = max(0, start_i - chunk_start_i)
            fallback_end = max(fallback_start + 1, end_i - chunk_start_i)
            return fallback_start, fallback_end

    def _calc_apply_highlight_span(local_start, local_end):
        js = r"""
        setTimeout(function() {
            const box = document.getElementById('calc-content-box');
            const el = document.getElementById('calc-content-text');
            if (!el) return;
            const text = el.textContent || '';
            const s = __START__;
            const e = __END__;
            if (s < 0 || e <= s || s >= text.length) return;
            function esc(v) {
                return v.replace(/[&<>"]/g, function(c) {
                    return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];
                });
            }
            const end = Math.min(e, text.length);
            el.innerHTML =
                esc(text.slice(0, s))
                + '<mark class="calc-match current" id="calc-current-match">'
                + esc(text.slice(s, end))
                + '</mark>'
                + esc(text.slice(end));
            const mark = document.getElementById('calc-current-match');
            if (!box || !mark) return;
            const boxRect = box.getBoundingClientRect();
            const markRect = mark.getBoundingClientRect();
            const delta = (markRect.top - boxRect.top) - (box.clientHeight / 2);
            window.__calcChunkIgnoreScrollUntil = Date.now() + 350;
            window.__calcChunkProgrammaticScroll = true;
            box.scrollTop += delta;
            setTimeout(function() {
                window.__calcChunkProgrammaticScroll = false;
            }, 120);
        }, 0);
        """
        js = js.replace('__START__', str(int(local_start))).replace('__END__', str(int(local_end)))
        _run_js(js)

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
        chunk_mode = bool(
            state.get('file_is_preview')
            and int(state.get('selected_file_size') or 0) > CALC_TEXT_FULL_READ_BYTES
        )
        top_spacer_html = ''
        bottom_spacer_html = ''
        if chunk_mode:
            total_size = max(1, int(state.get('selected_file_size') or 1))
            chunk_start = max(0, int(state.get('file_chunk_start') or 0))
            chunk_end = max(chunk_start, int(state.get('file_chunk_end') or chunk_start))
            virtual_h = max(12000, min(180000, int(total_size / 96)))
            top_px = int((chunk_start / total_size) * virtual_h)
            bottom_px = int((max(0, total_size - chunk_end) / total_size) * virtual_h)
            top_spacer_html = f"<div id='calc-chunk-top-spacer' style='height:{top_px}px;'></div>"
            bottom_spacer_html = f"<div id='calc-chunk-bottom-spacer' style='height:{bottom_px}px;'></div>"
        calc_content_area.value = (
            "<style>"
            ".calc-match { background: #fff59d; padding: 0 2px; }"
            ".calc-match.current { background: #ffcc80; }"
            "</style>"
            "<div id='calc-content-box' style='height:100%;"
            " overflow-y:auto; overflow-x:hidden; border:1px solid #ddd; padding:6px;"
            " background:#fafafa; width:100%; box-sizing:border-box;'>"
            f"{top_spacer_html}"
            "<div id='calc-content-text' style='white-space:pre-wrap; overflow-wrap:anywhere;"
            " word-break:break-word; font-family:monospace; font-size:12px; line-height:1.3;'>"
            f"{_html.escape(text)}"
            "</div>"
            f"{bottom_spacer_html}"
            "</div>"
        )
        if chunk_mode:
            file_size_js = int(state.get('selected_file_size') or 0)
            chunk_start_js = int(state.get('file_chunk_start') or 0)
            _run_js("""
            setTimeout(function() {
                const box = document.getElementById('calc-content-box');
                const topSpacer = document.getElementById('calc-chunk-top-spacer');
                const txt = document.getElementById('calc-content-text');
                if (!box || !topSpacer || !txt) return;
                const startInput = document.querySelector('.calc-chunk-start-input input');
                const ratioInput = document.querySelector('.calc-chunk-ratio-input input');
                if (!startInput || !ratioInput) return;

                const fileSize = __FILE_SIZE__;
                const chunkBytes = __CHUNK_BYTES__;
                const currentStart = __CURRENT_START__;
                if (fileSize <= 0 || chunkBytes <= 0) return;

                function emit(el, value) {
                    el.value = String(value);
                    el.dispatchEvent(new Event('input', {bubbles: true}));
                    el.dispatchEvent(new Event('change', {bubbles: true}));
                }

                function requestChunkForRatio(ratio) {
                    if (!isFinite(ratio)) return;
                    const clamped = Math.max(0, Math.min(1, ratio));
                    if ((window.__calcChunkIgnoreScrollUntil || 0) > Date.now()) {
                        window.__calcChunkPendingRatio = clamped;
                        return;
                    }
                    if (window.__calcChunkProgrammaticScroll) return;
                    if (window.__calcChunkBusy) {
                        window.__calcChunkPendingRatio = clamped;
                        return;
                    }
                    let desired = Math.floor((clamped * fileSize) / chunkBytes) * chunkBytes;
                    if (desired < 0) desired = 0;
                    if (desired >= fileSize) desired = Math.max(0, fileSize - chunkBytes);
                    desired = Math.floor(desired / chunkBytes) * chunkBytes;
                    if (desired === currentStart) return;
                    window.__calcChunkBusy = true;
                    window.__calcChunkPendingRatio = null;
                    window.__calcChunkRequestedRatio = clamped;
                    emit(startInput, desired);
                    emit(ratioInput, clamped);
                    setTimeout(function() {
                        if (window.__calcChunkBusy) window.__calcChunkBusy = false;
                    }, 2000);
                }

                if (!box.dataset.chunkBound) {
                    box.dataset.chunkBound = '1';
                    let timer = null;
                    box.addEventListener('scroll', function() {
                        if (window.__calcChunkProgrammaticScroll) return;
                        if ((window.__calcChunkIgnoreScrollUntil || 0) > Date.now()) return;
                        if (timer) clearTimeout(timer);
                        timer = setTimeout(function() {
                            const maxScroll = Math.max(1, box.scrollHeight - box.clientHeight);
                            const ratio = box.scrollTop / maxScroll;
                            requestChunkForRatio(ratio);
                        }, 80);
                    }, {passive: true});
                }

                const maxScroll = Math.max(0, box.scrollHeight - box.clientHeight);
                const requestedRatio = window.__calcChunkRequestedRatio;
                window.__calcChunkProgrammaticScroll = true;
                if (typeof requestedRatio === 'number' && isFinite(requestedRatio)) {
                    const clamped = Math.max(0, Math.min(1, requestedRatio));
                    box.scrollTop = Math.floor(clamped * maxScroll);
                } else {
                    box.scrollTop = Math.min(maxScroll, topSpacer.offsetHeight || 0);
                }
                window.__calcChunkRequestedRatio = null;
                setTimeout(function() {
                    window.__calcChunkProgrammaticScroll = false;
                }, 120);
                setTimeout(function() {
                    window.__calcChunkBusy = false;
                    const pending = window.__calcChunkPendingRatio;
                    if (typeof pending === 'number' && isFinite(pending)) {
                        window.__calcChunkPendingRatio = null;
                        requestChunkForRatio(pending);
                    }
                }, 0);
            }, 0);
            """.replace('__FILE_SIZE__', str(file_size_js))
              .replace('__CHUNK_BYTES__', str(CALC_TEXT_CHUNK_BYTES))
              .replace('__CURRENT_START__', str(chunk_start_js)))
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
        # input.txt may contain only a SMILES string (optionally with short comment lines).
        # Convert to XYZ so it can be visualized like regular coordinates.
        if is_smiles(text):
            smiles_line = lines[0].strip()
            xyz_string, num_atoms, _method, error = smiles_to_xyz(smiles_line)
            if not error and xyz_string:
                return f"{num_atoms}\n{title}\n{xyz_string}"
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
                _mol3d_counter[0] += 1
                viewer_id = f"calc_trj_viewer_{_mol3d_counter[0]}"
                wrapper_id = f"calc_mol_wrap_{_mol3d_counter[0]}"
                view_scope_json = json.dumps(state.get('current_path') or '/')
                with calc_mol_viewer:
                    html_content = f"""
                    <div id="{wrapper_id}" class="calc-mol-stage-wrapper" style="width:100%;">
                        <div id="{viewer_id}" style="width:100%;height:{CALC_MOL_SIZE}px;position:relative;"></div>
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
                            var mv = el ? el.closest('.calc-mol-viewer') : null;
                            if (!el || typeof $3Dmol === "undefined"
                                || !mv || mv.offsetParent === null) {{
                                tries += 1;
                                if (tries < 80) setTimeout(initViewer, 50);
                                return;
                            }}
                            /* Pre-size to correct dimensions */
                            var lft = document.querySelector('.calc-left');
                            if (lft) {{
                                var leftRect = lft.getBoundingClientRect();
                                var mvRect = mv.getBoundingClientRect();
                                if (mvRect.top > 0 || mvRect.height > 0) {{
                                    var availH = leftRect.bottom - mvRect.top - 6;
                                    var availW = mv.parentElement
                                        ? mv.parentElement.getBoundingClientRect().width - 4
                                        : mvRect.width;
                                    var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
                                    var w = Math.floor(Math.min(h * 1.2, availW));
                                    if (h >= 200 && w >= 240) {{
                                        mv.style.width = w + 'px';
                                        mv.style.height = h + 'px';
                                        el.style.width = w + 'px';
                                        el.style.height = h + 'px';
                                    }}
                                }}
                            }}
                            window._calcMolViewStateByScope = window._calcMolViewStateByScope || {{}};
                            var viewScope = {view_scope_json};
                            var previousViewer = window._calcMolViewer || window.calc_trj_viewer;
                            var previousScope = window._calcMolViewScopeKey || viewScope;
                            if (previousViewer && typeof previousViewer.getView === 'function') {{
                                try {{
                                    window._calcMolViewStateByScope[previousScope] = previousViewer.getView();
                                }} catch (_e) {{}}
                            }}
                            var savedView = window._calcMolViewStateByScope[viewScope] || null;
                            var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                            var xyz = `{full_xyz}`;
                            viewer.addModelsAsFrames(xyz, "xyz");
                            viewer.setStyle({{}}, {{stick: {{radius: 0.1}}, sphere: {{scale: 0.25}}}});
                            if (savedView && typeof viewer.setView === 'function') {{
                                try {{
                                    viewer.setView(savedView);
                                }} catch (_e) {{
                                    viewer.zoomTo();
                                    viewer.zoom({CALC_MOL_ZOOM});
                                }}
                            }} else {{
                                viewer.zoomTo();
                                viewer.zoom({CALC_MOL_ZOOM});
                            }}
                            viewer.setFrame({idx});
                            viewer.render();
                            window.calc_trj_viewer = viewer;
                            window._calcMolViewer = viewer;
                            window._calcMolViewScopeKey = viewScope;
                            var wrappers = document.querySelectorAll('.calc-mol-stage-wrapper');
                            wrappers.forEach(function(w) {{
                                if (w.id !== "{wrapper_id}") w.remove();
                            }});
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
        state['search_truncated'] = False
        state['selected_file_path'] = None
        state['selected_file_size'] = 0
        state['file_is_preview'] = False
        state['file_preview_note'] = ''
        state['file_chunk_start'] = 0
        state['file_chunk_end'] = 0
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
        _calc_hide_chunk_controls()
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
        notes = []
        if state.get('search_truncated'):
            notes.append(f'limited to {CALC_SEARCH_MAX_MATCHES} matches')
        note_html = ''
        if notes:
            note_html = f' <span style="color:#777;">({", ".join(notes)})</span>'
        if state['current_match'] < 0:
            calc_search_result.value = (
                f'<span style="color:green;">{len(state["search_spans"])} matches</span>{note_html}'
            )
            return
        start, _ = state['search_spans'][state['current_match']]
        if _calc_is_chunk_mode():
            path_str = state.get('selected_file_path')
            if path_str and Path(path_str).exists():
                line, col = _calc_line_col_from_file_offset(Path(path_str), start)
            else:
                line, col = 0, 0
        else:
            line, col = calc_pos_to_line_col(start)
        calc_search_result.value = (
            f'<b>{state["current_match"] + 1}/{len(state["search_spans"])}</b> '
            f'<span style="color:#555;">(line {line}, col {col})</span>{note_html}'
        )

    def calc_show_match():
        if not state['search_spans'] or state['current_match'] < 0:
            return
        calc_update_nav_buttons()
        calc_update_search_result()
        if _calc_is_chunk_mode():
            start, end = state['search_spans'][state['current_match']]
            chunk_start = int(state.get('file_chunk_start') or 0)
            chunk_end = int(state.get('file_chunk_end') or 0)
            if start < chunk_start or end > chunk_end or not state.get('file_content'):
                _run_js(
                    "window.__calcChunkRequestedRatio = null;"
                    "window.__calcChunkPendingRatio = null;"
                    "window.__calcChunkBusy = false;"
                    "window.__calcChunkIgnoreScrollUntil = Date.now() + 1200;"
                    "window.__calcChunkProgrammaticScroll = false;"
                )
                desired = max(0, int(start) - (CALC_TEXT_CHUNK_BYTES // 4))
                _calc_request_chunk_start(desired)
                chunk_start = int(state.get('file_chunk_start') or 0)
                chunk_end = int(state.get('file_chunk_end') or 0)
                if start < chunk_start or end > chunk_end:
                    _calc_request_chunk_start(start)
                    chunk_start = int(state.get('file_chunk_start') or 0)
            path_str = state.get('selected_file_path')
            path = Path(path_str) if path_str else None
            if path and path.exists():
                local_start, local_end = _calc_local_text_span_from_file_bytes(
                    path, chunk_start, start, end
                )
            else:
                local_start = max(0, int(start) - chunk_start)
                local_end = max(local_start + 1, int(end) - chunk_start)
            _calc_apply_highlight_span(local_start, local_end)
            return
        if _calc_should_highlight():
            calc_apply_highlight(calc_search_input.value.strip(), state['current_match'])
            calc_scroll_to('match')

    # -- options dropdown ---------------------------------------------------
    def calc_update_options_dropdown():
        selected = calc_file_list.value
        sel_lower = selected.lower() if selected else ''
        if selected:
            name = selected[2:].strip()
            if _calc_is_complete_mutation_csv(name):
                calc_options_dropdown.options = ['(Options)', 'Preselection', 'Visualize']
                calc_options_dropdown.value = '(Options)'
                calc_options_dropdown.layout.display = 'block'
                return
            if _calc_is_preselected_mutation_csv(name) or _calc_is_rejected_mutation_csv(name):
                calc_options_dropdown.options = ['(Options)', 'Visualize']
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
        state['search_truncated'] = False

        if not query or not state['file_content']:
            calc_search_result.value = ''
            calc_update_nav_buttons()
            if _calc_is_chunk_mode():
                _run_js("""
                (function() {
                    const el = document.getElementById('calc-content-text');
                    if (!el) return;
                    el.textContent = el.textContent || '';
                })();
                """)
            elif _calc_should_highlight():
                calc_apply_highlight('', -1)
            return

        if _calc_is_chunk_mode():
            _run_js(
                "window.__calcChunkRequestedRatio = null;"
                "window.__calcChunkPendingRatio = null;"
                "window.__calcChunkBusy = false;"
                "window.__calcChunkIgnoreScrollUntil = Date.now() + 1200;"
                "window.__calcChunkProgrammaticScroll = false;"
            )
            path_str = state.get('selected_file_path')
            path = Path(path_str) if path_str else None
            if not path or not path.exists():
                calc_search_result.value = '<span style="color:red;">0 matches</span>'
                calc_update_nav_buttons()
                return
            spans = _calc_find_global_matches(path, query, CALC_SEARCH_MAX_MATCHES + 1)
            if len(spans) > CALC_SEARCH_MAX_MATCHES:
                state['search_truncated'] = True
                spans = spans[:CALC_SEARCH_MAX_MATCHES]
            state['search_spans'] = spans
            if not state['search_spans']:
                calc_search_result.value = '<span style="color:red;">0 matches</span>'
                calc_update_nav_buttons()
                _run_js("""
                (function() {
                    const el = document.getElementById('calc-content-text');
                    if (!el) return;
                    el.textContent = el.textContent || '';
                })();
                """)
                return
            state['current_match'] = 0
            calc_show_match()
            return

        pattern = re.compile(re.escape(query), re.IGNORECASE)
        spans = []
        for m in pattern.finditer(state['file_content']):
            spans.append(m.span())
            if len(spans) >= CALC_SEARCH_MAX_MATCHES:
                state['search_truncated'] = True
                break
        state['search_spans'] = spans

        if not state['search_spans']:
            calc_search_result.value = '<span style="color:red;">0 matches</span>'
            calc_update_nav_buttons()
            if _calc_should_highlight():
                calc_apply_highlight('', -1)
            return

        state['current_match'] = 0
        calc_update_nav_buttons()
        calc_update_search_result()
        if _calc_should_highlight():
            calc_apply_highlight(query, state['current_match'])
            calc_scroll_to('match')

    def calc_on_suggest(change):
        value = change['new']
        if not value or value == '(Select)':
            return
        calc_search_input.value = value

    def calc_on_search_input_change(change):
        # For large text buffers, avoid running search on every input update.
        if len(state.get('file_content', '')) > CALC_HIGHLIGHT_MAX_CHARS:
            return
        calc_do_search()

    def _calc_request_chunk_start(requested):
        path_str = state.get('selected_file_path')
        size = int(state.get('selected_file_size') or 0)
        if not path_str or size <= CALC_TEXT_FULL_READ_BYTES:
            return False
        path = Path(path_str)
        if not path.exists():
            return False
        try:
            req = int(requested)
        except Exception:
            req = 0
        req = max(0, min(req, max(0, size - 1)))
        req = (req // CALC_TEXT_CHUNK_BYTES) * CALC_TEXT_CHUNK_BYTES
        current_start = int(state.get('file_chunk_start') or 0)
        if req == current_start:
            return False
        _calc_load_text_preview_chunk(path, size, req)
        return True

    def calc_prev_chunk(button):
        start = max(0, int(state.get('file_chunk_start') or 0) - CALC_TEXT_CHUNK_BYTES)
        _calc_request_chunk_start(start)

    def calc_next_chunk(button):
        size = int(state.get('selected_file_size') or 0)
        start = int(state.get('file_chunk_end') or 0)
        if start >= size:
            return
        _calc_request_chunk_start(start)

    def calc_on_chunk_request(button):
        _calc_request_chunk_start(calc_chunk_request_start.value)

    def calc_on_chunk_request_start_change(change):
        if change.get('name') != 'value':
            return
        _calc_request_chunk_start(change.get('new', 0))

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
            if not _calc_is_complete_mutation_csv(csv_path.name):
                calc_preselect_status.value = (
                    '<span style="color:#d32f2f;">Preselection only works with complete_mutation_space.csv.</span>'
                )
                _calc_preselect_show(True)
                return
            _calc_preselect_load(csv_path, mode='complete_preselection')
            _calc_preselect_show(True)
            _calc_preselect_render_current()
        elif change['new'] == 'Visualize':
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
            name = csv_path.name
            if _calc_is_complete_mutation_csv(name):
                mode = 'complete_visualize'
            elif _calc_is_preselected_mutation_csv(name):
                mode = 'preselected_visualize'
            elif _calc_is_rejected_mutation_csv(name):
                mode = 'rejected_visualize'
            else:
                calc_preselect_status.value = '<span style="color:#d32f2f;">No visualization mode for this file.</span>'
                _calc_preselect_show(True)
                return
            _calc_preselect_load(csv_path, mode=mode)
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
        state['search_truncated'] = False
        calc_search_result.value = ''
        calc_search_input.value = ''
        calc_update_nav_buttons()
        calc_delete_hide_confirm()
        calc_update_options_dropdown()

        next_suffix = full_path.suffix.lower()
        next_name_lower = full_path.name.lower()
        keep_previous_viewer_during_load = (
            next_name_lower == 'coord' or next_suffix in ['.xyz', '.cube', '.cub']
        )

        # Reset display (avoid flashing text area before viewer is ready)
        if not keep_previous_viewer_during_load:
            calc_mol_viewer.clear_output()
            calc_mol_container.layout.display = 'none'
        calc_content_area.layout.display = 'none'
        calc_content_label.layout.display = 'none'
        if not keep_previous_viewer_during_load:
            _set_view_toggle(False, True)
        calc_file_info.value = ''
        calc_path_display.value = ''
        state['file_content'] = ''
        state['selected_file_path'] = None
        state['selected_file_size'] = 0
        state['file_is_preview'] = False
        state['file_preview_note'] = ''
        state['file_chunk_start'] = 0
        state['file_chunk_end'] = 0
        _calc_hide_chunk_controls()
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
        state['selected_file_path'] = str(full_path)
        state['selected_file_size'] = int(size)

        # --- Turbomole coord file ---
        if name_lower == 'coord':
            _set_view_toggle(True, False)
            calc_mol_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
            calc_update_view()
            try:
                content = full_path.read_text()
                state['file_content'] = content
                state['file_is_preview'] = False
                state['file_preview_note'] = ''
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
                state['file_is_preview'] = False
                state['file_preview_note'] = ''
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
                    "<div style='height:100%; width:100%; border:1px solid #ddd; padding:6px;"
                    " background:#fafafa; box-sizing:border-box; display:flex;"
                    " align-items:center; justify-content:center; overflow:hidden;'>"
                    f"<img src='data:image/png;base64,{b64}' style='max-width:100%; max-height:100%;"
                    " width:auto; height:auto; object-fit:contain; display:block;' />"
                    "</div>"
                )
                calc_update_view()
            except Exception as e:
                calc_set_message(f'Error: {e}')
            return

        # --- Cube/Cub volumetric data ---
        if suffix in ['.cube', '.cub']:
            _set_view_toggle(True, False)
            calc_mol_container.layout.display = 'block'
            calc_content_area.layout.display = 'none'
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
                    def _cube_extra(view, raw):
                        view.addVolumetricData(raw, 'cube', {'isoval': 0.02, 'color': '#0026ff', 'opacity': 0.85})
                        view.addVolumetricData(raw, 'cube', {'isoval': -0.02, 'color': '#b00010', 'opacity': 0.85})
                    _render_3dmol(content, fmt='cube', extra_fn=_cube_extra)
                # Trigger dynamic resize after volumetric data is loaded
                _run_js("setTimeout(function(){ if(window.calcResizeMolViewer) window.calcResizeMolViewer(); }, 600);")
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
                state['file_is_preview'] = False
                state['file_preview_note'] = ''

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
            if size > CALC_TEXT_FULL_READ_BYTES:
                _calc_load_text_preview_chunk(full_path, size, 0)
                calc_copy_btn.disabled = False
                calc_copy_path_btn.disabled = False
                chunk_lines = state['file_content'].count('\n') + (1 if state['file_content'] else 0)
                calc_file_info.value = (
                    f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                    f' ({size_str}, chunked view, {chunk_lines} lines loaded)'
                )
                if suffix == '.inp':
                    state['selected_inp_path'] = full_path
                    state['selected_inp_base'] = full_path.stem
                    calc_recalc_btn.disabled = False
                calc_update_view()
                return

            content = full_path.read_text(errors='ignore')
            state['file_content'] = content
            state['file_is_preview'] = False
            state['file_preview_note'] = ''
            state['file_chunk_start'] = 0
            state['file_chunk_end'] = 0
            _calc_update_chunk_controls()
            calc_copy_btn.disabled = False
            calc_copy_path_btn.disabled = False
            line_count = content.count('\n') + (1 if content else 0)
            calc_file_info.value = (
                f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                f' ({size_str}, {line_count} lines)'
            )
            calc_render_content(scroll_to='top')

            if _calc_is_mutation_csv(name):
                _set_view_toggle(False, False)
                _calc_preselect_show(False)

            # input.txt/start.txt: show coordinates in molecule viewer
            if name in ('input.txt', 'start.txt'):
                # For SMILES-only input.txt, open the same side-by-side 2D/3D visual mode
                # used for mutation-space CSV browsing.
                if name == 'input.txt' and is_smiles(content):
                    _set_view_toggle(True, False)
                    _on_view_toggle()
                    return

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
    calc_search_input.observe(calc_on_search_input_change, names='value')
    calc_search_suggest.observe(calc_on_suggest, names='value')
    calc_chunk_prev_btn.on_click(calc_prev_chunk)
    calc_chunk_next_btn.on_click(calc_next_chunk)
    calc_chunk_request_btn.on_click(calc_on_chunk_request)
    calc_chunk_request_start.observe(calc_on_chunk_request_start_change, names='value')
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
    calc_preselect_prev.on_click(_calc_preselect_prev_entry)
    calc_preselect_next.on_click(_calc_preselect_next_entry)
    calc_preselect_move.on_click(_calc_preselect_move_current)
    calc_preselect_close.on_click(_calc_preselect_close_view)
    calc_preselect_new3d.on_click(_calc_preselect_new_3d_structure)
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
        ' border-radius:4px; display:block;'
        ' z-index:10; pointer-events:auto !important; position:relative; }'
        '.calc-splitter:hover { background:linear-gradient('
        'to right, #b0b0b0, #e0e0e0, #b0b0b0); }'
        '.calc-mol-viewer { overflow:hidden !important; padding:0 !important;'
        ' transition: height 0.15s ease-out; }'
        '.calc-mol-viewer .output_area, .calc-mol-viewer .output_subarea,'
        ' .calc-mol-viewer .output_wrapper, .calc-mol-viewer .jp-OutputArea-child,'
        ' .calc-mol-viewer .jp-OutputArea-output'
        ' { padding:0 !important; margin:0 !important; }'
        '.calc-preselect-viz { overflow:hidden !important; padding:0 !important; margin:0 !important; }'
        '.calc-preselect-viz-wrap { flex:1 1 220px !important; min-width:180px !important;'
        f' max-width:min({CALC_PRESELECT_VIZ_SIZE}px, 42vh) !important; width:100% !important; }}'
        '.calc-preselect-viz { width:100% !important; aspect-ratio:1 / 1 !important;'
        ' height:auto !important; max-width:100% !important; }'
        '.calc-preselect-viz .output_area, .calc-preselect-viz .output_subarea,'
        ' .calc-preselect-viz .output_wrapper, .calc-preselect-viz .jp-OutputArea-child,'
        ' .calc-preselect-viz .jp-OutputArea-output'
        ' { padding:0 !important; margin:0 !important; width:100% !important;'
        ' height:100% !important; overflow:hidden !important; }'
        '.calc-preselect-viz img { width:100% !important; height:100% !important; object-fit:contain !important; }'
        '.calc-preselect-viz [id^="3dmolviewer"], .calc-preselect-viz canvas {'
        ' width:100% !important; height:100% !important; }'
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
        calc_chunk_hidden_row,
        calc_override_status,
        calc_content_area,
        calc_edit_area,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='5px',
        overflow_x='hidden', overflow_y='hidden',
    ))

    calc_splitter = widgets.HTML(
        "<div class='calc-splitter' title='Drag to resize'></div>",
        layout=widgets.Layout(height='100%', width='10px', margin='0', overflow='visible'),
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

    # Enable draggable splitter + dynamic mol-viewer resize + $3Dmol patch
    # NOTE: must be a SINGLE _run_js call because run_js uses clear_output.
    _run_js(f"""
    setTimeout(function() {{
        /* --- Splitter drag logic --- */
        var left = document.querySelector('.calc-left');
        var right = document.querySelector('.calc-right');
        var splitter = document.querySelector('.calc-splitter');
        if (left && right && splitter && splitter.dataset.bound !== '1') {{
            splitter.dataset.bound = '1';
            var minW = {CALC_LEFT_MIN};
            var maxW = {CALC_LEFT_MAX};
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
                if (window.calcResizeMolViewer) setTimeout(window.calcResizeMolViewer, 50);
            }}
            splitter.addEventListener('mousedown', function(e) {{
                e.preventDefault();
                document.addEventListener('mousemove', onMove);
                document.addEventListener('mouseup', onUp);
            }});
        }}

        /* --- Monkey-patch $3Dmol.createViewer to store viewer instances --- */
        function patchCreateViewer() {{
            if (typeof $3Dmol !== 'undefined' && !$3Dmol._patched) {{
                var orig = $3Dmol.createViewer;
                $3Dmol.createViewer = function() {{
                    var v = orig.apply(this, arguments);
                    window._calcMolViewer = v;
                    setTimeout(function() {{
                        if (window.calcResizeMolViewer) window.calcResizeMolViewer();
                    }}, 300);
                    return v;
                }};
                $3Dmol._patched = true;
            }}
        }}
        patchCreateViewer();
        var patchInterval = setInterval(function() {{
            patchCreateViewer();
            if (typeof $3Dmol !== 'undefined' && $3Dmol._patched) clearInterval(patchInterval);
        }}, 500);

        /* --- Dynamic mol-viewer resize --- */
        window.calcResizeMolViewer = function() {{
            var lft = document.querySelector('.calc-left');
            var mv = document.querySelector('.calc-mol-viewer');
            if (!lft || !mv || mv.offsetParent === null) return;
            var container = mv.closest('.widget-vbox');
            if (!container || container.style.display === 'none') return;
            var leftRect = lft.getBoundingClientRect();
            var mvRect = mv.getBoundingClientRect();
            if (mvRect.top === 0 && mvRect.height === 0) return;
            var availH = leftRect.bottom - mvRect.top - 6;
            var availW = mv.parentElement.getBoundingClientRect().width - 4;
            var h = Math.floor(availH * {CALC_MOL_DYNAMIC_SCALE});
            var w = Math.floor(Math.min(h * 1.2, availW));
            if (h < 200) h = 200;
            if (w < 240) w = 240;
            mv.style.width = w + 'px';
            mv.style.height = h + 'px';
            var inner = mv.querySelector('[id^="3dmolviewer"], [id^="calc_trj"]');
            if (inner) {{ inner.style.width = w + 'px'; inner.style.height = h + 'px'; }}
            var v = window._calcMolViewer || window.calc_trj_viewer;
            if (v && typeof v.resize === 'function') {{ v.resize(); v.render(); }}
        }};
        window.addEventListener('resize', function() {{
            setTimeout(window.calcResizeMolViewer, 100);
        }});
        window.calcResizePreselect3D = function() {{
            var v = window._calcPreselect3dViewer;
            var elId = window._calcPreselect3dElId;
            if (!v || !elId) return;
            var el = document.getElementById(elId);
            var box = el ? el.closest('.calc-preselect-viz') : null;
            if (!el || !box || box.offsetParent === null) return;
            var rect = box.getBoundingClientRect();
            var side = Math.floor(Math.min(rect.width || 0, rect.height || 0));
            if (side >= 100) {{
                el.style.width = side + 'px';
                el.style.height = side + 'px';
            }}
            if (typeof v.resize === 'function') {{
                v.resize();
                v.render();
            }}
        }};
        window.addEventListener('resize', function() {{
            setTimeout(function() {{
                if (window.calcResizePreselect3D) window.calcResizePreselect3D();
            }}, 120);
        }});
        var tab = document.querySelector('.calc-tab');
        if (tab) {{
            new MutationObserver(function() {{
                setTimeout(window.calcResizeMolViewer, 200);
                setTimeout(function() {{
                    if (window.calcResizePreselect3D) window.calcResizePreselect3D();
                }}, 220);
            }}).observe(tab, {{attributes: true, subtree: true, attributeFilter: ['style']}});
        }}
    }}, 0);
    """)

    return tab_widget, {
        'calc_list_directory': calc_list_directory,
    }
