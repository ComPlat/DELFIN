"""Calculations Browser tab: file browser with viewer, search, recalc, and more."""

import base64
import html as _html
import json
import re
import shutil
import subprocess
import tempfile
import threading
from collections import Counter
from itertools import permutations, product
from pathlib import Path

import ipywidgets as widgets
import numpy as np
from IPython.display import HTML, clear_output, display

from .constants import CALC_SEARCH_OPTIONS
from .input_processing import (
    parse_inp_resources,
    parse_resource_settings,
    smiles_to_xyz,
    smiles_to_xyz_quick,
    contains_metal,
    is_smiles,
)
from .helpers import disable_spellcheck
from .molecule_viewer import coord_to_xyz, parse_xyz_frames, DEFAULT_3DMOL_STYLE_JS
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
    CALC_LEFT_DEFAULT = 360
    CALC_LEFT_MIN = 360
    CALC_LEFT_MAX = 520
    CALC_PRESELECT_VIZ_SIZE = 520
    CALC_RMSD_PANEL_HEIGHT = f'{CALC_MOL_SIZE}px'
    CALC_RMSD_COL_GAP = '18px'
    # Large-file guardrails
    CALC_TEXT_FULL_READ_BYTES = 2 * 1024 * 1024
    CALC_TEXT_CHUNK_BYTES = 8 * 1024 * 1024
    CALC_SEARCH_MAX_MATCHES = 2000
    CALC_HIGHLIGHT_MAX_CHARS = 400_000
    CALC_DOWNLOAD_MAX_BYTES = 25 * 1024 * 1024
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
        'select_mode': False,
        'file_is_preview': False,
        'file_preview_note': '',
        'file_chunk_start': 0,
        'file_chunk_end': 0,
        'chunk_dom_initialized': False,
        'xyz_frames': [],
        'xyz_current_frame': [0],
        'report_running': {},
        'traj_viewer_ready': False,
        'rmsd_available': False,
        'rmsd_mode_active': False,
        'rmsd_saved_display': {},
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
    calc_scope_id = f'calc-scope-{abs(id(state))}'
    calc_resize_mol_fn = f'calcResizeMolViewer_{abs(id(state))}'
    calc_resize_pre_fn = f'calcResizePreselect3D_{abs(id(state))}'

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
    calc_select_btn = widgets.Button(
        description='☑ Select',
        layout=widgets.Layout(width='76px', height='26px'),
    )

    # Detect whether we are inside the Archive tab (calc_dir == archiv_dir)
    _is_archive_tab = ctx.calc_dir.resolve() == ctx.archive_dir.resolve()

    # Move-to-Archive button (hidden when we are already inside the Archive)
    calc_move_archive_btn = widgets.Button(
        description='→ Archive', button_style='info',
        layout=widgets.Layout(
            width='90px', height='26px',
            display='none' if _is_archive_tab else 'inline-flex',
        ),
    )
    calc_move_archive_yes_btn = widgets.Button(
        description='Yes', button_style='warning',
        layout=widgets.Layout(width='60px', height='26px'),
    )
    calc_move_archive_no_btn = widgets.Button(
        description='No',
        layout=widgets.Layout(width='60px', height='26px'),
    )
    calc_move_archive_label = widgets.HTML('<b>Move to Archive?</b>')
    calc_move_archive_confirm = widgets.HBox(
        [calc_move_archive_label, calc_move_archive_yes_btn, calc_move_archive_no_btn],
        layout=widgets.Layout(display='none', gap='6px', align_items='center'),
    )
    calc_move_archive_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
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

    # Multi-select list (shown only in Select mode)
    calc_multi_select = widgets.SelectMultiple(
        options=[], rows=22,
        layout=widgets.Layout(
            width='100%', flex='1 1 0', min_height='0', margin='-4px 0 0 0',
            display='none',
        ),
    )
    calc_multi_select.add_class('delfin-multi-select')

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
        [
            widgets.HTML('<b>Frame:</b>'),
            calc_xyz_frame_input,
            calc_xyz_frame_total,
            calc_xyz_copy_btn,
        ],
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

    # RMSD controls (for single-frame XYZ files)
    calc_rmsd_ref_input = widgets.Textarea(
        value='',
        placeholder=(
            'Paste reference coordinates here (XYZ block or "Element x y z" lines).'
        ),
        layout=widgets.Layout(width='100%', height=CALC_RMSD_PANEL_HEIGHT),
    )
    calc_rmsd_ref_input.add_class('delfin-nospell')
    calc_rmsd_ref_input.add_class('calc-rmsd-ref-input')
    disable_spellcheck(ctx, class_name='calc-rmsd-ref-input')
    calc_rmsd_info_input = widgets.Text(
        value='',
        description='Info:',
        placeholder='Optional info for aligned comment line',
        layout=widgets.Layout(width='100%'),
        style={'description_width': '55px'},
    )
    calc_rmsd_info_input.add_class('delfin-nospell')
    calc_rmsd_info_input.add_class('calc-rmsd-info-input')
    disable_spellcheck(ctx, class_name='calc-rmsd-info-input')
    calc_rmsd_run_btn = widgets.Button(
        description='RMSD', button_style='success',
        layout=widgets.Layout(width='95px', height='32px'),
    )
    calc_rmsd_hide_btn = widgets.Button(
        description='Hide RMSD', button_style='warning',
        layout=widgets.Layout(width='120px', height='32px', display='none'),
    )
    calc_rmsd_hide_btn.add_class('calc-rmsd-trigger-btn')
    calc_rmsd_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    calc_rmsd_preview = widgets.Output(
        layout=widgets.Layout(
            width='100%',
            flex='1 1 0',
            min_width='0',
            height=CALC_RMSD_PANEL_HEIGHT,
            border='2px solid #1976d2',
            overflow='hidden',
            padding='0',
        ),
    )
    calc_rmsd_preview.add_class('calc-rmsd-preview')
    calc_rmsd_input_col = widgets.VBox(
        [
            calc_rmsd_info_input,
            widgets.HTML('<b>Reference Coordinates</b>'),
            calc_rmsd_ref_input,
            widgets.HBox(
                [calc_rmsd_run_btn, calc_rmsd_hide_btn],
                layout=widgets.Layout(gap='8px', align_items='center', margin='4px 0 0 0'),
            ),
            calc_rmsd_status,
        ],
        layout=widgets.Layout(
            width='100%',
            flex='1 1 0',
            min_width='0',
            margin='0',
        ),
    )
    calc_rmsd_controls = widgets.HBox(
        [calc_rmsd_input_col, calc_rmsd_preview],
        layout=widgets.Layout(
            display='none',
            width='100%',
            gap=CALC_RMSD_COL_GAP,
            overflow_x='hidden',
            flex_flow='row',
            align_items='stretch',
            margin='0 0 10px 0',
        ),
    )
    calc_rmsd_controls.add_class('calc-rmsd-controls')
    calc_rmsd_input_col.add_class('calc-rmsd-input-col')

    # Copy / path / report buttons
    calc_copy_path_btn = widgets.Button(
        description='PATH', button_style='',
        layout=widgets.Layout(width='70px', min_width='70px', height='26px'), disabled=True,
    )
    calc_copy_btn = widgets.Button(
        description='Copy', button_style='info',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px'), disabled=True,
    )
    calc_download_btn = widgets.Button(
        description='Download', button_style='',
        layout=widgets.Layout(width='130px', min_width='130px', height='26px'), disabled=True,
    )
    calc_report_btn = widgets.Button(
        description='Report', button_style='success',
        layout=widgets.Layout(width='80px', min_width='80px', height='26px'), disabled=True,
    )
    calc_download_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
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
        placeholder='Search in file...', continuous_update=True,
        layout=widgets.Layout(width='140px', min_width='140px', height='26px'),
    )
    calc_search_input.add_class('delfin-nospell')
    calc_search_input.add_class('calc-search-input')
    disable_spellcheck(ctx, class_name='calc-search-input')
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
        calc_rmsd_controls,
        calc_mol_viewer,
    ], layout=widgets.Layout(display='none', margin='0 0 10px 0', width='100%', align_items='stretch'))

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

    def _calc_label_to_name(label):
        if not label or label.startswith('('):
            return ''
        head, sep, tail = label.partition(' ')
        if sep:
            return tail.strip()
        return head.strip()

    def _calc_get_selected_path():
        selected = calc_file_list.value
        if not selected or selected.startswith('('):
            return None
        name = _calc_label_to_name(selected)
        if not name:
            return None
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
                        viewer.setStyle({{}}, {DEFAULT_3DMOL_STYLE_JS});
                        viewer.zoomTo();
                        viewer.center();
                        viewer.zoom(0.90);
                        viewer.render();
                        window._calcPreselect3dViewerByScope = window._calcPreselect3dViewerByScope || {{}};
                        window._calcPreselect3dElIdByScope = window._calcPreselect3dElIdByScope || {{}};
                        var scopeKey = {json.dumps(calc_scope_id)};
                        window._calcPreselect3dViewerByScope[scopeKey] = viewer;
                        window._calcPreselect3dElIdByScope[scopeKey] = "{viewer_id}";
                        if (window["{calc_resize_pre_fn}"]) {{
                            setTimeout(window["{calc_resize_pre_fn}"], 120);
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

    def _calc_clear_main_viewer_state(reset_view_state=False):
        calc_mol_viewer.clear_output()
        scope_key_json = json.dumps(calc_scope_id)
        clear_views_flag = 'true' if reset_view_state else 'false'
        _run_js(
            f"""
            (function() {{
                var scopeKey = {scope_key_json};
                var scopeRoot = document.querySelector('.{calc_scope_id}');
                var wrappers = scopeRoot
                    ? scopeRoot.querySelectorAll('.calc-mol-stage-wrapper')
                    : document.querySelectorAll('.calc-mol-stage-wrapper');
                wrappers.forEach(function(w) {{ w.remove(); }});
                if (window._calcMolViewerByScope) {{
                    delete window._calcMolViewerByScope[scopeKey];
                }}
                if (window._calcTrajViewerByScope) {{
                    delete window._calcTrajViewerByScope[scopeKey];
                }}
                if ({clear_views_flag} && window._calcMolViewStateByScope) {{
                    var prefix = scopeKey + ':';
                    Object.keys(window._calcMolViewStateByScope).forEach(function(k) {{
                        if (k.indexOf(prefix) === 0) {{
                            delete window._calcMolViewStateByScope[k];
                        }}
                    }});
                }}
            }})();
            """
        )

    _mol3d_counter = [0]

    def _render_3dmol(data, fmt='xyz', extra_fn=None):
        """Render a 3D molecule via JS with correct initial sizing."""
        import json
        _mol3d_counter[0] += 1
        viewer_id = f"mol3d_{_mol3d_counter[0]}"
        wrapper_id = f"calc_mol_wrap_{_mol3d_counter[0]}"
        data_json = json.dumps(data)
        view_scope_json = json.dumps(f"{calc_scope_id}:{state.get('current_path') or '/'}")
        scope_id_json = json.dumps(calc_scope_id)

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
                var scopeRoot = el ? el.closest('.{calc_scope_id}') : null;
                if (!scopeRoot) scopeRoot = document.querySelector('.{calc_scope_id}');
                if (!el || typeof $3Dmol === "undefined"
                    || !mv || mv.offsetParent === null) {{
                    tries += 1;
                    if (tries < 80) setTimeout(initViewer, 50);
                    return;
                }}
                /* Compute correct size from layout before creating viewer */
                var lft = scopeRoot ? scopeRoot.querySelector('.calc-left') : null;
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
                window._calcMolViewerByScope = window._calcMolViewerByScope || {{}};
                window._calcMolViewScopeKeyByScope = window._calcMolViewScopeKeyByScope || {{}};
                window._calcTrajViewerByScope = window._calcTrajViewerByScope || {{}};
                var scopeKey = {scope_id_json};
                var viewScope = {view_scope_json};
                var previousViewer =
                    window._calcMolViewerByScope[scopeKey]
                    || window._calcTrajViewerByScope[scopeKey]
                    || null;
                var previousScope =
                    window._calcMolViewScopeKeyByScope[scopeKey] || viewScope;
                if (previousViewer && typeof previousViewer.getView === 'function') {{
                    try {{
                        window._calcMolViewStateByScope[previousScope] = previousViewer.getView();
                    }} catch (_e) {{}}
                }}
                var savedView = window._calcMolViewStateByScope[viewScope] || null;
                var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
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
                        viewer.zoom({CALC_MOL_ZOOM});
                    }}
                }} else {{
                    viewer.zoomTo();
                    viewer.center();
                    viewer.zoom({CALC_MOL_ZOOM});
                }}
                viewer.render();
                window._calcMolViewerByScope[scopeKey] = viewer;
                window._calcMolViewScopeKeyByScope[scopeKey] = viewScope;
                var wrappers = scopeRoot
                    ? scopeRoot.querySelectorAll('.calc-mol-stage-wrapper')
                    : document.querySelectorAll('.calc-mol-stage-wrapper');
                wrappers.forEach(function(w) {{
                    if (w.id !== "{wrapper_id}") w.remove();
                }});
                if (window["{calc_resize_mol_fn}"]) {{
                    setTimeout(window["{calc_resize_mol_fn}"], 200);
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
        name = _calc_label_to_name(selected)
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

    def _calc_selected_item_path():
        selected = calc_file_list.value
        if not selected or selected.startswith('('):
            return None
        name = _calc_label_to_name(selected)
        if not name:
            return None
        full_path = (
            (_calc_dir() / state['current_path'] / name)
            if state['current_path']
            else (_calc_dir() / name)
        )
        return full_path if full_path.exists() else None

    def _calc_download_target_path():
        selected_path = _calc_selected_item_path()
        if selected_path is not None:
            return selected_path
        if state['current_path']:
            current_dir = _calc_dir() / state['current_path']
            if current_dir.exists():
                return current_dir
        return None

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

    def _calc_update_chunk_dom(content, size_bytes, chunk_start, chunk_end):
        total_size = max(1, int(size_bytes))
        c_start = max(0, int(chunk_start))
        c_end = max(c_start, int(chunk_end))
        virtual_h = max(12000, min(180000, int(total_size / 96)))
        top_px = int((c_start / total_size) * virtual_h)
        bottom_px = int((max(0, total_size - c_end) / total_size) * virtual_h)
        fallback_ratio = c_start / max(1, total_size - 1)
        js = r"""
        setTimeout(function() {
            const box = document.getElementById('calc-content-box');
            const topSpacer = document.getElementById('calc-chunk-top-spacer');
            const bottomSpacer = document.getElementById('calc-chunk-bottom-spacer');
            const txt = document.getElementById('calc-content-text');
            if (!box || !topSpacer || !bottomSpacer || !txt) return;
            const req = window.__calcChunkRequestedRatio;
            const fallback = __FALLBACK_RATIO__;
            const targetRatio = (typeof req === 'number' && isFinite(req))
                ? Math.max(0, Math.min(1, req))
                : Math.max(0, Math.min(1, fallback));
            txt.style.opacity = '0';
            topSpacer.style.height = '__TOP_PX__px';
            bottomSpacer.style.height = '__BOTTOM_PX__px';
            txt.textContent = __TEXT_JSON__;
            const maxScroll = Math.max(0, box.scrollHeight - box.clientHeight);
            window.__calcChunkProgrammaticScroll = true;
            box.scrollTop = Math.floor(targetRatio * maxScroll);
            window.__calcChunkRequestedRatio = null;
            setTimeout(function() {
                window.__calcChunkProgrammaticScroll = false;
                window.__calcChunkBusy = false;
                txt.style.opacity = '1';
                const fn = window.__calcChunkProcessPendingRatio;
                if (typeof fn === 'function') fn(false);
            }, 20);
        }, 0);
        """
        js = (
            js.replace('__TOP_PX__', str(top_px))
              .replace('__BOTTOM_PX__', str(bottom_px))
              .replace('__FALLBACK_RATIO__', f'{fallback_ratio:.12f}')
              .replace('__TEXT_JSON__', json.dumps(content))
        )
        _run_js(js)

    def _calc_load_text_preview_chunk(path, size_bytes, start_byte, rerun_search=False, scroll_to=None):
        content, chunk_start, chunk_end = _calc_read_text_chunk(path, size_bytes, start_byte)
        state['file_content'] = content
        state['file_is_preview'] = True
        state['file_chunk_start'] = chunk_start
        state['file_chunk_end'] = chunk_end
        state['file_preview_note'] = ''
        if state.get('chunk_dom_initialized'):
            _calc_update_chunk_dom(content, size_bytes, chunk_start, chunk_end)
            if scroll_to:
                calc_scroll_to(scroll_to)
        else:
            calc_render_content(scroll_to=scroll_to)
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

    def _calc_set_requested_ratio_for_offset(byte_offset):
        size = int(state.get('selected_file_size') or 0)
        if size <= 0:
            return
        try:
            off = int(byte_offset)
        except Exception:
            off = 0
        off = max(0, min(off, max(0, size - 1)))
        ratio = off / max(1, size - 1)
        _run_js(f"window.__calcChunkRequestedRatio = {ratio:.12f};")

    def _calc_find_global_matches(path, query, max_matches):
        q_bytes = query.encode('utf-8', errors='ignore')
        if not q_bytes:
            return []
        # Fast path: use ripgrep for large-file fixed-string matching.
        # Falls back to Python scanning if rg is unavailable.
        if '\n' not in query and '\r' not in query:
            try:
                cmd = [
                    'rg',
                    '-a',
                    '-i',
                    '-F',
                    '--byte-offset',
                    '--only-matching',
                    '--no-filename',
                    '--no-line-number',
                    '--color',
                    'never',
                    '--max-count',
                    str(int(max_matches)),
                    query,
                    str(path),
                ]
                proc = subprocess.run(
                    cmd, capture_output=True, text=True, check=False,
                )
                # 0: matches found, 1: no matches. Anything else means fallback.
                if proc.returncode in (0, 1):
                    spans = []
                    q_len = len(q_bytes)
                    for line in proc.stdout.splitlines():
                        if not line:
                            continue
                        off_txt, _, _rest = line.partition(':')
                        try:
                            abs_start = int(off_txt)
                        except Exception:
                            continue
                        spans.append((abs_start, abs_start + q_len))
                        if len(spans) >= max_matches:
                            break
                    return spans
            except Exception:
                pass
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

    def calc_update_download_btn():
        target = _calc_download_target_path()
        if target is None:
            calc_download_btn.disabled = True
            calc_download_btn.description = 'Download'
            return
        calc_download_btn.disabled = False
        calc_download_btn.description = 'Download ZIP' if target.is_dir() else 'Download'

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
        text_opacity = '1'
        if chunk_mode:
            total_size = max(1, int(state.get('selected_file_size') or 1))
            chunk_start = max(0, int(state.get('file_chunk_start') or 0))
            chunk_end = max(chunk_start, int(state.get('file_chunk_end') or chunk_start))
            virtual_h = max(12000, min(180000, int(total_size / 96)))
            top_px = int((chunk_start / total_size) * virtual_h)
            bottom_px = int((max(0, total_size - chunk_end) / total_size) * virtual_h)
            top_spacer_html = f"<div id='calc-chunk-top-spacer' style='height:{top_px}px;'></div>"
            bottom_spacer_html = f"<div id='calc-chunk-bottom-spacer' style='height:{bottom_px}px;'></div>"
            text_opacity = '0'
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
            f" word-break:break-word; font-family:monospace; font-size:12px; line-height:1.3; opacity:{text_opacity};'>"
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
                    if (window.__calcChunkBusy) {
                        window.__calcChunkPendingRatio = clamped;
                        return;
                    }
                    const tailStart = Math.max(0, fileSize - chunkBytes);
                    let desired;
                    if (clamped >= 0.9995) {
                        desired = tailStart;
                    } else {
                        desired = Math.floor((clamped * fileSize) / chunkBytes) * chunkBytes;
                        if (desired < 0) desired = 0;
                        if (desired > tailStart) desired = tailStart;
                    }
                    if (desired === currentStart) return;
                    window.__calcChunkBusy = true;
                    window.__calcChunkPendingRatio = null;
                    window.__calcChunkRequestedRatio = clamped;
                    emit(startInput, desired);
                    emit(ratioInput, clamped);
                    setTimeout(function() {
                        if (!window.__calcChunkBusy) return;
                        window.__calcChunkBusy = false;
                        processPendingRatio();
                    }, 2000);
                }

                function processPendingRatio(forceRun) {
                    if (window.__calcChunkBusy) return;
                    if (!forceRun && window.__calcChunkScrollbarDragActive) return;
                    const pending = window.__calcChunkPendingRatio;
                    if (typeof pending !== 'number' || !isFinite(pending)) return;
                    const now = Date.now();
                    const fastUntil = (window.__calcChunkFastDragUntil || 0);
                    if (!forceRun && fastUntil > now) {
                        const waitFast = Math.max(20, (fastUntil - now) + 5);
                        if (window.__calcChunkResumeTimer) clearTimeout(window.__calcChunkResumeTimer);
                        window.__calcChunkResumeTimer = setTimeout(function() {
                            window.__calcChunkResumeTimer = null;
                            processPendingRatio(false);
                        }, waitFast);
                        return;
                    }
                    const ignoreUntil = (window.__calcChunkIgnoreScrollUntil || 0);
                    if (!forceRun && ignoreUntil > now) {
                        const waitMs = Math.max(20, (ignoreUntil - now) + 5);
                        if (window.__calcChunkResumeTimer) clearTimeout(window.__calcChunkResumeTimer);
                        window.__calcChunkResumeTimer = setTimeout(function() {
                            window.__calcChunkResumeTimer = null;
                            processPendingRatio(false);
                        }, waitMs);
                        return;
                    }
                    window.__calcChunkPendingRatio = null;
                    requestChunkForRatio(pending);
                }
                window.__calcChunkProcessPendingRatio = processPendingRatio;

                if (!window.__calcChunkReleaseBound) {
                    window.__calcChunkReleaseBound = true;
                    const onRelease = function() {
                        if (!window.__calcChunkScrollbarDragActive) return;
                        window.__calcChunkScrollbarDragActive = false;
                        const fn = window.__calcChunkProcessPendingRatio;
                        if (typeof fn === 'function') fn(true);
                    };
                    window.addEventListener('pointerup', onRelease, {passive: true});
                    window.addEventListener('mouseup', onRelease, {passive: true});
                    window.addEventListener('pointercancel', onRelease, {passive: true});
                    window.addEventListener('blur', onRelease, {passive: true});
                }

                if (!box.dataset.chunkBound) {
                    box.dataset.chunkBound = '1';
                    let timer = null;
                    box.addEventListener('wheel', function() {
                        window.__calcChunkLastWheelTs = Date.now();
                    }, {passive: true});
                    box.addEventListener('pointerdown', function(ev) {
                        if (window.__calcChunkProgrammaticScroll) return;
                        if (!ev) return;
                        const rect = box.getBoundingClientRect();
                        if (ev.clientX >= (rect.right - 28)) {
                            window.__calcChunkScrollbarDragActive = true;
                        }
                    }, {passive: true});
                    box.addEventListener('scroll', function(ev) {
                        if (window.__calcChunkProgrammaticScroll) {
                            return;
                        }
                        const nowTs = Date.now();
                        const maxScroll = Math.max(1, box.scrollHeight - box.clientHeight);
                        const ratio = box.scrollTop / maxScroll;
                        const prevRatio = window.__calcChunkPrevScrollRatio;
                        const delta = (typeof prevRatio === 'number' && isFinite(prevRatio))
                            ? Math.abs(ratio - prevRatio)
                            : 0;
                        const wheelRecent = ((nowTs - (window.__calcChunkLastWheelTs || 0)) <= 90);
                        if (!wheelRecent && delta >= 0.02) {
                            window.__calcChunkScrollbarDragActive = true;
                        }
                        if (typeof prevRatio === 'number' && isFinite(prevRatio)) {
                            if (delta >= 0.08) {
                                window.__calcChunkFastDragUntil = nowTs + 180;
                            }
                        }
                        window.__calcChunkPrevScrollRatio = ratio;
                        const buttonsDown = !!(ev && typeof ev.buttons === 'number' && ev.buttons > 0);
                        if (buttonsDown) {
                            window.__calcChunkScrollbarDragActive = true;
                        }
                        if (window.__calcChunkScrollbarDragActive) {
                            if (window.__calcChunkDragIdleTimer) {
                                clearTimeout(window.__calcChunkDragIdleTimer);
                            }
                            window.__calcChunkDragIdleTimer = setTimeout(function() {
                                window.__calcChunkDragIdleTimer = null;
                                if (!window.__calcChunkScrollbarDragActive) return;
                                window.__calcChunkScrollbarDragActive = false;
                                processPendingRatio(true);
                            }, 130);
                        }
                        window.__calcChunkPendingRatio = ratio;
                        if (window.__calcChunkScrollbarDragActive && window.__calcChunkBusy) {
                            return;
                        }
                        if (timer) clearTimeout(timer);
                        timer = setTimeout(function() {
                            processPendingRatio(false);
                        }, 30);
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
                    window.__calcChunkBusy = false;
                    if (txt) txt.style.opacity = '1';
                    processPendingRatio(false);
                }, 40);
            }, 0);
            """.replace('__FILE_SIZE__', str(file_size_js))
              .replace('__CHUNK_BYTES__', str(CALC_TEXT_CHUNK_BYTES))
              .replace('__CURRENT_START__', str(chunk_start_js)))
            state['chunk_dom_initialized'] = True
        else:
            state['chunk_dom_initialized'] = False
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
            if (__INDEX__ >= 0) {
                setTimeout(function() {
                    const mark = document.getElementById('calc-current-match');
                    if (!mark) return;
                    const b = document.getElementById('calc-content-box');
                    if (!b) return;
                    const br = b.getBoundingClientRect();
                    const mr = mark.getBoundingClientRect();
                    b.scrollTop += (mr.top - br.top) - (b.clientHeight / 2);
                }, 0);
            }
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
            xyz_string, num_atoms, _method, error = smiles_to_xyz_quick(smiles_line)
            if not error and xyz_string:
                return f"{num_atoms}\n{title}\n{xyz_string}"
        atom_count = len(lines)
        return f"{atom_count}\n{title}\n" + '\n'.join(lines)

    def _calc_reset_rmsd_ui(clear_input=False):
        calc_rmsd_controls.layout.display = 'none'
        _calc_set_rmsd_mode(False)
        calc_rmsd_status.value = ''
        calc_rmsd_hide_btn.layout.display = 'none'
        calc_rmsd_preview.clear_output()
        if clear_input:
            calc_rmsd_ref_input.value = ''
            calc_rmsd_info_input.value = ''

    def _calc_set_rmsd_mode(active):
        if active:
            if state.get('rmsd_mode_active'):
                return
            saved = {}
            widgets_to_hide = [
                ('calc_mol_label', calc_mol_label),
                ('calc_xyz_frame_label', calc_xyz_frame_label),
                ('calc_xyz_controls', calc_xyz_controls),
                ('calc_coord_controls', calc_coord_controls),
                ('calc_mol_viewer', calc_mol_viewer),
            ]
            for key, widget in widgets_to_hide:
                saved[key] = widget.layout.display if widget.layout.display is not None else ''
                widget.layout.display = 'none'
            state['rmsd_saved_display'] = saved
            state['rmsd_mode_active'] = True
            return

        if not state.get('rmsd_mode_active'):
            return
        saved = state.get('rmsd_saved_display') or {}
        widgets_to_restore = [
            ('calc_mol_label', calc_mol_label),
            ('calc_xyz_frame_label', calc_xyz_frame_label),
            ('calc_xyz_controls', calc_xyz_controls),
            ('calc_coord_controls', calc_coord_controls),
            ('calc_mol_viewer', calc_mol_viewer),
        ]
        for key, widget in widgets_to_restore:
            widget.layout.display = saved.get(key, '')
        state['rmsd_saved_display'] = {}
        state['rmsd_mode_active'] = False

    def _calc_set_rmsd_available(enabled):
        state['rmsd_available'] = bool(enabled)
        if not enabled:
            _calc_reset_rmsd_ui(clear_input=False)
        calc_update_options_dropdown()

    def _calc_parse_xyz_line(line):
        parts = line.split()
        if len(parts) < 4:
            return None
        index_patterns = [(0, 1, 2, 3)]
        if len(parts) >= 5:
            index_patterns.append((1, 2, 3, 4))
        for sym_idx, x_idx, y_idx, z_idx in index_patterns:
            if z_idx >= len(parts):
                continue
            symbol_raw = re.sub(r'[^A-Za-z]', '', parts[sym_idx] or '')
            if not symbol_raw:
                continue
            symbol = symbol_raw[0].upper() + symbol_raw[1:].lower()
            try:
                x = float(parts[x_idx].replace('D', 'E').replace('d', 'e'))
                y = float(parts[y_idx].replace('D', 'E').replace('d', 'e'))
                z = float(parts[z_idx].replace('D', 'E').replace('d', 'e'))
            except (TypeError, ValueError):
                continue
            return symbol, np.array([x, y, z], dtype=float)
        return None

    def _calc_parse_xyz_lines(lines, expected_atoms=None, allow_skip=False):
        symbols = []
        coords = []
        for raw_line in lines:
            stripped = raw_line.strip()
            if not stripped:
                continue
            parsed = _calc_parse_xyz_line(stripped)
            if parsed is None:
                if allow_skip:
                    continue
                raise ValueError(f'Could not parse XYZ line: {stripped[:120]}')
            symbol, xyz = parsed
            symbols.append(symbol)
            coords.append(xyz)
        if not symbols:
            raise ValueError('No coordinate lines found.')
        if expected_atoms is not None and len(symbols) != int(expected_atoms):
            raise ValueError(
                f'Atom count mismatch: expected {expected_atoms}, parsed {len(symbols)}.'
            )
        return symbols, np.asarray(coords, dtype=float)

    def _calc_build_xyz_from_symbols_coords(symbols, coords, comment='Aligned reference'):
        body = []
        for symbol, (x, y, z) in zip(symbols, coords):
            body.append(f'{symbol:<2}  {x: .8f}  {y: .8f}  {z: .8f}')
        return f'{len(symbols)}\n{comment}\n' + '\n'.join(body) + '\n'

    def _calc_parse_reference_xyz_input(raw_text):
        text = (raw_text or '').strip()
        if not text:
            raise ValueError('No reference coordinates pasted.')
        frames = parse_xyz_frames(text)
        if frames:
            comment, xyz_block, n_atoms = frames[0]
            symbols, coords = _calc_parse_xyz_lines(
                xyz_block.splitlines(), expected_atoms=n_atoms, allow_skip=False
            )
            return symbols, coords, (comment or 'Reference')
        symbols, coords = _calc_parse_xyz_lines(text.splitlines(), allow_skip=True)
        return symbols, coords, 'Reference'

    def _calc_get_current_xyz_model():
        frames = state['xyz_frames']
        if not frames:
            raise ValueError('No XYZ frame is currently loaded.')
        idx = state['xyz_current_frame'][0]
        if idx < 0 or idx >= len(frames):
            raise ValueError('Invalid XYZ frame index.')
        comment, xyz_block, n_atoms = frames[idx]
        symbols, coords = _calc_parse_xyz_lines(
            xyz_block.splitlines(), expected_atoms=n_atoms, allow_skip=False
        )
        xyz_content = f'{n_atoms}\n{comment}\n{xyz_block}\n'
        return symbols, coords, comment, xyz_content

    def _calc_kabsch_align(reference_coords, target_coords, mapping=None, return_rotation=False):
        ref = np.asarray(reference_coords, dtype=float)
        target = np.asarray(target_coords, dtype=float)
        if ref.shape != target.shape:
            raise ValueError(
                f'Coordinate shape mismatch: {ref.shape} vs {target.shape}.'
            )

        if mapping is not None:
            map_idx = np.asarray(mapping, dtype=int)
            if map_idx.ndim != 1 or map_idx.size != ref.shape[0]:
                raise ValueError('Invalid mapping length for Kabsch alignment.')
            if np.any(map_idx < 0) or np.any(map_idx >= target.shape[0]):
                raise ValueError('Mapping indices out of bounds.')
            if np.unique(map_idx).size != map_idx.size:
                raise ValueError('Mapping must be one-to-one.')
            target = target[map_idx]

        ref_centroid = ref.mean(axis=0)
        target_centroid = target.mean(axis=0)
        ref_centered = ref - ref_centroid
        target_centered = target - target_centroid
        covariance = ref_centered.T @ target_centered
        u, _s, vt = np.linalg.svd(covariance)
        rotation = vt.T @ u.T
        if np.linalg.det(rotation) < 0:
            vt[-1, :] *= -1.0
            rotation = vt.T @ u.T
        aligned = ref_centered @ rotation + target_centroid
        diff = aligned - target
        rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
        if return_rotation:
            return aligned, rmsd, rotation
        return aligned, rmsd

    def _calc_sq_distance_matrix(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        diff = a[:, None, :] - b[None, :, :]
        return np.einsum('ijk,ijk->ij', diff, diff)

    def _calc_element_assignment_for_rotation(
        ref_symbols_lc,
        target_symbols_lc,
        ref_rot_centered,
        target_centered,
    ):
        try:
            from scipy.optimize import linear_sum_assignment
        except Exception as exc:
            raise RuntimeError(
                'SciPy is required for permutation-based RMSD alignment.'
            ) from exc

        n_atoms = len(ref_symbols_lc)
        mapping = np.full(n_atoms, -1, dtype=int)

        ref_groups = {}
        target_groups = {}
        for idx, symbol in enumerate(ref_symbols_lc):
            ref_groups.setdefault(symbol, []).append(idx)
        for idx, symbol in enumerate(target_symbols_lc):
            target_groups.setdefault(symbol, []).append(idx)

        if set(ref_groups.keys()) != set(target_groups.keys()):
            return None

        for symbol, ref_idx in ref_groups.items():
            target_idx = target_groups.get(symbol, [])
            if len(ref_idx) != len(target_idx):
                return None
            if len(ref_idx) == 1:
                mapping[ref_idx[0]] = target_idx[0]
                continue

            ref_block = ref_rot_centered[np.asarray(ref_idx, dtype=int)]
            target_block = target_centered[np.asarray(target_idx, dtype=int)]
            costs = _calc_sq_distance_matrix(ref_block, target_block)
            row_ind, col_ind = linear_sum_assignment(costs)
            for row_pos, col_pos in zip(row_ind, col_ind):
                mapping[int(ref_idx[row_pos])] = int(target_idx[col_pos])

        if np.any(mapping < 0):
            return None
        return mapping

    def _calc_generate_proper_axis_rotations():
        mats = []
        for perm in permutations((0, 1, 2)):
            for signs in product((-1.0, 1.0), repeat=3):
                mat = np.zeros((3, 3), dtype=float)
                for new_axis, old_axis in enumerate(perm):
                    mat[old_axis, new_axis] = signs[new_axis]
                if np.linalg.det(mat) > 0.0:
                    mats.append(mat)
        return mats

    def _calc_random_rotation_matrices(n_samples=64, seed=20260220):
        rng = np.random.default_rng(seed)
        mats = []
        for _ in range(max(0, int(n_samples))):
            q = rng.normal(size=4)
            norm = float(np.linalg.norm(q))
            if norm < 1e-12:
                continue
            q = q / norm
            w, x, y, z = q
            mat = np.array([
                [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
                [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
                [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
            ], dtype=float)
            if np.linalg.det(mat) > 0.0:
                mats.append(mat)
        return mats

    def _calc_rdkit_best_alignment_from_xyz(
        ref_symbols,
        ref_coords,
        target_symbols,
        target_coords,
        max_matches=20000,
    ):
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds, rdMolAlign
        except Exception:
            return None

        ref_xyz = _calc_build_xyz_from_symbols_coords(ref_symbols, ref_coords, comment='Reference')
        target_xyz = _calc_build_xyz_from_symbols_coords(target_symbols, target_coords, comment='Target')
        ref_mol = Chem.MolFromXYZBlock(ref_xyz)
        target_mol = Chem.MolFromXYZBlock(target_xyz)
        if ref_mol is None or target_mol is None:
            return None

        try:
            rdDetermineBonds.DetermineConnectivity(ref_mol)
            rdDetermineBonds.DetermineConnectivity(target_mol)
        except Exception:
            return None

        try:
            rmsd, transform, atom_map = rdMolAlign.GetBestAlignmentTransform(
                ref_mol,
                target_mol,
                maxMatches=int(max_matches),
                reflect=False,
            )
        except Exception:
            return None

        if not atom_map:
            return None
        n_atoms = len(ref_symbols)
        mapping = np.full(n_atoms, -1, dtype=int)
        for probe_idx, ref_idx in atom_map:
            pi = int(probe_idx)
            ri = int(ref_idx)
            if pi < 0 or pi >= n_atoms or ri < 0 or ri >= n_atoms:
                return None
            mapping[pi] = ri
        if np.any(mapping < 0) or np.unique(mapping).size != n_atoms:
            return None

        ref_seq = [s.lower() for s in ref_symbols]
        target_seq = [s.lower() for s in target_symbols]
        if not all(ref_seq[i] == target_seq[mapping[i]] for i in range(n_atoms)):
            return None

        transform = np.asarray(transform, dtype=float)
        if transform.shape != (4, 4):
            return None
        ref_arr = np.asarray(ref_coords, dtype=float)
        aligned = ref_arr @ transform[:3, :3].T + transform[:3, 3]
        return aligned, float(rmsd), mapping

    def _calc_topology_mappings_from_xyz(
        ref_symbols,
        ref_coords,
        target_symbols,
        target_coords,
        max_matches=4096,
    ):
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
        except Exception:
            return []

        ref_xyz = _calc_build_xyz_from_symbols_coords(ref_symbols, ref_coords, comment='Reference')
        target_xyz = _calc_build_xyz_from_symbols_coords(target_symbols, target_coords, comment='Target')
        ref_mol = Chem.MolFromXYZBlock(ref_xyz)
        target_mol = Chem.MolFromXYZBlock(target_xyz)
        if ref_mol is None or target_mol is None:
            return []

        try:
            rdDetermineBonds.DetermineConnectivity(ref_mol)
            rdDetermineBonds.DetermineConnectivity(target_mol)
        except Exception:
            return []

        if ref_mol.GetNumAtoms() != target_mol.GetNumAtoms():
            return []

        try:
            matches = target_mol.GetSubstructMatches(
                ref_mol,
                uniquify=True,
                useChirality=False,
                maxMatches=max_matches,
            )
        except TypeError:
            matches = target_mol.GetSubstructMatches(
                ref_mol,
                uniquify=True,
                useChirality=False,
            )

        ref_seq = [s.lower() for s in ref_symbols]
        target_seq = [s.lower() for s in target_symbols]
        uniq = set()
        mappings = []
        for match in matches:
            if len(match) != len(ref_seq):
                continue
            key = tuple(int(v) for v in match)
            if key in uniq:
                continue
            if not all(ref_seq[i] == target_seq[j] for i, j in enumerate(key)):
                continue
            uniq.add(key)
            mappings.append(np.asarray(key, dtype=int))
            if len(mappings) >= max_matches:
                break
        return mappings

    def _calc_invert_mapping(mapping, n_atoms):
        map_idx = np.asarray(mapping, dtype=int)
        if map_idx.ndim != 1 or map_idx.size != int(n_atoms):
            return None
        inv = np.full(int(n_atoms), -1, dtype=int)
        for src_idx, dst_idx in enumerate(map_idx.tolist()):
            dst = int(dst_idx)
            if dst < 0 or dst >= int(n_atoms):
                return None
            if inv[dst] != -1:
                return None
            inv[dst] = int(src_idx)
        if np.any(inv < 0):
            return None
        return inv

    def _calc_local_swap_optimize_mapping(
        ref_symbols_lc,
        ref_coords,
        target_coords,
        mapping,
        max_passes=2,
        max_group_size=12,
        max_total_swaps=2200,
    ):
        ref = np.asarray(ref_coords, dtype=float)
        target = np.asarray(target_coords, dtype=float)
        best_mapping = np.asarray(mapping, dtype=int).copy()
        best_aligned, best_rmsd = _calc_kabsch_align(ref, target, mapping=best_mapping)

        groups = {}
        for idx, symbol in enumerate(ref_symbols_lc):
            groups.setdefault(symbol, []).append(idx)

        swap_checks = 0
        for _ in range(max(0, int(max_passes))):
            improved = False
            for idx_list in groups.values():
                if len(idx_list) < 2 or len(idx_list) > int(max_group_size):
                    continue
                for i_pos in range(len(idx_list) - 1):
                    i = int(idx_list[i_pos])
                    for j_pos in range(i_pos + 1, len(idx_list)):
                        j = int(idx_list[j_pos])
                        trial_mapping = best_mapping.copy()
                        trial_mapping[i], trial_mapping[j] = trial_mapping[j], trial_mapping[i]
                        trial_aligned, trial_rmsd = _calc_kabsch_align(
                            ref, target, mapping=trial_mapping
                        )
                        swap_checks += 1
                        if trial_rmsd < best_rmsd - 1e-12:
                            best_mapping = trial_mapping
                            best_aligned = trial_aligned
                            best_rmsd = float(trial_rmsd)
                            improved = True
                        if swap_checks >= int(max_total_swaps):
                            return best_aligned, best_rmsd, best_mapping
            if not improved:
                break
        return best_aligned, best_rmsd, best_mapping

    def _calc_refine_alignment_with_mapping(
        ref_symbols_lc,
        target_symbols_lc,
        ref_coords,
        target_coords,
        initial_mapping,
        max_iter=14,
    ):
        ref = np.asarray(ref_coords, dtype=float)
        target = np.asarray(target_coords, dtype=float)
        ref_centered = ref - ref.mean(axis=0)
        target_centered = target - target.mean(axis=0)
        mapping = np.asarray(initial_mapping, dtype=int)

        visited = set()
        best = None
        for _ in range(max_iter):
            key = tuple(int(v) for v in mapping.tolist())
            if key in visited:
                break
            visited.add(key)

            aligned, rmsd, rotation = _calc_kabsch_align(
                ref, target, mapping=mapping, return_rotation=True
            )
            if best is None or rmsd < best[0] - 1e-12:
                best = (rmsd, aligned.copy(), mapping.copy())

            ref_rot_centered = ref_centered @ rotation
            new_mapping = _calc_element_assignment_for_rotation(
                ref_symbols_lc, target_symbols_lc, ref_rot_centered, target_centered
            )
            if new_mapping is None or np.array_equal(new_mapping, mapping):
                break
            mapping = new_mapping

        if best is None:
            aligned, rmsd = _calc_kabsch_align(ref, target, mapping=mapping)
            best_aligned, best_rmsd, best_mapping = aligned, float(rmsd), mapping.copy()
        else:
            best_aligned, best_rmsd, best_mapping = best[1], float(best[0]), best[2].copy()

        # Local remapping refinement: swap assignments of same-element atoms
        # and keep any swap that lowers the Kabsch RMSD.
        if ref.shape[0] <= 220:
            try:
                swap_aligned, swap_rmsd, swap_mapping = _calc_local_swap_optimize_mapping(
                    ref_symbols_lc,
                    ref,
                    target,
                    best_mapping,
                )
                if swap_rmsd < best_rmsd - 1e-12:
                    best_aligned, best_rmsd, best_mapping = (
                        swap_aligned,
                        float(swap_rmsd),
                        swap_mapping,
                    )
            except Exception:
                pass

        return best_aligned, best_rmsd, best_mapping

    def _calc_align_reference_to_target(
        ref_symbols,
        ref_coords,
        target_symbols,
        target_coords,
    ):
        if len(ref_symbols) != len(target_symbols):
            raise ValueError(
                f'Atom count mismatch: ref {len(ref_symbols)} vs target {len(target_symbols)}.'
            )

        ref_seq = [s.lower() for s in ref_symbols]
        target_seq = [s.lower() for s in target_symbols]
        ref_counts = Counter(ref_seq)
        target_counts = Counter(target_seq)
        if ref_counts != target_counts:
            missing = sorted((ref_counts - target_counts).items())
            extra = sorted((target_counts - ref_counts).items())
            raise ValueError(
                f'Element composition mismatch. Missing in target: {missing}; extra in target: {extra}.'
            )

        ref = np.asarray(ref_coords, dtype=float)
        target = np.asarray(target_coords, dtype=float)
        n_atoms = ref.shape[0]
        target_centered = target - target.mean(axis=0)

        def _solve_for_reference(ref_work, orientation_label):
            ref_work = np.asarray(ref_work, dtype=float)
            ref_centered = ref_work - ref_work.mean(axis=0)
            candidates = []
            candidate_index_by_key = {}

            def _register_candidate(rmsd, aligned, mapping, source):
                map_idx = np.asarray(mapping, dtype=int)
                if map_idx.size != n_atoms:
                    return
                key = tuple(int(v) for v in map_idx.tolist())
                row = (float(rmsd), np.asarray(aligned, dtype=float), map_idx.copy(), source)
                idx = candidate_index_by_key.get(key)
                if idx is None:
                    candidate_index_by_key[key] = len(candidates)
                    candidates.append(row)
                elif row[0] < candidates[idx][0] - 1e-12:
                    candidates[idx] = row

            def _add_candidate(mapping, source):
                map_idx = np.asarray(mapping, dtype=int)
                if map_idx.size != n_atoms:
                    return
                aligned, rmsd, refined = _calc_refine_alignment_with_mapping(
                    ref_seq, target_seq, ref_work, target, map_idx
                )
                _register_candidate(rmsd, aligned, refined, source)

            # Keep direct order candidate for identical atom order.
            if ref_seq == target_seq:
                _add_candidate(np.arange(n_atoms, dtype=int), 'direct')

            # Best topological alignment via RDKit (symmetry-aware).
            rdkit_best = _calc_rdkit_best_alignment_from_xyz(
                ref_symbols, ref_work, target_symbols, target
            )
            if rdkit_best is not None:
                rd_aligned, rd_rmsd, rd_mapping = rdkit_best
                _register_candidate(rd_rmsd, rd_aligned, rd_mapping, 'rdkit-topology')
            rdkit_best_reverse = _calc_rdkit_best_alignment_from_xyz(
                target_symbols, target, ref_symbols, ref_work
            )
            if rdkit_best_reverse is not None:
                _rd_aligned_rev, rd_rmsd_rev, rd_mapping_rev = rdkit_best_reverse
                inv_map = _calc_invert_mapping(rd_mapping_rev, n_atoms)
                if inv_map is not None:
                    _add_candidate(inv_map, 'rdkit-topology-reverse')

            # Topology-guided candidates from inferred connectivity.
            topo_maps = _calc_topology_mappings_from_xyz(
                ref_symbols, ref_work, target_symbols, target
            )
            if len(topo_maps) > 256:
                scored = []
                for mapping in topo_maps:
                    try:
                        _aligned_quick, quick_rmsd = _calc_kabsch_align(
                            ref_work, target, mapping=mapping
                        )
                        scored.append((float(quick_rmsd), mapping))
                    except Exception:
                        continue
                scored.sort(key=lambda item: item[0])
                topo_maps = [mapping for _rmsd, mapping in scored[:256]]
            for mapping in topo_maps:
                _add_candidate(mapping, 'topology')

            topo_maps_reverse_raw = _calc_topology_mappings_from_xyz(
                target_symbols, target, ref_symbols, ref_work
            )
            topo_maps_reverse = []
            for mapping_rev in topo_maps_reverse_raw:
                inv_map = _calc_invert_mapping(mapping_rev, n_atoms)
                if inv_map is not None:
                    topo_maps_reverse.append(inv_map)
            if len(topo_maps_reverse) > 256:
                scored_rev = []
                for mapping in topo_maps_reverse:
                    try:
                        _aligned_quick, quick_rmsd = _calc_kabsch_align(
                            ref_work, target, mapping=mapping
                        )
                        scored_rev.append((float(quick_rmsd), mapping))
                    except Exception:
                        continue
                scored_rev.sort(key=lambda item: item[0])
                topo_maps_reverse = [mapping for _rmsd, mapping in scored_rev[:256]]
            for mapping in topo_maps_reverse:
                _add_candidate(mapping, 'topology-reverse')

            # Global orientation + Hungarian assignment candidates.
            n_random = 120 if n_atoms <= 40 else (72 if n_atoms <= 120 else 36)
            rotation_guesses = _calc_generate_proper_axis_rotations()
            rotation_guesses.extend(_calc_random_rotation_matrices(n_samples=n_random))
            for rot_guess in rotation_guesses:
                try:
                    initial_map = _calc_element_assignment_for_rotation(
                        ref_seq,
                        target_seq,
                        ref_centered @ rot_guess,
                        target_centered,
                    )
                except RuntimeError:
                    initial_map = None
                if initial_map is not None:
                    _add_candidate(initial_map, 'global-permutation')

            if not candidates:
                aligned, rmsd = _calc_kabsch_align(ref_work, target)
                return aligned, rmsd, {
                    'method': 'direct-fallback',
                    'mapping': list(range(n_atoms)),
                    'orientation': orientation_label,
                }

            best_rmsd, best_aligned, best_mapping, best_source = min(
                candidates, key=lambda item: item[0]
            )
            return best_aligned, best_rmsd, {
                'method': best_source,
                'mapping': [int(v) for v in np.asarray(best_mapping, dtype=int).tolist()],
                'orientation': orientation_label,
            }

        normal_aligned, normal_rmsd, normal_meta = _solve_for_reference(ref, 'normal')
        ref_centroid = ref.mean(axis=0)
        inverted_ref = (2.0 * ref_centroid) - ref
        inverted_aligned, inverted_rmsd, inverted_meta = _solve_for_reference(
            inverted_ref, 'inverted'
        )

        if inverted_rmsd < normal_rmsd - 1e-12:
            meta = dict(inverted_meta)
            meta['rmsd_normal'] = float(normal_rmsd)
            meta['rmsd_inverted'] = float(inverted_rmsd)
            meta['orientation_choice'] = 'inverted'
            return inverted_aligned, inverted_rmsd, meta
        meta = dict(normal_meta)
        meta['rmsd_normal'] = float(normal_rmsd)
        meta['rmsd_inverted'] = float(inverted_rmsd)
        meta['orientation_choice'] = 'normal'
        return normal_aligned, normal_rmsd, meta

    def _render_3dmol_dual_xyz(reference_xyz, target_xyz):
        import json
        _mol3d_counter[0] += 1
        viewer_id = f'mol3d_rmsd_{_mol3d_counter[0]}'
        wrapper_id = f'calc_mol_wrap_{_mol3d_counter[0]}'
        ref_json = json.dumps(reference_xyz)
        target_json = json.dumps(target_xyz)
        view_scope_json = json.dumps(f"{calc_scope_id}:{state.get('current_path') or '/'}")
        scope_id_json = json.dumps(calc_scope_id)

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
                var scopeRoot = el ? el.closest('.{calc_scope_id}') : null;
                if (!scopeRoot) scopeRoot = document.querySelector('.{calc_scope_id}');
                if (!el || typeof $3Dmol === "undefined" || !mv || mv.offsetParent === null) {{
                    tries += 1;
                    if (tries < 80) setTimeout(initViewer, 50);
                    return;
                }}
                var lft = scopeRoot ? scopeRoot.querySelector('.calc-left') : null;
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
                window._calcMolViewerByScope = window._calcMolViewerByScope || {{}};
                window._calcMolViewScopeKeyByScope = window._calcMolViewScopeKeyByScope || {{}};
                window._calcTrajViewerByScope = window._calcTrajViewerByScope || {{}};
                var scopeKey = {scope_id_json};
                var viewScope = {view_scope_json};
                var previousViewer =
                    window._calcMolViewerByScope[scopeKey]
                    || window._calcTrajViewerByScope[scopeKey]
                    || null;
                var previousScope =
                    window._calcMolViewScopeKeyByScope[scopeKey] || viewScope;
                if (previousViewer && typeof previousViewer.getView === 'function') {{
                    try {{
                        window._calcMolViewStateByScope[previousScope] = previousViewer.getView();
                    }} catch (_e) {{}}
                }}
                var savedView = window._calcMolViewStateByScope[viewScope] || null;
                var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                var targetData = {target_json};
                var refData = {ref_json};
                viewer.addModel(targetData, "xyz");
                viewer.setStyle({{model:0}}, {{
                    stick: {{radius: 0.20, color: "#1f5fff"}},
                    sphere: {{scale: 0.24, color: "#1f5fff"}}
                }});
                viewer.addModel(refData, "xyz");
                viewer.setStyle({{model:1}}, {{
                    stick: {{radius: 0.20, color: "#d32f2f"}},
                    sphere: {{scale: 0.24, color: "#d32f2f"}}
                }});
                if (savedView && typeof viewer.setView === 'function') {{
                    try {{
                        viewer.setView(savedView);
                    }} catch (_e) {{
                        viewer.zoomTo();
                        viewer.center();
                        viewer.zoom({CALC_MOL_ZOOM});
                    }}
                }} else {{
                    viewer.zoomTo();
                    viewer.center();
                    viewer.zoom({CALC_MOL_ZOOM});
                }}
                viewer.render();
                window._calcMolViewerByScope[scopeKey] = viewer;
                window._calcMolViewScopeKeyByScope[scopeKey] = viewScope;
                var wrappers = scopeRoot
                    ? scopeRoot.querySelectorAll('.calc-mol-stage-wrapper')
                    : document.querySelectorAll('.calc-mol-stage-wrapper');
                wrappers.forEach(function(w) {{
                    if (w.id !== "{wrapper_id}") w.remove();
                }});
                if (window["{calc_resize_mol_fn}"]) {{
                    setTimeout(window["{calc_resize_mol_fn}"], 200);
                }}
            }}
            setTimeout(initViewer, 0);
        }})();
        </script>
        """
        with calc_mol_viewer:
            display(HTML(html))

    def _render_rmsd_preview_dual_xyz(
        reference_xyz,
        target_xyz,
        output_widget=None,
        viewer_height=None,
    ):
        if output_widget is None:
            output_widget = calc_rmsd_preview
        if not viewer_height:
            viewer_height = CALC_RMSD_PANEL_HEIGHT
        in_main_mol_viewer = output_widget is calc_mol_viewer
        _mol3d_counter[0] += 1
        if in_main_mol_viewer:
            viewer_id = f'3dmolviewer_rmsd_{_mol3d_counter[0]}'
            wrapper_id = f'calc_mol_wrap_rmsd_{_mol3d_counter[0]}'
            viewer_container = (
                f'<div id="{wrapper_id}" class="calc-mol-stage-wrapper" style="width:100%;">'
                f'<div id="{viewer_id}" style="width:100%;height:{viewer_height};position:relative;"></div>'
                '</div>'
            )
        else:
            viewer_id = f'rmsd_preview_{_mol3d_counter[0]}'
            wrapper_id = ''
            viewer_container = (
                f'<div id="{viewer_id}" style="width:100%;height:{viewer_height};position:relative;"></div>'
            )
        target_json = json.dumps(target_xyz)
        ref_json = json.dumps(reference_xyz)
        calc_scope_selector = _html.escape(calc_scope_id, quote=True)
        if in_main_mol_viewer:
            post_render_js = (
                f'var scopeRoot = document.querySelector(".{calc_scope_selector}");'
                'var wrappers = scopeRoot ? scopeRoot.querySelectorAll(".calc-mol-stage-wrapper")'
                ' : document.querySelectorAll(".calc-mol-stage-wrapper");'
                f'wrappers.forEach(function(w) {{ if (w.id !== "{wrapper_id}") w.remove(); }});'
                f'if (window["{calc_resize_mol_fn}"]) {{'
                f' setTimeout(window["{calc_resize_mol_fn}"], 120);'
                f' setTimeout(window["{calc_resize_mol_fn}"], 360);'
                f'}}'
            )
        else:
            post_render_js = ''
        html = f"""
        {viewer_container}
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
                if (!el || typeof $3Dmol === "undefined") {{
                    tries += 1;
                    if (tries < 80) setTimeout(initViewer, 50);
                    return;
                }}
                var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                var targetData = {target_json};
                var refData = {ref_json};
                viewer.addModel(targetData, "xyz");
                viewer.setStyle({{model:0}}, {{
                    stick: {{radius: 0.20, color: "#1f5fff"}},
                    sphere: {{scale: 0.24, color: "#1f5fff"}}
                }});
                viewer.addModel(refData, "xyz");
                viewer.setStyle({{model:1}}, {{
                    stick: {{radius: 0.20, color: "#d32f2f"}},
                    sphere: {{scale: 0.24, color: "#d32f2f"}}
                }});
                viewer.zoomTo();
                viewer.center();
                viewer.zoom({CALC_MOL_ZOOM});
                viewer.render();
                {post_render_js}
            }}
            setTimeout(initViewer, 0);
        }})();
        </script>
        """
        with output_widget:
            clear_output()
            display(HTML(html))

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
        selected_path = _calc_selected_item_path()
        rmsd_enabled = bool(
            selected_path
            and selected_path.suffix.lower() == '.xyz'
            and len(frames) == 1
        )
        _calc_set_rmsd_available(rmsd_enabled)

        # Trajectory: use JS viewer and keep orientation when switching frames
        if len(frames) > 1:
            if initial_load or not state['traj_viewer_ready']:
                full_xyz = ""
                for comm, block, natoms in frames:
                    full_xyz += f"{natoms}\n{comm}\n{block}\n"
                _mol3d_counter[0] += 1
                viewer_id = f"calc_trj_viewer_{_mol3d_counter[0]}"
                wrapper_id = f"calc_mol_wrap_{_mol3d_counter[0]}"
                view_scope_json = json.dumps(
                    f"{calc_scope_id}:{state.get('current_path') or '/'}"
                )
                scope_id_json = json.dumps(calc_scope_id)
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
                            var scopeRoot = el ? el.closest('.{calc_scope_id}') : null;
                            if (!scopeRoot) scopeRoot = document.querySelector('.{calc_scope_id}');
                            if (!el || typeof $3Dmol === "undefined"
                                || !mv || mv.offsetParent === null) {{
                                tries += 1;
                                if (tries < 80) setTimeout(initViewer, 50);
                                return;
                            }}
                            /* Pre-size to correct dimensions */
                            var lft = scopeRoot ? scopeRoot.querySelector('.calc-left') : null;
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
                            window._calcMolViewerByScope = window._calcMolViewerByScope || {{}};
                            window._calcMolViewScopeKeyByScope = window._calcMolViewScopeKeyByScope || {{}};
                            window._calcTrajViewerByScope = window._calcTrajViewerByScope || {{}};
                            var scopeKey = {scope_id_json};
                            var viewScope = {view_scope_json};
                            var previousViewer =
                                window._calcMolViewerByScope[scopeKey]
                                || window._calcTrajViewerByScope[scopeKey]
                                || null;
                            var previousScope =
                                window._calcMolViewScopeKeyByScope[scopeKey] || viewScope;
                            if (previousViewer && typeof previousViewer.getView === 'function') {{
                                try {{
                                    window._calcMolViewStateByScope[previousScope] = previousViewer.getView();
                                }} catch (_e) {{}}
                            }}
                            var savedView = window._calcMolViewStateByScope[viewScope] || null;
                            var viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
                            var xyz = `{full_xyz}`;
                            viewer.addModelsAsFrames(xyz, "xyz");
                            viewer.setStyle({{}}, {DEFAULT_3DMOL_STYLE_JS});
                            if (savedView && typeof viewer.setView === 'function') {{
                                try {{
                                    viewer.setView(savedView);
                                }} catch (_e) {{
                                    viewer.zoomTo();
                                    viewer.center();
                                    viewer.zoom({CALC_MOL_ZOOM});
                                }}
                            }} else {{
                                viewer.zoomTo();
                                viewer.center();
                                viewer.zoom({CALC_MOL_ZOOM});
                            }}
                            viewer.setFrame({idx});
                            viewer.render();
                            window._calcTrajViewerByScope[scopeKey] = viewer;
                            window._calcMolViewerByScope[scopeKey] = viewer;
                            window._calcMolViewScopeKeyByScope[scopeKey] = viewScope;
                            var wrappers = scopeRoot
                                ? scopeRoot.querySelectorAll('.calc-mol-stage-wrapper')
                                : document.querySelectorAll('.calc-mol-stage-wrapper');
                            wrappers.forEach(function(w) {{
                                if (w.id !== "{wrapper_id}") w.remove();
                            }});
                            if (window["{calc_resize_mol_fn}"]) {{
                                setTimeout(window["{calc_resize_mol_fn}"], 200);
                            }}
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
                    var scopeKey = {json.dumps(calc_scope_id)};
                    var v = window._calcTrajViewerByScope
                        ? window._calcTrajViewerByScope[scopeKey]
                        : null;
                    if (v) {{
                        v.setFrame({idx});
                        v.render();
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

    def calc_collect_local_job_status_dirs():
        """Return latest local queue status per top-level calc folder."""
        backend = getattr(ctx, 'backend', None)
        load_jobs = getattr(backend, '_load_jobs', None)
        if not callable(load_jobs):
            return {}

        try:
            data = load_jobs()
        except Exception:
            return {}
        if not isinstance(data, dict):
            return {}

        calc_root = _calc_dir()
        try:
            calc_root_resolved = calc_root.resolve()
        except Exception:
            calc_root_resolved = calc_root

        latest_by_folder = {}
        for job in data.get('jobs', []) or []:
            if not isinstance(job, dict):
                continue
            raw_dir = job.get('job_dir')
            raw_status = str(job.get('status', '')).upper().strip()
            if not raw_dir or not raw_status:
                continue
            try:
                job_dir = Path(raw_dir).resolve()
            except Exception:
                continue
            try:
                rel = job_dir.relative_to(calc_root_resolved)
            except Exception:
                continue
            if not rel.parts:
                continue

            folder_path = calc_root_resolved / rel.parts[0]
            try:
                job_id = int(job.get('job_id', 0))
            except Exception:
                job_id = 0

            prev = latest_by_folder.get(folder_path)
            if prev is None or job_id >= prev[0]:
                latest_by_folder[folder_path] = (job_id, raw_status)

        return {folder: status for folder, (_job_id, status) in latest_by_folder.items()}

    def calc_folder_has_completed_results(folder: Path):
        """Best-effort check whether a folder contains finished calculation outputs."""
        if not folder.is_dir():
            return False

        finished_markers = (
            'DELFIN.txt',
            'DELFIN.docx',
            'DELFIN_Data.json',
            'ESD.txt',
            'report.docx',
        )
        if any((folder / marker).exists() for marker in finished_markers):
            return True

        co2_markers = (
            folder / 'CO2_coordination' / 'orientation_scan' / 'orientation_scan.csv',
            folder / 'CO2_coordination' / 'relaxed_surface_scan' / 'scan.relaxscanact.dat',
            folder / 'CO2_coordination' / 'relaxed_surface_scan' / 'scan.out',
        )
        if any(marker.exists() for marker in co2_markers):
            return True

        for scan_dir in (
            folder,
            folder / 'builder',
            folder / 'CO2_coordination' / 'relaxed_surface_scan',
        ):
            try:
                if not scan_dir.exists():
                    continue
                for out_path in scan_dir.glob('*.out'):
                    if out_path.name.lower().startswith('delfin_'):
                        continue
                    if calc_orca_terminated_normally(out_path):
                        return True
            except Exception:
                continue

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
        state['chunk_dom_initialized'] = False
        calc_file_list.options = []
        state['traj_viewer_ready'] = False
        _calc_clear_main_viewer_state(reset_view_state=False)
        calc_mol_container.layout.display = 'none'
        _calc_set_rmsd_available(False)
        _calc_reset_rmsd_ui(clear_input=True)
        calc_file_info.value = ''
        calc_path_display.value = ''
        calc_download_status.value = ''
        calc_search_result.value = ''
        calc_search_suggest.value = '(Select)'
        calc_search_input.value = ''
        calc_folder_search.value = ''
        calc_prev_btn.disabled = True
        calc_next_btn.disabled = True
        _set_view_toggle(False, True)
        calc_copy_btn.disabled = True
        calc_copy_path_btn.disabled = True
        calc_download_btn.disabled = True
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
            try:
                is_top_level_calc_view = current_dir.resolve() == _calc_dir().resolve()
            except Exception:
                is_top_level_calc_view = (state.get('current_path') or '') == ''
            local_status_dirs = calc_collect_local_job_status_dirs()
            for entry in entries:
                if entry.is_dir():
                    if is_top_level_calc_view:
                        try:
                            entry_resolved = entry.resolve()
                        except Exception:
                            entry_resolved = entry
                        if entry_resolved in local_status_dirs:
                            status = local_status_dirs.get(entry_resolved, '')
                            folder_icon = '✅' if status == 'COMPLETED' else '📂'
                        elif calc_folder_has_completed_results(entry):
                            folder_icon = '✅'
                        else:
                            folder_icon = '📂'
                    else:
                        folder_icon = '📂'
                    items.append(f'{folder_icon} {entry.name}')
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
        calc_multi_select.options = state['all_items']
        calc_update_path_display()
        calc_update_report_btn()
        calc_update_download_btn()

    def calc_filter_file_list(change=None):
        query = calc_folder_search.value.strip().lower()
        if not query:
            calc_file_list.options = state['all_items']
            calc_multi_select.options = state['all_items']
        else:
            filtered = [item for item in state['all_items'] if query in item.lower()]
            result = filtered if filtered else ['(No matches)']
            calc_file_list.options = result
            calc_multi_select.options = result

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
                _calc_set_requested_ratio_for_offset(start)
                _run_js(
                    "window.__calcChunkPendingRatio = null;"
                    "window.__calcChunkBusy = false;"
                    "window.__calcChunkIgnoreScrollUntil = Date.now() + 120;"
                    "window.__calcChunkProgrammaticScroll = false;"
                )
                desired = max(0, int(start) - (CALC_TEXT_CHUNK_BYTES // 4))
                _calc_request_chunk_start(desired)
                chunk_start = int(state.get('file_chunk_start') or 0)
                chunk_end = int(state.get('file_chunk_end') or 0)
                if start < chunk_start or end > chunk_end:
                    _calc_set_requested_ratio_for_offset(start)
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

    # -- options dropdown ---------------------------------------------------
    def calc_update_options_dropdown():
        selected = calc_file_list.value
        sel_lower = selected.lower() if selected else ''
        rmsd_available = bool(state.get('rmsd_available'))
        if selected:
            name = _calc_label_to_name(selected)
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
        elif selected and rmsd_available:
            calc_options_dropdown.options = ['(Options)', 'RMSD']
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
            calc_file_list.value = None

    def calc_on_home(b):
        state['current_path'] = ''
        calc_list_directory()
        calc_file_list.value = None

    def calc_on_refresh(b):
        calc_list_directory()
        calc_file_list.value = None

    def calc_go_top(b):
        if not state['file_content']:
            return
        if _calc_is_chunk_mode():
            _run_js(
                "window.__calcChunkRequestedRatio = null;"
                "window.__calcChunkPendingRatio = null;"
                "window.__calcChunkBusy = false;"
                "window.__calcChunkProgrammaticScroll = false;"
            )
            _calc_request_chunk_start(0, scroll_to='top')
            return
        calc_scroll_to('top')

    def calc_go_bottom(b):
        if not state['file_content']:
            return
        if _calc_is_chunk_mode():
            size = int(state.get('selected_file_size') or 0)
            # Load the true file tail (last chunk-sized window up to EOF).
            target = max(0, size - CALC_TEXT_CHUNK_BYTES)
            _run_js(
                "window.__calcChunkRequestedRatio = null;"
                "window.__calcChunkPendingRatio = null;"
                "window.__calcChunkBusy = false;"
                "window.__calcChunkProgrammaticScroll = false;"
            )
            _calc_request_chunk_start(target, scroll_to='bottom', align_to_chunk=False)
            return
        calc_scroll_to('bottom')

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
                "window.__calcChunkIgnoreScrollUntil = Date.now() + 120;"
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

        q_lower = query.lower()
        q_len = len(query)
        spans = []
        hay = state['file_content'].lower()
        pos = hay.find(q_lower)
        while pos >= 0:
            spans.append((pos, pos + q_len))
            if len(spans) >= CALC_SEARCH_MAX_MATCHES:
                state['search_truncated'] = True
                break
            pos = hay.find(q_lower, pos + 1)
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

    def calc_on_search_submit(_widget):
        calc_do_search()

    def _calc_request_chunk_start(requested, scroll_to=None, align_to_chunk=True):
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
        if align_to_chunk:
            tail_start = max(0, size - CALC_TEXT_CHUNK_BYTES)
            if req >= tail_start:
                req = tail_start
            else:
                req = (req // CALC_TEXT_CHUNK_BYTES) * CALC_TEXT_CHUNK_BYTES
        else:
            req = min(req, max(0, size - CALC_TEXT_CHUNK_BYTES))
        current_start = int(state.get('file_chunk_start') or 0)
        if scroll_to == 'top':
            _run_js("window.__calcChunkRequestedRatio = 0.0;")
        elif scroll_to == 'bottom':
            _run_js("window.__calcChunkRequestedRatio = 1.0;")
        if req == current_start:
            if scroll_to:
                calc_scroll_to(scroll_to)
            return False
        _calc_load_text_preview_chunk(path, size, req, scroll_to=scroll_to)
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

    def calc_on_rmsd_toggle(button=None):
        if not state.get('rmsd_available'):
            return
        if calc_rmsd_controls.layout.display == 'none':
            _calc_set_rmsd_mode(True)
            calc_rmsd_controls.layout.display = 'flex'
            calc_rmsd_hide_btn.layout.display = 'inline-flex'
            calc_rmsd_status.value = (
                '<span style="color:#555;">Paste reference coordinates and click RMSD.</span>'
            )
            with calc_rmsd_preview:
                clear_output()
                display(HTML(
                    f"<div style='height:{CALC_RMSD_PANEL_HEIGHT};display:flex;align-items:center;justify-content:center;"
                    "color:#666;border:1px dashed #bbb;background:#fafafa;'>"
                    "Visualization appears here after RMSD calculation."
                    "</div>"
                ))
        else:
            _calc_reset_rmsd_ui(clear_input=False)

    def calc_on_rmsd_hide(button):
        _calc_reset_rmsd_ui(clear_input=False)

    def calc_on_rmsd_run(button):
        selected_path = _calc_selected_item_path()
        if not selected_path or selected_path.suffix.lower() != '.xyz':
            calc_rmsd_status.value = (
                '<span style="color:#d32f2f;">Select a single-frame .xyz file first.</span>'
            )
            return

        try:
            target_symbols, target_coords, target_comment, target_xyz = _calc_get_current_xyz_model()
            ref_symbols, ref_coords, ref_comment = _calc_parse_reference_xyz_input(
                calc_rmsd_ref_input.value
            )
        except Exception as exc:
            calc_rmsd_status.value = (
                f'<span style="color:#d32f2f;">Input error: {_html.escape(str(exc))}</span>'
            )
            return

        try:
            aligned_ref_coords, rmsd_value, align_meta = _calc_align_reference_to_target(
                ref_symbols, ref_coords, target_symbols, target_coords
            )
            align_method = align_meta.get('method', 'direct')
            align_orientation = align_meta.get('orientation', 'normal')
            rmsd_normal = align_meta.get('rmsd_normal', None)
            rmsd_inverted = align_meta.get('rmsd_inverted', None)

            def _one_line_comment(text, limit):
                cleaned = ' '.join((text or '').split()).strip()
                if len(cleaned) > limit:
                    cleaned = cleaned[:limit].rstrip()
                return cleaned

            target_comment_text = _one_line_comment(target_comment, 180)
            ref_comment_text = _one_line_comment(ref_comment, 160)
            info_text = ' '.join((calc_rmsd_info_input.value or '').split()).strip()
            if len(info_text) > 240:
                info_text = info_text[:240].rstrip()
            aligned_comment = (
                f'Aligned reference to {selected_path.name} '
                f'(RMSD={rmsd_value:.6f} Angstrom)'
            )
            aligned_comment = f'{aligned_comment} | Orientation: {align_orientation}'
            if rmsd_normal is not None and rmsd_inverted is not None:
                aligned_comment = (
                    f'{aligned_comment} | Compare(normal={float(rmsd_normal):.6f}, '
                    f'inverted={float(rmsd_inverted):.6f})'
                )
            if target_comment_text:
                aligned_comment = f'{aligned_comment} | TargetComment: {target_comment_text}'
            if ref_comment_text:
                aligned_comment = f'{aligned_comment} | RefComment: {ref_comment_text}'
            if align_method not in ('direct', 'direct-fallback'):
                aligned_comment = f'{aligned_comment} | Mapping: {align_method}'
            if info_text:
                aligned_comment = f'{aligned_comment} | Info: {info_text}'
            aligned_ref_xyz = _calc_build_xyz_from_symbols_coords(
                ref_symbols, aligned_ref_coords, comment=aligned_comment
            )
            base_stem = re.sub(r'(?i)_aligned_ref(?:_\d+)?$', '', selected_path.stem)
            aligned_base = f'{base_stem}_aligned_ref'
            output_path = selected_path.with_name(f'{aligned_base}.xyz')
            if output_path.exists():
                idx = 2
                while True:
                    candidate = selected_path.with_name(f'{aligned_base}_{idx}.xyz')
                    if not candidate.exists():
                        output_path = candidate
                        break
                    idx += 1
            output_path.write_text(aligned_ref_xyz, encoding='utf-8')
            _render_rmsd_preview_dual_xyz(aligned_ref_xyz, target_xyz)
            compare_text = ''
            if rmsd_normal is not None and rmsd_inverted is not None:
                compare_text = (
                    f'Compare(normal={float(rmsd_normal):.6f}, '
                    f'inverted={float(rmsd_inverted):.6f}). '
                )
            calc_rmsd_status.value = (
                '<span style="color:#2e7d32;">'
                f'RMSD = {rmsd_value:.6f} Angstrom. '
                f'Saved: {_html.escape(output_path.name)} '
                f'Orientation: {_html.escape(align_orientation)}. '
                f'Alignment: {_html.escape(align_method)}. '
                f'{compare_text}'
                f'(ref comment: {_html.escape(ref_comment[:80])}, '
                f'target: {_html.escape(target_comment[:80])}).'
                '</span>'
            )
        except Exception as exc:
            calc_rmsd_status.value = (
                f'<span style="color:#d32f2f;">RMSD failed: {_html.escape(str(exc))}</span>'
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
        full_path_obj = _calc_selected_item_path()
        if full_path_obj is None and state['current_path']:
            full_path_obj = _calc_dir() / state['current_path']
        if full_path_obj is None:
            return
        full_path = str(full_path_obj)
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

    def _calc_trigger_download(filename, payload, mime='application/octet-stream'):
        b64 = base64.b64encode(payload).decode('ascii')
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

    def calc_on_download(button):
        target = _calc_download_target_path()
        if target is None:
            calc_download_status.value = '<span style="color:#d32f2f;">No file/folder selected.</span>'
            return

        try:
            payload = b''
            filename = ''
            mime = 'application/octet-stream'

            if target.is_dir():
                calc_download_status.value = (
                    f'<span style="color:#1976d2;">Preparing ZIP for '
                    f'{_html.escape(target.name or str(target))}...</span>'
                )
                with tempfile.TemporaryDirectory(prefix='delfin_download_') as tmpdir:
                    archive_name = target.name or 'folder'
                    archive_base = Path(tmpdir) / archive_name
                    archive_path = Path(shutil.make_archive(
                        str(archive_base),
                        'zip',
                        root_dir=str(target.parent),
                        base_dir=target.name,
                    ))
                    payload = archive_path.read_bytes()
                    filename = archive_path.name
                    mime = 'application/zip'
            else:
                payload = target.read_bytes()
                filename = target.name

            size_bytes = len(payload)
            if size_bytes > CALC_DOWNLOAD_MAX_BYTES:
                calc_download_status.value = (
                    '<span style="color:#d32f2f;">'
                    f'Download too large ({_calc_format_bytes(size_bytes)}). '
                    f'Limit: {_calc_format_bytes(CALC_DOWNLOAD_MAX_BYTES)}.'
                    '</span>'
                )
                return

            _calc_trigger_download(filename, payload, mime=mime)
            calc_download_status.value = (
                f'<span style="color:#2e7d32;">Download started: '
                f'{_html.escape(filename)} ({_calc_format_bytes(size_bytes)})</span>'
            )
        except Exception as exc:
            calc_download_status.value = (
                f'<span style="color:#d32f2f;">Download failed: '
                f'{_html.escape(str(exc))}</span>'
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
                name = _calc_label_to_name(selected)
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
        elif change['new'] == 'RMSD':
            calc_override_input.layout.display = 'none'
            calc_override_time.layout.display = 'none'
            calc_override_btn.layout.display = 'none'
            calc_override_status.layout.display = 'none'
            calc_override_status.value = ''
            calc_edit_area.layout.display = 'none'
            _calc_preselect_show(False)
            if not state.get('rmsd_available'):
                try:
                    calc_options_dropdown.unobserve(calc_on_options_change, names='value')
                except Exception:
                    pass
                calc_options_dropdown.value = '(Options)'
                calc_options_dropdown.observe(calc_on_options_change, names='value')
                return
            _set_view_toggle(True, False)
            calc_update_view()
            calc_on_rmsd_toggle()
            try:
                calc_options_dropdown.unobserve(calc_on_options_change, names='value')
            except Exception:
                pass
            calc_options_dropdown.value = '(Options)'
            calc_options_dropdown.observe(calc_on_options_change, names='value')
        else:
            if state.get('rmsd_mode_active'):
                return
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
        selected_name = _calc_label_to_name(selected) if selected else ''
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

    def calc_on_select_toggle(button):
        state['select_mode'] = not state['select_mode']
        if state['select_mode']:
            calc_select_btn.button_style = 'success'
            calc_file_list.layout.display = 'none'
            calc_multi_select.layout.display = ''
            calc_multi_select.value = ()
        else:
            calc_select_btn.button_style = ''
            calc_multi_select.layout.display = 'none'
            calc_file_list.layout.display = ''
            calc_multi_select.value = ()
            calc_file_list.value = None
            calc_delete_confirm.layout.display = 'none'
            calc_move_archive_confirm.layout.display = 'none'

    def calc_on_delete_click(button):
        if state['select_mode']:
            real = [s for s in calc_multi_select.value if not s.startswith('(')]
            if not real:
                calc_delete_status.value = '<span style="color:#d32f2f;">No items selected.</span>'
                return
            names = [_calc_label_to_name(s) for s in real]
            label_str = ', '.join(f'<code>{_html.escape(n)}</code>' for n in names)
            calc_delete_label.value = f'<b>Delete {len(names)} item(s)?</b> {label_str}'
            state['delete_current'] = False
            calc_delete_confirm.layout.display = 'flex'
            return
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
        if state['select_mode']:
            real = [s for s in calc_multi_select.value if not s.startswith('(')]
            if not real:
                calc_delete_status.value = '<span style="color:#d32f2f;">Nothing to delete.</span>'
                calc_delete_confirm.layout.display = 'none'
                return
            errors, deleted = [], []
            for label in real:
                name = _calc_label_to_name(label)
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
                    deleted.append(name)
                except Exception as exc:
                    errors.append(f'{_html.escape(name)}: {_html.escape(str(exc))}')
            calc_delete_confirm.layout.display = 'none'
            if deleted and not errors:
                calc_delete_status.value = (
                    f'<span style="color:#2e7d32;">Deleted {len(deleted)} item(s).</span>'
                )
            elif deleted and errors:
                calc_delete_status.value = (
                    f'<span style="color:#e65100;">Deleted {len(deleted)}, '
                    f'failed: {"; ".join(errors)}</span>'
                )
            else:
                calc_delete_status.value = (
                    f'<span style="color:#d32f2f;">{"; ".join(errors)}</span>'
                )
            calc_multi_select.value = ()
            saved_filter = calc_folder_search.value
            calc_list_directory()
            if saved_filter:
                calc_folder_search.value = saved_filter
                calc_filter_file_list()
            return
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
            name = _calc_label_to_name(selected)
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

    def calc_on_move_archive_click(button):
        if _is_archive_tab:
            return
        if state['select_mode']:
            real = [s for s in calc_multi_select.value if not s.startswith('(')]
        else:
            sel = calc_file_list.value
            real = [sel] if sel and not sel.startswith('(') else []
        if not real:
            calc_move_archive_status.value = '<span style="color:#d32f2f;">No items selected.</span>'
            return
        names = [_calc_label_to_name(s) for s in real]
        label_str = ', '.join(f'<code>{_html.escape(n)}</code>' for n in names)
        calc_move_archive_label.value = (
            f'<b>Move {len(names)} item(s) to Archive?</b> {label_str}'
        )
        calc_move_archive_confirm.layout.display = 'flex'

    def calc_on_move_archive_yes(button):
        if state['select_mode']:
            real = [s for s in calc_multi_select.value if not s.startswith('(')]
        else:
            sel = calc_file_list.value
            real = [sel] if sel and not sel.startswith('(') else []
        calc_move_archive_confirm.layout.display = 'none'
        if not real:
            return
        dest_dir = ctx.archive_dir
        dest_dir.mkdir(parents=True, exist_ok=True)
        errors, moved = [], []
        for label in real:
            name = _calc_label_to_name(label)
            src = (
                (_calc_dir() / state['current_path'] / name)
                if state['current_path']
                else (_calc_dir() / name)
            )
            dst = dest_dir / name
            try:
                shutil.move(str(src), str(dst))
                moved.append(name)
            except Exception as exc:
                errors.append(f'{_html.escape(name)}: {_html.escape(str(exc))}')
        if moved and not errors:
            calc_move_archive_status.value = (
                f'<span style="color:#2e7d32;">Moved {len(moved)} item(s) to Archive.</span>'
            )
        elif moved and errors:
            calc_move_archive_status.value = (
                f'<span style="color:#e65100;">Moved {len(moved)}, '
                f'failed: {"; ".join(errors)}</span>'
            )
        else:
            calc_move_archive_status.value = (
                f'<span style="color:#d32f2f;">{"; ".join(errors)}</span>'
            )
        if state['select_mode']:
            calc_multi_select.value = ()
        else:
            calc_file_list.value = None
        saved_filter = calc_folder_search.value
        calc_list_directory()
        if saved_filter:
            calc_folder_search.value = saved_filter
            calc_filter_file_list()

    def calc_on_move_archive_no(button):
        calc_move_archive_confirm.layout.display = 'none'

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

        name = _calc_label_to_name(selected)
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
            calc_file_list.value = None
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
            next_name_lower == 'coord' or next_suffix in ['.cube', '.cub']
        )

        # Reset display (avoid flashing text area before viewer is ready)
        if not keep_previous_viewer_during_load:
            _calc_clear_main_viewer_state(reset_view_state=False)
            calc_mol_container.layout.display = 'none'
        calc_content_area.layout.display = 'none'
        calc_content_label.layout.display = 'none'
        if not keep_previous_viewer_during_load:
            _set_view_toggle(False, True)
        calc_file_info.value = ''
        calc_path_display.value = ''
        calc_download_status.value = ''
        state['file_content'] = ''
        state['selected_file_path'] = None
        state['selected_file_size'] = 0
        state['file_is_preview'] = False
        state['file_preview_note'] = ''
        state['file_chunk_start'] = 0
        state['file_chunk_end'] = 0
        state['chunk_dom_initialized'] = False
        _calc_hide_chunk_controls()
        calc_reset_recalc_state()
        _calc_preselect_show(False)

        # Reset xyz trajectory controls
        state['xyz_frames'].clear()
        state['xyz_current_frame'][0] = 0
        state['traj_viewer_ready'] = False
        calc_xyz_controls.layout.display = 'none'
        calc_coord_controls.layout.display = 'none'
        _calc_set_rmsd_available(False)
        _calc_reset_rmsd_ui(clear_input=True)
        calc_xyz_frame_label.value = ''
        calc_xyz_frame_input.value = 1
        calc_xyz_frame_input.max = 1
        calc_xyz_frame_total.value = '<b>/ 1</b>'

        if not full_path.exists():
            calc_update_download_btn()
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
        calc_update_download_btn()

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
                state['chunk_dom_initialized'] = False
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
                    _render_3dmol(xyz_data)
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
                state['chunk_dom_initialized'] = False
                calc_render_content(scroll_to='top')
                calc_copy_btn.disabled = False
                calc_copy_path_btn.disabled = False
                lines = content.split('\n')

                state['xyz_frames'].clear()
                state['xyz_frames'].extend(parse_xyz_frames(content))
                n_frames = len(state['xyz_frames'])
                is_aligned_overlay_file = bool(
                    re.search(r'(?i)_aligned_ref(?:_\d+)?\.xyz$', name)
                )

                if is_aligned_overlay_file and state['xyz_frames']:
                    def _frame_to_xyz(frame):
                        comment, xyz_block, n_atoms = frame
                        return f"{n_atoms}\n{comment}\n{xyz_block}\n"

                    aligned_xyz = _frame_to_xyz(state['xyz_frames'][-1])
                    original_xyz = None
                    original_name = re.sub(
                        r'(?i)_aligned_ref(?:_\d+)?\.xyz$',
                        '.xyz',
                        name,
                    )
                    original_path = full_path.with_name(original_name)
                    if original_path.exists():
                        try:
                            original_content = original_path.read_text()
                            original_frames = parse_xyz_frames(original_content)
                            if original_frames:
                                original_xyz = _frame_to_xyz(original_frames[0])
                        except Exception:
                            original_xyz = None
                    # Fallback for legacy aligned files that already contain both frames.
                    if original_xyz is None and len(state['xyz_frames']) > 1:
                        original_xyz = _frame_to_xyz(state['xyz_frames'][0])

                    if original_xyz is not None:
                        calc_file_info.value = (
                            f'<b><span style="word-break:break-all;">{_html.escape(name)}</span></b>'
                            f' ({size_str}, aligned overlay with '
                            f'<code>{_html.escape(original_name)}</code>)'
                        )
                        state['xyz_current_frame'][0] = 0
                        calc_xyz_controls.layout.display = 'none'
                        calc_coord_controls.layout.display = 'none'
                        _calc_set_rmsd_available(False)
                        calc_xyz_frame_label.value = ''
                        _render_rmsd_preview_dual_xyz(
                            aligned_xyz,
                            original_xyz,
                            output_widget=calc_mol_viewer,
                            viewer_height=f'{CALC_MOL_SIZE}px',
                        )
                        calc_render_content(scroll_to='top')
                        return

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
                    _calc_set_rmsd_available(False)
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
                        _calc_set_rmsd_available(True)
                    else:
                        calc_xyz_controls.layout.display = 'none'
                        _calc_set_rmsd_available(False)
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
                _run_js(
                    f"setTimeout(function(){{"
                    f" if(window['{calc_resize_mol_fn}']) window['{calc_resize_mol_fn}']();"
                    f"}}, 600);"
                )
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
                state['chunk_dom_initialized'] = False

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
            state['chunk_dom_initialized'] = False
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
    if hasattr(calc_search_input, 'on_submit'):
        try:
            calc_search_input.on_submit(calc_on_search_submit)
        except Exception:
            pass
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
    calc_select_btn.on_click(calc_on_select_toggle)
    calc_move_archive_btn.on_click(calc_on_move_archive_click)
    calc_move_archive_yes_btn.on_click(calc_on_move_archive_yes)
    calc_move_archive_no_btn.on_click(calc_on_move_archive_no)
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
    calc_rmsd_run_btn.on_click(calc_on_rmsd_run)
    calc_rmsd_hide_btn.on_click(calc_on_rmsd_hide)
    calc_coord_copy_btn.on_click(calc_on_coord_copy)
    calc_copy_btn.on_click(calc_on_content_copy)
    calc_copy_path_btn.on_click(calc_on_path_copy)
    calc_download_btn.on_click(calc_on_download)
    calc_report_btn.on_click(calc_on_report)
    disable_spellcheck(ctx, class_name='calc-search-input')

    # -- initialise ---------------------------------------------------------
    if ctx.calc_dir.exists():
        calc_list_directory()
    else:
        calc_file_list.options = ['(calc folder not found)']

    # Inject JS once: single-click toggles selection, Shift+click selects range.
    # A global flag prevents re-installation on every tab re-render.
    _multi_select_js = (
        '(function(){'
        'if(window._delfinMultiSelectReady)return;'
        'window._delfinMultiSelectReady=true;'
        'var lastIdx=-1;'
        'function install(sel){'
        'if(sel._delfinToggle)return;sel._delfinToggle=true;'
        'sel.addEventListener("mousedown",function(e){'
        'if(e.target.tagName!=="OPTION")return;'
        'e.preventDefault();'
        'var opts=Array.from(sel.options);'
        'var idx=opts.indexOf(e.target);'
        'if(e.shiftKey&&lastIdx>=0){'
        'var lo=Math.min(lastIdx,idx),hi=Math.max(lastIdx,idx);'
        'for(var i=lo;i<=hi;i++)opts[i].selected=true;'
        '}else{e.target.selected=!e.target.selected;lastIdx=idx;}'
        'sel.dispatchEvent(new Event("change",{bubbles:true}));'
        'sel.focus();'
        '},true);}'
        'function scanAndInstall(root){'
        '(root.querySelectorAll?root.querySelectorAll(".delfin-multi-select select"):[]).forEach(install);}'
        'scanAndInstall(document.body);'
        'new MutationObserver(function(ms){'
        'ms.forEach(function(m){'
        'm.addedNodes.forEach(function(n){if(n.nodeType===1)scanAndInstall(n);});'
        '});}).observe(document.body,{childList:true,subtree:true});'
        '})();'
    )
    from IPython.display import Javascript as _Javascript
    with ctx.js_output:
        display(_Javascript(_multi_select_js))

    # -- layout -------------------------------------------------------------
    _archive_row_children = [] if _is_archive_tab else [calc_move_archive_btn]
    calc_nav_bar = widgets.VBox([
        calc_path_label,
        widgets.HBox(
            [calc_back_btn, calc_home_btn, calc_refresh_btn, calc_delete_btn, calc_select_btn],
            layout=widgets.Layout(
                width='100%', overflow_x='hidden',
                justify_content='flex-start', gap='6px',
            ),
        ),
        widgets.HBox(
            _archive_row_children,
            layout=widgets.Layout(
                width='100%', overflow_x='hidden',
                justify_content='flex-start', gap='6px',
                display='none' if _is_archive_tab else 'flex',
            ),
        ),
    ], layout=widgets.Layout(width='100%', overflow_x='hidden'))

    calc_left = widgets.VBox([
        calc_nav_bar,
        calc_filter_sort_row,
        calc_file_list,
        calc_multi_select,
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
        '.calc-left .widget-select-multiple { flex:1 1 0 !important; min-height:0 !important; }'
        '.calc-left .widget-select-multiple select { height:100% !important; }'
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
        '.calc-rmsd-trigger-btn button { display:flex !important; align-items:center !important;'
        ' justify-content:center !important; text-align:center !important; height:32px !important; }'
        '.calc-rmsd-trigger-btn button > span { display:inline-flex !important; align-items:center !important;'
        ' justify-content:center !important; text-align:center !important; width:100% !important; line-height:1.2 !important; }'
        '.calc-rmsd-controls { width:100% !important;'
        ' flex-wrap:nowrap !important; align-items:stretch !important; }'
        '.calc-rmsd-controls > .widget-vbox, .calc-rmsd-controls > .widget-output'
        ' { flex:1 1 0 !important; min-width:0 !important; }'
        '.calc-rmsd-input-col .widget-textarea { flex:1 1 0 !important; min-height:0 !important; }'
        '.calc-rmsd-input-col .widget-textarea textarea { height:100% !important; min-height:220px !important; }'
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
        '.calc-rmsd-preview { overflow:hidden !important; padding:0 !important; }'
        '.calc-rmsd-preview .output_area, .calc-rmsd-preview .output_subarea,'
        ' .calc-rmsd-preview .output_wrapper, .calc-rmsd-preview .jp-OutputArea-child,'
        ' .calc-rmsd-preview .jp-OutputArea-output'
        ' { padding:0 !important; margin:0 !important; width:100% !important; height:100% !important; }'
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
                [calc_copy_path_btn, calc_copy_btn, calc_download_btn, calc_report_btn, calc_view_toggle],
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
        calc_download_status,
        calc_report_status,
        calc_mol_container,
        calc_content_label,
        calc_preselect_container,
        calc_delete_confirm,
        calc_delete_status,
        calc_move_archive_confirm,
        calc_move_archive_status,
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
    tab_widget.add_class(calc_scope_id)
    calc_left.add_class('calc-left')
    calc_right.add_class('calc-right')
    calc_content_area.add_class('calc-content-area')

    # Enable draggable splitter + dynamic mol-viewer resize + $3Dmol patch
    # Stored as a plain string; the CALLER (create_dashboard in __init__.py)
    # runs ALL tab init scripts in one ctx.run_js() call so that no tab's
    # clear_output() wipes another tab's init JS.
    _init_js = f"""
    (function() {{
        function initCalcScopeBind(attempt) {{
            var scopeKey = {json.dumps(calc_scope_id)};
            var root = document.querySelector('.{calc_scope_id}');
            if (!root) {{
                if ((attempt || 0) < 40) {{
                    setTimeout(function() {{
                        initCalcScopeBind((attempt || 0) + 1);
                    }}, 100);
                    return;
                }}
                // Fallback to old global behavior if scoped root is unavailable.
                root = document.querySelector('.calc-tab');
                if (!root) return;
            }}

        /* --- Splitter drag logic --- */
        var left = root.querySelector('.calc-left');
        var right = root.querySelector('.calc-right');
        var splitter = root.querySelector('.calc-splitter');
        if (left && right && splitter && splitter.dataset.bound !== scopeKey) {{
            splitter.dataset.bound = scopeKey;
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
                if (window["{calc_resize_mol_fn}"]) {{
                    setTimeout(window["{calc_resize_mol_fn}"], 50);
                }}
            }}
            splitter.addEventListener('mousedown', function(e) {{
                e.preventDefault();
                document.addEventListener('mousemove', onMove);
                document.addEventListener('mouseup', onUp);
            }});
        }}

        /* --- Monkey-patch $3Dmol.createViewer to trigger all scoped resizers --- */
        function patchCreateViewer() {{
            if (typeof $3Dmol !== 'undefined' && !$3Dmol._calcResizePatched) {{
                var orig = $3Dmol.createViewer;
                $3Dmol.createViewer = function() {{
                    var v = orig.apply(this, arguments);
                    setTimeout(function() {{
                        try {{
                            var fns = window._calcResizeMolViewerFns || {{}};
                            Object.keys(fns).forEach(function(k) {{
                                var fn = fns[k];
                                if (typeof fn === 'function') fn();
                            }});
                        }} catch (_e) {{}}
                    }}, 300);
                    return v;
                }};
                $3Dmol._calcResizePatched = true;
            }}
        }}
        patchCreateViewer();
        var patchInterval = setInterval(function() {{
            patchCreateViewer();
            if (typeof $3Dmol !== 'undefined' && $3Dmol._calcResizePatched) {{
                clearInterval(patchInterval);
            }}
        }}, 500);

        /* --- Dynamic mol-viewer resize --- */
        window["{calc_resize_mol_fn}"] = function() {{
            var scopeRoot = document.querySelector('.{calc_scope_id}');
            if (!scopeRoot || scopeRoot.offsetParent === null) return;
            var lft = scopeRoot.querySelector('.calc-left');
            var mv = scopeRoot.querySelector('.calc-mol-viewer');
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
            var v = (window._calcMolViewerByScope && window._calcMolViewerByScope[scopeKey])
                || (window._calcTrajViewerByScope && window._calcTrajViewerByScope[scopeKey])
                || null;
            if (v && typeof v.resize === 'function') {{ v.resize(); v.render(); }}
        }};
        window._calcResizeMolViewerFns = window._calcResizeMolViewerFns || {{}};
        window._calcResizeMolViewerFns[scopeKey] = window["{calc_resize_mol_fn}"];
        window.addEventListener('resize', function() {{
            if (window["{calc_resize_mol_fn}"]) {{
                setTimeout(window["{calc_resize_mol_fn}"], 100);
            }}
        }});
        window["{calc_resize_pre_fn}"] = function() {{
            var v = window._calcPreselect3dViewerByScope
                ? window._calcPreselect3dViewerByScope[scopeKey]
                : null;
            var elId = window._calcPreselect3dElIdByScope
                ? window._calcPreselect3dElIdByScope[scopeKey]
                : null;
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
                if (window["{calc_resize_pre_fn}"]) window["{calc_resize_pre_fn}"]();
            }}, 120);
        }});
        if (root) {{
            new MutationObserver(function() {{
                if (window["{calc_resize_mol_fn}"]) {{
                    setTimeout(window["{calc_resize_mol_fn}"], 200);
                }}
                setTimeout(function() {{
                    if (window["{calc_resize_pre_fn}"]) window["{calc_resize_pre_fn}"]();
                }}, 220);
            }}).observe(root, {{attributes: true, subtree: true, attributeFilter: ['style']}});
        }}
        }}
        initCalcScopeBind(0);
    }})();
    """

    def calc_set_root(root_dir):
        """Switch browser root directory and reset to top level."""
        try:
            new_root = Path(root_dir).expanduser()
        except Exception:
            new_root = Path(root_dir)
        new_root.mkdir(parents=True, exist_ok=True)
        ctx.calc_dir = new_root
        state['current_path'] = ''
        calc_list_directory()

    return tab_widget, {
        'calc_list_directory': calc_list_directory,
        'calc_set_root': calc_set_root,
        'init_js': _init_js,
    }
